/*
 * Bigmatrix.cc
 *
 *  Created on: Dec 24, 2014
 *      Author: mimitantono
 */

#include "cluster.h"

Bigmatrix::Bigmatrix() {
	targetampliconids = new unsigned long*[Property::threads];
	scanner = new class scanner[Property::threads];
	for (int i = 0; i < Property::threads; i++) {
		scanner[i].search_begin();
		targetampliconids[i] = new unsigned long[Property::db_data.sequences / 4];
	}
	matrix_x = new std::vector<unsigned long int>[Property::threads];
	matrix_y = new std::vector<unsigned long int>[Property::threads];
	total_match = 0;
	total_qgram = 0;
	total_scan = 0;
	temp_cleaned = 0;
	temp_written = 0;
	total_data = 0;
	pthread_mutex_init(&workmutex, NULL);
}

Bigmatrix::~Bigmatrix() {
	delete[] scanner;
	scanner = NULL;
	delete[] matrix_x;
	matrix_x = NULL;
	delete[] matrix_y;
	matrix_y = NULL;
	for (int i = 0; i < Property::threads; i++) {
		delete[] targetampliconids[i];
	}
	delete[] targetampliconids;

	pthread_mutex_destroy(&workmutex);
}

void Bigmatrix::calculate_partition(int thread_id, int total_thread) {
	if (thread_id == 0)
		progress_init("Calculating matrix :", Property::db_data.sequences);
	unsigned long int max_targetcount = 0;
	unsigned long int max_qgramcount = 0;
	for (unsigned long int row_id = thread_id; row_id < Property::db_data.sequences; row_id += total_thread) {
		std::vector<unsigned long int> temp_next;
		seqinfo_t * row_sequence = Property::db_data.get_seqinfo(row_id);
		bool using_reference = false;
		unsigned long int targetcount = 0;
		if (next_comparison.find(row_sequence->reference) == next_comparison.end()) {
			for (unsigned long col_id = row_id + 1; col_id < Property::db_data.sequences; col_id++) {
				unsigned int diff_length = abs(row_sequence->seqlen - Property::db_data.get_seqinfo(col_id)->seqlen);
				if (diff_length <= Property::resolution) {
					unsigned long qgramdiff = qgram_diff(Property::db_data.get_qgram_vector(row_id),
							Property::db_data.get_qgram_vector(col_id));
					if (qgramdiff <= Property::resolution) {
						targetampliconids[thread_id][targetcount] = col_id;
						targetcount++;
					} else if (qgramdiff <= Property::max_next) {
						temp_next.push_back(col_id);
					}
				} else if (diff_length <= Property::max_next) {
					temp_next.push_back(col_id);
				}
			}
		} else {
			using_reference = true;
			std::vector<unsigned long int> * references = &next_comparison[row_sequence->reference];

			for (unsigned int k = 0; k < references->size(); k++) {
				unsigned long int col_id = (*references)[k];
				if (col_id > row_id) {
					unsigned long qgramdiff = qgram_diff(Property::db_data.get_qgram_vector(row_id),
							Property::db_data.get_qgram_vector(col_id));
					if (qgramdiff <= Property::resolution) {
						targetampliconids[thread_id][targetcount] = col_id;
						targetcount++;
					} else if (qgramdiff <= Property::max_next) {
						temp_next.push_back(col_id);
					}
				}
			}
			int log_count = comparison_log[row_sequence->reference];
			if (log_count <= 1) {
				pthread_mutex_lock(&workmutex);
				next_comparison.erase(next_comparison.find(row_sequence->reference));
				pthread_mutex_unlock(&workmutex);
				temp_cleaned++;
#ifdef DEBUG
				fprintf(Property::dbdebug, "Clean next comparison %ld\n", row_sequence->reference);
#endif
			} else {
				pthread_mutex_lock(&workmutex);
				comparison_log[row_sequence->reference] = log_count - 1;
				pthread_mutex_unlock(&workmutex);
			}
#ifdef DEBUG
			fprintf(Property::dbdebug, "Log count for %ld is now %d\n", row_sequence->reference, log_count - 1);
#endif
		}
		total_scan += targetcount;

		if (targetcount > max_targetcount)
			max_targetcount = targetcount;

		scanner[thread_id].search_do(row_id, targetcount, targetampliconids[thread_id]);

		for (unsigned long j = 0; j < targetcount; j++) {
			if (scanner[thread_id].master_result[j].diff <= Property::resolution) {
				vector_put(thread_id, row_id, targetampliconids[thread_id][j]);
			} else if (scanner[thread_id].master_result[j].diff <= Property::max_next) {
				temp_next.push_back(targetampliconids[thread_id][j]);
			} else {
#ifdef DEBUG
				fprintf(Property::dbdebug, "%ld and %ld are far away by %ld\n", row_id, targetampliconids[j],
						scanner[thread_id].master_result[j].diff);
#endif
			}
		}
		if (!using_reference) {
			int log_count = 0;
			for (unsigned long j = 0; j < targetcount; j++) {
				if (scanner[thread_id].master_result[j].diff <= Property::resolution) {
					if (Property::db_data.get_seqinfo(targetampliconids[thread_id][j])->reference == targetampliconids[thread_id][j]
							&& targetampliconids[thread_id][j] > row_id + Property::threads * 2) {
						Property::db_data.get_seqinfo(targetampliconids[thread_id][j])->reference = row_id;
						log_count++;
					}
				}
			}
			if (log_count > 0) {
				total_data += temp_next.size();
				temp_written++;
				pthread_mutex_lock(&workmutex);
//				next_comparison[row_id] = temp_next;
//				comparison_log[row_id] = log_count;
				pthread_mutex_unlock(&workmutex);
#ifdef DEBUG
				fprintf(Property::dbdebug, "Create next comparison %ld Log count %d\n", row_id, log_count);
#endif
			}
		}
		if (thread_id == 0)
			progress_update(row_id);
	}
	if (thread_id == 0) {
		progress_done();
	}

	fprintf(Property::dbdebug, "Max qgram count\t\t: %ld\n", max_qgramcount);
	fprintf(Property::dbdebug, "Max target count\t\t: %ld\n", max_targetcount);
}

void Bigmatrix::form_clusters() {
	unsigned long int cluster_id = 1;
	for (int t = 0; t < Property::threads; t++) {
		while (matrix_x[t].size() > 0) {
			unsigned long int first = matrix_x[t].back();
			unsigned long int second = matrix_y[t].back();
			cluster_info * existing_first = result.find_member(first);
			cluster_info * existing_second = result.find_member(second);
			if (existing_first != NULL && existing_second == NULL) {
				result.add_member(existing_first, second);
#ifdef DEBUG
				fprintf(Property::dbdebug, "Add %ld to cluster %ld\n", second, existing_first->cluster_id);
#endif
			} else if (existing_first == NULL && existing_second != NULL) {
				result.add_member(existing_second, first);
#ifdef DEBUG
				fprintf(Property::dbdebug, "Add %ld to cluster %ld\n", first, existing_second->cluster_id);
#endif
			} else if (existing_first == NULL && existing_second == NULL) {
				cluster_info * added = result.new_cluster(cluster_id);
				result.add_member(added, first);
				result.add_member(added, second);
#ifdef DEBUG
				fprintf(Property::dbdebug, "Create cluster %ld for %ld and %ld\n", cluster_id, member1.id,
						member.id);
#endif
				cluster_id++;
			} else if (existing_first != NULL && existing_second != NULL) {
				if (existing_first->cluster_id != existing_second->cluster_id) {
#ifdef DEBUG
					fprintf(Property::dbdebug, "Merge cluster %ld with %ld\n", existing_first->cluster_id, existing_second->cluster_id);
#endif
					result.merge_cluster(existing_first, existing_second);
				}
			}
			matrix_x[t].pop_back();
			matrix_y[t].pop_back();
		}
	}
	std::vector<unsigned long int> singletons;
	for (unsigned long int i = 0; i < Property::db_data.sequences; i++) {
		if (result.find_member(i) == NULL) {
			singletons.push_back(i);
		}
	}
	for (unsigned long int i = 0; i < singletons.size(); i++) {
		cluster_info * added = result.new_cluster(cluster_id);
		result.add_member(added, singletons[i]);
		cluster_id++;
	}
}

void Bigmatrix::print_clusters() {
#ifdef DEBUG
	result.print(Property::outfile, true);
#else
	result.print(Property::outfile, false);
#endif
}

void Bigmatrix::print_debug() {
	fprintf(Property::dbdebug, "Total match\t\t: %ld\n", total_match);
	fprintf(Property::dbdebug, "Total estimate\t\t: %ld\n", total_qgram);
	fprintf(Property::dbdebug, "Total search\t\t: %ld\n", total_scan);
	fprintf(Property::dbdebug, "Temp written\t\t: %ld\n", temp_written);
	fprintf(Property::dbdebug, "Temp cleaned\t\t: %ld\n", temp_cleaned);
	fprintf(Property::dbdebug, "Total data\t\t: %ld\n", total_data);
	for (int t = 0; t < Property::threads; t++) {
		fprintf(Property::dbdebug, "Map [%d] size\t: %ld\n", t, matrix_x[t].size());
	}
}

void Bigmatrix::vector_put(int thread_id, unsigned long int row, unsigned long int col) {
	matrix_x[thread_id].push_back(row);
	matrix_y[thread_id].push_back(col);
#ifdef DEBUG
	fprintf(Property::dbdebug, "%ld and %ld are connected\n", row, col);
#endif
	total_match++;
}

