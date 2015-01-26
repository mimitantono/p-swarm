/*
 * Bigmatrix.cc
 *
 *  Created on: Dec 24, 2014
 *      Author: mimitantono
 */

#include "Bigmatrix.h"

Bigmatrix::Bigmatrix(Db_data * db) {
	this->db = db;
	scanner = new class scanner[Property::threads];
	for (int i = 0; i < Property::threads; i++) {
		scanner[i].set_db(db);
		scanner[i].search_begin();
	}
	matrix = new std::vector<unsigned long int *>[Property::threads];
	total_match = 0;
	total_qgram = 0;
	total_scan = 0;
	temp_cleaned = 0;
	temp_written = 0;
}

Bigmatrix::~Bigmatrix() {
	delete[] scanner;
	scanner = NULL;
	delete[] matrix;
	matrix = NULL;
}

void Bigmatrix::calculate_partition(int thread_id, int total_thread, pthread_mutex_t workmutex) {
	unsigned long int * qgramamps = new unsigned long int[db->sequences];
	unsigned long int * qgramdiffs = new unsigned long int[db->sequences];
	unsigned long int * targetampliconids = new unsigned long[db->sequences];
	if (thread_id == 0)
		progress_init("Calculating matrix :", db->sequences);
	for (unsigned long int row_id = thread_id; row_id < db->sequences; row_id += total_thread) {
		std::vector<unsigned long int> temp_next;
		unsigned long int qgram_count = 0;
		seqinfo_t * row_sequence = db->get_seqinfo(row_id);
		bool using_reference = false;
		if (next_comparison.find(row_sequence->reference) == next_comparison.end()) {
			for (unsigned long col_id = row_id + 1; col_id < db->sequences; col_id++) {
				unsigned int diff_length = abs(row_sequence->seqlen - db->get_seqinfo(col_id)->seqlen);
				if (diff_length <= Property::resolution) {
					qgramamps[qgram_count] = col_id;
					qgram_count++;
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
					qgramamps[qgram_count] = col_id;
					qgram_count++;
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
		total_qgram += qgram_count;
		qgram_work_diff(row_id, qgram_count, qgramamps, qgramdiffs, db);

		unsigned long int targetcount = 0;
		for (unsigned long j = 0; j < qgram_count; j++) {
			if (qgramdiffs[j] <= Property::resolution) {
				targetampliconids[targetcount] = qgramamps[j];
				targetcount++;
			} else if (qgramdiffs[j] <= Property::max_next) {
				temp_next.push_back(qgramamps[j]);
			}
#ifdef DEBUG
			fprintf(Property::dbdebug, "%ld and %ld has %ld different qgrams\n", row_id, qgramamps[j], qgramdiffs[j]);
#endif
		}
		total_scan += targetcount;
		scanner[thread_id].search_do(row_id, targetcount, targetampliconids);

		for (unsigned long j = 0; j < targetcount; j++) {
			if (scanner[thread_id].master_result[j].diff <= Property::resolution) {
				vector_put(&matrix[thread_id], row_id, targetampliconids[j]);
			} else if (scanner[thread_id].master_result[j].diff <= Property::max_next) {
				temp_next.push_back(targetampliconids[j]);
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
					if (db->get_seqinfo(targetampliconids[j])->reference == targetampliconids[j]
							&& targetampliconids[j] > row_id + Property::threads * 2) {
						db->get_seqinfo(targetampliconids[j])->reference = row_id;
						log_count++;
					}
				}
			}
			if (log_count > 0) {
				total_data += temp_next.size();
				temp_written++;
				pthread_mutex_lock(&workmutex);
				next_comparison[row_id] = temp_next;
				comparison_log[row_id] = log_count;
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
	delete[] qgramamps;
	delete[] qgramdiffs;
	delete[] targetampliconids;
}

void Bigmatrix::form_clusters() {
	unsigned long int cluster_id = 1;
	for (int t = 0; t < Property::threads; t++) {
		while (matrix[t].size() > 0) {
			unsigned long int * pair = matrix[t].back();
			unsigned long int first = pair[0];
			unsigned long int second = pair[1];
			cluster_info * existing_first = result.find_member(first);
			cluster_info * existing_second = result.find_member(second);
			if (existing_first != NULL && existing_second == NULL) {
				member_info member;
				member.sequence = *db->get_seqinfo(second);
				member.sequence.clusterid = second;
				member.generation = 0;
				result.add_member(existing_first, member);
#ifdef DEBUG
				fprintf(Property::dbdebug, "Add %ld to cluster %ld\n", second, existing_first->cluster_id);
#endif
			} else if (existing_first == NULL && existing_second != NULL) {
				member_info member;
				member.sequence = *db->get_seqinfo(first);
				member.sequence.clusterid = first;
				member.generation = 0;
				result.add_member(existing_second, member);
#ifdef DEBUG
				fprintf(Property::dbdebug, "Add %ld to cluster %ld\n", first, existing_second->cluster_id);
#endif
			} else if (existing_first == NULL && existing_second == NULL) {
				cluster_info * added = result.new_cluster(cluster_id);
				member_info member;
				member.sequence = *db->get_seqinfo(first);
				member.sequence.clusterid = first;
				member.generation = 0;
				result.add_member(added, member);
				member_info member1;
				member1.sequence = *db->get_seqinfo(second);
				member1.sequence.clusterid = second;
				member1.generation = 0;
				result.add_member(added, member1);
#ifdef DEBUG
				fprintf(Property::dbdebug, "Create cluster %ld for %ld and %ld\n", cluster_id, member1.sequence.clusterid,
						member.sequence.clusterid);
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
			delete[] matrix[t].back();
			matrix[t].pop_back();
		}
	}
	std::vector<unsigned long int> singletons;
	for (unsigned long int i = 0; i < db->sequences; i++) {
		if (result.find_member(i) == NULL) {
			singletons.push_back(i);
		}
	}
	for (unsigned long int i = 0; i < singletons.size(); i++) {
		cluster_info * added = result.new_cluster(cluster_id);
		member_info member;
		member.sequence = *db->get_seqinfo(singletons[i]);
		member.sequence.clusterid = singletons[i];
		member.generation = 0;
		result.add_member(added, member);
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
		fprintf(Property::dbdebug, "Map [%d] size\t: %ld\n", t, matrix[t].size());
	}
}

void Bigmatrix::vector_put(std::vector<unsigned long int *> * vector, unsigned long int row, unsigned long int col) {
	unsigned long int * item = new unsigned long int[2];
	item[0] = row;
	item[1] = col;
	vector->push_back(item);
#ifdef DEBUG
	fprintf(Property::dbdebug, "%ld and %ld are connected\n", row, col);
#endif
	total_match++;
}

