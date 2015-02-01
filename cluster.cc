/*
 * Bigmatrix.cc
 *
 *  Created on: Dec 24, 2014
 *      Author: mimitantono
 */

#include "cluster.h"

Cluster::Cluster() {
	targetampliconids = new unsigned long*[Property::threads];
	scanner = new class scanner[Property::threads];
	row_stat = new unsigned long int[Property::threads];
	for (int i = 0; i < Property::threads; i++) {
		scanner[i].search_begin();
		targetampliconids[i] = new unsigned long[Property::db_data.sequences / 4];
		row_stat[i] = 0;
	}
	matrix_x = new std::vector<unsigned long int>[Property::threads];
	matrix_y = new std::vector<unsigned long int>[Property::threads];
	next_step = new std::vector<unsigned long int>[Property::threads];
	next_comparison = new std::vector<unsigned long int>[Property::threads];
	total_match = 0;
	total_qgram = 0;
	total_scan = 0;
	row_reference = 0;
	row_full = 0;
	total_data = 0;
	pthread_mutex_init(&workmutex, NULL);
}

Cluster::~Cluster() {
	if (scanner)
		delete[] scanner;
	if (matrix_x)
		delete[] matrix_x;
	if (matrix_y)
		delete[] matrix_y;
	if (next_step)
		delete[] next_step;
	if (next_comparison)
		delete[] next_comparison;
	if (row_stat)
		delete[] row_stat;
	for (int i = 0; i < Property::threads; i++) {
		if (targetampliconids[i])
			delete[] targetampliconids[i];
	}
	if (targetampliconids)
		delete[] targetampliconids;

	pthread_mutex_destroy(&workmutex);
}

void Cluster::run_thread(int thread_id, int total_thread) {
	if (thread_id == 0)
		progress_init("Calculating matrix :", Property::db_data.sequences);
	for (unsigned long int row_id = thread_id; row_id < Property::db_data.sequences; row_id += total_thread) {
		process_row(true, thread_id, row_id);
		for (unsigned long int next_row = 0; next_row < next_step[thread_id].size(); next_row++) {
			process_row(false, thread_id, next_step[thread_id][next_row]);
		}
		std::vector<unsigned long int>().swap(next_comparison[thread_id]);
		std::vector<unsigned long int>().swap(next_step[thread_id]);
	}
	if (thread_id == 0) {
		progress_done();
	}
}

/**
 * Excluding sequences with distance more than twice of resolution for the next connection must be every each row
 * If write_reference is true it means that this step will calculate distance with all columns and write a list of sequences
 * to be considered next step.
 * If write_reference is false, it means that this step will use the list of sequences that was calculated before and do not
 * need to write information for next_step.
 */
void Cluster::process_row(bool write_reference, int thread_id, unsigned long int row_id) {
	seqinfo_t * row_sequence = Property::db_data.get_seqinfo(row_id);
	if (row_sequence->visited) {
		return;
	} else {
		row_stat[thread_id]++;
		row_sequence->visited = true;
	}
	unsigned long int targetcount = 0;
	if (write_reference) {
		row_full++;
		for (unsigned long col_id = row_id + 1; col_id < Property::db_data.sequences; col_id++) {
			unsigned int diff_length = abs(row_sequence->seqlen - Property::db_data.get_seqinfo(col_id)->seqlen);
			if (diff_length <= Property::resolution) {
				unsigned long qgramdiff = qgram_diff(Property::db_data.get_qgram_vector(row_id),
						Property::db_data.get_qgram_vector(col_id));
				total_qgram++;
				if (qgramdiff <= Property::resolution) {
					targetampliconids[thread_id][targetcount] = col_id;
					targetcount++;
				} else if (write_reference && qgramdiff <= Property::max_next) {
					next_comparison[thread_id].push_back(col_id);
				}
			} else if (write_reference && diff_length <= Property::max_next) {
				next_comparison[thread_id].push_back(col_id);
			}
		}
	} else {
		row_reference++;
		for (unsigned int k = 0; k < next_comparison[thread_id].size(); k++) {
			unsigned long int col_id = next_comparison[thread_id][k];
			if (col_id > row_id) {
				unsigned long qgramdiff = qgram_diff(Property::db_data.get_qgram_vector(row_id),
						Property::db_data.get_qgram_vector(col_id));
				total_qgram++;
				if (qgramdiff <= Property::resolution) {
					targetampliconids[thread_id][targetcount] = col_id;
					targetcount++;
				}
			}
		}
	}
	total_scan += targetcount;

	scanner[thread_id].search_do(row_id, targetcount, targetampliconids[thread_id]);

	for (unsigned long j = 0; j < targetcount; j++) {
		if (scanner[thread_id].master_result[j].diff <= Property::resolution) {
			vector_put(thread_id, row_id, targetampliconids[thread_id][j]);
			if (write_reference) {
				next_step[thread_id].push_back(targetampliconids[thread_id][j]);
			}
		} else if (write_reference && scanner[thread_id].master_result[j].diff <= Property::max_next) {
			next_comparison[thread_id].push_back(targetampliconids[thread_id][j]);
		} else {
#ifdef DEBUG
			fprintf(Property::dbdebug, "%ld and %ld are far away by %ld\n", row_id, targetampliconids[thread_id][j],
					scanner[thread_id].master_result[j].diff);
#endif
		}
	}
	if (write_reference) {
//			std::unique(next_comparison[thread_id].begin(), next_comparison[thread_id].end());
		if (thread_id == 0)
			progress_update(row_id);
	}
}

void Cluster::form_clusters() {
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
				fprintf(Property::dbdebug, "Create cluster %ld for %ld and %ld\n", cluster_id,first,second);
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

void Cluster::print_clusters() {
#ifdef DEBUG
	result.print(Property::outfile, true);
#else
	result.print(Property::outfile, false);
#endif
}

void Cluster::print_debug() {
	fprintf(Property::dbdebug, "Total match\t\t: %ld\n", total_match);
	fprintf(Property::dbdebug, "Total estimate\t\t: %ld\n", total_qgram);
	fprintf(Property::dbdebug, "Total search\t\t: %ld\n", total_scan);
	fprintf(Property::dbdebug, "Full calculation\t: %ld\n", row_full);
	fprintf(Property::dbdebug, "Referenced calculation\t: %ld\n", row_reference);
	fprintf(Property::dbdebug, "Total data\t\t: %ld\n", total_data);
	for (int t = 0; t < Property::threads; t++) {
		fprintf(Property::dbdebug, "Map [%d] size\t: %ld\n", t, matrix_x[t].size());
	}
	for (int t = 0; t < Property::threads; t++) {
		fprintf(Property::dbdebug, "Row stat [%d]\t: %ld\n", t, row_stat[t]);
	}
}

void Cluster::vector_put(int thread_id, unsigned long int row, unsigned long int col) {
	matrix_x[thread_id].push_back(row);
	matrix_y[thread_id].push_back(col);
#ifdef DEBUG
	fprintf(Property::dbdebug, "%ld and %ld are connected\n", row, col);
#endif
	total_match++;
}

