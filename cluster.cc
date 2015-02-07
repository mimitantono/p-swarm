/*
 * Bigmatrix.cc
 *
 *  Created on: Dec 24, 2014
 *      Author: mimitantono
 */

#include "cluster.h"

Cluster::Cluster() {
	targetampliconids = new std::vector<unsigned long int>[Property::threads];
	scanner = new class scanner[Property::threads];
	row_stat = new unsigned long int[Property::threads];
	for (int i = 0; i < Property::threads; i++) {
		scanner[i].search_begin();
		row_stat[i] = 0;
	}
	next_step = new std::vector<unsigned long int>[Property::threads];
	next_comparison = new std::vector<unsigned long int>[Property::threads];
	total_match = 0;
	total_qgram = 0;
	total_scan = 0;
	row_reference = 0;
	row_full = 0;
	row_id_status = 1;
	cluster_id = 1;
	pthread_mutex_init(&result_mutex, NULL);
	pthread_mutex_init(&row_id_mutex, NULL);
}

Cluster::~Cluster() {
	if (scanner)
		delete[] scanner;
	if (next_step)
		delete[] next_step;
	if (next_comparison)
		delete[] next_comparison;
	if (row_stat)
		delete[] row_stat;
	if (targetampliconids)
		delete[] targetampliconids;

	pthread_mutex_destroy(&result_mutex);
	pthread_mutex_destroy(&row_id_mutex);
}

unsigned long int Cluster::row_id_dispenser() {
	unsigned long int return_row = 0;
	pthread_mutex_lock(&row_id_mutex);
	if (row_id_status <= Property::db_data.sequences) {
		return_row = row_id_status;
		row_id_status++;
	}
	pthread_mutex_unlock(&row_id_mutex);
	return return_row;
}

void Cluster::run_thread(int thread_id, int total_thread) {
	if (thread_id == 0)
		progress_init("Calculating matrix :", Property::db_data.sequences);
	unsigned long int row_id = row_id_dispenser();
	while (row_id > 0) {
		row_id--; //0 means that loop should be finished
		if (Property::enable_flag) {
			process_row(true, false, thread_id, row_id);
			for (unsigned long int next_row = 0; next_row < next_step[thread_id].size(); next_row++) {
				process_row(false, true, thread_id, next_step[thread_id][next_row]);
			}
			std::vector<unsigned long int>().swap(next_comparison[thread_id]);
			std::vector<unsigned long int>().swap(next_step[thread_id]);
		} else {
			process_row(false, false, thread_id, row_id);
		}
		row_id = row_id_dispenser();
	}
	if (thread_id == 0) {
		progress_done();
	}
}

void Cluster::process_row(bool write_reference, bool use_reference, int thread_id, unsigned long int row_id) {
	seqinfo_t * row_sequence = Property::db_data.get_seqinfo(row_id);
	if (!use_reference && row_sequence->visited) {
		return;
	} else {
		row_stat[thread_id]++;
		row_sequence->visited = true;
	}
	if (!use_reference) {
		row_full++;
		for (unsigned long col_id = row_id + 1; col_id < Property::db_data.sequences; col_id++) {
			unsigned long qgramdiff = qgram_diff(row_sequence->qgram, Property::db_data.get_seqinfo(col_id)->qgram);
			if (qgramdiff <= Property::resolution) {
				targetampliconids[thread_id].push_back(col_id);
			} else if (write_reference && qgramdiff <= Property::max_next) {
				next_comparison[thread_id].push_back(col_id);
			}
		}
		total_qgram = total_qgram + Property::db_data.sequences - row_id - 1;
	} else {
		row_reference++;
		for (unsigned int k = 0; k < next_comparison[thread_id].size(); k++) {
			unsigned long int col_id = next_comparison[thread_id][k];
			if (col_id > row_id) {
				unsigned long qgramdiff = qgram_diff(row_sequence->qgram, Property::db_data.get_seqinfo(col_id)->qgram);
				if (qgramdiff <= Property::resolution) {
					targetampliconids[thread_id].push_back(col_id);
				}
			}
		}
		total_qgram = total_qgram + next_comparison[thread_id].size();
	}
	total_scan += targetampliconids[thread_id].size();

	scanner[thread_id].search_do(row_id, &targetampliconids[thread_id]);

	for (unsigned long j = 0; j < targetampliconids[thread_id].size(); j++) {
		unsigned long int col_id = targetampliconids[thread_id][j];
		if (scanner[thread_id].master_result[j].diff <= Property::resolution) {
			vector_put(thread_id, row_id, col_id);
			if (write_reference && !Property::db_data.get_seqinfo(col_id)->visited) {
				Property::db_data.get_seqinfo(col_id)->visited = true;
				next_step[thread_id].push_back(col_id);
			}
		} else if (write_reference && scanner[thread_id].master_result[j].diff <= Property::max_next) {
			next_comparison[thread_id].push_back(col_id);
		} else {
#ifdef DEBUG
			fprintf(Property::dbdebug, "%ld and %ld are far away by %ld\n", row_id, targetampliconids[thread_id][j],
					scanner[thread_id].master_result[j].diff);
#endif
		}
	}
	std::vector<unsigned long int>().swap(targetampliconids[thread_id]);
	if (thread_id == 0)
		progress_update(row_full + row_reference);
}

void Cluster::form_clusters() {
	for (unsigned long int i = 0; i < Property::db_data.sequences; i++) {
		if (result.find_member(i) == NULL) {
			cluster_info * added = result.new_cluster(cluster_id);
			result.add_member(added, i);
			cluster_id++;
		}
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
	for (int t = 0; t < Property::threads; t++) {
		fprintf(Property::dbdebug, "Row stat [%d]\t: %ld\n", t, row_stat[t]);
	}
}

void Cluster::vector_put(int thread_id, unsigned long int first, unsigned long int second) {
	pthread_mutex_lock(&result_mutex);
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
		fprintf(Property::dbdebug, "Create cluster %ld for %ld and %ld\n", cluster_id, first, second);
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
#ifdef DEBUG
	fprintf(Property::dbdebug, "%ld and %ld are connected\n", first, second);
#endif
	pthread_mutex_unlock(&result_mutex);
	total_match++;
}

