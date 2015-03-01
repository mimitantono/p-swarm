/*
 * Bigmatrix.cc
 *
 *  Created on: Dec 24, 2014
 *      Author: mimitantono
 */

#include "cluster.h"
#include "qgram.h"
#include "scan.h"
#include "db.h"
#include "util.h"
#include <locale.h>
#include "seqinfo.h"
#include "property.h"

Cluster::Cluster() {
	targetampliconids = new std::vector<unsigned long int>[Property::threads];
	scanner = new class scanner[Property::threads];
	row_stat_by_thread = new unsigned long int[Property::threads];
	row_stat_by_iteration = new unsigned long int[Property::threads];
	next_comparison = new std::vector<unsigned long int> * [Property::threads];
	for (int i = 0; i < Property::threads; i++) {
		scanner[i].search_begin();
		row_stat_by_thread[i] = 0;
		row_stat_by_iteration[i] = 0;
		next_comparison[i] = new std::vector<unsigned long int>[Property::max_next + 1];
	}
	next_step = new std::queue<unsigned long int>[Property::threads];
	next_step_level = new std::queue<unsigned int>[Property::threads];
	match_statistics = new std::map<unsigned long int, bool>[Property::threads];
	matches_found = 0;
	qgram_performed = 0;
	scan_performed = 0;
	row_reference = 0;
	row_full = 0;
	current_row_id = 1;
	current_cluster_id = 1;
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
	if (match_statistics)
		delete[] match_statistics;
	if (row_stat_by_iteration)
		delete[] row_stat_by_iteration;
	if (row_stat_by_thread)
		delete[] row_stat_by_thread;
	if (targetampliconids)
		delete[] targetampliconids;

	pthread_mutex_destroy(&result_mutex);
	pthread_mutex_destroy(&row_id_mutex);
}

unsigned long int Cluster::get_next_row_id(int thread_id) {
	unsigned long int return_row = 0;
	pthread_mutex_lock(&row_id_mutex);
	if (current_row_id <= Property::db_data.sequences) {
		return_row = current_row_id;
		current_row_id++;
	}
	pthread_mutex_unlock(&row_id_mutex);
	return return_row;
}

void Cluster::run_thread(int thread_id, int total_thread) {
	if (thread_id == 0)
		progress_init("Calculating matrix :", Property::db_data.sequences);
	unsigned long int row_id = get_next_row_id(thread_id);
	while (row_id > 0) {
		row_id--; //0 means that loop should be finished
		if (Property::enable_flag) {
			process_row(true, false, thread_id, row_id, 1);
			while (next_step[thread_id].size() > 0) {
				process_row(false, true, thread_id, next_step[thread_id].front(), next_step_level[thread_id].front());
				next_step[thread_id].pop();
				next_step_level[thread_id].pop();
			}
			for (unsigned int i = 0; i < Property::max_next + 1; i++) {
				std::vector<unsigned long int>().swap(next_comparison[thread_id][i]);
			}
			std::map<unsigned long int, bool>().swap(match_statistics[thread_id]);
		} else {
			process_row(false, false, thread_id, row_id, 0);
		}
		row_id = get_next_row_id(thread_id);
	}
	if (thread_id == 0) {
		progress_done();
	}
}

void Cluster::process_row(bool write_reference, bool use_reference, int thread_id, unsigned long int row_id, unsigned int iteration) {
	seqinfo_t * row_sequence = Property::db_data.get_seqinfo(row_id);
	if (!use_reference && row_sequence->is_visited()) {
		return;
	} else {
		row_stat_by_thread[thread_id]++;
		row_sequence->set_visited();
	}
	fprintf(Property::dbdebug, "Calculate row %ld iteration %d\n", row_id, iteration);
	if (!use_reference) {
		row_full++;
		for (unsigned long col_id = row_id + 1; col_id < Property::db_data.sequences; col_id++) {
			seqinfo_t * col_sequence = Property::db_data.get_seqinfo(col_id);
			unsigned long qgramdiff = qgram_diff(row_sequence->qgram, col_sequence->qgram);
			if (qgramdiff <= Property::resolution) {
				targetampliconids[thread_id].push_back(col_id);
			}
			if (write_reference && qgramdiff <= Property::max_next) {
				next_comparison[thread_id][qgramdiff].push_back(col_id);
#ifdef DEBUG
				fprintf(Property::dbdebug, "%ld and %ld are far away by %ld\n", row_id, col_id, qgramdiff);
#endif
			}
		}
		qgram_performed = qgram_performed + Property::db_data.sequences - row_id - 1;
	} else if (use_reference) {
		row_reference++;
		unsigned int max_next = iteration * Property::resolution;
		unsigned int min_next = (iteration - 1) * Property::resolution;
		for (unsigned int j = 0; j <= max_next; j++) {
			for (unsigned int k = 0; k < next_comparison[thread_id][j].size(); k++) {
				unsigned long int col_id = next_comparison[thread_id][j][k];
				if (col_id > row_id && match_statistics[thread_id].find(col_id) == match_statistics[thread_id].end()) {
					seqinfo_t * col_sequence = Property::db_data.get_seqinfo(col_id);
					unsigned long qgramdiff = qgram_diff(row_sequence->qgram, col_sequence->qgram);
					if (qgramdiff <= Property::resolution) {
						targetampliconids[thread_id].push_back(col_id);
					}
				}
			}
		}
	}
	scan_performed += targetampliconids[thread_id].size();

	scanner[thread_id].search_do(row_id, &targetampliconids[thread_id]);

	for (unsigned long j = 0; j < targetampliconids[thread_id].size(); j++) {
		unsigned long int col_id = targetampliconids[thread_id][j];
		unsigned long int diff = scanner[thread_id].master_result[j].diff;
		if (diff <= Property::resolution) {
			add_match_to_cluster(thread_id, row_id, col_id);
			match_statistics[thread_id][col_id] = true;
			if (Property::enable_flag && !Property::db_data.get_seqinfo(col_id)->is_visited() && iteration < Property::depth) {
				Property::db_data.get_seqinfo(col_id)->set_visited();
				next_step[thread_id].push(col_id);
				next_step_level[thread_id].push(iteration + 1);
			}
		}
	}
	std::vector<unsigned long int>().swap(targetampliconids[thread_id]);
	if (thread_id == 0)
		progress_update(row_full + row_reference);
}

void Cluster::find_and_add_singletons() {
	for (unsigned long int i = 0; i < Property::db_data.sequences; i++) {
		if (result.find_member(i) == NULL) {
			cluster_info * added = result.new_cluster(current_cluster_id);
			result.add_member(added, i);
			current_cluster_id++;
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
	fprintf(Property::dbdebug, "Total match\t\t: %13ld\n", matches_found);
	fprintf(Property::dbdebug, "Total estimate\t\t: %13ld\n", qgram_performed);
	fprintf(Property::dbdebug, "Total search\t\t: %13ld\n", scan_performed);
	fprintf(Property::dbdebug, "Full calculation\t: %13ld\n", row_full);
	fprintf(Property::dbdebug, "Referenced calculation\t: %13ld\n", row_reference);
	for (int t = 0; t < Property::threads; t++) {
		fprintf(Property::dbdebug, "Row stat [%d]\t\t: %13ld\n", t, row_stat_by_thread[t]);
	}
}

void Cluster::add_match_to_cluster(int thread_id, unsigned long int first, unsigned long int second) {
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
		cluster_info * added = result.new_cluster(current_cluster_id);
		result.add_member(added, first);
		result.add_member(added, second);
#ifdef DEBUG
		fprintf(Property::dbdebug, "Create cluster %ld for %ld and %ld\n", current_cluster_id, first, second);
#endif
		current_cluster_id++;
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
	matches_found++;
}

void Cluster::write_next_comparison(int thread_id, unsigned long int col_id, unsigned int distance) {
	if (distance <= Property::max_next)
		next_comparison[thread_id][distance].push_back(col_id);
}

