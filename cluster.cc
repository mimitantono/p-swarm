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
	current_row_id = 1;
	current_cluster_id = 1;
	pthread_mutex_init(&result_mutex, NULL);
	pthread_mutex_init(&row_id_mutex, NULL);
}

Cluster::~Cluster() {
	pthread_mutex_destroy(&result_mutex);
	pthread_mutex_destroy(&row_id_mutex);
}

unsigned long int Cluster::get_next_row_id() {
	unsigned long int return_row = 0;
	pthread_mutex_lock(&row_id_mutex);
	if (current_row_id <= Property::db_data.sequences) {
		return_row = current_row_id++;
	}
	pthread_mutex_unlock(&row_id_mutex);
	return return_row;
}

void Cluster::run_thread(cluster_data *data, int total_thread) {
	if (data->thread_id == 0)
		progress_init("Calculating matrix :", Property::db_data.sequences);
	unsigned long int row_id = get_next_row_id();
	while (row_id > 0) {
		row_id--; //0 means that loop should be finished
		if (Property::enable_flag) {
			process_row(true, false, data, row_id, 1);
			while (data->next_step.size() > 0) {
				process_row(false, true, data, data->next_step.front(), data->next_step_level.front());
				data->next_step.pop();
				data->next_step_level.pop();
			}
			data->reset();
		} else {
			process_row(false, false, data, row_id, 1);
		}
		row_id = get_next_row_id();
	}
	if (data->thread_id == 0) {
		progress_done();
	}
}

void Cluster::process_row(bool write_reference, bool use_reference, cluster_data * data, unsigned long int row_id, unsigned int iteration) {
	seqinfo_t * row_sequence = Property::db_data.get_seqinfo(row_id);
	if (row_sequence->is_visited()) {
		return;
	}
	row_sequence->set_visited();
	++data->row_stat;
	++data->iteration_stat[iteration];
	if (!use_reference) {
		for (unsigned long col_id = row_id + 1; col_id < Property::db_data.sequences; ++col_id) {
			seqinfo_t * col_sequence = Property::db_data.get_seqinfo(col_id);
			unsigned long qgramdiff = qgram_diff(row_sequence->qgram, col_sequence->qgram);
			if (qgramdiff <= Property::resolution) {
				data->targetampliconids.push_back(col_id);
			} else if (write_reference && qgramdiff <= Property::max_next) {
				data->next_comparison[Property::max_next_map[qgramdiff]].push_back(col_id);
			}
		}
		data->qgram_performed += Property::db_data.sequences - row_id - 1;
	} else if (use_reference) {
		for (unsigned int j = 0; j <= iteration; ++j) {
			std::vector<unsigned long int> new_comparison;
			for (unsigned int k = 0; k < data->next_comparison[j].size(); ++k) {
				bool push = true;
				unsigned long int col_id = data->next_comparison[j][k];
				if (col_id > row_id) {
					seqinfo_t * col_sequence = Property::db_data.get_seqinfo(col_id);
					unsigned long qgramdiff = qgram_diff(row_sequence->qgram, col_sequence->qgram);
					if (qgramdiff <= Property::resolution) {
						data->targetampliconids.push_back(col_id);
						push = false;
					}
				}
				if (push) {
					new_comparison.push_back(col_id);
				}
			}
			new_comparison.swap(data->next_comparison[j]);
			std::vector<unsigned long int>().swap(new_comparison);
		}

	}
	data->scan_performed += data->targetampliconids.size();

	data->scanner.search_do(row_id, &data->targetampliconids);

	for (unsigned long j = 0; j < data->targetampliconids.size(); ++j) {
		unsigned long int col_id = data->targetampliconids[j];
		unsigned long int diff = data->scanner.master_result[j];
		if (diff <= Property::resolution) {
			add_match_to_cluster(data, row_id, col_id);
			if (Property::enable_flag && !Property::db_data.get_seqinfo(col_id)->is_visited() && iteration < Property::depth) {
//				Property::db_data.get_seqinfo(col_id)->set_visited();
				data->next_step.push(col_id);
				data->next_step_level.push(iteration + 1);
			}
		} else if (write_reference || use_reference) {
			data->next_comparison[iteration].push_back(col_id);
		}
	}
	std::vector<unsigned long int>().swap(data->targetampliconids);
	if (data->thread_id == 0)
		progress_update(current_row_id);
}

void Cluster::find_and_add_singletons() {
	for (unsigned long int i = 0; i < Property::db_data.sequences; ++i) {
		if (result.find_member(i) == NULL) {
			cluster_info * added = result.new_cluster(current_cluster_id++);
			result.add_member(added, i);
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

void Cluster::print_debug(cluster_data ** data) {
	unsigned long int matches_found = 0;
	unsigned long int qgram_performed = 0;
	unsigned long int scan_performed = 0;
	unsigned long int * iteration_stat = new unsigned long int[Property::depth];
	for (unsigned int j = 0; j < Property::depth; j++) {
		iteration_stat[j] = 0;
	}
	for (int t = 0; t < Property::threads; t++) {
		fprintf(Property::dbdebug, "Row stat [%d]\t\t: %13ld\n", t, ((*data)[t]).row_stat);
		matches_found += (*data)[t].matches_found;
		qgram_performed += (*data)[t].qgram_performed;
		scan_performed += (*data)[t].scan_performed;
		for (unsigned int j = 0; j < Property::depth; j++) {
			iteration_stat[j] += (*data)[t].iteration_stat[j + 1];
		}
	}
	fprintf(Property::dbdebug, "Total match\t\t: %13ld\n", matches_found);
	fprintf(Property::dbdebug, "Total estimate\t\t: %13ld\n", qgram_performed);
	fprintf(Property::dbdebug, "Total search\t\t: %13ld\n", scan_performed);
	fprintf(stderr, "Total estimate     : %ld\n", qgram_performed);
	fprintf(stderr, "Total search       : %ld\n", scan_performed);
	for (unsigned int j = 0; j < Property::depth; j++) {
		fprintf(Property::dbdebug, "Iteration[%d]\t: %13ld\n", j + 1, iteration_stat[j]);
	}
}

void Cluster::add_match_to_cluster(cluster_data * data, unsigned long int first, unsigned long int second) {
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
		cluster_info * added = result.new_cluster(current_cluster_id++);
		result.add_member(added, first);
		result.add_member(added, second);
#ifdef DEBUG
		fprintf(Property::dbdebug, "Create cluster %ld for %ld and %ld\n", current_cluster_id, first, second);
#endif
	} else if (existing_first != NULL && existing_second != NULL) {
		if (existing_first->cluster_id != existing_second->cluster_id) {
#ifdef DEBUG
			fprintf(Property::dbdebug, "Merge cluster %ld with %ld\n", existing_first->cluster_id, existing_second->cluster_id);
#endif
			if (existing_first->cluster_members.size() > existing_second->cluster_members.size()) {
				result.merge_cluster(existing_first, existing_second);
			} else {
				result.merge_cluster(existing_second, existing_first);
			}
		}
	}
#ifdef DEBUG
	fprintf(Property::dbdebug, "%ld and %ld are connected\n", first, second);
#endif
	++data->matches_found;
	pthread_mutex_unlock(&result_mutex);
}

