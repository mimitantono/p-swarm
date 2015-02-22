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
	for (int thread_id = 0; thread_id < Property::threads; ++thread_id) {
		scanner[thread_id].search_begin();
	}

	row_stat_by_thread = new unsigned long int[Property::threads];
	for (int thread_id = 0; thread_id < Property::threads; ++thread_id) {
		row_stat_by_thread[thread_id] = 0;
	}

//	next_comparison = new std::vector<std::vector<unsigned long int> >[Property::threads];
	next_comparison = new std::vector<std::queue<unsigned long int> >[Property::threads];
	for (int thread_id = 0; thread_id < Property::threads; ++thread_id) {
		for (unsigned int i = 0; i < Property::depth; ++i) {
//			std::vector<unsigned long int> push;
			std::queue<unsigned long int> push;
			next_comparison[thread_id].push_back(push);
		}
	}

//	outdated = new bool*[Property::threads];
//	for (int thread_id = 0; thread_id < Property::threads; ++thread_id) {
//		outdated[thread_id] = new bool[Property::db_data.sequences];
//		memset(outdated[thread_id], false, Property::db_data.sequences);
//	}

	row_stat_by_iteration = new unsigned long int[Property::depth + 1];
	memset(row_stat_by_iteration, 0, (Property::depth + 1) * sizeof(unsigned long int));

	row_visited = new bool[Property::db_data.sequences];
	memset(row_visited, false, Property::db_data.sequences);

	next_step = new std::queue<unsigned long int>[Property::threads];
	next_step_level = new std::queue<unsigned int>[Property::threads];

	matches_found = 0;
	qgram_performed = 0;
	scan_performed = 0;
	row_reference = 0;
	row_full = 0;
	current_row_id = 1;
	current_cluster_id = 1;

	max_next = new unsigned long int[Property::depth + 1];
	for (unsigned int depth = 2; depth <= Property::depth; ++depth) {
		max_next[depth] = (depth + 1) * Property::resolution;
	}

	pthread_mutex_init(&result_mutex, NULL);
	pthread_mutex_init(&row_id_mutex, NULL);
}

Cluster::~Cluster() {
	if (scanner)
		delete[] scanner;
	if (next_step)
		delete[] next_step;
	if (next_step_level)
		delete[] next_step_level;
	if (next_comparison)
		delete[] next_comparison;
	if (row_stat_by_thread)
		delete[] row_stat_by_thread;
	if (row_stat_by_iteration)
		delete[] row_stat_by_iteration;
	if (targetampliconids)
		delete[] targetampliconids;
	if (row_visited)
		delete[] row_visited;
	if (max_next)
		delete[] max_next;

	pthread_mutex_destroy(&result_mutex);
	pthread_mutex_destroy(&row_id_mutex);
}

unsigned long int Cluster::get_next_row_id(int thread_id) {
	unsigned long int return_row = 0;
	pthread_mutex_lock(&row_id_mutex);
	while (current_row_id <= Property::db_data.sequences) {
		if (!row_visited[current_row_id - 1]) {
			return_row = current_row_id;
			row_visited[current_row_id - 1] = true;
			break;
		} else {
			++current_row_id;
		}
	}
	pthread_mutex_unlock(&row_id_mutex);
	return return_row;
}

void Cluster::reset_flags(int thread_id) {
	for (unsigned int i = 0; i < Property::depth; ++i) {
		while (!next_comparison[thread_id][i].empty())
			next_comparison[thread_id][i].pop();
	}
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
			reset_flags(thread_id);
		} else {
			process_row(false, false, thread_id, row_id, 0);
		}
		row_id = get_next_row_id(thread_id);
	}
	if (thread_id == 0) {
		progress_done();
	}
}

void Cluster::qgram_diff_full_row(unsigned long int row_id, int thread_id, bool write_reference) {
#ifdef DEBUG
	fprintf(Property::dbdebug, "Calculate row %ld full\n", row_id);
#endif
	++row_full;
	for (unsigned long col_id = row_id + 1; col_id < Property::db_data.sequences; ++col_id) {
		unsigned long qgramdiff = qgram_diff_by_id(row_id, col_id);
		if (qgramdiff <= Property::resolution) {
			targetampliconids[thread_id].push_back(col_id);
		} else {
			write_next_comparison(thread_id, col_id, qgramdiff);
		}
	}
	qgram_performed = qgram_performed + Property::db_data.sequences - row_id - 1;
}

void Cluster::prepare_alignment(unsigned long int col_id, unsigned long int row_id, int thread_id) {
	if (col_id > row_id) {
		if (qgram_diff_by_id(row_id, col_id) <= Property::resolution) {
			targetampliconids[thread_id].push_back(col_id);
		} else {
			next_comparison[thread_id][0].push(col_id);
#ifdef DEBUG
			fprintf(Property::dbdebug, "next_comparison[%d] push = %ld\n", 0, col_id);
#endif
		}
		++qgram_performed;
	} else {
		next_comparison[thread_id][0].push(col_id);
	}
}

void Cluster::walkthrough_row_by_reference(unsigned int iteration, int thread_id, unsigned long int row_id) {
#ifdef DEBUG
	fprintf(Property::dbdebug, "Calculate row %ld by reference iteration %d\n", row_id, iteration);
#endif
	++row_reference;
	for (unsigned int j = 0; j < iteration; j += iteration - 1) {
		unsigned long int size = next_comparison[thread_id][j].size();
		for (unsigned int k = 0; k < size; ++k) {
			unsigned long int col_id = next_comparison[thread_id][j].front();
#ifdef DEBUG
			fprintf(Property::dbdebug, "next_comparison[%d] pop = %ld\n", j, col_id);
#endif
			prepare_alignment(col_id, row_id, thread_id);
			next_comparison[thread_id][j].pop();
		}
	}
}

void Cluster::process_row(bool write_reference, bool use_reference, int thread_id, unsigned long int row_id, unsigned int iteration) {
	++row_stat_by_iteration[iteration];
	++row_stat_by_thread[thread_id];
	if (!use_reference) {
		qgram_diff_full_row(row_id, thread_id, write_reference);
	} else if (use_reference) {
		walkthrough_row_by_reference(iteration, thread_id, row_id);
	}

	scanner[thread_id].search_do(row_id, &targetampliconids[thread_id]);
	scan_performed += targetampliconids[thread_id].size();

	for (unsigned long j = 0; j < targetampliconids[thread_id].size(); ++j) {
		unsigned long int col_id = targetampliconids[thread_id][j];
		unsigned long int diff = scanner[thread_id].master_result[j].diff;
		if (diff <= Property::resolution) {
			add_match_to_cluster(thread_id, row_id, col_id);
			if (Property::enable_flag && !row_visited[col_id] && iteration < Property::depth) {
				row_visited[col_id] = true;
				next_step[thread_id].push(col_id);
				next_step_level[thread_id].push(iteration + 1);
			}
		} else if (write_reference || use_reference) {
			next_comparison[thread_id][0].push(col_id);
#ifdef DEBUG
			fprintf(Property::dbdebug, "next_comparison[%d] push = %ld\n", 0, col_id);
#endif
		}
#ifdef DEBUG
		fprintf(Property::dbdebug, "%ld and %ld are far away by edit distance %ld\n", row_id, col_id, diff);
#endif
	}
	std::vector<unsigned long int>().swap(targetampliconids[thread_id]);
	if (thread_id == 0)
		progress_update(row_full + row_reference);
}

void Cluster::find_and_add_singletons() {
	for (unsigned long int i = 0; i < Property::db_data.sequences; ++i) {
		if (result.find_member(i) == NULL) {
			cluster_info * added = result.new_cluster(current_cluster_id);
			result.add_member(added, i);
			++current_cluster_id;
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
	for (int t = 0; t < Property::threads; ++t) {
		fprintf(Property::dbdebug, "Row stat thread[%d]\t: %13ld\n", t, row_stat_by_thread[t]);
	}
	for (unsigned int i = 0; i < Property::depth; ++i) {
		fprintf(Property::dbdebug, "Row stat iteration[%d]\t: %13ld\n", i, row_stat_by_iteration[i]);
	}
}

void Cluster::add_match_to_cluster(int thread_id, unsigned long int row_id, unsigned long int col_id) {
#ifdef DEBUG
	fprintf(Property::dbdebug, "%ld and %ld are connected\n", row_id, col_id);
#endif
	pthread_mutex_lock(&result_mutex);
	cluster_info * row_cluster = result.find_member(row_id);
	cluster_info * col_cluster = result.find_member(col_id);
	if (row_cluster != NULL && col_cluster == NULL) {
		result.add_member(row_cluster, col_id);
#ifdef DEBUG
		fprintf(Property::dbdebug, "Add %ld to cluster %ld\n", col_id, row_cluster->cluster_id);
#endif
	} else if (row_cluster == NULL && col_cluster != NULL) {
		result.add_member(col_cluster, row_id);
#ifdef DEBUG
		fprintf(Property::dbdebug, "Add %ld to cluster %ld\n", row_id, col_cluster->cluster_id);
#endif
	} else if (row_cluster == NULL && col_cluster == NULL) {
		cluster_info * added = result.new_cluster(current_cluster_id);
		result.add_member(added, row_id);
		result.add_member(added, col_id);
#ifdef DEBUG
		fprintf(Property::dbdebug, "Create cluster %ld for %ld and %ld\n", current_cluster_id, row_id, col_id);
#endif
		++current_cluster_id;
	} else if (row_cluster != NULL && col_cluster != NULL) {
		if (row_cluster->cluster_id != col_cluster->cluster_id) {
#ifdef DEBUG
			fprintf(Property::dbdebug, "Merge cluster %ld with %ld\n", row_cluster->cluster_id, col_cluster->cluster_id);
#endif
			result.merge_cluster(row_cluster, col_cluster);
		}
	}
	pthread_mutex_unlock(&result_mutex);
	++matches_found;
}

inline void Cluster::write_next_comparison(int thread_id, unsigned long int col_id, unsigned int distance) {
	for (unsigned int i = 2; i <= Property::depth; ++i) {
		if (distance <= max_next[i]) {
#ifdef DEBUG
			fprintf(Property::dbdebug, "Next comparison distance %d [%d] = %ld\n", distance, i, col_id);
#endif
			next_comparison[thread_id][i - 1].push(col_id);
			break;
		}
	}
}

