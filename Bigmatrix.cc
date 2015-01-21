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
	total_skip = 0;
	total_qgram = 0;
	total_scan = 0;
}

Bigmatrix::~Bigmatrix() {
	delete[] scanner;
	scanner = NULL;
//	for (int i = 0; i < Property::threads; i++) {
//		for (unsigned int j = 0; j < matrix[i].size(); j++) {
//			delete[] matrix[i][j];
//		}
//		matrix[i].clear();
//	}
	delete[] matrix;
	matrix = NULL;
}

void Bigmatrix::init_partition(int thread_id, int total_thread) {
	if (thread_id == 0)
		progress_init("Initialize matrix  :", db->sequences);
	for (unsigned long int i = thread_id; i < db->sequences; i += total_thread) {
		if (thread_id == 0)
			progress_update(i);
	}
	if (thread_id == 0)
		progress_done();
}

void Bigmatrix::calculate_partition(int thread_id, int total_thread, pthread_mutex_t workmutex) {
	unsigned long int * qgramamps = new unsigned long int[db->sequences];
	unsigned long int * qgramdiffs = new unsigned long int[db->sequences];
	unsigned long int * targetampliconids = new unsigned long[db->sequences];
	if (thread_id == 0)
		progress_init("Calculating matrix :", db->sequences);
	for (unsigned long int i = thread_id; i < db->sequences; i += total_thread) {
		unsigned long int qgram_count = 0;
		for (unsigned long j = 0; j < db->sequences - i - 1; j++) {
			unsigned long int col = j + i + 1;
			if (abs(db->get_seqinfo(i)->seqlen - db->get_seqinfo(col)->seqlen) <= Property::resolution) {
				qgramamps[qgram_count] = col;
				qgram_count++;
			}
		}
		total_qgram += qgram_count;
		qgram_work_diff(i, qgram_count, qgramamps, qgramdiffs, db);

		unsigned long int targetcount = 0;
		for (unsigned long j = 0; j < qgram_count; j++) {
			if (qgramdiffs[j] <= Property::resolution) {
				targetampliconids[targetcount] = qgramamps[j];
				targetcount++;
//			} else {
//				if (Property::enable_debug)
//					fprintf(Property::dbdebug, "%ld and %ld has %ld different qgrams\n", i, qgramamps[j],
//							scanner[thread_id].master_result[j].diff);
			}
		}
		total_scan += targetcount;
		scanner[thread_id].search_do(i, targetcount, targetampliconids);

		for (unsigned long j = 0; j < targetcount; j++) {
			if (scanner[thread_id].master_result[j].diff <= Property::resolution) {
				vector_put(&matrix[thread_id], i, targetampliconids[j]);
//			} else {
//				if (Property::enable_debug)
//					fprintf(Property::dbdebug, "%ld and %ld are far away by %ld\n", i, targetampliconids[j],
//							scanner[thread_id].master_result[j].diff);
			}
		}
		if (thread_id == 0)
			progress_update(i);
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
				if (Property::enable_debug)
					fprintf(Property::dbdebug, "Add %ld to cluster %ld\n", second, existing_first->cluster_id);
			} else if (existing_first == NULL && existing_second != NULL) {
				member_info member;
				member.sequence = *db->get_seqinfo(first);
				member.sequence.clusterid = first;
				member.generation = 0;
				result.add_member(existing_second, member);
				if (Property::enable_debug)
					fprintf(Property::dbdebug, "Add %ld to cluster %ld\n", first, existing_second->cluster_id);
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
				if (Property::enable_debug)
					fprintf(Property::dbdebug, "Create cluster %ld for %ld and %ld\n", cluster_id, member1.sequence.clusterid,
							member.sequence.clusterid);
				cluster_id++;
			} else if (existing_first != NULL && existing_second != NULL) {
				if (existing_first->cluster_id != existing_second->cluster_id) {
					if (Property::enable_debug)
						fprintf(Property::dbdebug, "Merge cluster %ld with %ld\n", existing_first->cluster_id, existing_second->cluster_id);
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
	if (Property::enable_debug)
		result.print(Property::outfile, true);
	else
		result.print(Property::outfile, false);
}

void Bigmatrix::print_debug() {
	fprintf(Property::dbdebug, "Total match\t\t: %ld\n", total_match);
	fprintf(Property::dbdebug, "Total skip\t\t: %ld\n", total_skip);
	fprintf(Property::dbdebug, "Total estimate\t\t: %ld\n", total_qgram);
	fprintf(Property::dbdebug, "Total search\t\t: %ld\n", total_scan);
	for (int t = 0; t < Property::threads; t++) {
		fprintf(Property::dbdebug, "Map [%d] size\t: %ld\n", t, matrix[t].size());
	}
}

void Bigmatrix::vector_put(std::vector<unsigned long int *> * vector, unsigned long int row, unsigned long int col) {
	unsigned long int * item = new unsigned long int[2];
	item[0] = row;
	item[1] = col;
	vector->push_back(item);
//	if (Property::enable_debug)
//		fprintf(Property::dbdebug, "%ld and %ld are connected\n", row, col);
	total_match++;
}

