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
}

Bigmatrix::~Bigmatrix() {
	delete[] scanner;
	for (int i = 0; i < Property::threads; i++) {
		matrix[i].clear();
	}
	delete[] matrix;
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

void Bigmatrix::calculate_partition(int thread_id, int total_thread) {
	unsigned long int * qgramamps = new unsigned long int[db->sequences];
	unsigned long int * qgramdiffs = new unsigned long int[db->sequences];
	unsigned long int * targetampliconids = new unsigned long[db->sequences];
	if (thread_id == 0)
		progress_init("Calculating matrix :", db->sequences);
	for (unsigned long int i = thread_id; i < db->sequences; i += total_thread) {
//		std::vector<int> visited;

		int qgram_count = 0;
		for (unsigned long j = i + 1; j < db->sequences; j++) {
			qgramamps[qgram_count] = j;
			qgram_count++;
		}
		qgram_work_diff(i, qgram_count, qgramamps, qgramdiffs, db);

		unsigned long int targetcount = 0;
		for (unsigned long j = 0; j < db->sequences - i - 1; j++) {
			if ((long) qgramdiffs[j] <= Property::resolution) {
				targetampliconids[targetcount] = qgramamps[j];
				targetcount++;
			}
		}
		scanner[thread_id].search_do(i, targetcount, targetampliconids);

		for (unsigned long j = 0; j < targetcount; j++) {
			if (scanner[thread_id].master_result[j].diff <= Property::resolution) {
				vector_put(&matrix[thread_id], i, targetampliconids[j]);
//				vector_put(&matrix[thread_id], targetampliconids[j], i);
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

int Bigmatrix::organize_sequences(int cluster_id, std::unordered_map<unsigned long int, unsigned long int>& map) {
	for (int t = 0; t < Property::threads; t++) {
		for (unsigned long int i = 0; i < matrix[t].size(); i++) {
			unsigned long int first = matrix[t][i][0];
			unsigned long int second = matrix[t][i][1];
			unsigned long int existing_first = map.find(first) == map.end() ? 0 : map[first];
			unsigned long int existing_second = map.find(second) == map.end() ? 0 : map[second];
			if (existing_first > 0 && existing_second > 0 && existing_first != existing_second) {
				for (std::unordered_map<unsigned long int, unsigned long int>::const_iterator it = map.begin(); it != map.end(); ++it) {
					if (it->second == existing_second) {
						map[it->first] = existing_first;
					}
				}
			} else if (existing_first == 0 && existing_second > 0) {
				map[first] = map[second];
			} else if (existing_first > 0 && existing_second == 0) {
				map[second] = map[first];
			} else {
				map[first] = cluster_id;
				map[second] = cluster_id;
				cluster_id++;
			}
		}
	}
	return cluster_id;
}

void Bigmatrix::push_members(int cluster_id, std::unordered_map<unsigned long int, unsigned long int>& map) {
	for (unsigned long int i = 1; i < cluster_id; i++) {
		cluster_info* new_cluster = result.new_cluster(cluster_id);
		for (std::unordered_map<unsigned long int, unsigned long int>::const_iterator it = map.begin(); it != map.end(); ++it) {
			if (it->second == i) {
				member_info member;
				member.sequence = *db->get_seqinfo(it->first);
				member.generation = 0;
				new_cluster->cluster_members.push_back(member);
			}
		}
		if (new_cluster->cluster_members.size() == 0) {
			new_cluster->erased = true;
		}
	}
}

void Bigmatrix::form_clusters() {
	std::unordered_map<unsigned long int, unsigned long int> map;
	int cluster_id = 1;
	cluster_id = organize_sequences(cluster_id, map);
	push_members(cluster_id, map);
//	std::vector<unsigned long int> singletons;
//	for (unsigned long int i = 0; i < db->sequences; i++) {
//		if (map.find(i) == map.end()) {
//			singletons.push_back(i);
//		}
//	}
//	for (unsigned long int i = 0; i < singletons.size(); i++) {
//		cluster_info * new_cluster = result.new_cluster(cluster_id);
//		member_info member;
//		member.sequence = *db->get_seqinfo(singletons[i]);
//		member.generation = 0;
//		new_cluster->cluster_members.push_back(member);
//		cluster_id++;
//	}
}

bool Bigmatrix::has_match(unsigned long int row_id) {
//	for (unsigned long int j = row_id + 1; j < db->sequences; j++) {
//		if (vector_contains(&matrix, row_id, j)) {
//			return true;
//		}
//	}
	return false;
}

void Bigmatrix::crawl_row(bool ** row_guestbook, unsigned long int cluster_id, unsigned long int row_id, int generation) {
	(*row_guestbook)[row_id] = true;
	std::vector<member_info> members;
	for (unsigned long int i = 0; i < db->sequences; i++) {
		if (!(*row_guestbook)[i] && vector_contains(matrix, row_id, i)) {
			member_info member;
			member.sequence = *db->get_seqinfo(i);
			member.generation = generation;
			member.sequence.clusterid = i;
			members.push_back(member);
		}
	}
	if (generation == 0 && members.size() > 0) {
		result.new_cluster(cluster_id);
		member_info member;
		member.sequence = *db->get_seqinfo(row_id);
		member.generation = 0;
		result.clusters[cluster_id].cluster_members[row_id] = member;
	}
	for (int i = 0; i < members.size(); i++) {
		result.clusters[cluster_id].cluster_members[i] = members[i];
		crawl_row(row_guestbook, cluster_id, members[i].sequence.clusterid, generation + 1);
	}
}

void Bigmatrix::print_clusters() {
	result.print(Property::outfile);
}

void Bigmatrix::print_matrix() {
	fprintf(Property::dbdebug, "Total match\t: %ld\nTotal skip\t: %ld\n", total_match, total_skip);
	for (int t = 0; t < Property::threads; t++) {
		fprintf(Property::dbdebug, "Map [%d] size\t: %ld\n", t, matrix[t].size());
	}

}

bool Bigmatrix::vector_contains(std::vector<unsigned long int *> * vector, unsigned long int row, unsigned long int col) {
	return true;
}

void Bigmatrix::vector_put(std::vector<unsigned long int *> * vector, unsigned long int row, unsigned long int col) {
	unsigned long int * item = new unsigned long int[2];
	item[0] = row;
	item[1] = col;
	vector->push_back(item);
	fprintf(Property::dbdebug, "%s and %s are connected\n", db->get_seqinfo(row)->header, db->get_seqinfo(col)->header);
	total_match++;
}

