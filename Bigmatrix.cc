/*
 * Bigmatrix.cc
 *
 *  Created on: Dec 24, 2014
 *      Author: mimitantono
 */

#include "Bigmatrix.h"

Bigmatrix::Bigmatrix(Db_data * db) {
	this->db = db;
//	search_data = new struct search_data[Property::threads];
//	guestbook = new unsigned long int*[db->sequences];
	scanner = new class scanner[Property::threads];
	for (int i = 0; i < Property::threads; i++) {
		scanner[i].set_db(db);
		scanner[i].search_begin();
	}
	matrix = new std::vector<unsigned long int *>[Property::threads];
	total_match = 0;
	total_skip = 0;
	skip_by_guestbook = 0;
//	guestbook.set_empty_key(a_pair(-1,-1));
}

Bigmatrix::~Bigmatrix() {
//	delete guestbook;
//	delete matrix;
	delete[] scanner;
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
//		std::vector<unsigned long int> visited;

		unsigned long int qgram_count = 0;
		for (unsigned long j = 0; j < db->sequences - i - 1; j++) {
			bool include = true;
//			if (guestbook.find(a_pair(i, j + i + 1)) != guestbook.end()) {
//				include = false;
//				skip_by_guestbook++;
//			}
			if (include) {
				qgramamps[j] = j + i + 1;
				qgram_count++;
			}
		}
		qgram_work_diff(i, qgram_count, qgramamps, qgramdiffs, db);

		unsigned long int targetcount = 0;
		for (unsigned long j = 0; j < qgram_count; j++) {
			if (qgramdiffs[j] <= Property::resolution) {
				targetampliconids[targetcount] = qgramamps[j];
				targetcount++;
			}
//			if (qgramdiffs[j] == 0) {
//				visited.push_back(qgramamps[j]);
//			}
		}
		scanner[thread_id].search_do(i, targetcount, targetampliconids);

		for (unsigned long j = 0; j < targetcount; j++) {
			if (scanner[thread_id].master_result[j].diff <= Property::resolution) {
				vector_put(&matrix[thread_id], i, targetampliconids[j]);
//				visited.push_back(j);
			}
		}
//		if (visited.size() > 0) {
//			for (unsigned long int k = 0; k < visited.size() - 1; k++) {
//				for (unsigned long int l = k + 1; l < visited.size(); l++) {
//					if (visited[k] < visited[l] && visited[k] > i) {
//						pthread_mutex_lock(&workmutex);
//						guestbook[a_pair(visited[k], visited[l])] = true;
//						pthread_mutex_unlock(&workmutex);
//					} else if (visited[l] > i) {
//						pthread_mutex_lock(&workmutex);
//						guestbook[a_pair(visited[l], visited[k])] = true;
//						pthread_mutex_unlock(&workmutex);
//					}
//				}
//			}
//		}
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
		for (unsigned long int i = 0; i < matrix[t].size(); i++) {
			unsigned long int first = matrix[t][i][0];
			unsigned long int second = matrix[t][i][1];
			cluster_info * existing_first = result.find_member(first);
			cluster_info * existing_second = result.find_member(second);
			if (existing_first != NULL && existing_second == NULL) {
				member_info member;
				member.sequence = *db->get_seqinfo(second);
				member.sequence.clusterid = second;
				member.generation = 0;
				result.add_member(existing_first, member);
				fprintf(Property::dbdebug, "Add %s to cluster %ld\n", member.sequence.header, existing_first->cluster_id);
			} else if (existing_first == NULL && existing_second != NULL) {
				member_info member;
				member.sequence = *db->get_seqinfo(first);
				member.sequence.clusterid = first;
				member.generation = 0;
				result.add_member(existing_second, member);
				fprintf(Property::dbdebug, "Add %s to cluster %ld\n", member.sequence.header, existing_second->cluster_id);
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
				fprintf(Property::dbdebug, "Create cluster %ld for %s and %s\n", cluster_id, member1.sequence.header,
						member.sequence.header);
				cluster_id++;
			} else if (existing_first != NULL && existing_second != NULL) {
				if (existing_first->cluster_id != existing_second->cluster_id) {
					fprintf(Property::dbdebug, "Merge cluster %ld with %ld\n", existing_first->cluster_id, existing_second->cluster_id);
					result.merge_cluster(existing_first, existing_second);
				}
			}
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

//void Bigmatrix::form_clusters() {
//	progress_init("Forming clusters   :", db->sequences);
//	int cluster_id = 0;
//	bool * row_guestbook = new bool[db->sequences];
//	for (unsigned long int i = 0; i < db->sequences; i++) {
//		row_guestbook[i] = false;
//	}
//	for (unsigned long int i = 0; i < db->sequences; i++) {
//		if (row_guestbook[i] == false) {
//			if (has_match(i)) {
//				result.new_cluster(cluster_id);
//				member_info member;
//				member.sequence = *db->get_seqinfo(i);
//				member.generation = 0;
//				result.clusters[cluster_id].cluster_members.push_back(member);
//				crawl_row(&row_guestbook, cluster_id, i, 0);
//				cluster_id++;
//			}
//		}
//		progress_update(i);
//	}
//	for (unsigned long int i = 0; i < db->sequences; i++) {
//		if (row_guestbook[i] == false) {
//			result.new_cluster(cluster_id);
//			member_info member;
//			member.sequence = *db->get_seqinfo(i);
//			member.generation = 0;
//			result.clusters[cluster_id].cluster_members.push_back(member);
//			cluster_id++;
////			fprintf(stderr, "\n%s", member.sequence.header);
//		}
//	}
//	delete row_guestbook;
//	progress_done();
//}

bool Bigmatrix::has_match(unsigned long int row_id) {
	for (unsigned long int j = row_id + 1; j < db->sequences; j++) {
//		if (vector_contains(&matrix, row_id, j)) {
//			return true;
//		}
	}
	return false;
}

void Bigmatrix::crawl_row(bool ** row_guestbook, unsigned long int cluster_id, unsigned long int row_id, int generation) {
	(*row_guestbook)[row_id] = true;
	for (unsigned long int i = 0; i < db->sequences; i++) {
//		if (!(*row_guestbook)[i] && vector_contains(&matrix, row_id, i)) {
//			member_info member;
//			member.sequence = *db->get_seqinfo(i);
//			member.generation = generation;
//			result.clusters[cluster_id].cluster_members.push_back(member);
//			crawl_row(row_guestbook, cluster_id, i, generation + 1);
//		}
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
	fprintf(Property::dbdebug, "Guestbook size\t: %ld\n", guestbook.size());
	fprintf(Property::dbdebug, "Skipped by guestbook\t: %ld\n", skip_by_guestbook);
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

