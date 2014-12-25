/*
 * Bigmatrix.cc
 *
 *  Created on: Dec 24, 2014
 *      Author: mimitantono
 */

#include "Bigmatrix.h"

Bigmatrix::Bigmatrix(Db_data * db) {
	this->db = db;
	matrix = new int*[db->sequences];
	unsigned long int dirbuffersize = 8 * Property::longest * ((Property::longest + 3) / 4) * 4;
	search_data = new struct search_data[Property::threads];
	for (int i = 0; i < Property::threads; i++) {
		search_data[i].qtable = new BYTE*[Property::longest];
		search_data[i].qtable_w = new WORD*[Property::longest];
		search_data[i].dprofile = new BYTE[4 * 16 * 32];
		search_data[i].dprofile_w = new WORD[4 * 2 * 8 * 32];
		search_data[i].hearray = new BYTE[Property::longest * 32];
		search_data[i].dir_array = new unsigned long[dirbuffersize];
		search_data[i].target_count = 1;
		memset(search_data[i].hearray, 0, Property::longest * 32);
		memset(search_data[i].dir_array, 0, dirbuffersize);
	}
}

Bigmatrix::~Bigmatrix() {
	for (unsigned long int i = 0; i < db->sequences; i++) {
		if (matrix[i])
			delete matrix[i];
	}
	if (matrix)
		delete[] matrix;
}

void Bigmatrix::init_partition(int thread_id, int total_thread) {
	if (thread_id == 0)
		progress_init("Initialize matrix  :", db->sequences);
	for (unsigned long int i = thread_id; i < db->sequences; i += total_thread) {
		matrix[i] = new int[db->sequences];
		if (thread_id == 0)
			progress_update(i);
	}
	if (thread_id == 0)
		progress_done();
}

void Bigmatrix::calculate_partition(int thread_id, int total_thread) {
	if (thread_id == 0)
		progress_init("Calculating matrix :", db->sequences);
	for (unsigned long int i = thread_id; i < db->sequences; i += total_thread) {
		std::vector<int> guestbook;
		for (unsigned long int j = i + 1; j < db->sequences; j++) {
			if (matrix[i][j] != 3) {
				if (qgram_diff(db->get_qgram_vector(i), db->get_qgram_vector(j)) <= Property::resolution) {
					seqinfo_t* _query = db->get_seqinfo(i);
					seqinfo_t* _target = db->get_seqinfo(j);
					search_result _result = searcher::search_single(_query, _target, &search_data[thread_id]);
					if (_result.diff <= Property::resolution) {
						matrix[i][j] = 1;
						matrix[j][i] = 1;
						guestbook.push_back(j);
					} else {
						matrix[i][j] = 0;
					}
				} else {
					matrix[i][j] = 0;
				}
			}
		}
		if (guestbook.size() > 0) {
			for (int k = 0; k < guestbook.size() - 1; k++) {
				for (int l = k + 1; l < guestbook.size(); l++) {
					matrix[guestbook[k]][guestbook[l]] = 3;
					matrix[guestbook[l]][guestbook[k]] = 3;
				}
			}
		}
		if (thread_id == 0)
			progress_update(i);
	}
	if (thread_id == 0)
		progress_done();
}

void Bigmatrix::form_clusters() {
	progress_init("Forming clusters   :", db->sequences);
	int cluster_id = 0;
	bool * row_guestbook = new bool[db->sequences];
	for (unsigned long int i = 0; i < db->sequences; i++) {
		row_guestbook[i] = false;
	}
	for (unsigned long int i = 0; i < db->sequences; i++) {
		if (row_guestbook[i] == false) {
			if (has_match(i)) {
				result.new_cluster(cluster_id);
				member_info member;
				member.sequence = *db->get_seqinfo(i);
				member.generation = 0;
				result.clusters[cluster_id].cluster_members.push_back(member);
				crawl_row(&row_guestbook, cluster_id, i, 0);
				cluster_id++;
			}
		}
		progress_update(i);
	}
	for (unsigned long int i = 0; i < db->sequences; i++) {
		if (row_guestbook[i] == false) {
			result.new_cluster(cluster_id);
			member_info member;
			member.sequence = *db->get_seqinfo(i);
			member.generation = 0;
			result.clusters[cluster_id].cluster_members.push_back(member);
			cluster_id++;
//			fprintf(stderr, "\n%s", member.sequence.header);
		}
	}
	delete row_guestbook;
	progress_done();
}

bool Bigmatrix::has_match(int row_id) {
	for (unsigned long int j = row_id + 1; j < db->sequences; j++) {
		if (matrix[row_id][j] == 1) {
			return true;
		}
	}
	return false;
}

void Bigmatrix::crawl_row(bool ** row_guestbook, unsigned long int cluster_id, unsigned long int row_id, int generation) {
	(*row_guestbook)[row_id] = true;
	for (unsigned long int i = 0; i < db->sequences; i++) {
		if (!(*row_guestbook)[i] && matrix[row_id][i] == 1) {
			member_info member;
			member.sequence = *db->get_seqinfo(i);
			member.generation = generation;
			result.clusters[cluster_id].cluster_members.push_back(member);
			crawl_row(row_guestbook, cluster_id, i, generation + 1);
		}
	}
}

void Bigmatrix::print_clusters() {
	result.print(Property::outfile);
}

void Bigmatrix::print_matrix() {
	long int total_match = 0;
	long int total_skipped = 0;
	for (unsigned long int i = 0; i < db->sequences; i++) {
		for (unsigned long int j = 0; j < db->sequences; j++) {
			if (matrix[i][j] == 1) {
				fprintf(Property::dbdebug, "%s and %s are connected\n", db->get_seqinfo(i)->header, db->get_seqinfo(j)->header);
				total_match++;
			} else if (matrix[i][j] == 3) {
				fprintf(Property::dbdebug, "%s and %s are skipped\n", db->get_seqinfo(i)->header, db->get_seqinfo(j)->header);
				total_skipped++;
			}
		}
		fprintf(Property::dbdebug, "\n");
	}
	fprintf(Property::dbdebug, "Total match\t: %ld\nTotal skip\t: %ld", total_match, total_skipped);
}
