/*
 * merger.cc
 *
 *  Created on: Nov 16, 2014
 *      Author: mimitantono
 */

#include "merger.h"

merger::merger(cluster_result** cluster_results, int count) {
	this->cluster_results = cluster_results;
	this->result_count = count;
	merge_result.partition_id = -1;
}

merger::~merger() {
}

void merger::merge_groups() {
	//asigning first group as seed
	for (int h = 0; h < result_count; h++) {
		for (int i = h + 1; i < result_count; i++) {
			int recursive = 0;
			bool repeat = true;
			while (repeat) {
				repeat = false;
				for (long j = 0; j < (*cluster_results)[h].clusters.size(); j++) {
					if (recursive == 0 || (*cluster_results)[i].clusters[j].expired) {
						long merged = 0;
						for (long k = 0; k < (*cluster_results)[i].clusters.size(); k++) {
							if (merge_clusters(&(*cluster_results)[h].clusters[j], &(*cluster_results)[i].clusters[k])) {
								merged++;
							}
						}
						if (merged > 0) {
							(*cluster_results)[i].clusters[j].expired = true;
							repeat = true;
						} else {
							(*cluster_results)[i].clusters[j].expired = false;
						}
					}
				}
				recursive++;
			}
			fprintf(stderr, "\nRecursively compared first stage #%d with #%d for %d times", h, i, recursive);
		}
	}
//comparing seed with each groups
	for (int i = 0; i < result_count; i++) {
		for (long k = 0; k < (*cluster_results)[i].clusters.size(); k++) {
			if (!(*cluster_results)[i].clusters[k].erased) {
				merge_result.clusters.push_back((*cluster_results)[i].clusters[k]);
			}
		}
	}
}

void merger::final_merge() {
	bool repeat = true;
	int recursive = 0;
	while (repeat) {
		repeat = false;
		for (long i = 0; i < merge_result.clusters.size() - 1; i++) {
			if (recursive == 0 || merge_result.clusters[i].expired) {
				long merged = 0;
				for (long j = i + 1; j < merge_result.clusters.size(); j++) {
					if (merge_clusters(&merge_result.clusters[i], &merge_result.clusters[j])) {
						merged++;
					}
				}
				if (merged > 0) {
					merge_result.clusters[i].expired = true;
//					repeat = true;
				} else {
					merge_result.clusters[i].expired = false;
				}
			}
		}
		recursive++;
	}
	fprintf(stderr, "\nRecursively compared second stage for %d times", recursive);
}

bool merger::merge_clusters(cluster_info *cluster, cluster_info* other) {
	if (other->erased || cluster->erased)
		return false;
	int min_member = 100;
	if (cluster->cluster_members.size() > min_member && other->cluster_members.size() > min_member) {
		int max_seed_distance = (cluster->max_generation + other->max_generation) * Property::resolution;
		if (qgram_diff(cluster->cluster_members[0].qgrams, other->cluster_members[0].qgrams) > max_seed_distance) {
			return false;
		} else {
			seqinfo_t query = cluster->cluster_members[0].sequence;
			seqinfo_t target = other->cluster_members[0].sequence;
			search_result result = searcher::search_single(&query, &target, searcher::search_data);
			if (result.diff > (cluster->max_generation + other->max_generation) * Property::resolution) {
				return false;
			}
		}
	}
	for (int i = cluster->cluster_members.size() - 1; i >= 0; i--) {
		for (int j = other->cluster_members.size() - 1; j >= 0; j--) {
			if (qgram_diff(cluster->cluster_members[i].qgrams, other->cluster_members[j].qgrams) <= Property::resolution) {
				seqinfo_t _query = cluster->cluster_members[i].sequence;
				seqinfo_t _target = other->cluster_members[j].sequence;
				search_result _result = searcher::search_single(&_query, &_target, searcher::search_data);
//				fprintf(stderr, "Diff of %s and %s %ld\n", _query.header, _target.header, _result.diff);
				if (_result.diff <= Property::resolution) {
					merge_result.merge_cluster(cluster, other);
					other->erased = true;
					return true;
				}
			}
		}
	}
	return false;
}
