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
	for (int i = 0; i < (*cluster_results)[0].clusters.size(); i++) {
		merge_result.clusters.push_back((*cluster_results)[0].clusters[i]);
	}
	//comparing seed with each groups
	for (int i = 1; i < result_count; i++) {
		for (long j = 0; j < merge_result.clusters.size(); j++) {
			for (long k = 0; k < (*cluster_results)[i].clusters.size(); k++) {
				if (merge_clusters(&merge_result.clusters[j], &(*cluster_results)[i].clusters[k])) {
					merge_result.merge_cluster(&merge_result.clusters[j], &(*cluster_results)[i].clusters[k]);
					merge_result.clusters[j].expired = true;
				}
			}
		}
		for (long k = 0; k < (*cluster_results)[i].clusters.size(); k++) {
			if (!(*cluster_results)[i].clusters[k].erased) {
				merge_result.clusters.push_back((*cluster_results)[i].clusters[k]);
			}
		}
	}
}

void merger::final_merge() {
	bool repeat = true;
	while (repeat) {
		repeat = false;
		for (long i = 0; i < merge_result.clusters.size() - 1; i++) {
			long merged = 0;
			for (int j = i + 1; j < merge_result.clusters.size(); j++) {
				if (merge_result.clusters[i].expired) {
					if (merge_clusters(&merge_result.clusters[i], &merge_result.clusters[j])) {
						merge_result.merge_cluster(&merge_result.clusters[i], &merge_result.clusters[j]);
						merged++;
					}
				}
			}
			if (merged > 0) {
				merge_result.clusters[i].expired = true;
				repeat = true;
			} else {
				merge_result.clusters[i].expired = false;
			}
		}
	}
}

bool merger::merge_clusters(cluster_info *cluster, cluster_info* other) {
	if (other->erased)
		return false;
//	fprintf(stderr, "Test merge cluster\n");
	seqinfo_t query = cluster->cluster_members[0].sequence;
	seqinfo_t target = other->cluster_members[0].sequence;
	search_result result = searcher::search_single(&query, &target);
	if (result.diff > (cluster->max_generation + other->max_generation) * Property::resolution) {
		return false;
	}
	for (int i = 0; i < cluster->cluster_members.size(); i++) {
		for (int j = 0; j < other->cluster_members.size(); j++) {
			seqinfo_t _query = cluster->cluster_members[i].sequence;
			seqinfo_t _target = other->cluster_members[j].sequence;
			search_result _result = searcher::search_single(&_query, &_target);
//			fprintf(stderr, "Diff of %s and %s %ld\n", _query.header, _target.header, _result.diff);
			if (_result.diff <= Property::resolution) {
				other->erased = true;
				return true;
			}
		}
	}
	return false;
}
