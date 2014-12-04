/*
 * merger.cc
 *
 *  Created on: Nov 16, 2014
 *      Author: mimitantono
 */

#include "merger.h"

merger::merger(cluster_result** cluster_results, int count) {
	this->cluster_results = cluster_results;
	this->count = count;
	merge_result.partition_id = -1;
}

merger::~merger() {
}

void merger::merge_groups() {
	cluster_result temp;
	//asigning first group as seed
	for (int i = 0; i < (*cluster_results)[0].clusters.size(); i++) {
		temp.clusters.push_back((*cluster_results)[0].clusters[i]);
	}
	//comparing seed with each groups
	for (int i = 1; i < count; i++) {
		for (int j = 0; j < temp.clusters.size(); j++) {
			for (int k = 0; k < (*cluster_results)[i].clusters.size(); k++) {
//				fprintf(stderr, "Comparing cluster #0.%d with cluster #%d.%d\n", j, i, k);
				//comparing seed.j with i.k
				if (merge_clusters(temp.clusters[j], (*cluster_results)[i].clusters[k])) {
					(*cluster_results)[i].clusters.erase((*cluster_results)[i].clusters.begin() + k);
				}
			}
		}
		for (int k = 0; k < (*cluster_results)[i].clusters.size(); k++) {
			temp.clusters.push_back((*cluster_results)[i].clusters[k]);
		}
	}
	temp.print();
	while (!temp.clusters.empty()) {
		merge_result.clusters.push_back(temp.clusters.back());
		temp.clusters.pop_back();
		for (int i = 0; i < temp.clusters.size(); i++) {
			if (merge_clusters(merge_result.clusters.back(), temp.clusters[i])) {
				merge_result.merge_cluster(&merge_result.clusters.back(), &temp.clusters[i]);
				temp.clusters.erase(temp.clusters.begin() + i);
			}
		}
	}
//	merge_result.print();
}

bool merger::merge_clusters(cluster_info cluster, cluster_info other) {
	seqinfo_t query = cluster.cluster_members[0].sequence;
	seqinfo_t target = other.cluster_members[0].sequence;
	search_result result = searcher::search_single(&query, &target);
	//compare sum of max generations with distance of seeds
	//if less or equal, merge clusters
	if ((cluster.max_generation + other.max_generation) * Property::resolution <= result.diff) {
//		fprintf(stderr, "Diff of %s and %s %ld\n", query.header, target.header, result.diff);
		merge_result.merge_cluster(&cluster, &other);
		return true;
	}
	return false;

//otherwise proceed to search every member

//if any of cluster members are near with each other, merge clusters

//otherwise store the cluster as it was in the result set
}
