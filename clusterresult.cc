/*
 * clusterresult.cpp
 *
 *  Created on: Nov 7, 2014
 *      Author: mimitantono
 */

#include "clusterresult.h"

cluster_result::cluster_result() {
	partition_id = -1;
}

cluster_result::~cluster_result() {
}

cluster_info * cluster_result::new_cluster(long cluster_id) {
	cluster_info info;
	info.cluster_id = cluster_id;
	clusters.push_back(info);
	return &clusters.back();
}

void cluster_result::print() {
	long total = 0;
	fprintf(Property::outfile, "\n\nThere are %lu clusters produced for thread #%ld", clusters.size(), partition_id);
	for (long i = 0; i < clusters.size(); i++) {
		fprintf(Property::outfile, "\nCluster #%u max generation %d with %lu members", clusters[i].cluster_id,
				clusters[i].max_generation, clusters[i].cluster_members.size());
		for (long j = 0; j < clusters[i].cluster_members.size(); j++) {
			fprintf(Property::outfile, "\n%ld. %s", j + 1, clusters[i].cluster_members[j].sequence.header);
			total++;
		}
	}
	fprintf(Property::outfile, "\n\n In total we have #%ld sequences", total);
}

void cluster_result::merge_cluster(cluster_info* cluster, cluster_info* merge) {
	for (int i = 0; i < merge->cluster_members.size(); i++) {
		//need to determine max generation more precisely for the merged cluster otherwise it will be difficult
		//for second time merging (because we don't know actual max generation)
		cluster->cluster_members.push_back(merge->cluster_members[i]);
	}
	if (merge->max_generation > cluster->max_generation) {
		cluster->max_generation = merge->max_generation;
	}
}
