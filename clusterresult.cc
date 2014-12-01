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
		fprintf(Property::outfile, "\nCluster #%ld.%u: %lu members", partition_id, clusters[i].cluster_id,
				clusters[i].cluster_members.size());
		for (long j = 0; j < clusters[i].cluster_members.size(); j++) {
			fprintf(Property::outfile, "\n%ld. %s", j + 1, clusters[i].cluster_members[j].sequence.header);
			total++;
		}
	}
	fprintf(Property::outfile, "\n\n In total we have #%ld sequences", total);
}
