/*
 * clusterresult.cpp
 *
 *  Created on: Nov 7, 2014
 *      Author: mimitantono
 */

#include "clusterresult.h"

cluster_result::cluster_result() {
	partition_id = 0;
}

cluster_result::~cluster_result() {
}

cluster_info * cluster_result::new_cluster() {
	clusters.push_back(new cluster_info());
	return clusters.back();
}

void cluster_result::print() {
	for (int i = 0; i < clusters.size(); i++) {
		fprintf(stderr, "\nCluster #%d.%u", partition_id, clusters[i]->cluster_id);
	}
}
