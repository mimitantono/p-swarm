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
}

merger::~merger() {
	zap(cluster_results);
}

void merger::merge_results(cluster_result * merge_result) {
	for (int i = 1; i < count; i++) {
		for (int j = 0; j < cluster_results[0]->clusters.size(); j++) {
			//Comparing cluster 0,j with cluster i,k
			for (int k = 0; k < cluster_results[i]->clusters.size(); k++) {

				//compare sum of max generations with distance of seeds

				//if less or equal, merge clusters

				//otherwise proceed to search every member

				//if any of cluster members are near with each other, merge clusters

				//otherwise store the cluster as it was in the result set
			}
		}
	}
}
