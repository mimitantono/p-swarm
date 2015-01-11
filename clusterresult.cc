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
	info.max_generation = 1;
	info.erased = false;
	info.expired = true;
	clusters.push_back(info);
	return &clusters.back();
}

struct compare_cluster {
	inline bool operator()(const cluster_info& struct1, const cluster_info& struct2) {
		return (struct1.cluster_members[0].sequence.header < struct2.cluster_members[0].sequence.header);
	}
};

struct compare_member {
	inline bool operator()(const member_info & struct1, const member_info & struct2) {
		return (struct1.sequence.header < struct2.sequence.header);
	}
};

struct sort_erased {
	inline bool operator()(const cluster_info& struct1, const cluster_info& struct2) {
		return (struct1.erased < struct2.erased);
	}
};

/**
 * Need to print out consistent format (such as correct result will look exactly the same)
 * this will be an expensive method, turn off except for unit test
 */
void cluster_result::print(FILE * stream) {
	std::sort(clusters.begin(), clusters.end(), sort_erased());
	while (clusters.back().erased) {
		clusters.pop_back();
	}
	for (unsigned int i = 0; i < clusters.size(); i++) {
		std::sort(clusters[i].cluster_members.begin(), clusters[i].cluster_members.end(), compare_member());
	}
	std::sort(clusters.begin(), clusters.end(), compare_cluster());
	long total = 0;
	for (unsigned int i = 0; i < clusters.size(); i++) {
		for (unsigned int j = 0; j < clusters[i].cluster_members.size(); j++) {
			fprintf(stream, "\n%s", clusters[i].cluster_members[j].sequence.header);
			total++;
		}
		fprintf(stream, "\n");
	}
	fprintf(stream, "\n\n In total we have %ld clusters of %ld sequences", clusters.size(), total);
}

void cluster_result::merge_cluster(cluster_info* cluster, cluster_info* merge) {
	for (unsigned int i = 0; i < merge->cluster_members.size(); i++) {
		//need to determine max generation more precisely for the merged cluster otherwise it will be difficult
		//for second time merging (because we don't know actual max generation)
		cluster->cluster_members.push_back(merge->cluster_members[i]);
	}
	cluster->max_generation += merge->max_generation * 2;
}
