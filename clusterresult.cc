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
	clusters[cluster_id] = info;
	return &clusters[cluster_id];
}

struct compare_cluster {
	inline bool operator()(const std::pair<cluster_info, std::vector<member_info> > & struct1,
			const std::pair<cluster_info, std::vector<member_info> > & struct2) {
		return (struct1.second[0].sequence.header < struct2.second[0].sequence.header);
	}
};

struct compare_member {
	inline bool operator()(const member_info & struct1,
			const member_info & struct2) {
		return (struct1.sequence.header < struct2.sequence.header);
	}
};

/**
 * Need to print out consistent format (such as correct result will look exactly the same)
 * this will be an expensive method, turn off except for unit test
 */
void cluster_result::print(FILE * stream, bool sort) {
	long total = 0;
	long clust = 0;
	if (sort) {
		std::vector<std::pair<cluster_info, std::vector<member_info> > > vector_clusters;
		for (std::map<unsigned long int, cluster_info>::const_iterator cit = clusters.begin(); cit != clusters.end(); ++cit) {
			if (!cit->second.erased) {
				std::pair<cluster_info, std::vector<member_info> > pair;
				for (std::map<unsigned long int, member_info>::const_iterator it = cit->second.cluster_members.begin();
						it != cit->second.cluster_members.end(); ++it) {
					pair.second.push_back(it->second);
				}
				pair.first = cit->second;
				std::sort(pair.second.begin(), pair.second.end(), compare_member());
				vector_clusters.push_back(pair);
			}
		}
		std::sort(vector_clusters.begin(), vector_clusters.end(), compare_cluster());
		for (unsigned int i = 0; i < vector_clusters.size(); i++) {
			for (unsigned int j = 0; j < vector_clusters[i].second.size(); j++) {
				fprintf(stream, "\n%s", vector_clusters[i].second[j].sequence.header);
				total++;
			}
			fprintf(stream, "\n");
			clust++;
		}
	} else {
		for (std::map<unsigned long int, cluster_info>::const_iterator cit = clusters.begin(); cit != clusters.end(); ++cit) {
			if (!cit->second.erased) {
				for (std::map<unsigned long int, member_info>::const_iterator it = cit->second.cluster_members.begin();
						it != cit->second.cluster_members.end(); ++it) {
					fprintf(stream, "\n%s", it->second.sequence.header);
					total++;
				}
				fprintf(stream, "\n");
				clust++;
			}
		}
	}
	fprintf(stream, "\n\nTotal: %ld clusters of %ld sequences", clust, total);
}

void cluster_result::merge_cluster(cluster_info* cluster, cluster_info* merge) {
	for (std::map<unsigned long int, member_info>::const_iterator it = merge->cluster_members.begin(); it != merge->cluster_members.end();
			++it) {
		add_member(cluster, it->second);
	}
	merge->erased = true;
	cluster->max_generation += merge->max_generation * 2;
}

cluster_info * cluster_result::find_member(unsigned long int sequence_id) {
	if (member_stat.find(sequence_id) != member_stat.end()) {
		return &(clusters[member_stat[sequence_id]]);
	}
	return NULL;
}

void cluster_result::add_member(cluster_info* cluster, member_info member) {
	cluster->cluster_members[member.sequence.clusterid] = member;
	member_stat[member.sequence.clusterid] = cluster->cluster_id;
}
