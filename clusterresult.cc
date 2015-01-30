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
	info.erased = false;
	clusters[cluster_id] = info;
	return &clusters[cluster_id];
}

struct compare_cluster {
	inline bool operator()(const std::pair<cluster_info, std::vector<unsigned long int> > & struct1,
			const std::pair<cluster_info, std::vector<unsigned long int> > & struct2) {
		return (Property::db_data.get_seqinfo(struct1.second[0])->header < Property::db_data.get_seqinfo(struct2.second[0])->header);
	}
};

struct compare_member {
	inline bool operator()(const unsigned long int id1, unsigned long int id2) {
		return (Property::db_data.get_seqinfo(id1)->header < Property::db_data.get_seqinfo(id2)->header);
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
		std::vector<std::pair<cluster_info, std::vector<unsigned long int> > > vector_clusters;
		for (std::map<unsigned long int, cluster_info>::const_iterator cit = clusters.begin(); cit != clusters.end(); ++cit) {
			if (!cit->second.erased) {
				std::pair<cluster_info, std::vector<unsigned long int> > pair;
				for (unsigned long i = 0; i < cit->second.cluster_members.size(); i++) {
					pair.second.push_back(cit->second.cluster_members[i]);
				}
				pair.first = cit->second;
				std::sort(pair.second.begin(), pair.second.end(), compare_member());
				vector_clusters.push_back(pair);
			}
		}
		std::sort(vector_clusters.begin(), vector_clusters.end(), compare_cluster());
		for (unsigned int i = 0; i < vector_clusters.size(); i++) {
			for (unsigned int j = 0; j < vector_clusters[i].second.size(); j++) {
				fprintf(stream, "\n%s", Property::db_data.get_seqinfo(vector_clusters[i].second[j])->header);
				total++;
			}
			fprintf(stream, "\n");
			clust++;
		}
	} else {
		for (std::map<unsigned long int, cluster_info>::const_iterator cit = clusters.begin(); cit != clusters.end(); ++cit) {
			if (!cit->second.erased) {
				for (unsigned long int i = 0; i < cit->second.cluster_members.size(); i++) {
					fprintf(stream, "\n%s", Property::db_data.get_seqinfo(cit->second.cluster_members[i])->header);
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
	for (unsigned long int i = 0; i < merge->cluster_members.size(); i++) {
		add_member(cluster, merge->cluster_members[i]);
	}
	merge->erased = true;
}

cluster_info * cluster_result::find_member(unsigned long int sequence_id) {
	if (member_stat.find(sequence_id) != member_stat.end()) {
		return &(clusters[member_stat[sequence_id]]);
	}
	return NULL;
}

void cluster_result::add_member(cluster_info* cluster, unsigned long int id) {
	cluster->cluster_members.push_back(id);
	member_stat[id] = cluster->cluster_id;
}
