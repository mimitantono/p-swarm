/*
 * clusterresult.h
 *
 *  Created on: Nov 7, 2014
 *      Author: mimitantono
 */

#ifndef CLUSTERRESULT_H_
#define CLUSTERRESULT_H_

#include<vector>
#include<map>
#include<string>
#include<algorithm>
#include "db.h"

typedef struct member_info {
	seqinfo_t sequence;
	unsigned qgram_diff;
	unsigned generation;
	unsigned radius;
	unsigned char * qgrams;
} member_info;

typedef struct cluster_info {
	unsigned long int cluster_id;
	unsigned max_generation;
	std::map<unsigned long int, member_info> cluster_members;
	bool expired;
	bool erased;
} cluster_info;

class cluster_result {
public:
	cluster_result();
	virtual ~cluster_result();
	cluster_info * new_cluster(long cluster_id);
	long partition_id;
	void merge_cluster(cluster_info* cluster, cluster_info* merge);
	void print(FILE * stream);
	void add_member(cluster_info * cluster, member_info member);
	cluster_info * find_member(unsigned long int sequence_id);
private:
	std::map<unsigned long int, cluster_info> clusters;
	std::map<unsigned long int, unsigned long int> member_stat;
};

#endif /* CLUSTERRESULT_H_ */
