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

typedef struct cluster_info {
	unsigned long int cluster_id;
	std::vector<unsigned long int> cluster_members;
} cluster_info;

class cluster_result {
public:
	cluster_result();
	virtual ~cluster_result();
	cluster_info * new_cluster(unsigned long int cluster_id);
	long partition_id;
	void merge_cluster(cluster_info* cluster, cluster_info* merge);
	void print(FILE * stream, bool sort);
	void add_member(cluster_info * cluster, unsigned long int id);
	cluster_info * find_member(unsigned long int sequence_id);
private:
	std::map<unsigned long int, cluster_info> clusters;
	unsigned long int * member_stat;
};

#endif /* CLUSTERRESULT_H_ */
