/*
 * merger.h
 *
 *  Created on: Nov 16, 2014
 *      Author: mimitantono
 */

#ifndef MERGER_H_
#define MERGER_H_

#include "clusterresult.h"
#include "search.h"

class merger {
private:
	int count;
	cluster_result** cluster_results;
	cluster_result merge_result;
	bool merge_clusters(cluster_info cluster, cluster_info other);
public:
	merger(cluster_result** cluster_results, int count);
	virtual ~merger();
	void merge_groups();
};

#endif /* MERGER_H_ */
