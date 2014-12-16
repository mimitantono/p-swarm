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
	int result_count;
	cluster_result** cluster_results;
	bool merge_clusters(cluster_info * cluster, cluster_info * other);
public:
	merger(cluster_result** cluster_results, int count);
	virtual ~merger();
	void merge_groups();
	void final_merge();
	cluster_result merge_result;
};

#endif /* MERGER_H_ */
