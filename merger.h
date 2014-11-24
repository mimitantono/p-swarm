/*
 * merger.h
 *
 *  Created on: Nov 16, 2014
 *      Author: mimitantono
 */

#ifndef MERGER_H_
#define MERGER_H_

#include "clusterresult.h"

class merger {
private:
	int count;
	cluster_result** cluster_results;
public:
	merger(cluster_result** cluster_results, int count);
	virtual ~merger();
	void merge_results(cluster_result * merge_result);
};

#endif /* MERGER_H_ */
