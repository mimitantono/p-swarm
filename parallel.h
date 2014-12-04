/*
 * parallel.h
 *
 *  Created on: Nov 2, 2014
 *      Author: mimitantono
 */

#ifndef PARALLEL_H_
#define PARALLEL_H_

#include <stdio.h>
#include <pthread.h>
#include "cluster.h"
#include "property.h"
#include "clusterresult.h"
#include "db.h"
#include "util.h"
#include "merger.h"

class Parallel {
public:
	Parallel();
	virtual ~Parallel();
	void run();
	static cluster_result * results;
};

#endif /* PARALLEL_H_ */

