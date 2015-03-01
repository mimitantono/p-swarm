/*
 * clusterdata.h
 *
 *  Created on: Mar 1, 2015
 *      Author: mimitantono
 */

#ifndef CLUSTERDATA_H_
#define CLUSTERDATA_H_

#include "scan.h"
#include<queue>
#include<map>
#include<vector>

class cluster_data {
public:
	cluster_data();
	virtual ~cluster_data();

	class scanner scanner;

	int thread_id;
	unsigned long int matches_found;
	unsigned long int qgram_performed;
	unsigned long int scan_performed;
	unsigned long int row_full;
	unsigned long int row_reference;
	unsigned long int row_stat;

	std::vector<unsigned long int> targetampliconids;
	std::queue<unsigned long int> next_step;
	std::queue<unsigned int> next_step_level;
	std::vector<unsigned long int> * next_comparison;
	bool* match_statistics;

	void write_next_comparison(unsigned long int col, unsigned int distance);
	void reset();
};

#endif /* CLUSTERDATA_H_ */
