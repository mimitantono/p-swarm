/*
 * clusterdata.cc
 *
 *  Created on: Mar 1, 2015
 *      Author: mimitantono
 */

#include "clusterdata.h"
#include "property.h"
#include <string.h>

cluster_data::cluster_data() {
	thread_id = -1;
	next_comparison = new std::vector<unsigned long int>[Property::max_next + 1];
	matches_found = 0;
	qgram_performed = 0;
	scan_performed = 0;
	row_full = 0;
	row_reference = 0;
	row_stat = 0;
	match_statistics = new bool[Property::db_data.sequences];
	memset(match_statistics, 0, Property::db_data.sequences * sizeof(bool));
	scanner.search_begin();
}

cluster_data::~cluster_data() {
	delete[] next_comparison;
//	delete[] match_statistics;
}

void cluster_data::reset() {
	for (unsigned int i = 0; i < Property::max_next + 1; i++) {
		std::vector<unsigned long int>().swap(next_comparison[i]);
	}
//	std::map<unsigned long int, bool>().swap(match_statistics);
	memset(match_statistics, 0, Property::db_data.sequences * sizeof(bool));
}

void cluster_data::write_next_comparison(unsigned long int col_id, unsigned int distance) {
	if (distance <= Property::max_next)
		next_comparison[distance].push_back(col_id);
}

