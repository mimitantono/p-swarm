/*
 * search.cc
 *
 *  Created on: Nov 9, 2014
 *      Author: mimitantono
 */

#include "search.h"

searcher::searcher() {
	query = 0;
}

searcher::~searcher() {
}

void searcher::set_query(queryinfo_t *query) {
	this->query = query;
}
