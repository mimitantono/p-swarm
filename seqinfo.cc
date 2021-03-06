/*
 * seqinfo.cc
 *
 *  Created on: Feb 27, 2015
 *      Author: mimitantono
 */

#include "seqinfo.h"
#include <stddef.h>

//pthread_mutex_t seqinfo_t::mutex;

seqinfo_t::seqinfo_t() {
	header = 0;
	seq = 0;
	headerlen = 0;
	hdrhash = 0;
	seqlen = 0;
	abundance = 0;
	visited = false;
//	pthread_mutex_init(&mutex, NULL);
}

seqinfo_t::~seqinfo_t() {
//	pthread_mutex_destroy(&mutex);
}

bool seqinfo_t::is_visited() {
	bool result;
//	pthread_mutex_lock(&mutex);
	result = visited;
//	pthread_mutex_unlock(&mutex);
	return result;
}

void seqinfo_t::set_visited() {
//	pthread_mutex_lock(&mutex);
	visited = true;
//	pthread_mutex_unlock(&mutex);
}

