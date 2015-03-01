/*
 * seqinfo.h
 *
 *  Created on: Feb 27, 2015
 *      Author: mimitantono
 */

#ifndef SEQINFO_H_
#define SEQINFO_H_

#include<pthread.h>

#define QGRAMLENGTH 4
#define QGRAMVECTORBITS (1<<(2*QGRAMLENGTH))
#define QGRAMVECTORBYTES (QGRAMVECTORBITS/8)

typedef unsigned char qgramvector_t[QGRAMVECTORBYTES];

class seqinfo_t {
private:
	pthread_mutex_t mutex;
	bool visited;
public:
	seqinfo_t();
	virtual ~seqinfo_t();
	char * header;
	char * seq;
	unsigned int headerlen;
	unsigned int seqlen;
	unsigned int abundance;
	unsigned int hdrhash;
	qgramvector_t qgram;
	void set_visited();
	bool is_visited();
};

#endif /* SEQINFO_H_ */
