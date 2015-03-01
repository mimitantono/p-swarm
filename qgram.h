#ifndef QGRAM_H_
#define QGRAM_H_

void findqgrams(unsigned char * seq, unsigned long seqlen, unsigned char * qgramvector);
void qgram_work_diff(unsigned long seed, unsigned long listlen, unsigned long * amplist, unsigned long * difflist);
unsigned long qgram_diff_by_id(unsigned long int a, unsigned long int b);
unsigned long qgram_diff(unsigned char * a, unsigned char * b);
void printqgrams(unsigned char * qgramvector);

#endif /* QGRAM_H_ */
