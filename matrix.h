/*
 * matrix.h
 *
 *  Created on: Oct 24, 2014
 *      Author: mimitantono
 */

#ifndef MATRIX_H_
#define MATRIX_H_

extern char sym_nt[];

class Matrix {
private:
	static void score_matrix_read();
	static void score_matrix_dump();
public:
	static long * score_matrix_63;
	static unsigned char * score_matrix_8;
	static unsigned short * score_matrix_16;
	static void score_matrix_free();
	static void score_matrix_init();
};

#endif /* MATRIX_H_ */
