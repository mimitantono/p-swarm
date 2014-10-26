/*
 * main.cc
 *
 *  Created on: Oct 19, 2014
 *      Author: mimitantono
 */
#include "main.h"

using namespace std;

int main(int argc, char** argv) {
	init();
	cout << db_getqgramvector(1);
	destroy();
}

void init() {
	Matrix::score_matrix_init();
	CPU_Info::cpu_features_detect();
	CPU_Info::cpu_features_show();
	Property::init();
	Property::print();
	db_read("/Users/mimitantono/Downloads/bio/BioMarks.fas");
	algo_run();
}

void destroy() {
	Matrix::score_matrix_free();
	db_free();
	search_end();
}
