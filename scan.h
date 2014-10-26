/*
 * scan.h
 *
 *  Created on: Oct 26, 2014
 *      Author: mimitantono
 */

#ifndef SCAN_H_
#define SCAN_H_

#include <pthread.h>
#include "property.h"
#include "util.h"
#include "matrix.h"
#include "search.h"
#include "db.h"

void search_all(unsigned long query_no);
void search_do(unsigned long query_no,
               unsigned long listlength,
               unsigned long * targets,
               unsigned long * scores,
               unsigned long * diffs,
               unsigned long * alignlengths,
               long bits);
void search_begin();
void search_end();

#endif /* SCAN_H_ */
