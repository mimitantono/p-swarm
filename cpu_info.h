/*
 * cpu_info.h
 *
 *  Created on: Oct 25, 2014
 *      Author: mimitantono
 */

#ifndef CPU_INFO_H_
#define CPU_INFO_H_

#define cpuid(f1, f2, a, b, c, d)                                       \
  __asm__ __volatile__ ("cpuid"                                         \
                        : "=a" (a), "=b" (b), "=c" (c), "=d" (d)        \
                        : "a" (f1), "c" (f2));

class CPU_Info {
public:
	static long mmx_present;
	static long sse_present;
	static long sse2_present;
	static long sse3_present;
	static long ssse3_present;
	static long sse41_present;
	static long sse42_present;
	static long popcnt_present;
	static long avx_present;
	static long avx2_present;
	static void cpu_features_detect();
	static void cpu_features_show();
};

#endif /* CPU_INFO_H_ */
