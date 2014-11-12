#include "cpu_info.h"

long CPU_Info::mmx_present;
long CPU_Info::sse_present;
long CPU_Info::sse2_present;
long CPU_Info::sse3_present;
long CPU_Info::ssse3_present;
long CPU_Info::sse41_present;
long CPU_Info::sse42_present;
long CPU_Info::popcnt_present;
long CPU_Info::avx_present;
long CPU_Info::avx2_present;

void CPU_Info::cpu_features_detect() {
	unsigned int a, b, c, d;

	cpuid(0, 0, a, b, c, d);
	unsigned int maxlevel = a & 0xff;

	if (maxlevel >= 1) {
		cpuid(1, 0, a, b, c, d);
		mmx_present = (d >> 23) & 1;
		sse_present = (d >> 25) & 1;
		sse2_present = (d >> 26) & 1;
		sse3_present = (c >> 0) & 1;
		ssse3_present = (c >> 9) & 1;
		sse41_present = (c >> 19) & 1;
		sse42_present = (c >> 20) & 1;
		popcnt_present = (c >> 23) & 1;
		avx_present = (c >> 28) & 1;

		if (maxlevel >= 7) {
			cpuid(7, 0, a, b, c, d);
			avx2_present = (b >> 5) & 1;
		}
	}
}

void CPU_Info::cpu_features_show() {
	fprintf(stderr, "CPU features       : ");
	if (mmx_present)
		fprintf(stderr, " mmx");
	if (sse_present)
		fprintf(stderr, " sse");
	if (sse2_present)
		fprintf(stderr, " sse2");
	if (sse3_present)
		fprintf(stderr, " sse3");
	if (ssse3_present)
		fprintf(stderr, " ssse3");
	if (sse41_present)
		fprintf(stderr, " sse4.1");
	if (sse42_present)
		fprintf(stderr, " sse4.2");
	if (popcnt_present)
		fprintf(stderr, " popcnt");
	if (avx_present)
		fprintf(stderr, " avx");
	if (avx2_present)
		fprintf(stderr, " avx2");
	fprintf(stderr, "\n");
}
