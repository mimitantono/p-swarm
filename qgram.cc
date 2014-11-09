#include "qgram.h"

qgramvector_t * qgrams = new qgramvector_t[sequences];

void printqgrams(unsigned char * qgramvector) {
	/* print qgramvector */
	fprintf(stderr, "qgram vector:\n");
	for (int i = 0; i < QGRAMVECTORBYTES; i++) {
		fprintf(stderr, "%02x", qgramvector[i]);
		if ((i % 32) == 31)
			fprintf(stderr, "\n");
	}
}

void findqgrams(unsigned char * seq, unsigned long seqlen, unsigned char * qgramvector) {
	/* set qgram bit vector by xoring occurrences of qgrams in sequence */

	memset(qgramvector, 0, QGRAMVECTORBYTES);

	unsigned long qgram = 0;
	unsigned long i = 0;

	while ((i < QGRAMLENGTH - 1) && (i < seqlen)) {
		qgram = (qgram << 2) | (seq[i] - 1);
		i++;
	}

	while (i < seqlen) {
		qgram = (qgram << 2) | (seq[i] - 1);
		qgramvector[(qgram >> 3) & (QGRAMVECTORBYTES - 1)] ^= (1 << (qgram & 7));
		i++;
	}
}

/* 
 Unable to get the Mac gcc compiler v 4.2.1 produce the real
 popcnt instruction. Therefore resorting to assembly code.
 */

#define popcnt_asm(x,y)                                         \
  __asm__ __volatile__ ("popcnt %1,%0" : "=r"(y) : "r"(x));

inline unsigned long popcount(unsigned long x) {
	unsigned long y;
	popcnt_asm(x, y);
	return y;
}

unsigned long popcount_128(__m128i x) {
	__m128i mask1 = _mm_set_epi8(0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55);

	__m128i mask2 = _mm_set_epi8(0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33);

	__m128i mask4 = _mm_set_epi8(0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f);

	__m128i zero = _mm_setzero_si128();

	/* add together 2 bits: 0+1, 2+3, 3+4, ... 126+127 */

	__m128i a = _mm_srli_epi64(x, 1);
	__m128i b = _mm_and_si128(x, mask1);
	__m128i c = _mm_and_si128(a, mask1);
	__m128i d = _mm_add_epi64(b, c);

	/* add together 4 bits: (0+1)+(2+3), ... (124+125)+(126+127) */

	__m128i e = _mm_srli_epi64(d, 2);
	__m128i f = _mm_and_si128(d, mask2);
	__m128i g = _mm_and_si128(e, mask2);
	__m128i h = _mm_add_epi64(f, g);

	/* add together 8 bits: (0..3)+(4..7), ... (120..123)+(124..127) */

	__m128i i = _mm_srli_epi64(h, 4);
	__m128i j = _mm_add_epi64(h, i);
	__m128i k = _mm_and_si128(j, mask4);

	/* add together 8 bytes: (0..63) and (64..127) */

	__m128i l = _mm_sad_epu8(k, zero);

	/* add together 64-bit values into final 128 bit value */

	__m128i m = _mm_srli_si128(l, 8);
	__m128i n = _mm_add_epi64(m, l);

	/* return low 64 bits: return value is always in range 0 to 128 */

	unsigned long o = (unsigned long) _mm_movepi64_pi64(n);

	return o;
}

unsigned long compareqgramvectors_128(unsigned char * a, unsigned char * b) {
	/* Count number of different bits */
	/* Uses SSE2 but not POPCNT instruction */
	/* input MUST be 16-byte aligned */

	__m128i * ap = (__m128i *) a;
	__m128i * bp = (__m128i *) b;
	unsigned long count = 0;

	while ((unsigned char*) ap < a + QGRAMVECTORBYTES)
		count += popcount_128(_mm_xor_si128(*ap++, *bp++));

	return count;
}

unsigned long compareqgramvectors_64(unsigned char * a, unsigned char * b) {
	/* Count number of different bits */
	/* Uses the POPCNT instruction, requires CPU with this feature */

	unsigned long *ap = (unsigned long*) a;
	unsigned long *bp = (unsigned long*) b;
	unsigned long count = 0;

	while ((unsigned char*) ap < a + QGRAMVECTORBYTES)
		count += popcount(*ap++ ^ *bp++);

	return count;
}

unsigned long compareqgramvectors(unsigned char * a, unsigned char * b) {
	if (CPU_Info::popcnt_present)
		return compareqgramvectors_64(a, b);
	else
		return compareqgramvectors_128(a, b);
}

inline unsigned long qgram_diff(unsigned long a, unsigned long b) {
	unsigned long diffqgrams = compareqgramvectors(db_getqgramvector(a), db_getqgramvector(b));
	unsigned long mindiff = (diffqgrams + 2 * QGRAMLENGTH - 1) / (2 * QGRAMLENGTH);
	return mindiff;
}

void qgram_work_diff(unsigned long seed, unsigned long listlen, unsigned long * amplist, unsigned long * difflist) {
	for (unsigned long i = 0; i < listlen; i++)
		difflist[i] = qgram_diff(seed, amplist[i]);
}

