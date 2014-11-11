#include "db.h"

#define MEMCHUNK 1048576
#define LINEALLOC LINE_MAX

char map_nt[256] = { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, 1, -1, 2, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 4, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		1, -1, 2, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 4, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };

char map_hex[256] = { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, -1, -1, -1, -1, -1, -1,
		-1, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };


Db_data::Db_data() {
	longestheader = 0;
	qgrams = 0;
	longest = 0;
	sequences = 0;
	headerchars = 0;
	seqindex = 0;
	nucleotides = 0;
	threadid = 0;
}

unsigned char * Db_data::get_qgram_vector(unsigned long seq_no) {
	return (unsigned char*) (qgrams + seq_no);
}

void Db_data::showseq(char * seq) {
	char * p = seq;
	while (char c = *p++) {
		putchar(sym_nt[(unsigned int) c]);
	}
}

int compare_abundance(const void * a, const void * b) {
	seqinfo_t * x = (seqinfo_t *) a;
	seqinfo_t * y = (seqinfo_t *) b;

	if (x->abundance > y->abundance)
		return -1;
	else if (x->abundance < y->abundance)
		return +1;
	else if (x < y)
		return -1;
	else if (x > y)
		return +1;
	else
		return 0;
}

bool Db_data::process_line(long line) {
	if ((line % (Property::threads * 2)) == threadid * 2)
		return true;
	else if ((line % (Property::threads * 2)) == threadid * 2 + 1)
		return true;
	else
		return false;
}

void Db_data::read_file(string filename) {
	/* allocate space */

	unsigned long dataalloc = MEMCHUNK;
	char * datap = (char *) xmalloc(dataalloc);
	unsigned long datalen = 0;

	longest = 0;
	longestheader = 0;
	sequences = 0;
	nucleotides = 0;
	headerchars = 0;

	FILE * fp = NULL;
	if (filename.c_str()) {
		fp = fopen(filename.c_str(), "r");
		if (!fp)
			fatal("Error: Unable to open input data file (%s).", filename.c_str());
	} else
		fp = stdin;

	/* get file size */

	long filesize = 0;
	if (filename.c_str()) {
		if (fseek(fp, 0, SEEK_END))
			fatal("Error: Unable to seek in database file (%s)", filename.c_str());
		filesize = ftell(fp);
		rewind(fp);
	}

	char line[LINEALLOC];
	line[0] = 0;
	fgets(line, LINEALLOC, fp);

	long lineno = 1;

	progress_init("Reading database: ", filesize);
	while (line[0]) {
		/* read header */
		/* the header ends at a space character, a newline or a nul character */

		if (line[0] != '>')
			fatal("Illegal header line in fasta file.");

		long headerlen = 0;
		if (char * stop = strpbrk(line + 1, " \n"))
			headerlen = stop - (line + 1);
		else
			headerlen = strlen(line + 1);

		headerchars += headerlen;

		if (headerlen > longestheader)
			longestheader = headerlen;

		/* store the header */

		while (datalen + headerlen + 1 > dataalloc) {
			dataalloc += MEMCHUNK;
			datap = (char *) xrealloc(datap, dataalloc);
		}
		memcpy(datap + datalen, line + 1, headerlen);
		*(datap + datalen + headerlen) = 0;
		datalen += headerlen + 1;

		/* get next line */

		line[0] = 0;
		fgets(line, LINEALLOC, fp);
		lineno++;

		/* read sequence */

		unsigned long seqbegin = datalen;

		while (line[0] && (line[0] != '>')) {
			char m;
			char c;
			char * p = line;
			while ((c = *p++))
				if ((m = map_nt[(int) c]) >= 0) {
					while (datalen >= dataalloc) {
						if (process_line(lineno)) {
							dataalloc += MEMCHUNK;
							datap = (char *) xrealloc(datap, dataalloc);
						}
					}

					*(datap + datalen) = m;
					datalen++;
				} else if (c != '\n') {
					char msg[100];
					snprintf(msg, 100, "Illegal character '%c' in sequence on line %ld", c, lineno);
					fatal(msg);
				}
			line[0] = 0;
			fgets(line, LINEALLOC, fp);
			lineno++;
		}

		while (datalen >= dataalloc) {
			if (process_line(lineno)) {
				dataalloc += MEMCHUNK;
				datap = (char *) xrealloc(datap, dataalloc);
			}
		}

		long length = datalen - seqbegin;

		nucleotides += length;

		if (length > longest)
			longest = length;

		*(datap + datalen) = 0;
		datalen++;

		sequences++;

		if (filename.c_str())
			progress_update(ftell(fp));
	}
	progress_done();

	fclose(fp);

	/* set up hash to check for unique headers */

	unsigned long hdrhashsize = 2 * sequences;

	seqinfo_t * * hdrhashtable = new seqinfo_t*[hdrhashsize];
	memset(hdrhashtable, 0, hdrhashsize * sizeof(seqinfo_t *));

	unsigned long duplicatedidentifiers = 0;

	/* create indices */

	seqindex = new seqinfo_t[sequences];
	seqinfo_t * seqindex_p = seqindex;

	regex_t db_regexp;
	regmatch_t pmatch[4];

	if (regcomp(&db_regexp, "(_)([0-9]+)$", REG_EXTENDED))
		fatal("Regular expression compilation failed");

	long lastabundance = LONG_MAX;

	int presorted = 1;
	int missingabundance = 0;

	char * p = datap;
	for (unsigned long i = 0; i < sequences; i++) {
		seqindex_p->header = p;
		seqindex_p->headerlen = strlen(seqindex_p->header);
		seqindex_p->headeridlen = seqindex_p->headerlen;

		p += seqindex_p->headerlen + 1;

		/* get amplicon abundance */
		seqindex_p->abundance = 0;
		if (!regexec(&db_regexp, seqindex_p->header, 4, pmatch, 0)) {
			seqindex_p->abundance = atol(seqindex_p->header + pmatch[2].rm_so);
			seqindex_p->abundance_start = pmatch[0].rm_so;
			seqindex_p->abundance_end = pmatch[0].rm_eo;
		} else {
			seqindex_p->abundance_start = 0;
			seqindex_p->abundance_end = 0;
		}

		if (seqindex_p->abundance < 1)
			missingabundance++;

		if (seqindex_p->abundance > lastabundance)
			presorted = 0;

		lastabundance = seqindex_p->abundance;

		/* check hash, fatal error if found, otherwize insert new */
		unsigned long hdrhash = HASH((unsigned char*) seqindex_p->header, seqindex_p->headeridlen);
		seqindex_p->hdrhash = hdrhash;
		unsigned long hashindex = hdrhash % hdrhashsize;

		seqinfo_t * found;

		while ((found = hdrhashtable[hashindex])) {
			if ((found->hdrhash == hdrhash) && (found->headeridlen == seqindex_p->headeridlen)
					&& (strncmp(found->header, seqindex_p->header, found->headeridlen) == 0))
				break;
			hashindex = (hashindex + 1) % hdrhashsize;
		}

		if (found) {
			duplicatedidentifiers++;
			fprintf(stderr, "Duplicated sequence identifier: %s\n", seqindex_p->header);
		}

		hdrhashtable[hashindex] = seqindex_p;

		seqindex_p->seq = p;
		seqindex_p->seqlen = strlen(p);
		p += seqindex_p->seqlen + 1;

		seqindex_p++;
	}

	if (missingabundance) {
		char * msg;
		asprintf(&msg, "Abundance annotation not found for %d sequences", missingabundance);
		fatal(msg);
	}

	if (!presorted) {
		qsort(seqindex, sequences, sizeof(seqinfo_t), compare_abundance);
	}

	delete (hdrhashtable);

	if (duplicatedidentifiers)
		exit(1);

	qgrams_init();
	free(datap);
}

void Db_data::print_info() {
	fprintf(stderr, "Database info:     %lu nt", nucleotides);
	fprintf(stderr, " in %ld sequences,", sequences);
	fprintf(stderr, " longest %d nt\n", longest);
}

void Db_data::qgrams_init() {
	qgrams = new qgramvector_t[sequences];

	seqinfo_t * seqindex_p = seqindex;
	for (long i = 0; i < sequences; i++) {
		/* find qgrams */
		findqgrams((unsigned char*) seqindex_p->seq, seqindex_p->seqlen, qgrams[i]);
		seqindex_p++;
	}
}

seqinfo_t * Db_data::get_seqinfo(unsigned long seqno) {
	return seqindex + seqno;
}

void Db_data::get_sequence_and_length(unsigned long seqno, char ** address, long * length) {
	*address = seqindex[seqno].seq;
	*length = (long) (seqindex[seqno].seqlen);
}

void Db_data::put_seq(long seqno) {
	char * seq;
	long len;
	get_sequence_and_length(seqno, &seq, &len);
	for (int i = 0; i < len; i++)
		putchar(sym_nt[(int) (seq[i])]);
}

void Db_data::print_debug() {
	fprintf(Property::debugfile, "\nThis is DB #%d containing %lu sequences", threadid, sequences);
	for (long i = 0; i < sequences; i++) {
		fprintf(Property::debugfile, "\n%ld: %s", i, seqindex[i].header);
	}
}

Db_data::~Db_data() {
	if (seqindex)
		delete (seqindex);
	if (qgrams)
		delete[] (qgrams);
}

