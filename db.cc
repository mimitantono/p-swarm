#include "db.h"

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
	longest = 0;
	sequences = 0;
	nucleotides = 0;
	threadid = 0;
	datap = (char *) xmalloc(MEMCHUNK);
}

Db_data::~Db_data() {
	if (datap)
		free(datap);
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

bool Db_data::detect_duplicates() {
	/* set up hash to check for unique headers */

	unsigned long hdrhashsize = 0;

	hdrhashsize += sequences;
	hdrhashsize *= 2;

	seqinfo_t ** hdrhashtable = new seqinfo_t*[hdrhashsize];
	memset(hdrhashtable, 0, hdrhashsize * sizeof(seqinfo_t *));

	unsigned long duplicatedidentifiers = 0;
	/* check hash, fatal error if found, otherwize insert new */
	for (unsigned long j = 0; j < sequences; j++) {
		seqinfo_t *seqindex_p = &seqindex[j];
		unsigned long hdrhash = HASH((unsigned char*) seqindex_p->header, seqindex_p->headerlen);
		seqindex_p->hdrhash = hdrhash;
		unsigned long hashindex = hdrhash % hdrhashsize;

		seqinfo_t * found;

		while ((found = hdrhashtable[hashindex])) {
			if ((found->hdrhash == hdrhash) && (found->headerlen == seqindex_p->headerlen)
					&& (strncmp(found->header, seqindex_p->header, found->headerlen) == 0))
				break;
			hashindex = (hashindex + 1) % hdrhashsize;
			found = hdrhashtable[hashindex];
		}
		if (found) {
			duplicatedidentifiers++;
			fprintf(stderr, "Duplicated sequence identifier: %s\n", seqindex_p->header);
		}
		hdrhashtable[hashindex] = seqindex_p;
	}

	if (hdrhashtable)
		delete[] hdrhashtable;

	if (duplicatedidentifiers)
		return true;

	return false;
}

void Db_data::read_file() {
	/* allocate space */
	unsigned long dataalloc = MEMCHUNK;
	unsigned long datalen = 0;

	FILE * fp = NULL;
	if (Property::databasename.c_str()) {
		fp = fopen(Property::databasename.c_str(), "r");
		if (!fp)
			fatal("Error: Unable to open input data file (%s).", Property::databasename.c_str());
	} else
		fp = stdin;

	/* get file size */

	long filesize = 0;
	if (Property::databasename.c_str()) {
		if (fseek(fp, 0, SEEK_END))
			fatal("Error: Unable to seek in database file (%s)", Property::databasename.c_str());
		filesize = ftell(fp);
		rewind(fp);
	}

	char line[LINEALLOC];
	line[0] = 0;
	fgets(line, LINEALLOC, fp);

	long lineno = 1;

	progress_init("Reading database   :", filesize);
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

		while (line[0] && (line[0] != '>')) {
			char m;
			char c;
			char * p = line;
			while ((c = *p++))
				if ((m = map_nt[(int) c]) >= 0) {
					while (datalen >= dataalloc) {
						dataalloc += MEMCHUNK;
						datap = (char *) xrealloc(datap, dataalloc);
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
			dataalloc += MEMCHUNK;
			datap = (char *) xrealloc(datap, dataalloc);
		}

		*(datap + datalen) = 0;
		datalen++;

		sequences++;

		if (Property::databasename.c_str())
			progress_update(ftell(fp));
	}
	progress_done();

	fclose(fp);

	/* create indices */
	int missingabundance = 0;

	long lastabundance = -1;
	int presorted = 1;

	regex_t db_regexp;
	regmatch_t pmatch[4];

	if (regcomp(&db_regexp, "(_)([0-9]+)$", REG_EXTENDED))
		fatal("Regular expression compilation failed");

	char * p = datap;
	std::vector<seqinfo_t>(sequences).swap(seqindex);
	for (unsigned long i = 0; i < sequences; i++) {
		seqindex[i].header = p;
		seqindex[i].headerlen = strlen(seqindex[i].header);
		seqindex[i].visited = false;

		p += seqindex[i].headerlen + 1;

		/* get amplicon abundance */
		seqindex[i].abundance = 0;
		if (!regexec(&db_regexp, seqindex[i].header, 4, pmatch, 0)) {
			seqindex[i].abundance = atol(seqindex[i].header + pmatch[2].rm_so);
		}

		if (seqindex[i].abundance < 1)
			missingabundance++;

		if (lastabundance > 0 && seqindex[i].abundance > lastabundance)
			presorted = 0;

		lastabundance = seqindex[i].abundance;

		seqindex[i].seq = p;
		seqindex[i].seqlen = strlen(p);

		if (seqindex[i].seqlen > longest) {
			longest = seqindex[i].seqlen;
		}
		nucleotides += seqindex[i].seqlen;

		p += seqindex[i].seqlen + 1;
	}

	if (missingabundance) {
		char * msg;
		asprintf(&msg, "Abundance annotation not found for %d sequences", missingabundance);
		fatal(msg);
	}

	if (!presorted) {
		fatal("Input file was not presorted.\n");
	}
	qgrams_init();
	print_info();
	Property::longest = longest;

	if (detect_duplicates())
		fatal("There are duplicates in input file.");

	regfree(&db_regexp);
}

void Db_data::print_info() {
	fprintf(stderr, "Database info      : %lu nt", nucleotides);
	fprintf(stderr, " in %ld sequences,", sequences);
	fprintf(stderr, " longest %ld nt\n", longest);
}

void Db_data::qgrams_init() {
	for (unsigned long i = 0; i < sequences; i++) {
		/* find qgrams */
		findqgrams((unsigned char*) seqindex[i].seq, seqindex[i].seqlen, seqindex[i].qgram);
	}
}

seqinfo_t * Db_data::get_seqinfo(unsigned long seqno) {
	return &seqindex[seqno];
}

queryinfo_t Db_data::get_queryinfo(unsigned long seqno) {
	queryinfo_t query;
	query.seq = seqindex[seqno].seq;
	query.len = (long) (seqindex[seqno].seqlen);
	query.qno = seqno;
	return query;
}

void Db_data::put_seq(long seqno) {
	queryinfo_t query = get_queryinfo(seqno);
	for (int i = 0; i < query.len; i++)
		putchar(sym_nt[(int) (query.seq[i])]);
}

void Db_data::print_debug() {
	fprintf(Property::dbdebug, "\nThis is DB #%d containing %lu sequences", threadid, sequences);
	for (unsigned long i = 0; i < sequences; i++) {
		if (seqindex[i].header == NULL) {
			fatal("Sequence index header should not be null");
		} else if (seqindex[i].abundance == 0) {
			fatal("Sequence index abundance should not be null");
		} else if (seqindex[i].seq == NULL) {
			fatal("Sequence index content should not be null");
		} else {
			fprintf(Property::dbdebug, "\n %ld : %s [%d]\n%s", i, seqindex[i].header, seqindex[i].abundance, seqindex[i].seq);
		}
	}
	fclose(Property::dbdebug);
}

