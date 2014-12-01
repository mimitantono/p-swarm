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
	qgrams = 0;
	longest = 0;
	sequences = 0;
	nucleotides = 0;
	threadid = 0;
}

Db_data::~Db_data() {
	delete[] (qgrams);
	qgrams = NULL;
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

void Db_data::read_file(std::vector<Db_data*>& db, char* datap) {
	for (int i = 0; i < Property::threads; i++) {
		db.push_back(new Db_data());
	}
	/* allocate space */

	unsigned long dataalloc = MEMCHUNK;
	unsigned long datalen = 0;

	unsigned long sequences = 0;

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

	/* set up hash to check for unique headers */

	unsigned long hdrhashsize = 2 * sequences;

	seqinfo_t ** hdrhashtable = new seqinfo_t*[hdrhashsize];
	memset(hdrhashtable, 0, hdrhashsize * sizeof(seqinfo_t *));

	unsigned long duplicatedidentifiers = 0;

	/* create indices */
	int missingabundance = 0;

	int* longest_array = new int[Property::threads];
	unsigned long * sequences_array = new unsigned long[Property::threads];
	unsigned long * nucleotides_array = new unsigned long[Property::threads];
	long *lastabundance = new long[Property::threads];
	int *presorted = new int[Property::threads];

	int tail = sequences - (sequences / Property::threads) * Property::threads;
	for (int threadid = 0; threadid < Property::threads; threadid++) {
		if ((Property::threads - threadid) <= tail) {
//			db[threadid]->seqindex = new seqinfo_t[sequences / Property::threads];
		} else { //threadid is not tail, meaning that it will get more data
//			db[threadid]->seqindex = new seqinfo_t[sequences / Property::threads + 1];
		}
		longest_array[threadid] = 0;
		sequences_array[threadid] = 0;
		nucleotides_array[threadid] = 0;
		lastabundance[threadid] = -1;
		presorted[threadid] = 1;
	}

	regex_t db_regexp;
	regmatch_t pmatch[4];

	if (regcomp(&db_regexp, "(_)([0-9]+)$", REG_EXTENDED))
		fatal("Regular expression compilation failed");

	char * p = datap;
	for (unsigned long i = 0; i < sequences; i++) {
		int threadid = i % Property::threads;
		sequences_array[threadid]++;
		seqinfo_t seqinfo;
		db[threadid]->seqindex.push_back(seqinfo);
		seqinfo_t * seqindex_p = &(db[threadid]->seqindex.back());
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

		if (lastabundance[threadid] > 0 && seqindex_p->abundance > lastabundance[threadid])
			presorted[threadid] = 0;

		lastabundance[threadid] = seqindex_p->abundance;

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

		if (seqindex_p->seqlen > longest_array[threadid]) {
			longest_array[threadid] = seqindex_p->seqlen;
		}
		nucleotides_array[threadid] += seqindex_p->seqlen;

		p += seqindex_p->seqlen + 1;
	}

	if (missingabundance) {
		char * msg;
		asprintf(&msg, "Abundance annotation not found for %d sequences", missingabundance);
		fatal(msg);
	}

	delete[] (hdrhashtable);

	if (duplicatedidentifiers)
		exit(1);

	for (int threadid = 0; threadid < Property::threads; threadid++) {
//		if (!presorted[threadid]) {
//			fprintf(stderr, "Input file was not presorted, sorting now.\n");
//			qsort(db[threadid]->seqindex, db[threadid]->sequences, sizeof(seqinfo_t), compare_abundance);
//		}
		db[threadid]->nucleotides = nucleotides_array[threadid];
		db[threadid]->longest = longest_array[threadid];
		db[threadid]->sequences = sequences_array[threadid];
		db[threadid]->threadid = threadid;
		db[threadid]->qgrams_init();
//		db[threadid]->print_debug();
		db[threadid]->print_info();
	}

	delete[] longest_array;
	longest_array = NULL;
	delete[] presorted;
	presorted = NULL;
	delete[] sequences_array;
	sequences_array = NULL;
	delete[] nucleotides_array;
	nucleotides_array = NULL;
	delete[] lastabundance;
	lastabundance = NULL;
//		free(datap);
}

void Db_data::print_info() {
	fprintf(stderr, "Database info:     %lu nt", nucleotides);
	fprintf(stderr, " in %ld sequences,", sequences);
	fprintf(stderr, " longest %d nt\n", longest);
}

void Db_data::qgrams_init() {
	qgrams = new qgramvector_t[sequences];
	for (long i = 0; i < sequences; i++) {
		/* find qgrams */
		findqgrams((unsigned char*) seqindex[i].seq, seqindex[i].seqlen, qgrams[i]);
	}
}

seqinfo_t * Db_data::get_seqinfo(unsigned long seqno) {
	return &seqindex[seqno];
}

queryinfo_t Db_data::get_sequence_and_length(unsigned long seqno) {
	queryinfo_t query;
	query.seq = seqindex[seqno].seq;
	query.len = (long) (seqindex[seqno].seqlen);
	query.qno = seqno;
	return query;
}

void Db_data::put_seq(long seqno) {
	queryinfo_t query = get_sequence_and_length(seqno);
	for (int i = 0; i < query.len; i++)
		putchar(sym_nt[(int) (query.seq[i])]);
}

void Db_data::print_debug() {
	fprintf(Property::dbdebug, "\nThis is DB #%d containing %lu sequences", threadid, sequences);
	for (long i = 0; i < sequences; i++) {
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

