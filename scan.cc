#include "scan.h"

scanner::scanner() {
	sd = new struct search_data;
	master_alignlengths = 0;
	master_bits = 0;
	master_diffs = 0;
	master_length = 0;
	master_next = 0;
	master_scores = 0;
	master_targets = 0;
	remainingchunks = 0;
	dirbufferbytes = 0;
	db = 0;
}

scanner::~scanner() {
	delete (sd);
	delete[] (sd->qtable);
	delete[] (sd->qtable_w);
	delete[] (sd->dprofile);
	delete[] (sd->dprofile_w);
	delete[] (sd->hearray);
	delete[] (sd->dir_array);
	sd = NULL;
}

void scanner::search_alloc() {
	dirbufferbytes = 8 * db->longest * ((db->longest + 3) / 4) * 4;
	sd->qtable = new BYTE*[db->longest];
	sd->qtable_w = new WORD*[db->longest];
	sd->dprofile = new BYTE[4 * 16 * 32];
	sd->dprofile_w = new WORD[4 * 2 * 8 * 32];
	sd->hearray = new BYTE[db->longest * 32];
	sd->dir_array = new unsigned long[dirbufferbytes];

	memset(sd->hearray, 0, db->longest * 32);
	memset(sd->dir_array, 0, dirbufferbytes);
}

void scanner::search_init() {
	for (long i = 0; i < db->query.len; i++) {
		sd->qtable[i] = sd->dprofile + 64 * db->query.seq[i];
		sd->qtable_w[i] = sd->dprofile_w + 32 * db->query.seq[i];
	}
}

void scanner::search_chunk(long bits) {
	if (sd->target_count == 0)
		return;

	if (bits == 16)
		searcher.search16(sd->qtable_w, Property::penalty_gapopen, Property::penalty_gapextend, (WORD*) Matrix::score_matrix_16,
				sd->dprofile_w, (WORD*) sd->hearray, sd->target_count, master_targets + sd->target_index, master_scores + sd->target_index,
				master_diffs + sd->target_index, master_alignlengths + sd->target_index, db->query.len, dirbufferbytes / 8, sd->dir_array,
				db);
	else
		searcher.search8(sd->qtable, Property::penalty_gapopen, Property::penalty_gapextend, (BYTE*) Matrix::score_matrix_8, sd->dprofile,
				sd->hearray, sd->target_count, master_targets + sd->target_index, master_scores + sd->target_index,
				master_diffs + sd->target_index, master_alignlengths + sd->target_index, db->query.len, dirbufferbytes / 8, sd->dir_array,
				db);
}

int scanner::search_getwork(unsigned long * countref, unsigned long * firstref) {
	// * countref = how many sequences to search
	// * firstref = index into master_targets/scores/diffs where thread should start

	unsigned long status = 0;

	if (master_next < master_length) {
		unsigned long chunksize = ((master_length - master_next + remainingchunks - 1) / remainingchunks);

		*countref = chunksize;
		*firstref = master_next;

		master_next += chunksize;
		remainingchunks--;
		status = 1;
	}

	return status;
}

void scanner::master_dump() {
	printf("master_dump\n");
	printf("   i    t    s    d\n");
	for (unsigned long i = 0; i < 1403; i++) {
		printf("%4lu %4lu %4lu %4lu\n", i, master_targets[i], master_scores[i], master_diffs[i]);
	}
}

void scanner::search_worker_core() {
	search_init();
	while (search_getwork(&sd->target_count, &sd->target_index))
		search_chunk(master_bits);
}

void scanner::search_do(unsigned long query_no, unsigned long listlength, unsigned long * targets, unsigned long * scores,
		unsigned long * diffs, unsigned long * alignlengths, long bits) {
	db->query.qno = query_no;
	db->query = db->get_sequence_and_length(query_no);

	master_next = 0;
	master_length = listlength;
	master_targets = targets;
	master_scores = scores;
	master_diffs = diffs;
	master_alignlengths = alignlengths;
	master_bits = bits;

	//TODO Thread assumed 1 here, need refactoring
	remainingchunks = 1;

	search_worker_core();
}

void scanner::search_begin() {
	search_alloc();
}

void scanner::set_db(Db_data * db) {
	this->db = db;
}
