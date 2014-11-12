#include "scan.h"

scanner::scanner() {
	searcher.set_query(&query);
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
}

void scanner::search_alloc(struct search_data * sdp) {
	dirbufferbytes = 8 * db->longest * ((db->longest + 3) / 4) * 4;
	sdp->qtable = new BYTE*[db->longest];
	sdp->qtable_w = new WORD*[db->longest];
	sdp->dprofile = new BYTE[4 * 16 * 32];
	sdp->dprofile_w = new WORD[4 * 2 * 8 * 32];
	sdp->hearray = new BYTE[db->longest * 32];
	sdp->dir_array = new unsigned long[dirbufferbytes];

	memset(sdp->hearray, 0, db->longest * 32);
	memset(sdp->dir_array, 0, dirbufferbytes);
}

void scanner::search_delete(struct search_data * sdp) {
	delete (sdp->qtable);
	delete (sdp->qtable_w);
	delete (sdp->dprofile);
	delete (sdp->dprofile_w);
	delete (sdp->hearray);
	delete (sdp->dir_array);
}

void scanner::search_init(struct search_data * sdp) {
	for (long i = 0; i < query.len; i++) {
		sdp->qtable[i] = sdp->dprofile + 64 * query.seq[i];
		sdp->qtable_w[i] = sdp->dprofile_w + 32 * query.seq[i];
	}
}

void scanner::search_chunk(struct search_data * sdp, long bits) {
	if (sdp->target_count == 0)
		return;

	if (bits == 16)
		searcher.search16(sdp->qtable_w, Property::penalty_gapopen, Property::penalty_gapextend, (WORD*) Matrix::score_matrix_16,
				sdp->dprofile_w, (WORD*) sdp->hearray, sdp->target_count, master_targets + sdp->target_index,
				master_scores + sdp->target_index, master_diffs + sdp->target_index, master_alignlengths + sdp->target_index, query.len,
				dirbufferbytes / 8, sdp->dir_array, db);
	else
		searcher.search8(sdp->qtable, Property::penalty_gapopen, Property::penalty_gapextend, (BYTE*) Matrix::score_matrix_8, sdp->dprofile,
				sdp->hearray, sdp->target_count, master_targets + sdp->target_index, master_scores + sdp->target_index,
				master_diffs + sdp->target_index, master_alignlengths + sdp->target_index, query.len, dirbufferbytes / 8, sdp->dir_array,
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
	search_init(sd);
	while (search_getwork(&sd->target_count, &sd->target_index))
		search_chunk(sd, master_bits);
}

void scanner::search_do(unsigned long query_no, unsigned long listlength, unsigned long * targets, unsigned long * scores,
		unsigned long * diffs, unsigned long * alignlengths, long bits) {
	query.qno = query_no;
	db->get_sequence_and_length(query_no, &query.seq, &query.len);

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
	search_alloc(sd);
}

void scanner::set_db(Db_data * db) {
	this->db = db;
}
