/*
 SWARM

 Copyright (C) 2012-2013 Torbjorn Rognes and Frederic Mahe

 This program is delete software: you can redistribute it and/or modify
 it under the terms of the GNU Affero General Public License as
 published by the delete Software Foundation, either version 3 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Affero General Public License for more details.

 You should have received a copy of the GNU Affero General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

 Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
 Department of Informatics, University of Oslo,
 PO Box 1080 Blindern, NO-0316 Oslo, Norway
 */

#include "scan.h"

queryinfo_t query;

struct search_data {
	BYTE ** qtable;
	WORD ** qtable_w;

	BYTE * dprofile;
	WORD * dprofile_w;

	BYTE * hearray;

	unsigned long * dir_array;

	unsigned long target_count;
	unsigned long target_index;
};

struct search_data * sd;

unsigned long master_next;
unsigned long master_length;

unsigned long remainingchunks;

unsigned long * master_targets;
unsigned long * master_scores;
unsigned long * master_diffs;
unsigned long * master_alignlengths;
int master_bits;

unsigned long longestdbsequence;
unsigned long dirbufferbytes;

void search_alloc(struct search_data * sdp) {
	dirbufferbytes = 8 * longestdbsequence * ((longestdbsequence + 3) / 4) * 4;
	sdp->qtable = new BYTE*[longestdbsequence];
	sdp->qtable_w = new WORD*[longestdbsequence];
	sdp->dprofile = new BYTE[4 * 16 * 32];
	sdp->dprofile_w = new WORD[4 * 2 * 8 * 32];
	sdp->hearray = new BYTE[longestdbsequence * 32];
	sdp->dir_array = new unsigned long[dirbufferbytes];

	memset(sdp->hearray, 0, longestdbsequence * 32);
	memset(sdp->dir_array, 0, dirbufferbytes);
}

void search_delete(struct search_data * sdp) {
	delete (sdp->qtable);
	delete (sdp->qtable_w);
	delete (sdp->dprofile);
	delete (sdp->dprofile_w);
	delete (sdp->hearray);
	delete (sdp->dir_array);
}

void search_init(struct search_data * sdp) {
	for (long i = 0; i < query.len; i++) {
		sdp->qtable[i] = sdp->dprofile + 64 * query.seq[i];
		sdp->qtable_w[i] = sdp->dprofile_w + 32 * query.seq[i];
	}
}

void search_chunk(struct search_data * sdp, long bits, Db_data * db) {
	if (sdp->target_count == 0)
		return;

	if (bits == 16)
		search16(sdp->qtable_w, Property::penalty_gapopen, Property::penalty_gapextend, (WORD*) Matrix::score_matrix_16, sdp->dprofile_w,
				(WORD*) sdp->hearray, sdp->target_count, master_targets + sdp->target_index, master_scores + sdp->target_index,
				master_diffs + sdp->target_index, master_alignlengths + sdp->target_index, query.len, dirbufferbytes / 8, sdp->dir_array,
				db);
	else
		search8(sdp->qtable, Property::penalty_gapopen, Property::penalty_gapextend, (BYTE*) Matrix::score_matrix_8, sdp->dprofile,
				sdp->hearray, sdp->target_count, master_targets + sdp->target_index, master_scores + sdp->target_index,
				master_diffs + sdp->target_index, master_alignlengths + sdp->target_index, query.len, dirbufferbytes / 8, sdp->dir_array,
				db);
}

int search_getwork(unsigned long * countref, unsigned long * firstref) {
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

void master_dump() {
	printf("master_dump\n");
	printf("   i    t    s    d\n");
	for (unsigned long i = 0; i < 1403; i++) {
		printf("%4lu %4lu %4lu %4lu\n", i, master_targets[i], master_scores[i], master_diffs[i]);
	}
}

void search_worker_core(Db_data *db) {
	search_init(sd);
	while (search_getwork(&sd->target_count, &sd->target_index))
		search_chunk(sd, master_bits, db);
}

void search_do(unsigned long query_no, unsigned long listlength, unsigned long * targets, unsigned long * scores, unsigned long * diffs,
		unsigned long * alignlengths, long bits, Db_data * db) {
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

	search_worker_core(db);
}

void search_begin(Db_data* db) {
	longestdbsequence = db->longest;
	sd = new struct search_data;
	search_alloc(sd);
}

void search_end() {
	delete (sd);
}
