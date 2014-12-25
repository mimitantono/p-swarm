#include "search.h"

searcher::searcher() {
}

searcher::~searcher() {
}

search_result searcher::search_single(seqinfo_t * query, seqinfo_t * target, search_data * search_data) {
	std::vector<search_result> results;
	search_result sr;
	results.push_back(sr);
	std::vector<queryinfo_t> targets;
	queryinfo_t query_t;
	query_t.len = query->seqlen;
	query_t.seq = query->seq;
	queryinfo_t target_t;
	target_t.len = target->seqlen;
	target_t.seq = target->seq;
	targets.push_back(target_t);
	unsigned long int dirbuffersize = 8 * Property::longest * ((Property::longest + 3) / 4) * 4;
	for (long i = 0; i < query_t.len; i++) {
		search_data->qtable[i] = search_data->dprofile + 64 * query_t.seq[i];
		search_data->qtable_w[i] = search_data->dprofile_w + 32 * query_t.seq[i];
	}

//	fprintf(stderr, "Longest: %lu\n", Property::longest);
//	fprintf(stderr, "Comparing %s[%lu] with %s[%lu]\n", query_t.seq, query_t.len, target_t.seq, target_t.len);

	if (Property::bits == 8) {
		search8(search_data, targets, &results, &query_t, dirbuffersize / 8, Property::longest);
	} else {

		search16(search_data, targets, &results, &query_t, dirbuffersize / 8, Property::longest);
	}
	return results[0];
}
