#include "cluster.h"

#define SEPCHAR ' '
#define BITS 8

static unsigned long swarmed;
static unsigned long seeded;

cluster_result * algo_run(partition_info partition) {
	unsigned long count_comparisons_8 = 0;
	unsigned long count_comparisons_16 = 0;
	unsigned long amplicons = partition.end - partition.start;

	unsigned long targetcount;
	unsigned long * targetindices = new unsigned long[amplicons];
	unsigned long * targetampliconids = new unsigned long[amplicons];
	unsigned long * scores = new unsigned long[amplicons];
	unsigned long * diffs = new unsigned long[amplicons];
	unsigned long * alignlengths = new unsigned long[amplicons];
	unsigned long * qgramamps = new unsigned long[amplicons];
	unsigned long * qgramdiffs = new unsigned long[amplicons];
	unsigned long * qgramindices = new unsigned long[amplicons];
	unsigned long * hits = new unsigned long[amplicons];

	unsigned long searches = 0;
	unsigned long estimates = 0;

	unsigned long largestswarm = 0;
	unsigned long swarmsize = 0;

	unsigned long maxgenerations = 0;

	ampliconinfo_s * amps = new ampliconinfo_s[amplicons];

	/* set ampliconid for all */
	for (unsigned long i = partition.start; i < partition.end; i++) {
		amps[i].ampliconid = i;
	}

	/* always search in 8 bit mode unless resolution is very high */

	unsigned long bits;

	if ((unsigned long) Property::resolution <= Property::diff_saturation())
		bits = 8;
	else
		bits = 16;

	seeded = 0;
	swarmed = 0;

	unsigned long swarmid = 0;

	progress_init("Clustering:       ", amplicons);
	while (seeded < amplicons) {

		/* process each initial seed */

		swarmid++;

		unsigned long amplicons_copies = 0;
		unsigned long singletons = 0;
		unsigned long hitcount = 0;
		unsigned long diffsum = 0;
		unsigned long maxradius = 0;
		unsigned long maxgen = 0;
		unsigned long seedindex;

		seedindex = seeded;
		seeded++;

		amps[seedindex].swarmid = swarmid;
		amps[seedindex].generation = 0;
		amps[seedindex].radius = 0;

		unsigned long seedampliconid = amps[seedindex].ampliconid;
		hits[hitcount++] = seedampliconid;

		unsigned long abundance = db_getabundance(seedampliconid);
		amplicons_copies += abundance;
		if (abundance == 1)
			singletons++;

		swarmsize = 1;
		swarmed++;

		/* find diff estimates between seed and each amplicon in pool */

		targetcount = 0;

		unsigned long listlen = amplicons - swarmed;

		for (unsigned long i = 0; i < listlen; i++)
			qgramamps[i] = amps[swarmed + i].ampliconid;

		qgram_work_diff(seedampliconid, listlen, qgramamps, qgramdiffs);

		estimates += listlen;

		for (unsigned long i = 0; i < listlen; i++) {
			unsigned poolampliconid = qgramamps[i];
			long diff = qgramdiffs[i];
			amps[swarmed + i].diffestimate = diff;
			if (diff <= Property::resolution) {
				targetindices[targetcount] = swarmed + i;
				targetampliconids[targetcount] = poolampliconid;
				targetcount++;
			}
		}

		if (targetcount > 0) {
			search_do(seedampliconid, targetcount, targetampliconids, scores, diffs, alignlengths, bits);
			searches++;

			if (bits == 8)
				count_comparisons_8 += targetcount;
			else
				count_comparisons_16 += targetcount;

			for (unsigned long t = 0; t < targetcount; t++) {
				unsigned diff = diffs[t];

				if (diff <= (unsigned long) Property::resolution) {
					unsigned i = targetindices[t];

					/* move the target (i) to the position (swarmed)
					 of the first unswarmed amplicon in the pool */

					if (swarmed < i) {
						struct ampliconinfo_s temp = amps[i];
						for (unsigned j = i; j > swarmed; j--) {
							amps[j] = amps[j - 1];
						}
						amps[swarmed] = temp;
					}

					amps[swarmed].swarmid = swarmid;
					amps[swarmed].generation = 1;
					if (maxgen < 1)
						maxgen = 1;
					amps[swarmed].radius = diff;
					if (diff > maxradius)
						maxradius = diff;

					unsigned poolampliconid = amps[swarmed].ampliconid;
					hits[hitcount++] = poolampliconid;

					diffsum += diff;

					abundance = db_getabundance(poolampliconid);
					amplicons_copies += abundance;
					if (abundance == 1)
						singletons++;

					swarmsize++;

					swarmed++;
				}
			}

			while (seeded < swarmed) {

				/* process each subseed */

				unsigned subseedampliconid;
				unsigned subseedradius;

				unsigned long subseedindex;
				unsigned long subseedgeneration;

				subseedindex = seeded;
				subseedampliconid = amps[subseedindex].ampliconid;
				subseedradius = amps[subseedindex].radius;
				subseedgeneration = amps[subseedindex].generation;

				seeded++;

				targetcount = 0;

				unsigned long listlen = 0;
				for (unsigned long i = swarmed; i < amplicons; i++) {
					unsigned long targetampliconid = amps[i].ampliconid;
					if (amps[i].diffestimate <= subseedradius + Property::resolution) {
						qgramamps[listlen] = targetampliconid;
						qgramindices[listlen] = i;
						listlen++;
					}
				}

				qgram_work_diff(subseedampliconid, listlen, qgramamps, qgramdiffs);

				estimates += listlen;

				for (unsigned long i = 0; i < listlen; i++)
					if ((long) qgramdiffs[i] <= Property::resolution) {
						targetindices[targetcount] = qgramindices[i];
						targetampliconids[targetcount] = qgramamps[i];
						targetcount++;
					}

				if (targetcount > 0) {
					search_do(subseedampliconid, targetcount, targetampliconids, scores, diffs, alignlengths, bits);
					searches++;

					if (bits == 8)
						count_comparisons_8 += targetcount;
					else
						count_comparisons_16 += targetcount;

					for (unsigned long t = 0; t < targetcount; t++) {
						unsigned diff = diffs[t];

						if (diff <= (unsigned long) Property::resolution) {
							unsigned i = targetindices[t];

							/* find correct position in list */

							/* move the target (i) to the position (swarmed)
							 of the first unswarmed amplicon in the pool
							 then move the target further into the swarmed
							 but unseeded part of the list, so that the
							 swarmed amplicons are ordered by id */

							unsigned long targetampliconid = amps[i].ampliconid;
							unsigned pos = swarmed;

							while ((pos > seeded) && (amps[pos - 1].ampliconid > targetampliconid)
									&& (amps[pos - 1].generation > subseedgeneration))
								pos--;

							if (pos < i) {
								struct ampliconinfo_s temp = amps[i];
								for (unsigned j = i; j > pos; j--) {
									amps[j] = amps[j - 1];
								}
								amps[pos] = temp;
							}

							amps[pos].swarmid = swarmid;
							amps[pos].generation = subseedgeneration + 1;
							if (maxgen < amps[pos].generation)
								maxgen = amps[pos].generation;
							amps[pos].radius = subseedradius + diff;
							if (amps[pos].radius > maxradius)
								maxradius = amps[pos].radius;

							unsigned poolampliconid = amps[pos].ampliconid;
							hits[hitcount++] = poolampliconid;
							diffsum += diff;

							abundance = db_getabundance(poolampliconid);
							amplicons_copies += abundance;
							if (abundance == 1)
								singletons++;

							swarmsize++;

							swarmed++;
						}
					}
				}
			}
		}

		if (swarmsize > largestswarm)
			largestswarm = swarmsize;

		if (maxgen > maxgenerations)
			maxgenerations = maxgen;

		progress_update(seeded);
	}
	progress_done();

	/* output results */

	cluster_result * result = new cluster_result;
	result->partition_id = partition.threadid + 1;

	if (amplicons > 0) {
		char sep_amplicons;
		char sep_swarms;

		/* native swarm output */
		sep_amplicons = SEPCHAR; /* usually a space */
		sep_swarms = '\n';

		long previd = -1;

		for (unsigned long i = 0; i < amplicons; i++) {
			member_info * member = new member_info;
			member->sequence = seqindex[amps[i].ampliconid];
			member->generation = amps[i].generation;
			member->radius = amps[i].radius;
			member->qgram_diff = amps[i].diffestimate;
			if (amps[i].swarmid != previd) {
				result->new_cluster(amps[i].swarmid);
			}
			result->clusters.back()->cluster_members.push_back(member);
			previd = amps[i].swarmid;
		}
	}

	delete (qgramdiffs);
	delete (qgramamps);
	delete (qgramindices);
	delete (hits);
	delete (alignlengths);
	delete (diffs);
	delete (scores);
	delete (targetindices);
	delete (targetampliconids);
	delete (amps);

	return result;
}

