/// 786

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "const.h"
#include "qualities.h"
#include "arithmetic.h"

/* calculate phred score of phred char c */
double phred (char c, int offset) { 
	return pow (10, -(c - offset) / 10.0);
}

/* initialize quality mapping using statistical analysis
 * data will be obtained from fastq file f
 * this will also find length of single read in fastq file */
void quality_mapping_init (quality_mapping *q, buffered_file *f, int *read_length) { 
	char line[MAXLINE];
	int stat[128];
	memset (stat, 0, sizeof stat);

	// read sample data
	for (int i = 0; i < _quality_sample_lines; i++) {
		if (_interleave == 20) { /* skip first part in interleaved */
			f_gets (f, line, MAXLINE);
			f_gets (f, line, MAXLINE);
			f_gets (f, line, MAXLINE);
			f_gets (f, line, MAXLINE);
		}

		f_gets (f, line, MAXLINE);
		f_gets (f, line, MAXLINE);
		f_gets (f, line, MAXLINE);
		if (!f_gets (f, line, MAXLINE)) 
			break;
		int l = strlen (line);
		for (int j = 0; j < l; j++) {
			stat[line[j]]++;
		}
		(*read_length) = l;

		if (_interleave == 10) { /* skip second part in interleaved */
			f_gets (f, line, MAXLINE);
			f_gets (f, line, MAXLINE);
			f_gets (f, line, MAXLINE);
			f_gets (f, line, MAXLINE);
		}
	}

	q->offset = 64;
	for (int i = 33; i < 64; i++) if (stat[i]) {
		q->offset = 33;
		break;
	}

	// calculate replacements
	for (int c = 0; c < 128; c++)
		q->values[c] = c;

	// if percentage is 0% then keep original qualities
	if (!_quality_lossy_percentage)
		return;

	char already_assigned[128];
	memset (already_assigned, 0, sizeof (already_assigned));

	// Discard all elements with >30% error
	char hlimit;
	for (hlimit = q->offset; 100 * phred (hlimit, q->offset) > 30; hlimit++) {
		q->values[hlimit] = q->offset;
		already_assigned[hlimit] = 1;
	}

	// sort stats by maximum value
	struct { char c; int stat; } bmm[128], t;
	for (int i = 0; i < 128; i++)
		bmm[i].c = i, bmm[i].stat = stat[i];
	for (int i = 0; i < 128; i++) { // sort, slow but used only once (anyway qsort is not faster)
		for (int j = i + 1; j < 128; j++) 
			if (bmm[j].stat > bmm[i].stat || (bmm[j].stat == bmm[i].stat && bmm[j].c < bmm[i].c)) {
				t = bmm[i];
				bmm[i] = bmm[j];
				bmm[j] = t;
			}
	}

	// obtain replacement table
	double percentage = _quality_lossy_percentage / 100.0;
	for (int i=0; i < 128 && bmm[i].stat; i++) {
		if (!already_assigned[bmm[i].c] && stat[bmm[i].c] >= stat[bmm[i].c - 1] && stat[bmm[i].c] >= stat[bmm[i].c + 1]) {
			double er = phred (bmm[i].c, q->offset);

			double total = er;
			char left = bmm[i].c, right = bmm[i].c;
			for (int c = bmm[i].c - 1; c >= 0; c--) {
				total += phred (c, q->offset);
				if (already_assigned[c] || total / (bmm[i].c - c + 1) > er + (er * percentage)) {
					left = c + 1;
					break;
				}
			}
			total = er;
			for (int c = bmm[i].c + 1; c < 128; c++) {
				total += phred(c, q->offset);
				if (already_assigned[c] || total / (c - bmm[i].c + 1) < er - (er * percentage)) {
					right = c - 1;
					break;
				}
			}

			for (char c = left; c <= right; c++) { // MrgChrs
				q->values[c] = bmm[i].c;
				already_assigned[c] = 1;
			}	
		}
	}
}

int output_quality (char* line, char *read, quality_mapping *q, uint8_t *dest, int ZZ) {
	static uint32_t prev[2][3] = { { 500, 500, 500 }, 
		                            { 500, 500, 500 } };

	int bc = 0, l = 0;
	while (line[l]) {
		dest[bc] = (read[l] == 'N' ? q->offset : q->values[line[l]]) - q->offset;
		if (prev[ZZ][1] < 256) {
			ac_freq3[ZZ][prev[ZZ][1]*AC_DEPTH + dest[bc]]++;
			if (prev[ZZ][0] < 256) 
				ac_freq4[ZZ][(prev[ZZ][0]*AC_DEPTH + prev[ZZ][1])*AC_DEPTH+dest[bc]]++;
		}
		else { 
			for (int e = 0; e < AC_DEPTH; e++)
				for (int j = 0; j < AC_DEPTH; j++)
					for (int l = 0; l < AC_DEPTH; l++)
						ac_freq4[ZZ][(e*AC_DEPTH + j)*AC_DEPTH+l] = 1;
		}
		prev[ZZ][0] = prev[ZZ][1];
		prev[ZZ][1] = dest[bc];
		bc++;
		l++;
	}
	return bc;
}


