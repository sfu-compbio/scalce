/// 786

#include <time.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <omp.h>

#include "const.h"
#include "buffio.h"
#include "compress.h"
#include "reads.h"
#include "qualities.h"
#include "names.h"
#include "arithmetic.h"

buffered_file file_pool[MAXOPENFILES];  /* handles for files */
int32_t read_length[2];                 /* length of read in one fastq file */
int64_t reads_count = 0;                /* total read count in all input files */

char global_buffer[GLOBALBUFSZ];        /* global shared buffer for reading/writing */
uint8_t metadata[MAXMETA];              /* metadata file contents */

int merhamet_merge (int *fl, int flc, int idx) {
	const int MAXSEE = MAXOPENFILES / 2 - 1; // max

	char buffer[MAXLINE];

	DLOG ("Merging part %d for files + ", idx);
	for (int i = 0; i < flc; i++) {
		struct stat S;
		snprintf (buffer, MAXLINE, "%s/t_%03d_%d.tmp", _temp_directory, fl[i], idx);
		stat (buffer, &S);
		DLOG ("%02d (%7.1lfM)%s", fl[i], (S.st_size/1024.0)/1024.0, (i==flc-1?"":", "));
	}
	DLOG (" ... ");

	snprintf (buffer, MAXLINE, "%s/t_TMP_%d.tmp", _temp_directory, idx);
	f_open (file_pool + 0, buffer, IO_WRITE);
	for (int i = 0; i < flc; i++) {
		snprintf(buffer, MAXLINE, "%s/t_%03d_%d.tmp", _temp_directory, fl[i], 3); // meta
		f_open (file_pool + i + 1, buffer, IO_READ);
		snprintf(buffer, MAXLINE, "%s/t_%03d_%d.tmp", _temp_directory, fl[i], idx); // other
		f_open (file_pool + flc + i + 1, buffer, IO_READ);
	}

	int32_t bins[MAXMERGE];
	int32_t cores[MAXMERGE];
	int minV = MAXBIN, minC = MAXBIN;
	int minI = 0;
	for (int i = 0; i < flc; i++) { 
		f_read (file_pool + i + 1, &(bins[i]), sizeof(int32_t));
		f_read (file_pool + i + 1, &(cores[i]), sizeof(int32_t));
		if (bins[i] < minV) {
			minV = bins[i];
			minC = cores[i];
			minI = i;
		}
	}

	int     binex = flc, metadata_pos = 0;
	int64_t totalLen = 0;
	int32_t prevMinV, prevMinC;
	while (binex) {
		int64_t len;
		for (int i = 0; i <= idx-(idx>3); i++)
			f_read (file_pool + minI + 1, &len, sizeof(int64_t));
		for (int64_t s = 0; s < len; ) {
			int64_t l = f_read (file_pool + flc + minI + 1, global_buffer, MIN (GLOBALBUFSZ, len-s));
			f_write (file_pool + 0, global_buffer, l);
			s += l;
		}
		totalLen += len;
		for (int i = idx + (idx<=3); i < 3 + 2 * _use_second_file; i++)
			f_read (file_pool + minI + 1, &len, sizeof(int64_t));

		prevMinV = minV;
		prevMinC = minC;
		minV = MAXBIN;
		int rd = f_read (file_pool + minI + 1, &(bins[minI]), sizeof(int32_t));
		if (rd != sizeof(int32_t)) {
			bins[minI] = MAXBIN;
			binex--;
		}
		else f_read( file_pool + minI + 1, &(cores[minI]), sizeof(int32_t) );
		for (int i = 0; i < flc; i++) 
			if (bins[i] < minV) {
				minV = bins[i];
				minC = cores[i];
				minI = i;
			}
		if (minV != prevMinV) {
			memcpy (metadata + metadata_pos, &prevMinV, sizeof(int32_t));
			metadata_pos += sizeof(int32_t);
			memcpy (metadata + metadata_pos, &prevMinC, sizeof(int32_t));
			metadata_pos += sizeof(int32_t) + sizeof(int64_t) * (idx - (idx > 3));
			memcpy (metadata + metadata_pos, &totalLen, sizeof(int64_t));
			metadata_pos += sizeof(int64_t) *  (3 + 2 * _use_second_file - idx + (idx > 3));
			totalLen = 0;
		}
	}

	f_close (file_pool + 0);
	for (int i = 0; i < flc; i++) {
		f_close(file_pool + i + 1);
		snprintf(buffer, MAXLINE, "%s/t_%03d_%d.tmp", _temp_directory, fl[i], idx); 
		f_close(file_pool + flc + i + 1);
		unlink (buffer);
	}

	DLOG("OK\n");
	return metadata_pos;
}

void WW(uint32_t u,buffered_file *f) {
	do {
		uint8_t x=u&0x7f; u>>=7;
		if(u) x|=0x80;
		f_write(f,&x,1);
	}while(u);
}

int64_t combine_and_compress_with_split (int fl, const char *output, int64_t split_size, int64_t phred_offset) {
	LOG ("Generating %lld scalce file(s) from temp #%d ", reads_count / split_size + (reads_count % split_size != 0), fl);
	if (split_size < reads_count) 
		LOG ("with max. %lld reads, ", split_size);
	if (_use_second_file)
		LOG ("using paired files, ");

	char buffer[MAXLINE];
	for (int i = 0; i <= 3 + 2 * _use_second_file; i++) { 
		snprintf(buffer, MAXLINE, "%s/t_%03d_%d.tmp", _temp_directory, fl, i);
		f_open (file_pool + i, buffer, IO_READ);
	}
	
	buffered_file *fN = file_pool + 0,
	              *fR = file_pool + 1,
	              *fQ = file_pool + 2,
	              *fM = file_pool + 3;

	/*
	 * Structure:
	 * - Read size 		int32
	 * - Num reads			int64
	 * - Phred off     	int64
	 * - Qual off    		int64
	 * - Name off    		int64
	 * - Reads				*
	 * - Quals				*
	 * - NamesMode			char
	 *   ?1	Names 		*
	 *   ?0	StartIdx 	int64
	 *   		LibName		*
	 */

	buffered_file foQ, foR, foN;

	f_init (&foQ, IO_SYS);
	f_init (&foR, _compression_mode);
	f_init (&foN, _compression_mode);
	LOG ("compression: %s\n", 
			_compression_mode == IO_SYS 
				? "no compression"
				: (_compression_mode == IO_PGZIP
						? "pigz"
						: (_compression_mode == IO_GZIP
								? "gzip"
								: "bzip"
						)
				)
	);

	int64_t new_size = 0;

	uint8_t magic[] = { 's','c','a','l', 'c','e','2','0' };

	uint64_t size;
	for (int NP = 0; NP < _use_second_file + 1; NP++) {
		int32_t id, core;
		int64_t lR = 0, lQ = 0, lN = 0, temp;
		//	char alive = 1;
		//	for (int64_t F = 0; alive; F++) {
		//		int64_t limit = MIN (split_size, reads_count - F * split_size);

		snprintf (buffer, MAXLINE, "%s_%d.scalcen", output, NP+1);
		f_open (&foN, buffer, IO_WRITE);
		snprintf (buffer, MAXLINE, "%s_%d.scalceq", output, NP+1);
		f_open (&foQ, buffer, IO_WRITE);
		snprintf (buffer, MAXLINE, "%s_%d.scalcer", output, NP+1);
		f_open (&foR, buffer, IO_WRITE);

		f_write (&foR, magic, 8);
		f_write (&foR, read_length + NP, sizeof(int32_t));

		f_write (&foQ, magic, 8);
		f_write (&foQ, &phred_offset, sizeof(int64_t));
		for (int i = 0; i < AC_DEPTH; i++) // write stat stuff
			for (int j = 0; j < AC_DEPTH; j++) 
				for (int k = 0; k < AC_DEPTH; k++) 
					f_write(&foQ, &ac_freq4[NP][(i*AC_DEPTH + j)*AC_DEPTH+k], sizeof(uint32_t));	
		set_ac_stat(ac_freq3[NP], ac_freq4[NP]);
		uint64_t qualtotalsize = reads_count * read_length[NP];
		f_write(&foQ,&qualtotalsize,sizeof(uint64_t));

		f_write (&foN, magic, 8);
		f_write (&foN, &_use_names, 1);
		if (!_use_names) {
			int64_t name_idx = 0; // F * split_size;
			f_write (&foN, &name_idx, sizeof(int64_t));
			f_write (&foN, _library_name, strlen(_library_name));
		}
		else {
			f_write(&foN, &sanger, 1);
			if (sanger) {
				int l; char *cx;
			
				cx = prefix;
				l = strlen(cx); f_write(&foN, &l, sizeof(int)); f_write(&foN, cx, l);
				cx = machine_name;
				l = strlen(cx); f_write(&foN, &l, sizeof(int)); f_write(&foN, cx, l);
				cx = run_id;
				l = strlen(cx); f_write(&foN, &l, sizeof(int)); f_write(&foN, cx, l);
			}
		}


		/* read metadata */
		while(1){
			if (f_read (fM, &id, sizeof(int32_t)) != sizeof(int32_t))
				break;
			f_read (fM, &core, sizeof(int32_t));
			//LOG("%d\n",core);
			f_read (fM, &lN, sizeof(int64_t));
			f_read (fM, &lR, sizeof(int64_t));
			f_read (fM, &lQ, sizeof(int64_t));
			if (_use_second_file) {
				int64_t *lR2 = (NP ? &lR : &temp);
				int64_t *lQ2 = (NP ? &lQ : &temp);
				f_read (fM, lR2, sizeof(int64_t));
				f_read (fM, lQ2, sizeof(int64_t));
			}

			if (!NP) {
				f_write(&foR, &core, sizeof(int32_t));
				uint64_t size;
				if (core != MAXBIN - 1) // root or non root
					size = lR / ( SZ_READ( read_length[0] - strlen(patterns[core]) ) + sizeof(int16_t) );
				else
					size = lR / ( SZ_READ(read_length[0])+sizeof(int16_t) );
				f_write(&foR, &size, sizeof(int64_t));
			}
			for (int64_t i = 0; i < lR; i += GLOBALBUFSZ) {
				int64_t r = f_read (fR, global_buffer, MIN (lR - i, GLOBALBUFSZ));
				f_write (&foR, global_buffer, r);
			}
			for (int64_t i = 0; i < lQ; i += GLOBALBUFSZ) {
				int64_t r = f_read (fQ, global_buffer, MIN (lQ - i, GLOBALBUFSZ));
				ac_write (&foQ, (uint8_t*)global_buffer, r);
			}
			if (_use_names && sanger) {
				uint32_t n_id = 0, n_lane = 0;
				for (int64_t nn = 0; nn < lN; ) {
					uint8_t ll;
					f_read(fN, &ll, sizeof(uint8_t));
					f_read(fN, global_buffer, ll);
					global_buffer[ll]=0;
					nn += ll + 1;

					char cc[200], cp = 0, ep = 0;
					for (int i = 0; i <= ll; i++) if (i == ll || global_buffer[i] == ' ' || global_buffer[i] == ':' || global_buffer[i] == '.' || global_buffer[i] == '/') {
						cc[cp]=0;
						if (ep == 1) {
							WW(atoi(cc)-n_id,&foN);
							if ( core != MAXBIN-1 ) n_id=atoi(cc);
						}
						else if (ep == 4) {
							WW(atoi(cc)-n_lane,&foN);
							if(core != MAXBIN-1) n_lane=atoi(cc);
						}
						else if (ep == 5) {
							WW(atoi(cc),&foN);	
						}
						else if (ep == 6) {
							WW(atoi(cc),&foN);	
						}

						cp = 0;
						ep++;
					} else cc[cp++]=global_buffer[i];

					//int64_t r = f_read (fN, global_buffer, MIN (lN - i, GLOBALBUFSZ));
					//f_write (&foN, global_buffer, r);
				}
			}
			else if (_use_names) {
				for (int64_t i = 0; i < lN; i += GLOBALBUFSZ) {
					int64_t r = f_read (fN, global_buffer, MIN (lN - i, GLOBALBUFSZ));
					f_write (&foN, global_buffer, r);
				}
			}
			else;
		}

		ac_write(&foQ, 0, 0);
		f_close (&foR);
		f_close (&foN);
		f_close (&foQ);
		DLOG("\tData for paired-end %d with %lld reads in files (%s, %s, %s)\n", NP, reads_written, foN.file_name, foQ.file_name, foR.file_name);

		// ...
		struct stat s;
		stat (foR.file_name, &s); new_size += s.st_size;
		stat (foQ.file_name, &s); new_size += s.st_size;
		stat (foN.file_name, &s); new_size += s.st_size;

		if (_use_second_file && !NP) {
			fR = file_pool + 4;
			fQ = file_pool + 5;
			f_seek (fN, 0);
			f_seek (fM, 0);
		}
	}

	DLOG("\tCleaning ...\n");
	f_free (&foQ);
	f_free (&foN);
	f_free (&foR);
	for (int i = 0; i <= 3 + 2 * (_use_second_file); i++) { 
		snprintf(buffer, MAXLINE, "%s/t_%03d_%d.tmp", _temp_directory, fl, i);
		f_close (file_pool + i);
		unlink (buffer);
	}

	DLOG("\tDone!\n");
	return new_size;
}

void merge (int total, int param) {
	char buffer[MAXLINE], buffer2[MAXLINE];

	char *used = (char*) mallox (total + 5);
	memset (used, 0, total + 5);

	int found[50];
	for (;;) {
		int totalFound = 0;
		int fc = 0;
		for (int i = 0; i < total; i++) {
			if (!used[i]) {
				found[fc++] = i;
				used[i] = 1;
				totalFound++;
			}
			if (fc == param || (i == total - 1)) {
				if (fc > 1) {
					int mm = 0;
					for (int f = 0; f < 4 + _use_second_file * 2; f++) if (f != 3) {
						mm = merhamet_merge (found, fc, f);
						snprintf (buffer, MAXLINE, "%s/t_TMP_%d.tmp",  _temp_directory, f);
						snprintf (buffer2, MAXLINE, "%s/t_%03d_%d.tmp", _temp_directory, found[fc - 1], f);
						unlink(buffer2);
						rename(buffer,buffer2);
					}
					for (int f = 0; f < fc; f++) {
						snprintf(buffer, MAXLINE, "%s/t_%03d_%d.tmp", _temp_directory, found[f], 3); 
						unlink(buffer);
					}
					snprintf(buffer, MAXLINE, "%s/t_%03d_%d.tmp", _temp_directory, found[fc - 1], 3);
					f_open(file_pool + 3, buffer, IO_WRITE);
					f_write(file_pool + 3, metadata, mm);
					f_close(file_pool + 3);
				}
				used[found[fc - 1]] = 0;
				fc = 0;
			}
		}
		if (totalFound == 1) 
			break;
	}

	frex (used, total + 5);
}

void dump_trie (int fl, aho_trie *t) {
	char buffer[MAXLINE];
	for (int i = 0; i < 4 + (_use_second_file * 2); i++) {
		snprintf(buffer, MAXLINE, "%s/t_%03d_%d.tmp", _temp_directory, fl, i);
		f_open (file_pool + i, buffer, IO_WRITE);	
	}
	aho_output (t, file_pool);
	for (int i = 0; i < 4 + (_use_second_file * 2); i++) 
		f_close (file_pool + i);
	LOG("\tCreated temp file %d\n", fl);
}

void get_quality_stats (buffered_file *f, char *path, quality_mapping *qmap) {
	f_open (f, path, IO_READ);
	if (!f_alive (f)) ERROR ("Cannot open file %s!\n", path);

	if (_interleave) _interleave = 10;
	quality_mapping_init (qmap, f, read_length);
	if (_interleave) _interleave = 1;
	
	if (_interleave) {
		_interleave = 20;
		quality_mapping_init (qmap + 1, f, read_length + 1);
		_interleave = 1;
	}
	else if (_use_second_file) {
		f_close (f);
		f_open (f, get_second_file (path), IO_READ);
		if (!f_alive (f)) ERROR ("Cannot open file %s!\n", get_second_file (path));
		quality_mapping_init (qmap + 1, f, read_length + 1);
	}
	f_close (f);

	for (int i = 0; i < (_use_second_file | _interleave) + 1; i++)
		LOG("\tPaired end #%d, quality offset: %d\n"
			 "\t               read length: %d\n", i+1, qmap[i].offset, read_length[i]);

}

void compress (char **files, int nf, const char *output, const char *pattern_path) {
	char line[MAXLINE], read[2][MAXLINE];

	int start_time = TIME;

	LOG("Buffer size: %lldK, bucket storage size: %lldK\n", _file_buffer_size/1024, _max_bucket_set_size/1024);


	/* initialize the trie */
	DLOG ("Reading core strings ... ");
	aho_trie *trie = pattern_path[0] ? read_patterns_from_file(pattern_path) : read_patterns ();
	DLOG ("OK\n");


	/* initialize basic variables */
	LOG ("Preprocessing FASTQ files ...\n");
	uint64_t total_size = 0;     /* current working set size in bytes */
	int temp_file_count = 0;    
	int64_t original_size = 0;   /* size of input files, combined */

	char use_only_second_file = _use_second_file;        /* we need original _use... variable for simplification */
	if (_interleave)             /* all auxiliray functions will behave as in second file mode */
		_use_second_file = 1;

	aho_trie **buckets = (aho_trie**) mallox (_thread_count * sizeof(aho_trie*));
	char ****data = (char****) mallox (_thread_count * sizeof(char*));
	read_data *rdat = (read_data*) mallox (_thread_count * sizeof (read_data));
	for(int t = 0; t < _thread_count; t++) {
		rdat[t].data = (uint8_t*) mallox (5 * MAXLINE);
		data[t] = (char***) mallox (2 * sizeof(char*));
		data[t][0] = (char**) mallox (3 * sizeof (char*));
		data[t][1] = (char**) mallox (3 * sizeof (char*));
		for (int i=0; i<6; i++)
			data[t][i/3][i%3] = (char*) mallox (MAXLINE * sizeof (char));
	}

	/* initialize the files */
	for (int i = 0; i < MAXOPENFILES; i++)
		f_init (file_pool + i, IO_SYS);
	buffered_file f[2];
	f_init (f + 0, IO_GZIP);
	if (use_only_second_file)
		f_init (f + 1, IO_GZIP);

	/* obtain quality tables and stats */
	quality_mapping qmap[2];
	get_quality_stats (f, files[0], qmap);


	/* iterate through files */
	for (int i = 0; i < nf; i++) {
		/* open files */
		for (int fi = 0; fi < use_only_second_file + 1; fi++) {
			char *fn = !fi ? files[i] : get_second_file (files[i]); 

			f_open (f + fi, fn, IO_READ); 
			if (!f_alive (f + fi)) 
				ERROR ("Cannot read file %s\n", fn);
			
			struct stat s;
			stat (f[fi].file_name, &s);
			original_size += s.st_size;
		}

		int64_t file_reads = 0;  /* read count for this file */


		LOG("Memory alloc'd: %.2lf\n", getmemx());


		int d_idx = 0; 
		char done = 0;
		while (1) {
			if (!f_gets (f + 0, data[d_idx][0][0], MAXLINE)) {
				done = 1;
			}
			if (!done) {
				f_gets (f + 0, data[d_idx][0][1], MAXLINE);
				f_gets (f + 0, data[d_idx][0][2], MAXLINE);
				f_gets (f + 0, data[d_idx][0][2], MAXLINE);
				buffered_file *fn = (_interleave ? f : (_use_second_file ? f + 1 : 0));
				if (fn) {
					f_gets (fn, data[d_idx][1][0], MAXLINE);
					f_gets (fn, data[d_idx][1][1], MAXLINE);
					f_gets (fn, data[d_idx][1][2], MAXLINE);
					f_gets (fn, data[d_idx][1][2], MAXLINE);
				}

				d_idx++;
			}

			if (d_idx == _thread_count || done) {
				#pragma omp parallel for reduction(+:total_size) num_threads(_thread_count) 
				for (int t = 0; t < d_idx; t++) {
					uint8_t *out = rdat[t].data;

					int n = aho_search (data[t][0][1], trie, &(buckets[t]));

					int64_t sz = output_name (data[t][0][0], out);
					if (n != -1) {
						sz += output_read (data[t][0][1], out + sz, n-buckets[t]->level+1, buckets[t]->level);
						rdat[t].end = n-buckets[t]->level+1;
					}
					else {
						sz += output_read (data[t][0][1], out + sz, 0, 0);
						rdat[t].end = 0; //read_length[0];
					}
					sz += output_quality(data[t][0][2], data[t][0][1], qmap + 0, out + sz,0);
					int64_t of2 = sz;
					if (_use_second_file || _interleave) {
						sz += output_read    (data[t][1][1], out + sz, 0, 0);
						sz += output_quality (data[t][1][2], data[t][1][1], qmap + 1, out + sz,1);
					}

					rdat[t].of = of2;
					rdat[t].sz = sz;
					aho_trie_bucket (buckets[t], rdat + t);

					total_size += sz + sizeof(read_data) + sizeof(struct bin_node*);
				}
				
				reads_count += d_idx;
				file_reads  += d_idx;
				d_idx = 0;

				if (total_size >= _max_bucket_set_size) {
					total_size = 0;
			//LOG("> Memory alloc'd: %.2lf\n", getmemx());
					dump_trie (temp_file_count++, trie);
		//	LOG("< Memory alloc'd: %.2lf\n", getmemx());
				}
			}

			if (done) break;
		}

		for (int fi = 0; fi < use_only_second_file + 1; fi++)
			f_close (f + fi);
		LOG("\tDone with file %s, %lld reads found\n", files[i], file_reads);
		if (use_only_second_file)
			LOG("\t     also file %s with %lld reads\n", get_second_file (files[i]), file_reads); 
	}

	/* clean all stuff */
//	LOG("> Memory alloc'd: %.2lf\n", getmemx());
	if (total_size)
		dump_trie (temp_file_count++, trie);
//	LOG("< Memory alloc'd: %.2lf\n", getmemx());
	for (int fi = 0; fi < use_only_second_file + 1; fi++)
		f_free (f + fi);
	DLOG("\tDone!\n");



	/* merge */
	LOG("Merging results ...\n");
	merge (temp_file_count, 8);

	/* compress and output */
	LOG("Doing something crazy...\n");
	int scalce_preprocessing_time = TIME;
	int64_t new_size = combine_and_compress_with_split (temp_file_count-1, output, (_split_reads ? _split_reads : reads_count), qmap[0].offset);

	/* cleaning again */
	LOG("Cleaning ...\n");
	aho_trie_free (trie);
	remove(_temp_directory);
	for (int i = 0; i < MAXOPENFILES; i++)
		f_free (file_pool + i);

	/* print useful stats */
	LOG("Statistics:\n");
	LOG("\tTotal number of reads: %lld\n", reads_count);
	LOG("\tRead length: first end %d\n", read_length[0]);
	if (_use_second_file) 
		LOG("\t             second end %d\n", read_length[1]);
	LOG("\tUnbucketed reads count: %d, bucketed percentage %.2lf\n", unbuck(), 100.0 * ((double)(reads_count - unbuck()) / reads_count));
	LOG("\tLossy percentage: %d\n", _quality_lossy_percentage);
	LOG("\tTime elapsed: %02d:%02d:%02d\n", TIME/3600, (TIME/60)%60, TIME%60);
	LOG("\tCompression time: %02d:%02d:%02d\n", (TIME-scalce_preprocessing_time)/3600, ((TIME-scalce_preprocessing_time)/60)%60, (TIME-scalce_preprocessing_time)%60);
	LOG("\tOriginal size: %.2lfM, new size: %.2lfM, compression factor: %.2lf\n", original_size/(1024.0*1024.4), new_size/(1024.0*1024.0), (original_size/(double)new_size));
}

