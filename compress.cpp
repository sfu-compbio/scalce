/*
 * Copyright (c) 2011 - 2012, Simon Fraser University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, 
 * are permitted provided that the following conditions are met:
 *   
 * Redistributions of source code must retain the above copyright notice, this list
 * of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or other
 *   materials provided with the distribution.
 * - Neither the name of the Simon Fraser University nor the names of its contributors may be
 *   used to endorse or promote products derived from this software without specific
 *   prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * Author         : Ibrahim Numanagic
 * Email          : inumanag AT sfu DOT ca
 * Last Update    : 25. vii 2012.
 */

#include <time.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <pthread.h>

#include "const.h"
#include "buffio.h"
#include "compress.h"
#include "reads.h"
#include "qualities.h"
#include "names.h"
#include "arithmetic.h"

//buffered_file file_pool[MAXOPENFILES];  /* handles for files */
int32_t read_length[2];                 /* length of read in one fastq file */
int64_t reads_count = 0;                /* total read count in all input files */

char global_buffer[GLOBALBUFSZ];        /* global shared buffer for reading/writing */
uint8_t metadata[MAXMETA];              /* metadata file contents */

int merhamet_merge (int flc, int idx) {
	//const int MAXSEE = MAXOPENFILES / 2 - 1; // max
	char buffer[MAXLINE];

	DLOG ("Merging part %d for files + ", idx);
	for (int i = 0; i < flc; i++) {
		struct stat S;
		snprintf (buffer, MAXLINE, "%s/t_%03d_%d.tmp", _temp_directory, i, idx);
		stat (buffer, &S);
		DLOG ("%02d (%7.1lfM)%s", i, (S.st_size/1024.0)/1024.0, (i==flc-1?"":", "));
	}
	DLOG (" ... ");

	buffered_file *file_pool = (buffered_file*)mallox((2 * flc + 1) * sizeof(buffered_file));
	snprintf (buffer, MAXLINE, "%s/t_TMP_%d.tmp", _temp_directory, idx);
	f_init (file_pool + 0, IO_SYS);
	f_open (file_pool + 0, buffer, IO_WRITE);
	for (int i = 0; i < flc; i++) {
		snprintf(buffer, MAXLINE, "%s/t_%03d_%d.tmp", _temp_directory, i, 3); // meta
		f_init (file_pool + i + 1, IO_SYS);
		f_open (file_pool + i + 1, buffer, IO_READ);
		snprintf(buffer, MAXLINE, "%s/t_%03d_%d.tmp", _temp_directory, i, idx); // other
		f_init (file_pool + flc + i + 1, IO_SYS);
		f_open (file_pool + flc + i + 1, buffer, IO_READ);
	}

	int32_t *bins = (int32_t*)mallox(flc * sizeof(int32_t));
	int32_t *cores = (int32_t*)mallox(flc * sizeof(int32_t));
//	int32_t bins[MAXMERGE];
//	int32_t cores[MAXMERGE];
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
		snprintf(buffer, MAXLINE, "%s/t_%03d_%d.tmp", _temp_directory, i, idx); 
		f_close(file_pool + flc + i + 1);
		unlink (buffer);
	}
	free(bins);
	free(cores);
	free(file_pool);

	DLOG("OK\n");
	return metadata_pos;
}

int64_t combine_and_compress_with_split (int fl, const char *output, int64_t split_size, int64_t phred_offset) {
	LOG ("Generating %lld scalce file(s) from temp #%d ", reads_count / split_size + (reads_count % split_size != 0), fl);
	if (split_size < reads_count) 
		LOG ("with max. %lld reads, ", split_size);
	if (_use_second_file)
		LOG ("using paired files, ");

	buffered_file file_pool[6];
	char buffer[MAXLINE];
	for (int i = 0; i <= 3 + 2 * _use_second_file; i++) { 
		snprintf(buffer, MAXLINE, "%s/t_%03d_%d.tmp", _temp_directory, fl, i);
		f_init (file_pool + i, IO_SYS);
		f_open (file_pool + i, buffer, IO_READ);
	}
	
	buffered_file *fN = file_pool + 0,
	              *fR = file_pool + 1,
	              *fQ = file_pool + 2,
	              *fM = file_pool + 3;

	if (_compress_qualities && !_no_ac) 
		ac_init();	
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

	f_init (&foQ, _no_ac ? _compression_mode : IO_SYS);
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
	uint8_t magic[] = { 's','c','a','l', 'c','e','2','2' };

	int sz_meta = 1;
	if (read_length[0] > 255) sz_meta = 2;

	uint64_t size;
	for (int NP = 0; NP < _use_second_file + 1; NP++) {
		int32_t id, core;
		int64_t lR = 0, lQ = 0, lN = 0, temp;

		snprintf (buffer, MAXLINE, "%s_%d.scalcen", output, NP+1);
		f_open (&foN, buffer, IO_WRITE);
		snprintf (buffer, MAXLINE, "%s_%d.scalceq", output, NP+1);
		f_open (&foQ, buffer, IO_WRITE);
		snprintf (buffer, MAXLINE, "%s_%d.scalcer", output, NP+1);
		f_open (&foR, buffer, IO_WRITE);

		f_write (&foR, magic, 8);
		f_write (&foR, &_no_ac, sizeof(int));
		f_write (&foR, read_length + NP, sizeof(int32_t));

		f_write (&foQ, magic, 8);
		f_write (&foQ, &phred_offset, sizeof(int64_t));
		
		if (_compress_qualities && !_no_ac) {
			uint64_t _lmt_ = (1LL << 32) - 1;
			int factor = 1;
			factor += (reads_count * read_length[NP]) / _lmt_;
			LOG("shrink factor for %d: %d\n", NP, factor);

			for (int i = 0; i < AC_DEPTH; i++) // write stat stuff
				for (int j = 0; j < AC_DEPTH; j++) 
					for (int k = 0; k < AC_DEPTH; k++) {

						uint64_t *p = &ac_freq4[NP][(i*AC_DEPTH+j)*AC_DEPTH+k];
						*p = *p / factor;
						if (*p == 0) *p = 1;

						assert((*p)>0);
						assert((*p)<((1LL<<32)-1));

						uint32_t x = *p;
						f_write(&foQ, &x, sizeof(uint32_t));
					}
			set_ac_stat(ac_freq3[NP], ac_freq4[NP]);
		
			uint64_t qualtotalsize = reads_count * read_length[NP];
			f_write(&foQ,&qualtotalsize,sizeof(uint64_t));
		}

		f_write (&foN, magic, 8);
		f_write (&foN, &_use_names, 1);
		if (!_use_names) {
			int64_t name_idx = 0; // F * split_size;
			f_write (&foN, &name_idx, sizeof(int64_t));
			f_write (&foN, _library_name, strlen(_library_name));
		}
		

		/* read metadata */
		int __nC__=0;
		while (1) {
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
					size = lR / ( SZ_READ( read_length[0] - strlen(patterns[core]) ) + sz_meta );
				else
					size = lR / ( SZ_READ(read_length[0])+sz_meta );
				f_write(&foR, &size, sizeof(int64_t));
			}
			for (int64_t i = 0; i < lR; i += GLOBALBUFSZ) {
				int64_t r = f_read (fR, global_buffer, MIN (lR - i, GLOBALBUFSZ));
				f_write (&foR, global_buffer, r);
			}
			if (_compress_qualities)
			for (int64_t i = 0; i < lQ; i += GLOBALBUFSZ) {
				int64_t r = f_read (fQ, global_buffer, MIN (lQ - i, GLOBALBUFSZ));
				if (!_no_ac) 
					ac_write (&foQ, (uint8_t*)global_buffer, r);
				else
					f_write (&foQ, (uint8_t*)global_buffer, r);
			}
			if (_use_names) {
				for (int64_t i = 0; i < lN; i += GLOBALBUFSZ) {
					int64_t r = f_read (fN, global_buffer, MIN (lN - i, GLOBALBUFSZ));
					f_write (&foN, global_buffer, r);
				}
			}
			else;

			__nC__++;
		}
		printf("Written %d cores\n", __nC__);

		if (!_no_ac) 
			ac_write(&foQ, 0, 0);
		f_close (&foR);
		f_close (&foN);
		f_close (&foQ);
//		DLOG("\tData for paired-end %d with %lld reads in files (%s, %s, %s)\n", NP, reads_written, foN.file_name, foQ.file_name, foR.file_name);

		// ...
		struct stat s;
		stat (foR.file_name, &s); new_size += s.st_size;
		LOG("\tRead bit size:    paired end %d = %.2lf\n", NP, (s.st_size * 8.0) / double(reads_count * read_length[NP]));
		stat (foQ.file_name, &s); new_size += s.st_size;
		LOG("\tQuality bit size: paired end %d = %.2lf\n", NP, (s.st_size * 8.0) / double(reads_count * read_length[NP]));
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
	if (!_no_ac) 
		ac_finish();

	DLOG("\tDone!\n");
	return new_size;
}

void merge (int total) {
	char buffer[MAXLINE], buffer2[MAXLINE];
	int mm;
	for (int f = 0; f < 4 + _use_second_file * 2; f++) if (f != 3) {
		mm = merhamet_merge (total, f);
		snprintf (buffer, MAXLINE, "%s/t_TMP_%d.tmp",  _temp_directory, f);
		snprintf (buffer2, MAXLINE, "%s/t_%03d_%d.tmp", _temp_directory, total, f);
		unlink(buffer2);
		rename(buffer,buffer2);
	}
	for (int i = 0; i < total; i++) {
		snprintf(buffer, MAXLINE, "%s/t_%03d_%d.tmp", _temp_directory, i, 3); 
		unlink(buffer);
	}
	snprintf(buffer, MAXLINE, "%s/t_%03d_%d.tmp", _temp_directory, total, 3);
	buffered_file f;
	f_init(&f, IO_SYS);
	f_open(&f, buffer, IO_WRITE);
	f_write(&f, metadata, mm);
	f_close(&f);
}

void dump_trie (int fl, aho_trie *t) {
	char buffer[MAXLINE];

	buffered_file file_pool[6];
	for (int i = 0; i < 4 + (_use_second_file * 2); i++) {
		snprintf(buffer, MAXLINE, "%s/t_%03d_%d.tmp", _temp_directory, fl, i);
		f_init (&file_pool[i], IO_SYS);
		f_open (&file_pool[i], buffer, IO_WRITE);	
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
//		if (!f_alive (f)) ERROR ("Cannot open file %s!\n", get_second_file (path));
		quality_mapping_init (qmap + 1, f, read_length + 1);
	}
	f_close (f);

	for (int i = 0; i < (_use_second_file | _interleave) + 1; i++)
		LOG("\tPaired end #%d, quality offset: %d\n"
			 "\t               read length: %d\n", i+1, qmap[i].offset, read_length[i]);
}

/************************/

static pthread_t *threads;
static pthread_spinlock_t r_spin;
static pthread_spinlock_t w_spin;
static int *thread_index;

aho_trie *trie;
quality_mapping qmap[2];
buffered_file input[2];

int temp_file_count = 0;    
uint64_t total_size, file_reads;

void *thread (void *vt) {
	int t = *(int*)vt;

	char read[MAXLINE],
		  name[MAXLINE],
		  qual[MAXLINE],
		  read2[MAXLINE],
		  qual2[MAXLINE];

	aho_trie *bucket;
	read_data rd;
	uint8_t out[MAXLINE*5];
	rd.data = out;

	int txx;

	while (1) {
		pthread_spin_lock(&r_spin);
		if (!f_gets(input, name, MAXLINE)) {
			pthread_spin_unlock(&r_spin);
			return 0;
		}
		f_gets(input, read, MAXLINE);
		f_gets(input, qual, MAXLINE);
		f_gets(input, qual, MAXLINE);
		buffered_file *fn = (_interleave ? input : (_use_second_file ? input + 1 : 0));
		if (fn) {
			f_gets(fn, read2, MAXLINE);
			f_gets(fn, read2, MAXLINE);
			f_gets(fn, qual2, MAXLINE);
			f_gets(fn, qual2, MAXLINE);
		}	
		pthread_spin_unlock(&r_spin);

		int n = aho_search(read, trie, &bucket);
		
		rd.sz = output_name(name, rd.data);
		if (n != -1) {
			rd.sz += output_read(read, rd.data + rd.sz, n - bucket->level + 1, bucket->level);
			rd.end = n + 1;
		}
		else {
			rd.sz += output_read(read, rd.data + rd.sz, 0, 0);
			rd.end = 0;
		}
		
		pthread_spin_lock(&w_spin);
		if (_compress_qualities)
			rd.sz += output_quality(qual, read, qmap + 0, rd.data + rd.sz, 0);
		rd.of = rd.sz;
		if (_use_second_file || _interleave) {
			rd.sz += output_read(read2, rd.data + rd.sz, 0, 0);
			if (_compress_qualities)
				rd.sz += output_quality(qual2, read2, qmap + 1, rd.data + rd.sz, 1);
		}

		bin_node *bn = aho_trie_bucket (bucket, &rd);
		total_size += rd.sz + sizeof(bin_node);
		file_reads++;
		pthread_spin_unlock(&w_spin);

		memcpy (bn->data.data, rd.data, rd.sz);

		if (total_size >= _max_bucket_set_size) {
			pthread_spin_lock(&w_spin);
			if (total_size >= _max_bucket_set_size) {
				dump_trie(temp_file_count++, trie);
				total_size = 0;
			}
			pthread_spin_unlock(&w_spin);
		}
	}
}

/************************/

void compress (char **files, int nf, const char *output, const char *pattern_path) {
	char line[MAXLINE], read[2][MAXLINE];


	LOG("Buffer size: %lldK, bucket storage size: %lldK\n", _file_buffer_size/1024, _max_bucket_set_size/1024);


	/* initialize the trie */
	DLOG ("Reading core strings ... ");
	trie = pattern_path[0] ? read_patterns_from_file(pattern_path) : read_patterns ();
	DLOG ("OK\n");


	/* initialize basic variables */
	LOG ("Preprocessing FASTQ files ...\n");
	int64_t original_size = 0;   /* size of input files, combined */

	char use_only_second_file = _use_second_file;        /* we need original _use... variable for simplification */
	if (_interleave)             /* all auxiliray functions will behave as in second file mode */
		_use_second_file = 1;

	/* init threading */
	threads = (pthread_t*) mallox(sizeof(pthread_t) * _thread_count);
	thread_index = (int*) mallox(_thread_count * sizeof(int));
	for (int i = 0; i < _thread_count; i++)
		thread_index[i] = i;
	pthread_spin_init(&r_spin, 0);
	pthread_spin_init(&w_spin, 0);

	/* initialize the files */
	f_init (input + 0, IO_GZIP);
	if (use_only_second_file)
		f_init (input + 1, IO_GZIP);

	/* obtain quality tables and stats */
	get_quality_stats (input, files[0], qmap);
	LOG("OK\n");

	LOG("Using %d threads...\n", _thread_count);

	/* iterate through files */
	for (int F = 0; F < nf; F++) {
		/* open files */
		for (int fi = 0; fi < use_only_second_file + 1; fi++) {
			char *fn = !fi ? files[F] : get_second_file (files[F]); 

			f_open (input + fi, fn, IO_READ); 
			if (!f_alive (input + fi)) 
				ERROR ("Cannot read file %s\n", fn);
			
			struct stat s;
			stat (input[fi].file_name, &s);
			original_size += s.st_size;
		}

		file_reads = 0;
		for (int ti = 0; ti < _thread_count; ti++)
			pthread_create(&threads[ti], 0, thread, (void*)&thread_index[ti]);
		for (int ti = 0; ti < _thread_count; ti++)
			pthread_join (threads[ti], 0);
		reads_count += file_reads;

		for (int fi = 0; fi < use_only_second_file + 1; fi++)
			f_close (input + fi);
		LOG("\tDone with file %s, %lld reads found\n", files[F], file_reads);
		if (use_only_second_file)
			LOG("\t          file %s, %lld reads found\n", get_second_file (files[F]), file_reads);

		file_reads = 0;
	}

	/* clean all stuff */
	if (total_size) {
		dump_trie (temp_file_count++, trie);
	}
	for (int fi = 0; fi < use_only_second_file + 1; fi++)
		f_free (input + fi);
	DLOG("\tDone!\n");

	/* merge */
	LOG("Merging results ... %d\n", temp_file_count);
	if (temp_file_count == 1)
		temp_file_count = 0;
	else
		merge(temp_file_count);

	/* compress and output */
	int64_t scalce_preprocessing_time = TIME;
	int64_t new_size = combine_and_compress_with_split (temp_file_count, output, (_split_reads ? _split_reads : reads_count), qmap[0].offset);

	/* cleaning again */
	LOG("Cleaning ...\n");
	aho_trie_free (trie);
	free(threads);
	free(thread_index);
	remove(_temp_directory);
	pthread_spin_destroy(&r_spin);
	pthread_spin_destroy(&w_spin);

	/* print useful stats */
	LOG("Statistics:\n");
	LOG("\tTotal number of reads: %lld\n", reads_count);
	LOG("\tRead length: first end %d\n", read_length[0]);
	if (_use_second_file) 
		LOG("\t             second end %d\n", read_length[1]);
	LOG("\tUnbucketed reads count: %d, bucketed percentage %.2lf\n", unbuck(), 100.0 * ((double)(reads_count - unbuck()) / reads_count));
	LOG("\tLossy percentage: %d\n", _quality_lossy_percentage);
	_time_elapsed = (TIME-_time_elapsed)/1000000;
	LOG("\tTime elapsed:     %02lld:%02lld:%02lld\n", _time_elapsed/3600, (_time_elapsed/60)%60, _time_elapsed%60);
	int64_t compress_time = (TIME - scalce_preprocessing_time) / 1000000;
	LOG("\tCompression time: %02lld:%02lld:%02lld\n", compress_time/3600, (compress_time/60)%60, compress_time%60);
	LOG("\tOriginal size: %.2lfM, new size: %.2lfM, compression factor: %.2lf\n", original_size/(1024.0*1024.4), new_size/(1024.0*1024.0), (original_size/(double)new_size));
	
}

