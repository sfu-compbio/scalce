/// 786

#include <time.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "const.h"
#include "buffio.h"
#include "compress.h"
#include "reads.h"
#include "qualities.h"
#include "names.h"

buffered_file file_pool[MAXOPENFILES];  /* handles for files */
int32_t read_length[2];                 /* length of read in one fastq file */
int64_t reads_count = 0;                /* total read count in all input files */

char global_buffer[GLOBALBUFSZ];        /* global shared buffer for reading/writing */
uint8_t metadata[MAXMETA];              /* metadata file contents */

int merhamet_merge (int *fl, int flc, int idx) {
	const int MAXSEE = MAXOPENFILES / 2 - 1; // max

	char buffer[MAXLINE];

	LOG ("Merging part %d for files + ", idx);
	for (int i = 0; i < flc; i++) {
		struct stat S;
		snprintf (buffer, MAXLINE, "%s/t_%03d_%d.tmp", _temp_directory, fl[i], idx);
		stat (buffer, &S);
		LOG ("%02d (%7.1lfM)%s", fl[i], (S.st_size/1024.0)/1024.0, (i==flc-1?"":", "));
	}
	LOG (" ... ");

	snprintf (buffer, MAXLINE, "%s/t_TMP_%d.tmp", _temp_directory, idx);
	f_open (file_pool + 0, buffer, IO_WRITE);
	for (int i = 0; i < flc; i++) {
		snprintf(buffer, MAXLINE, "%s/t_%03d_%d.tmp", _temp_directory, fl[i], 3); // meta
		f_open (file_pool + i + 1, buffer, IO_READ);
		snprintf(buffer, MAXLINE, "%s/t_%03d_%d.tmp", _temp_directory, fl[i], idx); // other
		f_open (file_pool + flc + i + 1, buffer, IO_READ);
	}

	int32_t bins[MAXMERGE];
	int minV = MAXBIN;
	int minI = 0;
	for (int i = 0; i < flc; i++) { 
		f_read (file_pool + i + 1, &(bins[i]), sizeof(int32_t));
		if (bins[i] < minV) {
			minV = bins[i];
			minI = i;
		}
	}

	int     binex = flc, metadata_pos = 0;
	int64_t totalLen = 0;
	int32_t prevMinV;
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
		minV = MAXBIN;
		int rd = f_read (file_pool + minI + 1, &(bins[minI]), sizeof(int32_t));
		if (rd != sizeof(int32_t)) {
			bins[minI] = MAXBIN;
			binex--;
		}
		for (int i = 0; i < flc; i++) 
			if (bins[i] < minV) {
				minV = bins[i];
				minI = i;
			}
		if (minV != prevMinV) {
			memcpy (metadata + metadata_pos, &prevMinV, sizeof(int32_t));
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

	LOG("OK\n");
	return metadata_pos;
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
	f_init (&foQ, _compression_mode);
	f_init (&foR, _compression_mode);
	f_init (&foN, _compression_mode);
	LOG ("compression: %s\n", _compression_mode==IO_SYS?"no compression":(_compression_mode==IO_GZIP?"gzip":"bzip"));

	int64_t new_size = 0;

	for (int NP = 0; NP < _use_second_file + 1; NP++) {
		int64_t lR = 0, lQ = 0, lN = 0;
		char alive = 1;
		for (int64_t F = 0; alive; F++) {
			int64_t limit = MIN (split_size, reads_count - F * split_size);

			snprintf (buffer, MAXLINE, "%s.%02d_%d.scalcen", output, F, NP+1);
			f_open (&foN, buffer, IO_WRITE);
			snprintf (buffer, MAXLINE, "%s.%02d_%d.scalceq", output, F, NP+1);
			f_open (&foQ, buffer, IO_WRITE);
			snprintf (buffer, MAXLINE, "%s.%02d_%d.scalcer", output, F, NP+1);
			f_open (&foR, buffer, IO_WRITE);

			f_write (&foR, read_length + NP, sizeof(int32_t));
			f_write (&foR, &limit, sizeof(int64_t));
			f_write (&foQ, &phred_offset, sizeof(int64_t));
			f_write (&foN, &_use_names, 1);
			if (!_use_names) {
				int64_t name_idx = F * split_size;
				f_write (&foN, &name_idx, sizeof(int64_t));
				f_write (&foN, _library_name, strlen(_library_name));
			}

			int64_t reads_written = 0, temp;
			while (reads_written < split_size) {
				int64_t num_reads = lR / (SZ_READ (read_length[NP]));
				if (reads_written + num_reads <= split_size) {
					// write complete bucket
					for (int64_t i = 0; i < lR; i += GLOBALBUFSZ) {
						int64_t r = f_read (fR, global_buffer, MIN (lR - i, GLOBALBUFSZ));
						f_write (&foR, global_buffer, r);
					}
					for (int64_t i = 0; i < lQ; i += GLOBALBUFSZ) {
						int64_t r = f_read (fQ, global_buffer, MIN (lQ - i, GLOBALBUFSZ));
						f_write (&foQ, global_buffer, r);
					}
					if (_use_names) for (int64_t i = 0; i < lN; i += GLOBALBUFSZ) {
						int64_t r = f_read (fN, global_buffer, MIN (lN - i, GLOBALBUFSZ));
						f_write (&foN, global_buffer, r);
					}
					reads_written += num_reads;		
					if (f_read (fM, &lN, sizeof(int32_t)) < sizeof(int32_t)) {
						alive = 0;
						break;
					}
				
					f_read (fM, &lN, sizeof(int64_t));
					f_read (fM, &lR, sizeof(int64_t));
					f_read (fM, &lQ, sizeof(int64_t));
					if (_use_second_file) {
						int64_t *lR2 = (NP ? &lR : &temp);
						int64_t *lQ2 = (NP ? &lQ : &temp);
						f_read (fM, lR2, sizeof(int64_t));
						f_read (fM, lQ2, sizeof(int64_t));
					} 
				}
				else {
					int64_t left = split_size - reads_written;
					for (int64_t i = 0; i < left; i++) {
						int64_t r = f_read (fR, buffer, SZ_READ(read_length[NP]));
						f_write (&foR, buffer, r);
					}
					lR -= left * SZ_READ(read_length[NP]);
					for (int64_t i = 0; i < left; i++) { 
						int pos = 0, len = 0;
						while (len < read_length[NP]) {
							f_read (fQ, buffer + pos, 1);
							uint8_t t = buffer[pos];
							len += (t >= 128 ? t - 128 : 1);
							pos++;
						}
						f_write (&foQ, buffer, pos);
						lQ -= pos;
					}
					if (_use_names) for (int64_t i = 0; i < left; i++) {
						f_read (fN, buffer, 1);
						f_read (fN, buffer + 1, buffer[0]);
						f_write (&foN, buffer, buffer[0] + 1);
						lN -= (buffer[0] + 1);
					}
					reads_written += left;
				}
			}

			f_close (&foR);
			f_close (&foN);
			f_close (&foQ);
			LOG("\tCreated part %d (paired-end %d) with %lld reads, files (%s, %s, %s)\n", F, NP, reads_written, foN.file_name, foQ.file_name, foR.file_name);

			// ...
			struct stat s;
			stat (foR.file_name, &s); new_size += s.st_size;
			stat (foQ.file_name, &s); new_size += s.st_size;
			stat (foN.file_name, &s); new_size += s.st_size;

		}
		if (_use_second_file && !NP) {
			fR = file_pool + 4;
			fQ = file_pool + 5;
			f_seek (fN, 0);
			f_seek (fM, 0);
		}
	}

	LOG("\tCleaning ...\n");
	f_free (&foQ);
	f_free (&foN);
	f_free (&foR);
	for (int i = 0; i <= 3 + 2 * (_use_second_file); i++) { 
		snprintf(buffer, MAXLINE, "%s/t_%03d_%d.tmp", _temp_directory, fl, i);
		f_close (file_pool + i);
		unlink (buffer);
	}

	LOG("\tDone!\n");
	return new_size;
}

void merge (int total, int param) {
	char buffer[MAXLINE], buffer2[MAXLINE];

	char *used = calloc (total + 5, sizeof (char));

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

	free (used);
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

void compress (char **files, int nf, const char *output, const char *pattern_path) {
	char line[MAXLINE], read[2][MAXLINE];
	quality_mapping qmap[2];

	if (_interleave) _use_second_file = 1;

	LOG("Buffer size: %lldK, bucket storage size: %lldK\n", _file_buffer_size/1024, _max_bucket_set_size/1024);

	for (int i = 0; i < MAXOPENFILES; i++)
		f_init (file_pool + i, IO_SYS);

	LOG("Reading core strings ... ");
	aho_trie *trie = read_patterns ();
	LOG("OK\n");

	LOG("Preprocessing FASTQ files ...\n");
	uint64_t total_size = 0;
	int temp_file_count = 0;

	buffered_file f[2];
	f_init (f + 0, IO_GZIP);
	if (_use_second_file && !_interleave)
		f_init (f + 1, IO_GZIP);
	
	int64_t original_size = 0;
	for (int i = 0; i < nf; i++) {
		for (int F = 0; F < (_use_second_file&&!_interleave) + 1; F++) {
			if (F) {
				if (get_second_file (files[i]))
					f_open (f + F, get_second_file (files[i]), IO_READ);
				else
					ERROR ("Cannot find paired-end for file %s\n", files[i]);
			}
			else
				f_open (f + F, files[i], IO_READ);
			if (!f_alive(f + F))
				ERROR ("File %s does not exist\n", f[F].file_name);
			
			DLOG("\tAnalyzing file %s\n", files[i]);
			struct stat s;
			stat (f[F].file_name, &s);
			original_size += s.st_size;

			/* Get quality information from first file only */
			if (i == 0) {
				if(_interleave) _interleave=10;
				quality_mapping_init (qmap + F, f + F, read_length + F);
				if(_interleave) _interleave=1;
				f_seek (f + F, 0);
				LOG("\tPaired end #%d, quality offset: %d\n"
					 "\t               read length: %d\n", F+1,qmap[F].offset, read_length[F]);
			}
		}
		if (_interleave&&i==0) {
			_interleave=20;
			quality_mapping_init (qmap + 1, f, read_length + 1);
			_interleave=1;
			f_seek (f, 0);
			LOG("\tPaired end #%d, quality offset: %d\n"
				 "\t               read length: %d\n", 2,qmap[1].offset, read_length[1]);
		}

		
		int64_t file_reads = 0;
		read_data lo;
		lo.data = global_buffer;
		while (f_gets (f + 0, line, MAXLINE)) {
			int64_t sz = output_name (line, (uint8_t*)global_buffer);
			f_gets(f + 0, read[0], MAXLINE);
			sz += output_read (read[0], (uint8_t*)(global_buffer + sz));
			f_gets (f + 0, line, MAXLINE);
			f_gets (f + 0, line, MAXLINE);
			sz += output_quality (line, read[0], qmap + 0, (uint8_t*)(global_buffer + sz));
			int64_t of2 = sz;
			if (_use_second_file&&!_interleave) {
				f_gets(f + 1, line, MAXLINE);
				f_gets(f + 1, read[1], MAXLINE);
				sz += output_read (read[1], (uint8_t*)(global_buffer + sz));
				f_gets(f + 1, line, MAXLINE);
				f_gets(f + 1, line, MAXLINE);
				sz += output_quality (line, read[1], qmap + 1, (uint8_t*)(global_buffer + sz));
			}
			if (_interleave) {
				f_gets(f, line, MAXLINE);
				f_gets(f, read[1], MAXLINE);
				sz += output_read (read[1], (uint8_t*)(global_buffer + sz));
				f_gets(f, line, MAXLINE);
				f_gets(f, line, MAXLINE);
				sz += output_quality (line, read[1], qmap + 1, (uint8_t*)(global_buffer + sz));
			}
		
			lo.of = of2;
			lo.sz = sz;
			aho_search (read[0], trie, &lo);

			total_size += sz + 2 * sizeof(int64_t) + sizeof(struct bin_node*);
			reads_count++;
			file_reads++;
			if (total_size >= _max_bucket_set_size) {
				total_size = 0;
				dump_trie (temp_file_count++, trie);
			}
		}
		f_close (f + 0);
		if (_use_second_file&&!_interleave)
			f_close (f + 1);
		LOG("\tDone with file %s with %lld reads\n", files[i], file_reads);
		if (_use_second_file&&!_interleave) 
			LOG("\t     also file %s with %lld reads\n", get_second_file (files[i]), file_reads); 
	}
	if (total_size)
		dump_trie (temp_file_count++, trie);

	f_free (f + 0);
	if (_use_second_file&&!_interleave)
		f_free (f + 1);
	LOG("\tDone!\n");

	LOG("Merging results\n");
	merge (temp_file_count, 8);

	int tm = TIME;

	int64_t new_size = combine_and_compress_with_split (temp_file_count-1, output, (_split_reads ? _split_reads : reads_count), qmap[0].offset);

	LOG("Cleaning\n");
	aho_trie_free (trie);
	remove(_temp_directory);
	for (int i = 0; i < MAXOPENFILES; i++)
		f_free (file_pool + i);

	LOG("Statistics\n");
	LOG("\tTotal number of reads: %lld\n", reads_count);
	LOG("\tRead length: first end %d\n", read_length[0]);
	if (_use_second_file) 
		LOG("\t             second end %d\n", read_length[1]);
	LOG("\tUnbucketed reads count: %d, bucketed percentage %.2lf\n", unbuck(), 100.0 * ((double)(reads_count - unbuck()) / reads_count));
	LOG("\tLossy percentage: %d\n", _quality_lossy_percentage);
	LOG("\tTime elapsed: %02d:%02d:%02d\n", TIME/3600, (TIME/60)%60, TIME%60);
	LOG("\tCompression time: %02d:%02d:%02d\n", (TIME-tm)/3600, ((TIME-tm)/60)%60, (TIME-tm)%60);
	LOG("\tOriginal size: %.2lfM, new size: %.2lfM, compression factor: %.2lf\n", original_size/(1024.0*1024.4), new_size/(1024.0*1024.0), (original_size/(double)new_size));
}

