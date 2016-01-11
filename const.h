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

#ifndef CONST_H__
#define CONST_H__

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <sys/types.h>
#include <sys/sysinfo.h>
#include <unistd.h>
#include <sys/stat.h>

void MEM(char *A);

int64_t _TIME_();

#define GLOBALBUFSZ  5*MAXLINE
#define MIN(x,y)  (((x)<(y))?(x):(y))
#define MAX(x,y)  (((x)>(y))?(x):(y))

#define TIME	(_TIME_())
#define SZ_READ(l)\
	(((l) / 4) + ((l) % 4 > 0))
#define SZ_QUAL(l)\
	( l )
//3 * ( (l / 4) + (l % 4 > 0) ) )
#ifdef DEBUG
#define LOG(c,...)\
	fprintf(stderr, "[%02d:%02d:%02d]: "c, TIME/3600, (TIME/60)%60, TIME%60, ##__VA_ARGS__)
#define DLOG(c,...)\
	fprintf(stderr, "[%02d:%02d:%02d]: "c, TIME/3600, (TIME/60)%60, TIME%60, ##__VA_ARGS__)
#else 
#define LOG(c,...) fprintf(stderr,c,##__VA_ARGS__)
#define DLOG(c,...) //fprintf(stderr,c,##__VA_ARGS__)
#endif
#define ERROR(c,...)\
	{ fprintf (stderr, "(ERROR) "c, ##__VA_ARGS__); exit (1); }


// Max. characters per line
#ifdef PACBIO
	#define MAXLINE 	100000
#else
	#define MAXLINE 	2500
#endif
	
// Max. files in file pool
#define MAXOPENFILES 20
// Max. concurrent files for merging
#define MAXMERGE     (MAXOPENFILES / 2 - 1)
// Max. bin ID
#define MAXBIN       (1 << 30)
// Max. metadata file size, 200MB, will 
#define MAXMETA      (300 * 1024 * 1024)

// Global program arguments
extern int       _quality_sample_lines;
extern int       _quality_lossy_percentage;
extern char      _use_second_file;
extern char      _use_names;
extern uint64_t  _file_buffer_size;
extern uint64_t  _max_bucket_set_size;
extern char      _temp_directory[MAXLINE];
extern char      _library_name[MAXLINE];
extern char      _pattern_path[MAXLINE];  
extern int64_t       _time_elapsed;
extern int       _split_reads;
extern int       _compression_mode;
extern int       _thread_count;
extern char      _interleave;
extern int       _decompress;
extern char      _is_fasta;
extern int       _no_ac ;
extern int     _compress_qualities;
extern int32_t read_length[2];                 /* length of read in one fastq file */
extern int64_t reads_count;                /* total read count in all input files */

char *get_second_file (const char *c);
void *mallox (size_t size);
void frex (void *ptr, size_t size);
double getmemx ();

extern int _tbl[];// ={ 0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0 };
#define getval(c) (_tbl[c-'A'])
#endif // CONST_H__
