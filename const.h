/// 786

#ifndef CONST_H__
#define CONST_H__

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <sys/types.h>
#include <sys/sysinfo.h>

void MEM(char *A);

int64_t _TIME_();

#define GLOBALBUFSZ  5*MAXLINE
#define MIN(x,y)  (((x)<(y))?(x):(y))
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
#define MAXLINE 		2500
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

extern int       _no_ac ;
extern int32_t read_length[2];                 /* length of read in one fastq file */
extern int64_t reads_count;                /* total read count in all input files */

char *get_second_file (const char *c);
void *mallox (size_t size);
void frex (void *ptr, size_t size);
double getmemx ();

extern int _tbl[];// ={ 0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0 };
#define getval(c) (_tbl[c-'A'])
#endif // CONST_H__
