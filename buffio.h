/// 786

#ifndef BUFFIO_H__
#define BUFFIO_H__

#include "const.h"

typedef int64_t (*fn_rw)    (void*, void*, int64_t);
typedef int     (*fn_close) (void*);
typedef int64_t (*fn_open)  (const char*, int, char*);  
typedef int64_t (*fn_seek)  (void*, int64_t);

typedef struct {
	char file_name[MAXLINE];
	char *buffer;
	int mode;
	int64_t bufsz, bufpos, bufend;
	int64_t handle;

	// Function pointers
	fn_rw    my_read,
				my_write;
	fn_open  my_open;
	fn_close my_close;
	fn_seek  my_seek;
} buffered_file;

void     f_init  (buffered_file *f, int algorithm);
void     f_free  (buffered_file *f);
char     f_alive (buffered_file *f);
int64_t  f_write (buffered_file *f, void *c, int64_t sz);	
int64_t  f_read  (buffered_file *f, void *c, int64_t sz);	
void     f_seek  (buffered_file *f, int64_t pos);
void     f_open  (buffered_file *f, const char *c, int m);
void     f_close (buffered_file *f);
char*    f_gets  (buffered_file *f, char *c, int64_t maxsz);

#endif // BUFFIO_H
