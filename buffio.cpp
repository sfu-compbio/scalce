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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include "zlib-1.2.5/zlib.h"
#include "bzip2-1.0.6/bzlib.h"
#include "const.h"
#include "buffio.h"

int64_t Zopen (const char *path, int mode, char *c) {
	gzFile handle;
	if (mode) {
		char bf[MAXLINE];
		snprintf (bf, MAXLINE, "%s", path); 
		handle = gzopen64 (bf, "wb");
		strncpy (c, bf, MAXLINE);
	}
	else {
		handle = gzopen64 (path, "rb");
		strncpy (c, path, MAXLINE);
	}
	return (int64_t)handle;
}

int     Zclose (void* handle)                      { return gzclose (handle); }
int64_t Zwrite (void* handle, void *c, int64_t sz) { return gzwrite (handle, c, sz); }
int64_t Zread  (void* handle, void *c, int64_t sz) { return gzread (handle, c, sz); }
int64_t Zseek  (void* handle, int64_t l)           { return gzseek64 (handle, l, SEEK_SET); }
char   *Zgets  (void* handle, char *c, int64_t l)  { return gzgets (handle, c, l); }

int64_t Bopen (const char *path, int mode, char *c) {
	BZFILE *handle;
	if (mode) {
		char bf[MAXLINE];
		snprintf (bf, MAXLINE, "%s", path); 
		handle = BZ2_bzopen (bf, "wb");
		strncpy (c, bf, MAXLINE);
	}
	else {
		handle = BZ2_bzopen (path, "rb");
		strncpy (c, path, MAXLINE);
	}
	return (int64_t)handle;
}

int     Bclose (void* handle)                      { BZ2_bzclose (handle); return 0; }
int64_t Bwrite (void* handle, void *c, int64_t sz) { return BZ2_bzwrite (handle, c, sz); }	
int64_t Bread  (void* handle, void *c, int64_t sz) { return BZ2_bzread (handle, c, sz); }
int64_t Bseek  (void* handle, int64_t l)           { ERROR("Seek not supported with bzlib"); }
char   *Bgets  (void* handle, char *c, int64_t l)  { ERROR("Gets not supported with bzlib"); }

int64_t Sopen (const char *path, int mode, char *c) {
	int64_t handle;

	//int handle;
	if (!strcmp (path, "-")) handle = 1; // (int)stdout;
	else handle = (int64_t)fopen(path, mode ? "wb" : "rb"
/*	else handle = open64 (path, //O_DIRECT | O_SYNC |
			(mode ? (O_CREAT | O_WRONLY | O_TRUNC) 
			 : (O_RDONLY)), 
			S_IRUSR | S_IWUSR*/
			);
	strncpy (c, path, MAXLINE);
	if (handle <= 0) handle = 0;
	return handle;
}

int     Sclose (void* handle)                      { return fclose ((FILE*)handle); }
int64_t Swrite (void* handle, void *c, int64_t sz) { return fwrite (c, 1, sz, (FILE*)handle); }	
int64_t Sread  (void* handle, void *c, int64_t sz) { return fread (c, 1, sz, (FILE*)handle); }
int64_t Sseek  (void* handle, int64_t l)           { return fseek ((FILE*)handle, l, SEEK_SET); }
char   *Sgets  (void* handle, char *c, int64_t l)  { return fgets (c, l, (FILE*)handle); }
	

int64_t PZopen (const char *path, int mode, char *c) {
	if (!mode) 
		ERROR("Parallel support available only in write-mode");

	char bf[MAXLINE];
	snprintf (bf, MAXLINE, "%s", path); 
	strncpy (c, bf, MAXLINE);

	int fd[2];
	if (pipe (fd) == -1)
		ERROR ("Pipe failed");
	int pid = fork ();
	if (pid == -1)
		ERROR ("Fork failed");

	if (pid == 0) { // child, call pigz
		close (fd[1]);

		int fo = open64 (bf, O_CREAT | O_WRONLY | O_TRUNC, S_IRUSR | S_IWUSR);
		dup2 (fo, 1); // set stdout to file
		dup2 (fd[0], 0); // set stdin to pipe
		close (fd[0]);

		char* pl[5]; char tmp[20];
		pl[0] = (char*)"pigz";
		pl[1] = (char*)"-c";
		pl[2] = tmp;
		sprintf(pl[2], "-p%d", _thread_count);
		pl[3] = (mode ? 0 : (char*)"-d");
		pl[4] = 0;
		execvp ("pigz", pl);
		ERROR ("pigz exec failed");
		abort();
	}
	else { // parent
		close (fd[0]);
		return fd[1];
	}

	return 0;
}

int     PZclose (void* handle)                      { return close ((int64_t)handle); }
int64_t PZwrite (void* handle, void *c, int64_t sz) { return write ((int64_t)handle, c, sz); }	
int64_t PZread  (void* handle, void *c, int64_t sz) { return read ((int64_t)handle, c, sz); }
int64_t PZseek  (void* handle, int64_t l)           { ERROR("Seek not ssupported with pigz"); }
char   *PZgets  (void* handle, char *c, int64_t l)  { ERROR("Gets not supported with pigz"); }

char f_alive (buffered_file *f) {
	return (f->handle != 0);
}

int64_t f_write (buffered_file *f, void *c, int64_t sz) {
//	if (f->mode == 0 || !f_alive (f))
//		ERROR ("Invalid file handle (closed or opened in read mode)");

	return f->my_write((void*)f->handle, c, sz);

/*	int64_t wp = 0;
	while (sz >= f->bufsz - f->bufpos) {
		memcpy (f->buffer + f->bufpos, (char*)c + wp, f->bufsz - f->bufpos);
		sz -= (f->bufsz - f->bufpos);

		int64_t n = (f->my_write) ((void*)(f->handle), f->buffer, f->bufsz);
		wp += (f->bufsz - f->bufpos);
		if (n != f->bufsz)
			return wp;
		f->bufpos = 0;
	}
	memcpy (f->buffer + f->bufpos, (char*)c + wp, sz);
	f->bufpos += sz;
	return wp + sz;*/
}

int64_t f_read (buffered_file *f, void *c, int64_t sz) {
//	if (f->mode == 1 || !f_alive (f))
//		ERROR ("Invalid file handle (closed or opened in write mode)");

	return f->my_read((void*)f->handle, c, sz);

	/*

	int64_t rp = 0;
	while (sz - rp >= f->bufend - f->bufpos) {
		memcpy ((char*)c + rp, f->buffer + f->bufpos, f->bufend - f->bufpos);
		rp += f->bufend - f->bufpos;

		f->bufend = (f->my_read) ((void*)(f->handle), f->buffer, f->bufsz);
		f->bufpos = 0;
		if (f->bufend == 0) 
			return rp;
	}
	memcpy ((char*)c + rp, f->buffer + f->bufpos, sz - rp);
	f->bufpos += sz - rp;
	return sz;*/
}

void f_seek (buffered_file *f,  int64_t pos) {
	f->my_seek((void*)f->handle, pos);
/*
	if (!f_alive (f)) 
		ERROR("Invalid file handle (probably closed)");
	if (f->mode == 1 && f->bufpos > 0)
		(f->my_write) ((void*)(f->handle), f->buffer, f->bufpos);
	int64_t x;
	if ((x = (f->my_seek) ((void*)(f->handle), pos)) != pos) {
		ERROR("Seek failed! %lld, got %lld", pos, x);
		abort();
	}
	f->bufpos = 0;
	f->bufend = f->bufsz;
	if (!f->mode)
		f->bufend = (f->my_read) ((void*)(f->handle), f->buffer, f->bufsz);*/
}

void f_open (buffered_file *f, const char *c, int m) {
	if (f_alive (f))
		f_close (f);
	f->mode = m;
	f->handle = (f->my_open) (c, m, f->file_name);

	/*f->bufpos = 0;
	f->bufend = f->bufsz;

	if (m == 0 && f_alive (f)) {
		f->bufend = (f->my_read) ((void*)(f->handle), f->buffer, f->bufsz);
	}*/
}

void f_close (buffered_file *f) {
//	if (f->mode == 1 && f_alive (f) && f->bufpos > 0)
//		(f->my_write) ((void*)(f->handle), f->buffer, f->bufpos);
	if (f_alive (f))
		(f->my_close) ((void*)(f->handle));
	f->handle = 0;
}

char* f_gets (buffered_file *f, char *c, int64_t maxsz) {
	return f->my_gets((void*)f->handle, c, maxsz);
	/*
	int64_t i, n;
	for (i = 0; i < maxsz - 1; i++) {
		n = f_read (f, c + i, 1);
		if (n == 0 || c[i] == '\n')
			break;
	}
	if (i) c[i] = 0;
	return (i ? c : 0);*/
}


void f_init (buffered_file *f, int algorithm) { 
//	f->bufsz = _file_buffer_size;
//	f->buffer = (char*) mallox (_file_buffer_size);
	f->handle = 0;
	f->mode = 0; //f->bufpos = f->bufend = 0;

	switch (algorithm) {
		case IO_SYS:
			f->my_open  = Sopen;
			f->my_close = Sclose;
			f->my_seek  = Sseek;
			f->my_read  = Sread;
			f->my_write = Swrite;
			f->my_gets  = Sgets;
			break;
		case IO_GZIP:
			f->my_open  = Zopen;
			f->my_close = Zclose;
			f->my_seek  = Zseek;
			f->my_read  = Zread;
			f->my_write = Zwrite;
			f->my_gets  = Zgets;
			break;
		case IO_PGZIP:
			f->my_open  = PZopen;
			f->my_close = PZclose;
			f->my_seek  = PZseek;
			f->my_read  = PZread;
			f->my_write = PZwrite;
			f->my_gets  = PZgets;
			break;
		case IO_BZIP:
			f->my_open  = Bopen;
			f->my_close = Bclose;
			f->my_seek  = Bseek;
			f->my_read  = Bread;
			f->my_write = Bwrite;
			f->my_gets = Bgets;
			break;
		default:
			ERROR("Unknown compression mode");
	}
}

void f_free (buffered_file *f) {
	f_close (f);
//	frex (f->buffer, _file_buffer_size);
}

