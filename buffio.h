/*
 * Copyright (c) 2011 - 2012, Simon Fraser University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list
 * of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice,
 * this
 *   list of conditions and the following disclaimer in the documentation and/or
 * other
 *   materials provided with the distribution.
 * - Neither the name of the Simon Fraser University nor the names of its
 * contributors may be
 *   used to endorse or promote products derived from this software without
 * specific
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

#ifndef BUFFIO_H__
#define BUFFIO_H__

#include "const.h"

#define IO_READ 0
#define IO_WRITE 1
#define IO_SYS 10
#define IO_GZIP 20
#define IO_BZIP 30
#define IO_PGZIP 40

typedef int64_t (*fn_rw)(void *, void *, int64_t);
typedef int (*fn_close)(void *);
typedef int64_t (*fn_open)(const char *, int, char *);
typedef int64_t (*fn_seek)(void *, int64_t);
typedef char *(*fn_gets)(void *, char *, int64_t);

typedef struct {
  char file_name[MAXLINE];
  char *buffer;
  int mode;
  int64_t bufsz, bufpos, bufend;
  int64_t handle;

  // Function pointers
  fn_rw my_read, my_write;
  fn_open my_open;
  fn_close my_close;
  fn_seek my_seek;
  fn_gets my_gets;
} buffered_file;

void f_init(buffered_file *f, int algorithm);
void f_free(buffered_file *f);
char f_alive(buffered_file *f);
int64_t f_write(buffered_file *f, void *c, int64_t sz);
int64_t f_read(buffered_file *f, void *c, int64_t sz);
void f_seek(buffered_file *f, int64_t pos);
void f_open(buffered_file *f, const char *c, int m);
void f_close(buffered_file *f);
char *f_gets(buffered_file *f, char *c, int64_t maxsz);

#endif // BUFFIO_H
