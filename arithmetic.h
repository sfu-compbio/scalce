/// 786

#ifndef ARITHMETIC_H__
#define ARITHMETIC_H__

#include "buffio.h"

extern int Z;
extern
uint32_t	ac_freq3[2][46*46],
		   ac_freq4[2][46*46][46],
			ac_sz;
void a_init_coder (buffered_file *f);
void a_init_decoder (buffered_file *f);
void a_read (buffered_file *f, uint8_t *arr, int sz);
void a_write (buffered_file *f, uint8_t *arr, int sz);
void a_finish_coder (buffered_file *f);

#endif // ARITHMETIC_H__
