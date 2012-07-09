/// 786

#include <assert.h>
#include <inttypes.h>
#include <string.h>
#include <pthread.h>
#include "arithmetic.h"

const int buffer_size = 10 * 1024 * 1024;
uint64_t ac_freq3[2][AC_DEPTH*AC_DEPTH],
		   ac_freq4[2][AC_DEPTH*AC_DEPTH*AC_DEPTH];

/*******************************************/

ac_stat::ac_stat (uint64_t *a3, uint64_t *a4) {
	acfq3 = a3;
	acfq4 = a4;

	for (int i = 0; i < AC_DEPTH; i++) {
		for (int j = 0; j < AC_DEPTH; j++) {
			hi_fq[i*AC_DEPTH + j][0] = acfq4[(i*AC_DEPTH + j)*AC_DEPTH+0];
			for (int l = 1; l < AC_DEPTH; l++) 
				hi_fq[i*AC_DEPTH + j][l] = hi_fq[i*AC_DEPTH + j][l - 1] + acfq4[(i*AC_DEPTH + j)*AC_DEPTH+l];
			acfq3[i*AC_DEPTH + j] = hi_fq[i*AC_DEPTH + j][AC_DEPTH - 1];
			assert(acfq3[i*AC_DEPTH+j]<((1LL<<32)-1));
			int p = -1;
			for (int l = 0; l < AC_DEPTH; l++) if (acfq4[(i*AC_DEPTH + j)*AC_DEPTH+l]) {
				if (p != -1) 
					lo_fq[i*AC_DEPTH + j][l] = hi_fq[i*AC_DEPTH + j][p];
				p = l;
			}
		}
	}
}

/****************************************/

ac_coder::ac_coder (uint8_t *ou, ac_stat *as)
	: ou(ou), as(as), os(ou) 
{
	reset();
}
void ac_coder::reset () {
	lo = 0;
	hi = -1;
	underflow = 0;
	prev[0] = prev[1] = prev[2] = 500; // molim?!	
	ou = os;
	buf_pos = 7;
}

void ac_coder::O (int t) {
	*ou <<= 1;
	if (t) *ou |= 1;
	if (buf_pos == 0) {
		ou++;
		assert(ou-os < buffer_size);
		buf_pos = 7;
		*ou = 0;
	}
	else --buf_pos;
}

void ac_coder::write (uint8_t *arr, int sz) {
	int i = 0;
	if (prev[0] > 255) { 
		*ou = arr[0]; ou++;
		*ou = arr[1]; ou++;
		prev[0] = arr[0];
		prev[1] = arr[1];
		i = 2;
	}
	for (; i < sz; i++) {
		uint32_t c_lo = as->lo_fq[prev[0]*AC_DEPTH + prev[1]][arr[i]],
					c_hi = as->hi_fq[prev[0]*AC_DEPTH + prev[1]][arr[i]],
					p_sz = as->acfq3[prev[0]*AC_DEPTH + prev[1]];
		uint64_t range = hi - lo + 1LL;
		hi = lo + (range * c_hi) / p_sz - 1;
		lo = lo + (range * c_lo) / p_sz;

		while (1) {
			if ((hi & 0x80000000) == (lo & 0x80000000)) {
				O (hi & 0x80000000);
				while (underflow) {
					O (~hi & 0x80000000); 
					underflow--;
				}
			}
			else if (!(hi & 0x40000000) && (lo & 0x40000000)) { // underflow ante portas
				underflow++;
				lo &= 0x3fffffff; // set 2nd highest digit to 0 -- 1st will be removed after
				hi |= 0x40000000; // set 2nd highest digit to 1
			}
			else break;

			lo <<= 1;
			hi <<= 1; hi |= 1;
		}

		prev[0] = prev[1];
		prev[1] = arr[i];
	}
}

void ac_coder::flush () {
	O (lo & 0x40000000);
	underflow++;
	while (underflow) {
		O (~lo & 0x40000000);
		underflow--;
	}
	while (buf_pos != 7)
		O (0);
}

/***************************/

ac_decoder::ac_decoder(ac_stat *as, uint8_t *in)
	: as(as), in(in), is(in)
{
	reset();
}

uint8_t ac_decoder::I() {
	uint8_t c = ((*in) >> buf_pos) & 1;
	if (!buf_pos) {
		in++;
		assert(in-is < buffer_size);
		buf_pos = 7;
	}
	else buf_pos--;
	return (c ? 1 : 0);
}

void ac_decoder::reset() {
	lo = 0;
	hi = -1;
	prev[0] = prev[1] = prev[2] = 500;
	in = is;
}

uint8_t ac_decoder::read_single() {
	uint8_t c = 0;
	uint64_t range = hi - lo + 1LL;
	uint32_t count = ((code - lo + 1LL) * as->acfq3[prev[0]*AC_DEPTH + prev[1]] - 1LL) / range;

	uint32_t i;
	for (i = 0; i < AC_DEPTH; i++) 
		if (as->acfq4[(prev[0]*AC_DEPTH + prev[1]) * AC_DEPTH + i] 
				&& count >= as->lo_fq[prev[0]*AC_DEPTH + prev[1]][i] 
				&& count <  as->hi_fq[prev[0]*AC_DEPTH + prev[1]][i]) 
			break;
	assert(i<AC_DEPTH);

	hi = lo + (range * as->hi_fq[prev[0]*AC_DEPTH + prev[1]][i]) / as->acfq3[prev[0]*AC_DEPTH + prev[1]] - 1;
	lo = lo + (range * as->lo_fq[prev[0]*AC_DEPTH + prev[1]][i]) / as->acfq3[prev[0]*AC_DEPTH + prev[1]];
	while (1) {
		if ((hi & 0x80000000) == (lo & 0x80000000)) ;
		else if (!(hi & 0x40000000) && (lo & 0x40000000)) { // underflow ante portas
			code ^= 0x40000000;
			lo   &= 0x3fffffff; 
			hi   |= 0x40000000; 
		}
		else break;
		lo <<= 1;
		hi <<= 1; hi |= 1;
		code <<= 1; code |= I();
	}

	prev[0] = prev[1];
	prev[1] = i;
	return i;
}

void ac_decoder::read(uint8_t *ar, int sz) {
	int i = 0;
	if (prev[0] > 255) {
		ar[0] = *in; in++;
		ar[1] = *in; in++;
		prev[0] = ar[0];
		prev[1] = ar[1];
		i = 2;

		// f_read(f, &buffer, sizeof(uint8_t));
		buf_pos = 7;
		code = 0;
		for (int j = 0; j < 32; j++) {
			code <<= 1;
			code |= I();
		}
	} 

	for (; i < sz; i++)
		ar[i] = read_single();
}

/****************************/

static pthread_t *threads;
static uint8_t  *input, *output;
static uint32_t *input_size, *output_size;
static int *thread_index;
ac_stat as;

void set_ac_stat(uint64_t *a3, uint64_t *a4) {
	as = ac_stat(a3,a4);
}

void *thread_c (void *tv) {
	int tx = *((int*)tv);
	ac_coder ax(output + tx * buffer_size, &as);
	ax.write(input + tx * buffer_size, input_size[tx]);
	ax.flush();
	output_size[tx] = ax.output() - (output + tx * buffer_size);
	return 0;
}

void *thread_d (void *tv) {
	int tx = *((int*)tv);
	ac_decoder ad(&as, input + tx * buffer_size);
	ad.read(output + tx * buffer_size, output_size[tx]);
	return 0;
}

/****************************/

void ac_init () {
	threads = (pthread_t*)malloc(sizeof(pthread_t) * _thread_count);
	input = (uint8_t*)malloc(buffer_size * _thread_count);
	output = (uint8_t*)malloc(buffer_size * _thread_count);
	output_size = (uint32_t*)malloc(sizeof(uint32_t) * _thread_count);
	input_size = (uint32_t*)malloc(sizeof(uint32_t) * _thread_count);
	thread_index = (int*)malloc(_thread_count * sizeof(int));
	for (int i = 0; i < _thread_count; i++)
		thread_index[i] = i;
}

void ac_finish() {
	free(threads);
	free(input);
	free(output);
	free(output_size);
	free(input_size);
	free(thread_index);
}

void ac_write (buffered_file *fo, uint8_t *data, int size) {
	static int current_position = 0;
	if (!size) { // flush
		int l_th = current_position / buffer_size;
		for (int i = 0; i < l_th; i++) {
			input_size[i] = buffer_size;
			pthread_create(&threads[i], 0, thread_c, (void*) &thread_index[i]);
		}
		if (current_position % buffer_size != 0) {
			input_size[l_th] = current_position % buffer_size;
			pthread_create(&threads[l_th], 0, thread_c, (void*) &thread_index[l_th]);
			l_th++;
		}

		for (int i = 0; i < l_th; i++) {
			pthread_join (threads[i], 0);
			f_write(fo, &output_size[i], sizeof(uint32_t));
			f_write(fo, output + i * buffer_size, output_size[i]);
		}
		current_position = 0;
		return;
	}

	int offset = 0;
	while (size != offset) {
		int w = _thread_count * buffer_size - current_position;
		if (size-offset < w) w = size-offset;
		memcpy(input + current_position, data + offset, w);
		if (current_position + w == buffer_size * _thread_count) {
			for (int i = 0; i < _thread_count; i++) { 
				input_size[i] = buffer_size;
				pthread_create(&threads[i], 0, thread_c, (void*) &thread_index[i]);
			}
			for (int i = 0; i < _thread_count; i++) {
				pthread_join (threads[i], 0);
				f_write(fo, &output_size[i], sizeof(uint32_t));
				f_write(fo, output + i * buffer_size, output_size[i]);
			}
			current_position = 0;
		}
		else current_position += w;
		offset += w;
	}
}

void ac_read (buffered_file *fi, uint8_t *data, int size) {
	static uint64_t f_sz = 0;
	static int current_position = 0;

	int offset = 0, w;
	if (!size) goto x;
	while (offset != size) {
	  	w = _thread_count * buffer_size - current_position;
		if (size-offset < w) w = size-offset;
		memcpy(data + offset, output + current_position, w);
		if (current_position + w == buffer_size * _thread_count) {
x:
			if (!size) f_read(fi, &f_sz, sizeof(uint64_t));
			int ct = 0;
			for (; ct < _thread_count; ct++) {
				if (!f_read(fi, &input_size[ct], sizeof(int)))
					break;
				f_read(fi, input + ct * buffer_size, input_size[ct]);
				output_size[ct] = buffer_size < f_sz ? buffer_size : f_sz;
				f_sz -= output_size[ct];
			}
			for (int i = 0; i < ct; i++)
				pthread_create(&threads[i], 0, thread_d, (void*) &thread_index[i]);
			for (int i = 0; i < ct; i++) 
				pthread_join (threads[i], 0);
			current_position = 0;
			if (!size) return;
		}
		else current_position += w;
		offset += w;
	}
}
