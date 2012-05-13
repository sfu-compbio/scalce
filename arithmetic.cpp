/// 786

#include <assert.h>
#include <inttypes.h>
#include "arithmetic.h"

static uint32_t hi_fq[2][46*46][46],   // 760k 
					 lo_fq[2][46*46][46];   // 760k
static uint32_t lo[2], 
					 hi[2], 
					 underflow[2],
					 code[2],
					 prev[2][3];

/* bit i/o routines */
static uint8_t  buffer[2];
static int      buf_pos[2];

int Z=0;

uint8_t I (buffered_file *f) {
	uint8_t c = (buffer[Z] >> buf_pos[Z]) & 1;
	if (!buf_pos[Z]) {
		if (f_read(f, &buffer[Z], 1) != 1)
			buffer[Z] = 0;
		buf_pos[Z] = 7;
	}
	else buf_pos[Z]--;
	return (c ? 1 : 0);
}

void O (buffered_file *f, int t) {
	buffer[Z] <<= 1;
	if (t) buffer[Z] |= 1;
	if (buf_pos[Z] == 0) {
		f_write (f, &buffer[Z], sizeof(uint8_t));
		buf_pos[Z] = 7;
		buffer[Z] = 0;
	}
	else --buf_pos[Z];
}

/* initialization */
void a_init_coder (buffered_file *f) {
//	double ratio = 1;
//	LOG("%llu --",ac_sz);
//	if (ac_sz >= (1ll << 32))
//		ratio = ((1ll << 32) - 256*256 - 5) / (double)(ac_sz);
//	ac_sz *= ratio;
/*
	for (int i = 0; i < 256; i++) {
//		ac_freq1[i] *= ratio;
//		ac_freq1[i] += 256;
//		ac_sz += 256;
		ac_freq1[i] = 0;
		for (int j = 0; j < 256; j++) {
			ac_freq2[i * 256 + j] *= ratio;
			ac_freq2[i * 256 + j] ++;
			ac_freq1[i] += ac_freq2[i*256+j];
		}
	}
*/
	for (int i = 0; i < 46; i++) {
		for (int j = 0; j < 46; j++) {
				hi_fq[Z][i*46 + j][0] = ac_freq4[Z][i*46 + j][0];
				for (int l = 1; l < 46; l++) 
					hi_fq[Z][i*46 + j][l] = hi_fq[Z][i*46 + j][l - 1] + ac_freq4[Z][i*46 + j][l];
				ac_freq3[Z][i*46 + j] = hi_fq[Z][i*46 + j][45];
				int p = -1;
				for (int l = 0; l < 46; l++) if (ac_freq4[Z][i*46 + j][l]) {
					if (p != -1) 
						lo_fq[Z][i*46 + j][l] = hi_fq[Z][i*46 + j][p];
					p = l;
				}
		}
	}

	for (int i = 0; i < 46; i++) 
		for (int j = 0; j < 46; j++) 
			for (int k = 0; k < 46; k++) 
				f_write(f, &ac_freq4[Z][i*46 + j][k], sizeof(uint32_t));	

	lo[Z] = 0;
	hi[Z] = -1;
	underflow[Z] = 0;
	prev[Z][0] = prev[Z][1] = prev[Z][2] = 500;	
	buffer[Z] = 0;
	buf_pos[Z] = 7;
	//LOG("----");
}

void a_write (buffered_file *f, uint8_t *arr, int sz) {
	int i = 0;
	if (prev[Z][0] > 255) { // first 3 char hack
		f_write (f, &arr[0], sizeof(uint8_t));
		f_write (f, &arr[1], sizeof(uint8_t));
//		f_write (f, &arr[2], sizeof(uint8_t));
		prev[Z][0] = arr[0];
		prev[Z][1] = arr[1];
//		prev[Z][2] = arr[2];
		i = 2;
	}
	for (; i < sz; i++) {
		uint32_t c_lo = lo_fq[Z][prev[Z][0]*46 + prev[Z][1]][arr[i]],
					c_hi = hi_fq[Z][prev[Z][0]*46 + prev[Z][1]][arr[i]],
					p_sz = ac_freq3[Z][prev[Z][0]*46 + prev[Z][1]];
		uint64_t range = hi[Z] - lo[Z] + 1LL;
		hi[Z] = lo[Z] + (range * c_hi) / p_sz - 1;
		lo[Z] = lo[Z] + (range * c_lo) / p_sz;

		while (1) {
			if ((hi[Z] & 0x80000000) == (lo[Z] & 0x80000000)) {
				O (f, hi[Z] & 0x80000000);
				while (underflow[Z]) {
					O (f, ~hi[Z] & 0x80000000); 
					underflow[Z]--;
				}
			}
			else if (!(hi[Z] & 0x40000000) && (lo[Z] & 0x40000000)) { // underflow[Z] ante portas
				underflow[Z]++;
				lo[Z] &= 0x3fffffff; // set 2nd highest digit to 0 -- 1st will be removed after
				hi[Z] |= 0x40000000; // set 2nd highest digit to 1
			}
			else break;

			lo[Z] <<= 1;
			hi[Z] <<= 1; hi[Z] |= 1;
		}
	
		prev[Z][0] = prev[Z][1];
	//	prev[Z][1] = prev[Z][2];
		prev[Z][1] = arr[i];
	}
}

void a_finish_coder (buffered_file *f) {
	O (f, lo[Z] & 0x40000000);
	underflow[Z]++;
	while (underflow[Z]) {
		O (f, ~lo[Z] & 0x40000000);
		underflow[Z]--;
	}
	while (buf_pos[Z] != 7)
		O (f, 0);
}


void a_init_decoder (buffered_file *f) {
	for (int i = 0; i < 46; i++) 
		for (int j = 0; j < 46; j++) 
			for (int k = 0; k < 46; k++) 
				f_read(f, &ac_freq4[Z][i*46 + j][k], sizeof(uint32_t));	

	for (int i = 0; i < 46; i++) {
		for (int j = 0; j < 46; j++) {
				hi_fq[Z][i*46 + j][0] = ac_freq4[Z][i*46 + j][0];
				for (int l = 1; l < 46; l++) 
					hi_fq[Z][i*46 + j][l] = hi_fq[Z][i*46 + j][l - 1] + ac_freq4[Z][i*46 + j][l];
				ac_freq3[Z][i*46 + j] = hi_fq[Z][i*46 + j][45];
				int p = -1;
				for (int l = 0; l < 46; l++) if (ac_freq4[Z][i*46 + j][l]) {
					if (p != -1) 
						lo_fq[Z][i*46 + j][l] = hi_fq[Z][i*46 + j][p];
					p = l;
				}
			
		}
	}

	//f_read(f, &prev[Z], sizeof(uint8_t));
	prev[Z][0] = prev[Z][1] = prev[Z][2] = 500;	
	lo[Z] = 0;
	hi[Z] = -1;
}

uint8_t a_read_single (buffered_file *f) {
	uint8_t c = 0;

	uint64_t range = hi[Z] - lo[Z] + 1LL;
	uint32_t count = ((code[Z] - lo[Z] + 1LL) * 
			ac_freq3[Z][prev[Z][0]*46 + prev[Z][1]] - 1LL) / range;

	uint32_t i;
	for (i = 0; i < 46; i++) 
		if (ac_freq4[Z][prev[Z][0]*46 + prev[Z][1]][i] 
				&& count >= lo_fq[Z][prev[Z][0]*46 + prev[Z][1]][i] 
				&& count <  hi_fq[Z][prev[Z][0]*46 + prev[Z][1]][i]) 
			break;
//	if(i>=46) {
	//	LOG("\n\n %u %u[%c%c%c] ", count, hi_fq[Z][prev[Z][0]*46*46+prev[Z][1]*46+prev[Z][2]][45],prev[Z][0]+'!',prev[Z][1]+'!',prev[Z][2]+'!');
//	}
	assert(i<46);

	hi[Z] = lo[Z] + (range * hi_fq[Z][prev[Z][0]*46 + prev[Z][1]][i]) / ac_freq3[Z][prev[Z][0]*46 + prev[Z][1]] - 1;
	lo[Z] = lo[Z] + (range * lo_fq[Z][prev[Z][0]*46 + prev[Z][1]][i]) / ac_freq3[Z][prev[Z][0]*46 + prev[Z][1]];
	while (1) {
		if ((hi[Z] & 0x80000000) == (lo[Z] & 0x80000000)) ;
		else if (!(hi[Z] & 0x40000000) && (lo[Z] & 0x40000000)) { // underflow[Z] ante portas
			code[Z] ^= 0x40000000;
			lo[Z]   &= 0x3fffffff; 
			hi[Z]   |= 0x40000000; 
		}
		else break;

		lo[Z] <<= 1;
		hi[Z] <<= 1; hi[Z] |= 1;
		code[Z] <<= 1; code[Z] |= I(f);
	}

	prev[Z][0] = prev[Z][1];
//	prev[Z][1] = prev[Z][2];
	prev[Z][1] = i;

//	LOG("%c",i+'!');
	return i;
}

void a_read (buffered_file *f, uint8_t *ar, int sz) {
	int i = 0;
	if (prev[Z][0] > 255) {
		f_read(f, &ar[0], sizeof(uint8_t));
		f_read(f, &ar[1], sizeof(uint8_t));
		//f_read(f, &ar[2], sizeof(uint8_t));
		prev[Z][0] = ar[0];
		prev[Z][1] = ar[1];
		//prev[Z][2] = ar[2];
		i = 2;

		f_read(f, &buffer[Z], sizeof(uint8_t));
		buf_pos[Z] = 7;
		code[Z] = 0;
		for (int j = 0; j < 32; j++) {
			code[Z] <<= 1;
			code[Z] |= I(f);
		}
	} 

	for (; i < sz; i++)
		ar[i] = a_read_single(f);
}

