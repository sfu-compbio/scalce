/// 786

#ifndef ARITHMETIC_H__
#define ARITHMETIC_H__

#include "buffio.h"

#define AC_DEPTH 50
extern uint64_t 
			ac_freq3[2][AC_DEPTH*AC_DEPTH],    // 16K
		   ac_freq4[2][AC_DEPTH*AC_DEPTH*AC_DEPTH]; // 760k
void a_init_decoder (buffered_file *f);
void a_read (buffered_file *f, uint8_t *arr, int sz);

struct ac_stat {
	uint64_t *acfq3, *acfq4;
	uint32_t hi_fq[AC_DEPTH * AC_DEPTH][AC_DEPTH];
	uint32_t lo_fq[AC_DEPTH * AC_DEPTH][AC_DEPTH];
	ac_stat(){};
	ac_stat(uint64_t*, uint64_t*);
};	
class ac_coder {
private:
	ac_stat *as;
	uint32_t lo, hi, underflow, code;
	uint32_t prev[3]; /* TODO once ...*/
	uint8_t  *ou, *os;
	int      buf_pos;

public:
	ac_coder(uint8_t*, ac_stat*);
	void O(int t);
	void init();
	void write(uint8_t *arr, int sz);
	void flush();

	void reset();
	uint8_t *output() { return ou; }
};
class ac_decoder {
private:
	ac_stat *as;
	uint32_t lo, hi, code;
	uint32_t prev[3];
	uint8_t  *in, *is;
	int      buf_pos;
public:
	ac_decoder(ac_stat*, uint8_t*);
	uint8_t I();
	uint8_t read_single();
	void read(uint8_t*, int);
	void reset();
};

void ac_write (buffered_file *fo, uint8_t *data, int size);
void ac_read (buffered_file *fi, uint8_t *data, int size);
void set_ac_stat(uint64_t *a3, uint64_t *a4);
void ac_init();
void ac_finish();
#endif // ARITHMETIC_H__
