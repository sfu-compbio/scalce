/// 786

#include "util.h"

void order0 (uint64_t *freq, uint64_t sz) {
	int alsz=0;
	for (int i = 0; i < 256; i++) if (freq[i]) alsz++;

	double entropy = log(alsz)/log(2);
	printf("O0: entropy: %.2lf b/c, total: %.0lf (%.2lfM)\n", entropy, sz * entropy / 8, sz * entropy / (1024 * 1024 * 8));
}

void order1 (uint64_t *freq, uint64_t sz) {
	double entropy = 0;
	for (int i = 0; i < 256; i++) if (freq[i]) {
		double p = double(freq[i])/sz;
		entropy += log(p) * p;
	}
	entropy *= -1;
	entropy /= log(2);

	printf("O1: entropy: %.2lf b/c, total: %.0lf (%.2lfM)\n", entropy, sz * entropy / 8, sz * entropy / (1024 * 1024 * 8));
}

void order2 (uint64_t *freq, uint64_t **cond, uint64_t sz) {
	double entropy = 0;
	for (int i = 0; i < 256; i++) if (freq[i]) {
		double p = double(freq[i]) / sz;

		double s = 0;
		for (int j = 0; j < 256; j++) if (cond[i][j]) {
			double c = double(cond[i][j]) / freq[i];
			s += log(c) * c;
		}

		entropy += s * p;
	}
	entropy *= -1;
	entropy /= log(2);

	printf("O2: entropy: %.2lf b/c, total: %.0lf (%.2lfM)\n", entropy, sz * entropy / 8, sz * entropy / (1024 * 1024 * 8));
}

void order3 (uint64_t *freq, uint64_t **cond, uint64_t ***markoff, uint64_t sz) {
	double entropy = 0;
	for (int i = 0; i < 256; i++) if (freq[i]) {
		double p = double(freq[i]) / sz;
		double s = 0;
		for (int j = 0; j < 256; j++) if (cond[i][j]) {
			double c = double(cond[i][j]) / freq[i];
			double si = 0;
			for (int k = 0; k < 256; k++) if (markoff[i][j][k]) {
				double h = double(markoff[i][j][k]) / cond[i][j];
				si += h * log(h);
			}
			s += c * si;
		}
		entropy += s * p;
	}
	entropy *= -1;
	entropy /= log(2);

	printf("O3: entropy: %.2lf b/c, total: %.0lf (%.2lfM)\n", entropy, sz * entropy / 8, sz * entropy / (1024 * 1024 * 8));
}

int main (int argc, char **argv) {
	E("entropy calculator, first argument [input file]\n");
	
	uint64_t *freq1, **freq2, ***freq3;
	freq1 = new uint64_t[256];
	memset(freq1, 0, 256*sizeof(uint64_t));
	freq2 = new uint64_t*[256];
	freq3 = new uint64_t**[256];
	for (int i = 0; i < 256; i++) {
		freq2[i] = new uint64_t[256];
		memset(freq2[i], 0, 256*sizeof(uint64_t));
		freq3[i] = new uint64_t*[256];
		for (int j = 0; j < 256; j++) {
			freq3[i][j] = new uint64_t[256];
			memset(freq3[i][j], 0, 256*sizeof(uint64_t));
		}
	}

	FILE *fi = fopen (argv[1], "r");
	uint8_t c, pc = 0, ppc = 0;
	uint64_t sz = 0;
	while (fread (&c, 1, 1, fi) == 1) {
		freq1[c]++;
		if (sz > 0) freq2[pc][c]++;
		if (sz > 1) freq3[ppc][pc][c]++;
		ppc = pc;
		pc = c;
		sz++;
	}

/*	for(int i=0;i<256;i++) printf("%llu,",freq[i]);printf("\n");
	for (int i = 0; i < 256; i++) { printf("{ ");for (int j = 0; j < 256; j++) {
		if(j)printf(" ,");printf("%llu",cond[i][j]);
	} printf(" }, "); }
*/
	printf("total chars: %llu\n", sz);

	order0(freq1,sz);
	order1(freq1,sz);
	order2(freq1,freq2,sz);
	order3(freq1,freq2,freq3,sz);

	return 0;
}

