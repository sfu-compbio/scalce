/// 786

#include "util.h"

void order0 (uint64_t *freq, uint64_t sz, int limit) {
	int alsz=0;
	for (int i = 0; i < limit; i++) if (freq[i]) alsz++;

	double entropy = log(alsz)/log(2);
	printf("O0: entropy: %.2lf b/c, total: %.0lf (%.2lfM)\n", entropy, sz * entropy / 8, sz * entropy / (1024 * 1024 * 8));
}

void order1 (uint64_t *freq, uint64_t sz, int limit) {
	double entropy = 0;
	for (int i = 0; i < limit; i++) if (freq[i]) {
		double p = double(freq[i])/sz;
		entropy += log(p) * p;
	}
	entropy *= -1;
	entropy /= log(2);

	printf("O1: entropy: %.2lf b/c, total: %.0lf (%.2lfM)\n", entropy, sz * entropy / 8, sz * entropy / (1024 * 1024 * 8));
}

void order2 (uint64_t *freq, uint64_t **cond, uint64_t sz, int limit) {
	double entropy = 0;
	for (int i = 0; i < limit; i++) if (freq[i]) {
		double p = double(freq[i]) / sz;

		double s = 0;
		for (int j = 0; j < limit; j++) if (cond[i][j]) {
			double c = double(cond[i][j]) / freq[i];
			s += log(c) * c;
		}

		entropy += s * p;
	}
	entropy *= -1;
	entropy /= log(2);

	printf("O2: entropy: %.2lf b/c, total: %.0lf (%.2lfM)\n", entropy, sz * entropy / 8, sz * entropy / (1024 * 1024 * 8));
}

void order3 (uint64_t *freq, uint64_t **cond, uint64_t ***markoff, uint64_t sz, int limit) {
	double entropy = 0;
	for (int i = 0; i < limit; i++) if (freq[i]) {
		double p = double(freq[i]) / sz;
		double s = 0;
		for (int j = 0; j < limit; j++) if (cond[i][j]) {
			double c = double(cond[i][j]) / freq[i];
			double si = 0;
			for (int k = 0; k < limit; k++) if (markoff[i][j][k]) {
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


void order4 (uint64_t *freq, uint64_t **cond, uint64_t ***markoff, uint64_t ****f4, uint64_t sz, int limit) {
	double entropy = 0;
	for (int i = 0; i < limit; i++) if (freq[i]) {
		double p = double(freq[i]) / sz;
		double s = 0;
		for (int j = 0; j < limit; j++) if (cond[i][j]) {
			double c = double(cond[i][j]) / freq[i];
			double si = 0;
			for (int k = 0; k < limit; k++) if (markoff[i][j][k]) {
				double h = double(markoff[i][j][k]) / cond[i][j];
				double sii = 0;
				for (int l = 0; l < limit; l++) if (f4[i][j][k][l]) {
					double g = double(f4[i][j][k][l]) / markoff[i][j][k];
					sii += g * log(g);
				}
				si += h * sii;
			}
			s += c * si;
		}
		entropy += s * p;
	}
	entropy *= -1;
	entropy /= log(2);

	printf("O4: entropy: %.2lf b/c, total: %.0lf (%.2lfM)\n", entropy, sz * entropy / 8, sz * entropy / (1024 * 1024 * 8));
}

int main (int argc, char **argv) {
	E("entropy calculator, first argument [input file]\n");
	srand(time(0));
	
	E("preparing...\n");
	FILE *fi = fopen (argv[1], "r");
	/*uint8_t c, min=255,max=0;
	while (fread (&c, 1, 1, fi) == 1) {
		min=std::min(c,min);
		max=std::max(c,max);
	}
	int limit = max-min+1;
	E("[%d,%d] lim=%d\n",min,max,limit);*/
	int limit=256;

	E("allocating...\n");
	uint64_t *freq1, **freq2, ***freq3, ****freq4;
	freq1 = new uint64_t[limit];
	memset(freq1, 0, limit*sizeof(uint64_t));
	freq2 = new uint64_t*[limit];
	freq3 = new uint64_t**[limit];
//	freq4 = new uint64_t***[limit];
	
	for (int i = 0; i < limit; i++) {
		freq2[i] = new uint64_t[limit];
		memset(freq2[i], 0, limit*sizeof(uint64_t));
		freq3[i] = new uint64_t*[limit];
//		freq4[i] = new uint64_t**[limit];
		for (int j = 0; j < limit; j++) {
			freq3[i][j] = new uint64_t[limit];
			memset(freq3[i][j], 0, limit*sizeof(uint64_t));
//			freq4[i][j] = new uint64_t*[limit];
//			for (int k = 0; k < limit; k++) {
//				freq4[i][j][k] = new uint64_t[limit];
//				memset(freq4[i][j][k], 0, limit*sizeof(uint64_t));
//			}
		}
	}
		
	E("obtaining stats...\n");
	rewind(fi);
	uint8_t pc = 0, ppc = 0, pppc = 0;
	uint64_t sz = 0;
//	while (fread (&c, 1, 1, fi) == 1) {
		
	uint8_t cx[ 100 ],c;
	while ( fread(cx, 1, 100, fi) == 100 ) {
		for (int i = 0; i < 100; i++) {
			c = cx[i]-'!';
//			c -= min;
			freq1[c]++;
			if (sz > 0) freq2[pc][c]++;
			if (sz > 1) freq3[ppc][pc][c]++;

//				if (sz > 2) freq4[pppc][ppc][pc][c]++;

			pppc = ppc;
			ppc = pc;
			pc = c;
			sz++;
		}
	}

/*	for(int i=0;i<limit;i++) printf("%llu,",freq[i]);printf("\n");
	for (int i = 0; i < limit; i++) { printf("{ ");for (int j = 0; j < limit; j++) {
		if(j)printf(" ,");printf("%llu",cond[i][j]);
	} printf(" }, "); }
*/
	E("entropy calc...\n");
	printf("total chars: %llu\n", sz);

	order0(freq1,sz,limit);
	order1(freq1,sz,limit);
	order2(freq1,freq2,sz,limit);
	order3(freq1,freq2,freq3,sz,limit);
//	order4(freq1,freq2,freq3,freq4,sz,limit);

	return 0;
}

