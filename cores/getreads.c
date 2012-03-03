/// 786

#include <ctime>
#include <cstdio>
#include <string>
#include <cstdlib>
#include <iostream>
#include <algorithm>
using namespace std;

#include <inttypes.h>

#define LOG(c,...)   fprintf(stderr,"[%6d]: "c"\n",clock()/CLOCKS_PER_SEC,##__VA_ARGS__)
#define ERROR(c,...) { fprintf(stderr,"Error! "c"\n",##__VA_ARGS__); exit(1); }
const int MAXLINE=1000;
void uncompress (const char *path, const char *out) {
	int len; int64_t oQ, oN, nl;

	FILE *fR = fopen(path, "rb");
	fread (&len, sizeof(int), 1, fR);
	fread (&oQ, sizeof(int64_t), 1, fR); // just hack
	fread (&oQ, sizeof(int64_t), 1, fR);
	fread (&oN, sizeof(int64_t), 1, fR);
	fread (&nl, sizeof(int64_t), 1, fR);

	FILE *fo = fopen(out,"wb");

	uint8_t buffer[MAXLINE];
	int bc = 0;

	char alphabet[] = "ACGT";
	uint8_t l[MAXLINE], o[MAXLINE], chr, off;
	
	for (long K = 0; K < nl; K++) {
		int W = fread (l, 1, len/4 + (len%4 > 0), fR);
		fread (o, 1, len/8 + (len%8 > 0), fR);
		bc = 0;
		for(int i = 0; i < len; i++) {
			chr = (l[i/4] >> ((3 - i%4) * 2)) & 3;
			off = (o[i/8] >> (7 - i%8)) & 1;

			buffer[bc++] = (off ? 'N' : alphabet[chr]);
		}
		buffer[bc++] = '\n';
		fwrite (buffer, 1, bc, fo);
	}

	LOG("NUMLINES: %lld",nl);
	fclose(fo);
}


int main (int argc, char **argv) {
	uncompress(argv[1],argv[2]);

	LOG("Done!");
	return 0;
}

