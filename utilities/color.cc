/// 786

#include "util.h"

int getval (char c) {
	return (c == 'C' ? 1 : (c == 'G' ? 2 : (c == 'T' ? 3 : 0)));
}

int colorspace_table[4][4] = {
	{ 0, 1, 2, 3 },
	{ 1, 0, 3, 2 },
	{ 2, 3, 0, 1 },
	{ 3, 2, 1, 0 }
}; 

int main (int argc, char **argv) {
	E("color-space encoder, no newlines, first arg [read len] (stdin/out)\n");

	int rlen = atoi(argv[1]);

	const int maxsz = 1050;
	char buffer[maxsz];
	int l;
	while ((l = fread (buffer, 1, rlen, stdin)) == rlen) {
		buffer[0] = getval(buffer[0]);
		for(int i=1; i<l; i++)
			buffer[i]=colorspace_table[ buffer[i-1] ][ getval(buffer[i]) ];
		for(int i=0; i<l; i++)
			putc(buffer[i]+'0', stdout);
	}
	

	return 0;
}

