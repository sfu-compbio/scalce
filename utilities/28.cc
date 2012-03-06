/// 786

#include "util.h"

int getval (char c) {
	return 
		((c == 'C' || c == '1') 
			? 1 
			: ((c == 'G' || c == '2') 
				? 2 
				: ((c == 'T' || c == '3') ? 3 : 0)
			)
		);
}

int main (int argc, char **argv) {
	E("2/8 encoder; input must be uppercase ACTG w/o newlines (stdin/out)\n");

	uint8_t buffer = 0;
	uint8_t ch;
	int cnt=0;
	while (fread(&ch,1,1,stdin) == 1) {
		ch=getval(ch);
		buffer=(buffer<<2)|ch;
		cnt++;

		if(cnt%4==0) {
			fwrite(&buffer,1,1,stdout);
			buffer=0;
		}
		else;
	}

	if (cnt%4 != 0) {
		while (cnt%4!=0) { buffer<<=2; cnt++; }
		fwrite(&buffer,1,1,stdout);
	}

	return 0;
}

