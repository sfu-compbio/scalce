/// 786

#include "util.h"

int main (int argc, char **argv) {
	E("quality 4qual in 3byte encoder, no newlines allowed (stdin/out)\n");

	uint8_t buffer[3] = {0,0,0};
	uint8_t ch[4];
	int cnt=0;
	while (fread(&ch[0],1,1,stdin) == 1) {
		if (fread(&ch[1],1,1,stdin) != 1) ch[1]=0;
		if (fread(&ch[2],1,1,stdin) != 1) ch[2]=0;
		if (fread(&ch[3],1,1,stdin) != 1) ch[3]=0;

		buffer[0] = (ch[0] << 2) | (ch[1] >> 6);
		buffer[1] = (ch[1] << 4) | (ch[2] >> 2);
		buffer[2] = (ch[1] << 6) | (ch[3]);

		fwrite(buffer,3,1,stdout);
	}

	return 0;
}

