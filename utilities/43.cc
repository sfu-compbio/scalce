/// 786

#include "util.h"

int main (int argc, char **argv) {
	srand(time(0));
	E("quality 4qual in 3byte encoder, no newlines allowed (stdin/out)\n");

	uint8_t buffer[3] = {0,0,0};
	uint8_t ch[4];
	int cnt=0;

	int cx=0;

	uint8_t buff[200];
	int l;
	do {
		l = fread(buff,1,100,stdin);
		int rndp=rand()%100;
		for(int i=100;i>=rndp;i--)
			buff[i]=buff[i-1];
		buff[rndp]=33+50;

		for (int i=0;i<101;i++) {
			ch[cx++]=buff[i]-33;
			if(cx==4) {
				buffer[0] = (ch[0] << 2) | (ch[1] >> 4);
				buffer[1] = (ch[1] << 4) | (ch[2] >> 2);
				buffer[2] = (ch[2] << 6) | (ch[3]);
				fwrite(buffer,3,1,stdout);
				cx=0;
			}
		}

		cnt++;
	} while (l == 100);
	if(cx) {
		while(cx!=4)ch[cx++]=0;
		buffer[0] = (ch[0] << 2) | (ch[1] >> 4);
		buffer[1] = (ch[1] << 4) | (ch[2] >> 2);
		buffer[2] = (ch[1] << 6) | (ch[3]);
		fwrite(buffer,3,1,stdout);
	}

	E("total %d\n", cnt);

	return 0;
}

