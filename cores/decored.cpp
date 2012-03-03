/// 786

#include <cstdio>
#include <cstring>
#include <inttypes.h>
using namespace std;

int main(void) {
	char qqq[]="ACGT";
	while(1) {
		int16_t ln;
		if(fread(&ln,1,2,stdin)<2) break;
		int32_t cnt;
		fread(&cnt,1,4,stdin);

		int sz=ln/4+(ln%4!=0);
		fprintf(stderr,"%d\n",sz);
		int64_t x;
		for(int i=0;i<cnt;i++) {
			fread(&x,1,sz,stdin);
			for(int j=ln-1;j>=0;j--) {
				putchar(qqq[  (x>>(2*j))&3   ]);
			}
			putchar('\n');
		}
	}
	return 0;
}

