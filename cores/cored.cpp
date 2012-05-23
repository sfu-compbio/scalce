/// 786

#include <iostream>
#include <cstdio>
#include <vector>
#include <cstring>
#include <inttypes.h>
using namespace std;

vector<uint64_t> v;

void dox(int l) {
	fwrite(&l,1,2,stdout);
	int32_t x=v.size();
	fwrite(&x,1,4,stdout);

	int sz=l/4+(l%4!=0);
	for(int i=0;i<x;i++)
		fwrite(&(v[i]),1,sz,stdout);
	v.clear();
}

int64_t cnt(char *c){
	int64_t r=0;
	for(int i=0;i<strlen(c)-1;i++) {
		r=(r<<2)| ( c[i]=='C'?1:(c[i]=='G'?2:(c[i]=='T'?3:0)));
	}
	return r;
}

int main (void) {
	int pl=8;

	uint64_t l, cr;
	while (scanf("%llu %llu", &l, &cr) != EOF){
		if(pl!=l) {
			dox(pl);
			pl=l;
		}
		v.push_back( cr );
	}
	dox(pl);

	return 0;
}
