/// 786

#include <iostream>
#include <cstdio>
#include <vector>
using namespace std;

vector<int64_t> v;

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
	char l[50];
	int pl=0;
	bool f=false;
	while (fgets(l,50,stdin)){
		int xl=strlen(l)-1;

		if(f&&pl!=xl) {
			dox(pl);
			pl=xl;
		}
		v.push_back(cnt(l));
		if(!f) {
			f=1;
			pl=xl;
		}
	}
	dox(pl);

	return 0;
}
