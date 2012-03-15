/// 786

#include <set>
#include <list>
#include "util.h"

struct Ps {
	string first;
	string second;
	short pos;
	Ps(const string &a,const string &b,short i):first(a),second(b),pos(i){}
};

//typedef P Ps;

inline int getval (char c) {
	return (c == 'C' ? 1 : (c == 'G' ? 2 : (c == 'T' ? 3 : 0)));
}

const int maxsz = 105;

string shift (const string &s) {
	static int low_p[maxsz];
	int low_p_sz = 0;
	uint64_t low = -1ll;

	uint64_t current = 0;
	for (int i = 0; i < s.size() + 32; i++) {
		current = (current << 2) | getval(s[i % s.size()]);
		if (i >= 31 && current < low) {
		//	if (current < low) low_p_sz = 0;
			low = current;
			low_p[0/*low_p_sz++*/] = i - 31;
		}
	}

//	if (low_p_sz == 1) {
	return s.substr(low_p[0]) + s.substr(0,low_p[0]);
//	}
//	else {
//		throw "yay!";
//	}
}

string shift (const string &s, int p, int l) {
	return s.substr(p + l) + s.substr(0, p);
}

/** fast radix sorting for dna strings  *
  * params: pos   - in-string position (default 0)
  *         start - start position of partition
  *			size  - partition size
  *			v     - vector to be sorted (vector of pointers to string references!)
  */
void rs (int pos, int start, int size, vector<Ps*> &v) {
	if (size <= 1)  /* trivial case: subarray of size 1 */
		return;
	if (pos >= v[0]->first.size()) /* trivial case: exceeded string length */
		return;

	/* obtain character stats in count array */
	int count[4] = {0}; 
	vector<Ps*> original(size);
	for (int i = start; i < start + size; i++) {
		count[ getval(v[i]->first.at(pos)) ]++;
		original[i - start] = v[i]; /* keep original ordering of vector safe */
	}

	/* get cumulative sums for count array */
	int cumulative[5] = {0};
	for (int i = 1; i < 5; i++)
		cumulative[i] = cumulative[i - 1] + count[i - 1]; 

	/* partition the original array according to the starting character */
	for (int i = start; i < start + size; i++) {
		char c = getval(original[i - start]->first.at(pos));
		v[start + cumulative[c]] = original[i - start];
		cumulative[c]++; /* here, "fake" update cumulative sum in order to keep track of partition size */
	}

	/* update cumulatove sums */
	cumulative[0] = 0;
	for (int i = 1; i < 5; i++)
		cumulative[i] = cumulative[i - 1] + count[i - 1];

	/* sort partitions */
	for (int i = 0; i < 4; i++)
		rs(pos + 1, start + cumulative[i], cumulative[i + 1] - cumulative[i], v);
}

vector<Ps*> radix (vector<Ps> &v) {
	vector<Ps*> vx;
	for (int i = 0; i < v.size(); i++)
		vx.push_back(&v[i]);
	rs(0, 0, vx.size(), vx);
	return vx;
}

void magic28 (const string &s, FILE *f) {
	uint8_t buffer = 0, bsz = 0;

/*	if (!s.size()) {
		while (bsz != 4) {
			buffer <<= 2;
			bsz++;
		}
		fwrite(&buffer, 1, 1, stdout);
		return;
	}*/

	for (int i = 0; i < s.size(); i++) {
//		if(s[i]=='@') { x = i++; continue; }
//		if(s[i]=='#') { y = i++; continue; }
		buffer = (buffer<<2) | getval(s[i]);
		bsz++;
		if (bsz == 4) { 
			fwrite(&buffer, 1, 1, f);
			buffer = bsz = 0;
		}
	}
	if (bsz) {
		while (bsz != 4) { buffer <<= 2; bsz++; }
		fwrite(&buffer, 1, 1, f);
	}
}

void magicQl (const string &s, FILE *f) {
	for(int i=0;i<s.size();i++){
		uint8_t c=s[i]-'!';
		fwrite(&c,1,1,f);
	}
//	for(int i=0; i<s.size(); i++)
}


int main (int argc, char **argv) {
	E("color-space encoder, no newlines, first arg [read len] (stdin/out)\n");

	FILE *fr = fopen(argv[1], "r"),
		  *fq = fopen(argv[2], "r"),
		  *fx = fopen("_R_", "wb"),
		  *fy = fopen("_Q_", "wb"),
		  *ft = fopen("_T_", "wb");

	vector<Ps> reads;

	char c_core[maxsz], c_buffer[maxsz];
	int  ns, pos;
	while (fscanf(fr, "%s %d", c_core, &ns) != EOF) {
		int l = strlen(c_core);
		for (int k = 0; k < ns; k++) {
			fscanf(fr, "%d %s", &pos, c_buffer);
			string s = string(c_buffer);
			fread(c_buffer, 1, 100, fq);
			for (int i=0;i<100; i++) {
				if(s[i]=='N') { c_buffer[i]='!'; s[i]='A'; }
			}
			if (l > 1) {
				s = shift(s, pos, l);
				assert(s.size() == 100 - l);
			}
			//	s = s.substr(0, pos) + "#" + s.substr(pos + l) + "@";
			else // root case
				;
			//	s += "@";
			short px=(l>1)?(100-pos-l):(100);
		//	E("%d\n",px);
			reads.push_back(Ps(s, string(c_buffer), px));

		}

		vector<Ps*> r = radix(reads);
//		RLE(r);
		fwrite(&ns, sizeof(int), 1, fx);   // bucketsz
		char _l=l; fwrite(&_l, 1, 1, fx);  // corelen
		for (int k = 0; k < r.size(); k++) {
			magic28(r[k]->first,  fx);
			magicQl(r[k]->second, fy);
			fwrite(&r[k]->pos, sizeof(short), 1, ft);
		}
//			printf("%s\n", r[k]->c_str());
//		magic28("");
		reads.clear();
	}

	fclose(fq); fclose(fr); fclose(fx); fclose(fy); fclose(ft);

	return 0;
}

