/// 786

#include <set>
#include <list>
#include "util.h"

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

/** fast radix sorting for dna strings  *
  * params: pos   - in-string position (default 0)
  *         start - start position of partition
  *			size  - partition size
  *			v     - vector to be sorted (vector of pointers to string references!)
  */
void rs (int pos, int start, int size, vector<string*> &v) {
	if (size <= 1)  /* trivial case: subarray of size 1 */
		return;
	if (pos >= v[0]->size()) /* trivial case: exceeded string length */
		return;

	/* obtain character stats in count array */
	int count[4] = {0}; 
	vector<string*> original(size);
	for (int i = start; i < start + size; i++) {
		count[ getval(v[i]->at(pos)) ]++;
		original[i - start] = v[i]; /* keep original ordering of vector safe */
	}

	/* get cumulative sums fot count array */
	int cumulative[5] = {0};
	for (int i = 1; i < 5; i++)
		cumulative[i] = cumulative[i - 1] + count[i - 1]; 

	/* partition the original array according to the starting character */
	for (int i = start; i < start + size; i++) {
		char c = getval(original[i - start]->at(pos));
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

vector<string*> radix (vector<string> &v) {
	vector<string*> vx;
	for (int i = 0; i < v.size(); i++)
		vx.push_back(&v[i]);
	rs(0, 0, vx.size(), vx);
	return vx;
}

void RLE (const vector<string*> &v, int block_size = 1000) {
//	fwrite(&block_size, sizeof(int), 1, stdout);
	for (int k = 0; k < v.size(); k += block_size) {
		for (int p = 0; p < v[0]->size(); p++) { 
			uint8_t prev = v[k]->at(p);
			uint8_t count = 1;
			for (int i = 1; i <= block_size && i + k <= v.size(); i++) {
				if (i < block_size && i + k < v.size() && v[k + i]->at(p) == prev && count < 127)
					count++;
				else {
					fwrite(&prev, 1, 1, stdout);
					if (count > 1) {
						count += 128;
						fwrite(&count, 1, 1, stdout);
					}
					if (k + i < v.size())
						prev = v[k + i]->at(p);
					count = 0;
				}
			}
		}
	}
}

void magic28 (const string &s) {
	uint8_t buffer = 0;
	int i = 0;

	uint8_t x, y;
	while (i < s.size()) {
		if(s[i]=='@') { x = i++; continue; }
		if(s[i]=='#') { y = i++; continue; }

		buffer = (buffer<<2) | getval(s[i]);

		if (i % 4 == 0) 
			fwrite(&buffer, 1, 1, stdout);
		i++;
	}

	if (i % 4 != 0) {
		while (i % 4 != 0) { buffer <<= 2; i++; }
		fwrite(&buffer, 1, 1, stdout);
		fwrite(&x, 1, 1, stdout);
		fwrite(&y, 1, 1, stdout);
	}
}


int main (int argc, char **argv) {
	E("color-space encoder, no newlines, first arg [read len] (stdin/out)\n");

	vector<string> reads;

	char c_core[maxsz], c_buffer[maxsz];
	int  ns, pos;
	while (scanf("%s %d", c_core, &ns) != EOF) {
		int l = strlen(c_core);
		for (int k = 0; k < ns; k++) {
			scanf("%d %s", &pos, c_buffer);
			string s = string(c_buffer);
			if (l > 1)
				s = s.substr(0, pos) + "#" + s.substr(pos + l) + "@";
			else
				s += "@";
			reads.push_back(shift(s));
		}

		vector<string*> r = radix(reads);
		RLE(r);
//		for (int k = 0; k < r.size(); k++)
//			magic28(*r[k]);
//			printf("%s\n", r[k]->c_str());
		reads.clear();
	}

	return 0;
}

