/// 786

#include <map>
#include <set>
#include <list>
#include "util.h"

struct fq_rec {
	string read, qual;
};
map<string, fq_rec> dict;

int main (int argc, char **argv) {
	E("FASTQ reordering checker\n");

	FILE *fr = fopen(argv[1], "r"),
		  *fq = fopen(argv[2], "r");

	int frc=0;
	char buff[500];
	while (fgets(buff, 500, fr)) {
		string name = buff;
		fq_rec f;
		fgets(buff, 500, fr); f.read = buff;
		fgets(buff, 500, fr);
		fgets(buff, 500, fr); f.qual = buff;
		dict[name] = f;
		frc++;
	}
	fclose(fr);

	int fqc=0;
	while (fgets(buff, 500, fq)) {
		string name = buff;
		fq_rec f;
		fgets(buff, 500, fq); f.read = buff;
		fgets(buff, 500, fq);
		fgets(buff, 500, fq); f.qual = buff;

		map<string, fq_rec>::iterator i = dict.find(name);
		if (i == dict.end() || (i->second.read != f.read && i->second.qual != f.qual)) {
			E("%s\t%s\t%s\t%s\t%s", name.c_str(), i->second.read.c_str(), i->second.qual.c_str(), f.read.c_str(), f.qual.c_str());
			exit(0);
		}
		fqc++;
	}
	fclose(fq);
	if(frc!=fqc) E("A %d B %d\n", frc, fqc);
	else E("Match!\n");

	return 0;
}

