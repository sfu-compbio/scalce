/// 786

#include "buffio.h"
#include <cstdio>
#include <inttypes.h>
#include <clocale>
#include <cstring>
using namespace std;

const int blocksz = 4096;
char block[blocksz];

void read_stdio() {
	FILE *fi = fopen("test", "r");
	
	uint64_t sz = 0, p;
	while (p = fread(block, 1, blocksz, fi)) sz += p;
	printf("read %'llu bytes\n", sz);
	fclose(fi);
}


buffered_file xfi;
void read_myio() {
	f_open(&xfi, "test", IO_READ);
	uint64_t sz = 0, p;
	while (p = f_read(&xfi, block, blocksz)) sz += p;
	printf("read %'llu bytes\n", sz);
	f_close(&xfi);
}


void gets_stdio() {
	FILE *fi = fopen("test", "r");
	uint64_t sz = 0;
	while (fgets(block, blocksz, fi)) sz += strlen(block);
	printf("read %'llu bytes\n", sz);
	fclose(fi);
}


void gets_myio() {
	f_open(&xfi, "test", IO_READ);
	uint64_t sz = 0;
	while (f_gets(&xfi, block, blocksz)) sz += strlen(block)+1;
	printf("read %'llu bytes\n", sz);
	f_close(&xfi);
}

const uint64_t write_sz = 6ll * 1024 * 1024 * 1024;
void write_stdio() {
	FILE *fi = fopen("testw", "wb");
	uint64_t sz = write_sz;
	while (sz) sz -= fwrite(block, 1, blocksz, fi);
	fclose(fi);
}
void write_myio() {
	f_open(&xfi, "testw", IO_WRITE);
	uint64_t sz = write_sz;
	while (sz) sz -= f_write(&xfi, block, blocksz);
	f_close(&xfi);
}

uint64_t _file_buffer_size = 20 * 1024 * 1024;

int main(int argc, char **argv) {
	setlocale(LC_ALL,"");
	int t;

	_file_buffer_size = 1024 * 1024 * atof(argv[1]);
	printf("myion buffer size %'llu\n", _file_buffer_size);

	f_init(&xfi, IO_SYS);

	printf("read\n");
	t = TIME;
	read_stdio();
	printf("stdio: %'12d\n", TIME-t);
	t = TIME;
	read_myio();
	printf("myio:  %'12d\n", TIME-t);

	printf("gets\n");
	t = TIME;
	gets_stdio();
	printf("stdio: %'12d\n", TIME-t);
	t = TIME;
	//gets_myio();
	printf("myio:  %'12d\n", TIME-t);

	printf("write %'llu bytes\n", write_sz);
	t = TIME;
	write_stdio();
	printf("stdio: %'12d\n", TIME-t);
	t = TIME;
	write_myio();
	printf("myio:  %'12d\n", TIME-t);


	return 0;
}

