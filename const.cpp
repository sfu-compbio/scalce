/// 786

#include <time.h>
#include <sys/time.h>
#include <string.h>
#include "const.h"

int _tbl[] ={ 0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0 };

char *get_second_file (const char *c) {
	static char buf[MAXLINE];
	strncpy (buf, c, MAXLINE);
	for (int i = strlen(buf) - 1; i >= 0; i--)
		if (buf[i] == '1') { 
			buf[i] = '2';
			break;
		}

	if (!strcmp (buf, c))
		return 0;

	return buf;
}
int _TIME_()
{
	struct timeval t;
	gettimeofday(&t,0);
	return (t.tv_sec*1000000+t.tv_usec);
}


int parseLine(char* line){
	int i = strlen(line);
	while (*line < '0' || *line > '9') line++;
	line[i-3] = '\0';
	i = atoi(line);
	return i;
}

void MEM(char *A){
/*	FILE* file = fopen("/proc/self/status", "r");
	int result = -1;
	char line[128];
	while (fgets(line, 128, file) != NULL){
		if (strncmp(line, "VmSize:", 7) == 0) result = parseLine(line);
		break;
	}
	fclose(file);
	
	
	LOG("%s: %d\n", A, result);*/
}

static uint64_t mem_usg = 0;
void *mallox (size_t size) {
	mem_usg += size;
	void *v = malloc(size);

	if (!v) { ERROR("mallox failed whoa whoa whoa!\n");  return 0; }
	else return v;
}
void frex (void *ptr, size_t size) {
	mem_usg -= size;
	free(ptr);
}
double getmemx () {
	return mem_usg / (1024*1024.0);
}


