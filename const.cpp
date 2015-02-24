/*
 * Copyright (c) 2011 - 2012, Simon Fraser University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, 
 * are permitted provided that the following conditions are met:
 *   
 * Redistributions of source code must retain the above copyright notice, this list
 * of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or other
 *   materials provided with the distribution.
 * - Neither the name of the Simon Fraser University nor the names of its contributors may be
 *   used to endorse or promote products derived from this software without specific
 *   prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * Author         : Ibrahim Numanagic
 * Email          : inumanag AT sfu DOT ca
 * Last Update    : 25. vii 2012.
 */

#include <time.h>
#include <sys/time.h>
#include <string.h>
#include "const.h"

int _tbl[] ={ 
	0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0,
	0,0,0,0,0,0,
	0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0 
};

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
int64_t _TIME_()
{
	struct timeval t;
	gettimeofday(&t,0);
	return (t.tv_sec*1000000ll+t.tv_usec);
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


