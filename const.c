/// 786

#include <string.h>
#include "const.h"

int tbl[] ={ 0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0 };

int getval (char c) {
	return tbl[c-'A'];
//	return (c == 'C' ? 1 : (c == 'G' ? 2 : (c == 'T' ? 3 : 0)));
}

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
