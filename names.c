/// 786

#include <stdlib.h>
#include <string.h>
#include "const.h"
#include "names.h"

char first = 0;
char sanger = 1;
char prefix[200], machine_name[200], run_id[200];

/* encode read name. discards everything after first space */
int output_name (char* name, uint8_t *dest) {
	if (!_use_names)
		return 0;
	#pragma omp critical
	{
	if (!first) {
		char cc[200], cp = 0, ep = 0;
		for (int i = 0; ; i++) 
			if (!name[i] || name[i] == '\n' || name[i] == ' ' || name[i] == ':' || name[i] == '.' || name[i] == '/') {
				cc[cp]=0;
				if (ep == 0)
					strncpy(prefix, cc, cp);
				else if (ep == 1) {
					if (!atoi(cc)) { sanger = 0; break; }
				}
				else if (ep == 2)
					strncpy(machine_name, cc, cp);
				else if (ep == 3)
					strncpy(run_id, cc, cp);
				else if (ep == 4) {
					if (!atoi(cc)) { sanger = 0; break; }
				}
				else if (ep == 5) {
					if (!atoi(cc)) { sanger = 0; break; }
				}
				else if (ep == 6) {
					if (!atoi(cc)) { sanger = 0; break; }
				}
				else {
					if (cc[0] != '1' && cc[0] != '2') { sanger = 0; break; }
				}
				cp = 0;
				ep++;
				if (!name[i] || name[i] == '\n')
					break;
			} 
			else cc[cp++] = name[i];

		if (ep != 8) sanger = 0;
		if (sanger) LOG("Sanger names detected\n");
		first = 1;
	}
	}

	int i;
	for (i = 1; name[i] != 0 && name[i] != '\n' /*&& name[i] != ' '*/; i++)
		dest[i] = name[i];
	dest[0] = i - 1;
	if (_interleave && i>=2 && dest[i-2] == '/' &&  (dest[i-1] == '1'||dest[i-1]=='2'))
		dest[0] -= 2;
	return dest[0] + 1;
}

