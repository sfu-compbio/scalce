/// 786

#include <stdlib.h>
#include <string.h>
#include "const.h"
#include "names.h"

/* encode read name. discards everything after first space */
int output_name (char* name, uint8_t *dest) {
	if (!_use_names) {
		dest[0] = 0;
		return 1;
	}

	int i;
	for (i = 1; name[i] != '\n' && name[i] != ' '; i++)
		dest[i] = name[i];
	dest[0] = i - 1;
	if (_interleave && i>=2 && dest[i-2] == '/' &&  (dest[i-1] == '1'||dest[i-1]=='2'))
		dest[0] -= 2;
	return dest[0] + 1;
}

