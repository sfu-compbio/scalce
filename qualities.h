/// 786

#ifndef QUALITIES_H__
#define QUALITIES_H__

#include "buffio.h"

/* quality_mapping - keeps lossy transformation information for qualities */
typedef struct {
	int offset;        /* phred quality offset - 33 or 64 */
	int values[128];   /* replacement table */
} quality_mapping;

void quality_mapping_init (quality_mapping *q, buffered_file *f, int *read_length); 
int output_quality (char *line, char *read, quality_mapping *q, uint8_t *dest);

#endif // QUALITIES_H__

