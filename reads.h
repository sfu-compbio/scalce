/// 786

#ifndef READS_H__
#define READS_H__

#include "const.h"
#include "buffio.h"

/* read_data - contains all neccessary data for one read (name, quality, read and eventually paired read data */
typedef struct {
	uint8_t *data;   /* read data (with name, quality and so on) */
	int32_t sz;      /* size of data */
	int32_t of;      /* offset in data for paired read */
} read_data;

/* bin_node - linked list node for read_data */
typedef struct bin_node bin_node;
struct bin_node {
	read_data data;  /* read data */
	bin_node  *next; /* neighbor node */
};

/* bin - linked list wrapper */
typedef struct bin bin;
struct bin {
	bin_node *first, *last;  /* pointers to first/last element */
	int32_t size;            /* number of elements */
};

/* aho_trie - trie and aho automaton for core strings */
typedef struct aho_trie aho_trie;
struct aho_trie {
	aho_trie *child[4];        /* trie nodes (A, C, T, G) */
	aho_trie *fail;            /* aho automaton fail pointer */
	aho_trie *next_to_output;  /* pointer to "output" node */
	int level;                 /* node level in trie */
	int id;                    /* node identificator */
	uint64_t bin_size;         /* total number of reads bucketed so far */
	char output;               /* is this node "output" node? */
	bin bin;                   /* bin linked list */
};

aho_trie *read_patterns ();
void aho_search (char *text, aho_trie *root, aho_trie **bucket);
void aho_trie_bucket (aho_trie *t, read_data *d);
int output_read (char *line, uint8_t *dest);
void aho_output (aho_trie *root, buffered_file *f);
void aho_trie_free (aho_trie *t);
int unbuck (void);
#endif // READS_H__

