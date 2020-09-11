/*
 * Copyright (c) 2011 - 2012, Simon Fraser University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list
 * of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice,
 * this
 *   list of conditions and the following disclaimer in the documentation and/or
 * other
 *   materials provided with the distribution.
 * - Neither the name of the Simon Fraser University nor the names of its
 * contributors may be
 *   used to endorse or promote products derived from this software without
 * specific
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

#ifndef READS_H__
#define READS_H__

#include "buffio.h"
#include "const.h"

extern char **patterns;

/* read_data - contains all neccessary data for one read (name, quality, read
 * and eventually paired read data */
typedef struct {
  uint8_t *data; /* read data (with name, quality and so on) */
  int32_t sz;    /* size of data */
  int32_t of;    /* offset in data for paired read */
  int16_t end;   /* end marker */
  int32_t read_length, read_length_2;
} read_data;

/* bin_node - linked list node for read_data */
typedef struct bin_node bin_node;
struct bin_node {
  read_data data; /* read data */
  bin_node *next; /* neighbor node */
};

/* bin - linked list wrapper */
typedef struct bin bin;
struct bin {
  bin_node *first, *last; /* pointers to first/last element */
  int32_t size;           /* number of elements */
};

/* aho_trie - trie and aho automaton for core strings */
typedef struct aho_trie aho_trie;
struct aho_trie {
  aho_trie *child[4];       /* trie nodes (A, C, T, G) */
  aho_trie *fail;           /* aho automaton fail pointer */
  aho_trie *next_to_output; /* pointer to "output" node */
  int level;                /* node level in trie */
  int id;                   /* pattern identificator */
  uint64_t bin_size;        /* total number of reads bucketed so far */
  int32_t output;           /* is this node "output" node? */
  struct bin bin;           /* bin linked list */
};

aho_trie *read_patterns();
aho_trie *read_patterns_from_file(const char *f);
int aho_search(char *text, aho_trie *root, aho_trie **bucket);
bin_node *aho_trie_bucket(aho_trie *t, read_data *d);
int output_read(char *line, uint8_t *dest, int n, int l);
void aho_output(aho_trie *root, buffered_file *f);
void aho_trie_free(aho_trie *t);
int unbuck(void);
#endif // READS_H__
