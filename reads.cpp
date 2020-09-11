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

#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "const.h"
#include "reads.h"

char **patterns;
static int pattern_c = 0;

void bin_prepare(aho_trie *t);

uint8_t *_pool;
uint64_t _pool_pos;
void *getpool(uint64_t sz) {
  void *p = _pool + _pool_pos;
  _pool_pos += sz;
  if (_pool_pos >= _max_bucket_set_size + 500 * MAXLINE)
    ERROR("pool inconsistent");
  return p;
}
void clrpool() { _pool_pos = 0; }

/* initialize bin */
void bin_init(bin *b) {
  b->first = b->last = 0;
  b->size = 0;
}

/* insert read_data into bin (linked list) */
void bin_insert(bin *b, read_data *d) {}

/* safely dispose bin */
void bin_free(bin *b) {
  //	if (b->size) for (bin_node *n = b->first; n; ) {
  //		bin_node *x = n->next;
  //		frex(n->data.data, sizeof(uint8_t)*n->data.sz);
  //		frex(n, sizeof(struct bin_node));
  //		n = x;
  //	}
  b->first = b->last = 0;
  b->size = 0;
}

/* dump bin contents to file
 * f must contain 5 pointers, for names, reads and qualities files
 * (and reads and qualities for paired end) */
void bin_dump(aho_trie *t, buffered_file *f) {
  int64_t tN = 0, tQ = 0, tR = 0, tR2 = 0, tQ2 = 0, tL = 0;
  int64_t nN = 0, nR = 0, nQ = 0, nR2 = 0, nQ2 = 0, nL = 0;

  size_t lenCore = 0;
  if (t->output >= 0)
    lenCore = strlen(patterns[t->output]);
#ifndef PACBIO
  int sz_read;
  if (t->output >= 0)
    sz_read = SZ_READ(read_length[0] - lenCore);
  else
    sz_read = SZ_READ(read_length[0]);
#endif

  int sz_meta = 1;
  if (read_length[0] > 255)
    sz_meta = 2;

#ifdef PACBIO
  sz_meta = 4;
#endif

  for (bin_node *n = t->bin.first; n; n = n->next) {
    read_data *r = &(n->data);
    if (_use_names)
      nN = f_write(f + 0, r->data, r->data[0] + 1);
    else
      nN = 1;

#ifdef PACBIO
    size_t sz_read;
    if (t->output >= 0)
      sz_read = SZ_READ(r->read_length - lenCore);
    else
      sz_read = SZ_READ(r->read_length);
#endif
    nR = f_write(f + 1, r->data + nN, sz_read);
    // write the end marker
    f_write(f + 1, &(r->end), sz_meta);

    nQ = f_write(f + 2, r->data + nN + nR, r->of - nN - nR);

    if (_use_second_file) {
#ifdef PACBIO
      sz_read = SZ_READ(r->read_length_2);
      nR2 = f_write(f + 4, r->data + r->of, sz_read);
#else
      nR2 = f_write(f + 4, r->data + r->of, SZ_READ(read_length[1]));
#endif
      nQ2 = f_write(f + 5, r->data + r->of + nR2, r->sz - (r->of + nR2));
    }

#ifdef PACBIO
    nL = f_write(f + 6, &r->read_length, sizeof(uint32_t));
    if (_use_second_file)
      nL += f_write(f + 6, &r->read_length_2, sizeof(uint32_t));
    tL += nL;
#endif

    if (_use_names)
      tN += nN;
    tQ += nQ;
    tR += nR + sz_meta;
    tR2 += nR2;
    tQ2 += nQ2;
  }
  bin_free(&(t->bin));

  // META
  if (t->output == -1) {
    int32_t x = MAXBIN - 1;
    f_write(f + 3, &x, sizeof(int32_t)); // which core
    f_write(f + 3, &x, sizeof(int32_t)); // which core?
  } else {
    f_write(f + 3, &(t->id), sizeof(int32_t));     // which core?
    f_write(f + 3, &(t->output), sizeof(int32_t)); // which core?
  }

  f_write(f + 3, &tN, sizeof(int64_t));
  f_write(f + 3, &tR, sizeof(int64_t));
  f_write(f + 3, &tQ, sizeof(int64_t));
  if (_use_second_file) {
    f_write(f + 3, &tR2, sizeof(int64_t));
    f_write(f + 3, &tQ2, sizeof(int64_t));
  }
#ifdef PACBIO
  f_write(f + 3, &tL, sizeof(int64_t));
#endif
}

/* trie_queue - simple queue for aho_trie (for bfs traversal) */
typedef struct trie_queue trie_queue;
struct trie_queue {
  aho_trie **data;
  int sz, start, end, size;
};

/* queue initialization */
void trie_queue_init(trie_queue *q, int sz) {
  q->data = (aho_trie **)mallox(sz * sizeof(struct aho_trie *));
  q->sz = sz;
  q->start = q->end = q->size = 0;
}

/* queue push operation */
void trie_queue_push(trie_queue *q, aho_trie *t) {
  q->data[q->end] = t;
  q->size++;
  q->end = (q->end + 1) % q->sz;
}

/* and queue pop operation */
aho_trie *trie_queue_pop(trie_queue *q) {
  aho_trie *t = q->data[q->start];
  q->start = (q->start + 1) % q->sz;
  q->size--;
  return t;
}

/* safely dispose queue */
void trie_queue_free(trie_queue *q) {
  frex(q->data, q->sz * sizeof(struct aho_trie *));
}

char *visited = 0;   /* visited flags for bfs traversal */
int nodes_count = 0; /* total number of nodes in tree */
int unbucketed = 0;  /* number of unbucketed reads */

/* initialize aho trie */
void aho_trie_init(aho_trie *t) {
  t->fail = 0;
  t->output = -1;
  t->bin_size = 0;
  t->next_to_output = 0;
  t->level = 0;
  t->id = 0;
  t->child[0] = t->child[1] = t->child[2] = t->child[3] = 0;
  bin_init(&(t->bin));
}

/* bucket one read */
bin_node *aho_trie_bucket(aho_trie *t, read_data *d) {
  bin_node *n = (bin_node *)getpool(sizeof(struct bin_node));
  n->data = *d;
  n->data.data = (uint8_t *)getpool(d->sz * sizeof(uint8_t));
  n->next = 0;

  if (t->bin.size == 0)
    t->bin.first = t->bin.last = n;
  else {
    t->bin.last->next = n;
    t->bin.last = n;
  }
  t->bin.size++;
  t->bin_size++;

  //	memcpy (n->data.data, d->data, d->sz);
  return n;
}

/* insert core string into trie */
void pattern_insert(char *c, aho_trie *n, int level, int id) {
  if (*c != 0 && *c != '\n') {
    char cx = getval(*c);
    if (n->child[cx] == 0) {
      n->child[cx] = (aho_trie *)mallox(sizeof(struct aho_trie));
      aho_trie_init(n->child[cx]);
      n->child[cx]->level = level + 1;
      nodes_count++;
    }
    pattern_insert(c + 1, n->child[cx], level + 1, id);
  } else {
    n->output = id;
    //		n-> = 1;
  }
}

/* initialize aho automaton */
void prepare_aho_automata(aho_trie *root) {
  trie_queue q;
  trie_queue_init(&q, nodes_count + 1);
  root->fail = root;

  for (int i = 0; i < 4; i++)
    if (root->child[i]) {
      root->child[i]->fail = root;
      trie_queue_push(&q, root->child[i]);
    }

  int traversed = 0;
  while (q.size) {
    aho_trie *cur = trie_queue_pop(&q);
    for (int i = 0; i < 4; i++) {
      aho_trie *t = cur->child[i];
      if (t) {
        aho_trie *f = cur->fail;
        while (f != root && !f->child[i])
          f = f->fail;
        t->fail = (f->child[i] ? f->child[i] : root);
        trie_queue_push(&q, t);
        t->next_to_output =
            (t->fail->output >= 0) ? t->fail : t->fail->next_to_output;
      }
    }
    cur->id = ++traversed;
  }
  q.start = q.end = q.size = 0;
  trie_queue_push(&q, root);
  while (q.size) {
    aho_trie *cur = trie_queue_pop(&q);
    for (int i = 0; i < 4; i++) {
      if (cur->child[i]) {
        trie_queue_push(&q, cur->child[i]);
      }
      aho_trie *c = cur;
      while (c != root && !c->child[i])
        c = c->fail;
      cur->child[i] = (c->child[i] ? c->child[i] : root);
      c = cur;
      while (c && c->output == -1)
        c = c->next_to_output;
      cur->next_to_output = c;
    }
  }

  trie_queue_free(&q);
  visited = (char *)mallox((nodes_count + 1) * sizeof(char));

  LOG("Allocating pool... ");
  _pool_pos = 0;
  _pool = (uint8_t *)mallox(_max_bucket_set_size + 100 * 5 * MAXLINE);
  LOG("OK!\n");
}

/* read core strings and create trie and aho automaton */
extern char _binary_patterns_bin_start;
extern char _binary_patterns_bin_end;
extern char _binary_patterns_bin_size;
aho_trie *read_patterns() {
  aho_trie *root = (aho_trie *)mallox(sizeof(struct aho_trie));
  aho_trie_init(root);
  char line[MAXLINE];
  char alphabet[] = "ACGT";

  patterns = (char **)mallox(5000000 * sizeof(char *));

  char *data = &_binary_patterns_bin_start;
  int pos = 0;
  int64_t size =
      (int64_t)(&_binary_patterns_bin_end - &_binary_patterns_bin_start);
  while (1) {
    int16_t ln;
    if (pos == size)
      break;
    memcpy(&ln, data + pos, sizeof(int16_t));
    pos += sizeof(int16_t);

    int32_t cnt;
    memcpy(&cnt, data + pos, sizeof(int32_t));
    pos += sizeof(int32_t);

    int sz = ln / 4 + (ln % 4 != 0);
    int64_t x;
    for (int i = 0; i < cnt; i++) {
      int nl = 0;
      memcpy(&x, data + pos, sz);
      assert(sz <= 8);
      pos += sz;

      patterns[pattern_c] = (char *)mallox(ln + 1);
      for (int j = ln - 1; j >= 0; j--)
        patterns[pattern_c][nl++] = alphabet[(x >> (2 * j)) & 3];
      patterns[pattern_c][nl++] = 0;
      // printf("%d %s\n", pattern_c, patterns[pattern_c]);
      pattern_insert(patterns[pattern_c], root, 0, pattern_c);
      pattern_c++;
    }
  }

  if (_decompress)
    return 0;
  else {
    prepare_aho_automata(root);
    return root;
  }
}

aho_trie *read_patterns_from_file(const char *path) {
  aho_trie *root = (aho_trie *)mallox(sizeof(struct aho_trie));

  LOG("Pattern file %s ... ", path);

  aho_trie_init(root);
  patterns = (char **)mallox(10000000 * sizeof(char *));
  memset(patterns, 0, sizeof(char *) * 10000000);

  FILE *f = fopen(path, "r");
  while (fscanf(f, "%ms", &(patterns[pattern_c])) != EOF) {
    //		int len = strlen(pt);
    //		printf("-- %d --\n", len);
    pattern_insert(patterns[pattern_c], root, 0, pattern_c);
    pattern_c++;
  }

  LOG("read %'d patterns\n", pattern_c);

  /*	strcpy(patterns[pattern_c], "AAAAAAAAAAAAAAAA");
          pattern_insert (patterns[pattern_c], root, 0, pattern_c); pattern_c++;
          strcpy(patterns[pattern_c], "CCCCCCCCCCCCCCCC");
          pattern_insert (patterns[pattern_c], root, 0, pattern_c); pattern_c++;
          strcpy(patterns[pattern_c], "GGGGGGGGGGGGGGGG");
          pattern_insert (patterns[pattern_c], root, 0, pattern_c); pattern_c++;
          strcpy(patterns[pattern_c], "TTTTTTTTTTTTTTTT");
          pattern_insert (patterns[pattern_c], root, 0, pattern_c); pattern_c++;
  */

  prepare_aho_automata(root);
  return root;
}

/* search for core strings in the read */
int aho_search(char *text, aho_trie *root, aho_trie **bucket) {
  aho_trie *cur = root, *largest = 0, *x;
  int bestpos = -1;
  for (int i = 0; text[i] != '\n'; i++) {
    cur = cur->child[getval(text[i])];
    if (cur->next_to_output) {
      x = cur->next_to_output;
      if (!largest || largest->level < x->level ||
          (largest->level == x->level && largest->bin_size < x->bin_size)) {
        bestpos = i;
        largest = x;
      }
    }
  }
  *bucket = (largest ? largest : root);
  return bestpos;
}

/* 2/8 read encoding */
int output_read(char *line, uint8_t *dest, int n, int l) {
  int bc = 0;

  uint8_t ca = 0, cc = 0;
  for (int i = n + l; line[i] != '\n'; i++) {
    ca = (ca << 2) | (line[i] ? getval(line[i]) : 0);
    cc++;
    if (cc == 4) {
      dest[bc++] = ca;
      cc = 0;
    }
  }
  for (int i = 0; i < n; i++) {
    ca = (ca << 2) | (line[i] ? getval(line[i]) : 0);
    cc++;
    if (cc == 4) {
      dest[bc++] = ca;
      cc = 0;
    }
  }
  if (cc) {
    while (cc != 4) {
      ca <<= 2;
      cc++;
    }
    dest[bc++] = ca;
  }

  return bc;
}

/* flush (output) all contents of trie in the file
 * unbucketed reads are flushed at the end
 * see dump_trie for more information */
void aho_output(aho_trie *root, buffered_file *f) {
  memset(visited, 0, sizeof(char) * (nodes_count + 1));
  visited[root->id] = 1;

  trie_queue q;
  trie_queue_init(&q, nodes_count + 1);

  for (int i = 0; i < 4; i++) {
    trie_queue_push(&q, root->child[i]);
    visited[root->child[i]->id] = 1;
  }

  while (q.size) {
    aho_trie *cur = trie_queue_pop(&q);
    assert(cur->id <= nodes_count);
    if (cur->bin.size) {
      bin_prepare(cur);
      bin_dump(cur, f);
    }
    for (int i = 0; i < 4; i++)
      if (cur->child[i] && !visited[cur->child[i]->id]) {
        trie_queue_push(&q, cur->child[i]);
        visited[cur->child[i]->id] = 1;
      }
  }
  if (root->bin.size) {
    unbucketed += root->bin.size;
    bin_prepare(root);
    bin_dump(root, f);
  }

  clrpool();
  trie_queue_free(&q);
}

/* get number fo unbucketed reads */
int unbuck(void) { return unbucketed; }

/* safely dispose trie */
void aho_trie_free(aho_trie *t) {
  trie_queue q;
  trie_queue_init(&q, nodes_count + 1);
  trie_queue_push(&q, t);

  aho_trie **p =
      (aho_trie **)mallox((nodes_count + 1) * sizeof(struct aho_trie *));
  memset(p, 0, (1 + nodes_count) * sizeof(struct aho_trie *));
  p[t->id] = t;

  while (q.size) {
    aho_trie *cur = trie_queue_pop(&q);
    for (int i = 0; i < 4; i++)
      if (cur->child[i] && !p[cur->child[i]->id]) {
        trie_queue_push(&q, cur->child[i]);
        p[cur->child[i]->id] = cur->child[i];
      }
  }
  for (int i = 0; i < nodes_count; i++)
    if (p[i]) {
      bin_free(&(p[i]->bin));
      if (p[i]->output >= 0 && patterns[p[i]->output]) {
        frex(patterns[p[i]->output], 15);
        patterns[p[i]->output] = 0;
      }
      frex(p[i], sizeof(struct aho_trie));
    }
  frex(visited, nodes_count);
  frex(p, nodes_count * sizeof(struct aho_trie *));
  frex(patterns, sizeof(char *) * 5000000);
}

/** fast radix sorting for dna strings  *
  * params: pos   - in-string position (default 0)
  *         start - start position of partition
  *			size  - partition size
  *			v     - vector to be sorted (vector of pointers to string
  *references!)
  */
/*#define _POS(n,x,i) \
        (( (n[x]->data.data[ n[x]->data.data[0] + 1 + ((i) / 4)]) >> ((3 - ((i)
   % 4)) * 2) ) & 3)*/
#define _POSX(n, x, i)                                                         \
  (((n[x]->data.data[n[x]->data.data[0] + 1 + ((i) / 4)]) >>                   \
    ((3 - ((i) % 4)) * 2)) &                                                   \
   3)
#ifdef PACBIO
#define _POS(n, x, i)                                                          \
  ((i + n[x]->data.end < n[x]->data.read_length)                               \
       ? _POSX(n, x, i)                                                        \
       : 0) // n,x,iarray, pos in array, pos in item
#else
#define _POS(n, x, i)                                                          \
  ((i + n[x]->data.end < read_length[0]) ? _POSX(n, x, i) : 0)
#endif

bin_node **_nodes = 0, **_temp = 0;
int _bin_sz = 0, _allocd = 0, _limit = 0;
void _radix_sort(int pos, int64_t start, int64_t size) {
  if (size <= 1) /* trivial case: subarray of size 1 */
    return;
  if (pos >= _limit) /* trivial case: exceeded string length */
    return;

  /* obtain character stats in count array */
  int count[4] = {0};
  for (int i = start; i < start + size; i++) {
    count[_POS(_nodes, i, pos)]++;
    _temp[i] = _nodes[i]; /* keep original ordering of vector safe */
  }

  /* get cumulative sums for count array */
  int cumulative[5] = {0};
  for (int i = 1; i < 5; i++)
    cumulative[i] = cumulative[i - 1] + count[i - 1];

  /* partition the original array according to the starting character */
  for (int i = start; i < start + size; i++) {
    char c = _POS(_temp, i, pos);
    _nodes[start + cumulative[c]] = _temp[i];
    cumulative[c]++; /* here, "fake" update cumulative sum in order to keep
                        track of partition size */
  }

  /* update cumulatove sums */
  cumulative[0] = 0;
  for (int i = 1; i < 5; i++)
    cumulative[i] = cumulative[i - 1] + count[i - 1];

  /* sort partitions */
  for (int i = 0; i < 4; i++)
    _radix_sort(pos + 1, start + cumulative[i],
                cumulative[i + 1] - cumulative[i]);
}

void bin_prepare(aho_trie *t) {
  bin *b = &(t->bin);
  // LOG("core %10d    sz %10d orig len %2d\n",t->id,b->size,t->level);

  if (_allocd == 0 && b->size < 10000) {
    _allocd = 10000;
    _temp = (bin_node **)mallox(_allocd * sizeof(bin_node *));
    _nodes = (bin_node **)mallox(_allocd * sizeof(bin_node *));
  } else if (b->size > _allocd) {
    frex(_temp, _allocd * sizeof(bin_node *));
    frex(_nodes, _allocd * sizeof(bin_node *));
    _temp = (bin_node **)mallox(b->size * sizeof(bin_node *));
    _nodes = (bin_node **)mallox(b->size * sizeof(bin_node *));
    _allocd = b->size;
  }
  _bin_sz = 0;

  _limit = 0;
  for (bin_node *n = b->first; n; n = n->next) {
    _nodes[_bin_sz++] = n;
#ifdef PACBIO
    _limit = MAX(_limit, n->data.read_length - t->level);
#endif
  }
#ifndef PACBIO
  _limit = read_length[0] - t->level;
#endif
  _radix_sort(0, 0, _bin_sz);

  b->first = _nodes[0];
  b->last = _nodes[_bin_sz - 1];
  for (int i = 0; i < _bin_sz - 1; i++)
    _nodes[i]->next = _nodes[i + 1];
  _nodes[_bin_sz - 1]->next = 0;
}
