/// 786

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

#include "reads.h"
#include "const.h"

/* initialize bin */
void bin_init (bin *b) {
	b->first = b->last = 0;
	b->size = 0;
}

/* insert read_data into bin (linked list) */
void bin_insert (bin *b, read_data *d) {
	bin_node *n  = malloc (sizeof (struct bin_node));
	n->data.sz   = d->sz;
	n->next      = 0;
	n->data.of   = d->of;
	n->data.data = malloc (d->sz * sizeof (uint8_t));
	memcpy (n->data.data, d->data, d->sz);
	if (b->size == 0)
		b->first = b->last = n;
	else {
		b->last->next = n;
		b->last = n;
	}
	b->size++;
}

/* safely dispose bin */
void bin_free (bin *b) {
	if (b->size) for (bin_node *n = b->first; n; ) {
		bin_node *x = n->next;
		free (n);
		n = x;
	}
	b->first = b->last = 0;
	b->size = 0;
}

/* dump bin contents to file 
 * f must contain 5 pointers, for names, reads and qualities files
 * (and reads and qualities for paired end) */
void bin_dump (bin *b, buffered_file *f, int32_t id) {
	int64_t tN = 0, tQ = 0, tR = 0, tR2 = 0, tQ2 = 0;
	int64_t nN = 0, nR = 0, nQ = 0, nR2 = 0, nQ2 = 0;
	
	for (bin_node *n = b->first; n; n = n->next) {
		read_data *r = &(n->data);
		if (_use_names)
			nN = f_write (f + 0, r->data, r->data[0] + 1);
		nR = f_write (f + 1, r->data + nN, SZ_READ (read_length[0]));
		nQ = f_write (f + 2, r->data + nN + nR, r->of - nN - nR);

		if (_use_second_file) {
			nR2 = f_write (f + 4, r->data + r->of, SZ_READ (read_length[1]));
			nQ2 = f_write (f + 5, r->data + r->of + nR2, r->sz - (r->of + nR2));
		}
		tN += nN; tQ += nQ; tR += nR; tR2 += nR2; tQ2 += nQ2;
		free (r->data);
	}
	bin_free (b);

	f_write (f + 3, &id, sizeof(int32_t));
	f_write (f + 3, &tN, sizeof(int64_t));
	f_write (f + 3, &tR, sizeof(int64_t));
	f_write (f + 3, &tQ, sizeof(int64_t));
	if (_use_second_file) {
		f_write (f + 3, &tR2, sizeof(int64_t));
		f_write (f + 3, &tQ2, sizeof(int64_t));
	}
}

/* trie_queue - simple queue for aho_trie (for bfs traversal) */
typedef struct trie_queue trie_queue;
struct trie_queue {
	aho_trie **data;
	int sz, start, end, size;
};

/* queue initialization */
void trie_queue_init (trie_queue *q, int sz) {
	q->data = malloc (sz * sizeof (struct aho_trie*));
	q->sz = sz;
	q->start = q->end = q->size = 0;
}

/* queue push operation */
void trie_queue_push (trie_queue *q, aho_trie *t) {
	q->data[q->end] = t;
	q->size++;
	q->end = (q->end + 1) % q->sz;
}

/* and queue pop operation */
aho_trie *trie_queue_pop (trie_queue *q) {
	aho_trie *t = q->data[q->start];
	q->start = (q->start + 1) % q->sz;
	q->size--;
	return t;
}

/* safely dispose queue */
void trie_queue_free (trie_queue *q) {
	free (q->data);
}

char *visited    = 0;  /* visited flags for bfs traversal */
int  nodes_count = 0;  /* total number of nodes in tree */
int  unbucketed  = 0;  /* number of unbucketed reads */

/* initialize aho trie */
void aho_trie_init (aho_trie *t) {
	t->fail = 0;
	t->output = 0;
	t->bin_size = 0;
	t->next_to_output = 0;
	t->level = 0;
	t->id = 0;
	t->child[0] = t->child[1] = t->child[2] = t->child[3] = 0;
	bin_init (&(t->bin));
}

/* bucket one read */
void aho_trie_bucket (aho_trie *t, read_data *d) {
	bin_insert (&(t->bin), d);
	t->bin_size++;
}
/* insert core string into trie */
void pattern_insert (char *c, aho_trie *n, int level) {
	if (*c != 0 && *c != '\n') {
		*c = getval (*c);
		if (n->child[*c] == 0) {
			n->child[*c] = malloc (sizeof (struct aho_trie));
			aho_trie_init (n->child[*c]);
			n->child[*c]->level = level + 1;
			nodes_count++;
		}
		pattern_insert (c + 1, n->child[*c], level + 1); 
	}
	else {
		n->output = 1;
	}
}

/* initialize aho automaton */
void prepare_aho_automata (aho_trie *root) {
	trie_queue q;
	trie_queue_init (&q, nodes_count + 1);
	root->fail = root;

	for (int i = 0; i < 4; i++) 
		if (root->child[i]) {
			root->child[i]->fail = root;
			trie_queue_push (&q, root->child[i]);
		}

	int traversed = 0;
	while (q.size) {
		aho_trie *cur = trie_queue_pop (&q);
		for (int i = 0; i < 4; i++) {
			aho_trie *t = cur->child[i];
			if (t) {
				aho_trie *f = cur->fail;
				while (f != root && !f->child[i])
					f = f->fail;
				t->fail = (f->child[i] ? f->child[i] : root);
				trie_queue_push (&q, t);
				t->next_to_output = t->fail->output ? t->fail : t->fail->next_to_output;
			}
		}
		cur->id = ++traversed;
	}
	q.start = q.end = q.size = 0;
	trie_queue_push (&q, root);
	while (q.size) {
		aho_trie *cur = trie_queue_pop (&q);
		for (int i = 0; i < 4; i++) {
			if (cur->child[i]) {
				trie_queue_push (&q, cur->child[i]);
			}
			aho_trie *c = cur;
			while (c != root && !c->child[i]) 
				c = c->fail;
			cur->child[i] = (c->child[i] ? c->child[i] : root);
			c = cur;
			while (c && !c->output)
				c = c->next_to_output;
			cur->next_to_output = c;
		}
	}

	trie_queue_free (&q);
	visited = malloc (nodes_count * sizeof (char));
}

/* read core strings and create trie and aho automaton */
extern char _binary_patterns_bin_start;
extern char _binary_patterns_bin_end;
extern char _binary_patterns_bin_size;
aho_trie *read_patterns () {
	aho_trie *root = malloc (sizeof (struct aho_trie));
	aho_trie_init (root);
	char line[MAXLINE];
	char alphabet[] = "ACGT";

	char *data = &_binary_patterns_bin_start;
	int pos = 0;
	int64_t size = (int64_t)(&_binary_patterns_bin_size);
	while (1) {
		int16_t ln;
		if (pos == size)
			break;
		memcpy (&ln, data + pos, sizeof(int16_t));
		pos += sizeof(int16_t);
		
		int32_t cnt;
		memcpy (&cnt, data + pos, sizeof(int32_t));
		pos += sizeof(int32_t);
		
		int sz = ln / 4 + (ln % 4 != 0);
		int64_t x;
		for(int i = 0; i < cnt; i++) {
			int nl = 0;
			memcpy (&x, data + pos, sz);
			pos += sz;
			for(int j = ln - 1; j >= 0; j--)
				line[nl++] = alphabet[ (x >> (2 * j)) & 3 ];
			line[nl++] = 0;
			pattern_insert (line, root, 0);
		}
	}
	prepare_aho_automata (root);
	return root;
}

/* search for core strings in the read */
void aho_search (char *text, aho_trie *root, read_data *pos) {
	aho_trie *cur = root, *largest = 0, *x;
	for (int i = 0; text[i] != 0 && text[i] != '\n'; i++) {
		if (text[i] == 'N')
			cur = root;
		else {
			cur = cur->child[getval (text[i])];
			if (cur->next_to_output) {
				x = cur->next_to_output;
				if (!largest || 
						largest->bin_size < x->bin_size || 
						(largest->bin_size == x->bin_size && largest->level < x->level)) 
					largest = x;
			}
		} 
	}
	aho_trie_bucket (largest ? largest : root, pos);
}

/* 2/8 read encoding */
int output_read (char *line, uint8_t *dest) {
	int bc = 0;
	int l = strlen (line);
	for (int i = 0; i < l; i += 4) {
		uint8_t ca = 0;
		for (int j = 0; j < 4; j++)
			ca = (ca << 2) | ((i+j<l)?getval(line[i+j]):0);
		dest[bc++] = ca;
	}

	return bc;
}

/* flush (output) all contents of trie in the file
 * unbucketed reads are flushed at the end
 * see dump_trie for more information */
void aho_output (aho_trie *root, buffered_file *f) { 
	memset (visited, 0, sizeof(char) * nodes_count);
	visited[root->id] = 1;

	trie_queue q;
	trie_queue_init (&q, nodes_count + 1);

	for (int i = 0; i < 4; i++) {
		trie_queue_push (&q, root->child[i]);
		visited[root->child[i]->id] = 1;
	}

	while (q.size) {
		aho_trie *cur = trie_queue_pop (&q);
		if (cur->bin.size)
			bin_dump (&(cur->bin), f, cur->id);
		for (int i = 0; i < 4; i++) 
			if (cur->child[i] && !visited[cur->child[i]->id]) {
				trie_queue_push (&q, cur->child[i]);
				visited[cur->child[i]->id] = 1;
			}
	}
	if (root->bin.size) {
		unbucketed += root->bin.size;
		bin_dump (&(root->bin), f, MAXBIN - 1);
	}

	trie_queue_free (&q);
}

/* get number fo unbucketed reads */
int unbuck (void) {
	return unbucketed;
}

/* safely dispose trie */
void aho_trie_free (aho_trie *t) {
	trie_queue q;
	trie_queue_init (&q, nodes_count + 1);
	trie_queue_push (&q, t);

	aho_trie **p = calloc (nodes_count, sizeof (struct aho_trie*));
	p[t->id] = t;

	while (q.size) {
		aho_trie *cur = trie_queue_pop (&q);
		for (int i = 0; i < 4; i++)
			if (cur->child[i] && !p[cur->child[i]->id]) {
				trie_queue_push (&q, cur->child[i]);
				p[cur->child[i]->id] = cur->child[i];
			}
	}
	for (int i = 0; i < nodes_count; i++) if (p[i]) {
		bin_free (&(p[i]->bin));
		free (p[i]);
	}
	free (visited);
	free (p);
}
