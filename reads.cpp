/// 786

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

#include "reads.h"
#include "const.h"

char **patterns;
static int pattern_c = 0;

void bin_prepare (aho_trie *t);

/* initialize bin */
void bin_init (bin *b) {
	b->first = b->last = 0;
	b->size = 0;
}

/* insert read_data into bin (linked list) */
void bin_insert (bin *b, read_data *d) {
	bin_node *n  = (bin_node*) mallox (sizeof (struct bin_node));
	n->data = *d;
	n->data.data = (uint8_t*) mallox (d->sz * sizeof (uint8_t));
	memcpy (n->data.data, d->data, d->sz);
	n->next = 0;

	#pragma omp critical 
	{
		if (b->size == 0)
			b->first = b->last = n;
		else {
			b->last->next = n;
			b->last = n;
		}
		b->size++;
	}
}

/* safely dispose bin */
void bin_free (bin *b) {
	if (b->size) for (bin_node *n = b->first; n; ) {
		bin_node *x = n->next;

		frex(n->data.data, sizeof(uint8_t)*n->data.sz);
		frex(n, sizeof(struct bin_node));

		n = x;
	}
	b->first = b->last = 0;
	b->size = 0;
}

/* dump bin contents to file 
 * f must contain 5 pointers, for names, reads and qualities files
 * (and reads and qualities for paired end) */
void bin_dump (aho_trie *t, buffered_file *f) {
	int64_t tN = 0, tQ = 0, tR = 0, tR2 = 0, tQ2 = 0;
	int64_t nN = 0, nR = 0, nQ = 0, nR2 = 0, nQ2 = 0;
	
	int sz_read;
	if (t->output >= 0) sz_read = SZ_READ( read_length[0] - strlen(patterns[t->output]) );
	else sz_read = SZ_READ (read_length[0]);

	int sz_meta = 1;
	if (read_length[0] > 255) sz_meta = 2;

	for (bin_node *n = t->bin.first; n; n = n->next) {
		read_data *r = &(n->data);
		if (_use_names)
			nN = f_write (f + 0, r->data, r->data[0] + 1);
		nR = f_write (f + 1, r->data + nN, sz_read);
		// write the end marker
		f_write(f+1, &(r->end), sz_meta);

		nQ = f_write (f + 2, r->data + nN + nR, r->of - nN - nR);

		if (_use_second_file) {
			nR2 = f_write (f + 4, r->data + r->of, SZ_READ (read_length[1]));
			nQ2 = f_write (f + 5, r->data + r->of + nR2, r->sz - (r->of + nR2));
		}
		tN += nN; tQ += nQ; tR += nR + sz_meta; tR2 += nR2; tQ2 += nQ2;
	}
	bin_free (&(t->bin));

	
	if (t->output == -1) {
		int32_t x = MAXBIN - 1;
		f_write (f + 3, &x, sizeof(int32_t)); // which core
		f_write (f + 3, &x, sizeof(int32_t)); // which core?
	}
	else {
		f_write (f + 3, &(t->id),     sizeof(int32_t)); // which core?
		f_write (f + 3, &(t->output), sizeof(int32_t)); // which core?
	}

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
	q->data = (aho_trie**) mallox (sz * sizeof (struct aho_trie*));
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
	frex(q->data, q->sz*sizeof(struct aho_trie*));
}

char *visited    = 0;  /* visited flags for bfs traversal */
int  nodes_count = 0;  /* total number of nodes in tree */
int  unbucketed  = 0;  /* number of unbucketed reads */

/* initialize aho trie */
void aho_trie_init (aho_trie *t) {
	t->fail = 0;
	t->output = -1;
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
	#pragma omp critical 
	{
		t->bin_size++;
	}
}

/* insert core string into trie */
void pattern_insert (char *c, aho_trie *n, int level, int id) {
	if (*c != 0 && *c != '\n') {
		printf("%c", *c);
		char cx = getval (*c);
		if (n->child[cx] == 0) {
			n->child[cx] = (aho_trie*) mallox (sizeof (struct aho_trie));
			aho_trie_init (n->child[cx]);
			n->child[cx]->level = level + 1;
			nodes_count++;
		}
		pattern_insert (c + 1, n->child[cx], level + 1, id); 
	}
	else {
		n->output = id;
		printf("\n");
//		n-> = 1;
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
				t->next_to_output = (t->fail->output >= 0) ? t->fail : t->fail->next_to_output;
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
			while (c && c->output == -1)
				c = c->next_to_output;
			cur->next_to_output = c;
		}
	}

	trie_queue_free (&q);
	visited = (char*) mallox ((nodes_count + 1) * sizeof (char));
}

/* read core strings and create trie and aho automaton */
extern char _binary_patterns_bin_start;
extern char _binary_patterns_bin_end;
extern char _binary_patterns_bin_size;
aho_trie *read_patterns () {
	aho_trie *root = (aho_trie*) mallox (sizeof (struct aho_trie));
	aho_trie_init (root);
	char line[MAXLINE];
	char alphabet[] = "ACGT";

	patterns = (char**) mallox(5000000*sizeof(char*));
	for(int i=0;i<5000000;i++)
		patterns[i]=(char*) mallox(15);


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
			assert(sz<=8);
			pos += sz;
			
//			patterns[pattern_c] = mallox( ln );
			for(int j = ln - 1; j >= 0; j--)
				patterns[pattern_c][nl++] = alphabet[ (x >> (2 * j)) & 3 ];
			patterns[pattern_c][nl++] = 0;
			pattern_insert (patterns[pattern_c], root, 0, pattern_c);
			pattern_c++;
		}
	}

	printf("Read %d patterns\n", pattern_c);

	prepare_aho_automata (root);
	return root;
}

aho_trie *read_patterns_from_file (const char *path) {
	aho_trie *root = (aho_trie*) mallox (sizeof (struct aho_trie));
	aho_trie_init (root);
	char line[MAXLINE];
	char alphabet[] = "ACGT";
	patterns = (char**) mallox(5000000*sizeof(char*));
	for(int i=0;i<5000000;i++)
		patterns[i]=(char*) mallox(15);


	int len;
	uint64_t val;

	FILE *f = fopen(path, "r");
	LOG("reading patterns from %s\n", path);
	while (fscanf(f, "%d %llu", &len, &val) != EOF) {
		if (len <= 4) {
			LOG("\tskipping (%llu,%d)\n", val, len);
			continue;
		}
		for (int i = 0; i < len; i++) {
			patterns[pattern_c][len - i - 1] = alphabet[val & 3];
			val >>= 2;
		}
		patterns[pattern_c][len] = 0;
		pattern_insert (patterns[pattern_c], root, 0, pattern_c);
		pattern_c++;
	}
	printf("Read %d patterns\n", pattern_c);
	
	prepare_aho_automata (root);
	return root;
}

/* search for core strings in the read */
int aho_search (char *text, aho_trie *root, aho_trie **bucket) {
	aho_trie *cur = root, *largest = 0, *x;

	int bestpos = -1;
	for (int i = 0; text[i] != 0 && text[i] != '\n'; i++) {
		cur = cur->child[getval (text[i])];
		if (cur->next_to_output) {
			x = cur->next_to_output;
			if (!largest ||
					largest->level < x->level ||
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
int output_read (char *line, uint8_t *dest, int n, int l) {
	int bc = 0;

	uint8_t ca = 0, cc = 0;
	for (int i = n + l; line[i]; i++) {
		ca = (ca << 2) | (line[i] ? getval(line[i]): 0);
		cc++;
		if (cc == 4) {
			dest[bc++] = ca;
			cc = 0;
		}
	}
	for (int i = 0; i < n; i++) {
		ca = (ca << 2) | (line[i] ? getval(line[i]): 0);
		cc++;
		if (cc == 4) {
			dest[bc++] = ca;
			cc = 0;
		}
	}
	if (cc) {
		while (cc != 4) { ca <<= 2; cc++; }
		dest[bc++] = ca;
	}

	return bc;
}


/* flush (output) all contents of trie in the file
 * unbucketed reads are flushed at the end
 * see dump_trie for more information */
void aho_output (aho_trie *root, buffered_file *f) { 
	memset (visited, 0, sizeof(char) * (nodes_count + 1));
	visited[root->id] = 1;

	trie_queue q;
	trie_queue_init (&q, nodes_count + 1);


	for (int i = 0; i < 4; i++) {
		trie_queue_push (&q, root->child[i]);
		visited[root->child[i]->id] = 1;
	}

	while (q.size) {
		aho_trie *cur = trie_queue_pop (&q);
		assert(cur->id <= nodes_count);
		if (cur->bin.size) {
			bin_prepare (cur);
			bin_dump (cur, f);
		}
		for (int i = 0; i < 4; i++) 
			if (cur->child[i] && !visited[cur->child[i]->id]) {
				trie_queue_push (&q, cur->child[i]);
				visited[cur->child[i]->id] = 1;
			}
	}
	if (root->bin.size) {
		unbucketed += root->bin.size;
		bin_prepare (root);
		bin_dump (root, f);
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

	aho_trie **p = (aho_trie**) mallox ((nodes_count+1) * sizeof (struct aho_trie*));
	memset(p, 0, (1+nodes_count)*sizeof(struct aho_trie*));
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
		if (p[i]->output >= 0 && patterns[p[i]->output]) {
			frex(patterns[p[i]->output], 15);
			patterns[p[i]->output] = 0;
		}
		frex(p[i], sizeof(struct aho_trie));
	}
	frex(visited, nodes_count);
	frex(p, nodes_count * sizeof(struct aho_trie*));
	frex(patterns, sizeof(char*)*5000000);
}


/** fast radix sorting for dna strings  *
  * params: pos   - in-string position (default 0)
  *         start - start position of partition
  *			size  - partition size
  *			v     - vector to be sorted (vector of pointers to string references!)
  */
/*#define _POS(n,x,i) \
	(( (n[x]->data.data[ n[x]->data.data[0] + 1 + ((i) / 4)]) >> ((3 - ((i) % 4)) * 2) ) & 3)*/
#define _POSX(n,x,i) \
	(( (n[x]->data.data[ n[x]->data.data[0] + 1 + ((i) / 4)]) >> ((3 - ((i) % 4)) * 2) ) & 3)
#define _POS(n,x,i) \
	((i + n[x]->data.end < read_length[0]) ? _POSX(n,x,i) : 0)

bin_node **_nodes = 0,
			**_temp  = 0;
int _bin_sz = 0, 
	 _allocd = 0, 
	 _limit  = 0;
void _radix_sort (int pos, int64_t start, int64_t size) {
	if (size <= 1)  /* trivial case: subarray of size 1 */
		return;
	if (pos >= _limit) /* trivial case: exceeded string length */
		return;

	/* obtain character stats in count array */
	int count[4] = {0}; 
	for (int i = start; i < start + size; i++) {
		count[ _POS(_nodes, i, pos) ]++;
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
		cumulative[c]++; /* here, "fake" update cumulative sum in order to keep track of partition size */
	}

	/* update cumulatove sums */
	cumulative[0] = 0;
	for (int i = 1; i < 5; i++)
		cumulative[i] = cumulative[i - 1] + count[i - 1];

	/* sort partitions */
	for (int i = 0; i < 4; i++)
		_radix_sort(pos + 1, start + cumulative[i], cumulative[i + 1] - cumulative[i]);
}

void bin_prepare (aho_trie *t) {
	bin *b = &(t->bin);
	//LOG("core %10d    sz %10d orig len %2d\n",t->id,b->size,t->level);

	if (_allocd == 0 && b->size < 10000) {
		_allocd = 10000;
		_temp = (bin_node**) mallox (_allocd * sizeof(bin_node*));
		_nodes = (bin_node**) mallox (_allocd * sizeof(bin_node*));
	}
	else if (b->size > _allocd) {
		frex(_temp, _allocd*sizeof(bin_node*));
		frex(_nodes, _allocd*sizeof(bin_node*));
		_temp = (bin_node**) mallox (b->size * sizeof(bin_node*));
		_nodes = (bin_node**) mallox (b->size * sizeof(bin_node*));
		_allocd=b->size;
	}
	_bin_sz = 0;
	for (bin_node *n = b->first; n; n = n->next) 
		_nodes[_bin_sz++] = n;
	_limit=read_length[0] - t->level;
	_radix_sort (0, 0, _bin_sz);
	
	b->first = _nodes[0];
	b->last = _nodes[_bin_sz - 1];
	for (int i = 0; i < _bin_sz - 1; i++)
		_nodes[i]->next = _nodes[i + 1];
	_nodes[_bin_sz - 1]->next = 0;
}

