/// 786

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

#include "buffio.h"
#include "const.h"
#include "decompress.h"
#include "arithmetic.h"
#include "reads.h"
#include "names.h"

void get_file_name (char *name, char c) {
	char *p;
	if (p = strstr (name, ".scalce")) {
		*(p + 7) = c;
	}
}

uint32_t WI(buffered_file *f) {
	uint32_t r=0,i=0;
	while(1) {
		uint8_t x; f_read(f,&x,1);

		r|=((x&0x7f)<<(i*7));
		if(x&0x80) i++;
		else break;
	}
	return r;
}

void decompress (const char *path, const char *out) {
	uint8_t buffer[MAXLINE];

	aho_trie *trie = read_patterns ();

	int32_t len[2]; int64_t nl[2], phredOff[2];

	int mode[2] = { IO_SYS, IO_SYS }; uint8_t c[2];

	buffered_file fR[2], fQ[2], fN[2];

	strncpy (buffer, path, MAXLINE);
	for (int i = 0; i < 1 + (_use_second_file|_interleave); i++) {
	/* file type detection */
		get_file_name ((char*)buffer, 'r');          
		FILE *ft = fopen ((char*)buffer, "rb");
		fread (c, 1, 2, ft);
		int mode = IO_SYS;
		if (c[0] == 0x1f && c[1] == 0x8b)
			mode = IO_GZIP;
		if (c[0] == 0x42 && c[1] == 0x5a)
			mode = IO_BZIP;
		fclose (ft);
		LOG("Paired-end %d, file format detected: %s\n",i,(mode==IO_GZIP?"gzip":(mode==IO_BZIP?"bzip":"plain")));
		
	/* initialization */
		f_init (fR + i, mode);
		f_init (fQ + i, IO_SYS);
		f_init (fN + i, mode);

	/* open read, name and quality file */
		get_file_name ((char*)buffer, 'r');
		f_open (fR + i, (char*)buffer, IO_READ);
		if (!f_alive (fR + i)) 
			ERROR("Cannot find read file %s!\n", buffer);
		get_file_name ((char*)buffer, 'q');
		f_open (fQ + i, (char*)buffer, IO_READ);
		if (!f_alive (fQ + i)) 
			ERROR("Cannot find quality file %s!\n", buffer);
		get_file_name ((char*)buffer, 'n');
		f_open (fN + i, (char*)buffer, IO_READ);
		if (!f_alive (fN + i) && _use_names) 
			ERROR("Cannot find name file %s! Use -n parameter if you want to skip name file lookup.\n", buffer);

	/* read magic */
		f_read (fR+i, buffer, 8);
		f_read (fQ+i, buffer, 8);
		f_read (fN+i, buffer, 8);
	/* read metadata */
		f_read (fR+i, len+i, sizeof(int32_t));
		f_read (fQ+i, phredOff+i, sizeof(int64_t));
		Z=i; a_init_decoder (fQ + i);

		if (!i) strncpy (buffer, get_second_file (path), MAXLINE);
	}
	if (!_use_names)
		LOG ("Using library name %s\n", _library_name);
	if ((_interleave|_use_second_file) && nl[1] != nl[0])
		ERROR("Paired-ends do not have the same number of reads.\n");

	/* prepare output files */
	buffered_file fo[2], *pfo[2];
	int file = 0;
	for (int i = 0; i < 1 + _use_second_file; i++) {
		f_init (fo + i, IO_SYS);
		pfo[i] = fo + i;
	
	//	if (_split_reads)
	//		snprintf((char*)buffer, MAXLINE, "%s_%d.part%02d.fastq", out,i+1, file++);
		if (_interleave && !strcmp("-",out))
			snprintf((char*)buffer, MAXLINE, "-");
		else
			snprintf((char*)buffer, MAXLINE, "%s_%d.fastq", out,i+1);
		f_open (fo + i, (char*)buffer, IO_WRITE);
	}
	if (_interleave)
		pfo[1] = fo;

	int bc[2]; bc[0] = SZ_QUAL(len[0]); bc[1] = SZ_QUAL(len[1]);
//	LOG("--->%d\n",bc[0]);
	char alphabet[] = "ACGT";
	uint8_t l[MAXLINE], o[MAXLINE], ox[MAXLINE], chr, off;
	char names = 0; int64_t nameidx = 0;
	int64_t limit = (_split_reads ? _split_reads : nl[0]);
	char library[MAXLINE+1];
	if (!_use_names) {
		strncpy (library, _library_name, MAXLINE);
		nameidx = 0;
	}
	else { // both ends should use same library name and starting number
		f_read(fN, &names, 1);
		if (_use_second_file|_interleave) /* repeat */ 
			f_read(fN+1, &names, 1);
		if (!names) {
			f_read(fN, &nameidx, sizeof(int64_t));	
			int t = f_read (fN, library, MAXLINE);
			library[t] = 0;

			if (_use_second_file|_interleave) { /* repeat */ 
				f_read(fN+1, &nameidx, sizeof(int64_t));	
				int t = f_read (fN+1, library, MAXLINE);
				library[t] = 0;
			}
		}
		else {
			for(int ee=0;ee<(_use_second_file|_interleave)+1;ee++) {
			f_read(fN+ee, &sanger, 1);
			if (sanger) {
				int l;
				f_read(fN+ee, &l, sizeof(int)); f_read(fN+ee,prefix,l);prefix[l]=0;
				f_read(fN+ee, &l, sizeof(int)); f_read(fN+ee,machine_name,l); machine_name[l]=0;
				f_read(fN+ee, &l, sizeof(int)); f_read(fN+ee,run_id,l); run_id[l]=0;
			}
			}
		}
	}

	int read_next_read_info = 0;
	int32_t core, corlen;
	uint32_t n_id = 0, n_lane = 0, n_i, n_l, n_x, n_y;
	for (int64_t K = 0; ; K++) {
		for (int F = 0; F < 1 + (_interleave|_use_second_file); F++) {
		/* names */
			if (K == read_next_read_info && !F) {
				uint64_t len = 0;
				if (f_read(fR+F, &core, sizeof(int32_t))!=sizeof(int32_t)) 
					goto end;
				f_read(fR+F, &len, sizeof(int64_t));
				read_next_read_info += len;

				corlen=(core==MAXBIN-1)?0:strlen(patterns[core]);
				n_id=n_lane=0;
			//	if(corlen)LOG("i=%d, next=%d, core=%s[%d]\n", K, read_next_read_info, patterns[core],core);
			}

			if (names) {
				if (sanger) {
					n_i=WI(fN+F);n_id+=n_i;
					n_l=WI(fN+F);n_lane+=n_l;
					n_x=WI(fN+F);
					n_y=WI(fN+F);

					sprintf(buffer,"%s.%u %s:%s:%u:%u:%u/%d",prefix,n_id,machine_name,run_id,n_lane,n_x,n_y,F+1);
					f_write(pfo[F], buffer, strlen(buffer));
					if(core==MAXBIN-1){n_id=n_lane=0;}
				}
				else {
					f_read(fN+F, &chr, 1);
					buffer[0] = '@'; 
					f_read(fN+F, buffer+1, chr);

					if((_interleave|_use_second_file) && chr > 0 && buffer[chr - 1] == '/')
						buffer[chr] = F + 1 + '0';

					f_write(pfo[F],buffer, chr+1);
				}
			}
			else {
				snprintf((char*)buffer, MAXLINE, "@%s.%lld", library, nameidx);
				f_write(pfo[F],buffer,strlen((char*)buffer));
			}
			if(_interleave && (!names ||  (chr > 0 && buffer[chr] != F + 1 - '0' && buffer[chr - 1] != '/'))) {
				snprintf((char*)buffer,3,"/%d",F+1);
				f_write(pfo[F],buffer, 2);
			}
			chr='\n';
			f_write(pfo[F],&chr,1);
			

		/* quals */
			Z=F; a_read(fQ + F, buffer, bc[F]);
			int qc = 0;
			for (int i = 0; i < bc[F]; i ++) buffer[i]+=phredOff[F]; 
				/*
				buffer[qc++] = phredOff[F] + (o[i] >> 2);
				buffer[qc++] = phredOff[F] + ( ((o[i + 0] << 4) | (o[i + 1] >> 4)) & 0x3f );
				buffer[qc++] = phredOff[F] + ( ((o[i + 1] << 2) | (o[i + 2] >> 6)) & 0x3f );
				buffer[qc++] = phredOff[F] + (o[i + 2] & 0x3f);
				LOG("%c%c%c%c [%d]\n", buffer[qc-4],buffer[qc-3],buffer[qc-2],buffer[qc-1],__gc); __gc+=3;
			}*/
			buffer[qc = bc[F]] = '\n';
			//bc = SZ_QUAL(len[F]);
			/*bc = 0;
			while (bc < len[F]) {
				f_read (fQ+F,&chr, 1);
				if (chr < 128) {
					buffer[bc++] = chr;
				}
				else {
					off = buffer[bc - 1];
					for (int i = 0; i < chr - 128; i++)
						buffer[bc++] = off;
				}
			}
			buffer[bc++] = '\n';*/

		/* reads */
					
			if (!F) {
				int r = f_read (fR+F, o, SZ_READ(len[F]-corlen));
				int16_t end; 
				f_read(fR+F, &end, sizeof(int16_t));
			//	LOG("[%d,%d,%d]\n",corlen,len[F],end);
				for(int i=0; i<len[F]-corlen; i++)
					ox[i] = alphabet[ (o[i/4] >> ((3 - i%4) * 2)) & 3 ];

				int le=0;
				if (end) {
					memcpy(l, ox+(len[F]-end-corlen), end); le+=end;
				}
				if (corlen) {
					memcpy(l+le, patterns[core], corlen); le+= corlen;
				}
				memcpy(l+le, ox, len[F]-end-corlen);
			}
			else {
				f_read(fR+F,o,SZ_READ(len[F]));
				for(int i=0; i<len[F]; i++)
					l[i] = alphabet[ (o[i/4] >> ((3 - i%4) * 2)) & 3 ];
			}
			for(int i = 0; i < len[F]; i++) {
				chr = (buffer[i] == phredOff[F]
						? 'N'
						: l[i] ); // alphabet[ (l[i/4] >> ((3 - i%4) * 2)) & 3 ]);
//				if (buffer[i] == '!') buffer[i] = phredOff[F];
				f_write (pfo[F], &chr, 1);
			}
			chr = '\n'; f_write (pfo[F],&chr, 1);
			chr = '+'; f_write (pfo[F],&chr, 1);
			chr = '\n'; f_write (pfo[F],&chr, 1);
			f_write (pfo[F], buffer, qc+1);

		}
		nameidx++;
	/*	limit--;
		if (_split_reads && !limit && K < nl[0] - 1) {
			f_close (pfo[0]);
			limit = _split_reads;
			LOG("Created part %d as %s with %d reads\n", file, pfo[0]->file_name,_split_reads);
			if (_use_second_file)
				LOG("Created part %d as %s with %d reads\n", file, pfo[1]->file_name,_split_reads);

			for (int i = 0; i < 1 + _use_second_file; i++) {
				f_close (pfo[i]);
				snprintf((char*)buffer, MAXLINE, "%s_%d.part%02d.fastq", out,i+1, file++);
				f_open (fo + i, (char*)buffer, IO_WRITE);
			}
		}*/
	}
end:;
//	int lR = (_split_reads? _split_reads-limit : nl[0]-limit);
	LOG("Created part %d as %s with %d reads\n", file, pfo[0]->file_name,read_next_read_info);
	if (_use_second_file)
		LOG("Created part %d as %s with %d reads\n", file, pfo[1]->file_name,read_next_read_info);

	f_free(fR);
	f_free(fQ);
	f_free(fN);
	f_free(fo);
	if (_interleave|_use_second_file) {
		f_free(fR+1);
		f_free(fN+1);
		f_free(fQ+1);
	}
	if(_use_second_file)
		f_free(fo+1);
	LOG("Time elapsed: %02d:%02d:%02d\n", TIME/3600, (TIME/60)%60, TIME%60);
}


