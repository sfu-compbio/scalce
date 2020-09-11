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

#include <getopt.h>
#include <locale.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/sysctl.h>
#include <sys/sysinfo.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include "buffio.h"
#include "compress.h"
#include "const.h"
#include "decompress.h"

//#define MIN(a,b) ((a)<(b)?(a):(b))

// Argument globals
int _quality_sample_lines = 100000;
int _quality_lossy_percentage = 0;
char _use_second_file = 0;
char _is_fasta = 0;
char _use_names = 1;
uint64_t _file_buffer_size = 128 * 1024 * 1024;
uint64_t _max_bucket_set_size = 4LL * 1024LL * 1024LL * 1024LL;
char _temp_directory[MAXLINE] = "__temp__";
char _output_path[MAXLINE] = "";
char _library_name[MAXLINE] = "";
char _pattern_path[MAXLINE];
int _split_reads = 0;
int _compression_mode = IO_GZIP;
char _interleave = 0;
int64_t _time_elapsed = 0;
int _thread_count = 1;
int _decompress = 0;
int _no_ac = 0;
int _compress_qualities = 1;

extern char _binary_HELP_start;
extern char _binary_HELP_end;
void help() {
  for (char *c = &_binary_HELP_start; c != &_binary_HELP_end; c++)
    putchar(*c);
  exit(0);
}

// based on
// http://stackoverflow.com/questions/12765010/c-test-for-file-existence-before-calling-execvp
int is_file(const char *path) {
  struct stat buf;
  stat(path, &buf);
  return S_ISREG(buf.st_mode);
}
char check_pigz() {
  char *path = getenv("PATH");
  char *item = NULL;
  int found = 0;

  if (!path)
    return 0;

  path = strdup(path);
  char real_path[MAXLINE]; // or PATH_MAX or something smarter
  for (char *item = strtok(path, ":"); !found && item;
       item = strtok(NULL, ":")) {
    sprintf(real_path, "%s/pigz", item);
    // printf("Testing %s\n", real_path);
    if (is_file(real_path) &&
        !(access(real_path, F_OK) || access(real_path, X_OK)))
      // check if the file exists and is executable
      found = 1;
  }
  free(path);
  return found;
}

void check_arguments(char **files, int length, int mode) {
  if (_compression_mode == IO_PGZIP && !check_pigz()) {
    LOG("** pigz not found in PATH. Using zlib compression instead. **\n");
    _compression_mode = IO_GZIP;
  }

  if (_interleave && _use_second_file)
    ERROR("Interleaved option (-i) cannot be used with paired-end option "
          "(-r).\n");
  if (strlen(_output_path) == 0)
    ERROR("No output file specified.\n");
  if (!_use_names && strlen(_library_name) == 0)
    ERROR("No library name specified.\n");
  if (mode == 1 && length > 1)
    ERROR("Too many files specified (decompression only supports one file).\n");
  if (_quality_lossy_percentage < 0 || _quality_lossy_percentage > 100)
    ERROR("Percentage must be in range [0,100].\n");
  if (!strcmp(_output_path, "-") && (_split_reads || _use_second_file))
    ERROR("stdout can be only used with single-end file decompression. It "
          "cannot be used with --split-reads option!\n");

  struct stat s;
  if (!mode) {
    if (stat(_temp_directory, &s) == 0) {
      if (!S_ISDIR(s.st_mode))
        ERROR("%s exists, but it should be a directory.\n", _temp_directory);
    } else if (mkdir(_temp_directory, 0777)) {
      ERROR("Cannot create directory %s.\n", _temp_directory);
    }
  }

  for (int i = 0; i < length; i++) {
    if (stat(files[i], &s) != 0)
      ERROR("File %s does not exist or it is not accessible.\n", files[i]);
    if (_use_second_file || (_interleave && mode)) {
      if (!get_second_file(files[i]))
        ERROR("Cannot get file name for paired end for file %s. File should "
              "contain character 1.\n",
              files[i]);
      if (stat(get_second_file(files[i]), &s) != 0)
        ERROR("File %s does not exist or it is not accessible.\n",
              get_second_file(files[i]));
    }
  }
}

int main(int argc, char **argv) {
  setlocale(LC_ALL, "");
  _time_elapsed = TIME;

  // set default number of threads
  _thread_count = MIN(4, sysconf(_SC_NPROCESSORS_ONLN) - 1);

  LOG("SCALCE %s [pthreads; available cores=%d]\n", SCALCE_VERSION,
      _thread_count + 1);
#ifdef PACBIO
  LOG("!! This is (very) experimental version which (supposedly) supports "
      "variable length long reads !!\n");
  LOG("!! If using Illumina, compile without -DPACBIO -- faster and safer "
      "!!\n");
#endif
  if (_thread_count > 1)
    _compression_mode = IO_PGZIP;
  else
    _compression_mode = IO_GZIP;

  int mode = 0, opt; // default -  compress
  struct option long_opt[] = {{"help", 0, NULL, 'h'},
                              {"lossy-percentage", 1, NULL, 'p'},
                              //	{ "file-buffer-size",  1, NULL, 'b' },
                              {"decompress", 0, NULL, 'd'},
                              {"compression", 1, NULL, 'c'},
                              {"output", 1, NULL, 'o'},
                              {"sample-size", 1, NULL, 's'},
                              {"no-qualities", 0, NULL, 'Q'},
                              {"patterns", 1, NULL, 'P'},
                              //	{ "interleave",        0, NULL, 'i' },
                              {"temp-directory", 1, NULL, 't'},
                              {"bucket-set-size", 1, NULL, 'B'},
                              {"paired-end", 0, NULL, 'r'},
                              {"skip-names", 1, NULL, 'n'},
                              {"split-reads", 1, NULL, 'S'},
                              {"fasta", 0, NULL, 'f'},
                              {"threads", 1, NULL, 'T'},
                              {"version", 0, NULL, 'v'},
                              {"no-arithmetic", 0, NULL, 'A'},
                              {NULL, 0, NULL, 0}};
  do {
    opt = getopt_long(argc, argv, "vhp:T:dc:o:fs:t:B:rQAn:P:S:" /*"i"*/,
                      long_opt, NULL);
    switch (opt) {
    case 'v':
      //	LOG("%s\n",SCALCE_VERSION);
      exit(0);
    case 'h':
      help();
      break;
    case 'A':
      _no_ac = 1;
      break;
    case 'f':
      _is_fasta = 1;
      _compress_qualities = 0;
      break;
    //	case 'i':
    //		_interleave = 1;
    //		break;
    case 'c':
      if (!strcmp(optarg, "bz"))
        _compression_mode = IO_BZIP;
      else if (!strcmp(optarg, "gz"))
        _compression_mode = IO_GZIP;
      else if (!strcmp(optarg, "pigz"))
        _compression_mode = IO_PGZIP;
      else if (!strcmp(optarg, "no"))
        _compression_mode = IO_SYS;
      else
        ERROR("Unknown compression mode. See help for details.\n");
      break;
    case 'B': {
      int l = strlen(optarg);
      char al = optarg[l - 1];
      uint64_t sz = 1024 * 1024LL;
      if (al == 'G')
        sz *= 1024;
      else if (al == 'M')
        ;
      else {
        ERROR("Size parameter must be ended with G or M.\n");
      }
      optarg[l - 1] = 0;
      if (opt == 'B')
        _max_bucket_set_size = sz * atoi(optarg);
      else
        _file_buffer_size = sz * atoi(optarg);
    } break;
    case 'p':
      _quality_lossy_percentage = atoi(optarg);
      break;
    case 'T':
      _thread_count = atoi(optarg);
      if (_thread_count == 1 && _compression_mode == IO_PGZIP)
        _compression_mode = IO_GZIP;
      break;
    case 'S':
      _split_reads = atoi(optarg);
      break;
    case 'r':
      _use_second_file = 1;
      break;
    case 'Q':
      _compress_qualities = 0;
      break;
    case 's':
      _quality_sample_lines = atoi(optarg);
      break;
    case 'd':
      mode = 1;
      break;
    case 'P':
      strncpy(_pattern_path, optarg, MAXLINE);
      break;
    case 't':
      strncpy(_temp_directory, optarg, MAXLINE);
      break;
    case 'o':
      strncpy(_output_path, optarg, MAXLINE);
      break;
    case 'n':
      _use_names = 0;
      strncpy(_library_name, optarg, MAXLINE);
      break;
    case -1:
      break;
    default: { help(); }
    }
  } while (opt != -1);
  check_arguments(argv + optind, argc - optind, mode);
  if (mode == 1) {
    _decompress = 1;
    decompress(argv[argc - 1], _output_path);
  } else {
    compress(argv + optind, argc - optind, _output_path, _pattern_path);
  }

  LOG("Done!\n");
  return 0;
}
