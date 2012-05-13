/// 786

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#include <sys/sysinfo.h>

#include "const.h"
#include "compress.h"
#include "decompress.h"

// Argument globals
int       _quality_sample_lines     = 100000;
int       _quality_lossy_percentage = 0;
char      _use_second_file          = 0;
char      _use_names						= 1;
uint64_t	 _file_buffer_size         = 128 * 1024 * 1024;
uint64_t  _max_bucket_set_size      = 4LL * 1024LL * 1024LL * 1024LL;
char      _temp_directory[MAXLINE]  = "__temp__";
char      _output_path[MAXLINE]     = "";
char      _library_name[MAXLINE]		= "";
char      _pattern_path[MAXLINE];  
int       _split_reads              = 0;
int		 _compression_mode         = IO_GZIP;
char      _interleave               = 0;
int       _time_elapsed 				= 0;
int       _thread_count 				= 1;

extern char _binary_HELP_start;
extern char _binary_HELP_end;
void help () {
	for (char *c = &_binary_HELP_start; c != &_binary_HELP_end; c++)
		putchar (*c);
	exit(0);
}

void check_arguments (char **files, int length, int mode) {
	if (_interleave && _use_second_file)
		ERROR("Interleaved option (-i) cannot be used with paired-end option (-r).\n");
	if (strlen (_output_path) == 0)
		ERROR("No output file specified.\n");
	if (!_use_names && strlen (_library_name) == 0)
		ERROR("No library name specified.\n");
	if (mode == 1 && length > 1) 
		ERROR ("Too many files specified (decompression onyl supports one file).\n");
	if (_quality_lossy_percentage < 0 || _quality_lossy_percentage > 100)
		ERROR ("Percentage must be in range [0,100].\n");
	if (!strcmp (_output_path, "-") && !(_interleave && mode))
		ERROR ("stdout can be only used with interleaved files in decompression mode.\n");
	if (_thread_count > 1 && _compression_mode == IO_GZIP)
		_compression_mode = IO_PGZIP;

	struct stat s;
	if (!mode) {
		if (stat (_temp_directory, &s) == 0) {
			if (!S_ISDIR (s.st_mode)) ERROR ("%s exists, but it should be a directory.\n", _temp_directory); 
		}
		else if (mkdir (_temp_directory, 0777)) {
			ERROR("Cannot create directory %s.\n", _temp_directory);
		}
	}

	for (int i = 0; i < length; i++) {
		if (stat (files[i], &s) != 0)
			ERROR ("File %s does not exist or it is not accessible.\n", files[i]);
		if (_use_second_file || (_interleave && mode)) {
			if (!get_second_file (files[i]))
				ERROR ("Cannot get file name for paired end for file %s. File should contain character 1.\n", files[i]);
			if (stat (get_second_file (files[i]), &s) != 0)
				ERROR ("File %s does not exist or it is not accessible.\n", get_second_file (files[i]));
		}
	}
}

int main (int argc, char **argv) {
	_time_elapsed = time(0);

	// set default number of threads
	_thread_count = sysconf( _SC_NPROCESSORS_ONLN ) - 1;

	LOG("SCALCE %s [OpenMP; available cores=%d]\n", SCALCE_VERSION, _thread_count+1);

	int mode = 0, opt; // default -  compress
	struct option long_opt[] = {
		{ "help",              0, NULL, 'h' },
		{ "lossy-percentage",  1, NULL, 'p' },
		{ "file-buffer-size",  1, NULL, 'b' },
		{ "decompress",        0, NULL, 'd' },
		{ "compression",       1, NULL, 'c' },
		{ "output",            1, NULL, 'o' },
		{ "sample-size",       1, NULL, 's' },
		{ "patterns",          1, NULL, 'P' },
		{ "interleave",        0, NULL, 'i' },
		{ "temp-directory",    1, NULL, 't' },
		{ "bucket-set-size",   1, NULL, 'B' },
		{ "paired-end",        0, NULL, 'r' },
		{ "skip-names",        1, NULL, 'n' },
		{ "split-reads",       1, NULL, 'S' },
		{ "threads",           1, NULL, 'T' },
		{ "version",		     0, NULL, 'v' },
		{ NULL,                0, NULL,  0  }
	};
	do {
		opt = getopt_long (argc, argv, "vhp:b:T:dc:o:s:P:t:B:rin:S:", long_opt, NULL);
		switch (opt) {
			case 'v':
			//	LOG("%s\n",SCALCE_VERSION);
				exit(0);
			case 'h':
				help ();
				break;
			case 'i':
				_interleave = 1;
				break;
			case 'c':
				if (!strcmp (optarg, "bz"))
					_compression_mode = IO_BZIP;
				else if (!strcmp (optarg, "gz"))
					_compression_mode = IO_GZIP;
				else if (!strcmp (optarg, "pigz"))
					_compression_mode = IO_PGZIP;
				else if (!strcmp (optarg, "no"))
					_compression_mode = IO_SYS;
				else
					ERROR ("Unknown compression mode. See help for details.\n");
				break;
			case 'b':
			case 'B': {
				int l = strlen (optarg);
				char al = optarg[l - 1];
				uint64_t sz = 1024 * 1024LL;
				if (al == 'G') 
					sz *= 1024;
				else if (al == 'M') ;
				else {
					ERROR ("Size parameter must be ended with G or M.\n");
				}
				optarg[l - 1] = 0;
				if (opt == 'B') 
					_max_bucket_set_size = sz * atoi (optarg);
				else 
					_file_buffer_size = sz * atoi (optarg);
			}
			break;
			case 'p':
				_quality_lossy_percentage = atoi (optarg);
				break;
			case 'T':
				_thread_count = atoi (optarg);
				break;
			case 'S':
				_split_reads = atoi (optarg);
				break;
			case 'r':
				_use_second_file = 1;
				break;
			case 's':
				_quality_sample_lines = atoi (optarg);
				break;
			case 'd':
				mode = 1;
				break;
			case 'P':
				strncpy (_pattern_path, optarg, MAXLINE);
				break; 
			case 't':
				strncpy (_temp_directory, optarg, MAXLINE);
				break;
			case 'o':
				strncpy (_output_path, optarg, MAXLINE);
				break;
			case 'n':
				_use_names = 0;
				strncpy (_library_name, optarg, MAXLINE);
				break;
			case -1:
				break;
			default: {
				help ();
			}
		}
	} while (opt != -1);
	check_arguments (argv + optind, argc - optind, mode);
	if (mode == 1) { 
		decompress(argv[argc-1],_output_path);
	}
	else {
		compress (argv + optind, argc - optind, _output_path, _pattern_path);
	}

	LOG("Done!\n");
	return 0;
}

