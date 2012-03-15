# 786

CC=gcc
CFLAGS= -c -O2 -std=c99 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -DSCALCE_VERSION=\"2.0_very_beta\" -fopenmp
LDFLAGS= -s -static -lm -Wl,-Bstatic zlib-1.2.5/libz.a bzip2-1.0.6/libbz2.a -fopenmp 
SOURCES=const.c buffio.c main.c names.c qualities.c reads.c compress.c decompress.c arithmetic.c
OBJECTS=HELP.o patterns.o $(SOURCES:.c=.o)
EXECUTABLE=scalce

all: $(SOURCES) $(EXECUTABLE) 

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

HELP.o:
	ld -r -b binary -o HELP.o HELP

patterns.o:
	ld -r -b binary -o patterns.o patterns.bin

clean:
	make clean -C zlib-1.2.5 
	make clean -C bzip2-1.0.6 
	rm -f *.o scalce

libs:
	make clean -C zlib-1.2.5 
	make -C zlib-1.2.5 
	make clean -C bzip2-1.0.6 
	make -C bzip2-1.0.6 
