# 786

CC=g++
CFLAGS= -c -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -DSCALCE_VERSION=\"2.1\" -O2
LDFLAGS=  -lm  bzip2-1.0.6/libbz2.a zlib-1.2.5/libz.a -lpthread -O2 #-g -lprofiler -ltcmalloc
SOURCES=const.cpp buffio.cpp arithmetic.cpp main.cpp names.cpp qualities.cpp reads.cpp compress.cpp decompress.cpp
OBJECTS=HELP.o patterns.o $(SOURCES:.cpp=.o)
EXECUTABLE=scalce

all: $(SOURCES) $(EXECUTABLE) 

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

HELP.o:
	ld -r -b binary -o HELP.o HELP

patterns.o:
	ld -r -b binary -o patterns.o patterns.bin

clean:
	rm -f *.o scalce

clean-all: clean
	make clean -C zlib-1.2.5 
	make clean -C bzip2-1.0.6 

libs:
	make clean -C zlib-1.2.5 
	make -C zlib-1.2.5 
	make clean -C bzip2-1.0.6 
	make -C bzip2-1.0.6 
