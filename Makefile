# 786

CC=g++
CFLAGS=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -DSCALCE_VERSION=\"2.8\"
LDFLAGS=-lm -lpthread 

BZIP_VER=1.0.6
ZLIB_VER=1.2.10

CFLAGS+= -I zlib-$(ZLIB_VER) -I bzip2-$(BZIP_VER)
LDFLAGS+= zlib-$(ZLIB_VER)/libz.a bzip2-$(BZIP_VER)/libbz2.a

OBJECTS=HELP.o patterns.o $(SOURCES:.cpp=.o)
EXECUTABLE=scalce
SOURCES=const.cpp buffio.cpp arithmetic.cpp main.cpp names.cpp qualities.cpp reads.cpp compress.cpp decompress.cpp

all: CFLAGS+= -O3 -DNDEBUG
all: zlib bzlib $(SOURCES) $(EXECUTABLE) 

debug: CFLAGS+= -g
debug: $(SOURCES) $(EXECUTABLE) 

pacbio: CFLAGS+=-DPACBIO -O3 -DNDEBUG
pacbio: $(SOURCES) scalce-pacbio

scalce-pacbio: $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) -c $(CFLAGS) $< -o $@

HELP.o:
	ld -r -b binary -o HELP.o HELP

patterns.o:
	ld -r -b binary -o patterns.o patterns.bin

clean:
	cd zlib-$(ZLIB_VER) && make clean 
	cd bzip2-$(BZIP_VER) && make clean
	rm -f *.o scalce

download:
	curl -L http://www.bzip.org/$(BZIP_VER)/bzip2-$(BZIP_VER).tar.gz | tar zxvf -
	curl -L http://zlib.net/zlib-$(ZLIB_VER).tar.gz | tar zxvf -

zlib:
	cd zlib-$(ZLIB_VER) && ./configure && make -j

bzlib:
	cd bzip2-$(BZIP_VER) && make -j

