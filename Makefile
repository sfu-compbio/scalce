# 786

CC=g++
CFLAGS=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -DSCALCE_VERSION=\"2.8\"
LDFLAGS=-lm -lpthread -lbz2 -lz

OBJECTS=HELP.o patterns.o $(SOURCES:.cpp=.o)
EXECUTABLE=scalce
SOURCES=const.cpp buffio.cpp arithmetic.cpp main.cpp names.cpp qualities.cpp reads.cpp compress.cpp decompress.cpp

all: CFLAGS+= -O3 -DNDEBUG
all: $(SOURCES) $(EXECUTABLE) 

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
	rm -f *.o scalce

