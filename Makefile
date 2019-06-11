GPP=gcc
CCFLAGS=-msse -Wall
CCFLAGS+=-g
CCFLAGS+=-O3
LDFLAGS=-lm 

OBJECTS=lib/bfio.o \

all: fse

fse: fse.c 
	$(GPP) -pg $(CCFLAGS) $(OBJECTS) $^ -o $@ $(LDFLAGS)

clean: fse 
	rm fse
	rm lib/*.o

#fse: fse_template.c Makefile
#	$(CC) -o fse -Wall -O3 -lm fse_template.c
#	$(CC) -o debug_fse -Wall -g -lm fse_template.c

