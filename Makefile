#
# Makefile for GFMKPS
#

CC = g++

all: mkps

gfchar.o: gfchar.cpp gfchar.h
	$(CC) -c gfchar.cpp

mkps.o: mkps.cpp mkps.h
	$(CC) -c mkps.cpp

mkps: mkps.o gfchar.o

clean: 
	rm -f *~ *.o mkps
