# Makefile for matrix library

MACHTYPE=$(shell uname -m)

ifeq ($(MACHTYPE), x86_64)
P=x64
else
P=ia32
endif

CC=gcc
OPTFLAGS=-O3 -ffast-math -Wall -mfpmath=sse -msse2 -funroll-loops -march=native -mtune=native
#OPTFLAGS=-g2
OTHERFLAGS=-Wall -Wno-unused-result -D__NO_MKL__

INCLUDE_PATH = -I../imagelib -I../../include
CFLAGS = $(OPTFLAGS) $(OTHERFLAGS) $(INCLUDE_PATH)

TARGET = libmatrix.a
OBJS = matrix.o vector.o svd.o

%.o : %.c
	$(CC) -c -o $@ $(CFLAGS) $<

all: $(TARGET)

$(TARGET): $(OBJS)
	ar r $@ $(OBJS)
	cp $@ ..

clean:
	rm -f *.o $(TARGET) *~
