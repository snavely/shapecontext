# Makefile for imagelib

MACHTYPE=$(shell uname -m)

GCC			= g++

CC=gcc
# OPTFLAGS=-g2
OPTFLAGS=-O3 -ffast-math -Wall -mfpmath=sse -msse2 -funroll-loops -march=core2
OTHERFLAGS=-Wall -Wno-unused-result

IMAGELIB_OBJS= affine.o alignicp.o alignlm.o anniface.o anntool.o bmp.o canny.o circle.o color.o cspnd.o dmap.o dmap-io.o	\
	error_metric.o fileio.o filter.o fit.o fmatrix.o histogram.o	\
	homography.o horn.o image.o lerp.o		\
	morphology.o pgm.o poly.o pyramid.o pyramid-io.o qsort.o  \
	ransac.o resample.o tps.o transform.o transform-opt.o triangulate.o	\
	util.o kruppa.o vanish.o

INCLUDE_PATH=-I../matrix -I../ann_L1_0.2/include -I../../include

CFLAGS=$(OTHERFLAGS) $(OPTFLAGS) $(INCLUDE_PATH)
CPPFLAGS=$(OTHERFLAGS) $(OPTFLAGS) $(INCLUDE_PATH)

IMAGELIB=libimage.a

%.o : %.c
	$(CC) -c -o $@ $(CFLAGS) $<

all: $(IMAGELIB)

$(IMAGELIB): $(IMAGELIB_OBJS)
	ar r $(IMAGELIB) $(IMAGELIB_OBJS)
	cp $(IMAGELIB) ..

clean:
	rm -f *.o *~ $(IMAGELIB)
