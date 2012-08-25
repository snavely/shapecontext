MACHTYPE=$(shell uname -m)

ifeq ($(MACHTYPE), x86_64)
P=x64
else
P=ia32
endif

OBJS=dpmpar.o \
	enorm.o \
	fdjac2.o \
	lmdif1.o \
	lmdif.o \
	lmpar.o \
	qrfac.o \
	qrsolv.o

AR=ar
FC=gcc
FFLAGS=-O2 -ffast-math

all: libminpack.a

libminpack.a: $(OBJS)
	$(AR) rs $@ $(OBJS)
	cp $@ ../

%.o : %.f
	$(FC) -c -o $@ $(FFLAGS) $<

clean:
	rm -f *.o *~ libminpack.a
