# HOSTNAME=jimantha

CC=g++
CPP=g++
OBJS=shapecontext.o match.o emd.o register.o

MKL_PATH=/hydra/S1/local/intel/Compiler/11.0/083/mkl/lib/em64t
INTEL_PATH=/hydra/S1/local/intel/Compiler/11.0/083/lib/intel64

INCLUDE_PATH=-I../lib/imagelib -I../lib/matrix -I../lib/ann_L1_1.1/include
LIB_PATH=-L../lib/imagelib -L../lib/matrix -L../lib/minpack -L../lib/ann_L1_0.2/lib -L$(MKL_PATH) -L$(INTEL_PATH)
# LIBS=-limage -lm -llapack -lminpack -lblas -lf2c -ldummy -lANN

LA_LIBS=-lmkl_intel_thread -lmkl_core -lmkl_intel_lp64 -lmkl_lapack \
	-liomp5 -lpthread 
LIBS=-limage -lmatrix -lm $(LA_LIBS) -lminpack -lANN

CFLAGS=-Wall -fopenmp -O3 -ffast-math -Wall -mfpmath=sse -msse2 -mssse3 -funroll-loops -march=core2 $(INCLUDE_PATH)
CPPFLAGS=-Wall -fopenmp -O3 -ffast-math -Wall -mfpmath=sse -msse2 -mssse3 -funroll-loops -march=core2 $(INCLUDE_PATH)

all: sc3 sc_test pts2dmap drawpts morphimg2

sc2: $(OBJS) sc2.o
	$(CPP) $(CFLAGS) -o $@ $(OBJS) sc2.o $(LIB_PATH) $(LIBS)
#	cp $@ ../bin

sc3: $(OBJS) sc3.o
	$(CPP) $(CFLAGS) -o $@ $(OBJS) sc3.o $(LIB_PATH) $(LIBS)
#	cp $@ ../bin

sc_test: $(OBJS) sc_test.o
	$(CPP) $(CFLAGS) -o $@ $(OBJS) sc_test.o $(LIB_PATH) $(LIBS)
#	cp $@ ../bin

pts2dmap: pts2dmap.o
	$(CPP) $(CFLAGS) -o $@ pts2dmap.o $(LIB_PATH) $(LIBS)
#	cp $@ ../bin

drawpts: drawpts.o
	$(CPP) $(CFLAGS) -o $@ drawpts.o $(LIB_PATH) $(LIBS)
#	cp $@ ../bin

morphimg2: morphimg2.o morph2.o
	$(CPP) $(CFLAGS) -o $@ $^ $(LIB_PATH) $(LIBS)

clean:
	rm -f *.o *~ sc
