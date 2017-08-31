CC=g++
CPP=g++
OBJS=shapecontext.o match.o emd.o register.o

INCLUDE_PATH=-I../lib/imagelib -I../lib/matrix -I../lib/ann_L1_1.1/include
LIB_PATH=-L../lib/imagelib -L../lib/matrix -L../lib/minpack -L../lib/ann_L1_0.2/lib

LA_LIBS=-lblas -llapack
LIBS=-limage -lmatrix -lm $(LA_LIBS) -lminpack -lANN

CFLAGS=-Wall -fopenmp -O3 $(INCLUDE_PATH)
CPPFLAGS=-Wall -fopenmp -O3 $(INCLUDE_PATH)

all: sc2 sc3

sc2: $(OBJS) sc2.o
	$(CPP) $(CFLAGS) -o $@ $(OBJS) sc2.o $(LIB_PATH) $(LIBS)

sc3: $(OBJS) sc3.o
	$(CPP) $(CFLAGS) -o $@ $(OBJS) sc3.o $(LIB_PATH) $(LIBS)

sc_test: $(OBJS) sc_test.o
	$(CPP) $(CFLAGS) -o $@ $(OBJS) sc_test.o $(LIB_PATH) $(LIBS)

pts2dmap: pts2dmap.o
	$(CPP) $(CFLAGS) -o $@ pts2dmap.o $(LIB_PATH) $(LIBS)

drawpts: drawpts.o
	$(CPP) $(CFLAGS) -o $@ drawpts.o $(LIB_PATH) $(LIBS)

morphimg2: morphimg2.o morph2.o
	$(CPP) $(CFLAGS) -o $@ $^ $(LIB_PATH) $(LIBS)

clean:
	rm -f *.o *~ sc
