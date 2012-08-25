#-----------------------------------------------------------------------------
# Makefile for the test and evaluation program
#
# ANN: Approximate Nearest Neighbors
# Version: 0.1 (Beta release)
#-----------------------------------------------------------------------------
# Revision 0.1  03/04/98
#	Initial release
#-----------------------------------------------------------------------------
# Note: For full performance measurements, it is assumed that the library
# and this program have both been compiled with the -DPERF flag.  See the
# Makefile in the ANN base directory for this flag.
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Basic definitions
#		BASEDIR		where include, src, lib, ... are
#		INCDIR		include directory
#		LIBDIR		library directory
#		BINDIR		bin directory
#		LDFLAGS		loader flags
#		ANNLIBS		ANN library
#		OTHERLIBS	other libraries
#-----------------------------------------------------------------------------

BASEDIR	= ..
INCDIR	= $(BASEDIR)/include
LIBDIR	= $(BASEDIR)/lib
BINDIR	= $(BASEDIR)/bin
LDFLAGS	= -L$(LIBDIR)
ANNLIBS	= -lANN
OTHERLIBS = -lm

#-----------------------------------------------------------------------------
# Some more definitions
#		ANNTEST		name of test program
#		PURIFTY		your purify (for memory leaks)
#-----------------------------------------------------------------------------

ANNTEST = ann_test
PURIFY  = purify

HEADERS = rand.h
TESTSOURCES = ann_test.cc rand.cc
TESTOBJECTS = $(TESTSOURCES:.cc=.o)

#-----------------------------------------------------------------------------
# Make the program
#-----------------------------------------------------------------------------

default: 
	@echo "Specify a target configuration"

targets: $(BINDIR)/$(ANNTEST)

$(BINDIR)/$(ANNTEST): $(TESTOBJECTS) $(LIBDIR)/$(ANNLIB)
	$(C++) $(TESTOBJECTS) -o $(ANNTEST) $(LDFLAGS) $(ANNLIBS) $(OTHERLIBS)
	mv $(ANNTEST) $(BINDIR)

purify: $(TESTOBJECTS)
	$(PURIFY) $(C++) $(TESTOBJECTS) -o $(ANNTEST) $(LDFLAGS) $(ANNLIBS) $(OTHERLIBS)

#-----------------------------------------------------------------------------
# configuration definitions
#-----------------------------------------------------------------------------

include ../Make-config

#-----------------------------------------------------------------------------
# Objects
#-----------------------------------------------------------------------------

ann_test.o: ann_test.cc
	$(C++) -c -I$(INCDIR) $(CFLAGS) $(CFLAGS2) ann_test.cc

rand.o: rand.cc
	$(C++) -c -I$(INCDIR) $(CFLAGS) $(CFLAGS2) rand.cc

#-----------------------------------------------------------------------------
# Cleaning
#-----------------------------------------------------------------------------

clean:
	-rm -f *.o core

realclean: clean
	-rm -f $(BINDIR)/$(ANNTEST)
