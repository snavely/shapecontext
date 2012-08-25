#-----------------------------------------------------------------------------
# Makefile for the sample program
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
#		ANNLIBS		libraries
#-----------------------------------------------------------------------------

BASEDIR	= ..
INCDIR	= $(BASEDIR)/include
LIBDIR	= $(BASEDIR)/lib
BINDIR	= $(BASEDIR)/bin
LDFLAGS	= -L$(LIBDIR)
ANNLIBS	= -lANN -lm

#-----------------------------------------------------------------------------
# Some more definitions
#		ANNSAMP		name of sample program
#-----------------------------------------------------------------------------

ANNSAMP = ann_sample

SAMPSOURCES = ann_sample.cc
SAMPOBJECTS = $(SAMPSOURCES:.cc=.o)

#-----------------------------------------------------------------------------
# Make the program
#-----------------------------------------------------------------------------

default: 
	@echo "Specify a target configuration"

targets: $(BINDIR)/$(ANNSAMP)

$(BINDIR)/$(ANNSAMP): $(SAMPOBJECTS) $(LIBDIR)/$(ANNLIB)
	$(C++) $(SAMPOBJECTS) -o $(ANNSAMP) $(LDFLAGS) $(ANNLIBS)
	mv $(ANNSAMP) $(BINDIR)

#-----------------------------------------------------------------------------
# configuration definitions
#-----------------------------------------------------------------------------

include ../Make-config

#-----------------------------------------------------------------------------
# Objects
#-----------------------------------------------------------------------------

ann_sample.o: ann_sample.cc
	$(C++) -c -I$(INCDIR) $(CFLAGS2) ann_sample.cc

#-----------------------------------------------------------------------------
# Cleaning
#-----------------------------------------------------------------------------

clean:
	-rm -f *.o core

realclean: clean
	-rm -f $(BINDIR)/$(ANNSAMP)
