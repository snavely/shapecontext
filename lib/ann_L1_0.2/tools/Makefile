#-----------------------------------------------------------------------------
# Makefile for tools
#
# ANN: Approximate Nearest Neighbors
# Version: 0.1 (Beta release)
#-----------------------------------------------------------------------------
# Revision 0.1  03/04/98
#	Initial release
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Basic definitions
#		BASEDIR		where include, src, lib, ... are
#		INCDIR		include directory
#		LIBDIR		library directory
#		BINDIR		bin directory
#		LDFLAGS		loader flags
#		ANNLIB		ANN library
#		OTHERLIB	other libraries
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
#-----------------------------------------------------------------------------

ANN2FIG = ann2fig			# name of figure drawing program
SOURCES = ann2fig.cc			# source files
OBJECTS = $(SOURCES:.cc=.o)		# object files

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

.cc.o:
	$(C++) -c -I$(INCDIR) $(CFLAGS) $(CFLAGS2) $<

#-----------------------------------------------------------------------------
# Make the program
#-----------------------------------------------------------------------------

default: 
	@echo "Specify a target configuration"

targets: $(BINDIR)/$(ANN2FIG)

$(BINDIR)/$(ANN2FIG): $(OBJECTS)
	$(C++) $(OBJECTS) -o $(ANN2FIG) $(LDFLAGS) $(ANNLIBS) $(OTHERLIBS)
	mv $(ANN2FIG) $(BINDIR)

#-----------------------------------------------------------------------------
# configuration definitions
#-----------------------------------------------------------------------------

include ../Make-config

#-----------------------------------------------------------------------------
# Objects
#-----------------------------------------------------------------------------

ann2fig.o: ann2fig.cc
	$(C++) -c -I$(INCDIR) $(CFLAGS2) ann2fig.cc

#-----------------------------------------------------------------------------
# Cleaning
#-----------------------------------------------------------------------------

clean:
	-rm -f *.o core

realclean: clean
	-rm -f $(BINDIR)/$(ANN2FIG)
