#-----------------------------------------------------------------------------
# Makefile for core library
#
# ANN: Approximate Nearest Neighbors
# Version: 0.1 (Beta release)
#-----------------------------------------------------------------------------
# Revision 0.1  03/04/98
#	Initial release
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Some basic definitions:
#		BASEDIR		where include, src, lib, ... are
#		INCLIB		include directory
#		LIBLIB		library directory
#-----------------------------------------------------------------------------
BASEDIR	= ..
INCDIR	= $(BASEDIR)/include
LIBDIR	= $(BASEDIR)/lib

SOURCES = ANN.cc brute.cc kd_tree.cc kd_util.cc kd_split.cc \
	kd_search.cc kd_pr_search.cc bd_tree.cc bd_search.cc \
	bd_pr_search.cc perf.cc

HEADERS = kd_tree.h kd_split.h kd_util.h kd_search.h \
	kd_pr_search.h perf.h pr_queue.h pr_queue_k.h

OBJECTS = $(SOURCES:.cc=.o)

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

.cc.o:
	$(C++) -c -I$(INCDIR) $(CFLAGS) $(CFLAGS2) $<

#-----------------------------------------------------------------------------
# Make the library
#-----------------------------------------------------------------------------

default:
	@echo "Specify a target configuration"

targets: $(LIBDIR)/$(ANNLIB)

$(LIBDIR)/$(ANNLIB): $(OBJECTS)
	$(MAKELIB) $(ANNLIB) $(OBJECTS)
	$(RANLIB) $(ANNLIB)
	mv $(ANNLIB) $(LIBDIR)

#-----------------------------------------------------------------------------
# Make object files
#-----------------------------------------------------------------------------

ANN.o: ANN.cc
	$(C++) -c -I$(INCDIR) $(CFLAGS) $(CFLAGS2) ANN.cc

brute.o: brute.cc
	$(C++) -c -I$(INCDIR) $(CFLAGS) $(CFLAGS2) brute.cc

kd_tree.o: kd_tree.cc
	$(C++) -c -I$(INCDIR) $(CFLAGS) $(CFLAGS2) kd_tree.cc

kd_util.o: kd_util.cc
	$(C++) -c -I$(INCDIR) $(CFLAGS) $(CFLAGS2) kd_util.cc

kd_split.o: kd_split.cc
	$(C++) -c -I$(INCDIR) $(CFLAGS) $(CFLAGS2) kd_split.cc

kd_search.o: kd_search.cc
	$(C++) -c -I$(INCDIR) $(CFLAGS) $(CFLAGS2) kd_search.cc

kd_pr_search.o: kd_pr_search.cc
	$(C++) -c -I$(INCDIR) $(CFLAGS) $(CFLAGS2) kd_pr_search.cc

bd_tree.o: bd_tree.cc
	$(C++) -c -I$(INCDIR) $(CFLAGS) $(CFLAGS2) bd_tree.cc

bd_search.o: bd_search.cc
	$(C++) -c -I$(INCDIR) $(CFLAGS) $(CFLAGS2) bd_search.cc

bd_pr_search.o: bd_pr_search.cc
	$(C++) -c -I$(INCDIR) $(CFLAGS) $(CFLAGS2) bd_pr_search.cc

perf.o: perf.cc
	$(C++) -c -I$(INCDIR) $(CFLAGS) $(CFLAGS2) perf.cc

#-----------------------------------------------------------------------------
# Configuration definitions
#-----------------------------------------------------------------------------

include ../Make-config

#-----------------------------------------------------------------------------
# Cleaning
#-----------------------------------------------------------------------------

clean:
	-rm -f *.o core

realclean: clean

