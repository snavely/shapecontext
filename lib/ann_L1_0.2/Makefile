#-----------------------------------------------------------------------------
# Top-level Makefile for ANN.
#
# ANN: Approximate Nearest Neighbors
# Version: 0.1 (Beta release)
#-----------------------------------------------------------------------------
# Copyright (c) 1997 University of Maryland and Sunil Arya and David Mount
# All Rights Reserved.
# 
# This software and related documentation is part of the 
# Approximate Nearest Neighbor Library (ANN).
# 
# Permission to use, copy, and distribute this software and its 
# documentation is hereby granted free of charge, provided that 
# (1) it is not a component of a commercial product, and 
# (2) this notice appears in all copies of the software and
#     related documentation. 
# 
# The University of Maryland and the authors make no representations
# about the suitability or fitness of this software for any purpose.
# It is provided "as is" without express or implied warranty.
#-----------------------------------------------------------------------------
# Revision 0.1  09/06/97
#	Initial release
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# default: list the options
#-----------------------------------------------------------------------------
default:
	@echo "Enter one of the following:"
	@echo "  make sunos4               for Sun with SunOS 4.x"
	@echo "  make sunos4-g++           for Sun with SunOS 4.x and g++"
	@echo "  make sunos5               for Sun with SunOS 5.x"
	@echo "  make sunos5-sl            for Sun with SunOS 5.x, make shared libs"
	@echo "  make sunos5-g++           for Sun with SunOS 5.x and g++"
	@echo "  make sunos5-g++-sl        for Sun with SunOS 5.x, g++, make shared libs"
	@echo "  make sgi                  for SGI with CC"
	@echo "  make alpha-g++            for DEC Alpha and g++"
	@echo "  make clean                remove .o files"
	@echo "  make realclean            remove .o, library and executable files"
	@echo " "
	@echo "See file Makefile for other compilation options, such as disabling"
	@echo "performance measurement code."

#-----------------------------------------------------------------------------
# main make entry point
#-----------------------------------------------------------------------------
sgi sunos4 sunos4-g++ sunos5 sunos5-g++ sunos5-g++-sl alpha-g++ authors-debug:
	cd src ; $(MAKE) $@ "CFLAGS2 = $(CFLAGS2)"
	cd test ; $(MAKE) $@ "CFLAGS2 = $(CFLAGS2)"
	cd tools ; $(MAKE) $@
#	cd samples ; $(MAKE) $@

#-----------------------------------------------------------------------------
# Remove .o files and core files
#-----------------------------------------------------------------------------
clean:
	cd src ; $(MAKE) clean
	cd test ; $(MAKE) clean
	cd tools ; $(MAKE) clean
	cd samples ; $(MAKE) clean

#-----------------------------------------------------------------------------
# Remove everthing that can be remade
#-----------------------------------------------------------------------------
realclean: 
	-rm -f lib/*.a lib/*.so
	cd src ; $(MAKE) realclean
	cd test ; $(MAKE) realclean
	cd tools ; $(MAKE) realclean
	cd samples ; $(MAKE) realclean

#-----------------------------------------------------------------------------
# Make distribution package (for use by authors only)
#-----------------------------------------------------------------------------
DISTR = ann_0.1

distribution: realclean 
	cd doc ; $(MAKE) 
	cd .. ; mkdir $(DISTR)
	cp Make-config ../$(DISTR)
	cp Makefile ../$(DISTR)
	cp README ../$(DISTR)
	cp -r bin ../$(DISTR)
	cp -r include ../$(DISTR)
	cp -r lib ../$(DISTR)
	cp -r samples ../$(DISTR)
	cp -r src ../$(DISTR)
	cp -r test ../$(DISTR)
	cp -r tools ../$(DISTR)
	cd ../$(DISTR); mkdir doc
	cp doc/*.ps ../$(DISTR)/doc
	cd .. ; tar -cfv $(DISTR).tar $(DISTR) ; gzip $(DISTR).tar
