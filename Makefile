#-----------------------------------------------------------------------------
# Top level makefile for Bundler
#
# Bundler: Structure from Motion for Unordered Photo Collections
# Version: 0.1 08/10/2008
#    http://phototour.cs.washington.edu/bundler/
#-----------------------------------------------------------------------------
# Copyright (c) 2008 University of Washington, Cornell University 
#     and Noah Snavely
# All Rights Reserved.
#-----------------------------------------------------------------------------

default:
# Make libraries
	cd lib/ann_L1_0.2; $(MAKE) alpha-g++ #linux-g++-shared
	cd lib/imagelib; $(MAKE)
	cd lib/matrix; $(MAKE)
	cd lib/minpack; $(MAKE)
# Make program
	cd src; $(MAKE)


clean:
	cd lib/ann_L1_0.2; $(MAKE) clean
	cd lib/imagelib; $(MAKE) clean
	cd lib/matrix; $(MAKE) clean
	cd lib/minpack; $(MAKE) clean
	cd src; $(MAKE) clean
