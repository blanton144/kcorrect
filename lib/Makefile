###############################################################################
# Sloan Digital Sky Survey (SDSS)
#
# M. Blanton
###############################################################################

#
# IDL support utilities for spectro2d and the fibermapper
#
SHELL = /bin/sh
#
# Install things in their proper places in $(KCORRECT_DIR)
#
all:

install :
	- cp *.dylib $(KCORRECT_DIR)/lib
	- cp *.so $(KCORRECT_DIR)/lib
	- cp *.a $(KCORRECT_DIR)/lib
	- cp Makefile $(KCORRECT_DIR)/lib
#
clean :
	- /bin/rm -f *~ core *.$(SO_EXT)
