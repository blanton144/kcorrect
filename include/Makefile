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
install :
	- cp *.h $(KCORRECT_DIR)/include
	- cp Makefile $(KCORRECT_DIR)/include

clean :
	- /bin/rm -f *~ core
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) clean ); \
	done
