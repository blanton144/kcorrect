###############################################################################
# Sloan Digital Sky Survey (SDSS) -- Value-added galaxy catalog
# M. Blanton & D. Hogg
###############################################################################

SHELL = /bin/sh
#
SUBDIRS = fit seds test utils

all :;

#
# Install things in their proper places in $(VAGC_DIR)
#
install :
	@ for f in $(SUBDIRS); do \
		(mkdir $(VAGC_DIR)/pro/$$f; \
		echo Copying $$f; \
		cp $$f/* $(VAGC_DIR)/pro/$$f ); \
	done
	-@ cp Makefile $(VAGC_DIR)/pro

clean :
	- /bin/rm -f *~ core
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; /bin/rm -f core *~ ); \
	done
