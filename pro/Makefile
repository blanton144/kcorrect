###############################################################################
# Sloan Digital Sky Survey (SDSS) 
# M. Blanton & D. Hogg
###############################################################################

SHELL = /bin/sh
#
SUBDIRS = fit lf seds test utils samples

all :;

#
# Install things in their proper places in $(KCORRECT_DIR)
#
install :
	@ for f in $(SUBDIRS); do \
		(mkdir $(KCORRECT_DIR)/pro/$$f; \
		echo Copying $$f; \
		cp $$f/* $(KCORRECT_DIR)/pro/$$f ); \
	done
	-@ cp Makefile $(KCORRECT_DIR)/pro

clean :
	- /bin/rm -f *~ core
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; /bin/rm -f core *~ ); \
	done
