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
.c.o :
	$(CC) -c $(CCCHK) $(CFLAGS) $*.c
#
CFLAGS  = $(SDSS_CFLAGS) -DCHECK_LEAKS -I../include

SUBDIRS = pro src 

all :
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) all ); \
	done

#
# Install things in their proper places in $(IDLUTILS_DIR)
#
install :
	- cp *.so $(IDLUTILS_DIR)/lib
	- cp *.a $(IDLUTILS_DIR)/lib
	- cp Makefile $(IDLUTILS_DIR)/lib

clean :
	- /bin/rm -f *~ core
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) clean ); \
	done
