###############################################################################
# Sloan Digital Sky Survey (SDSS)
#
# M. Blanton
###############################################################################

SHELL = /bin/sh
#
.c.o :
	$(CC) -c $(CCCHK) $(CFLAGS) $*.c
#
#

SUBDIRS = pro src docs

all :
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) all ); \
	done

#
# Install things in their proper places in $(IDLUTILS_DIR)
#
install : lib/libkcorrect.so
	cp lib/libkcorrect.so $KCORRECT_LD_LIB

clean :
	- /bin/rm -f *~ core
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) clean ); \
	done
