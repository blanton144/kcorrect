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

SUBDIRS = src pro data docs lib test src ups include

install:
	@echo "You should be sure to have updated before doing this."
	@echo ""
	@if [ "$(VAGC_DIR)" = "" ]; then \
		echo You have not specified a destination directory >&2; \
		exit 1; \
	fi 
	@if [ -e $(VAGC_DIR) ]; then \
		echo The destination directory already exists >&2; \
		exit 1; \
	fi 
	@echo ""
	@echo "You will be installing in \$$VAGC_DIR=$$VAGC_DIR"
	@echo "I'll give you 5 seconds to think about it"
	@sleep 5
	@echo ""
	@ rm -rf $(VAGC_DIR)
	@ mkdir $(VAGC_DIR)
	@ for f in $(SUBDIRS); do \
		(mkdir $(VAGC_DIR)/$$f; cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) install ); \
	done
	/bin/cp Makefile $(VAGC_DIR)

all :
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) all ); \
	done

clean :
	- /bin/rm -f *~ core
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) clean ); \
	done
