###############################################################################
# kcorrect
#
# M. Blanton
###############################################################################

SHELL = /bin/sh
#
.c.o :
	$(CC) -c $(CCCHK) $(CFLAGS) $*.c
#
#
SUBDIRS = src pro data docs lib test ups include bin

all :
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) all ); \
	done

install: 
	@echo "You should be sure to have updated before doing this."
	@echo ""
	@if [ "$(KCORRECT_DIR)" = "" ]; then \
		echo You have not specified a destination directory >&2; \
		exit 1; \
	fi 
	@if [ -e $(KCORRECT_DIR) ]; then \
		echo The destination directory already exists >&2; \
		exit 1; \
	fi 
	@echo ""
	@echo "You will be installing in \$$KCORRECT_DIR=$$KCORRECT_DIR"
	@echo "I'll give you 5 seconds to think about it"
	@sleep 5
	@echo ""
	@ rm -rf $(KCORRECT_DIR)
	@ mkdir $(KCORRECT_DIR)
	@ for f in $(SUBDIRS); do \
		(mkdir $(KCORRECT_DIR)/$$f; cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) install ); \
	done
	/bin/cp Makefile $(KCORRECT_DIR)

clean :
	- /bin/rm -f *~ core
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) clean ); \
	done
