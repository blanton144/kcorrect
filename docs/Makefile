###############################################################################
# Sloan Digital Sky Survey (SDSS)
# kcorrect
# M. Blanton
###############################################################################

SHELL = /bin/sh
#
SUBDIRS = paper
#
all :
	idl mk_doc
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) all ); \
	done
	cp -f paper/paper.ps html
	ln -sf ../../data/filters html/filters

#
# Install things in their proper places in $(KCORRECT_DIR)
#
install :
	- cp -r mk_doc v_11_msg html paper $(KCORRECT_DIR)/docs

clean :
	- /bin/rm -f *~ core www/*~ www/core


