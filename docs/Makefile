###############################################################################
# Sloan Digital Sky Survey (SDSS)
# IDL support code for products: kcorrect
#
# S. Burles & D. Schlegel
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
# Install things in their proper places in $(VAGC_DIR)
#
install :

clean :
	- /bin/rm -f *~ core www/*~ www/core


