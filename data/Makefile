###############################################################################
# kcorrect
# M. Blanton
###############################################################################

SHELL = /bin/sh
#
all :;
#
# Install things in their proper places in $(VAGC_DIR)
#
install :
	- cp -r basel $(VAGC_DIR)/data
	- cp -r dustmodels $(VAGC_DIR)/data
	- cp -r filters $(VAGC_DIR)/data
	- cp -r redshifts $(VAGC_DIR)/data
	- cp -r seds $(VAGC_DIR)/data
	- cp -r templates $(VAGC_DIR)/data

clean :
	- /bin/rm -f *~ core
