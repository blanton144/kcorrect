###############################################################################
# kcorrect
# M. Blanton
###############################################################################

SHELL = /bin/sh
#
all :;
#
# Install things in their proper places in $(KCORRECT_DIR)
#
install :
	- cp -r basel $(KCORRECT_DIR)/data
	- cp -r dustmodels $(KCORRECT_DIR)/data
	- cp -r filters $(KCORRECT_DIR)/data
	- cp -r redshifts $(KCORRECT_DIR)/data
	- cp -r seds $(KCORRECT_DIR)/data
	- cp -r templates $(KCORRECT_DIR)/data
	- cp -r test $(KCORRECT_DIR)/data

clean :
	- /bin/rm -f *~ core
