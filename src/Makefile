###############################################################################
# Sloan Digital Sky Survey (SDSS)
# Routines for calculating K-corrections
#
# M.R. Blanton
###############################################################################

SHELL = /bin/sh
#
CFLAGS  = -I$(IDL_DIR)/external -I$(KCORRECT_DIR)/include \
	-I$(KCORRECT_DIR)/src $(X_CFLAGS) $(SDSS_CFLAGS) -DCHECK_LEAKS
#
.c.o : $(KCORRECT_DIR)/include/kcorrect.h export.h
	$(CC) -c $(CCCHK) $(CFLAGS) $(IDL_FLAG) $*.c
.f.o :
	$(FC) $(X_FFLAGS) -c -o $*.o $*.f

LIB = $(KCORRECT_DIR)/lib
BIN = $(KCORRECT_DIR)/bin

RECIPES = 
OBJECTS = \
  idl_k_projection_table.o \
  k_projection_table.o \
  idl_k_fit_nonneg.o \
  k_fit_nonneg.o \
  idl_k_fit_photoz.o \
  k_fit_photoz.o \
  idl_k_reconstruct_maggies.o \
  k_reconstruct_maggies.o \
  k_nonneg_solve.o \
	k_interpolate.o \
	k_utils.o \
	k_qromo.o \
	k_midpnt.o \
	k_choldc.o \
	k_cholsl.o \
	k_locate.o \
	k_brent.o \
	k_polint.o 

$(LIB)/libkcorrect.$(SO_EXT): $(OBJECTS) 
	$(LD) -shared $(X_LD_FLAGS) -o $(LIB)/libkcorrect.$(SO_EXT) $(OBJECTS) 
#	nm -s ../lib/libkcorrect.$(SO_EXT)

#

#
# Install things in their proper places in $(BASEDIR)
#
install : 

clean :
	- /bin/rm -f *~ core *.o so_locations
