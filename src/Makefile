###############################################################################
# Sloan Digital Sky Survey (SDSS)
# Routines for calculating K-corrections
#
# M.R. Blanton
###############################################################################

SHELL = /bin/sh
#
CFLAGS  = -I$(IDL_DIR)/external $(SDSS_CFLAGS) -DCHECK_LEAKS
#
.c.o : kcorrect.h export.h
	$(CC) -g -c $(CCCHK) $(CFLAGS) $(X_CFLAGS) $(IDL_FLAG) $*.c
.f.o :
	f77 -ffixed-line-length-none -c -o $*.o $*.f

LIB = $(KCORRECT_DIR)/lib
BIN = $(KCORRECT_DIR)/bin

RECIPES = 
OBJECTS = \
	k_fit_coeffs.o \
	idl_k_fit_coeffs.o \
  k_create_r.o \
  idl_k_create_r.o \
	k_fit_templates.o \
	idl_k_fit_templates.o \
	k_model_fluxes.o \
	idl_k_model_fluxes.o \
	k_fit_photoz.o \
	idl_k_fit_photoz.o \
	k_fileopen.o \
	k_load_ascii_table.o \
	k_write_ascii_table.o \
	k_load_filters.o \
	k_interpolate.o \
	k_utils.o \
	k_qromo.o \
	k_zbrent.o \
	k_brent.o \
	k_polint.o \
	k_midpnt.o \
	k_locate.o \
	k_ilaenv.o \
	k_lsame.o \
	k_dposv.o  \
	k_ddot.o    k_dpotf2.o  \
	k_dgemm.o   k_dpotrf.o  k_dtrsm.o \
	k_dgemv.o   k_dpotrs.o  k_dsyrk.o \
	k_dscal.o   k_xerbla.o 

all : $(LIB)/libkcorrect.$(SO_EXT) $(BIN)/fit_coeffs $(BIN)/fit_photoz

$(BIN)/fit_coeffs: $(OBJECTS) $(LIB)/libkcorrect.$(SO_EXT) fit_coeffs.c 
	$(CC) -g -o $(BIN)/fit_coeffs fit_coeffs.c -lm -L$(LIB) -lkcorrect

$(BIN)/fit_photoz: $(OBJECTS) $(LIB)/libkcorrect.$(SO_EXT) fit_photoz.c 
	$(CC) -g -o $(BIN)/fit_photoz fit_photoz.c -lm -L$(LIB) -lkcorrect

$(BIN)/model_fluxes: $(OBJECTS) $(LIB)/libkcorrect.$(SO_EXT) model_fluxes.c 
	$(CC) -g -o $(BIN)/model_fluxes model_fluxes.c -lm -L$(LIB) -lkcorrect

$(LIB)/libkcorrect.$(SO_EXT): $(OBJECTS) 
	$(LD) -g -shared $(X_LD_FLAGS) -o $(LIB)/libkcorrect.$(SO_EXT) $(OBJECTS) \
		$(MAKE_FTNLIB) -L/usr/lib/gcc-lib/i386-redhat-linux/2.96
#	nm -s ../lib/libkcorrect.$(SO_EXT)

#

#
# Install things in their proper places in $(BASEDIR)
#
install : 

clean :
	- /bin/rm -f *~ core *.o so_locations
