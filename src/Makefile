###############################################################################
# Sloan Digital Sky Survey (SDSS)
# Routines for calculating K-corrections
#
# M.R. Blanton
###############################################################################

SHELL = /bin/sh
#
CFLAGS  = -I$(IDL_DIR)/external -I$(KCORRECT_DIR)/include \
	$(X_CFLAGS) $(SDSS_CFLAGS) -DCHECK_LEAKS
#
.c.o : $(KCORRECT_DIR)/include/kcorrect.h export.h
	$(CC) -c $(CCCHK) $(CFLAGS) $(IDL_FLAG) $*.c
.f.o :
	f77 $(X_FFLAGS) -c -o $*.o $*.f

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
	k_reconstruct_maggies.o \
	idl_k_reconstruct_maggies.o \
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

all : $(LIB)/libkcorrect.$(SO_EXT) $(BIN)/fit_coeffs $(BIN)/fit_photoz \
	$(BIN)/reconstruct_maggies $(BIN)/color_interp

$(BIN)/fit_coeffs: $(OBJECTS) $(LIB)/libkcorrect.$(SO_EXT) fit_coeffs.c 
	$(CC) $(CFLAGS) -o $(BIN)/fit_coeffs fit_coeffs.c -lm -L$(LIB) -lkcorrect

$(BIN)/fit_photoz: $(OBJECTS) $(LIB)/libkcorrect.$(SO_EXT) fit_photoz.c 
	$(CC) $(CFLAGS) -o $(BIN)/fit_photoz fit_photoz.c -lm -L$(LIB) -lkcorrect

$(BIN)/color_interp: color_interp.c $(LIB)/libkcorrect.$(SO_EXT) 
	$(CC) $(CFLAGS) -o $(BIN)/color_interp color_interp.c -lm -L$(LIB) -lkcorrect

$(BIN)/reconstruct_maggies: $(OBJECTS) $(LIB)/libkcorrect.$(SO_EXT) \
		reconstruct_maggies.c 
	$(CC) $(CFLAGS) -o $(BIN)/reconstruct_maggies reconstruct_maggies.c -lm -L$(LIB) -lkcorrect

$(LIB)/libkcorrect.$(SO_EXT): $(OBJECTS) 
	$(LD) -shared $(X_LD_FLAGS) -o $(LIB)/libkcorrect.$(SO_EXT) $(OBJECTS) \
		$(MAKE_FTNLIB) -L/usr/lib/gcc-lib/i386-redhat-linux/2.96
#	nm -s ../lib/libkcorrect.$(SO_EXT)

#

#
# Install things in their proper places in $(BASEDIR)
#
install : 

clean :
	- /bin/rm -f *~ core *.o so_locations
