###############################################################################
# Sloan Digital Sky Survey (SDSS)
# Routines for calculating K-corrections
#
# M.R. Blanton
###############################################################################

SHELL = /bin/sh
#
BLAH_CFLAGS  = -I$(KCORRECT_DIR)/include \
	-I$(KCORRECT_DIR)/src 
CFLAGS  = -I$(KCORRECT_DIR)/include \
	-I$(KCORRECT_DIR)/src $(X_CFLAGS) $(SDSS_CFLAGS) -DCHECK_LEAKS
#
.c.o : $(KCORRECT_DIR)/include/kcorrect.h export.h
	$(CC) -O -c $(CCCHK) $(CFLAGS) $(IDL_FLAG) $*.c

LIB = $(KCORRECT_DIR)/lib
BIN = $(KCORRECT_DIR)/bin

RECIPES = 
OBJECTS = \
  idl_k_projection_table.o \
  k_projection_table.o \
  idl_k_fit_spec.o \
  k_fit_spec.o \
  idl_k_fit_spec_linear.o \
  k_fit_spec_linear.o \
  idl_k_fit_nonneg.o \
  k_fit_nonneg.o \
  idl_k_fit_photoz.o \
  k_fit_photoz.o \
  idl_k_reconstruct_maggies.o \
  k_reconstruct_maggies.o \
  k_strparse.o \
  k_yanny_readone.o \
  k_filter_struct.o \
  k_fileopen.o \
  k_load_filters.o \
  k_read_ascii_table.o \
  k_nonneg_solve.o \
	k_interpolate.o \
	idl_k_binspec.o \
	k_binspec.o \
	k_utils.o \
	k_qromo.o \
	k_midpnt.o \
	k_choldc.o \
	k_cholsl.o \
	k_locate.o \
	k_brent.o \
	k_evolve.o \
	k_zbrent.o \
	k_polint.o \
	idl_lf_calc_vmax.o \
	idl_ztoV.o \
	idl_Vtoz.o \
	idl_z2dm.o \
	idl_z2t.o \
	idl_t2z.o \
	idl_lf_eep.o \
	idl_lf_select_eep.o \
	ztransform.o \
	iterate_lf.o \
	philike.o \
	phierrors_lf.o \
	gaussj.o \
	lf_select_eep.o \
	lf_calc_vmax.o \
	lf_eep.o \
	lf_eepfit.o \
	lf_WH_interp.o \
	lf_set_AB.o \
	lf_sum_AB.o 

all: \
  $(BIN)/k_test_filter $(BIN)/fit_coeffs \
	$(BIN)/reconstruct_maggies $(BIN)/fit_photoz \
	$(LIB)/libkcorrect.$(SO_EXT) 
#
$(LIB)/libkcorrect.$(SO_EXT): $(OBJECTS)
	$(LD) $(X_LD_FLAGS) -o $(LIB)/libkcorrect.$(SO_EXT) $(OBJECTS) 

$(BIN)/k_test_filter: $(OBJECTS) $(LIB)/libkcorrect.$(SO_EXT) k_test_filter.c 
	$(CC) $(CFLAGS) -o $(BIN)/k_test_filter k_test_filter.c -lm -lkcorrect \
  -L$(LIB) 

$(BIN)/fit_coeffs: $(OBJECTS) $(LIB)/libkcorrect.$(SO_EXT) fit_coeffs.c 
	$(CC) $(CFLAGS) -o $(BIN)/fit_coeffs fit_coeffs.c -lm -L$(LIB) \
	-lkcorrect

$(BIN)/fit_photoz: $(OBJECTS) $(LIB)/libkcorrect.$(SO_EXT) fit_photoz.c 
	$(CC)  $(CFLAGS) -o $(BIN)/fit_photoz fit_photoz.c -lm -L$(LIB) \
	-lkcorrect 

$(BIN)/reconstruct_maggies: $(OBJECTS) $(LIB)/libkcorrect.$(SO_EXT) \
	reconstruct_maggies.c 
	$(CC)  $(CFLAGS) -o $(BIN)/reconstruct_maggies reconstruct_maggies.c \
	-lm -L$(LIB) -lkcorrect

#
#
# Install things in their proper places in $(BASEDIR)
#
install : 
	cp *.c $(KCORRECT_DIR)/src 
	-@ cp Makefile $(KCORRECT_DIR)/src

clean :
	- /bin/rm -f *~ core *.o so_locations $(BIN)/fit_coeffs $(BIN)/k_test_filter $(BIN)/fit_photoz $(BIN)/reconstruct_maggies $(LIB)/libkcorrect.$(SO_EXT) 
