;+
; NAME:
;   kphotoz
;
; PURPOSE:
;   Given a set of AB magnitudes and errors, guess the redshift.
;   (If no redshift is specified to K-correct to, just return the 
;   set of reconstructed magnitudes based on the templates used). 
;
; CALLING SEQUENCE:
;   kphotoz, galaxy_mag, galaxy_magerr, photoz, [coeffs=, version=, $
;       vpath=, /maggies, rmatrix=, zvals=, ematrix=] 
;      
; INPUTS:
;   galaxy_mag    - AB magnitudes of galaxies [N_band, N_gal]
;   galaxy_magerr - errors in magnitudes  [N_band, N_gal]
;
; OPTIONAL INPUTS:
;   version       - version of templates to use (default 'default')
;   vpath   - path to templates (default $KCORRECT_DIR/data/etemplates)
;   maggies       - set if input and output in 10^{-0.4*mag}
;
; OUTPUTS:
;   photoz        - photometric redshifts
;
; OPTIONAL INPUT/OUTPUTS:
;   ematrix       - eigentemplates [N_dim, N_template]
;   rmatrix       - look up table for bmatrix and filter information 
;                   [N_z, N_dim, N_band]
;   zvals         - look up table for rmatrix [N_z]
;
; COMMENTS:
;   This program has a large amount of overhead. So use it on long
;   lists for best results (ie. don't call this 1,000,000 times to get
;   photo-z's for 1,000,000 objects). 
;
;   Be careful when sending SDSS "photo" outputs directly into this program.
;   Eg. occasionally the magnitudes or the errors have crazy values,
;   such as -9999. The normal "garbage in, garbage out" rules apply.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   k_fit_photoz
;
; REVISION HISTORY:
;   04-Jan-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro kphotoz, galaxy_mag, galaxy_magerr, photoz, coeffs=coeffs, version=version, vpath=vpath, maggies=maggies, rmatrix=rmatrix, zvals=zvals, ematrix=ematrix

; Need at least 6 parameters
if (N_params() LT 3) then begin
    print, 'Syntax - kphotoz, galaxy_mag, galaxy_magerr, photoz [, coeffs=, $'
    print, '         version=, vpath=, /maggies, rmatrix=, zvals=, ematrix=]'
    return
endif

ngalaxy=long((size(galaxy_mag))[2])
nk=long(n_elements(galaxy_mag))/ngalaxy

if(NOT keyword_set(filterpath)) then $
  filterpath=getenv('KCORRECT_DIR')+'/data/filters'
if(NOT keyword_set(vpath)) then $
  vpath=getenv('KCORRECT_DIR')+'/data/etemplates'
if(NOT keyword_set(version)) then $
  version='defaultpz'

; Calculate maggies if necessary
if(NOT keyword_set(maggies)) then begin
    galaxy_flux=10.^(-0.4*galaxy_mag)
    galaxy_invvar=1./(galaxy_flux*0.4*alog(10.)*galaxy_magerr)^2
endif else begin
    galaxy_flux=galaxy_mag
    galaxy_invvar=1./galaxy_magerr^2
endelse 

; Calculate coeffs
if(n_elements(rmatrix) gt 0 AND n_elements(zvals) gt 0 AND $
   n_elements(ematrix) gt 0) then begin
    k_fit_photoz,galaxy_flux,galaxy_invvar,photoz,coeffs, $
      filterpath=filterpath,rmatrix=rmatrix,zvals=zvals, $
      ematrix=ematrix
endif else begin
    k_fit_photoz,galaxy_flux,galaxy_invvar,photoz,coeffs,version=version, $
      vpath=vpath,filterpath=filterpath,rmatrix=rmatrix,zvals=zvals, $
      ematrix=ematrix
endelse

end
