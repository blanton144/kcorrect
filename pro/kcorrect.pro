;+
; NAME:
;   kcorrect
;
; PURPOSE:
;   Given a set of AB magnitudes, K-correct the magnitudes to a given 
;   redshift. (If no redshift is specified to K-correct to, just
;   return the set of reconstructed magnitudes based on the templates
;   used). 
;
; CALLING SEQUENCE:
;   kcorrect, galaxy_mag, galaxy_magerr, galaxy_z, galaxy_mag0, 
;      [kcorrectz=, version=, vpath=, /maggies] 
;      
; INPUTS:
;   galaxy_mag    - AB magnitudes of galaxies
;   galaxy_magerr - errors in magnitudes 
;   galaxy_z      - redshifts of galaxies
;
; OPTIONAL INPUTS:
;   kcorrectz     - redshift (or redshifts) to K-correct magnitudes to
;   version       - version of templates to use (default 'default')
;   vpath   - path to templates (default $KCORRECT_DIR/data/etemplates)
;   maggies       - set if input and output in 10^{-0.4*mag}
;
; OUTPUTS:
;   galaxy_mag0   - K-corrected AB magnitudes
;
; OPTIONAL INPUT/OUTPUTS:
;
; COMMENTS:
;   This program has a large amount of overhead. So use it on long
;   lists for best results (ie. don't call this 1,000,000 times to get
;   K-corrections for 1,000,000 objects). If it is *unavoidable* to
;   do this, then you should figure out how to use k_fit_coeffs and
;   k_model_fluxes to this more efficiently (by calculating rmatrix
;   only once). You can ask me for help if this is too hard. Or you
;   can convince me to let you do this through kcorrect.
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
;   k_fit_coeffs
;   k_model_fluxes
;
; REVISION HISTORY:
;   24-Jan-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro kcorrect, galaxy_mag, galaxy_magerr, galaxy_z, galaxy_mag0, kcorrectz=kcorrectz, version=version, vpath=vpath, maggies=maggies

; Need at least 6 parameters
if (N_params() LT 4) then begin
    print, 'Syntax - kcorrect, galaxy_mag, galaxy_magerr, galaxy_z, galaxy_mag0, $'
    print, '        [kcorrectz=, version=, vpath=, maggies=]'
    return
endif

ngalaxy=long(n_elements(galaxy_z))
nk=long(n_elements(galaxy_mag))/ngalaxy

if(NOT keyword_set(filterpath)) then $
  filterpath=getenv('KCORRECT_DIR')+'/data/filters'
if(NOT keyword_set(vpath)) then $
  vpath=getenv('KCORRECT_DIR')+'/data/etemplates'
if(NOT keyword_set(version)) then $
  version='default'

; Calculate maggies if necessary
if(NOT keyword_set(maggies)) then begin
    galaxy_flux=10.^(-0.4*galaxy_mag)
    galaxy_invvar=1./(galaxy_flux*0.4*alog(10.)*galaxy_magerr)^2
endif else begin
    galaxy_flux=galaxy_mag
    galaxy_invvar=1./galaxy_magerr^2
endelse 

; Calculate coeffs
k_fit_coeffs,galaxy_flux,galaxy_invvar,galaxy_z,coeff,version=version, $
  vpath=vpath,filterpath=filterpath,rmatrix=rmatrix,zvals=zvals, $
  ematrix=ematrix

; Calculate model fluxes
if(n_elements(kcorrectz) eq 0) then begin 
    correct_z=galaxy_z
endif else begin
    if(n_elements(kcorrectz) eq 1) then begin
        correct_z=replicate(kcorrectz,n_elements(galaxy_z))
    endif else begin
        correct_z=kcorrectz
    endelse
endelse
indx=where(correct_z lt zvals[0],count)
if(count gt 0) then correct_z[indx]=0.5*(zvals[0]+zvals[1])
k_model_fluxes,coeff,correct_z,model_flux,rmatrix=rmatrix,zvals=zvals, $
  ematrix=ematrix

; Calculate magnitudes if necessary
galaxy_mag0=dblarr(nk,ngalaxy)
if(NOT keyword_set(maggies)) then begin
    negindx=where(model_flux le 0.d,count)
    if(count gt 0) then galaxy_mag0[negindx]=0.d
    goodindx=where(model_flux gt 0.d,count)
    if(count gt 0) then $
      galaxy_mag0[goodindx]=-2.5*alog10(model_flux[goodindx])
endif else begin
    galaxy_mag0=model_flux
endelse

end
