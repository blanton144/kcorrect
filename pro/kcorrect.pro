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
;      [kcorrectz=, version=, vpath=, /maggies, rmatrix=, zvals=,
;       ematrix=] 
;      
; INPUTS:
;   galaxy_mag    - AB magnitudes of galaxies [N_band, N_gal]
;   galaxy_magerr - errors in magnitudes  [N_band, N_gal]
;   galaxy_z      - redshifts of galaxies [N_gal]
;
; OPTIONAL INPUTS:
;   kcorrectz     - redshift (or redshifts) to K-correct magnitudes to
;   version       - version of templates to use (default 'default')
;   vpath   - path to templates (default $KCORRECT_DIR/data/etemplates)
;   maggies       - set if input and output in 10^{-0.4*mag}
;   invvar        - if maggies is set, this means that magerr is
;                   actually invvar
;   sdssfix       - uses k_sdssfix to "fix" the SDSS magnitudes
;   vconstraint   - use the standard constraints for the version 
;   constraints_amp - amplitude of user-specified constraint
;   constraints_mean - mean of user-specified constraint
;   constraints_var - variance matrix of user-specified constraint
;
; OUTPUTS:
;   galaxy_mag0   - K-corrected AB magnitudes
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
;   K-corrections for 1,000,000 objects). If it is *unavoidable* to
;   do this, then you should figure out how to use k_fit_coeffs and
;   k_reconstruct_maggies to this more efficiently (by calculating rmatrix
;   only once). You can ask me for help if this is too hard. Or you
;   can convince me to let you do this through kcorrect.
;
;   Be careful when sending SDSS "photo" outputs directly into this program.
;   Eg. occasionally the magnitudes or the errors have crazy values,
;   such as -9999. The normal "garbage in, garbage out" rules
;   apply. However, I have supplied the /sdssfix flag, which does a
;   reasonable job in most cases of identifying problem cases and
;   doing something OK about it. If you really care about the u
;   magnitudes of every single object, you will have to think harder,
;   but if all you want is to K-correct gri magnitudes well, this 
;   will probably do just fine. /sdssfix also adds photometric
;   zeropoint errors (0.05,0.02,0.02,0.02,0.03) to the fit.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   k_fit_coeffs
;   k_reconstruct_maggies
;
; REVISION HISTORY:
;   24-Jan-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro kcorrect, galaxy_mag, galaxy_magerr, galaxy_z, galaxy_mag0, kcorrectz=kcorrectz, version=version, vpath=vpath, maggies=maggies, rmatrix=rmatrix, zvals=zvals, ematrix=ematrix, coeff=coeff, sdssfix=sdssfix, addgrgap=addgrgap, invvar=invvar, vconstraint=vconstraint, constraints_amp=constraints_amp, constraints_mean=constraints_mean, constraints_var=constraints_var

; Need at least 6 parameters
if (N_params() LT 4) then begin
    print, 'Syntax - kcorrect, galaxy_mag, galaxy_magerr, galaxy_z, galaxy_mag0, $'
    print, '        [kcorrectz=, version=, vpath=, /maggies, rmatrix=, zvals=, ematrix=,$'
    print, '         coeff=, /sdssfix, /invvar, /vconstraint]'
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

; Fix SDSS mags if desired
tmp_galaxy_mag=galaxy_mag
tmp_galaxy_magerr=galaxy_magerr
if(keyword_set(sdssfix)) then begin
   k_sdssfix,tmp_galaxy_mag,tmp_galaxy_magerr,galaxy_z=galaxy_z, $
      maggies=maggies, invvar=invvar, /errorsonly
endif

; set constraints for this version if necessary
if(keyword_set(vconstraint)) then begin
   constraints_amp=1.d
   k_load_ascii_table,constraints_var,vpath+'/scaledvar.'+version+'.dat'
   k_load_ascii_table,constraints_mean,vpath+'/scaledmean.'+version+'.dat'
endif
  
; Calculate maggies if necessary
if(NOT keyword_set(maggies)) then begin
    galaxy_maggies=10.^(-0.4*tmp_galaxy_mag)
    galaxy_invvar=1./(galaxy_maggies*0.4*alog(10.)*tmp_galaxy_magerr)^2
endif else begin
    galaxy_maggies=tmp_galaxy_mag
    if(NOT keyword_set(invvar)) then $
      galaxy_invvar=1./tmp_galaxy_magerr^2
    if(keyword_set(invvar)) then $
      galaxy_invvar=tmp_galaxy_magerr
endelse 
tmp_galaxy_mag=0l
tmp_galaxy_magerr=0l

; Set the grgap maggie if appropriate
useversion=version
if(keyword_set(addgrgap)) then begin
    useversion='addgrgap'
    tmp_galaxy_maggies=dblarr(nk+1,ngalaxy)
    tmp_galaxy_invvar=dblarr(nk+1,ngalaxy)
    tmp_galaxy_maggies[0:nk-1,*]=galaxy_maggies
    tmp_galaxy_invvar[0:nk-1,*]=galaxy_invvar
    tmp_galaxy_maggies[nk,*]=0.5*(galaxy_maggies[1,*]+galaxy_maggies[2,*])
    tmp_galaxy_invvar[nk,*]= $
      1.d/(0.25*tmp_galaxy_maggies[nk,*])^2
    galaxy_maggies=tmp_galaxy_maggies
    galaxy_invvar=tmp_galaxy_invvar
    tmp_galaxy_maggies=0l
    tmp_galaxy_invvar=0l
endif

; Calculate coeffs
if(n_elements(rmatrix) gt 0 AND n_elements(zvals) gt 0 AND $
   n_elements(ematrix) gt 0) then begin
    k_fit_coeffs,galaxy_maggies,galaxy_invvar,galaxy_z,coeff, $
      filterpath=filterpath,rmatrix=rmatrix,zvals=zvals, $
      ematrix=ematrix, constraints_amp=constraints_amp, $
      constraints_var=constraints_var, constraints_mean=constraints_mean
endif else begin
    k_fit_coeffs,galaxy_maggies,galaxy_invvar,galaxy_z,coeff, $
      version=useversion,vpath=vpath,filterpath=filterpath,rmatrix=rmatrix, $
      zvals=zvals,ematrix=ematrix, constraints_amp=constraints_amp, $
      constraints_var=constraints_var, constraints_mean=constraints_mean
endelse
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
k_reconstruct_maggies,coeff,correct_z,reconstruct_maggies,rmatrix=rmatrix, $
  zvals=zvals,ematrix=ematrix

; Remove the grgap maggie if appropriate
if(keyword_set(addgrgap)) then begin
    tmp_reconstruct_maggies=dblarr(nk,ngalaxy)
    tmp_reconstruct_maggies=reconstruct_maggies[0:nk-1,*]
    reconstruct_maggies=tmp_reconstruct_maggies
    tmp_reconstruct_maggies=0l
endif

; Calculate magnitudes if necessary
galaxy_mag0=dblarr(nk,ngalaxy)
if(NOT keyword_set(maggies)) then begin
    negindx=where(reconstruct_maggies le 0.d,count)
    if(count gt 0) then galaxy_mag0[negindx]=1000.d
    goodindx=where(reconstruct_maggies gt 0.d,count)
    if(count gt 0) then $
      galaxy_mag0[goodindx]=-2.5*alog10(reconstruct_maggies[goodindx])
endif else begin
    galaxy_mag0=reconstruct_maggies
endelse

end
