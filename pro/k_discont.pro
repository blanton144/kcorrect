;+
; NAME:
;   k_discont
;
; PURPOSE:
;   Calculate the mean coefficients in two redshift regimes.
;   Calculate the mean colors as well.
;   Now reconstruct the colors in each regime from the other regime.
;   Find the set of color offsets which minimizes the color
;   differences between the mean and reconsructed colors in both
;   regimes simultaneously.
;
; CALLING SEQUENCE:
;   k_discont,
;      
; INPUTS:
;
; OPTIONAL INPUTS:
;   version       - version of templates to use (default 'default')
;   vpath   - path to templates (default $KCORRECT_DIR/data/etemplates)
;
; OUTPUTS:
;
; OPTIONAL INPUT/OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   28-Jan-2002  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_discont, galaxy_mag, galaxy_magerr, galaxy_z, zregions, offsets, version=version, vpath=vpath, rmatrix=rmatrix, zvals=zvals, ematrix=ematrix

; Need at least 6 parameters
if (N_params() LT 4) then begin
    print, 'Syntax - k_discont, galaxy_mag, galaxy_magerr, galaxy_z, zregions, $'
    print, '        offsets [, version=, vpath=, rmatrix=, zvals=, ematrix=]'
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

; find two regions
indx1=where(galaxy_z gt zregions[0,0] and galaxy_z lt zregions[1,0])
meanz1=0.5*(zregions[1,0]-zregions[0,0])
indx2=where(galaxy_z gt zregions[0,1] and galaxy_z lt zregions[1,1])
meanz2=0.5*(zregions[1,1]-zregions[0,1])

; Calculate coeffs
galaxy_z1=galaxy_z[indx1]
galaxy_flux1=10.^(-0.4*galaxy_mag[indx1])
galaxy_invvar1=1./(galaxy_flux1*0.4*alog(10.)*galaxy_magerr[indx1])^2
galaxy_z2=galaxy_z[indx2]
galaxy_flux2=10.^(-0.4*galaxy_mag[indx2])
galaxy_invvar2=1./(galaxy_flux2*0.4*alog(10.)*galaxy_magerr[indx2])^2
if(n_elements(rmatrix) gt 0 AND n_elements(zvals) gt 0 AND $
   n_elements(ematrix) gt 0) then begin
    k_fit_coeffs,galaxy_flux1,galaxy_invvar1,galaxy_z1,coeff1, $
      filterpath=filterpath,rmatrix=rmatrix,zvals=zvals, $
      ematrix=ematrix
    k_fit_coeffs,galaxy_flux2,galaxy_invvar2,galaxy_z2,coeff2, $
      filterpath=filterpath,rmatrix=rmatrix,zvals=zvals, $
      ematrix=ematrix
endif else begin
    k_fit_coeffs,galaxy_flux1,galaxy_invvar1,galaxy_z1,coeff1, $
      version=version,vpath=vpath,filterpath=filterpath,rmatrix=rmatrix, $
      zvals=zvals,ematrix=ematrix
    k_fit_coeffs,galaxy_flux2,galaxy_invvar2,galaxy_z2,coeff2, $
      version=version,vpath=vpath,filterpath=filterpath,rmatrix=rmatrix, $
      zvals=zvals,ematrix=ematrix
endelse

; Find mean scaled coeffs
meancoeff1=dblarr(nt)
meancoeff1[0]=1.
for i=1, nt-1 do $
  meancoeff1[i]=djs_avsigclip(coeff1[i,*]/coeff1[0,*])
meancoeff2=dblarr(nt)
meancoeff2[0]=1.
for i=1, nt-1 do $
  meancoeff2[i]=djs_avsigclip(coeff2[i,*]/coeff2[0,*])

; Calculate model fluxes in each regime
k_model_fluxes,meancoeff1,meanz2,recmag2,rmatrix=rmatrix,zvals=zvals, $
  ematrix=ematrix
k_model_fluxes,meancoeff2,meanz1,recmag1,rmatrix=rmatrix,zvals=zvals, $
  ematrix=ematrix
reccolor1=dblarr(nk-1)
for k=0, nk-2 do $
  reccolor1=-2.5*alog10(recmag1[k]/recmag1[k+1])
reccolor2=dblarr(nk-1)
for k=0, nk-2 do $
  reccolor2=-2.5*alog10(recmag2[k]/recmag2[k+1])

; Calculate average colors in each regime
meancolor1=dblarr(nk-1)
for k=0, nk-2 do $
  meancolor1=djs_avsigclip(galaxy_mag[k,indx1]-galaxy_mag[k+1,indx1])
meancolor2=dblarr(nk-1)
for k=0, nk-2 do $
  meancolor2=djs_avsigclip(galaxy_mag[k,indx2]-galaxy_mag[k+1,indx2])
  
klog,reccolor1
klog,meancolor1
klog,reccolor2
klog,meancolor2

end
