;+
; NAME:
;   sfh_reconstruct
; PURPOSE:
;   get SFH parameters given coefficients
; CALLING SEQUENCE:
;   sfh_reconstruct, coeffs [, vname=, sfr=, dage=, metallicity=, age= ]
; INPUTS:
;   coeffs - [Nt, Ngals] coefficients for the galaxy in question
; OPTIONAL INPUTS:
;   vname - name of set of templates to use (default 'default')
; OUTPUTS:
;   age - [Nage] ages of bursts in model (yrs)
;   dage - [Nage] width of age bin (yrs)
;   sfr - [Nage, Ngals] rate of star formation in bin (M_solar/yr)
;   metallicity - [Nage, Ngals] metallicity in each bin
; REVISION HISTORY:
;   2005-08-15 MRB, NYU
;-
;------------------------------------------------------------------------------
pro sfh_reconstruct, coeffs, vname=vname, sfr=sfr, metallicity=metallicity, $
                     age=age, dage=dage

if (NOT keyword_set(vname)) then $
  vname='default'

metallicities=[0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05]
hdr=headfits(getenv('KCORRECT_DIR')+ $
             '/data/templates/k_nmf_mmatrix.'+vname+'.fits')
dust=mrdfits(getenv('KCORRECT_DIR')+ $
             '/data/templates/k_nmf_mmatrix.'+vname+'.fits',2)
met=mrdfits(getenv('KCORRECT_DIR')+ $
             '/data/templates/k_nmf_mmatrix.'+vname+'.fits',3)
age=mrdfits(getenv('KCORRECT_DIR')+ $
             '/data/templates/k_nmf_mmatrix.'+vname+'.fits',4)
ndusts=long(sxpar(hdr, 'NDUST'))
nmets=long(sxpar(hdr, 'NMET'))
nages=long(sxpar(hdr, 'NAGE'))
templates=mrdfits(getenv('KCORRECT_DIR')+ $
                  '/data/templates/k_nmf_soln.'+vname+'.fits',0)
if((size(templates))[0] eq 1) then $
  nt=1 $
else $
  nt=(size(templates, /dim))[1]

dage=fltarr(nages)
dage[1:nages-2L]=0.5*(age[2:nages-1L]-age[0:nages-3L])
dage[0]=age[1]-age[0]
dage[nages-1L]=age[nages-1L]-age[nages-2L]
sfh_tot=fltarr(nages, nt)
sfh_met=fltarr(nages, nt)
for i=0L, nt-1L do begin
    currsfh=reform(templates[0:nages*nmets*ndusts-1L,i], nages, nmets, ndusts)
    dust=reform(dust, nages, nmets, ndusts)
    for j=0L, nages-1L do $
      sfh_tot[j,i]=total(currsfh[j,*,*])/dage[j]
    for j=0L, nages-1L do $
      sfh_met[j,i]=total(currsfh[j,*,*]*metallicities[met[j,*,*]])/dage[j]
endfor
sfh_met=sfh_met/sfh_tot

sfr=sfh_tot#coeffs
metallicity=((sfh_tot*sfh_met)#coeffs)/sfr

end
