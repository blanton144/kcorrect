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

dfile=getenv('KCORRECT_DIR')+'/data/templates/k_nmf_derived.'+vname+'.fits'

age=mrdfits(dfile, 14)
dage=mrdfits(dfile, 15)
sfh_tot=mrdfits(dfile, 12)
sfh_met=mrdfits(dfile, 13)

sfr=sfh_tot#coeffs
metallicity=((sfh_tot*sfh_met)#coeffs)/sfr

end
