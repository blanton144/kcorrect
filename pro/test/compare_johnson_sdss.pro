;+
; NAME:
;   compare_johnson_sdss
; PURPOSE:
;   given input sdss galaxies, return rest frame johnson and sdss mags
; CALLING SEQUENCE:
;   compare_johnson_sdss, maggies, maggies_ivar, z, johnson, sdss [ $ 
;           , sdss_band_shift= ]
; INPUTS:
;   maggies - [5, N] input maggeis
;   maggies_ivar - [5, N] input maggeis inverse variance
;   z - [N] input redshifts
; OPTIONAL INPUTS:
;   sdss_band_shift    - blueshift of bandpasses to apply for sdss 
;                        (to get ^{z}b type bands) [default 0]
; OUTPUTS:
;   johnson - [5,N] UBVRI mags (rest)
;   sdss - [5,N] ugriz mags (rest)
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   02-Jun-2003  Updated to new v3_0 standards MRB, NYU
;-
;------------------------------------------------------------------------------
pro compare_johnson_sdss, in_maggies, in_maggies_ivar, extinction, redshift, $
                          johnson, sdss, sdss_model, kcorrect=kcorrect, $
                          sdss_band_shift=sdss_band_shift, chi2=chi2, $
                          ab=ab

if(NOT keyword_set(ab)) then $
  maggies=sdssflux2ab(in_maggies) $
else $
  maggies=in_maggies
if(NOT keyword_set(ab)) then $
  maggies_ivar=sdssflux2ab(in_maggies_ivar, /ivar) $
else $
  maggies_ivar=in_maggies_ivar
maggies=maggies*10.^(0.4*extinction)
maggies_ivar=maggies_ivar*10.^(-0.8*extinction)
ngalaxy=long(n_elements(redshift))
kcorrect, maggies, maggies_ivar, redshift, kcorrect, coeffs=coeffs, $
  band_shift=sdss_band_shift, chi2=chi2, filterlist=filterlist
sdss=22.5-2.5*alog10(maggies)-kcorrect
k_load_vmatrix, vmatrix, lambda, vfile=vfile, lfile=lfile,vpath=vpath
k_reconstruct_maggies, coeffs, fltarr(ngalaxy), smaggies, $
  filterlist=filterlist, lambda=lambda, vmatrix=vmatrix, $
  band_shift=sdss_band_shift

k_reconstruct_maggies, coeffs, fltarr(ngalaxy), jmaggies, vmatrix=vmatrix, $
  filterlist=['bessell_U.par', 'bessell_B.par', 'bessell_V.par', $
              'bessell_R.par', 'bessell_I.par'], lambda=lambda
johnson=22.5-2.5*alog10(jmaggies)
sdss_model=22.5-2.5*alog10(smaggies)
jab2vega=-k_vega2ab(filterlist=['bessell_U.par', 'bessell_B.par', $
                                'bessell_V.par', 'bessell_R.par', $
                                'bessell_I.par'], /hayes)
for i=0, 4 do $
  johnson[i,*]=johnson[i,*]+jab2vega[i]

end
