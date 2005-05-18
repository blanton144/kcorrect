;+
; NAME:
;   k_run_nmf
; PURPOSE:
;   run the nmf fitting code to get basis spectra
; CALLING SEQUENCE:
;   k_run_nmf [, nt=, niter= ]
; OPTIONAL INPUTS:
;   nt - number of templates to fit for (default 6)
;   niter - number of iterations of fit (default 1000)
; OPTIONAL KEYWORDS:
;   /qa - make qa plots at end
;   /reset - ignore k_nmf_soln file if it exists
; COMMENTS:
;   Requires k_nmf_mmatrix and k_nmf_spdata to have been run. 
;   Uses nmf_sparse for fitting. Puts results in k_nmf_soln.fits.
;   Checks for existing k_nmf_soln.fits to use as starting
;   point. Otherwise nmf_sparse chooses random starting point.
; REVISION HISTORY:
;   29-Nov-2004  Michael Blanton (NYU)
;-
;------------------------------------------------------------------------------
pro k_run_nmf, niter=niter, nt=nt, qa=qa, reset=reset

if(NOT keyword_set(niter)) then niter=1000L
if(NOT keyword_set(nt)) then nt=4

mmatrix=mrdfits('k_nmf_mmatrix.fits')
dust=mrdfits('k_nmf_mmatrix.fits',2)
age=mrdfits('k_nmf_mmatrix.fits',4)
datastr=mrdfits('k_nmf_spdata.fits',1)
vals=mrdfits('k_nmf_spdata.fits',2)
ivar=mrdfits('k_nmf_spdata.fits',3)
xx=mrdfits('k_nmf_spdata.fits',4)

data=create_struct(datastr, $
                   'val', fltarr(n_elements(vals)), $
                   'x', fltarr(n_elements(vals)))
data.val=vals
data.x=xx
data_ivar=create_struct(datastr, $
                        'val', fltarr(n_elements(vals)), $
                        'x', fltarr(n_elements(vals)))
data_ivar.val=ivar
data_ivar.x=xx
ngals=n_elements(data.rowstart)

ilez=where(data.val le 0. or data_ivar.val gt 1.d+30, nlez)
if(nlez gt 0) then begin
    data.val[ilez]=1.
    data_ivar.val[ilez]=0.
endif

if(file_test('k_nmf_soln.fits') eq 1 and $
   keyword_set(reset) eq 0) then begin
    templates=mrdfits('k_nmf_soln.fits',0)
    coeffs=mrdfits('k_nmf_soln.fits',1)
    if(keyword_set(coeffs)) then $
      if((size(coeffs,/dim))[1] ne ngals) then coeffs=0
endif else begin
    templates=0.5+randomu(seed, (size(mmatrix, /dim))[1], nt)
    for i=0L, nt-1L do $
      templates[0:n_elements(age)-1L,i]= $
      max(sqrt(age))*templates[0:n_elements(age)-1L,i]/sqrt(age)
    ii=where(dust.tauv gt 0., nii)
    if(nii gt 0) then $
      for i=0L, nt-1L do $
      templates[ii, i]=templates[ii,i]*0.001
endelse
nmf_sparse, data, data_ivar, nt, mmatrix, niter, coeffs=coeffs, $
  templates=templates

mwrfits, templates, 'k_nmf_soln.fits', /create
mwrfits, coeffs, 'k_nmf_soln.fits'

if(keyword_set(qa)) then $
  k_qa_nmf

end
