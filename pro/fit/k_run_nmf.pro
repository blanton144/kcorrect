;+
; NAME:
;   k_run_nmf
; PURPOSE:
;   run the nmf fitting code to get basis spectra
; CALLING SEQUENCE:
;   k_run_nmf
; REVISION HISTORY:
;   29-Nov-2004  Michael Blanton (NYU)
;-
;------------------------------------------------------------------------------
pro k_run_nmf

mmatrix=mrdfits('k_nmf_mmatrix.fits')
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

ilez=where(data.val le 0., nlez)
if(nlez gt 0) then begin
    data.val[ilez]=1.
    data_ivar.val[ilez]=0.
endif

if(file_test('k_nmf_soln.fits')) then begin
    templates=mrdfits('k_nmf_soln.fits',0)
    coeffs=mrdfits('k_nmf_soln.fits',1)
endif 
nmf_sparse, data, data_ivar, 6, mmatrix, 5000L, coeffs=coeffs, $
  templates=templates

mwrfits, templates, 'k_nmf_soln.fits', /create
mwrfits, coeffs, 'k_nmf_soln.fits'

end
