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
data=mrdfits('k_nmf_data.fits',1)
ivar=mrdfits('k_nmf_data.fits',2)

ilez=where(data le 0., nlez)
if(nlez gt 0) then begin
    data[ilez]=1.
    ivar[ilez]=0.
endif

if(file_test('k_nmf_soln.fits')) then begin
    ww0=mrdfits('k_nmf_soln.fits',0)
    hh0=mrdfits('k_nmf_soln.fits',1)
    nmf_sq_m_err, data, 6, mmatrix, ivar, 5000L, ww0, hh0, ww=ww, hh=hh
endif else begin
    nmf_sq_m_err, data, 6, mmatrix, ivar, 5000L, ww=ww, hh=hh
endelse

mwrfits, ww, 'k_nmf_soln.fits', /create
mwrfits, hh, 'k_nmf_soln.fits'

ww=mrdfits('k_nmf_soln.fits')
hh=mrdfits('k_nmf_soln.fits',1)
model=mmatrix#ww#hh



end
