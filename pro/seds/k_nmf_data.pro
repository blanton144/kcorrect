;+
; NAME:
;   k_nmf_data
; PURPOSE:
;   put together data for nmf fitting
; CALLING SEQUENCE:
;   k_nmf_data
; BUGS:
;   AB corrections need updated
;   what about filter curves?
;   smooth spectra to constant vdisp?
; REVISION HISTORY:
;   23-Nov-2004  Michael Blanton (NYU)
;-
;------------------------------------------------------------------------------
pro k_nmf_data, mmatrix=mmatrix, sample=sample

if(NOT keyword_set(mmatrix)) then mmatrix='k_nmf_mmatrix.fits'
if(NOT keyword_set(sample)) then sample='drtwo14'
if(NOT keyword_set(nsub)) then nsub=1000L
if(NOT keyword_set(seed1)) then seed1=1000L

hdr=headfits(mmatrix)
nspec=long(sxpar(hdr, 'NSPEC'))
nzf=long(sxpar(hdr, 'NZ'))
nfilter=long(sxpar(hdr, 'NFILTER'))
ndusts=long(sxpar(hdr, 'NDUST'))
nmets=long(sxpar(hdr, 'NMET'))
nages=long(sxpar(hdr, 'NAGE'))
lambda=mrdfits(mmatrix, 1)
avloglam=alog10(lambda[0:nspec-1])
filterlist=strtrim(mrdfits(mmatrix, 5), 2)
zf=mrdfits(mmatrix, 6)

lowz=mrdfits(getenv('VAGC_REDUX')+'/lowz/lowz_catalog.'+sample+'.fits',1)
zhelio=lg_to_helio(lowz.zlg, lowz.ra, lowz.dec)

seed=seed1
indx=shuffle_indx(n_elements(lowz), num_sub=nsub, seed=seed)
lowz=lowz[indx]

;;sdss_spec_block, lowz.plate, lowz.fiberid, lowz.mjd, $
;;  block_flux=block_flux, block_ivar=block_ivar, block_lambda=block_lambda, $
;;  avloglam=avloglam
;; convert!!

iz=long(floor((nzf-1.)*(zhelio-zf[0])/(zf[nzf-1]-zf[0])+0.5))

data=fltarr(n_elements(lambda), nsub)
ivar=fltarr(n_elements(lambda), nsub)
;;data[0:nspec-1, *]= block_flux
;;ivar[0:nspec-1, *]= block_ivar
;; ?? skip galex for the moment
for ifilter=2L, nfilter-1L do begin & $
    igood=where(lowz.absmag_ivar[ifilter-2] gt 0, ngood) & $
    if(ngood gt 0) then begin & $
        data[iz[igood]+ifilter*nzf+nspec, igood]= $
          10.^(-0.4*lowz[igood].absmag[ifilter-2]) & $
        ivar[iz[igood]+ifilter*nzf+nspec, igood]= $
  data[iz[igood]+ifilter*nzf+nspec, igood]*0.4*alog(10.)* $
  lowz[igood].absmag_ivar[ifilter-2]) & $
    endif & $
endfor

stop

end
