;+
; NAME:
;   k_nmf_data
; PURPOSE:
;   put together data for nmf fitting
; CALLING SEQUENCE:
;   k_nmf_data
; BUGS:
;   AB corrections need updated
;   tweak spectra to match fluxes??
;   what about filter curves?
;   smooth spectra to constant vdisp?
;   galactic reddening correction to spectra
; REVISION HISTORY:
;   23-Nov-2004  Michael Blanton (NYU)
;-
;------------------------------------------------------------------------------
pro k_nmf_data, mmatrix=mmatrix, sample=sample

if(NOT keyword_set(mmatrix)) then mmatrix='k_nmf_mmatrix.fits'
if(NOT keyword_set(outfile)) then outfile='k_nmf_data.fits'
if(NOT keyword_set(sample)) then sample='drtwo14'
if(NOT keyword_set(nsdss_photo)) then nsdss_photo=100L
if(NOT keyword_set(nsdss_spec)) then nsdss_spec=100L
if(NOT keyword_set(seed1)) then seed1=1000L
if(NOT keyword_set(omega0)) then omega0=0.3
if(NOT keyword_set(omegal0)) then omegal0=0.7

;; figure out what form we need the data in
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

;; read in the low-z catalog
lowz=mrdfits(getenv('VAGC_REDUX')+'/lowz/lowz_catalog.'+sample+'.fits',1)
seed=seed1

;; create full matrix
ntotal=nsdss_photo+nsdss_spec
data=fltarr(n_elements(lambda), ntotal)
ivar=fltarr(n_elements(lambda), ntotal)
redshifts=fltarr(ntotal)
isdss_spec=lindgen(nsdss_spec)
isdss_photo=nsdss_spec+lindgen(nsdss_photo)

;; collect some spectra
indx_spec=shuffle_indx(n_elements(lowz), num_sub=nsdss_spec, seed=seed)
lowz_spec=lowz[indx_spec]
sdss_spec_block, lowz_spec.plate, lowz_spec.fiberid, lowz_spec.mjd, $
  block_flux=block_flux, block_ivar=block_ivar, block_lambda=block_lambda, $
  avloglam=avloglam
absrc=3.631*2.99792*10.^(15-2.*avloglam)
dm=lf_distmod(lowz_spec.zdist)
for i=0L, n_elements(lowz_spec)-1L do $
  block_flux[*,i]=block_flux[*,i]/absrc*10.^(0.4*dm[i])
for i=0L, n_elements(lowz_spec)-1L do $
  block_ivar[*,i]=block_ivar[*,i]*absrc^2*10.^(-0.8*dm[i])
data[0:nspec-1, 0:nsdss_spec-1L]=block_flux
ivar[0:nspec-1, 0:nsdss_spec-1L]=block_ivar
redshifts[isdss_spec]=lowz_spec.zdist

;; collect some sdss photometry
indx_photo=shuffle_indx(n_elements(lowz), num_sub=nsdss_photo, seed=seed)
lowz_photo=lowz[indx_photo]
redshifts[isdss_photo]=lowz_photo.zdist
zhelio=lg_to_helio(lowz_photo.zlg, lowz_photo.ra, lowz_photo.dec)
iz=long(floor((nzf-1.)*(zhelio-zf[0])/(zf[nzf-1]-zf[0])+0.5))
for ifilter=2L, nfilter-1L do begin 
    igood=where(lowz_photo.absmag_ivar[ifilter-2] gt 0, ngood) 
    if(ngood gt 0) then begin 
        data[iz[igood]+ifilter*nzf+nspec, isdss_photo[igood]]= $
          10.^(-0.4*(lowz_photo[igood].absmag[ifilter-2]+ $
                     lowz_photo[igood].kcorrect[ifilter-2])) 
        ivar[iz[igood]+ifilter*nzf+nspec, isdss_photo[igood]]= $
          data[iz[igood]+ifilter*nzf+nspec,isdss_photo[igood]]*0.4*alog(10.)* $
          lowz_photo[igood].absmag_ivar[ifilter-2] 
    endif 
endfor

hdr=['']
sxaddpar, hdr, 'NSPEC', nspec, 'number of points in spectrum'
sxaddpar, hdr, 'NZ', nzf, 'number of redshifts'
sxaddpar, hdr, 'NFILTER', n_elements(filterlist), 'number of filters'
sxaddpar, hdr, 'NDUST', ndusts, 'number of dusts'
sxaddpar, hdr, 'NMET', nmets, 'number of metallicities'
sxaddpar, hdr, 'NAGE', nages, 'number of ages'
mwrfits, 0, outfile, hdr, /create
mwrfits, data, outfile
mwrfits, ivar, outfile
mwrfits, lambda, outfile
mwrfits, redshifts, outfile

stop

end
