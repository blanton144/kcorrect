;+
; NAME:
;   k_nmf_data
; PURPOSE:
;   put together data for nmf fitting
; CALLING SEQUENCE:
;   k_nmf_data
; COMMENTS:
;   Takes data from sources:
;     1. SDSS photometric survey: currently the low-z catalog
;             ugrizJHK magnitudes
;     2. SDSS spectroscopic survey: currently the low-z catalog
;             optical spectrum smoothed to 300 km/s resolution
;     3. GALEX-SDSS matrix
;             FNugriz magnitudes
;   Creates a file k_nmf_data.fits of the following form:
;     HDU0: [Ndata, Ngals]; data in absolute maggies, where
;         Ndata= nlambda+nfilter*nz
;         [0:nlambda-1L, *] is the spectral information
;         [nlambda:Ndata-1, *] is through the filters at each redshift
;     HDU1: [Ndata, Ngals]; inverse variance of data
;     HDU2: [Ndata]; effective central wavelength for each data poin5
;     HDU3: [Ngals]; distance (in redshift units) to each object
;     HDU4: [Ngals]; heliocentric redshift of each object
; BUGS:
;   AB corrections need updated
;   tweak spectra to match fluxes??
;   what about filter curves?
;   smooth spectra to constant vdisp?
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
if(NOT keyword_set(ngalex)) then ngalex=160L
if(NOT keyword_set(seed1)) then seed1=1000L
if(NOT keyword_set(omega0)) then omega0=0.3
if(NOT keyword_set(omegal0)) then omegal0=0.7
if(NOT keyword_set(velmodtype)) then velmodtype='sigv150'

;; relative weights
galex_weight=5.
sdss_spec_weight=0.3
sdss_photo_weight=1.

;; figure out what form we need the data in
hdr=headfits(mmatrix)
nspec=long(sxpar(hdr, 'NSPEC'))
nzf=long(sxpar(hdr, 'NZ'))
nfilter=long(sxpar(hdr, 'NFILTER'))
ndusts=long(sxpar(hdr, 'NDUST'))
nmets=long(sxpar(hdr, 'NMET'))
nages=long(sxpar(hdr, 'NAGE'))
vdisp=float(sxpar(hdr, 'VDISP'))
lambda=mrdfits(mmatrix, 1)
avloglam=alog10(lambda[0:nspec-1])
filterlist=strtrim(mrdfits(mmatrix, 5), 2)
zf=mrdfits(mmatrix, 6)

;; read in the low-z catalog
lowz=mrdfits(getenv('VAGC_REDUX')+'/lowz/lowz_catalog.'+sample+'.fits',1)
seed=seed1

;; create full matrix
ntotal=nsdss_photo+nsdss_spec+ngalex
data=fltarr(n_elements(lambda), ntotal)
ivar=fltarr(n_elements(lambda), ntotal)
zdist=fltarr(ntotal)
zhelio=fltarr(ntotal)
isdss_spec=lindgen(nsdss_spec)
isdss_photo=nsdss_spec+lindgen(nsdss_photo)
igalex=nsdss_spec+nsdss_photo+lindgen(ngalex)

;; collect some spectra
indx_spec=shuffle_indx(n_elements(lowz), num_sub=nsdss_spec, seed=seed)
lowz_spec=lowz[indx_spec]
sdss_spec_block, lowz_spec.plate, lowz_spec.fiberid, lowz_spec.mjd, $
  block_flux=block_flux, block_ivar=block_ivar, block_lambda=block_lambda, $
  avloglam=avloglam, /deextinct, vdisp=vdisp
absrc=3.631*2.99792*10.^(15-2.*avloglam)
dm=lf_distmod(lowz_spec.zdist)
for i=0L, n_elements(lowz_spec)-1L do $
  block_flux[*,i]=block_flux[*,i]/absrc*10.^(0.4*dm[i])
for i=0L, n_elements(lowz_spec)-1L do $
  block_ivar[*,i]=block_ivar[*,i]*absrc^2*10.^(-0.8*dm[i])
data[0:nspec-1, isdss_spec]=block_flux
ivar[0:nspec-1, isdss_spec]=block_ivar
zdist[isdss_spec]=lowz_spec.zdist
zhelio[isdss_spec]=lg_to_helio(lowz_spec.zlg, lowz_spec.ra, lowz_spec.dec)
ivar[*,isdss_spec]=ivar[*,isdss_spec]*sdss_spec_weight

;; collect some sdss photometry
indx_photo=shuffle_indx(n_elements(lowz), num_sub=nsdss_photo, seed=seed)
lowz_photo=lowz[indx_photo]
zdist[isdss_photo]=lowz_photo.zdist
zhelio[isdss_photo]=lg_to_helio(lowz_photo.zlg, lowz_photo.ra, lowz_photo.dec)
iz=long(floor((nzf-1.)*(zhelio[isdss_photo]-zf[0])/(zf[nzf-1]-zf[0])+0.5))
for ifilter=2L, nfilter-1L do begin 
    igood=where(lowz_photo.absmag_ivar[ifilter-2] gt 0, ngood) 
    if(ngood gt 0) then begin 
        data[iz[igood]+ifilter*nzf+nspec, isdss_photo[igood]]= $
          10.^(-0.4*(lowz_photo[igood].absmag[ifilter-2]+ $
                     lowz_photo[igood].kcorrect[ifilter-2])) 
        ivar[iz[igood]+ifilter*nzf+nspec, isdss_photo[igood]]= $
          lowz_photo[igood].absmag_ivar[ifilter-2]/ $
          (data[iz[igood]+ifilter*nzf+nspec,isdss_photo[igood]]* $
           0.4*alog(10.))^2
    endif 
endfor
ivar[*,isdss_photo]=ivar[*,isdss_photo]*sdss_photo_weight

;; collect some galex photometry with SDSS too
galex_objects=mrdfits(getenv('VAGC_REDUX')+'/galex/galex_objects.fits',1)
galex=mrdfits(getenv('VAGC_REDUX')+'/galex/galex_catalog.fits',1)
imatched=where(galex_objects.object_tag ge 0)
galex_objects=galex_objects[imatched]
galex=galex[imatched]
galex_lss=mrdfits(vagc_name('lss_index', sample=sample), 1, $
                  row=galex_objects.object_tag)
igot=where(galex_lss.ztype gt 0 and $
           galex_lss.z gt 0.01 and $
           galex_lss.z lt 0.30)
galex_lss=galex_lss[igot]
galex_objects=galex_objects[igot]
galex=galex[igot]
indx_galex=shuffle_indx(n_elements(galex), num_sub=ngalex, seed=seed)
galex=galex[indx_galex]
galex_objects=galex_objects[indx_galex]
galex_lss=galex_lss[indx_galex]
galex_vmod=mrdfits(getenv('VAGC_REDUX')+'/velmod_distance/distance_'+ $
                   velmodtype+'.fits',1, row=galex_objects.object_tag)
galex_dm=lf_distmod(galex_vmod.zdist)
galex_sdss=mrdfits(vagc_name('object_sdss_imaging'), 1, $
                   row=galex_objects.object_tag)
zdist[igalex]=galex_vmod.zdist
zhelio[igalex]=galex_vmod.zact
iz=long(floor((nzf-1.)*(zhelio[igalex]-zf[0])/(zf[nzf-1]-zf[0])+0.5))
ifilter=0
igood=where(galex.fuv_mag ne -99., ngood)
if(ngood gt 0) then begin
    data[iz[igood]+ifilter*nzf+nspec, igalex[igood]]= $
      10.^(-0.4*(galex[igood].fuv_mag-galex[igood].fuv_extinction- $
                 galex_dm[igood]))
    ivar[iz[igood]+ifilter*nzf+nspec, igalex[igood]]= $
      1./(galex[igood].fuv_magerr^2* $
          (data[iz[igood]+ifilter*nzf+nspec,igalex[igood]]* $
           0.4*alog(10.))^2)
endif
ifilter=1
igood=where(galex.nuv_mag ne -99., ngood)
if(ngood gt 0) then begin
    data[iz[igood]+ifilter*nzf+nspec, igalex[igood]]= $
      10.^(-0.4*(galex[igood].nuv_mag-galex[igood].nuv_extinction- $
                 galex_dm[igood]))
    ivar[iz[igood]+ifilter*nzf+nspec, igalex[igood]]= $
      1./(galex[igood].nuv_magerr^2* $
          (data[iz[igood]+ifilter*nzf+nspec,igalex[igood]]* $
           0.4*alog(10.))^2)
endif
for ifilter=2L, 6L do begin 
    igood=where(galex_sdss.petroflux_ivar[ifilter-2] gt 0, ngood) 
    if(ngood gt 0) then begin 
        data[iz[igood]+ifilter*nzf+nspec, igalex[igood]]= $
          galex_sdss[igood].petroflux[ifilter-2]* $
          10.^(-9.+ $
               0.4*(galex_sdss[igood].extinction[ifilter-2]+galex_dm[igood]))
        ivar[iz[igood]+ifilter*nzf+nspec, igalex[igood]]= $
          galex_sdss[igood].petroflux_ivar[ifilter-2]* $
          10.^(18.- $
               0.8*(galex_sdss[igood].extinction[ifilter-2]+galex_dm[igood]))
    endif 
endfor
ivar[*,igalex]=ivar[*,igalex]*galex_weight

hdr=['']
sxaddpar, hdr, 'NSPEC', nspec, 'number of points in spectrum'
sxaddpar, hdr, 'NZ', nzf, 'number of redshifts'
sxaddpar, hdr, 'NFILTER', n_elements(filterlist), 'number of filters'
sxaddpar, hdr, 'NDUST', ndusts, 'number of dusts'
sxaddpar, hdr, 'NMET', nmets, 'number of metallicities'
sxaddpar, hdr, 'NAGE', nages, 'number of ages'
sxaddpar, hdr, 'VDISP', vdisp, 'smoothed to this velocity dispersion (km/s)'
mwrfits, 0, outfile, hdr, /create
mwrfits, data, outfile
mwrfits, ivar, outfile
mwrfits, lambda, outfile
mwrfits, zdist, outfile
mwrfits, zhelio, outfile

stop

end
