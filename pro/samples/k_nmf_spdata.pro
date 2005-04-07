;+
; NAME:
;   k_nmf_spdata
; PURPOSE:
;   put together sparse data for nmf fitting
; CALLING SEQUENCE:
;   k_nmf_spdata
; COMMENTS:
;   Assumes k_nmf_mmatrix has been run.
;
;   Takes data from sources:
;     1. SDSS photometric survey
;             ugrizJHK magnitudes
;     2. SDSS spectroscopic survey
;             optical spectrum smoothed to 300 km/s resolution
;     3. GALEX-SDSS matrix
;             FNugriz magnitudes
; 
;   Creates a file k_nmf_spdata.fits of the following form:
;     HDU0: NULL
;     HDU1: single-element structure of the form:
;              .NX - number of possible observations for each galaxy
;              .NY - number of galaxies
;              .ROWSTART[NY] - starting position in data block for
;                              each galaxy
;              .NXROW[NY] - number of observations for each galaxy
;     HDU2: [NVAL] data block of floating-point values in absolute
;           maggies. If you constructed this sparse-sampling into a 
;           [NX, NY] matrix it would look like (ROWSTART and NXROW
;           tell you how to do so):
;            NY=Ngals
;            NX=Ndata= nlambda+nfilter*nz
;            [0:nlambda-1L, *] is the spectral information
;            [nlambda:Ndata-1, *] is through the filters at each redshift
;     HDU3: [NVAL] data block of floating-point inverse variances
;     HDU4: [NVAL] block of column positions for each row (ie. which
;           observations exist for each galaxy (ROWSTART and NXROW
;           tell you to index this list)
;     HDU5: [Ndata]; effective central wavelength for each data point
;     HDU6: [Ngals]; distance (in redshift units) to each object
;     HDU7: [Ngals]; heliocentric redshift of each object
;   This file is readable by k_nmf_run.pro
;   
;   assumes kcorrect files used:
;     correction = [-0.042,  0.036,  0.015, 0.013, -0.002]
;   and we want:
;     correction = [-0.036,  0.012,  0.010, 0.028,  0.040] 
;   from daniel's WDs
;   so we apply:
;                  [ 0.006, -0.024, -0.005, 0.015,  0.042]
; REVISION HISTORY:
;   23-Nov-2004  Michael Blanton (NYU)
;-
;------------------------------------------------------------------------------
pro k_nmf_spdata, mmatrix=mmatrix, sample=sample, flux=flux

if(NOT keyword_set(mmatrix)) then mmatrix='k_nmf_mmatrix.fits'
if(NOT keyword_set(outfile)) then outfile='k_nmf_spdata.fits'
if(NOT keyword_set(sample)) then sample='sample15'
if(NOT keyword_set(flux)) then flux='model'
if(NOT keyword_set(nsdss_photo)) then nsdss_photo=150L
if(NOT keyword_set(nsdss_spec)) then nsdss_spec=150L
if(NOT keyword_set(ngalex)) then ngalex=150L
if(NOT keyword_set(ndeep)) then ndeep=150L
if(NOT keyword_set(seed1)) then seed1=1000L
if(NOT keyword_set(omega0)) then omega0=0.3
if(NOT keyword_set(omegal0)) then omegal0=0.7
if(NOT keyword_set(velmodtype)) then velmodtype='sigv150'
; min errors in FNugrizJHKBRI
minerrors=[0.05, 0.05, 0.05, 0.02, 0.02, 0.02, 0.03, $
           0.05, 0.05, 0.05, 0.02, 0.02, 0.02]
kc2ab=[ 0.006, -0.024, -0.005, 0.015,  0.042, 0., 0., 0., 0., 0., 0.]
seed=seed1

;; relative weights
galex_weight=5.
sdss_spec_weight=0.3
sdss_photo_weight=1.
deep_weight=1.

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

;; create full matrix
ntotal=nsdss_photo+nsdss_spec+ngalex+ndeep
datastr={nx:n_elements(lambda), $
         ny:ntotal, $
         rowstart:lonarr(ntotal), $
         nxrow:lonarr(ntotal)}
currx=0L
data=0
ivar=0
xx=0
zdist=fltarr(ntotal)
zhelio=fltarr(ntotal)
isdss_spec=lindgen(nsdss_spec)
isdss_photo=nsdss_spec+lindgen(nsdss_photo)
igalex=nsdss_spec+nsdss_photo+lindgen(ngalex)
ideep=nsdss_spec+nsdss_photo+ngalex+lindgen(ndeep)

;; sdss spectra
postcat=hogg_mrdfits(vagc_name('post_catalog', sample=sample, letter='bsafe', $
                               post='1'), 1, nrow=28800)
indx_spec=shuffle_indx(n_elements(postcat), num_sub=nsdss_spec, seed=seed)
postcat=postcat[indx_spec]
kc=mrdfits(vagc_name('kcorrect', flux=flux, collision_type='none', $
                     band_shift='0.10'),1,row=postcat.object_position)
sp=mrdfits(vagc_name('object_sdss_spectro'),1,row=postcat.object_position)
vmod=mrdfits(getenv('VAGC_REDUX')+'/velmod_distance/distance_'+velmodtype+ $
             '.fits',1, row=postcat.object_position)
sdss_spec_block, sp.plate, sp.fiberid, sp.mjd, $
  block_flux=block_flux, block_ivar=block_ivar, block_lambda=block_lambda, $
  avloglam=avloglam, /deextinct, vdisp=vdisp
absrc=3.631*2.99792*10.^(15-2.*avloglam)
dm=lf_distmod(vmod.zdist)
for i=0L, n_elements(postcat)-1L do $
  block_flux[*,i]=block_flux[*,i]/absrc*10.^(0.4*dm[i])
for i=0L, n_elements(postcat)-1L do $
  block_ivar[*,i]=block_ivar[*,i]*absrc^2*10.^(-0.8*dm[i])
for i=0L, n_elements(postcat)-1L do begin
    igood=where(block_ivar[*,i] gt 0., ngood)
    if(ngood gt 0) then begin
        if(keyword_set(data)) then begin
            data=[data, block_flux[igood,i]]
            ivar=[ivar, block_ivar[igood,i]*sdss_spec_weight]
            xx=[xx, igood]
        endif else begin
            data=block_flux[igood,i]
            ivar=block_ivar[igood,i]*sdss_spec_weight
            xx=igood
        endelse
        datastr.rowstart[isdss_spec[i]]=currx
        datastr.nxrow[isdss_spec[i]]=ngood
        currx=currx+ngood
    endif
endfor
zdist[isdss_spec]=vmod.zdist
zhelio[isdss_spec]=kc.z

;; sdss photometry
postcat=hogg_mrdfits(vagc_name('post_catalog', sample=sample, letter='bsafe', $
                               post='1'), 1, nrow=28800)
indx_photo=shuffle_indx(n_elements(postcat), num_sub=nsdss_photo, seed=seed)
postcat=postcat[indx_photo]
kc=mrdfits(vagc_name('kcorrect', flux=flux, collision_type='none', $
                     band_shift='0.10'),1,row=postcat.object_position)
sp=mrdfits(vagc_name('object_sdss_spectro'),1,row=postcat.object_position)
vmod=mrdfits(getenv('VAGC_REDUX')+'/velmod_distance/distance_'+velmodtype+ $
             '.fits',1, row=postcat.object_position)
zdist[isdss_photo]=vmod.zdist
dm=lf_distmod(zdist[isdss_photo])
zhelio[isdss_photo]=sp.z
iz=long(floor((nzf-1.)*(zhelio[isdss_photo]-zf[0])/(zf[nzf-1]-zf[0])+0.5))
for i=0L, n_elements(postcat)-1L do begin
    igood=where(kc[i].abmaggies gt 0. and kc[i].abmaggies_ivar gt 0, ngood) 
    if(ngood gt 0) then begin
        new_data=1.e-9*kc[i].abmaggies[igood]*10.^(0.4*dm[i])* $
          10.^(-0.4*kc2ab[igood])
        new_ivar=1.e+18/ $
          (1./kc[i].abmaggies_ivar[igood]+ $
           (0.4*alog(10.)*kc[i].abmaggies[igood]*minerrors[igood+2])^2) $
          *10.^(-0.8*dm[i])*10.^(0.8*kc2ab[igood])*sdss_photo_weight
        new_xx=iz[i]+(igood+2L)*nzf+nspec
        if(keyword_set(data)) then begin
            data=[data, new_data]
            ivar=[ivar, new_ivar]
            xx=[xx, new_xx]
        endif else begin
            data=new_data
            ivar=new_ivar
            xx=new_xx
        endelse
        datastr.rowstart[isdss_photo[i]]=currx
    endif
    datastr.nxrow[isdss_photo[i]]=ngood
    currx=currx+ngood
endfor

;; collect some galex photometry with SDSS too
galex_objects=mrdfits(getenv('VAGC_REDUX')+'/galex/galex_objects.fits',1)
galex=mrdfits(getenv('VAGC_REDUX')+'/galex/galex_catalog.fits',1)
imatched=where(galex_objects.object_tag ge 0)
galex_objects=galex_objects[imatched]
galex=galex[imatched]
galex_lss=mrdfits(vagc_name('lss_index', sample=sample), 1, $
                  row=galex_objects.object_position)
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
                   velmodtype+'.fits',1, row=galex_objects.object_position)
galex_dm=lf_distmod(galex_vmod.zdist)
galex_sdss=mrdfits(vagc_name('object_sdss_imaging'), 1, $
                   row=galex_objects.object_position)
galex_kc=mrdfits(vagc_name('kcorrect', flux=flux, collision_type='none', $
                           band_shift='0.10'),1, $
                 row=galex_objects.object_position)
zdist[igalex]=galex_vmod.zdist
zhelio[igalex]=galex_vmod.zact
iz=long(floor((nzf-1.)*(zhelio[igalex]-zf[0])/(zf[nzf-1]-zf[0])+0.5))
for i=0L, n_elements(galex)-1L do begin
    datastr.rowstart[igalex[i]]=currx

    if(galex[i].fuv_mag ne -99.) then begin
        datastr.nxrow[igalex[i]]=datastr.nxrow[igalex[i]]+1L
        maggies=10.^(-0.4*(galex[i].fuv_mag-galex[i].fuv_extinction- $
                               galex_dm[i]))
        maggies_ivar= $
          1./((galex[i].fuv_magerr^2+minerrors[0]^2)* $
              (0.4*alog(10.)*maggies)^2)
        new_xx=iz[i]+(0L)*nzf+nspec
        if(keyword_set(data)) then begin
            data=[data, maggies]
            ivar=[ivar, maggies_ivar*galex_weight]
            xx=[xx, new_xx]
        endif else begin
            data=maggies
            ivar=maggies_ivar*galex_weight
            xx=new_xx
        endelse
        currx=currx+1L
    endif 

    if(galex[i].nuv_mag ne -99.) then begin
        datastr.nxrow[igalex[i]]=datastr.nxrow[igalex[i]]+1L
        maggies=10.^(-0.4*(galex[i].nuv_mag-galex[i].nuv_extinction- $
                               galex_dm[i]))
        maggies_ivar= $
          1./((galex[i].nuv_magerr^2+minerrors[1]^2)* $
              (0.4*alog(10.)*maggies)^2)
        new_xx=iz[i]+(1L)*nzf+nspec
        if(keyword_set(data)) then begin
            data=[data, maggies]
            ivar=[ivar, maggies_ivar*galex_weight]
            xx=[xx, new_xx]
        endif else begin
            data=maggies
            ivar=maggies_ivar*galex_weight
            xx=new_xx
        endelse
        currx=currx+1L
    endif 

    igood=where(galex_kc[i].abmaggies gt 0., ngood) 
    new_data=1.e-9*galex_kc[i].abmaggies[igood]*10.^(0.4*galex_dm[i])* $
      10.^(-0.4*kc2ab[igood])
    new_ivar=1.e+18/ $
      (1./galex_kc[i].abmaggies_ivar[igood]+ $
       (0.4*alog(10.)*galex_kc[i].abmaggies[igood]*minerrors[igood+2])^2)* $
      galex_weight*10.^(-0.8*galex_dm[i])*10.^(0.8*kc2ab[igood])
    new_xx=iz[i]+(igood+2L)*nzf+nspec
    if(keyword_set(data)) then begin
        data=[data, new_data]
        ivar=[ivar, new_ivar]
        xx=[xx, new_xx]
    endif else begin
        data=new_data
        ivar=new_ivar
        xx=new_xx
    endelse
    datastr.nxrow[igalex[i]]=datastr.nxrow[igalex[i]]+ngood
    currx=currx+ngood
endfor

;; collect some deep photometry 
deep=mrdfits(getenv('KCORRECT_DIR')+ $
             '/data/redshifts/deep/zcat.dr1.uniq.fits.gz',1)
igot=where(deep.zquality ge 3 and deep.zhelio gt 0.01 and deep.zhelio lt 2.) 
deep=deep[igot]
indx_deep=shuffle_indx(n_elements(deep), num_sub=ndeep, seed=seed)
deep=deep[indx_deep]
deep_dm=lf_distmod(deep.zhelio)
zdist[ideep]=deep.zhelio
zhelio[ideep]=deep.zhelio
iz=long(floor((nzf-1.)*(zhelio[ideep]-zf[0])/(zf[nzf-1]-zf[0])+0.5))
for i=0L, n_elements(deep)-1L do begin
    datastr.rowstart[ideep[i]]=currx

    ngood=0L

    iband=10L
    if(deep[i].magb lt 25.) then begin
        datastr.nxrow[ideep[i]]=datastr.nxrow[ideep[i]]+1L
        extinction=deep[i].sfd_ebv*4.32
        maggies=10.^(-0.4*(deep[i].magb-extinction-deep_dm[i]))
        maggies_ivar=1./((0.03^2)*(0.4*alog(10.)*maggies)^2)
        new_xx=iz[i]+(iband)*nzf+nspec
        if(keyword_set(data)) then begin
            data=[data, maggies]
            ivar=[ivar, maggies_ivar*deep_weight]
            xx=[xx, new_xx]
        endif else begin
            data=maggies
            ivar=maggies_ivar*deep_weight
            xx=new_xx
        endelse
        currx=currx+1L
        ngood=ngood+1L
    endif 

    iband=11L
    if(deep[i].magr lt 25.) then begin
        datastr.nxrow[ideep[i]]=datastr.nxrow[ideep[i]]+1L
        extinction=deep[i].sfd_ebv*2.63
        maggies=10.^(-0.4*(deep[i].magr-extinction-deep_dm[i]))
        maggies_ivar=1./((0.03^2)*(0.4*alog(10.)*maggies)^2)
        new_xx=iz[i]+(iband)*nzf+nspec
        if(keyword_set(data)) then begin
            data=[data, maggies]
            ivar=[ivar, maggies_ivar*deep_weight]
            xx=[xx, new_xx]
        endif else begin
            data=maggies
            ivar=maggies_ivar*deep_weight
            xx=new_xx
        endelse
        currx=currx+1L
        ngood=ngood+1L
    endif 

    iband=12L
    if(deep[i].magi lt 25.) then begin
        datastr.nxrow[ideep[i]]=datastr.nxrow[ideep[i]]+1L
        extinction=deep[i].sfd_ebv*1.96
        maggies=10.^(-0.4*(deep[i].magi-extinction-deep_dm[i]))
        maggies_ivar=1./((0.03^2)*(0.4*alog(10.)*maggies)^2)
        new_xx=iz[i]+(iband)*nzf+nspec
        if(keyword_set(data)) then begin
            data=[data, maggies]
            ivar=[ivar, maggies_ivar*deep_weight]
            xx=[xx, new_xx]
        endif else begin
            data=maggies
            ivar=maggies_ivar*deep_weight
            xx=new_xx
        endelse
        currx=currx+1L
        ngood=ngood+1L
    endif 
endfor

hdr=['']
sxaddpar, hdr, 'NSPEC', nspec, 'number of points in spectrum'
sxaddpar, hdr, 'NZ', nzf, 'number of redshifts'
sxaddpar, hdr, 'NFILTER', n_elements(filterlist), 'number of filters'
sxaddpar, hdr, 'NDUST', ndusts, 'number of dusts'
sxaddpar, hdr, 'NMET', nmets, 'number of metallicities'
sxaddpar, hdr, 'NAGE', nages, 'number of ages'
sxaddpar, hdr, 'VDISP', vdisp, 'smoothed to this velocity dispersion (km/s)'
mwrfits, 0, outfile, hdr, /create
mwrfits, datastr, outfile
mwrfits, data, outfile
mwrfits, ivar, outfile
mwrfits, xx, outfile
mwrfits, lambda, outfile
mwrfits, zdist, outfile
mwrfits, zhelio, outfile

end
