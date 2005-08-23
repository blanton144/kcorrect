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
;     4. DEEP data
;     5. GOODS redshift sample data
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
; REVISION HISTORY:
;   23-Nov-2004  Michael Blanton (NYU)
;-
;------------------------------------------------------------------------------
pro k_nmf_spdata_goods, mmatrix=mmatrix, sample=sample, flux=flux, $
                        galexblue=galexblue, nsdss_spec=nsdss_spec, $
                        nlrg_spec=nlrg_spec, spchop=spchop, few=few, $
                        seed=seed1

spchop=1
if(n_elements(mmatrix) eq 0) then mmatrix='k_nmf_mmatrix.fits'
if(n_elements(outfile) eq 0) then outfile='k_nmf_spdata.fits'
if(n_elements(sample) eq 0) then sample='dr4'
if(n_elements(flux) eq 0) then flux='petro'
if(n_elements(nlrg_photo) eq 0) then nlrg_photo=0L
if(n_elements(nlrg_spec) eq 0) then nlrg_spec=0L
if(n_elements(nsdss_photo) eq 0) then nsdss_photo=0L
if(n_elements(nsdss_spec) eq 0) then nsdss_spec=100L
if(n_elements(ngalex) eq 0) then ngalex=0L
if(n_elements(ndeep) eq 0) then ndeep=0L
if(n_elements(ngoods) eq 0) then ngoods=1000L
if(n_elements(nswire) eq 0) then nswire=0L
if(n_elements(seed1) eq 0) then seed1=1002L
if(n_elements(omega0) eq 0) then omega0=0.3
if(n_elements(omegal0) eq 0) then omegal0=0.7
seed=seed1

;; relative weights
galex_weight=1.0
sdss_spec_weight=0.003
sdss_photo_weight=[1.,1.,1.,1.,1.]
lrg_spec_weight=0.003
lrg_photo_weight=1.0
deep_weight=1.0
goods_weight=1.0
swire_weight=0.001
twomass_weight=1.0

spfull=sdss_spectro_matched(columns=['plate', $
                                    'fiberid', $
                                    'mjd', $
                                    'z', $
                                    'primtarget'])
imfull=hogg_mrdfits(vagc_name('object_sdss_imaging'),1,nrow=28800)

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

keep=bytarr(n_elements(avloglam))+1
if(keyword_set(spchop)) then begin
    wave=10.^avloglam
    keep=bytarr(n_elements(avloglam))
    ikeep=where((wave gt 3700. and wave lt 4300.) OR $
                (wave gt 6450. and wave lt 6800.) OR $
                (wave gt 4700. and wave lt 5100.), nkeep)
    if(nkeep gt 0) then $
      keep[ikeep]=1
endif

;; set up indices
ntotal=0L
if(nsdss_photo gt 0) then $
  isdss_photo  =ntotal+lindgen(nsdss_photo)
ntotal=ntotal+nsdss_photo
if(nlrg_photo gt 0) then $
  ilrg_photo   =ntotal+lindgen(nlrg_photo)
ntotal=ntotal+nlrg_photo
if(ngoods gt 0) then $
  igoods       =ntotal+lindgen(ngoods)
ntotal=ntotal+ngoods
if(nswire gt 0) then $
  iswire       =ntotal+lindgen(nswire)
ntotal=ntotal+nswire
if(ngalex gt 0) then $
  igalex       =ntotal+lindgen(ngalex)
ntotal=ntotal+ngalex
if(ndeep gt 0) then $
  ideep        =ntotal+lindgen(ndeep)
ntotal=ntotal+ndeep
if(nsdss_spec gt 0) then $
  isdss_spec   =ntotal+lindgen(nsdss_spec)
ntotal=ntotal+nsdss_spec
if(nlrg_spec gt 0) then $
  ilrg_spec    =ntotal+lindgen(nlrg_spec)
ntotal=ntotal+nlrg_spec

;; create full matrix
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

if(ngalex gt 0) then begin
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
    galex_dm=lf_distmod(galex_lss.z)
    galex_sdss=imfull[galex_objects.object_position]
    sdss_to_maggies, sdss_maggies, sdss_ivar, calibobj=galex_sdss, flux=flux
    zdist[igalex]=galex_lss.z
    zhelio[igalex]=galex_lss.z
    iz=long(floor((nzf-1.)*(zhelio[igalex]-zf[0])/(zf[nzf-1]-zf[0])+0.5))
    for i=0L, n_elements(galex)-1L do begin
        datastr.rowstart[igalex[i]]=currx

        if(galex[i].fuv_mag ne -999. AND galex[i].fuv_mag ne -99.) then begin
            datastr.nxrow[igalex[i]]=datastr.nxrow[igalex[i]]+1L
            maggies=10.^(-0.4*(galex[i].fuv_mag-galex[i].fuv_extinction- $
                               galex_dm[i]))
            maggies_ivar= $
              1./((galex[i].fuv_magerr^2+0.02^2)* $
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

        if(galex[i].nuv_mag ne -999. AND galex[i].nuv_mag ne -99.) then begin
            datastr.nxrow[igalex[i]]=datastr.nxrow[igalex[i]]+1L
            maggies=10.^(-0.4*(galex[i].nuv_mag-galex[i].nuv_extinction- $
                               galex_dm[i]))
            maggies_ivar= $
              1./((galex[i].nuv_magerr^2+0.02^2)* $
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

        if(keyword_set(galexblue)) then $
          igood=where(sdss_maggies[0:1,i] gt 0., ngood) $
        else $
          igood=where(sdss_maggies[*,i] gt 0., ngood) 
        new_data=sdss_maggies[igood,i]*10.^(0.4*galex_dm[i])
        new_ivar=sdss_ivar[igood,i]*sdss_photo_weight[igood]*10.^(-0.8*galex_dm[i])
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
endif

;; collect some goods photometry
if(ngoods gt 0) then begin
    gphoto=rsex(getenv('KCORRECT_DIR')+ $
                '/data/redshifts/goods/mb_cdfs_isaac_ks_photz_c1.1_d3.0k.cat')
    gz=rsex(getenv('KCORRECT_DIR')+ $
            '/data/redshifts/goods/mb_z.cdfs.c1.1z.20050218.cat')
    igot=where(gz.z lt 2. and gz.z gt 0.1)
    gz=gz[igot]
    gphoto=gphoto[igot]
    indx_g=shuffle_indx(n_elements(gz), num_sub=ngoods, seed=seed)
    gz=gz[indx_g]
    gphoto=gphoto[indx_g]
    goods_dm=lf_distmod(gz.z)
    zdist[igoods]=gz.z
    zhelio[igoods]=gz.z
    iz=long(floor((nzf-1.)*(zhelio[igoods]-zf[0])/(zf[nzf-1]-zf[0])+0.5))
    goods_to_maggies, gphoto, goods_maggies, goods_ivar
    
    gbands=[16,17,18,19,13,14,15]
    for i=0L, n_elements(gz)-1L do begin
        datastr.rowstart[igoods[i]]=currx
        for j=0L, 6L do begin
            iband=gbands[j]
            if(goods_ivar[j,i] gt 0.) then begin
                datastr.nxrow[igoods[i]]=datastr.nxrow[igoods[i]]+1L
                maggies=goods_maggies[j,i]*10.^(0.4*goods_dm[i])
                maggies_ivar=goods_ivar[j,i]*10.^(-0.8*goods_dm[i])
                new_xx=iz[i]+(iband)*nzf+nspec
                if(keyword_set(data)) then begin
                    data=[data, maggies]
                    ivar=[ivar, maggies_ivar*goods_weight]
                    xx=[xx, new_xx]
                endif else begin
                    data=maggies
                    ivar=maggies_ivar*goods_weight
                    xx=new_xx
                endelse
                currx=currx+1L
            endif
        endfor 
    endfor
endif

if(nsdss_spec gt 0) then begin
;; sdss spectra
    postcat=hogg_mrdfits(vagc_name('post_catalog', sample=sample, letter='bsafe', $
                                   post='1'), 1, nrow=28800)
    indx_spec=shuffle_indx(n_elements(postcat), num_sub=nsdss_spec, seed=seed)
    postcat=postcat[indx_spec]
    sp=spfull[postcat.object_position]
    sdss_spec_block, sp.plate, sp.fiberid, sp.mjd, $
      block_flux=block_flux, block_ivar=block_ivar, block_lambda=block_lambda, $
      avloglam=avloglam, /deextinct, vdisp=vdisp
    absrc=3.631*2.99792*10.^(15-2.*avloglam)
    dm=lf_distmod(sp.z)
    for i=0L, n_elements(postcat)-1L do $
      block_flux[*,i]=block_flux[*,i]/absrc*10.^(0.4*dm[i])
    for i=0L, n_elements(postcat)-1L do $
      block_ivar[*,i]=block_ivar[*,i]*absrc^2*10.^(-0.8*dm[i])
    for i=0L, n_elements(postcat)-1L do begin
        igood=where(block_ivar[*,i] gt 0. and keep gt 0, ngood)
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
    zdist[isdss_spec]=sp.z
    zhelio[isdss_spec]=sp.z
endif

if(nsdss_photo gt 0) then begin
;; sdss photometry
    postcat=hogg_mrdfits(vagc_name('post_catalog', sample=sample, letter='bsafe', $
                                   post='1'), 1, nrow=28800)
    indx_photo=shuffle_indx(n_elements(postcat), num_sub=nsdss_photo, seed=seed)
    postcat=postcat[indx_photo]
    im=imfull[postcat.object_position]
    twomass=mrdfits(vagc_name('object_twomass'),1,row=postcat.object_position)
    sdss_to_maggies, sdss_maggies, sdss_ivar, calibobj=im, flux=flux
    twomass_to_maggies, twomass, twomass_maggies, twomass_ivar
    maggies=fltarr(8, n_elements(twomass)) 
    maggies_ivar=fltarr(8, n_elements(twomass)) 
    maggies[0:4, *]= sdss_maggies
    maggies_ivar[0:4, *]= sdss_ivar
    maggies[5:7, *]= twomass_maggies
    maggies_ivar[5:7, *]= twomass_ivar
    sp=spfull[postcat.object_position]
    zdist[isdss_photo]=sp.z
    dm=lf_distmod(zdist[isdss_photo])
    zhelio[isdss_photo]=sp.z
    iz=long(floor((nzf-1.)*(zhelio[isdss_photo]-zf[0])/(zf[nzf-1]-zf[0])+0.5))
    weight=[sdss_photo_weight, replicate(twomass_weight,3)]
    for i=0L, n_elements(postcat)-1L do begin
        igood=where(maggies[*,i] gt 0. and maggies_ivar[*,i] gt 0, ngood) 
        if(ngood gt 0) then begin
            new_data=maggies[igood,i]*10.^(0.4*dm[i])
            new_ivar=maggies_ivar[igood,i]*10.^(-0.8*dm[i])*weight[igood]
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
endif

if(nlrg_photo gt 0) then begin
;; lrg photometry
    ilrg=where((spfull.primtarget AND 32) gt 0 AND $
               spfull.z gt 0.2 and spfull.z lt 0.6)
    sp=spfull[ilrg]
    im=imfull[ilrg]
    indx_lrg_photo=shuffle_indx(n_elements(sp), num_sub=nlrg_photo, seed=seed)
    sp=sp[indx_lrg_photo]
    im=im[indx_lrg_photo]
    sdss_to_maggies, maggies, maggies_ivar, calibobj=im, flux=flux
    maggies_ivar[0,*]=0. ;; ignore u-band for LRGs
    zdist[ilrg_photo]=sp.z
    dm=lf_distmod(zdist[ilrg_photo])
    zhelio[ilrg_photo]=sp.z
    iz=long(floor((nzf-1.)*(zhelio[ilrg_photo]-zf[0])/(zf[nzf-1]-zf[0])+0.5))
    for i=0L, n_elements(im)-1L do begin
        datastr.rowstart[ilrg_photo[i]]=currx

        igood=where(maggies[*,i] gt 0. and maggies_ivar[*,i] gt 0, ngood) 
        if(ngood gt 0) then begin
            new_data=maggies[igood,i]*10.^(0.4*dm[i])
            new_ivar=maggies_ivar[igood,i]*10.^(-0.8*dm[i])*lrg_photo_weight
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
        endif
        datastr.nxrow[ilrg_photo[i]]=ngood
        currx=currx+ngood
    endfor
endif

if(nlrg_spec gt 0) then begin
;; lrg spectroscopy
    ilrg=where((spfull.primtarget AND 32) gt 0 AND $
               spfull.z gt 0.2 and spfull.z lt 0.6)
    sp=spfull[ilrg]
    im=imfull[ilrg]
    indx_lrg_spec=shuffle_indx(n_elements(sp), num_sub=nlrg_spec, seed=seed)
    sp=sp[indx_lrg_spec]
    im=im[indx_lrg_spec]
    zdist[ilrg_spec]=sp.z
    zhelio[ilrg_spec]=sp.z
    dm=lf_distmod(zdist[ilrg_spec])
    iz=long(floor((nzf-1.)*(zhelio[ilrg_spec]-zf[0])/(zf[nzf-1]-zf[0])+0.5))
    for i=0L, n_elements(im)-1L do begin
        datastr.rowstart[ilrg_spec[i]]=currx

        sdss_spec_block, sp[i].plate, sp[i].fiberid, sp[i].mjd, $
          block_flux=block_flux, block_ivar=block_ivar, $
          block_lambda=block_lambda, avloglam=avloglam, /deextinct, vdisp=vdisp
        absrc=3.631*2.99792*10.^(15-2.*avloglam)
        block_flux=block_flux/absrc*10.^(0.4*dm[i])
        block_ivar=block_ivar*absrc^2*10.^(-0.8*dm[i])
        igood=where(block_ivar gt 0. and keep gt 0, ngood)
        if(ngood gt 0) then begin
            if(keyword_set(data)) then begin
                data=[data, block_flux[igood]]
                ivar=[ivar, block_ivar[igood]*lrg_spec_weight]
                xx=[xx, igood]
            endif else begin
                data=block_flux[igood]
                ivar=block_ivar[igood]*lrg_spec_weight
                xx=igood
            endelse
            datastr.nxrow[ilrg_spec[i]]=ngood
            currx=currx+ngood
        endif
    endfor
endif

if(nswire gt 0) then begin
;; collect some swire photometry 
    swire=mrdfits(getenv('VAGC_REDUX')+'/spitzer/swire_catalog.fits',1)
    objects=mrdfits(getenv('VAGC_REDUX')+'/spitzer/swire_objects.fits',1)
    iin=where(objects.object_position ge 0)
    sp=spfull[objects[iin].object_position]
    im=imfull[objects[iin].object_position]
    igal=where(sp.z gt 0.003 and sp.z lt 0.5)
    sp=sp[igal]
    im=im[igal]
    swire=swire[iin[igal]]
    indx_s=shuffle_indx(n_elements(swire), num_sub=nswire, seed=seed)
    swire=swire[indx_s]
    sp=sp[indx_s]
    im=im[indx_s]
    swire_to_maggies, swire, swire_maggies, swire_ivar
    sdss_to_maggies, sdss_maggies, sdss_ivar, calibobj=im, flux=flux
    swire_maggies[0:4,*]=sdss_maggies
    swire_ivar[0:4,*]=sdss_ivar
    swire_dm=lf_distmod(sp.z)
    zdist[iswire]=sp.z
    zhelio[iswire]=sp.z
    iz=long(floor((nzf-1.)*(zhelio[iswire]-zf[0])/(zf[nzf-1]-zf[0])+0.5))
    iband=[2,3,4,5,6,20,21,22,23,24]
    for i=0L, n_elements(swire)-1L do begin
        datastr.rowstart[iswire[i]]=currx
        igood=where(swire_ivar[*,i] gt 0., ngood)
        datastr.nxrow[iswire[i]]=ngood
        new_xx=iz[i]+(iband[igood])*nzf+nspec
        if(keyword_set(data)) then begin
            data=[data, swire_maggies[igood,i]*10.^(0.4*swire_dm[i])]
            ivar=[ivar, swire_ivar[igood,i]*swire_weight* $
                  10.^(-0.8*swire_dm[i])]
            xx=[xx, new_xx]
        endif else begin
            data=swire_maggies[igood,i]*10.^(0.4*swire_dm[i])
            ivar=swire_ivar[igood,i]*swire_weight*10.^(-0.8*swire_dm[i])
            xx=new_xx
        endelse
        currx=currx+ngood
    endfor
endif


if(ndeep gt 0) then begin
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
    sigmbase=0.05
    for i=0L, n_elements(deep)-1L do begin
        datastr.rowstart[ideep[i]]=currx

        iband=10L
        mbase=24.0
        if(deep[i].magb lt 25.) then begin
            datastr.nxrow[ideep[i]]=datastr.nxrow[ideep[i]]+1L
            extinction=deep[i].sfd_ebv*4.32
            deep_err2=0.02^2+(sigmbase*10.^(-0.4*(mbase-deep[i].magb)))^2
            maggies=10.^(-0.4*(deep[i].magb-extinction-deep_dm[i]))
            maggies_ivar=1./(deep_err2*(0.4*alog(10.)*maggies)^2)
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
        endif 

        iband=11L
        mbase=24.0
        if(deep[i].magr lt 25.) then begin
            datastr.nxrow[ideep[i]]=datastr.nxrow[ideep[i]]+1L
            extinction=deep[i].sfd_ebv*2.63
            maggies=10.^(-0.4*(deep[i].magr-extinction-deep_dm[i]))
            deep_err2=0.02^2+(sigmbase*10.^(-0.4*(mbase-deep[i].magr)))^2
            maggies_ivar=1./(deep_err2*(0.4*alog(10.)*maggies)^2)
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
        endif 

        iband=12L
        mbase=24.0
        if(deep[i].magi lt 25.) then begin
            datastr.nxrow[ideep[i]]=datastr.nxrow[ideep[i]]+1L
            extinction=deep[i].sfd_ebv*1.96
            maggies=10.^(-0.4*(deep[i].magi-extinction-deep_dm[i]))
            deep_err2=0.02^2+(sigmbase*10.^(-0.4*(mbase-deep[i].magi)))^2
            maggies_ivar=1./(deep_err2*(0.4*alog(10.)*maggies)^2)
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
        endif 
    endfor
endif

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
