;+
; NAME:
;   k_nmf_spdata_spec
; PURPOSE:
;   put together spectroscopic data for nmf fitting
; CALLING SEQUENCE:
;   k_nmf_spdata_spec
; COMMENTS:
;   Assumes k_nmf_mmatrix has been run.
;
;   Takes data from sources:
;     1. SDSS spectroscopic survey
;             optical spectrum smoothed to 300 km/s resolution
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
pro k_nmf_spdata_spec, mmatrix=mmatrix, sample=sample, flux=flux, $
                  galexblue=galexblue

if(NOT keyword_set(mmatrix)) then mmatrix='k_nmf_mmatrix.fits'
if(NOT keyword_set(outfile)) then outfile='k_nmf_spdata.fits'
if(NOT keyword_set(sample)) then sample='sample15'
if(NOT keyword_set(flux)) then flux='model'
if(NOT keyword_set(nlrg_photo)) then nlrg_photo=50L
if(NOT keyword_set(nlrg_spec)) then nlrg_spec=10L
if(NOT keyword_set(nsdss_photo)) then nsdss_photo=100L
if(NOT keyword_set(nsdss_spec)) then nsdss_spec=1000L
if(NOT keyword_set(ngalex)) then ngalex=100L
if(NOT keyword_set(ndeep)) then ndeep=100L
if(NOT keyword_set(ngoods)) then ngoods=100L
if(NOT keyword_set(seed1)) then seed1=1000L
if(NOT keyword_set(omega0)) then omega0=0.3
if(NOT keyword_set(omegal0)) then omegal0=0.7
if(NOT keyword_set(velmodtype)) then velmodtype='sigv150'
; min errors in FNugrizJHKBRIBV
minerrors=[0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, $
           0.05, 0.05, 0.05, 0.02, 0.02, 0.02, 0.02, 0.02]
kc2ab=[ 0.006, -0.024, -0.005, 0.015,  0.042, 0., 0., 0., 0., 0., 0., 0., 0.]
seed=seed1

;; relative weights
galex_weight=1.0
sdss_spec_weight=1.
sdss_photo_weight=1.0
lrg_spec_weight=1.
lrg_photo_weight=1.0
deep_weight=1.0
goods_weight=1.0

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

;; set up indices
ntotal=0L
isdss_spec   =ntotal+lindgen(nsdss_spec)
ntotal=ntotal+nsdss_spec
isdss_photo  =ntotal+lindgen(nsdss_photo)
ntotal=ntotal+nsdss_photo
ilrg_photo   =ntotal+lindgen(nlrg_photo)
ntotal=ntotal+nlrg_photo
ilrg_spec    =ntotal+lindgen(nlrg_spec)
ntotal=ntotal+nlrg_spec
igoods       =ntotal+lindgen(ngoods)
ntotal=ntotal+ngoods
igalex       =ntotal+lindgen(ngalex)
ntotal=ntotal+ngalex
ideep        =ntotal+lindgen(ndeep)
ntotal=ntotal+ndeep

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

;; sdss spectra
postcat=hogg_mrdfits(vagc_name('post_catalog', sample=sample, letter='bsafe', $
                               post='1'), 1, nrow=28800)
gmr=postcat.absm[1]-postcat.absm[2]
isort=sort(gmr)
nbins=11
bins=findgen(nbins+1L)/10.
for i=0L, nbins-1L do begin
    ii=where(gmr gt bins[i] and $
             gmr lt bins[i+1])
    tmp_indx_spec= $
      ii[shuffle_indx(n_elements(ii), num_sub=nsdss_spec/nbins, seed=seed)]
    if(i eq 0) then $
      indx_spec=tmp_indx_spec $
    else $
      indx_spec=[indx_spec, tmp_indx_spec]
endfor
nleft=nsdss_spec-n_elements(indx_spec)
if(nleft gt 0) then $
  indx_spec=[indx_spec, shuffle_indx(n_elements(postcat), num_sub=nleft, $
                                     seed=seed)]
postcat=postcat[indx_spec]
kc=mrdfits(vagc_name('kcorrect', flux=flux, collision_type='none', $
                     band_shift='0.10'),1,row=postcat.object_position)
sp=sdss_spectro_matched()
sp=sp[postcat.object_position]
vmod=mrdfits(getenv('VAGC_REDUX')+'/velmod_distance/distance_'+velmodtype+ $
             '.fits',1, row=postcat.object_position)
sdss_spec_block, sp.plate, sp.fiberid, sp.mjd, $
  block_flux=block_flux, block_ivar=block_ivar, block_lambda=block_lambda, $
  avloglam=avloglam, /deextinct, vdisp=vdisp, minerror=0.03
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
sp=sdss_spectro_matched(columns=['z'])
sp=sp[postcat.object_position]
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

;; lrg photometry
sp=sdss_spectro_matched(nrow=28800,columns=['z', 'primtarget'] )
im=hogg_mrdfits(vagc_name('object_sdss_imaging'),1,nrow=28800)
ilrg=where((sp.primtarget AND 32) gt 0 AND $
           sp.z gt 0.2 and sp.z lt 0.6)
sp=sp[ilrg]
im=im[ilrg]
indx_lrg_photo=shuffle_indx(n_elements(sp), num_sub=nlrg_photo, seed=seed)
sp=sp[indx_lrg_photo]
im=im[indx_lrg_photo]
sdss_to_maggies, maggies, maggies_ivar, calibobj=im
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
        new_ivar=1./ $
          (1./maggies_ivar[igood,i]+ $
           (0.4*alog(10.)*maggies[igood,i]*minerrors[igood+2])^2) $
          *10.^(-0.8*dm[i])*lrg_photo_weight
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

;; lrg spectroscopy
sp=sdss_spectro_matched(nrow=28800,columns=['z', 'primtarget', 'plate', $
	'fiberid', 'mjd'] )
im=hogg_mrdfits(vagc_name('object_sdss_imaging'),1,nrow=28800)
ilrg=where((sp.primtarget AND 32) gt 0 AND $
           sp.z gt 0.2 and sp.z lt 0.6)
sp=sp[ilrg]
im=im[ilrg]
indx_lrg_spec=shuffle_indx(n_elements(sp), num_sub=nlrg_spec, seed=seed)
sp=sp[indx_lrg_spec]
im=im[indx_lrg_spec]
sdss_to_maggies, maggies, maggies_ivar, calibobj=im
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
    igood=where(block_ivar gt 0., ngood)
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

;; collect some goods photometry 
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
glactc, gphoto.ra, gphoto.dec, 2000., gl, gb, 1, /deg
ebv=dust_getval(gl, gb, /noloop)
dfactors=[4.32, 3.32, 2.00, 1.54, 0.90, 0.58, 0.37]
for i=0L, n_elements(gz)-1L do begin
    datastr.rowstart[igoods[i]]=currx

    iband=16L
    if(gphoto[i].bmag_magauto ne -99.) then begin
        datastr.nxrow[igoods[i]]=datastr.nxrow[igoods[i]]+1L
        extinction=ebv[i]*dfactors[0]
        err2=0.02^2+gphoto[i].bmagerr_magauto^2
        maggies=10.^(-0.4*(gphoto[i].bmag_magauto-extinction-goods_dm[i]))
        maggies_ivar=1./(err2*(0.4*alog(10.)*maggies)^2)
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

    iband=17L
    if(gphoto[i].vmag_magauto ne -99.) then begin
        datastr.nxrow[igoods[i]]=datastr.nxrow[igoods[i]]+1L
        extinction=ebv[i]*dfactors[1]
        err2=0.02^2+gphoto[i].vmagerr_magauto^2
        maggies=10.^(-0.4*(gphoto[i].vmag_magauto-extinction-goods_dm[i]))
        maggies_ivar=1./(err2*(0.4*alog(10.)*maggies)^2)
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

    iband=18L
    if(gphoto[i].imag_magauto ne -99.) then begin
        datastr.nxrow[igoods[i]]=datastr.nxrow[igoods[i]]+1L
        extinction=ebv[i]*dfactors[2]
        err2=0.02^2+gphoto[i].imagerr_magauto^2
        maggies=10.^(-0.4*(gphoto[i].imag_magauto-extinction-goods_dm[i]))
        maggies_ivar=1./(err2*(0.4*alog(10.)*maggies)^2)
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

    iband=19L
    if(gphoto[i].zmag_magauto ne -99.) then begin
        datastr.nxrow[igoods[i]]=datastr.nxrow[igoods[i]]+1L
        extinction=ebv[i]*dfactors[3]
        err2=0.02^2+gphoto[i].zmagerr_magauto^2
        maggies=10.^(-0.4*(gphoto[i].zmag_magauto-extinction-goods_dm[i]))
        maggies_ivar=1./(err2*(0.4*alog(10.)*maggies)^2)
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

    iband=13L
    if(gphoto[i].jmag_magauto ne -99.) then begin
        datastr.nxrow[igoods[i]]=datastr.nxrow[igoods[i]]+1L
        extinction=ebv[i]*dfactors[4]
        err2=0.02^2+gphoto[i].jmagerr_magauto^2
        maggies=10.^(-0.4*(gphoto[i].jmag_magauto-extinction-goods_dm[i]))
        maggies_ivar=1./(err2*(0.4*alog(10.)*maggies)^2)
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

    if(0) then begin
        iband=14L
        if(gphoto[i].hmag_magauto ne -99.) then begin
            datastr.nxrow[igoods[i]]=datastr.nxrow[igoods[i]]+1L
            extinction=ebv[i]*dfactors[5]
            err2=0.02^2+gphoto[i].hmagerr_magauto^2
            maggies=10.^(-0.4*(gphoto[i].hmag_magauto-extinction-goods_dm[i]))
            maggies_ivar=1./(err2*(0.4*alog(10.)*maggies)^2)
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
    endif

    iband=15L
    if(gphoto[i].kmag_magauto ne -99.) then begin
        datastr.nxrow[igoods[i]]=datastr.nxrow[igoods[i]]+1L
        extinction=ebv[i]*dfactors[6]
        err2=0.02^2+gphoto[i].kmagerr_magauto^2
        maggies=10.^(-0.4*(gphoto[i].kmag_magauto-extinction-goods_dm[i]))
        maggies_ivar=1./(err2*(0.4*alog(10.)*maggies)^2)
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

    if(galex[i].fuv_mag ne -999. AND galex[i].fuv_mag ne -99.) then begin
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

    if(galex[i].nuv_mag ne -999. AND galex[i].nuv_mag ne -99.) then begin
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

    if(keyword_set(galexblue)) then $
      igood=where(galex_kc[i].abmaggies[0:1] gt 0., ngood) $
    else $
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
