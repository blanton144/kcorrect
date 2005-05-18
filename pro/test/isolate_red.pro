;+
; NAME:
;   isolate_red
; PURPOSE:
;   isolate concentrated, high luminosity galaxies and try to fit colors
; CALLING SEQUENCE:
;   isolate_red
; REVISION HISTORY:
;   04-May-2005  Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro isolate_red

savfile='data_isolate_red.sav'
if(NOT file_test(savfile)) then begin
    im=hogg_mrdfits(vagc_name('object_sdss_imaging'),1,nrow=28800)
    twomass=hogg_mrdfits(getenv('VAGC_REDUX')+'/object_twomass.fits',1, $
                         nrow=28800)
    sp=sdss_spectro_matched(columns='z')
    kc=hogg_mrdfits(vagc_name('kcorrect', flux='model', band_shift='0.10', $
                              coll='none'),1,nrow=28800)
    sersic=hogg_mrdfits(getenv('VAGC_REDUX')+'/sersic/sersic_catalog.fits', $
                        1, nrow=28800)
    
    umg=kc.absmag[0]-kc.absmag[1]
    gmr=kc.absmag[1]-kc.absmag[2]
    
    ii=where(sp.z gt 0.13 and sp.z lt 0.17 and $
             sersic.sersic_n[2] gt 3. and sersic.sersic_n[2] lt 5. and $
             twomass.twomass_tag gt 0 and $
             gmr gt 0.9 and gmr lt 1.1 and $
             umg gt 1.0 and umg lt 3.1 and $
             randomu(seed, n_elements(gmr)) lt 0.25, nii)

    im=im[ii]
    twomass=twomass[ii]
    sp=sp[ii]
    sersic=sersic[ii]
    kc=kc[ii]
    umg=umg[ii]
    gmr=gmr[ii]
    save, filename=savfile
endif else begin
    restore, savfile
endelse 

sdss_to_maggies, sdss_maggies,sdss_ivar, calibobj=im, flux='petro'
twomass_to_maggies, twomass, twomass_maggies,twomass_ivar

maggies=fltarr(8,n_elements(im))
maggies[0:4,*]=sdss_maggies
maggies[5:7,*]=twomass_maggies
maggies_ivar=fltarr(8,n_elements(im))
maggies_ivar[0:4,*]=sdss_ivar
maggies_ivar[5:7,*]=twomass_ivar

if(NOT file_test('k_nmf_mmatrix.fits')) then $
  k_nmf_mmatrix, filterlist=['sdss_u0.par', $
                             'sdss_g0.par', $
                             'sdss_r0.par', $
                             'sdss_i0.par', $
                             'sdss_z0.par', $
                             'twomass_J.par', $
                             'twomass_H.par', $
                             'twomass_Ks.par'], $
  minzf=0.12, maxzf=0.18, nzf=100, nagesmax=40, /noel, $
  dusts={geometry:'dusty', dust:'MW', structure:'c', tauv:0.}

;; figure out what form we need the data in
if(NOT keyword_set(mmatrix)) then mmatrix='k_nmf_mmatrix.fits'
if(NOT keyword_set(outfile)) then outfile='k_nmf_spdata.fits'
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
ntotal=n_elements(im)
datastr={nx:n_elements(lambda), $
         ny:ntotal, $
         rowstart:lonarr(ntotal), $
         nxrow:lonarr(ntotal)}
zdist=sp.z
zhelio=sp.z
iz=long(floor((nzf-1.)*(zhelio-zf[0])/(zf[nzf-1]-zf[0])+0.5))
xx=0
data=0
ivar=0
for i=0L, ntotal-1L do begin & $
    new_xx=iz[i]+[0,1,2,3,4,5,6,7]*nzf+nspec & $
    new_data= maggies[*,i] & $
    new_ivar= maggies_ivar[*,i] & $
    datastr.rowstart[i]=i*8L & $
    datastr.nxrow[i]=8 & $
    if(NOT keyword_set(xx)) then $
      xx=new_xx $
    else $
      xx=[xx, new_xx] & $
    if(NOT keyword_set(data)) then $
      data=new_data $
    else $
      data=[data, new_data] & $
    if(NOT keyword_set(ivar)) then $
      ivar=new_ivar $
    else $
      ivar=[ivar, new_ivar] & $
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

k_run_nmf, niter=1000, nt=1, /qa, /reset

end
