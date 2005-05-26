;+
; NAME:
;   k_nmf_spdata_swire_blue
; PURPOSE:
;   put together blue SWIRE catalog to fit
; CALLING SEQUENCE:
;   k_nmf_spdata_swire_blue
; COMMENTS:
;   Assumes k_nmf_mmatrix has been run.
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
pro k_nmf_spdata_swire_blue, mmatrix=mmatrix, sample=sample, flux=flux
                             

if(NOT keyword_set(mmatrix)) then mmatrix='k_nmf_mmatrix.fits'
if(NOT keyword_set(outfile)) then outfile='k_nmf_spdata.fits'
if(NOT keyword_set(sample)) then sample='sample15'
if(NOT keyword_set(flux)) then flux='petro'
if(NOT keyword_set(nswire)) then nswire=65L
if(NOT keyword_set(seed1)) then seed1=1000L
seed=seed1

spfull=sdss_spectro_matched(columns=['plate', $
                                    'fiberid', $
                                    'mjd', $
                                    'z', $
                                    'primtarget'])

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
if(nswire gt 0) then $
  iswire       =ntotal+lindgen(nswire)
ntotal=ntotal+nswire

swire_weight=1.

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

if(nswire gt 0) then begin
;; collect some swire photometry 
    swire=mrdfits(getenv('VAGC_REDUX')+'/spitzer/swire_catalog.fits',1)
    objects=mrdfits(getenv('VAGC_REDUX')+'/spitzer/swire_objects.fits',1)
    iin=where(objects.object_position ge 0)
    sp=spfull[objects[iin].object_position]

    ;; get galaxies
    igal=where(sp.z gt 0.003 and sp.z lt 0.5)
    sp=sp[igal]
    swire=swire[iin[igal]]
    objects=objects[iin[igal]]

    ;; get blue
    swire_to_maggies, swire, swire_maggies, swire_ivar
    gmr=-2.5*alog10(swire_maggies[1,*]/swire_maggies[2,*])
    iblue=where(gmr lt sp.z*3.+0.33 and gmr gt -0.5)
    gmr=gmr[iblue]
    swire_maggies=swire_maggies[*,iblue]
    swire_ivar=swire_ivar[*,iblue]
    swire=swire[iblue]
    objects=objects[iblue]
    sp=sp[iblue]

    ;; limit number 
    indx_s=shuffle_indx(n_elements(swire), num_sub=nswire, seed=seed)
    swire=swire[indx_s]
    sp=sp[indx_s]
    objects=objects[indx_s]
    swire_to_maggies, swire, swire_maggies, swire_ivar
    swire_ivar[0,*]=0.
    swire_ivar[4,*]=0.
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
