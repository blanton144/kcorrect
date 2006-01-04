;+
; NAME:
;   k_nmf_spdata_combo17
; PURPOSE:
;   put together data for nmf fitting to combo17
; CALLING SEQUENCE:
;   k_nmf_spdata_combo17
; COMMENTS:
;   Assumes k_nmf_mmatrix_combo17 has been run.
;
;   Takes data from sources:
;      $KCORRECT_DIR/data/redshifts/combo17
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
; REVISION HISTORY:
;   23-Nov-2004  Michael Blanton (NYU)
;-
;------------------------------------------------------------------------------
pro k_nmf_spdata_combo17, mmatrix=mmatrix, sample=sample, flux=flux, $
                          seed=seed1

spchop=1
if(n_elements(mmatrix) eq 0) then mmatrix='k_nmf_mmatrix.fits'
if(n_elements(outfile) eq 0) then outfile='k_nmf_spdata.fits'


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

c17=mrdfits(getenv('KCORRECT_DIR')+ $
                      '/data/redshifts/combo17/COMBO17_table3.fits',1)
ipick=where(c17.mc_z gt 0.3 and c17.mc_z lt 1.0 and $
            strmid(c17.mc_class, 0, 6) eq 'Galaxy' and $
            c17.e_mc_z gt 0. and c17.e_mc_z lt 0.1 and $ 
            c17.ap_rmag eq c17.ap_rmag, npick)
c17=c17[ipick]
c17_to_maggies, c17, maggies, ivar

;; create full matrix
ntotal=npick
datastr={nx:n_elements(lambda), $
         ny:ntotal, $
         rowstart:lonarr(ntotal), $
         nxrow:lonarr(ntotal)}
currx=0L
xx=0
zdist=fltarr(ntotal)
zhelio=fltarr(ntotal)

ic17=lindgen(ntotal)

zdist[ic17]=c17.mc_z
zhelio[ic17]=c17.mc_z
iz=long(floor((nzf-1.)*(zhelio[ic17]-zf[0])/(zf[nzf-1]-zf[0])+0.5))
for i=0L, n_elements(ic17)-1L do begin
    datastr.rowstart[ic17[i]]=currx
    datastr.nxrow[ic17[i]]=17L
    new_xx=iz[i]+(lindgen(17))*nzf+nspec
    if(keyword_set(xx)) then begin
        xx=[xx, new_xx]
    endif else begin
        xx=new_xx
    endelse
    currx=currx+17L
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
mwrfits, maggies, outfile
mwrfits, ivar, outfile
mwrfits, xx, outfile
mwrfits, lambda, outfile
mwrfits, zdist, outfile
mwrfits, zhelio, outfile

end
