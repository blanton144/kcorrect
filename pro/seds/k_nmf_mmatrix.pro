;+
; NAME:
;   k_nmf_mmatrix
; PURPOSE:
;   make grid of BC03 model spectra for fitting
; CALLING SEQUENCE:
;   k_nmf_mmatrix [, outfile= ]
; COMMENTS:
;   makes [Nobs, Nsfh] matrix where 
;        Nobs = nlambda + nlines + nz*nfilters
; BUGS:
;   needs to output *unsmoothed* spectra 
;   should have better constrained emission lines
; REVISION HISTORY:
;   29-Jul-2004  Michael Blanton (NYU)
;-
;------------------------------------------------------------------------------
pro k_nmf_mmatrix, outfile=outfile

norm_lmin=1500.
norm_lmax=23000.
lmin=600.
lmax=30000.
navloglam=8000L
nagesmax=10
vdisp=300.
minzf=0.
maxzf=1.
nzf=100
;;nmets=6
;;mets=[0,1,2,3,4,5]
nmets=2
mets=[3,4]
sigma=vdisp/(2.99792e+5*alog(10.))

;; 1. make stellar pops

;;     a. find minimum set of bursts to use
bc03= k_im_read_bc03()
iwave=where(bc03.wave gt norm_lmin and bc03.wave lt norm_lmax, nwave)
nages=3000L + nagesmax
tol=5.
while(nages gt nagesmax) do begin
    iuse=n_elements(bc03.age)-1L 
    i=iuse
    while i ge 1 do begin
        j=iuse[0]-1L
        scalefact=total(bc03.flux[iwave,i]*bc03.flux[iwave,j])/ $
          total(bc03.flux[iwave,i]*bc03.flux[iwave,i])
        diff=total(((bc03.flux[iwave,i]*scalefact- $
                     bc03.flux[iwave,j])/ $
                    bc03.flux[iwave,j])^2,/double) 
        while(diff lt tol AND j gt 0L) do begin
            j=j-1L
            scalefact=total(bc03.flux[iwave,i]*bc03.flux[iwave,j])/ $
              total(bc03.flux[iwave,i]*bc03.flux[iwave,i])
            diff=total(((bc03.flux[iwave,i]*scalefact- $
                         bc03.flux[iwave,j])/ $
                        bc03.flux[iwave,j])^2,/double) 
        endwhile
        if(j ge 0 and diff ge tol) then begin
            iuse=[j,iuse]
        endif
        i=j
    endwhile
    tol=tol/0.85
    ages=bc03.age[iuse]
    nages=n_elements(ages)
    help,tol, nages
endwhile

;;     b. now make the full grid for the desired ages
tmp_bc03= k_im_read_bc03(age=1.)
nl=n_elements(tmp_bc03.flux)
wave=tmp_bc03.wave
loglam=alog10(wave)
grid=fltarr(nl, nages, nmets)
for im= 0L, nmets-1L do $
  grid[*,*,im]= (k_im_read_bc03(met=mets[im])).flux[*,iuse]

;;     c. put flux units into absolute maggies
dfact=3.826/3.086*1.e-5
absrc=3.631*2.99792*1.e-2/wave^2
grid=grid*dfact
for ia=0L, nages-1L do $
  for im= 0L, nmets-1L do $
  grid[*,ia,im]=grid[*,ia,im]*dfact/absrc

;;     d. interpolate grid onto flux grid
avloglam=double(alog10(lmin)+(alog10(lmax)-alog10(lmin))* $
                (dindgen(navloglam)+0.5)/float(navloglam))
sfgrid=fltarr(navloglam, nages*nmets)
ninterloglam=20000L
interloglam=double(alog10(lmin)+(alog10(lmax)-alog10(lmin))* $
                   (dindgen(ninterloglam)+0.5)/float(ninterloglam))
for im= 0L, nmets-1L do begin
    for ia= 0L, nages-1L do begin 
        splog, string(im)+string(ia) 
        intergrid=interpol(grid[*,ia,im], loglam, interloglam) 
        combine1fiber, interloglam, intergrid, fltarr(ninterloglam)+1., $
          newloglam=avloglam, newflux=tmp1, maxiter=0 
        sfgrid[*,ia+im*nages]=tmp1 
    endfor
endfor
  
;;     e. smooth to desired vdisp
for im= 0L, nmets-1L do begin
    for ia= 0L, nages-1L do begin 
        sfgrid[*,ia+im*nages]= $
          gauss_smooth(avloglam, sfgrid[*,ia+im*nages], sigma, avloglam)
    endfor
endfor

;; 2. make the dusty grid
ndusts=3L
dusts1={geometry:'', dust:'', structure:'', tauv:0.}
dusts=replicate(dusts1,3)
dusts.geometry=['dusty', $
               'dusty','dusty']
dusts.dust=['MW', $
           'MW','MW']
dusts.structure=['c', $
                'c','c']
dusts.tauv=[0.,1.,3.]
dustygrid=fltarr(navloglam, nages*nmets, ndusts)
dustfact=fltarr(navloglam, ndusts)
for i=0L, ndusts-1L do $
  dustfact[*,i]=exp(-witt_ext(dusts[i],dusts[i].tauv,10.^(avloglam)))
for i=0L, ndusts-1L do $
  for j=0L, nages*nmets-1L do $
  dustygrid[*,j,i]=sfgrid[*,j]*dustfact[*,i]

;; 3. make lookup table for properties
spgrid=reform(dustygrid,navloglam,nages*nmets*ndusts)
dust=replicate(dusts[0], nages,nmets,ndusts)
for i=0L, ndusts-1L do dust[*,*,i]=dusts[i]
met=fltarr(nages,nmets,ndusts)
for i=0L, nmets-1L do met[*,i,*]=mets[i]
age=fltarr(nages,nmets,ndusts)
for i=0L, nages-1L do age[i,*,*]=ages[i]

;; now make emission lines 
nel=0
if(1) then begin
    readcol, getenv('KCORRECT_DIR')+'/data/templates/linelist.txt', $
      f='(a,f,a)', comment=';', elname, lambda, type
    ii=where(type eq 'em' OR type eq 'both' and $
             lambda gt 3700. and $
             lambda lt 9000.)
    elname=elname[ii]
    lambda=lambda[ii]
    nel=n_elements(elname)
    emgrid=fltarr(navloglam, nel)
    absrc=3.631*2.99792*1.e-2/10.^(avloglam*2.)
    for i=0L, nel-1L do $
      emgrid[*, i]= 1.e-14*exp(-(avloglam-alog10(lambda[i]))^2/(sigma)^2)/ $
      (sqrt(2.*!DPI)*sigma)/absrc
    tmp_spgrid=spgrid
    spgrid=fltarr(navloglam,nages*nmets*ndusts+nel)
    spgrid[*,0L:nages*nmets*ndusts-1L]=tmp_spgrid
    spgrid[*,nages*nmets*ndusts:nages*nmets*ndusts+nel-1L]=emgrid
endif

;; 3. now make all filters at all redshifts
filterlist=['galex_FUV.par', $
            'galex_NUV.par', $
            'sdss_u0.par', $
            'sdss_g0.par', $
            'sdss_r0.par', $
            'sdss_i0.par', $
            'sdss_z0.par', $
            'twomass_J.par', $
            'twomass_H.par', $
            'twomass_Ks.par']
lambda=fltarr(navloglam+1L)
davloglam=avloglam[1]-avloglam[0]
lambda[0:navloglam-1L]=10.^(avloglam-davloglam)
lambda[1:navloglam]=10.^(avloglam+davloglam)
absrc=3.631*2.99792*1.e-2/10.^(avloglam*2.)
pgrid=spgrid
for i=0L, nages*ndusts*nmets+nel-1L do $
  pgrid[*,i]=pgrid[*,i]*absrc
k_projection_table, rmatrix, pgrid, lambda, zf, filterlist, zmin=minzf, $
  zmax=maxzf, nz=nzf
rmatrix=rmatrix > 0.
  
;; 4. prepare output
if(NOT keyword_set(outfile)) then outfile='k_nmf_mmatrix.fits'

outgrid=fltarr(navloglam+nzf*n_elements(filterlist),nages*nmets*ndusts+nel)
outgrid[0:navloglam-1L,*]=spgrid
for i=0L, (nages*nmets*ndusts)+nel-1L do $
  outgrid[navloglam:navloglam+nzf*n_elements(filterlist)-1L,i]= $
  rmatrix[*,i,*]

outlambda=fltarr(navloglam+nzf*n_elements(filterlist))
outlambda[0:navloglam-1]=10.^(avloglam)
outlambda[navloglam:navloglam+nzf*n_elements(filterlist)-1]= $
  (1./(1.+zf))#k_lambda_eff(filterlist=filterlist)

hdr=['']
sxaddpar, hdr, 'NSPEC', navloglam, 'number of points in spectrum'
sxaddpar, hdr, 'NZ', nzf, 'number of redshifts'
sxaddpar, hdr, 'NFILTER', n_elements(filterlist), 'number of filters'
sxaddpar, hdr, 'NDUST', ndusts, 'number of dusts'
sxaddpar, hdr, 'NMET', nmets, 'number of metallicities'
sxaddpar, hdr, 'NAGE', nages, 'number of ages'
sxaddpar, hdr, 'NEL', nel, 'number of emission lines'
sxaddpar, hdr, 'VDISP', vdisp, 'smoothed to this velocity dispersion (km/s)'
mwrfits, outgrid, outfile, hdr, /create
mwrfits, outlambda,outfile 
mwrfits, dust,outfile 
mwrfits, met, outfile
mwrfits, age, outfile
mwrfits, filterlist, outfile
mwrfits, zf, outfile

end