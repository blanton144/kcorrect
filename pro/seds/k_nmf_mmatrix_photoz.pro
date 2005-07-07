;+
; NAME:
;   k_nmf_mmatrix_photoz
; PURPOSE:
;   make grid of BC03 model spectra for fitting
; CALLING SEQUENCE:
;   k_nmf_mmatrix [, prefix=, back=, lmin=, lmax=, dusts= ]
; OPTIONAL INPUTS:
;   prefix - prefix to use for output files (default 'k_nmf')
;   back - # of Gyrs in the past to use for 'early' file
;   lmin, lmax - limits of spectrum output (Angstroms; default 600, 30000)
;   navloglam - number of wavelengths in final output
;   nagesmax - maximum number of instantantaneous burst ages to use
;   vdisp - smooth models to this velocity dispersion, in km/s (default 300)
;   minzf, maxzf - minimum and maximum redshifts of observation
;                  (default 0., 1.)
;   nzf - number of redshifts of observations
;   dusts - dust models (Witt & Gordon) to use for
;           extinction. Structure with elements:
;              .GEOMETRY - geometry from WG (default 'dusty', 'dusty', 'dusty')
;              .DUST - type of dust from WG (default 'MW', 'MW', 'MW')
;              .STRUCTURE - structure of dust from WG (default 'c', 'c', 'c')
;              .TAUV - optical depth at V (default 0., 1., 3.)
; COMMENTS:
;   makes [Nobs, Nsfh] matrix where 
;        Nobs = nlambda + nz*nfilters
;   units are absolute maggies per solar mass
; BUGS:
;   needs to output *unsmoothed* spectra 
;   should have better constrained emission lines
;   make sure DEEP filters are fixed 
;   include Draine & Li templates
; REVISION HISTORY:
;   29-Jul-2004  Michael Blanton (NYU)
;-
;------------------------------------------------------------------------------
pro k_nmf_mmatrix_photoz, prefix=prefix, back=back, lmin=lmin, lmax=lmax, $
                          dusts=dusts , filterlist=filterlist, minzf=minzf, $
                          maxzf=maxzf, nzf=nzf, nagesmax=nagesmax, noel=noel, $
                          isolib=isolib, nodust=nodust

nodust=1
noel=1
if(NOT keyword_set(isolib)) then isolib='Padova1994'
if(NOT keyword_set(back)) then back=0.5  ;; how many Gyrs earlier?
if(NOT keyword_set(lmin)) then lmin=600.
if(NOT keyword_set(lmax)) then lmax=3200000.
if(NOT keyword_set(prefix)) then prefix='k_nmf'
if(NOT keyword_set(navloglam)) then navloglam=10000L
if(NOT keyword_set(nagesmax)) then nagesmax=5
if(NOT keyword_set(vdisp)) then vdisp=300.
if(NOT keyword_set(minzf)) then minzf=0.
if(NOT keyword_set(maxzf)) then maxzf=2.
if(NOT keyword_set(nzf)) then nzf=300
if(NOT keyword_set(narbitrary)) then $
  narbitrary=long((alog10(lmax)-alog10(lmin))/0.06)
if(NOT keyword_set(filterlist)) then $
  filterlist=['galex_FUV.par', $
              'galex_NUV.par', $
              'sdss_u0.par', $
              'sdss_g0.par', $
              'sdss_r0.par', $
              'sdss_i0.par', $
              'sdss_z0.par', $
              'twomass_J.par', $
              'twomass_H.par', $
              'twomass_Ks.par', $
              'deep_B.par', $
              'deep_R.par', $
              'deep_I.par', $
              'goods_J_isaac_etc.par', $
              'goods_H_isaac_etc.par', $
              'goods_Ks_isaac_etc.par', $
              'goods_acs_f435w.par', $
              'goods_acs_f606w.par', $
              'goods_acs_f775w.par', $
              'goods_acs_f850lp.par', $
              'spitzer_irac_ch1.par', $
              'spitzer_irac_ch2.par', $
              'spitzer_irac_ch3.par', $
              'spitzer_irac_ch4.par', $
              'spitzer_mips_24.par', $
              'spitzer_mips_70.par', $
              'spitzer_mips_160.par']

if(n_tags(dusts) eq 0) then begin
    dusts1={geometry:'', dust:'', structure:'', tauv:0.}
    dusts=replicate(dusts1,3)
    dusts.geometry=['shell','shell','shell']
    dusts.dust=['MW','MW','SMC']
    dusts.structure=['h','h','h']
    dusts.tauv=[0.,3.0,3.0]
endif 
if(keyword_set(nodust)) then begin
    dusts={geometry:'dusty', dust:'MW', structure:'c', tauv:0.}
endif
ndusts=n_elements(dusts)

norm_lmin=800. ;; limits when testing whether two BC03 models are similar
norm_lmax=23000.
;; metallicities to use
mets=[0,1,2,3,4,5]
mets=[4]
nmets=n_elements(mets)   
sigma=vdisp/(2.99792e+5*alog(10.))  ;; smoothing sigma in log lambda
outfile=prefix+'_mmatrix.fits'  ;; output files
earlyfile=prefix+'_early.fits'
rawfile=prefix+'_rawspec.fits'

;; 1. make stellar pops

;;     a. find minimum set of bursts to use
bc03= k_im_read_bc03(isolib=isolib, /vac)
imaxage=max(where(bc03.age lt 14.e+9))
iwave=where(bc03.wave gt norm_lmin and bc03.wave lt norm_lmax, nwave)
nages=3000L + nagesmax
tol=5.
while(nages gt nagesmax) do begin
    iuse=imaxage 
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
tmp_bc03= k_im_read_bc03(age=1.,isolib=isolib, /vac)
nl=n_elements(tmp_bc03.flux)
wave=tmp_bc03.wave
loglam=alog10(wave)
grid=fltarr(nl, nages, nmets)
for im= 0L, nmets-1L do $
  grid[*,*,im]= (k_im_read_bc03(met=mets[im],isolib=isolib, /vac)).flux[*,iuse]
earlygrid=fltarr(nl, nages, nmets)
early=ages-back*1.e+9
iearly=where(early gt 0., nearly)
for ia= 0L, nearly-1L do $
  for im= 0L, nmets-1L do $
  earlygrid[*,iearly[ia],im]= $
  (k_im_read_bc03(met=mets[im], $
                  age=early[iearly[ia]]/1.e+9,isolib=isolib, /vac)).flux

;;     d. interpolate grid onto flux grid
avloglam=double(alog10(lmin)+(alog10(lmax)-alog10(lmin))* $
                (dindgen(navloglam)+0.5)/float(navloglam))
sfgrid=fltarr(navloglam, nages*nmets)
earlysfgrid=fltarr(navloglam, nages*nmets)
ninterloglam=20000L
interloglam=double(alog10(lmin)+(alog10(lmax)-alog10(lmin))* $
                   (dindgen(ninterloglam)+0.5)/float(ninterloglam))
for im= 0L, nmets-1L do begin
    for ia= 0L, nages-1L do begin 
        splog, string(im)+string(ia) 
        dloglam=loglam[1]-loglam[0]
        nloglam=n_elements(loglam)
        intergrid=interpol([0., 0., grid[*,ia,im], 0., 0.], $
                           [loglam[0]-dloglam, loglam[0], loglam, $
                            loglam[nloglam-1L], loglam[nloglam-1]+dloglam], $
                           interloglam) 
        combine1fiber, interloglam, intergrid, fltarr(ninterloglam)+1., $
          newloglam=avloglam, newflux=tmp1, maxiter=0 
        sfgrid[*,ia+im*nages]=tmp1 
        if(early[ia] gt 0.) then begin
            intergrid=interpol(earlygrid[*,ia,im], loglam, interloglam) 
            combine1fiber, interloglam, intergrid, fltarr(ninterloglam)+1., $
              newloglam=avloglam, newflux=tmp1, maxiter=0 
            earlysfgrid[*,ia+im*nages]=tmp1 
        endif
    endfor
endfor

;;     c. put flux units into absolute maggies
;; from solar luminosities to ergs s-1 cm-2 A-1
;;  (don't forget that 10 pc != 1 pc)
dfact=3.826/(4.*!DPI)/(3.086)^2*1.e-5 
absrc=3.631*2.99792*1.e-2/(10.^avloglam)^2
sfgrid=sfgrid*dfact
earlysfgrid=earlysfgrid*dfact
rawgrid=sfgrid
for ia=0L, nages-1L do $
  for im= 0L, nmets-1L do $
  sfgrid[*,ia+im*nages]=sfgrid[*,ia+im*nages]/absrc
for ia=0L, nages-1L do $
  for im= 0L, nmets-1L do $
  earlysfgrid[*,ia+im*nages]=earlysfgrid[*,ia+im*nages]/absrc

;;     e. smooth to desired vdisp
for im= 0L, nmets-1L do begin
    for ia= 0L, nages-1L do begin 
        sfgrid[*,ia+im*nages]= $
          k_smooth(avloglam, sfgrid[*,ia+im*nages], vdisp)
        earlysfgrid[*,ia+im*nages]= $
          k_smooth(avloglam, earlysfgrid[*,ia+im*nages], vdisp)
    endfor
endfor

;; 2. make the dusty grid
dustygrid=fltarr(navloglam, nages*nmets, ndusts)
earlydustygrid=fltarr(navloglam, nages*nmets, ndusts)
rawdustygrid=fltarr(navloglam, nages*nmets, ndusts)
dustfact=fltarr(navloglam, ndusts)
for i=0L, ndusts-1L do $
  dustfact[*,i]=exp(-witt_ext(dusts[i],dusts[i].tauv,10.^(avloglam)))
for i=0L, ndusts-1L do $
  for j=0L, nages*nmets-1L do $
  dustygrid[*,j,i]=sfgrid[*,j]*dustfact[*,i]
for i=0L, ndusts-1L do $
  for j=0L, nages*nmets-1L do $
  earlydustygrid[*,j,i]=earlysfgrid[*,j]*dustfact[*,i]
for i=0L, ndusts-1L do $
  for j=0L, nages*nmets-1L do $
  rawdustygrid[*,j,i]=rawgrid[*,j]*dustfact[*,i]

;; 3. make lookup table for properties
spgrid=reform(dustygrid,navloglam,nages*nmets*ndusts)
earlyspgrid=reform(earlydustygrid,navloglam,nages*nmets*ndusts)
rawspgrid=reform(rawdustygrid,navloglam,nages*nmets*ndusts)
dust=replicate(dusts[0], nages,nmets,ndusts)
for i=0L, ndusts-1L do dust[*,*,i]=dusts[i]
met=fltarr(nages,nmets,ndusts)
for i=0L, nmets-1L do met[*,i,*]=mets[i]
age=fltarr(nages,nmets,ndusts)
for i=0L, nages-1L do age[i,*,*]=ages[i]

;; now make emission lines 
nextra=0
nel=0
if(NOT keyword_set(noel)) then begin
    if(0) then begin
        readcol, getenv('KCORRECT_DIR')+'/data/templates/linelist.txt', $
          f='(a,f,a)', comment=';', elname, lambda, type
        ii=where(type eq 'em' OR type eq 'both' and $
             lambda gt 3700. and $
                 lambda lt 9000.)
        elname=elname[ii]
        lambda=lambda[ii]
        nel=n_elements(elname)
        emgrid=fltarr(navloglam, nel)
        for i=0L, nel-1L do $
          emgrid[*, i]= 1.e-7*exp(-(avloglam-alog10(lambda[i]))^2/(sigma)^2)/ $
          (sqrt(2.*!DPI)*sigma)/absrc
        tmp_spgrid=spgrid
        spgrid=fltarr(navloglam,nages*nmets*ndusts+nel)
        spgrid[*,0L:nages*nmets*ndusts-1L]=tmp_spgrid
        spgrid[*,nages*nmets*ndusts:nages*nmets*ndusts+nel-1L]=emgrid
        tmp_earlyspgrid=earlyspgrid
        earlyspgrid=fltarr(navloglam,nages*nmets*ndusts+nel)
        earlyspgrid[*,0L:nages*nmets*ndusts-1L]=tmp_earlyspgrid
        earlyspgrid[*,nages*nmets*ndusts:nages*nmets*ndusts+nel-1L]=emgrid
    endif
    
    atm='K'
    agestr='8Myr'
    model='SB99' 
    sfh='cont_n10'
    gasmets=['.05', '0.2', '0.4', '1.0', '2.0']
    qpars=[ '5.0e6', '1.0e7', '2.0e7', '4.0e7', '8.0e7', '1.5e8', '3.0e8'] 
    nel=n_elements(gasmets)*n_elements(qpars)
    emgrid=fltarr(navloglam, nel)
    for i=0L, n_elements(gasmets)-1L do begin
        for j=0L, n_elements(qpars)-1L do begin
            gasmet=gasmets[i]
            qpar=qpars[j]
            filename=getenv('KCORRECT_DIR')+'/data/seds/mappings/'+ $
              model+'_'+sfh+'/spec_Z'+gasmet+'_'+agestr+'_q'+qpar+'_'+model+ $
              '_'+atm+'.ph4'
            mappings=read_mappings(filename, /vac)
            for k=0L, n_elements(mappings.lambda)-1L do $
              emgrid[*, i*n_elements(qpars)+j]= $
              emgrid[*, i*n_elements(qpars)+j]+1.e-10* $
              mappings.flux[k]*exp(-(avloglam-alog10(mappings.lambda[k]))^2/ $
                                   (sigma)^2)/(sqrt(2.*!DPI)*sigma)/absrc
        endfor
    endfor
    tmp_spgrid=spgrid
    spgrid=fltarr(navloglam,nages*nmets*ndusts+nel)
    spgrid[*,0L:nages*nmets*ndusts-1L]=tmp_spgrid
    spgrid[*,nages*nmets*ndusts:nages*nmets*ndusts+nel-1L]=emgrid
    tmp_earlyspgrid=earlyspgrid
    earlyspgrid=fltarr(navloglam,nages*nmets*ndusts+nel)
    earlyspgrid[*,0L:nages*nmets*ndusts-1L]=tmp_earlyspgrid
    earlyspgrid[*,nages*nmets*ndusts:nages*nmets*ndusts+nel-1L]=emgrid
endif
nextra=nextra+nel

;; 3.5 dust from draine
ndraine=0L 
if(0) then begin
drainefiles=['spec_2.2.dat','spec_2.5.dat','spec_2.8.dat']
ndraine=n_elements(drainefiles)
drainegrid=fltarr(navloglam, ndraine)
for i=0L, n_elements(drainefiles)-1L do begin
    dst=read_draine(getenv('KCORRECT_DIR')+'/data/seds/draine/'+ $
                    drainefiles[i])
    lambda=[dst.lambda[0]-1., dst.lambda[0], dst.lambda, $
            dst.lambda[n_elements(dst.lambda)-1L], $
            dst.lambda[n_elements(dst.lambda)-1L]+1. ]
    flux=[0., 0., dst.flux, 0., 0.]
    drainegrid[*,i]= 1.e+18*interpol(flux, alog10(lambda), avloglam)/absrc
endfor
tmp_spgrid=spgrid
spgrid=fltarr(navloglam,nages*nmets*ndusts+nextra+ndraine)
spgrid[*,0L:nages*nmets*ndusts+nextra-1L]=tmp_spgrid
spgrid[*,nages*nmets*ndusts+nextra: $
       nages*nmets*ndusts+nextra+ndraine-1L]=drainegrid
tmp_earlyspgrid=earlyspgrid
earlyspgrid=fltarr(navloglam,nages*nmets*ndusts+nextra+ndraine)
earlyspgrid[*,0L:nages*nmets*ndusts+nextra-1L]=tmp_earlyspgrid
earlyspgrid[*,nages*nmets*ndusts+nextra: $
            nages*nmets*ndusts+nextra+ndraine-1L]= drainegrid
endif
nextra=nextra+ndraine

;; 3.7 add totally arbitrary components!
sarb=0.005
arblo=1200.
arbhi=24000.
narb=long((alog10(arbhi)-alog10(arblo))/sarb)
carb=alog10(arblo)+(alog10(arbhi)-alog10(arblo))*(findgen(narb)+0.5)/ $
  float(narb)
rsarb=carb[1]-carb[0]
arbgrid=fltarr(navloglam, narb)
for i=0L, narb-1L do $
  arbgrid[*,i]=k_bspline2((avloglam-carb[i])/rsarb)
tmp_spgrid=spgrid
spgrid=fltarr(navloglam,nages*nmets*ndusts+nextra+narb)
spgrid[*,0L:nages*nmets*ndusts+nextra-1L]=tmp_spgrid
spgrid[*,nages*nmets*ndusts+nextra: $
       nages*nmets*ndusts+nextra+narb-1L]=arbgrid
tmp_earlyspgrid=earlyspgrid
earlyspgrid=fltarr(navloglam,nages*nmets*ndusts+nextra+narb)
earlyspgrid[*,0L:nages*nmets*ndusts+nextra-1L]=tmp_earlyspgrid
earlyspgrid[*,nages*nmets*ndusts+nextra: $
            nages*nmets*ndusts+nextra+narb-1L]= arbgrid
nextra=nextra+narb

;; 3. now make all filters at all redshifts
lambda=fltarr(navloglam+1L)
davloglam=avloglam[1]-avloglam[0]
lambda[0:navloglam-1L]=10.^(avloglam-davloglam)
lambda[1:navloglam]=10.^(avloglam+davloglam)
pgrid=spgrid
for i=0L, nages*ndusts*nmets+nextra-1L do $
  pgrid[*,i]=pgrid[*,i]*absrc
k_projection_table, rmatrix, pgrid, lambda, zf, filterlist, zmin=minzf, $
  zmax=maxzf, nz=nzf
rmatrix=rmatrix > 0.
earlypgrid=earlyspgrid
for i=0L, nages*ndusts*nmets+nextra-1L do $
  earlypgrid[*,i]=earlypgrid[*,i]*absrc
k_projection_table, earlyrmatrix, earlypgrid, lambda, zf, filterlist, $
  zmin=minzf, zmax=maxzf, nz=nzf
earlyrmatrix=earlyrmatrix > 0.
  
;; 4. prepare output

outgrid=fltarr(navloglam+nzf*n_elements(filterlist),nages*nmets*ndusts+nextra)
outgrid[0:navloglam-1L,*]=spgrid
for i=0L, (nages*nmets*ndusts)+nextra-1L do $
  outgrid[navloglam:navloglam+nzf*n_elements(filterlist)-1L,i]= $
  rmatrix[*,i,*]

earlyoutgrid=fltarr(navloglam+nzf*n_elements(filterlist),nages*nmets*ndusts+ $
                    nextra)
earlyoutgrid[0:navloglam-1L,*]=earlyspgrid
for i=0L, (nages*nmets*ndusts)+nextra-1L do $
  earlyoutgrid[navloglam:navloglam+nzf*n_elements(filterlist)-1L,i]= $
  earlyrmatrix[*,i,*]

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
sxaddpar, hdr, 'NEXTRA', nextra, 'number of extra templates'
sxaddpar, hdr, 'NEL', nel, 'number of emission lines'
sxaddpar, hdr, 'NDRAINE', ndraine, 'number of draine models'
sxaddpar, hdr, 'ISOLIB', isolib, 'isochrone library used'
sxaddpar, hdr, 'VDISP', vdisp, 'smoothed to this velocity dispersion (km/s)'
mwrfits, outgrid, outfile, hdr, /create
mwrfits, outlambda,outfile 
mwrfits, dust,outfile 
mwrfits, met, outfile
mwrfits, age, outfile
mwrfits, filterlist, outfile
mwrfits, zf, outfile
mwrfits, gasmets, outfile
mwrfits, qpars, outfile
if(ndraine gt 0) then $
mwrfits, drainefiles, outfile

hdr=['']
sxaddpar, hdr, 'NSPEC', navloglam, 'number of points in spectrum'
sxaddpar, hdr, 'NZ', nzf, 'number of redshifts'
sxaddpar, hdr, 'NFILTER', n_elements(filterlist), 'number of filters'
sxaddpar, hdr, 'NDUST', ndusts, 'number of dusts'
sxaddpar, hdr, 'NMET', nmets, 'number of metallicities'
sxaddpar, hdr, 'NAGE', nages, 'number of ages'
sxaddpar, hdr, 'NEXTRA', nextra, 'number of extra templates'
sxaddpar, hdr, 'NEL', nel, 'number of emission lines'
sxaddpar, hdr, 'NDRAINE', ndraine, 'number of draine models'
sxaddpar, hdr, 'VDISP', vdisp, 'smoothed to this velocity dispersion (km/s)'
sxaddpar, hdr, 'BACK', back, 'Gyrs previous to mmatrix'
mwrfits, earlyoutgrid, earlyfile, hdr, /create

hdr=['']
sxaddpar, hdr, 'NSPEC', navloglam, 'number of points in spectrum'
sxaddpar, hdr, 'NDUST', ndusts, 'number of dusts'
sxaddpar, hdr, 'NMET', nmets, 'number of metallicities'
sxaddpar, hdr, 'NAGE', nages, 'number of ages'
sxaddpar, hdr, 'NEXTRA', nextra, 'number of extra templates'
sxaddpar, hdr, 'NEL', nel, 'number of emission lines'
sxaddpar, hdr, 'NDRAINE', ndraine, 'number of draine models'
mwrfits, rawspgrid, rawfile, hdr, /create

end
