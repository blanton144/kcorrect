;+
; NAME:
;   k_qa_nmf
; COMMENTS:
;   Creates k_qa_nmf.ps, which has a number of QA plots on the NMF 
;     fitting results, including colors as a function of redshift,
;     color-color plots, spectra vs. model spectra, star-formation 
;     histories, etc.
;   Also determines:
;     mass-weighted metallicity of each template
;     mass-weighted age of each template 
;     passive evolution estimates for each template, assuming constant
;       SFR into the future, based on looking 0.5 Gyr younger
;     makes standard vmatrix/lambda files 
; REVISION HISTORY:
;   30-Nov-2004  Michael Blanton (NYU)
;-
;------------------------------------------------------------------------------
pro k_qa_nmf

version='test'
nsubsp=130

;; read in the template basics
mmatrix=mrdfits('k_nmf_mmatrix.fits',0,hdr)
lambda=mrdfits('k_nmf_mmatrix.fits',1)
dust=mrdfits('k_nmf_mmatrix.fits',2)
met=mrdfits('k_nmf_mmatrix.fits',3)
age=mrdfits('k_nmf_mmatrix.fits',4)
rawspec=mrdfits('k_nmf_rawspec.fits',0)
filterlist=string(mrdfits('k_nmf_mmatrix.fits',5))
early=mrdfits('k_nmf_early.fits',0,earlyhdr)
zf=mrdfits('k_nmf_mmatrix.fits',6)
nspec=long(sxpar(hdr, 'NSPEC'))
back=long(sxpar(hdr, 'BACK'))
nel=long(sxpar(hdr, 'NEL'))
nzf=long(sxpar(hdr, 'NZ'))
nfilter=long(sxpar(hdr, 'NFILTER'))
ndusts=long(sxpar(hdr, 'NDUST'))
nmets=long(sxpar(hdr, 'NMET'))
nages=long(sxpar(hdr, 'NAGE'))
filternames=['F', 'N', 'u', 'g', 'r', 'i', 'z', 'J', 'H', 'K_s']
metallicities=[0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05]

;; read in the data
data=mrdfits('k_nmf_data.fits',1)
ivar=mrdfits('k_nmf_data.fits',2)
zhelio=mrdfits('k_nmf_data.fits',5)
iz=long(floor((nzf-1.)*(zhelio-zf[0])/(zf[nzf-1]-zf[0])+0.5))

;; read in the results
templates=mrdfits('k_nmf_soln.fits')
coeffs=mrdfits('k_nmf_soln.fits',1)
nt=(size(templates, /dim))[1]
model=mmatrix#templates#coeffs
emodel=early#templates#coeffs

;; determine basic derived quantities from the templates
t_mass=total(templates[0:nages*nmets*ndusts-1L,*], 1)
splog, t_mass
t_metallicity= $
  total((reform(metallicities[met], nages*nmets*ndusts, 1)# $
         replicate(1., nt))* $
        templates[0:nages*nmets*ndusts-1L,*], 1)/t_mass
splog, t_metallicity
t_age=total((reform(age, nages*nmets*ndusts, 1)# $ 
             replicate(1., nt))* $
            templates[0:nages*nmets*ndusts-1L,*], 1)/t_mass
splog, t_age
mwrfits, t_mass, 'k_nmf_derived.fits', hdr, /create
mwrfits, t_metallicity, 'k_nmf_derived.fits'
mwrfits, t_age, 'k_nmf_derived.fits'

;; make vmatrix and lambda
outvmatrix=mmatrix[0:nspec-1L,*]#templates
absrc=3.631*2.99792*1.e-2/(lambda[0:nspec-1L])^2
for i=0L, nt-1L do $
  outvmatrix[*,i]=outvmatrix[*,i]*absrc
outlambda=fltarr(nspec+1L)
dlg10l=alog10(lambda[1])-alog10(lambda[0])
outlambda[0:nspec-1L]= 10.^(alog10(lambda[0:nspec-1L])-0.5*dlg10l)
outlambda[nspec]= 10.^(alog10(lambda[nspec-1L])+0.5*dlg10l)
k_write_ascii_table,outvmatrix,'vmatrix.'+version+'.dat'
k_write_ascii_table,outlambda,'lambda.'+version+'.dat'

;; create broadband model fluxes + data
mfluxes=fltarr(nfilter, n_elements(zhelio))
fluxes=fltarr(nfilter, n_elements(zhelio))
fluxes_ivar=fltarr(nfilter, n_elements(zhelio))
for ifilter=0L, nfilter-1L do $
  mfluxes[ifilter,*]= model[iz+ifilter*nzf+nspec, $
                            lindgen(n_elements(zhelio))]
for ifilter=0L, nfilter-1L do $
  fluxes[ifilter,*]= data[iz+ifilter*nzf+nspec, $
                          lindgen(n_elements(zhelio))]
for ifilter=0L, nfilter-1L do $
  fluxes_ivar[ifilter,*]= ivar[iz+ifilter*nzf+nspec, $
                               lindgen(n_elements(zhelio))]

set_print, filename='k_qa_nmf.ps'

;; for each template, show the spectrum (full and just optical)
!P.MULTI=[0, 1, 2]
tspec=mmatrix#templates
for i=0L, nt-1L do begin
    djs_plot, lambda[0:nspec-1], tspec[*,i], /xlog, /ylog
    il=where(lambda[0:nspec-1] gt 3000. and lambda[0:nspec-1] lt 8500., nl)
    djs_plot, lambda[il], tspec[il, i], /ylog
endfor

;; for each template, the star-formation history and (mean metallicity
;; and dustiness of same) 
dage=fltarr(nages)
dage[1:nages-2L]=0.5*(age[2:nages-1L]-age[0:nages-3L])
dage[0]=age[1]-age[0]
dage[nages-1L]=age[nages-1L]-age[nages-2L]
!P.MULTI=[0, 1, 3]
for i=0L, nt-1L do begin
    sfh=reform(templates[0:nages*nmets*ndusts-1L,i], nages, nmets, ndusts)
    dust=reform(dust, nages, nmets, ndusts)
    sfh_tot=fltarr(nages)
    for j=0L, nages-1L do sfh_tot[j]=total(sfh[j,*,*])
    sfh_met=fltarr(nages)
    for j=0L, nages-1L do sfh_met[j]=total(sfh[j,*,*]* $
                                           metallicities[met[j,*,*]])
    sfh_met=sfh_met/sfh_tot
    sfh_dust=fltarr(nages)
    for j=0L, nages-1L do sfh_dust[j]=total(sfh[j,*,*]* $
                                            dust[j,*,*].tauv)
    sfh_dust=sfh_dust/sfh_tot
    djs_plot, age[*,0,0], sfh_tot, ytitle='SFR', /xlog, /ylog
    djs_oplot, age[*,0,0], sfh_tot, psym=4
    djs_plot, age[*,0,0], sfh_met, ytitle='metallicity', /xlog
    djs_oplot, age[*,0,0], sfh_met, psym=4
    djs_plot, age[*,0,0], sfh_dust, ytitle='avg \tau_V', /xlog
    djs_oplot, age[*,0,0], sfh_dust, psym=4
endfor

;; colors vs. redshift
!P.MULTI=[0,1,2]
for ifilter=0L, nfilter-2L do begin
    igood=where(fluxes[ifilter,*] gt 0. and $
                fluxes[ifilter+1L,*] gt 0.)
    mcolor=-2.5*alog10(mfluxes[ifilter,igood]/mfluxes[ifilter+1,igood])
    color=-2.5*alog10(fluxes[ifilter,igood]/fluxes[ifilter+1,igood])
    colorerr=2.5/alog(10.)*sqrt(1./fluxes_ivar[ifilter,igood]/ $
                                fluxes[ifilter,igood]^2+ $
                                1./fluxes_ivar[ifilter+1,igood]/ $
                                fluxes[ifilter+1,igood]^2)
    djs_plot, zhelio[igood], mcolor, psym=4, $
      ytitle=filternames[ifilter]+'-'+filternames[ifilter+1]
    djs_oplot, zhelio[igood], color, psym=6, color='red'
    djs_oploterr, zhelio[igood], color, yerr=colorerr, color='red'
    djs_plot, zhelio[igood], color-mcolor, psym=6, $
      ytitle=filternames[ifilter]+'-'+filternames[ifilter+1]+' residual', $
      xtitle='redshift z', color='red'
    djs_oplot, minmax(zhelio), [0., 0.], color='blue', linestyle=1
    djs_oploterr, zhelio[igood], color-mcolor, yerr=colorerr, color='red'
endfor

;; colors vs. color
!P.MULTI=[0,1,1]
for ifilter=0L, nfilter-3L do begin
    igood=where(fluxes[ifilter,*] gt 0. and $
                fluxes[ifilter+1L,*] gt 0.)
    mcolor0=-2.5*alog10(mfluxes[ifilter,igood]/mfluxes[ifilter+1,igood])
    color0=-2.5*alog10(fluxes[ifilter,igood]/fluxes[ifilter+1,igood])
    colorerr0=2.5/alog(10.)*sqrt(1./fluxes_ivar[ifilter,igood]/ $
                                 fluxes[ifilter,igood]^2+ $
                                 1./fluxes_ivar[ifilter+1,igood]/ $
                                 fluxes[ifilter+1,igood]^2)
    mcolor1=-2.5*alog10(mfluxes[ifilter+1,igood]/mfluxes[ifilter+2,igood])
    color1=-2.5*alog10(fluxes[ifilter+1,igood]/fluxes[ifilter+2,igood])
    colorerr1=2.5/alog(10.)*sqrt(1./fluxes_ivar[ifilter+1,igood]/ $
                                 fluxes[ifilter+1,igood]^2+ $
                                 1./fluxes_ivar[ifilter+2,igood]/ $
                                 fluxes[ifilter+2,igood]^2)
    djs_plot, color0, color1, psym=4, $
      xtitle=filternames[ifilter]+'-'+filternames[ifilter+1], $
      ytitle=filternames[ifilter+1]+'-'+filternames[ifilter+2]
    djs_oploterr, color0, color1, xerr=colorerr0, yerr=colorerr1
    djs_oplot, mcolor0, mcolor1, psym=4, color='red'
endfor

;; look at K-corrections for sanity check
kcorrect, fluxes, fluxes_ivar, zhelio, kcorrect, $
  vfile='vmatrix.test.dat', lfile='lambda.test.dat', $
  filterlist=filterlist
for i=0L, n_elements(filterlist)-1L do begin
    djs_plot, zhelio, kcorrect[i,*], psym=4
endfor

;; random set of spectra
iran=shuffle_indx(n_elements(zhelio), num_sub=nsubsp)
for i=0L, nsubsp-1L do begin
    djs_plot, lambda[0:nspec-1], model[0:nspec-1,iran[i]], color='red', $
      xra=[3000., 10000.]
    djs_oplot, lambda[0:nspec-1], data[0:nspec-1,iran[i]], $
      xra=[3000., 10000.]
endfor

end_print

end
