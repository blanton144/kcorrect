;+
; NAME:
;   k_plot_templates
; COMMENTS:
;   Creates:
;     templates_spec.ps - summary plot of the spectra of each template
;     templates_sfh.ps - summary plot of SFH of each template
; REVISION HISTORY:
;   30-Nov-2004  Michael Blanton (NYU)
;-
;------------------------------------------------------------------------------
pro k_plot_templates

;; read in the basis set info
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
nextra=long(sxpar(hdr, 'NEL'))
nzf=long(sxpar(hdr, 'NZ'))
nfilter=long(sxpar(hdr, 'NFILTER'))
ndusts=long(sxpar(hdr, 'NDUST'))
nmets=long(sxpar(hdr, 'NMET'))
nages=long(sxpar(hdr, 'NAGE'))
filternames=['F', 'N', 'u', 'g', 'r', 'i', 'z', 'J', 'H', 'K_s', 'B', 'R', 'I']
metallicities=[0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05]

;; read in the results
templates=mrdfits('k_nmf_soln.fits')
coeffs=mrdfits('k_nmf_soln.fits',1)
nbasis=n_elements(mmatrix)/n_elements(lambda)
nt=n_elements(templates)/nbasis

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

tspec=mmatrix#templates
;;dfact=3.826/(4.*!DPI)/(3.086)^2*1.e-5 
absrc=3.631*2.99792*1.e-2/lambda^2
for i=0L, nt-1L do $
  tspec[*,i]=tspec[*,i]*absrc

colors=['blue', 'green', 'magenta', 'red']

;; for each template, show the spectrum (full and just optical)
k_print, filename='templates_spec.ps', axis_char_scale=1.4
!P.MULTI=[0, 1, 2]
!Y.MARGIN=[1,0]
yra=[0.8, 1.2]*minmax(tspec[where(tspec gt 0.)])
xra=minmax(lambda)*[0.95,1.05]
djs_plot, [0.1], [0.1], /nodata, /xlog, /ylog, xra=xra, yra=yra, $
  ytitle='flux (erg/cm^{2}/s/\AA/M_{sun})',xch=0.001
for i=0L, nt-1L do $
  djs_oplot, lambda[0:nspec-1], tspec[0:nspec-1,i], /xlog, /ylog, th=3, $
  color=colors[i]
axis,!X.CRANGE[0],10.^(!Y.CRANGE[1]),xaxis=1, $
  xtitle=textoidl('!8\lambda !6(\AA)'), xcharsize=1.4, /ylog
xra=[3001.,9999.]
il=where(lambda[0:nspec-1] gt 3000. and lambda[0:nspec-1] lt 10000., nl)
yra=[0.8, 1.2]*minmax(tspec[il,*])
djs_plot, [1], [1], /nodata, /ylog, xra=xra, yra=yra, $
  xtitle='!8\lambda !6(\AA)', $
  ytitle='flux (erg/cm^{2}/s/\AA/M_{sun})'
ilab=where(lambda[0:nspec-2] lt 9200. and lambda[1:nspec-1] gt 9200.)
for i=0L, nt-1L do $
  djs_oplot, lambda[il], tspec[il, i], /ylog, th=3, color=colors[i] 
for i=0L, nt-1L do $
  djs_xyouts, lambda[ilab], tspec[ilab,i]*1.2, strtrim(string(i),2), $
  align=1., color=colors[i], charsize=1.5
k_end_print

sfh=reform(templates[0:nages*nmets*ndusts-1L,*], nages, nmets, ndusts,nt)
dust=reform(dust, nages, nmets, ndusts)
sfh_tot=fltarr(nages,nt)
for i=0L, nt-1L do $
  for j=0L, nages-1L do sfh_tot[j,i]=total(sfh[j,*,*,i])
sfh_met=fltarr(nages,nt)
for i=0L, nt-1L do $
  for j=0L, nages-1L do sfh_met[j,i]=total(sfh[j,*,*,i]* $
                                           metallicities[met[j,*,*]])
sfh_met=sfh_met/sfh_tot
sfh_dust=fltarr(nages,nt)
for i=0L, nt-1L do $
  for j=0L, nages-1L do sfh_dust[j,i]=total(sfh[j,*,*,i]* $
                                            dust[j,*,*].tauv)
sfh_dust=sfh_dust/sfh_tot

;; for each template, the star-formation history and (mean metallicity
;; and dustiness of same) 
k_print, filename='templates_sfh.ps', axis_char_scale=1.4
dage=fltarr(nages)
dage[1:nages-2L]=0.5*(age[2:nages-1L]-age[0:nages-3L])
dage[0]=age[1]-age[0]
dage[nages-1L]=age[nages-1L]-age[nages-2L]
!P.MULTI=[0, 1, 3]
!P.CHARSIZE=1.9
!X.CHARSIZE=1.9
!Y.CHARSIZE=1.9
yra=max(sfh_tot[where(sfh_tot gt 0)])*[1.e-5,1.5]
xra=[6.e+5, 2.2e+10]
djs_plot, [1], [1], ytitle='SFR (M_{sun}/yr', $
  /xlog, /ylog, /nodata, yra=yra, xra=xra, xch=0.001
for i=0L, nt-1L do $
  djs_oplot, age[*,0,0], sfh_tot[*,i], /xlog, /ylog, color=colors[i], $
  th=4
for i=0L, nt-1L do $
  djs_oplot, age[*,0,0], sfh_tot[*,i], psym=4, /xlog, /ylog, color=colors[i], $
  th=4
yra=[-1.e-2, 0.08]
djs_plot, [1], [1], ytitle='metallicity', $
  /xlog, /nodata, yra=yra, xra=xra, xch=0.001
for i=0L, nt-1L do $
  djs_oplot, age[*,0,0], sfh_met[*,i], /xlog, /ylog, color=colors[i], th=4
for i=0L, nt-1L do $
  djs_oplot, age[*,0,0], sfh_met[*,i], psym=4, /xlog, /ylog, color=colors[i], $
  th=4
yra=[-0.2, 6.1]
djs_plot, [1], [1], ytitle='avg \tau_V', $
  /xlog, /nodata, yra=yra, xra=xra
for i=0L, nt-1L do $
  djs_oplot, age[*,0,0], sfh_dust[*,i], /xlog, color=colors[i], th=4
for i=0L, nt-1L do $
  djs_oplot, age[*,0,0], sfh_dust[*,i], /xlog, color=colors[i], th=4, psym=4
k_end_print 

end
