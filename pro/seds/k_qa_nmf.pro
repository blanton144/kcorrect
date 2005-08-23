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
nsubsp=300

hogg_usersym, 10, /fill

;; read in the template basics
mmatrix=mrdfits('k_nmf_mmatrix.fits',0,hdr)
lambda=mrdfits('k_nmf_mmatrix.fits',1)
dust=mrdfits('k_nmf_mmatrix.fits',2)
met=mrdfits('k_nmf_mmatrix.fits',3)
age=mrdfits('k_nmf_mmatrix.fits',4)
rawspec=mrdfits('k_nmf_rawspec.fits',0)
filterlist=string(mrdfits('k_nmf_mmatrix.fits',5))
emmatrix=mrdfits('k_nmf_early.fits',0,earlyhdr)
lmmatrix=mrdfits('k_nmf_late.fits',0,latehdr)
zf=mrdfits('k_nmf_mmatrix.fits',6)
nspec=long(sxpar(hdr, 'NSPEC'))
back=long(sxpar(hdr, 'BACK'))
nextra=long(sxpar(hdr, 'NEXTRA'))
nzf=long(sxpar(hdr, 'NZ'))
nfilter=long(sxpar(hdr, 'NFILTER'))
ndusts=long(sxpar(hdr, 'NDUST'))
nmets=long(sxpar(hdr, 'NMET'))
nages=long(sxpar(hdr, 'NAGE'))
filternames=['F', 'N', 'u', 'g', 'r', 'i', 'z', 'J', 'H', 'K_s', 'B', 'R', $
             'I', 'J', "H", 'K_s', 'B', 'V', 'i', 'z', '[3.6]', '[4.5]', $
             '[5.8]', '[8.0]', '[24]' ]
metallicities=[0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05]

;; read in the data
datastr=mrdfits('k_nmf_spdata.fits',1)
vals=mrdfits('k_nmf_spdata.fits',2)
ivar=mrdfits('k_nmf_spdata.fits',3)
xx=mrdfits('k_nmf_spdata.fits',4)
data=create_struct(datastr, $
                   'val', fltarr(n_elements(vals)), $
                   'x', fltarr(n_elements(vals)))
data.val=vals
data.x=xx
data_ivar=create_struct(datastr, $
                        'val', fltarr(n_elements(vals)), $
                        'x', fltarr(n_elements(vals)))
data_ivar.val=ivar
data_ivar.x=xx
zhelio=mrdfits('k_nmf_spdata.fits',7)
iz=long(floor((nzf-1.)*(zhelio-zf[0])/(zf[nzf-1]-zf[0])+0.5))

;; read in the results
templates=mrdfits('k_nmf_soln.fits')
coeffs=mrdfits('k_nmf_soln.fits',1)
nt=n_elements(coeffs)/n_elements(zhelio)
model=data
early=data
late=data
mcoeffs=templates#coeffs
mmeval, model, transpose(mmatrix), mcoeffs

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

;; make vmatrix and lambda for early
outevmatrix=emmatrix[0:nspec-1L,*]#templates
absrc=3.631*2.99792*1.e-2/(lambda[0:nspec-1L])^2
for i=0L, nt-1L do $
  outevmatrix[*,i]=outevmatrix[*,i]*absrc
outlambda=fltarr(nspec+1L)
dlg10l=alog10(lambda[1])-alog10(lambda[0])
outlambda[0:nspec-1L]= 10.^(alog10(lambda[0:nspec-1L])-0.5*dlg10l)
outlambda[nspec]= 10.^(alog10(lambda[nspec-1L])+0.5*dlg10l)
k_write_ascii_table,outevmatrix,'vmatrix.'+version+'early.dat'
k_write_ascii_table,outlambda,'lambda.'+version+'early.dat'

;outlvmatrix=lmmatrix[0:nspec-1L,*]#templates
;absrc=3.631*2.99792*1.e-2/(lambda[0:nspec-1L])^2
;for i=0L, nt-1L do $
;  outlvmatrix[*,i]=outlvmatrix[*,i]*absrc
;outlambda=fltarr(nspec+1L)
;dlg10l=alog10(lambda[1])-alog10(lambda[0])
;outlambda[0:nspec-1L]= 10.^(alog10(lambda[0:nspec-1L])-0.5*dlg10l)
;outlambda[nspec]= 10.^(alog10(lambda[nspec-1L])+0.5*dlg10l)
;k_write_ascii_table,outlvmatrix,'vmatrix.'+version+'late.dat'
;k_write_ascii_table,outlambda,'lambda.'+version+'late.dat'

set_print, filename='k_qa_nmf.ps'

;; for each template, show the spectrum (full and just optical)
!P.MULTI=[0, 1, 2]
tspec=outvmatrix
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

;; derived quantities
stmass=transpose(t_mass#coeffs)
age=((t_mass*t_age)#coeffs)/stmass
metallicity=((t_mass*t_metallicity)#coeffs)/stmass
!P.MULTI=[0,1,2]
djs_plot, stmass, metallicity, psym=4, xtitle='stellar mass (M_\odot)', $
  ytitle='metallicity',/xlog, xcharsize=0.0001
djs_plot, stmass, age, psym=4, xtitle='stellar mass (M_\odot)', $
  ytitle='mean stellar age', /xlog, /ylog
djs_plot, zhelio, metallicity, psym=4, xtitle='redshift' , $
  ytitle='metallicity', xcharsize=0.0001
djs_plot, zhelio, age, psym=4, xtitle='redshift' , $
  ytitle='mean stellar age', /ylog

;; color and residuals vs. redshift
!P.MULTI=[0,1,2]
ngals=n_elements(data.rowstart)
col=fltarr(nfilter-1L,ngals)
colerr=fltarr(nfilter-1L,ngals)
mcol=fltarr(nfilter-1L,ngals)
for ifilter=0L, nfilter-2L do begin
    f0start=nspec+ifilter*nzf
    f0end=nspec+(ifilter+1L)*nzf-1L
    f1start=nspec+(ifilter+1L)*nzf
    f1end=nspec+(ifilter+2L)*nzf-1L
    for j=0L, ngals-1L do begin
        currx=data.rowstart[j]+lindgen(data.nxrow[j])
        i0=where(data.x[currx] ge f0start AND data.x[currx] le f0end, n0)
        i1=where(data.x[currx] ge f1start AND data.x[currx] le f1end, n1)
        if(n1 gt 0 and n0 gt 0) then begin
            mcol[ifilter,j]=-2.5*alog10(model.val[currx[i0[0]]]/ $
                             model.val[currx[i1[0]]])
            col[ifilter,j]=-2.5*alog10(data.val[currx[i0[0]]]/ $
                            data.val[currx[i1[0]]])
            colerr[ifilter,j]= $
              2.5/alog(10.)*sqrt(1./data_ivar.val[currx[i0[0]]]/ $
                                 data.val[currx[i0[0]]]^2+ $
                                 1./data_ivar.val[currx[i1[0]]]/ $
                                 data.val[currx[i1[0]]]^2)
        endif
    endfor
    igood=where(colerr[ifilter,*] ne 0., ngood)

    if(ngood gt 0) then begin
        djs_plot, zhelio[igood], col[ifilter,igood], psym=8, color='red', $
         symsize=0.2
        djs_oplot, zhelio[igood], mcol[ifilter,igood], psym=8, $
          ytitle=filternames[ifilter]+'-'+filternames[ifilter+1], symsize=0.2
        ;;djs_oploterr, zhelio[igood], col[ifilter,igood], $
        ;;  yerr=colerr[ifilter,igood], color='red', symsize=0.2
        yra=minmax(col[ifilter,igood]-mcol[ifilter,igood]) 
        yra[0]=yra[0] > (-0.9)
        yra[1]=yra[1] < (0.9)
        hogg_scatterplot, zhelio[igood], $
	        col[ifilter,igood]-mcol[ifilter,igood], yra=yra, $
          ytitle=filternames[ifilter]+'-'+filternames[ifilter+1]+' residual', $
          xtitle='redshift z', /cond, xnpix=24, ynpix=20, exp=0.5, $
	satfrac=0.001
        djs_oplot, minmax(zhelio), [0., 0.], color='blue', linestyle=1
    endif
endfor


;; colors vs. color
!P.MULTI=[0,1,1]
for ifilter=0L, nfilter-3L do begin
    igood=where(colerr[ifilter,*] gt 0. and $
                colerr[ifilter+1L,*] gt 0., ngood)
    if(ngood gt 0) then begin
        djs_plot, col[ifilter,igood], col[ifilter+1L,igood], psym=8, $
          color='red', symsize=0.2, $
          xtitle=filternames[ifilter]+'-'+filternames[ifilter+1], $
          ytitle=filternames[ifilter+1]+'-'+filternames[ifilter+2]
        djs_oplot, mcol[ifilter,igood], mcol[ifilter+1,igood], psym=8, $
	symsize=0.2, $
          xtitle=filternames[ifilter]+'-'+filternames[ifilter+1], $
          ytitle=filternames[ifilter+1]+'-'+filternames[ifilter+2]
        ;djs_oploterr, col[ifilter,igood], col[ifilter+1,igood], $
        ;  xerr=colerr[ifilter,igood], yerr=colerr[ifilter+1,igood], $
        ;  color='red'
    endif
endfor

;; random set of spectra
nsubsp=nsubsp < n_elements(zhelio)
iran=shuffle_indx(n_elements(zhelio), num_sub=nsubsp)
for i=0L, nsubsp-1L do begin
    icurr=data.rowstart[iran[i]]+lindgen(data.nxrow[iran[i]])
    xcurr=data.x[icurr]
    yra=minmax(data.val[icurr])*[0.7,1.3]
    djs_plot, lambda[xcurr], model.val[icurr], color='red', $
      xra=[1000., 10000.], /xlog, yra=yra
    djs_oplot, lambda[xcurr], data.val[icurr], psym=8, /xlog, symsize=0.2
    djs_oplot, lambda[xcurr], data.val[icurr], /xlog
endfor

end_print

end
