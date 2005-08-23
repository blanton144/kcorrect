;+
; NAME:
;   k_fit_templates
; PURPOSE:
;   run the template fitting code, with priming
; CALLING SEQUENCE:
;   k_fit_templates [, nt= ]
; OPTIONAL INPUTS:
;   nt - maximum number of templates to fit for (default 4)
; REVISION HISTORY:
;   09-Apr-2005  Michael Blanton (NYU)
;-
;------------------------------------------------------------------------------
pro k_fit_templates, nt=nt, nprime=nprime

if(NOT keyword_set(nt)) then nt=4

if(0) then begin
spawn, 'mkdir -p photo1'
cd, 'photo1'
k_nmf_mmatrix, /nodust, /noel
k_nmf_spdata, nsdss_spec=0L, nlrg_spec=0L, /few, flux='model', sample='dr4'
k_run_nmf, nt=1L, niter=1000L, /reset, /qa
cd, '../'

spawn, 'mkdir -p photo2'
cd, 'photo2'
spawn, 'cp ../photo1/k_nmf_mmatrix.fits .'
spawn, 'cp ../photo1/k_nmf_early.fits .'
spawn, 'cp ../photo1/k_nmf_spdata.fits .'
spawn, 'cp ../photo1/k_nmf_rawspec.fits .'
spawn, 'cp ../photo1/k_nmf_soln.fits .'
templates1=mrdfits('k_nmf_soln.fits',0)
coeffs1=mrdfits('k_nmf_soln.fits',1)
nb=(size(templates1,/dim))[0]
ng=(size(coeffs1,/dim))[1]
templates=fltarr(nb,2)
templates[*,0]=templates1*(0.95+0.1*randomu(seed, nb))+0.01/float(nb)
templates[*,1]=templates1*(0.95+0.1*randomu(seed, nb))+0.01/float(nb)
mwrfits, templates, 'k_nmf_soln.fits', /create
k_run_nmf, nt=2L, niter=10000L, /qa
cd, '../'

spawn, 'mkdir -p photo3'
cd, 'photo3'
spawn, 'cp ../photo2/k_nmf_mmatrix.fits .'
spawn, 'cp ../photo2/k_nmf_early.fits .'
spawn, 'cp ../photo2/k_nmf_spdata.fits .'
spawn, 'cp ../photo2/k_nmf_rawspec.fits .'
spawn, 'cp ../photo2/k_nmf_soln.fits .'
templates2=mrdfits('k_nmf_soln.fits',0)
coeffs2=mrdfits('k_nmf_soln.fits',1)
nb=(size(templates2,/dim))[0]
ng=(size(coeffs2,/dim))[1]
templates=fltarr(nb,3)
templates[*,0]=templates2[*,0]*(0.95+0.1*randomu(seed, nb))+0.01/float(nb)
templates[*,1]=templates2[*,1]*(0.95+0.1*randomu(seed, nb))+0.01/float(nb)
templates[*,2]=0.5*(templates2[*,0]+templates2[*,1])* $
  (0.5+1.0*randomu(seed, nb))+0.01/float(nb)
mwrfits, templates, 'k_nmf_soln.fits', /create
k_run_nmf, nt=3L, niter=20000L, /qa
cd, '../'

spawn, 'mkdir -p photo4'
cd, 'photo4'
spawn, 'cp ../photo3/k_nmf_mmatrix.fits .'
spawn, 'cp ../photo3/k_nmf_early.fits .'
spawn, 'cp ../photo3/k_nmf_spdata.fits .'
spawn, 'cp ../photo3/k_nmf_rawspec.fits .'
spawn, 'cp ../photo3/k_nmf_soln.fits .'
templates3=mrdfits('k_nmf_soln.fits',0)
coeffs3=mrdfits('k_nmf_soln.fits',1)
nb=(size(templates3,/dim))[0]
ng=(size(coeffs3,/dim))[1]
templates=fltarr(nb,4)
templates[*,0]=templates3[*,0]*(0.95+0.1*randomu(seed, nb))+0.01/float(nb)
templates[*,1]=templates3[*,1]*(0.95+0.1*randomu(seed, nb))+0.01/float(nb)
templates[*,2]=templates3[*,2]*(0.95+0.1*randomu(seed, nb))+0.01/float(nb)
templates[*,3]=0.333*(templates3[*,0]+templates3[*,1]+templates3[*,2])* $
  (0.1+1.8*randomu(seed, nb))+0.1/float(nb)
mwrfits, templates, 'k_nmf_soln.fits', /create
k_run_nmf, nt=4L, niter=40000L, /qa
cd, '../'

spawn, 'mkdir -p photo2m'
cd, 'photo2m'
spawn, 'cp ../photo2/k_nmf_mmatrix.fits .'
spawn, 'cp ../photo2/k_nmf_early.fits .'
spawn, 'cp ../photo2/k_nmf_rawspec.fits .'
spawn, 'cp ../photo2/k_nmf_soln.fits .'
k_nmf_spdata, nsdss_spec=0L, nlrg_spec=0L, flux='model', sample='dr4'
k_run_nmf, nt=2L, niter=5000L, /qa
templates2=mrdfits('k_nmf_soln.fits',0)
coeffs2=mrdfits('k_nmf_soln.fits',1)
nb=(size(templates2,/dim))[0]
ng=(size(coeffs2,/dim))[1]
templates=fltarr(nb,2)
templates[*,0]=templates2[*,0]*(0.99+0.02*randomu(seed, nb))+0.001/float(nb)
templates[*,1]=templates2[*,1]*(0.99+0.02*randomu(seed, nb))+0.001/float(nb)
mwrfits, templates, 'k_nmf_soln.fits', /create
k_run_nmf, nt=2L, niter=10000L, /qa
cd, '../'

spawn, 'mkdir -p photo3m'
cd, 'photo3m'
spawn, 'cp ../photo3/k_nmf_mmatrix.fits .'
spawn, 'cp ../photo3/k_nmf_early.fits .'
spawn, 'cp ../photo3/k_nmf_rawspec.fits .'
spawn, 'cp ../photo3/k_nmf_soln.fits .'
spawn, 'cp ../photo2m/k_nmf_spdata.fits .'
k_run_nmf, nt=3L, niter=5000L, /qa
templates3=mrdfits('k_nmf_soln.fits',0)
coeffs3=mrdfits('k_nmf_soln.fits',1)
nb=(size(templates3,/dim))[0]
ng=(size(coeffs3,/dim))[1]
templates=fltarr(nb,3)
templates[*,0]=templates3[*,0]*(0.99+0.02*randomu(seed, nb))+0.001/float(nb)
templates[*,1]=templates3[*,1]*(0.99+0.02*randomu(seed, nb))+0.001/float(nb)
templates[*,2]=templates3[*,2]*(0.99+0.02*randomu(seed, nb))+0.001/float(nb)
mwrfits, templates, 'k_nmf_soln.fits', /create
k_run_nmf, nt=3L, niter=40000L, /qa
cd, '../'

spawn, 'mkdir -p photo4m'
cd, 'photo4m'
spawn, 'cp ../photo4/k_nmf_mmatrix.fits .'
spawn, 'cp ../photo4/k_nmf_early.fits .'
spawn, 'cp ../photo4/k_nmf_rawspec.fits .'
spawn, 'cp ../photo4/k_nmf_soln.fits .'
spawn, 'cp ../photo2m/k_nmf_spdata.fits .'
k_run_nmf, nt=4L, niter=5000L, /qa
templates4=mrdfits('k_nmf_soln.fits',0)
coeffs4=mrdfits('k_nmf_soln.fits',1)
nb=(size(templates4,/dim))[0]
ng=(size(coeffs4,/dim))[1]
templates=fltarr(nb,4)
templates[*,0]=templates4[*,0]*(0.99+0.02*randomu(seed, nb))+0.001/float(nb)
templates[*,1]=templates4[*,1]*(0.99+0.02*randomu(seed, nb))+0.001/float(nb)
templates[*,2]=templates4[*,2]*(0.99+0.02*randomu(seed, nb))+0.001/float(nb)
templates[*,3]=templates4[*,3]*(0.99+0.02*randomu(seed, nb))+0.001/float(nb)
mwrfits, templates, 'k_nmf_soln.fits', /create
k_run_nmf, nt=4L, niter=40000L, /qa
cd, '../'

spawn, 'mkdir -p dust2m'
cd, 'dust2m'
k_nmf_mmatrix, /noel
hdr=headfits('k_nmf_mmatrix.fits')
ndusts=long(sxpar(hdr, 'NDUST'))
nmets=long(sxpar(hdr, 'NMET'))
nages=long(sxpar(hdr, 'NAGE'))
spawn, 'cp ../photo2m/k_nmf_soln.fits .'
spawn, 'cp ../photo2m/k_nmf_spdata.fits .'
templates2=mrdfits('k_nmf_soln.fits',0)
coeffs2=mrdfits('k_nmf_soln.fits',1)
nb=(size(templates2,/dim))[0]
ng=(size(coeffs2,/dim))[1]
templates=fltarr(nb*ndusts,2)
scale=fltarr(ndusts)+0.01
scale[0]=1.
for j=0L, 2L-1L do $
  for i=0L, ndusts-1L do $
  templates[i*nmets*nages:(i+1)*nmets*nages-1L,j]= $
  scale[i]*(templates2[*,j]*(1.00+0.00*randomu(seed, nb))+0.0001/float(nb))
mwrfits, templates, 'k_nmf_soln.fits', /create
mwrfits, coeffs2, 'k_nmf_soln.fits'
k_run_nmf, nt=2L, niter=10000L, /qa
cd, '../'

spawn, 'mkdir -p dust4m'
cd, 'dust4m'
spawn, 'cp ../dust2m/k_nmf_mmatrix.fits .'
spawn, 'cp ../dust2m/k_nmf_early.fits .'
hdr=headfits('k_nmf_mmatrix.fits')
ndusts=long(sxpar(hdr, 'NDUST'))
nmets=long(sxpar(hdr, 'NMET'))
nages=long(sxpar(hdr, 'NAGE'))
spawn, 'cp ../photo4m/k_nmf_soln.fits .'
spawn, 'cp ../photo4m/k_nmf_spdata.fits .'
templates4=mrdfits('k_nmf_soln.fits',0)
coeffs4=mrdfits('k_nmf_soln.fits',1)
nb=(size(templates4,/dim))[0]
ng=(size(coeffs4,/dim))[1]
templates=fltarr(nb*ndusts,4)
scale=fltarr(ndusts)+0.01
scale[0]=1.
for j=0L, 4L-1L do $
  for i=0L, ndusts-1L do $
  templates[i*nmets*nages:(i+1)*nmets*nages-1L,j]= $
  scale[i]*(templates4[*,j]*(1.00+0.00*randomu(seed, nb))+0.0001/float(nb))
mwrfits, templates, 'k_nmf_soln.fits', /create
mwrfits, coeffs4, 'k_nmf_soln.fits'
k_run_nmf, nt=4L, niter=20000L, /qa
cd, '../'


spawn, 'mkdir -p dust5m'
cd, 'dust5m'
spawn, 'cp ../dust4m/k_nmf_mmatrix.fits .'
spawn, 'cp ../dust4m/k_nmf_early.fits .'
spawn, 'cp ../dust4m/k_nmf_rawspec.fits .'
hdr=headfits('k_nmf_mmatrix.fits')
ndusts=long(sxpar(hdr, 'NDUST'))
nmets=long(sxpar(hdr, 'NMET'))
nages=long(sxpar(hdr, 'NAGE'))
spawn, 'cp ../dust4m/k_nmf_spdata.fits .'
spawn, 'cp ../dust4m/k_nmf_soln.fits .'
templates4=mrdfits('k_nmf_soln.fits',0)
coeffs4=mrdfits('k_nmf_soln.fits',1)
nb=(size(templates4,/dim))[0]
ng=(size(coeffs4,/dim))[1]
templates=fltarr(nb,5)
coeffs=fltarr(5,ng)
templates[*,0]=templates4[*,0]*(0.97+0.06*randomu(seed, nb))+0.0001/float(nb)
templates[*,1]=templates4[*,1]*(0.97+0.06*randomu(seed, nb))+0.0001/float(nb)
templates[*,2]=templates4[*,2]*(0.97+0.06*randomu(seed, nb))+0.0001/float(nb)
templates[*,3]=templates4[*,3]*(0.97+0.06*randomu(seed, nb))+0.0001/float(nb)
templates[*,4]=(total(templates4, 2))/4.*(0.95+0.1*randomu(seed, nb))+ $
  0.0001/float(nb)
coeffs[0:3, *]=coeffs4
coeffs[4,*]=0.005*total(coeffs4,1)
mwrfits, templates, 'k_nmf_soln.fits', /create
mwrfits, coeffs, 'k_nmf_soln.fits'
for j=0, 9 do $
k_run_nmf, nt=5L, niter=10000L, /qa
cd, '../'

cd, 'dust5m'
for j=0, 9 do $
k_run_nmf, nt=5L, niter=10000L, /qa
cd, '../'

if(0) then begin
spawn, 'mkdir -p spec4m'
cd, 'spec4m'
k_nmf_mmatrix
k_nmf_spdata, /spchop, flux='model', sample='dr4'
data=mrdfits('k_nmf_spdata.fits',1)
ngals=n_elements(data.rowstart)
hdr=headfits('k_nmf_mmatrix.fits')
nextra=long(sxpar(hdr, 'NEXTRA'))
spawn, 'cp ../dust4m/k_nmf_soln.fits .'
templates4=mrdfits('k_nmf_soln.fits',0)
coeffs4=mrdfits('k_nmf_soln.fits',1)
nb=(size(templates4,/dim))[0]
ng4=(size(coeffs4,/dim))[1]
templates=fltarr(nb+nextra,4)
templates[0:nb-1, 0:3]= templates4
templates[nb:nb+nextra-1, 0:3]=0.1/float(nb)
coeffs=fltarr(4,ngals)
coeffs[0:3,0:ng4-1L]=coeffs4
coeffs[0:3,ng4:ngals-1L]=0.3
mwrfits, templates, 'k_nmf_soln.fits', /create
mwrfits, coeffs4, 'k_nmf_soln.fits'
k_run_nmf, nt=4L, niter=50000L, /qa
cd, '../'
endif
endif

spawn, 'mkdir -p spec5m'
cd, 'spec5m'
if(0) then begin
k_nmf_mmatrix
k_nmf_spdata, /spchop, flux='petro', sample='dr4'
data=mrdfits('k_nmf_spdata.fits',1)
ngals=n_elements(data.rowstart)
hdr=headfits('k_nmf_mmatrix.fits')
nextra=long(sxpar(hdr, 'NEXTRA'))
spawn, 'cp ../dust5m/k_nmf_soln.fits .'
templates5=mrdfits('k_nmf_soln.fits',0)
coeffs5=mrdfits('k_nmf_soln.fits',1)
nb=(size(templates5,/dim))[0]
ng5=(size(coeffs5,/dim))[1]
templates=fltarr(nb+nextra,5)
templates[0:nb-1, 0:4]= templates5+1.e-5
templates[nb:nb+nextra-1, 0:4]=0.01/float(nb)
coeffs=fltarr(5,ngals)
coeffs[0:4,0:ng5-1L]=coeffs5+1.e-5
coeffs[0:4,ng5:ngals-1L]=0.2
mwrfits, templates, 'k_nmf_soln.fits', /create
mwrfits, coeffs, 'k_nmf_soln.fits'
for i=0L, 99L do $
k_run_nmf, nt=5L, niter=500L, /qa
endif
cd, '../'

spawn, 'mkdir -p spec5mb'
cd, 'spec5mb'
spawn, 'cp ../spec5m/k_nmf_mmatrix.fits .'
spawn, 'cp ../spec5m/k_nmf_early.fits .'
spawn, 'cp ../spec5m/k_nmf_late.fits .'
spawn, 'cp ../spec5m/k_nmf_rawspec.fits .'
spawn, 'cp ../spec5m/k_nmf_soln.fits .'
k_nmf_spdata, flux='petro', sample='dr4', seed=500
for i=0L, 999L do $
k_run_nmf, nt=5L, niter=500L, /qa
cd, '../'

end
