;+
; NAME:
;   k_derived
; PURPOSE:
;   create derived quantities file 
; CALLING SEQUENCE:
;   k_derived [, vname= ]
; OPTIONAL INPUTS:
;   vname - name of set of templates to use (default 'default')
; COMMENTS:
;   Create k_nmf_derived.vname.fits file, with the HDUs:
;       HDU0: [nbasis, nt] coefficients of each template
;       HDU1: [nspec, nt] spectrum of each template
;       HDU2: [nspec, nt] spectrum of each template (no lines)
;       HDU3: [nspec, nt] spectrum of each template (no dust)
;       HDU4: [nspec, nt] spectrum of each template (no lines or dust)
;       HDU5: [nspec, nt] spectrum of each template (unsmoothed)
;       HDU6: [nspec, nt] spectrum of each template (no lines; unsmoothed)
;       HDU7: [nspec, nt] spectrum of each template (no dust; unsmoothed)
;       HDU8: [nspec, nt] spectrum of each template (no lines or dust;
;                          unsmoothed)
;       HDU9: [nspec, nt] emission line spectrum of each template 
;       HDU10: [nspec, nt] dust extinction factor of each template
;       HDU11: [nspec] wavelengths for each spectrum
;       HDU12: [nage, nt] SFR for each template
;       HDU13: [nage, nt] metallicity for each template
;       HDU14: [nage] ages at which SFR given
;       HDU15: [nage] age differential (multiply SFR to get number of stars) 
;       HDU16: [nt] total mass formed in each template
;       HDU17: [nt] total current mass in each template
;       HDU18: [nt] metallicity of current stars in each template 
;       HDU19: [nt] total mass formed in each template in past 300 Myrs
;       HDU20: [nt] total mass formed in each template in past 1 Gyrs
;       HDU21: ages of each basis vector [Nages, Nmets, Ndust]
;       HDU22: metallicities of each basis vector [Nages, Nmets, Ndust]
;       HDU23: dust properties of each basis vector [Nages, Nmets, Ndust]
;       HDU24: fraction of original stellar mass surviving at age 
;              of this basis vector [Nages, Nmets, Ndust]
; REVISION HISTORY:
;   2005-08-15 MRB, NYU
;-
;------------------------------------------------------------------------------
pro k_derived, vname=vname

if (NOT keyword_set(vname)) then $
  vname='default'

metallicities=[0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05]

hdr=headfits(getenv('KCORRECT_DIR')+ $
             '/data/templates/k_nmf_mmatrix.'+vname+'.fits')
nspec=long(sxpar(hdr,'NSPEC'))
ndusts=long(sxpar(hdr, 'NDUST'))
nmets=long(sxpar(hdr, 'NMET'))
nages=long(sxpar(hdr, 'NAGE'))
nextra=long(sxpar(hdr,'NEXTRA'))
nel=long(sxpar(hdr,'NEL'))

spec=mrdfits(getenv('KCORRECT_DIR')+'/data/templates/k_nmf_mmatrix.'+ $
             vname+'.fits', 0, hdr)
spec=spec[0:nspec-1L, *]
lambda=mrdfits(getenv('KCORRECT_DIR')+ $
               '/data/templates/k_nmf_mmatrix.'+ $
               vname+'.fits', 1)
lambda=lambda[0:nspec-1L]
dust=mrdfits(getenv('KCORRECT_DIR')+ $
             '/data/templates/k_nmf_mmatrix.'+vname+'.fits',2)
met=mrdfits(getenv('KCORRECT_DIR')+ $
             '/data/templates/k_nmf_mmatrix.'+vname+'.fits',3)
met=metallicities[met]
mremain=mrdfits(getenv('KCORRECT_DIR')+ $
                '/data/templates/k_nmf_mmatrix.'+vname+'.fits',9)
age=mrdfits(getenv('KCORRECT_DIR')+ $
             '/data/templates/k_nmf_mmatrix.'+vname+'.fits',4)
rawspec=mrdfits(getenv('KCORRECT_DIR')+ $
                '/data/templates/k_nmf_rawspec.'+ $
                vname+'.fits', 0)
templates=mrdfits(getenv('KCORRECT_DIR')+ $
                  '/data/templates/k_nmf_soln.'+ $
                  vname+'.fits', 0)
if((size(templates))[0] eq 1) then $
  nt=1 $
else $
  nt=(size(templates, /dim))[1]
nl=n_elements(lambda)
nb=n_elements(spec)/nl
ns=n_elements(dust)

;; make regular spectrum
loglam=alog10(lambda[0:nspec-1L])
absrc=3.631*2.99792*1.e-2/lambda^2
for i=0L, nb-1L do spec[*,i]=spec[*,i]*absrc
tspec_v300=spec[0L:nspec-1L,*]#templates
tspec_v300_nl=spec[*,0:n_elements(dust)-1L]#templates[0:n_elements(dust)-1L,*]
tspec_v0_nl=rawspec[*,*]#templates[0:n_elements(dust)-1L,*]

;; get stellar pops spectrum
ii=where(dust.tauv eq 0, nii)
nd=ns/nii
templates_nd= $
  total(reform(templates[0:n_elements(dust)-1L,*], nii, nd,nt), 2)
tspec_v0_nl_nd=rawspec[*,0:nii-1L]#templates_nd
tspec_v300_nl_nd=spec[*,0:nii-1L]#templates_nd
extinction=tspec_v300_nl/tspec_v300_nl_nd

;; get emission line spectrum
if(nel gt 0) then begin
    lspec_v300=spec[*,n_elements(dust):n_elements(dust)+nel-1L]# $
      templates[n_elements(dust):n_elements(dust)+nel-1L, *]
endif else begin
    lspec_v300=fltarr(nspec,1)
endelse
tspec_v0=tspec_v0_nl+lspec_v300
tspec_v0_nd=tspec_v0_nl_nd+lspec_v300
tspec_v300_nd=tspec_v300_nl_nd+lspec_v300

ages=age[0:nages-1L]
dage=fltarr(nages)
dage[1:nages-2L]=0.5*(ages[2:nages-1L]-ages[0:nages-3L])
dage[0]=ages[1]-ages[0]
dage[nages-1L]=ages[nages-1L]-ages[nages-2L]
sfh_tot=fltarr(nages, nt)
sfh_met=fltarr(nages, nt)
for i=0L, nt-1L do begin
    currsfh=reform(templates[0:nages*nmets*ndusts-1L,i], nages, nmets, ndusts)
    dust=reform(dust, nages, nmets, ndusts)
    for j=0L, nages-1L do $
      sfh_tot[j,i]=total(currsfh[j,*,*])/dage[j]
    for j=0L, nages-1L do $
      sfh_met[j,i]=total(currsfh[j,*,*]*met[j,*,*])/dage[j]
endfor
sfh_met=sfh_met/sfh_tot

tmass=total(templates[0:nages*nmets*ndusts-1L,*], 1)
tmremain=fltarr(nt)
for i=0L, nt-1L do $
  tmremain[i]=total(templates_nd[0:nages*nmets-1L,i]*mremain)
i300=where(age lt 3.e+8, n300)
tmass300=total(templates[i300, *], 1)
i1000=where(age lt 1.e+9, n1000)
tmass1000=total(templates[i1000, *], 1)
tmetallicity=fltarr(nt)
for i=0L, nt-1L do $
  tmetallicity[i]=total(met[*,*,0]*templates_nd[0:nages*nmets-1L,i]*mremain)/ $
  tmremain[i]

outfile=getenv('KCORRECT_DIR')+'/data/templates/k_nmf_derived.'+vname+'.fits'
sxaddpar, hdr, 'NT', nt, 'number of templates'
mwrfits, float(templates), outfile, hdr, /create
mwrfits, float(tspec_v300), outfile
mwrfits, float(tspec_v300_nl), outfile
mwrfits, float(tspec_v300_nd), outfile
mwrfits, float(tspec_v300_nl_nd), outfile
mwrfits, float(tspec_v0), outfile
mwrfits, float(tspec_v0_nl), outfile
mwrfits, float(tspec_v0_nd), outfile
mwrfits, float(tspec_v0_nl_nd), outfile
mwrfits, float(lspec_v300), outfile
mwrfits, float(extinction), outfile
mwrfits, float(lambda), outfile
mwrfits, float(sfh_tot), outfile
mwrfits, float(sfh_met), outfile
mwrfits, float(ages), outfile
mwrfits, float(dage), outfile
mwrfits, float(tmass), outfile
mwrfits, float(tmremain), outfile
mwrfits, float(tmetallicity), outfile
mwrfits, float(tmass300), outfile
mwrfits, float(tmass1000), outfile
mwrfits, float(age), outfile
mwrfits, float(met), outfile
mwrfits, dust, outfile
mwrfits, float(mremain), outfile

end
