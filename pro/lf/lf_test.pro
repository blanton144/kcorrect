;+
; NAME:
;   lf_test
; PURPOSE:
;   test luminosity function fit
; REVISION HISTORY:
;   2002-5-23  written - Blanton
;-
function lf_test_kcorrect,z
return,2.5*alog10(1.+z)-2.5*alog10(1.1)
end
;
pro lf_test,rhozfit=rhozfit

; settings
nrandom=500000L
omega0=0.3
omegal0=0.7
dh=2.99792D+5/100.
sample_absmmin=-23.5
sample_absmmax=-15.5
sigabsmag=0.22
mkspacing=1.*sigabsmag
sample_zmin=0.02
sample_zmax=0.20
nmvals=500
nzvals=150
errmean=0.01
errspread=1.0
q=2.0
p=0.0
zzero=0.1
buffer=2.
qaschechter={schstr, phistar:1., mstar:-20.5, alpha:-1.}

; some other simple calculations
zvals=sample_zmin+(sample_zmax-sample_zmin)* $
  (dindgen(nzvals)+0.5)/double(nzvals)
vvals=dh^3*comdis(zvals,omega0,omegal0)^3/3.
sample_vmin=dh^3*comdis(sample_zmin,omega0,omegal0)^3/3.
sample_vmax=dh^3*comdis(sample_zmax,omega0,omegal0)^3/3.
volume=sample_vmax-sample_vmin

if(keyword_set(rhozfit)) then begin
    nrhozfunc=26

    openr,unit,getenv('BUILDDIR')+'/data/modes/modes4mike.dat',/get_lun
    modesarr=dblarr(28,179)
    readf,unit,modesarr
    close,unit
    free_lun,unit
    rmodes=modesarr[0,*]
    nbarmodes=modesarr[1,*]*sqrt(160000.)
    psimodes=modesarr[2:27,*]
    rhozfuncmodes=dblarr(26,n_elements(nbarmodes))
    for i=0, 25 do $
      rhozfuncmodes[i,*]=psimodes[i,*]/nbarmodes
    rvals=dh*comdis(zvals,omega0,omegal0)
    rhozfunc=dblarr(nzvals,nrhozfunc)
    for i=0, nrhozfunc-1 do $
        rhozfunc[*,i]=interpol(rhozfuncmodes[i,*],rmodes,rvals)
    modesarr=0L
        
    openr,unit,getenv('BUILDDIR')+'/data/modes/covar4mike.dat',/get_lun
    covararr=dblarr(27,27)
    readf,unit,covararr
    close,unit
    free_lun,unit
    rhozcovar=covararr[1:nrhozfunc,1:nrhozfunc]
    covararr=0L
endif

; find distribution of actual absolute mags
splog,'absolute mags'
absmag=sample_absmmin-0.5*buffer+randomu(seed,nrandom)* $
  (sample_absmmax-sample_absmmin+buffer)
schechter_test=lf_schechter(absmag,1.,qaschechter.mstar,qaschechter.alpha)/ $
  max(lf_schechter(absmag,1.,qaschechter.mstar,qaschechter.alpha))
chances=randomu(seed,nrandom)
indx=where(chances lt schechter_test)
absmag=absmag[indx]
nrandom=n_elements(absmag)

; what should phistar be? (remember to account for num den evolution)
splog,'get phistar'
amvals=sample_absmmin-0.5*buffer+(sample_absmmax-sample_absmmin+buffer)* $
  (dindgen(nmvals)+0.5)/double(nmvals)
schvals=lf_schechter(amvals,qaschechter.phistar,qaschechter.mstar, $
                     qaschechter.alpha)
totalsch=total(schvals*(amvals[1]-amvals[0]),/double)
numden=double(nrandom)/volume
qaschechter.phistar=numden/totalsch*10.^(-0.4*(sample_zmax-zzero))

; assign redshifts
splog,'redshifts'
vv=sample_vmin+randomu(seed,nrandom)*(sample_vmax-sample_vmin)
zz=dblarr(nrandom)
for i=0L, nrandom-1L do begin
    indx=where(vvals lt vv[i])
    j=indx[n_elements(indx)-1]
    if(j eq -1) then j=0
    if(j eq nzvals-1) then j=nzvals-2
    jp1=j+1
    sv=(vv[i]-vvals[j])/(vvals[jp1]-vvals[j])
    zz[i]=zvals[j]+sv*(zvals[jp1]-zvals[j])
endfor
zzchances=10.^(0.4*p*(zz-zzero))
zzchances=zzchances/max(zzchances)
indx=where(randomu(seed,n_elements(zzchances)) lt zzchances)
absmag=absmag[indx]
zz=zz[indx]
vv=vv[indx]

; add errors and evolution to absmag distribution
splog,'add errors'
absmerr=errmean*exp(errspread*randomn(seed,nrandom))
absmag=absmag+absmerr*randomn(seed,nrandom)-q*(zz-zzero)

; figure out distance modulus stuff
splog,'distance modulus'
dc=dh*comdis2(zz,omega0,omegal0)
dl=(1.+zz)*dc
dm=5.*alog10(dl)+25.
kcorrect=lf_test_kcorrect(zz)

; apply flux limits
splog,'flux limits'
mmin=13.5
mmax=17.5+0.*randomu(seed,nrandom)
absmmin=mmin-dm-kcorrect
absmmax=mmax-dm-kcorrect
indx=where(absmmin lt sample_absmmin,count)
if(count gt 0) then absmmin[indx]=sample_absmmin
indx=where(absmmin gt sample_absmmax,count)
if(count gt 0) then absmmin[indx]=sample_absmmax
indx=where(absmmax lt sample_absmmin,count)
if(count gt 0) then absmmax[indx]=sample_absmmin
indx=where(absmmax gt sample_absmmax,count)
if(count gt 0) then absmmax[indx]=sample_absmmax
indx=where(absmag gt absmmin and absmag lt absmmax)
ngals=n_elements(indx)
absmag=absmag[indx]
absmmin=absmmin[indx]
absmmax=absmmax[indx]
mmin=mmin[indx]
mmax=mmax[indx]
absmerr=absmerr[indx]
kcorrect=kcorrect[indx]
zz=zz[indx]

; put together k-correction info
splog,'kcorrections'
uniqkcorrect=dblarr(nzvals,1)
uniqkcorrect[*,0]=lf_test_kcorrect(zvals)
ikcorrect=lonarr(ngals)

; call npgauss
splog,'npgauss'
plot,zz,absmag,psym=3
lf_npgauss,zz,absmag,absmerr,kcorrect,ikcorrect,uniqkcorrect, $
  nzvals,mmin,mmax,sigabsmag,mkspacing,sample_absmmin,sample_absmmax, $
  sample_zmin,sample_zmax,absmk,phi,q,p,rhozpars,omega0=omega0, $
  omegal0=omegal0, iabsmerr=iabsmerr,immin=immin,immax=immax, $
  uniqabsmerr=uniqabsmerr,uniqmmin=uniqmmin,uniqmmax=uniqmmax, $
  qaplot='tmplf.ps',qaschechter=qaschechter,maxiter=10, $
  rhozfunc=rhozfunc,rhozcovar=rhozcovar,fixp=0.
save,filename='test.sav' 

splog,'selfunc'
sf=lf_selfunc_npgauss(zz,kcorrect,mmin,mmax,sample_absmmin,sample_absmmax, $
                      phi,absmk,sigabsmag,q,p,omega0=omega0,omegal0=omegal0)

splog,'mean_density'
lf_mean_density,zz,sf,[1.],iabsmerr,immin,immax,ikcorrect,uniqabsmerr, $
  uniqmmin,uniqmmax,uniqkcorrect,sample_absmmin,sample_absmmax,sample_zmin, $
  sample_zmax,phi,absmk,sigabsmag,q,p,n1,n1err,nzvals=nzvals,omega0=omega0, $
  omegal0=omegal0,/simple

lf_plotps,absmk,n1*phi,sigabsmag,sample_absmmin,sample_absmmax, $
  filename='lf_test.ps',schechter=qaschechter

save,filename='test.sav' 

end
