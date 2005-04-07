;+
; NAME:
;   lf_npgauss
; PURPOSE:
;   calculate luminosity function in presence of uncertainties and
;   errors, plus evolution and large scale structure
; USAGE:
;   lf_npgauss,zz,absmag,absmerr,kcorrect,ikcorrect,uniqkcorrect, $
;              nzvals,mmin,mmax,sigabsmag,mkspacing, $
;              sample_absmmin,sample_absmmax,sample_zmin,sample_zmax, $
;              absmk,phi,q,p [,rhozpars, omega0=, omegal0=, zzero=, $
;              limit_precision=, error_precision=, error_min=, $
;              error_max=, iabsmerr=, immin=, immax=, uniqabsmerr=, $
;              uniqmmin=, uniqmmax=, maxiter=, fixp=, fixq=, $
;              nborder=, qaplot=, qaschechter=, rhozfunc=, rhozcovar=, $
;              lssfix=, lnlike=, init_power= ]
; INPUTS:
;   zz               [N] redshifts
;   absmag           [N] absolute magnitudes
;   absmerr          [N] absolute magnitude errors
;   kcorrect         [N] K-correction for each galaxy
;   ikcorrect        [N] index of K-correction function for each gal
;   uniqkcorrect     [nzvals,N_k] K(z) for each K-correction function
;   nzvals           number of redshift slices in each K-correction fn
;   mmin             [N] minimum apparent mag for each galaxy
;   mmax             [N] maximum apparent mag for each galaxy
;   sigabsmag        width of each gaussian in fit
;   mkspacing        distance btw each gaussian in fit (approximate)
;   sample_absmmin   observed absolute magnitude minimum of sample
;   sample_absmmax   observed absolute magnitude maximum of sample
;   sample_zmin      minimum redshift of sample
;   sample_zmax      maximum redshift of sample
; OPTIONAL KEYWORDS:
;   /lssfix          fit for rhoz functions
; OPTIONAL INPUTS:
;   omega0           omega_matter to use (default: 0.3)
;   omegal0          omega_lambda to use (default: 0.7)
;   zzero            redshift to center fit on (default: 0.1)
;   limit_precision  limit t
;   rhozfunc, rhozcovar - radial density functions
; OUTPUTS:
;   absmk            peak of each gaussian
;   phi              amplitude of each gaussian
;   q                evolution of lum
;   p                evolution of num
;   rhozpars         parameters multiplying the radial density functions
;   lnlike           output likelihood
; REVISION HISTORY:
;   2002-5-22  written - Blanton
;-
;
; Calculate likelihood function and its derivs
function lf_npgauss_func, params, dparams

; common block 
common com_lf_npgauss, lf_absmk, lf_zdiff, lf_absmag, $
  lf_nphi, lf_sigabsmag, lf_dvdz, lf_nzvals, lf_zvals, $
  lf_zzero, lf_uniqkcorrect, lf_ikcorrect, lf_omega0, lf_omegal0, $
  lf_sample_absmmin, lf_sample_absmmax, lf_immin, lf_immax, $
  lf_iabsmerr, lf_uniqabsmerr, lf_uniqmmax, lf_uniqmmin, $
  lf_nborder, lf_qaplot, lf_qaschechter, lf_rhozfunc, lf_rhozcovar, $
  lf_nrhozfunc, lf_zz, lf_sample_zmin, lf_sample_zmax

; settings and constants
ngals=n_elements(lf_zdiff)
dh=2.99792D+5/100.
pi=3.14159265358979D+0
dz=lf_zvals[1]-lf_zvals[0]

; convert parameters to named parameters for convenience
; note: params is [ln(phi_k), q, p, [rhozpar] ]
lf_paramstophi,params,phi
q=params[lf_nphi-lf_nborder+0L]
p=params[lf_nphi-lf_nborder+1L]
if(lf_nrhozfunc gt 0) then $
  rhozpars=params[lf_nphi-lf_nborder+2:lf_nphi-lf_nborder+2+lf_nrhozfunc-1]

; set gaussian variances for each object
var=lf_uniqabsmerr[lf_iabsmerr]^2+lf_sigabsmag^2

;;;;;;;;;;;;;;;;;;;;;
; set aik and derivs
;

; calculate rhoz for each galaxy if desired
rhoz=dblarr(n_elements(lf_dvdz))+1.
if(lf_nrhozfunc gt 0) then begin
    rhoz_int=1.+lf_rhozfunc#rhozpars 
    indx=where(rhoz_int le 0.,count)
    if(count gt 0) then rhoz_int[indx]=1.e-3
    rhoz=interpolate(rhoz_int,(double(lf_nzvals))*(lf_zz-lf_sample_zmin)/ $
                     (lf_sample_zmax-lf_sample_zmin)-0.5)
    indx=where(rhoz le 0.,count)
    if(count gt 0) then rhoz[indx]=1.e-3
    rhozfunc=dblarr(lf_nrhozfunc,n_elements(rhoz))
    for i=0, lf_nrhozfunc-1 do begin
        rhozfunc[i,*]=1.+interpolate(lf_rhozfunc[*,i],double(lf_nzvals)* $
                                     (lf_zz-lf_sample_zmin)/ $
                                     (lf_sample_zmax-lf_sample_zmin))
    endfor
endif else begin
    rhoz_int=1.
    rhoz=1.
endelse

; set aik 
splog,'doing aik'
invvar=1./var
invvarfactor=1./sqrt(2.*pi*var)
aik=dblarr(ngals,lf_nphi)
pdiff=invvarfactor*lf_dvdz*10.^(0.4*p*lf_zdiff)
for k=0L, lf_nphi-1L do $
  aik[*,k]=pdiff*exp(-0.5*(lf_absmag-lf_absmk[k]+q*lf_zdiff)^2*invvar)

; set derivs of aik
daikdp=dblarr(ngals,lf_nphi)
daikdq=dblarr(ngals,lf_nphi)
splog,'doing aik derivs'
for k=0L, lf_nphi-1L do begin
    daikdq[*,k]=-(lf_absmag-lf_absmk[k]+q*lf_zdiff)*lf_zdiff*aik[*,k]*invvar
    daikdp[*,k]=0.4*alog(10.)*lf_zdiff*aik[*,k]
endfor

; calculate aik sums
splog,'doing aik sums'
fullphi=transpose(phi#replicate(1.,ngals))
sumaik=total(fullphi*aik,2)
sumdaikdq=total(fullphi*daikdq,2)
sumdaikdp=total(fullphi*daikdp,2)

;;;;;;;;;;;;;;;;;;;;;;;;;;;
; set bik and derivs
;

; construct dm for the lf_zvals
splog,'doing dm'
dcvals=dh*comdis(lf_zvals,lf_omega0,lf_omegal0) 
dlvals=(1.+lf_zvals)*dcvals
dmvals=5.*alog10(dlvals)+25.

; construct volume integral values
splog,'doing vol integral'
dvdzvals=dh^3*dcomvoldz(lf_zvals,lf_omega0,lf_omegal0)
int_dvdz_p=10.^(0.4*p*(lf_zvals-lf_zzero))*dvdzvals*dz

; construct interpolation arrays to speed up erfs 
splog,'doing erf arrays'
xerfarr=-2.5+(dindgen(300)+0.5)*5./300.
yerfarr=errorf(xerfarr)
yerfarr[0]=yerfarr[1]
yerfarr[298]=yerfarr[299]

; make the table of possible bik values to pick from
; based on the mmin, mmax, errors, and kcorrection for each
; object
splog,'make bik table'

; how big is table
nmmin=n_elements(lf_uniqmmin)
nmmax=n_elements(lf_uniqmmax)
nabsmerr=n_elements(lf_uniqabsmerr)
nkcorrect=n_elements(lf_uniqkcorrect)/lf_nzvals
help,nabsmerr,nmmin,nmmax,nkcorrect,lf_nphi

; allocate table (and sum tables)
biktable=dblarr(lf_nphi,nabsmerr,nmmin,nmmax,nkcorrect)
;dbikdqtable=dblarr(lf_nphi,nabsmerr,nmmin,nmmax,nkcorrect)
;dbikdptable=dblarr(lf_nphi,nabsmerr,nmmin,nmmax,nkcorrect)
sumbiktable=dblarr(nabsmerr,nmmin,nmmax,nkcorrect)
sumdbikdqtable=dblarr(nabsmerr,nmmin,nmmax,nkcorrect)
sumdbikdptable=dblarr(nabsmerr,nmmin,nmmax,nkcorrect)
if(lf_nrhozfunc gt 0) then begin
    sumdbikdrhozpartable=dblarr(lf_nrhozfunc,nabsmerr,nmmin,nmmax, $
                                nkcorrect)
endif

; now make it
for s=0L, nabsmerr-1L do begin

; calculate factors in front for this uncertainty
    currvarfactor=1./sqrt(2.*pi*(lf_uniqabsmerr[s]^2+lf_sigabsmag^2))
    svfactor=sqrt(0.5/(lf_uniqabsmerr[s]^2+lf_sigabsmag^2))

; now look at mmin,mmax, and kcorrect
    for immin=0L, nmmin-1L do begin
        for immax=0L, nmmax-1L do begin
            for ikcorrect=0L, nkcorrect-1L do begin
            
; construct zabsmmax and zabsmmin
                zabsmmax=lf_uniqmmax[immax]-lf_uniqkcorrect[*,ikcorrect]-dmvals
                zabsmmin=lf_uniqmmin[immin]-lf_uniqkcorrect[*,ikcorrect]-dmvals
                
; check that they are not out of bounds
                indx=where(zabsmmax gt lf_sample_absmmax,count)
                if(count gt 0) then zabsmmax[indx]=lf_sample_absmmax
                indx=where(zabsmmax lt lf_sample_absmmin,count)
                if(count gt 0) then zabsmmax[indx]=lf_sample_absmmin
                indx=where(zabsmmin lt lf_sample_absmmin,count)
                if(count gt 0) then zabsmmin[indx]=lf_sample_absmmin
                indx=where(zabsmmin gt lf_sample_absmmax,count)
                if(count gt 0) then zabsmmin[indx]=lf_sample_absmmax

; check that they are greater than zero
                indx=where(zabsmmin gt zabsmmax,count)
                if(count gt 0) then zabsmmin=zabsmmax

; calculate bik etc for this case -- see lf_npgauss.tex for 
; mathematical derivation of this
                dbikdrhoz_inttot=dblarr(lf_nzvals)
                dbikdp_inttot=dblarr(lf_nzvals)
                dbikdq_inttot=dblarr(lf_nzvals)
                for k=0L, lf_nphi-1L do begin
                    xmax=svfactor*(zabsmmax-lf_absmk[k]+q*(lf_zvals-lf_zzero))
                    xmin=svfactor*(zabsmmin-lf_absmk[k]+q*(lf_zvals-lf_zzero))
                    bik_int=interpolate(yerfarr,60.*(xmax+2.5))- $
                      interpolate(yerfarr,60.*(xmin+2.5))
                    ;dbikdp_int=bik_int*(lf_zvals-lf_zzero)
                    ;dbikdq_int=(lf_zvals-lf_zzero)*(exp(-xmax^2)-exp(-xmin^2))
                    ;dbikdptable[k,s,immin,immax,ikcorrect]=0.5*0.4*alog(10.)* $
                      ;total(dbikdp_int*int_dvdz_p*rhoz_int,/double)
                    ;dbikdqtable[k,s,immin,immax,ikcorrect]= $
                      ;currvarfactor*total(dbikdq_int*int_dvdz_p*rhoz_int, $
                                          ;/double)
                    dbikdq_inttot=dbikdq_inttot+ $
                      (exp(-xmax^2)-exp(-xmin^2))*phi[k]
                    dbikdp_inttot=dbikdp_inttot+bik_int*phi[k]
                    if(lf_nrhozfunc gt 0) then $
                      dbikdrhoz_inttot=dbikdrhoz_inttot+ $
                      0.5*bik_int*int_dvdz_p*phi[k]
                    biktable[k,s,immin,immax,ikcorrect]=0.5* $
                      total(bik_int*int_dvdz_p*rhoz_int)
                endfor

; tables of sums
                sumbiktable[s,immin,immax,ikcorrect]= $
                  total(biktable[*,s,immin,immax,ikcorrect]*phi,/double)
                ;sumdbikdptable[s,immin,immax,ikcorrect]= $
                  ;total(dbikdptable[*,s,immin,immax,ikcorrect]*phi,/double)
                ;sumdbikdqtable[s,immin,immax,ikcorrect]= $
                  ;total(dbikdqtable[*,s,immin,immax,ikcorrect]*phi,/double)
                sumdbikdptable[s,immin,immax,ikcorrect]=0.5*0.4*alog(10.)* $
                  total((lf_zvals-lf_zzero)*int_dvdz_p*dbikdp_inttot* $
                        rhoz_int,/double)
                sumdbikdqtable[s,immin,immax,ikcorrect]= $
                  currvarfactor*total((lf_zvals-lf_zzero)*int_dvdz_p* $
                                      dbikdq_inttot*rhoz_int,/double)
                if(lf_nrhozfunc gt 0) then begin
                    for irhoz=0L, lf_nrhozfunc-1L do begin
                        sumdbikdrhozpartable[irhoz,s,immin,immax,ikcorrect]= $
                          total(dbikdrhoz_inttot*lf_rhozfunc[*,irhoz], $
                                /double)
                    endfor
                endif
            endfor
        endfor
    endfor
endfor

; now use the table to make bik and its sums
splog,'make bik and derivs'
bik=dblarr(ngals,lf_nphi)
for k=0L, lf_nphi-1L do begin
    bik[*,k]=biktable[replicate(k,ngals),lf_iabsmerr,lf_immin, $
                      lf_immax,lf_ikcorrect]
endfor
sumbik=sumbiktable[lf_iabsmerr,lf_immin,lf_immax,lf_ikcorrect]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; now construct E and derivs 
;
splog,'constructing e and de'

; construct e
e=-total(alog(sumaik*rhoz/sumbik),/double)
if(lf_nrhozfunc gt 0) then begin
    chi2=rhozpars##transpose(invert(lf_rhozcovar)#rhozpars)
    e=e+0.5*chi2
endif

; construct de/dparams 
splog,'constructing de'
dparams=dblarr(lf_nphi-lf_nborder+2+lf_nrhozfunc)

; de/d(lnphi)
for k=0L, lf_nphi-lf_nborder-1L do $
  dparams[k]=phi[k]*total(bik[*,k]/sumbik-aik[*,k]/sumaik)

; contribution of border to de/d(lnphi)
for i=0L, lf_nborder-1L do begin
    dparams[lf_nphi-lf_nborder-1]=dparams[lf_nphi-lf_nborder-1]+ $
      phi[lf_nphi-lf_nborder-1]*total(bik[*,lf_nphi-lf_nborder+i]/sumbik- $
                                      aik[*,lf_nphi-lf_nborder+i]/sumaik)
endfor

; de/dq
dparams[lf_nphi-lf_nborder+0]= $
  total(sumdbikdqtable[lf_iabsmerr,lf_immin,lf_immax,lf_ikcorrect]/ $
        sumbik-sumdaikdq/sumaik)

; de/dp
dparams[lf_nphi-lf_nborder+1]= $
  total(sumdbikdptable[lf_iabsmerr,lf_immin,lf_immax,lf_ikcorrect]/ $
        sumbik-sumdaikdp/sumaik)

; de/da_i
if(lf_nrhozfunc gt 0) then begin
    for j=0L, lf_nrhozfunc-1L do begin
        tmpsum=sumdbikdrhozpartable[j,*,*,*,*]
        dparams[j+lf_nphi-lf_nborder+2]= $
          total(tmpsum[lf_iabsmerr,lf_immin,lf_immax,lf_ikcorrect]/sumbik- $
                (rhozfunc[j,*]-1.)/rhoz,/double)
    endfor

; a_i prior
    dparams[lf_nphi-lf_nborder+2:lf_nphi-lf_nborder+2+lf_nrhozfunc-1L]= $
      dparams[lf_nphi-lf_nborder+2:lf_nphi-lf_nborder+2+lf_nrhozfunc-1L]+ $
      invert(lf_rhozcovar)#rhozpars
endif

;splog,params
;splog,dparams
splog,e
;if(e[0] ne e[0]) then stop
;if(dparams[0] ne dparams[0]) then stop
return,e

end
;
; enforce normalization and make qa plot --- normalization falls
; out of likelihood calculation, so renormalizing should not affect
; maximization 
pro lf_npgauss_fixsum, myfunct, params, iter, fnorm, functargs=fcnargs, $
                       parinfo=parinfo,quiet=quiet, _extra=extra

common com_lf_npgauss

if(n_elements(fnorm) eq 0) then fnorm=lf_npgauss_func(params,dparams)

; generate full phi, normalize it, and tranlate it back to p
lf_paramstophi,params,phi
q=params[lf_nphi-lf_nborder+0L]
p=params[lf_nphi-lf_nborder+1L]
if(lf_nrhozfunc gt 0) then $
  rhozpars=params[lf_nphi-lf_nborder+2:lf_nphi-lf_nborder+2+lf_nrhozfunc-1]
lf_normalize,lf_absmk,phi,lf_sample_absmmin,lf_sample_absmmax,lf_sigabsmag
lf_phitoparams,phi,params

; write it
splog,iter
splog,'ln(phi) = '
splog,string(params[0:lf_nphi-lf_nborder-1])
splog,'rhozpars = '
if(lf_nrhozfunc gt 0) then $
  splog,string(params[lf_nphi-lf_nborder+2-1: $
                      lf_nphi-lf_nborder+2+lf_nrhozfunc-1])
splog,'q = '+string(params[lf_nphi-lf_nborder+0])
splog,'p = '+string(params[lf_nphi-lf_nborder+1])
splog,string(fnorm)

; make the plot 
if(keyword_set(lf_qaplot)) then $
  lf_plotps,lf_absmk,phi,lf_sigabsmag,lf_sample_absmmin,lf_sample_absmmax, $
  schechter=lf_qaschechter,filename=lf_qaplot

if(lf_nrhozfunc gt 0) then begin
    rhoz_int=1.+lf_rhozfunc#rhozpars 
    set_print,filename='tmp_rhozfunc.ps'
    djs_plot,lf_zvals,rhoz_int,thick=3
    end_print
endif

end
;
; Translate p to phi
pro lf_paramstophi,params,phi

common com_lf_npgauss

phi=dblarr(lf_nphi)  
phi[0:lf_nphi-lf_nborder-1]=exp(params[0L:lf_nphi-lf_nborder-1L])
if(lf_nborder gt 0) then $
  phi[lf_nphi-lf_nborder+lindgen(lf_nborder)]=phi[lf_nphi-lf_nborder-1]

end
;
; Translate phi to p
pro lf_phitoparams,phi,params

common com_lf_npgauss

if(n_elements(params) eq 0) then $
  params=dblarr(lf_nphi-lf_nborder+2+lf_nrhozfunc)
params[0:lf_nphi-lf_nborder-1]=alog(phi[0L:lf_nphi-lf_nborder-1L])

end
;
; minimization wrapper
pro lf_npgauss,zz,absmag,absmerr,kcorrect,ikcorrect,uniqkcorrect, $
               nzvals,mmin,mmax,sigabsmag,mkspacing, $
               sample_absmmin,sample_absmmax,sample_zmin,sample_zmax, $
               absmk,phi,q,p,rhozpars,omega0=omega0,omegal0=omegal0, $
               zzero=zzero,limit_precision=limit_precision, $
               error_precision=error_precision,error_min=error_min, $
               error_max=error_max,iabsmerr=iabsmerr,immin=immin, $
               immax=immax,uniqabsmerr=uniqabsmerr,uniqmmin=uniqmmin, $
               uniqmmax=uniqmmax,maxiter=maxiter,fixp=fixp,fixq=fixq, $
               nborder=nborder,qaplot=qaplot,qaschechter=qaschechter, $
               rhozfunc=rhozfunc, rhozcovar=rhozcovar, lssfix=lssfix, $
               lnlike=lnlike, init_power=init_power

common com_lf_npgauss

; set defaults
if(n_elements(omega0) eq 0) then omega0=0.3
if(n_elements(omegal0) eq 0) then omegal0=0.7
if(n_elements(zzero) eq 0) then zzero=0.1
if(NOT keyword_set(limit_precision)) then limit_precision=0.01
if(NOT keyword_set(error_precision)) then error_precision=0.02
if(NOT keyword_set(error_min)) then error_min=0.02
if(NOT keyword_set(error_max)) then error_max=1.0
if(NOT keyword_set(init_power)) then init_power=0.0
if(n_elements(maxiter) eq 0) then maxiter=50
if(n_elements(nborder) eq 0) then nborder=4

; settings
pi=3.14159265358979D
dh=2.99792D+5/100.D

; set common block variables
nphi=long((sample_absmmax-sample_absmmin)/mkspacing)+nborder
lf_nborder=nborder
lf_absmk=dblarr(nphi)
lf_absmk[0:nphi-nborder-1]=sample_absmmin+(sample_absmmax-sample_absmmin)* $
  (dindgen(nphi-nborder)+0.5)/double(nphi-nborder)
if(nborder gt 0) then $
  lf_absmk[nphi-nborder+lindgen(nborder)]=lf_absmk[nphi-nborder-1]+ $
  (dindgen(nborder)+1.)*(lf_absmk[1]-lf_absmk[0])
lf_zzero=zzero
lf_zz=zz
lf_zdiff=zz-lf_zzero
lf_absmag=absmag
lf_nphi=nphi
lf_sigabsmag=sigabsmag
lf_sample_absmmin=sample_absmmin
lf_sample_absmmax=sample_absmmax
lf_sample_zmin=sample_zmin
lf_sample_zmax=sample_zmax
lf_omega0=omega0
lf_omegal0=omegal0
lf_nzvals=nzvals
lf_uniqkcorrect=uniqkcorrect
lf_ikcorrect=ikcorrect
if(keyword_set(qaschechter)) then lf_qaschechter=qaschechter
if(keyword_set(qaplot)) then lf_qaplot=qaplot
lf_nrhozfunc=0L
if(keyword_set(rhozfunc)) then begin
    lf_rhozfunc=rhozfunc
    lf_nrhozfunc=n_elements(lf_rhozfunc)/nzvals
    if(keyword_set(rhozcovar)) then lf_rhozcovar=rhozcovar
endif

; construct zvals 
lf_zvals=sample_zmin+(sample_zmax-sample_zmin)* $
  (dindgen(lf_nzvals)+0.5)/double(lf_nzvals)

; find number of galaxies
ngals=n_elements(lf_absmag)

; calculate dvdz and the distance modulus for each object
dc=dh*comdis2(zz,lf_omega0,lf_omegal0) 
lf_dvdz=dh^3*dcomvoldz(zz,lf_omega0,lf_omegal0)
dl=(1.+zz)*dc
lf_dm=5.*alog10(dl)+25.

; find the set of limits, to the precision specificied; warn if there 
; galaxies outside of the limits, and set them to the limit
immin=long(mmin/limit_precision+0.5d)
sortimmin=immin[sort(immin)]
uniqimmin=sortimmin[uniq(sortimmin)]
lf_uniqmmin=double(uniqimmin)*limit_precision
lf_immin=lonarr(ngals)
for i=0L, n_elements(uniqimmin)-1 do begin
    indx=where(immin eq uniqimmin[i],count)
    if(count eq 0) then begin
        splog,'inconsistency in lf_npgauss (1)'
        return
    endif
    lf_immin[indx]=i
endfor
indx=where(lf_absmag+lf_dm+kcorrect lt lf_uniqmmin[lf_immin],count)
if(count gt 0) then begin
    splog,strtrim(string(count))+' galaxies outside of their mmin limit'
    lf_absmag[indx]=lf_uniqmmin[lf_immin[indx]]-lf_dm[indx]-kcorrect[indx]
endif
immin=lf_immin
uniqmmin=lf_uniqmmin

immax=long(mmax/limit_precision+0.5d)
sortimmax=immax[sort(immax)]
uniqimmax=sortimmax[uniq(sortimmax)]
lf_uniqmmax=double(uniqimmax)*limit_precision
lf_immax=lonarr(ngals)
for i=0L, n_elements(uniqimmax)-1 do begin
    indx=where(immax eq uniqimmax[i],count)
    if(count eq 0) then begin
        splog,'inconsistency in lf_npgauss (1)'
        return
    endif
    lf_immax[indx]=i
endfor
indx=where(lf_absmag+lf_dm+kcorrect gt lf_uniqmmax[lf_immax],count)
if(count gt 0) then begin
    splog,strtrim(string(count))+' galaxies outside of their mmax limit'
    lf_absmag[indx]=lf_uniqmmax[lf_immax[indx]]-lf_dm[indx]-kcorrect[indx]
endif 
immax=lf_immax
uniqmmax=lf_uniqmmax

; now find the set of uncertainties
tmpabsmerr=absmerr
indx=where(abs(tmpabsmerr) gt error_max,count)
if(count gt 0) then tmpabsmerr[indx]=error_max
indx=where(abs(tmpabsmerr) lt error_min,count)
if(count gt 0) then tmpabsmerr[indx]=error_min
tmpiabsmerr=long((tmpabsmerr-error_min)/error_precision)
absmerrvals=error_min+dindgen(max(tmpiabsmerr)+1L)*error_precision
sorttmpiabsmerr=tmpiabsmerr[sort(tmpiabsmerr)]
uniqiabsmerr=sorttmpiabsmerr[uniq(sorttmpiabsmerr)]
lf_uniqabsmerr=absmerrvals[uniqiabsmerr]
lf_iabsmerr=lonarr(ngals)
for i=0L, n_elements(uniqiabsmerr)-1 do begin
    indx=where(tmpiabsmerr eq uniqiabsmerr[i],count)
    if(count eq 0) then begin
        splog,'inconsistency in lf_npgauss (1)'
        return
    endif
    lf_iabsmerr[indx]=i
endfor
iabsmerr=lf_iabsmerr
uniqabsmerr=lf_uniqabsmerr

; set initial conditions, unless they have been fully specified
if(n_elements(q) eq 0 OR n_elements(p) eq 0 $ 
   OR (n_elements(phi) ne lf_nphi)) then begin
    phi=dblarr(lf_nphi)+10.^(init_power*(lf_absmk-min(lf_absmk)))
    q=0.
    p=0.
    if(n_elements(fixq) gt 0) then q=fixq
    if(n_elements(fixp) gt 0) then p=fixp
endif
lf_normalize,lf_absmk,phi,lf_sample_absmmin,lf_sample_absmmax,lf_sigabsmag
lf_phitoparams,phi,start
start[lf_nphi-lf_nborder+0]=q
start[lf_nphi-lf_nborder+1]=p
if(lf_nrhozfunc gt 0) then begin
    if(n_elements(rhozpars) gt 0) then $
      start[lf_nphi-lf_nborder+2:lf_nphi-lf_nborder+2+lf_nrhozfunc-1]= $
      rhozpars $
    else $
      start[lf_nphi-lf_nborder+2:lf_nphi-lf_nborder+2+lf_nrhozfunc-1]=0.
endif

; set up minimization
parinfo=replicate({value:0., $
                   fixed:0, $
                   limited:[0,0], $
                   limits:[0.,0], $
                   tied:''}, lf_nphi-lf_nborder+2+lf_nrhozfunc)
parinfo[lf_nphi-lf_nborder+0].limited=1
parinfo[lf_nphi-lf_nborder+1].limited=1
parinfo[lf_nphi-lf_nborder+0].limits=[-8.,8.]
parinfo[lf_nphi-lf_nborder+1].limits=[-8.,8.]
parinfo[lf_nphi-lf_nborder+0].fixed=n_elements(fixq)
parinfo[lf_nphi-lf_nborder+1].fixed=n_elements(fixp)
if(keyword_set(lssfix) AND lf_nrhozfunc gt 0) then begin
    parinfo[lf_nphi-lf_nborder+2:lf_nphi-lf_nborder+2+lf_nrhozfunc-1].fixed=1 
endif

; run minimization
if(maxiter gt 0) then begin
    params=tnmin('lf_npgauss_func',start,maxiter=maxiter,autoderivative=0, $
                 parinfo=parinfo,iterproc='lf_npgauss_fixsum')
endif else begin
    params=start
endelse

; last output and normalization
lf_npgauss_fixsum, 'lf_npgauss_func', params, 'last', fnorm

; set output values
lnlike=fnorm
absmk=lf_absmk
lf_paramstophi,params,phi
q=params[lf_nphi-lf_nborder+0]
p=params[lf_nphi-lf_nborder+1]
if(lf_nrhozfunc gt 0) then $
  rhozpars=params[lf_nphi-lf_nborder+2:lf_nphi-lf_nborder+2+lf_nrhozfunc-1]

end
