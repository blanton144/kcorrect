;+
; NAME:
;   k_fitnn_coeff
;
; PURPOSE:
;   Given a set of maggies, redshifts, and spectra to fit to, find the 
;   nonnegative combination of the inputs which produces the outputs.
;
; CALLING SEQUENCE:
;   k_fit_nn, galaxy_maggies, galaxy_invvar, galaxy_z, coeff, $
;      [ematrix=, zvals=, filterlist=, bmatrix=, lambda=, rmatrix=]
;
; INPUTS:
;   galaxy_maggies   - maggies in each band for each galaxy [N_band, N_gal]
;   galaxy_invvar - errors in each band for each galaxy [N_band, N_gal]
;   galaxy_z      - redshift for each galaxy [N_gal]
;   ematrix       - eigentemplates [N_dim, N_template]
;   filterlist    - list of files with filter information [N_band]
;   bmatrix       - orthogonal templates spanning SED space [N_lambda, N_dim]
;   lambda        - wavelengths for orthogonal templates [N_lambda]
;   rmatrix       - look up table for bmatrix and filter information 
;                   [N_z, N_dim, N_band]
;   zvals         - look up table for rmatrix [N_z]
;
; OPTIONAL INPUTS:
;
; KEYWORDS:
;   covar - galaxy_invvar is actually a set of covariance matrices
;
; OUTPUTS:
;   coeff    - output coefficients [N_template, N_gal]
;
; OPTIONAL INPUT/OUTPUTS:
;
; COMMENTS:
;   galaxy_maggies is in maggies (f=10.^{-0.4*mag}). galaxy_invvar is in 
;   maggies^{-2} ((f*0.4*ln(10)*magerr)^{-2}).
;
;   If filterlist, bmatrix, and lambda are specified, rmatrix is created 
;   (and replaced if asked by the user)
;
; EXAMPLES:
;   To get the coefficients using the standard templates:
; 
;   IDL> k_fit_nn,galaxy_maggies,galaxy_invvar,galaxy_z, coeff
; 
;   Then you can pass "coeff" into k_reconstruct_nn
;
; BUGS:
;   Will fail if N_dim is unity.
;
; PROCEDURES CALLED:
;   k_load_ascii_table
;   k_create_r
;   Dynamic link to idl_k_fit_nn.c in libkcorrect.so
;
; REVISION HISTORY:
;   04-Jancb2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_fitnn_coeffs, galaxy_maggies, galaxy_invvar, galaxy_z, coeff, $
                    ematrix=ematrix, bmatrix=bmatrix, lambda=lambda, $
                    zvals=zvals, filterlist=filterlist, rmatrix=rmatrix, $
                    version=version, vpath=vpath, filterpath=filterpath, $
                    covar=covar, qarun=qarun, qastop=qastop

; Need at least 6 parameters
if (N_params() LT 4) then begin
    print, 'Syntax - k_fitnn_coeff, galaxy_maggies, galaxy_invvar, galaxy_z, coeff, $'  
    print, '        [ematrix=, bmatrix=, lambda=, zvals=, filterlist=, rmatrix=, $'
    print, '         version=, vpath=, filterpath=, /covar]'
    return
endif

if(NOT keyword_set(filterpath)) then $
  filterpath=getenv('KCORRECT_DIR')+'/data/filters'

; Get bmatrix and stuff from files if necessary
if(keyword_set(version)) then begin 
    if(NOT keyword_set(vpath)) then $
      vpath=getenv('KCORRECT_DIR')+'/data/etemplates'
    k_load_ascii_table,ematrix,vpath+'/ematrix.'+version+'.dat'
    k_load_ascii_table,bmatrix,vpath+'/bmatrix.'+version+'.dat'
    k_load_ascii_table,lambda,vpath+'/lambda.'+version+'.dat'
    filtfile=vpath+'/filterlist.'+version+'.dat'
    spawn,'cat '+filtfile+' | wc -l',nfilters
    nk=long(nfilters[0])-1l
    filterlist=strarr(nk)
    openr,unit,vpath+'/filterlist.'+version+'.dat',/get_lun
    readf,unit,nk
    readf,unit,filterlist
    close,unit
    free_lun,unit
endif

; Set source object name
soname=filepath('libkcorrect.so', root_dir=getenv('KCORRECT_DIR'), $
                subdirectory='lib')

; Determine dimensions
ngalaxy=long(n_elements(galaxy_z))
nk=long(n_elements(galaxy_maggies)/ngalaxy)
if (keyword_set(bmatrix) AND keyword_set(filterlist)AND keyword_set(lambda)) $
  then begin
    k_create_r,rmatrix,bmatrix,lambda,zvals,filterlist,filterpath=filterpath
endif else begin
    klog,'Using precomputed rmatrix.'
    if (NOT keyword_set(rmatrix)) then begin
        klog, 'need to specify rmatrix or bmatrix and filterlist'
    endif
endelse 
nz=long(n_elements(zvals))
nb=long(n_elements(rmatrix)/(nz*nk)) 	 
nt=long(n_elements(ematrix))/nb

; set inverse covariance matrices
galaxy_invcovar=dblarr(nk,nk,ngalaxy)
if(NOT keyword_set(covar)) then begin 
    for k=0, nk-1 do $
      galaxy_invcovar[k,k,*]=galaxy_invvar[k,*]
endif else begin
    for i=0L, ngalaxy-1L do $
      galaxy_invcovar[*,*,i]=invert(galaxy_invvar[*,*,i])
endelse

; for each galaxy, fit coeffs
nconstraints=nb
neconstraints=0L
nconstraints_max=nconstraints+1L
nvar=nt
nvar_max=nt
amatrix=dblarr(nconstraints_max,nvar)
amatrix[0:nb-1,*]=ematrix
bbmatrix=dblarr(nconstraints_max)
xl=replicate(-1.d+30,nvar)
xu=replicate(1.d+30,nvar)
lagrange=dblarr(nconstraints+nvar*2)
ifail=-1L
iprint=1L
coeff=dblarr(nvar,ngalaxy)
qascale=3.631*2.99792e-2
for i=0L, ngalaxy-1L do begin

; first, construct rjk(z).C^{-1}kk'.rj'k'(z)
    rlocal=dblarr(nk,nb)
    zpos=double(nz)*(galaxy_z[i]-zvals[0])/(zvals[nz-1]-zvals[0])
    for k=0L, nk-1L do $
      for j=0L, nb-1L do begin
        rlocal[k,j]=interpolate(rmatrix[*,j,k],zpos)
    endfor
    ermatrix=ematrix##rlocal
    cmatrix=ermatrix##galaxy_invcovar[*,*,i]##transpose(ermatrix)       

; next, construct rjk(z).C^{-1}kk'.lk'
    dmatrix=-ermatrix##galaxy_invcovar[*,*,i]##galaxy_maggies[*,i]
    
; fit using qld routine
    tmpcoeff=dblarr(nvar)
    retval=call_external(soname, 'idl_k_qld', long(nconstraints), $
                         long(neconstraints), long(nconstraints_max), $
                         long(nvar), long(nvar_max), double(cmatrix), $
                         double(dmatrix), double(amatrix), double(bbmatrix), $
                         double(xl), double(xu), tmpcoeff, lagrange, ifail, $
                         iprint)
    coeff[*,i]=tmpcoeff

    if(keyword_set(qarun)) then begin
        bandctrs=[3220.,4340.,5660.,6940.,8300.]*1.1
        bandvals=qascale*galaxy_maggies[*,i]/bandctrs^2
        spec=(bmatrix#ematrix)#tmpcoeff
        plot,lambda*(1.+galaxy_z[i]),spec/(1.+galaxy_z[i]), $
		          xra=[2000.,10000.],/ylog
        oplot,bandctrs,bandvals,psym=4,color=255
        splog,'ifail = '+string(ifail)
        if(keyword_set(qastop)) then stop
    endif
endfor

end
;------------------------------------------------------------------------------
