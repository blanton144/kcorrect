;+
; NAME:
;   k_make_ktable
;
; PURPOSE:
;   Makes a ktable as used by the LSS sample
;
; CALLING SEQUENCE:
;   k_make_table, ktablebase, to_z [, version=version, vpath=vpath]
;
; INPUTS:
;   ktablebase - output file for table
;   to_z - redshift to correct to
;
; OPTIONAL INPUTS:
;   version   - version of eigentemplates to use
;   vpath   - version of eigentemplates to use
;
; OUTPUTS:
;   reconstruct_maggies  - maggies in each band
;
; OPTIONAL INPUT/OUTPUTS:
;   ematrix       - eigentemplates [N_dim, N_template]
;   filterlist    - list of files with filter information [N_band]
;   bmatrix       - orthogonal templates spanning SED space [N_lambda, N_dim]
;   lambda        - wavelengths for orthogonal templates [N_lambda]
;   rmatrix       - look up table for bmatrix and filter information 
;                   [N_z, N_dim, N_band]
;   zvals         - look up table for rmatrix [N_z]
;
; COMMENTS:
; 
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   16-Feb-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_make_ktable,ktablebase,to_z,ztable=ztable,ematrix=ematrix,zvals=zvals, $
                  version=version,vpath=vpath,filterpath=filterpath, $
                  filterlist=filterlist, $
                  rmatrix=rmatrix,bmatrix=bmatrix,lambda=lambda, $
                  rematrix=rematrix,rzvals=rzvals, $
                  rversion=rversion,rvpath=rvpath,rfilterpath=rfilterpath, $
                  rfilterlist=rfilterlist, $
                  rrmatrix=rrmatrix,rbmatrix=rbmatrix,rlambda=rlambda, $
                  qematrix=qematrix,qzvals=qzvals, $
                  qrmatrix=qrmatrix,qbmatrix=qbmatrix,qlambda=qlambda, $
                  qversion=qversion,qvpath=qvpath,qfilterpath=qfilterpath, $
                  qfilterlist=qfilterlist

; Need at least 1 parameters
if (N_params() LT 2) then begin
    print, 'Syntax - k_make_ktable, ktablebase, to_z, [, ematrix=, zvals=, rmatrix=, '
    print, '   bmatrix=, lambda=, version=, vpath=, filterlist=, filterpath=]'
    return
endif

if(keyword_set(ematrix)) then rematrix=ematrix
if(keyword_set(zvals)) then rzvals=zvals
if(keyword_set(version)) then rversion=version
if(keyword_set(vpath)) then rvpath=vpath
if(keyword_set(filterpath)) then rfilterpath=filterpath
if(keyword_set(filterlist)) then rfilterlist=filterlist
if(keyword_set(rmatrix)) then rrmatrix=rmatrix
if(keyword_set(bmatrix)) then rbmatrix=bmatrix
if(keyword_set(lambda)) then rlambda=lambda

if(NOT keyword_set(rfilterpath)) then $
  rfilterpath=getenv('KCORRECT_DIR')+'/data/filters'

; Get bmatrix and stuff from files if necessary
if(keyword_set(rversion)) then begin 
    if(NOT keyword_set(rvpath)) then $
      rvpath=getenv('KCORRECT_DIR')+'/data/etemplates'
    k_load_ascii_table,rematrix,rvpath+'/ematrix.'+rversion+'.dat'
    nt=(size(rematrix))[2]
    nb=(size(rematrix))[1]
    k_load_ascii_table,rbmatrix,rvpath+'/bmatrix.'+rversion+'.dat'
    nl=(size(rbmatrix))[1]
    k_load_ascii_table,rlambda,rvpath+'/lambda.'+rversion+'.dat'
    rfiltfile=rvpath+'/filterlist.'+rversion+'.dat'
    spawn,'cat '+rfiltfile+' | wc -l',nfilters
    nk=long(nfilters[0])-1l
    rfilterlist=strarr(nk)
    openr,unit,rvpath+'/filterlist.'+rversion+'.dat',/get_lun
    readf,unit,nk
    readf,unit,rfilterlist
    close,unit
    free_lun,unit
endif

if (keyword_set(rbmatrix) AND keyword_set(rfilterlist)  $
    AND keyword_set(rlambda)) then begin
    nk=long(n_elements(rfilterlist))
    k_create_r,rrmatrix,rbmatrix,rlambda,rzvals,rfilterlist, $
      filterpath=rfilterpath
endif else begin
    if (NOT keyword_set(rrmatrix) OR NOT keyword_set(rzvals)) then begin
        klog, 'need to specify rmatrix or bmatrix and filterlist'
    endif
endelse 
nz=long(n_elements(rzvals))

; Get bmatrix and stuff from files if necessary
if(NOT keyword_set(qversion)) then qversion=rversion
if(NOT keyword_set(qvpath)) then qvpath=rvpath
if(NOT keyword_set(qfilterpath)) then qfilterpath=rfilterpath
if(keyword_set(qversion)) then begin 
    if(NOT keyword_set(qvpath)) then $
      qvpath=getenv('KCORRECT_DIR')+'/data/etemplates'
    k_load_ascii_table,qematrix,qvpath+'/ematrix.'+qversion+'.dat'
    nt=(size(qematrix))[2]
    nb=(size(qematrix))[1]
    k_load_ascii_table,qbmatrix,qvpath+'/bmatrix.'+qversion+'.dat'
    nl=(size(qbmatrix))[1]
    k_load_ascii_table,qlambda,qvpath+'/lambda.'+qversion+'.dat'
    qfiltfile=qvpath+'/filterlist.'+qversion+'.dat'
    spawn,'cat '+qfiltfile+' | wc -l',nfilters
    nk=long(nfilters[0])-1l
    qfilterlist=strarr(nk)
    openr,unit,qvpath+'/filterlist.'+qversion+'.dat',/get_lun
    readf,unit,nk
    readf,unit,qfilterlist
    close,unit
    free_lun,unit
endif

if(NOT keyword_set(qematrix)) then qematrix=rematrix
if(NOT keyword_set(qbmatrix)) then qbmatrix=rbmatrix
if(NOT keyword_set(qlambda)) then qlambda=rlambda
if(NOT keyword_set(qfilterlist)) then qfilterlist=rfilterlist

if (keyword_set(qbmatrix) AND keyword_set(qfilterlist)  $
    AND keyword_set(qlambda)) then begin
    nk=long(n_elements(qfilterlist))
    k_create_r,qrmatrix,qbmatrix,qlambda,qzvals,qfilterlist, $
      filterpath=qfilterpath
endif else begin
    if (NOT keyword_set(qrmatrix) OR NOT keyword_set(qzvals)) then begin
        klog, 'need to specify rmatrix or bmatrix and filterlist'
    endif
endelse 
nz=long(n_elements(qzvals))

if(NOT keyword_set(qzvals)) then qzvals=rzvals
if(NOT keyword_set(qrmatrix)) then qrmatrix=rrmatrix

if(NOT keyword_set(ztable)) then ztable=rzvals

obs_maggies=dblarr(nk,nt,nz)
fix_maggies=dblarr(nk,nt,nz)

zero_z=replicate(0.,nz)
fixed_z=replicate(to_z,nz)
for t=0, nt-1 do begin
    coeff=dblarr(nt,nz)
    coeff[t,*]=1.d
    k_reconstruct_maggies,coeff,ztable,obs_maggies_temp,rmatrix=rrmatrix, $
      zvals=rzvals,ematrix=rematrix,band_shift=zero_z
    k_reconstruct_maggies,coeff,zero_z,fix_maggies_temp,rmatrix=qrmatrix, $
      zvals=qzvals,ematrix=qematrix,band_shift=fixed_z
    obs_maggies[*,t,*]=obs_maggies_temp
    fix_maggies[*,t,*]=fix_maggies_temp
endfor

for k=0, nk-1 do begin
    openw,unit,ktablebase+'.'+strtrim(string(k),2)+'.dat',/get_lun
    outarr=dblarr(1+nt*2,nz)
    outarr[0,*]=ztable
    outarr[1+2*lindgen(nt),*]=fix_maggies[k,*,*]
    outarr[2+2*lindgen(nt),*]=obs_maggies[k,*,*]
    printf,unit,format='('+strtrim(string(1+nt*2),2)+'e14.6)',outarr
    close,unit
    free_lun,unit
endfor

end
