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
pro k_make_ktable,ktablebase,to_z,ematrix=ematrix,zvals=zvals,rmatrix=rmatrix,bmatrix=bmatrix,lambda=lambda,version=version,vpath=vpath,filterpath=filterpath,filterlist=filterlist

; Need at least 1 parameters
if (N_params() LT 2) then begin
    print, 'Syntax - k_make_ktable, ktablebase, to_z, [, ematrix=, zvals=, rmatrix=, '
    print, '   bmatrix=, lambda=, version=, vpath=, filterlist=, filterpath=]'
    return
endif

if(NOT keyword_set(filterpath)) then $
  filterpath=getenv('KCORRECT_DIR')+'/data/filters'

; Get bmatrix and stuff from files if necessary
if(keyword_set(version)) then begin 
    if(NOT keyword_set(vpath)) then $
      vpath=getenv('KCORRECT_DIR')+'/data/etemplates'
    k_load_ascii_table,ematrix,vpath+'/ematrix.'+version+'.dat'
    nt=(size(ematrix))[2]
    nb=(size(ematrix))[1]
    k_load_ascii_table,bmatrix,vpath+'/bmatrix.'+version+'.dat'
    nl=(size(bmatrix))[1]
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

if (keyword_set(bmatrix) AND keyword_set(filterlist)  $
    AND keyword_set(lambda)) then begin
    nk=long(n_elements(filterlist))
    k_create_r,rmatrix,bmatrix,lambda,zvals,filterlist,filterpath=filterpath
endif else begin
    if (NOT keyword_set(rmatrix) OR NOT keyword_set(zvals)) then begin
        klog, 'need to specify rmatrix or bmatrix and filterlist'
    endif
endelse 
nz=long(n_elements(zvals))

obs_maggies=dblarr(nk,nt,nz)
fix_maggies=dblarr(nk,nt,nz)

fixed_z=replicate(to_z,nz)
for t=0, nt-1 do begin
    coeff=dblarr(nt,nz)
    coeff[t,*]=1.d
    k_reconstruct_maggies,coeff,zvals,obs_maggies_temp,rmatrix=rmatrix, $
      zvals=zvals,ematrix=ematrix
    k_reconstruct_maggies,coeff,fixed_z,fix_maggies_temp,rmatrix=rmatrix, $
      zvals=zvals,ematrix=ematrix
    obs_maggies[*,t,*]=obs_maggies_temp
    fix_maggies[*,t,*]=fix_maggies_temp
endfor

for k=0, nk-1 do begin
    openw,unit,ktablebase+'.'+strtrim(string(k),2)+'.dat',/get_lun
    outarr=dblarr(1+nt*2,nz)
    outarr[0,*]=zvals
    outarr[1+2*lindgen(nt),*]=fix_maggies[k,*,*]
    outarr[2+2*lindgen(nt),*]=obs_maggies[k,*,*]
    printf,unit,format='('+strtrim(string(1+nt*2),2)+'e14.6)',outarr
    close,unit
    free_lun,unit
endfor

end
