;+
; NAME:
;   k_mkvmatrix_gissel
;
; PURPOSE:
;   Read in GISSEL models and create vmatrix models
;
; CALLING SEQUENCE:
;
; INPUTS:
; OPTIONAL INPUTS:
; KEYWORDS:
; OUTPUTS:
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
;   k_load_ascii_table
; REVISION HISTORY:
;   25-Jul-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_mkvmatrix_gissel, vmatrix, lambda, metallicity, dust, age, agesigma, $
                        gisselpath=gisselpath, nl=nl

pi=3.14159265358979d

; Need at least 6 parameters
if (N_params() LT 5) then begin
    print, 'Syntax - k_mkvmatrix, vmatrix, metallicity, dust, age, agesigma $'
    print, '         [, gisselpath=]'
    return
endif

if(NOT keyword_set(gisselpath)) then $
  gisselpath='/data/specmodels/gissel/data'
if(NOT keyword_set(nl)) then $
  nl=3000L
if(NOT keyword_set(lmin)) then $
  lmin=1250.
if(NOT keyword_set(lmax)) then $
  lmax=33333.

nv=n_elements(metallicity)*n_elements(age)*n_elements(dust)
vmatrix=dblarr(nl,nv)
for m = 0, n_elements(metallicity)-1L do begin
    tsteps=dblarr(221)
    openr,unit,gisselpath+'/ssp_salp_z'+metallicity[m]+'.tsteps.dat',/get_lun
    readf,unit,tsteps
    free_lun,unit
    metspec=dblarr(nl,n_elements(tsteps))
    help,tsteps
    for i=0, n_elements(tsteps)-1 do begin
        help,i
        k_load_templates,gisselpath+'/ssp_salp_z'+metallicity[m]+ $
          '/ssp_salp_z'+metallicity[m]+'.flux.'+ $
          string(i,format='(i4.4)')+'.dat',tmpspec,lambda, $
          nl=nl,lmin=lmin,lmax=lmax
        metspec[*,i]=tmpspec
    endfor
    for a = 0, n_elements(age)-1L do begin
        scale=(1./(sqrt(2.*pi)*agesigma[a]))
        contrib=scale*exp(-0.5*((age[a]-tsteps)/agesigma[a])^2)
        indx=where(contrib gt 1.d-4*scale, count)
        if(count gt 0) then begin
            spec=dblarr(nl)
            for i = 0, count-1 do spec=spec+metspec[*,indx[i]]*contrib[indx[i]]
            for d = 0, n_elements(dust)-1L do begin
                vindx=m*n_elements(age)*n_elements(dust)+a*n_elements(dust)+d
                vmatrix[*,vindx]=spec*10.^(-0.4*dust[d]* $
                                           ext_ccm(lambda[0:nl-1],3.5))
            endfor
        endif
    endfor
endfor

end
;------------------------------------------------------------------------------
