;+
; NAME:
;   k_project_filters
;
; PURPOSE:
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL INPUT/OUTPUTS:
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
;   17-Jan-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function k_project_filters,lambda,flux,filterlist=filterlist, $
                           filterpath=filterpath,band_shift=band_shift, $
                           qa=qa

if(NOT keyword_set(filterpath)) then $
  filterpath=getenv('KCORRECT_DIR')+'/data/filters'
if(NOT keyword_set(filterlist)) then $
  filterlist=['sdss_u0','sdss_g0','sdss_r0','sdss_i0','sdss_z0']
if(NOT keyword_set(band_shift)) then band_shift=0.

; load in the filters
filterlist=filterpath+'/'+filterlist+'.dat'
k_load_filters,filterlist,filter_n,filter_lambda,filter_pass
filter_lambda=filter_lambda/(1.+band_shift)

; calculate differential for each point
nlambda=n_elements(lambda)
nspectra=n_elements(flux)/(nlambda-1)
dlambda1=dblarr(nlambda)
dlambda1[0]=0.5*(lambda[1]-lambda[0])
dlambda1[nlambda-1]=0.5*(lambda[nlambda-1]-lambda[nlambda-2])
dlambda1[1:nlambda-2]=0.5*(lambda[2:nlambda-1]-lambda[0:nlambda-3])
dlambda=dlambda1#replicate(1.,nspectra)

; flat spectrum
cspeed=2.99792e+18 ; ang per sec
abspec=3.631d-20
flat=cspeed*abspec/lambda^2

maggies=dblarr(n_elements(filterlist),nspectra)
for i=0L, n_elements(filterlist)-1L do begin
    ; interpolate filter onto lambda
    isort=sort(filter_lambda[0:filter_n[i]-1L,i])
    intfilter1=interpol((filter_pass[0:filter_n[i]-1L,i])[isort], $
                        (filter_lambda[0:filter_n[i]-1L,i])[isort], $
                        lambda)
    filtlammin=filter_lambda[0,i] < filter_lambda[filter_n[i]-1,i]
    filtlammax=filter_lambda[0,i] > filter_lambda[filter_n[i]-1L,i]
    indx=where(intfilter1 lt 0. or $
               lambda gt filtlammax or $
               lambda lt filtlammin, count)
    if(count gt 0) then intfilter1[indx]=0.
    intfilter=intfilter1#replicate(1.,nspectra)

    ; now total to get abfluxes
    lambdafull=lambda#replicate(1.,nspectra)
    abfluxes=total(lambdafull*dlambda*intfilter*flux,1,/double)

    ; normalize to our fav source
    absource=total(lambda*dlambda1*intfilter1*flat,/double)

    ; now make maggies
    maggies[i,*]=abfluxes/absource

    if(keyword_set(qa)) then begin
        plot,lambda,intfilter1,xra=[filtlammin,filtlammax]
        plot,lambda,flux,xra=[filtlammin,filtlammax]
    endif
endfor

return,maggies

end
;------------------------------------------------------------------------------

