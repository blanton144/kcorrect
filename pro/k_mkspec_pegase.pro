;+
; NAME:
;   k_mkspec_pegase
; PURPOSE:
;   Read in PEGASE.2 models and create vmatrix; requires that 
;   mrb_spectra.[metallicity].dat files be existent in the PEGASE 
;   directory
; CALLING SEQUENCE:
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
;   read_peg (in eplusa)
; REVISION HISTORY:
;   25-Jul-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_mkspec_pegase, vmatrix, lambda, metallicity, dust, sfhtype, sfhpars, $
                     pegasepath=pegasepath, attime=attime, maxage=maxage, $
                     minage=minage, lmin=lmin, lmax=lmax, nl=nl, $
                     nolines=nolines

; Need at least 6 parameters
if (N_params() LT 5) then begin
    print, 'Syntax - k_mkspec_pegase, vmatrix, metallicity, dust, sfhtype, sfhpars, $'
    print, '         [, pegasepath=]'
    return
endif

; defaults
if(NOT keyword_set(pegasepath)) then $
  pegasepath=getenv('DATA')+'/specmodels/PEGASE.2'
if(NOT keyword_set(attime)) then attime=0.d
if(NOT keyword_set(ntsteps)) then ntsteps=5000L
if(NOT keyword_set(nl)) then nl=5000L
if(NOT keyword_set(lmin)) then lmin=1250.
if(NOT keyword_set(lmax)) then lmax=33333.

nv=n_elements(metallicity)*n_elements(sfhpars)*n_elements(dust)
vmatrix=dblarr(nl,nv)
for m = 0, n_elements(metallicity)-1L do begin
;   read in the pegase file
    pegfile=pegasepath+'/mrb_spectra.0.'+metallicity[m]+'.dat'
    read_peg,pegfile,peg=peg
    k_spec_pegase,peg[0],dumspec,lambda,nl=nl,lmin=lmin,lmax=lmax, $
      nolines=nolines

;   set up times
    if(NOT file_test('data_metspec.'+metallicity[m]+'.sav')) then begin
        time=(min(peg.arr1)+(max(peg.arr1)-min(peg.arr1))*dindgen(ntsteps)/ $
              double(ntsteps-1L))*1.D+6
        dt=time[1]-time[0]
        metspec=dblarr(nl,ntsteps)
        for i=0, ntsteps-1 do $
          metspec[*,i]=k_interp_pegase(peg,time[i]*1.D-6,nl=nl,lmin=lmin, $
                                       lmax=lmax,nolines=nolines)
        save,time,dt,metspec,filename='data_metspec.'+metallicity[m]+'.sav'
    endif else begin
        restore,'data_metspec.'+metallicity[m]+'.sav'
    endelse

    for a = 0, n_elements(sfhpars)-1L do begin
        contrib=dt*exp(-0.5*((sfhpars[a].age-attime-time)/ $
                             sfhpars[a].agesigma)^2)/ $
          sqrt(2.*!DPI*sfhpars[a].agesigma^2)
        if(keyword_set(maxage)) then begin
            indx=where(time gt maxage-attime,count)
            if(count gt 0) then contrib[indx]=0.d
        endif
        if(keyword_set(minage)) then begin
            indx=where(time lt minage-attime,count)
            if(count gt 0) then contrib[indx]=0.d
        endif
        spec=dblarr(nl)
        for i = 0, n_elements(time)-1 do $
          spec=spec+metspec[*,i]*contrib[i]
        for d = 0, n_elements(dust)-1L do begin
            vindx=m*n_elements(sfhpars)*n_elements(dust)+ $
              d*n_elements(sfhpars)+a
            vmatrix[*,vindx]=spec* $
              exp(-witt_ext(dust[d],dust[d].tauv,lambda[0:nl-1]))
        endfor
    endfor
endfor

end
;------------------------------------------------------------------------------
