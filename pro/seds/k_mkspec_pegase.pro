;+
; NAME:
;   k_mkspec_pegase
; PURPOSE:
;   read in PEGASE.2 models, add dust, convolve with SFH
; CALLING SEQUENCE:
;   k_mkspec_pegase, vmatrix, lambda, metallicity, dust, sfhpars [, $
;      pegasepath=, attime=, maxage=, minage=, lmin=, lmax=, nl=, /nolines]
; INPUTS:
;   metallicity - string with metallicity name
;   dust - witt dust name 
;   sfhpars - structure describing SFH
; OPTIONAL INPUTS:
;   pegasepath - directory with PEGASE results (default
;                $DATA/specmodels/PEGASE.2)
;   attime - observe at look back time of this (default 0.)
;   [min|max]age - minimum and maximum ages of populations (default none)
;   l[min|max] - wavelength limits (default 1250., 33333.)
;   nl - number of pixels (default 5000)
;   pversion - prefix for spectra names (default 'mrb')
; KEYWORDS:
;   /nolines - don't use the lines
; OUTPUTS:
;   vmatrix - [nl, n_spectra] output spectra
;   lambda - [nl+1] pixel edges
; COMMENTS:
;   Requires you have created PEGASE.2 models in the files:
;     $DATA/specmodels/PEGASE.2/[pversion]_spectra.0.[metallicity].dat
; PROCEDURES CALLED:
;   k_load_ascii_table
;   k_read_peg 
; REVISION HISTORY:
;   25-Jul-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_mkspec_pegase, vmatrix, lambda, metallicity, dust, sfhpars, $
                     pegasepath=pegasepath, attime=attime, maxage=maxage, $
                     minage=minage, lmin=lmin, lmax=lmax, nl=nl, $
                     nolines=nolines, pversion=pversion

; Need at least 6 parameters
if (N_params() LT 5) then begin
    print, 'Syntax - k_mkspec_pegase, vmatrix, metallicity, dust, sfhpars, $'
    print, '         [, pegasepath=, pversion=]'
    return
endif

; defaults
if(NOT keyword_set(pegasepath)) then $
  pegasepath=getenv('DATA')+'/specmodels/PEGASE.2'
if(NOT keyword_set(attime)) then attime=0.d
if(NOT keyword_set(ntsteps)) then ntsteps=20000L
if(NOT keyword_set(nl)) then nl=5000L
if(NOT keyword_set(lmin)) then lmin=1250.
if(NOT keyword_set(lmax)) then lmax=33333.
if(NOT keyword_set(pversion)) then pversion='mrb'

nv=n_elements(metallicity)*n_elements(sfhpars)*n_elements(dust)* $
  n_elements(pversion)
vmatrix=fltarr(nl,nv)
for p = 0, n_elements(pversion)-1L do begin
    for m = 0, n_elements(metallicity)-1L do begin
;   read in the pegase file
        pegfile=pegasepath+'/'+pversion[p]+'_spectra.0.'+metallicity[m]+ $
          '.dat'
        k_read_peg,pegfile,peg=peg
        k_spec_pegase,peg[0],dumspec,lambda,nl=nl,lmin=lmin,lmax=lmax, $
          nolines=nolines

;   set up times
        if(NOT file_test('data_metspec'+pversion[p]+'.'+ $
                         metallicity[m]+'.sav')) then begin
            time=(min(peg.arr1)+(max(peg.arr1)-min(peg.arr1))* $
                  findgen(ntsteps)/float(ntsteps-1L))*1.e+6
            dt=time[1]-time[0]
            metspec=fltarr(nl,ntsteps)
            for i=0, ntsteps-1 do $
              metspec[*,i]=k_interp_pegase(peg,time[i]*1.e-6,nl=nl,lmin=lmin, $
                                           lmax=lmax,nolines=nolines)
            save,time,dt,metspec,filename='data_metspec'+pversion[p]+'.'+ $
              metallicity[m]+'.sav'
        endif else begin
            restore,'data_metspec'+pversion[p]+'.'+metallicity[m]+'.sav'
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
            spec=fltarr(nl)
            for i = 0, n_elements(time)-1 do $
              spec=spec+metspec[*,i]*contrib[i]
            for d = 0, n_elements(dust)-1L do begin
                vindx=p*n_elements(metallicity)*n_elements(sfhpars)* $
                  n_elements(dust)+ $
                  m*n_elements(sfhpars)*n_elements(dust)+ $
                  d*n_elements(sfhpars)+a
                vmatrix[*,vindx]=spec* $
                  exp(-witt_ext(dust[d],dust[d].tauv,lambda[0:nl-1]))
            endfor
        endfor
    endfor
endfor

end
;------------------------------------------------------------------------------
