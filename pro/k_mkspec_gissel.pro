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
pro k_mkspec_gissel, vmatrix, lambda, metallicity, dust, sfhtype, sfhpars, $
                     gisselpath=gisselpath, attime=attime, maxage=maxage, $
                     minage=minage

; Need at least 6 parameters
if (N_params() LT 5) then begin
    print, 'Syntax - k_mkspec_gissel, vmatrix, metallicity, dust, sfhtype, sfhpars, $'
    print, '         [, gisselpath=]'
    return
endif

pi=3.14159265358979D
littleh=0.7
km2Mpc=1./3.086d+19
s2yr=1./3600./24./365.25
hubble_time=s2yr/(littleh*100.*km2Mpc)

if(NOT keyword_set(gisselpath)) then $
  gisselpath='/data/specmodels/gissel/data'
if(NOT keyword_set(attime)) then attime=0.d
if(NOT keyword_set(nl)) then nl=3000L
if(NOT keyword_set(lmin)) then lmin=1250.
if(NOT keyword_set(lmax)) then lmax=33333.

nv=n_elements(metallicity)*n_elements(sfhpars)*n_elements(dust)
vmatrix=dblarr(nl,nv)
for m = 0, n_elements(metallicity)-1L do begin
    tsteps=dblarr(221)
    openr,unit,gisselpath+'/ssp_salp_z'+metallicity[m]+'.tsteps.dat',/get_lun
    readf,unit,tsteps
    free_lun,unit
    dt=dblarr(n_elements(tsteps))
    dt[0]=0.
    dt[1:n_elements(dt)-1]=tsteps[1:n_elements(dt)-1]- $
      tsteps[0:n_elements(dt)-2]
    help,tsteps
    if(NOT file_test('data_metspec.'+metallicity[m]+'.sav')) then begin
        metspec=dblarr(nl,n_elements(tsteps))
        for i=0, n_elements(tsteps)-1 do begin
            help,i
            k_load_templates,gisselpath+'/ssp_salp_z'+metallicity[m]+ $
              '/ssp_salp_z'+metallicity[m]+'.flux.'+ $
              string(i,format='(i4.4)')+'.dat',tmpspec,lambda, $
              nl=nl,lmin=lmin,lmax=lmax
            metspec[*,i]=tmpspec
        endfor
        save,metspec,filename='data_metspec.'+metallicity[m]+'.sav'
    endif else begin
        i=0
        k_load_templates,gisselpath+'/ssp_salp_z'+metallicity[m]+ $
          '/ssp_salp_z'+metallicity[m]+'.flux.'+ $
          string(i,format='(i4.4)')+'.dat',tmpspec,lambda, $
          nl=nl,lmin=lmin,lmax=lmax
        restore,'data_metspec.'+metallicity[m]+'.sav'
    endelse
    if(sfhtype eq 'gaussian') then begin
        for a = 0, n_elements(sfhpars)-1L do begin
            contrib=dt*exp(-0.5*((sfhpars[a].age-attime-tsteps)/ $
                                 sfhpars[a].agesigma)^2)/ $
              sqrt(2.*pi*sfhpars[a].agesigma^2)
            if(keyword_set(maxage)) then begin
                indx=where(tsteps gt maxage-attime,count)
                if(count gt 0) then contrib[indx]=0.d
            endif
            if(keyword_set(minage)) then begin
                indx=where(tsteps lt minage-attime,count)
                if(count gt 0) then contrib[indx]=0.d
            endif
            spec=dblarr(nl)
            for i = 0, n_elements(tsteps)-1 do $
              spec=spec+metspec[*,i]*contrib[i]
            for d = 0, n_elements(dust)-1L do begin
                vindx=m*n_elements(sfhpars)*n_elements(dust)+ $
                  a*n_elements(dust)+d
                vmatrix[*,vindx]=spec* $
                  exp(-witt_ext(dust[d],dust[d].tauv,lambda[0:nl-1]))
            endfor
        endfor
    endif else if(sfhtype eq 'zgaussian') then begin
        zarr=max(sfhpars.z)*(dindgen(1000)+0.5)/1000.
        tarr=hubble_time*lookback(zarr,0.3,0.7)
        ztsteps=interpol(zarr,tarr,tsteps)
        for a = 0, n_elements(sfhpars)-1L do begin
            contrib=dt*exp(-0.5*((sfhpars[a].z-ztsteps)/ $
                                 sfhpars[a].zsigma)^2)/ $
              sqrt(2.*pi*sfhpars[a].zsigma^2)
            spec=dblarr(nl)
            for i = 0, n_elements(tsteps)-1 do $
              spec=spec+metspec[*,i]*contrib[i]
            for d = 0, n_elements(dust)-1L do begin
                vindx=m*n_elements(sfhpars)*n_elements(dust)+ $
                  a*n_elements(dust)+d
                vmatrix[*,vindx]=spec* $
                  exp(-witt_ext(dust[d],dust[d].tauv,lambda[0:nl-1]))
            endfor
        endfor
    endif else if (sfhtype eq 'exponential') then begin
        for a = 0, n_elements(sfhpars)-1L do begin
            contrib=dblarr(n_elements(tsteps))
            indx=where(tsteps lt sfhpars[a].maxage)
            contrib[indx]=exp(tsteps/sfhpars[a].agescale)
            contrib=contrib/total(contrib,/double)
            indx=where(contrib gt 1.d-5, count)
            if(count gt 0) then begin
                spec=dblarr(nl)
                for i = 0, count-1 do $
                  spec=spec+metspec[*,indx[i]]*contrib[indx[i]]
                for d = 0, n_elements(dust)-1L do begin
                    vindx=m*n_elements(sfhpars)*n_elements(dust)+ $
                      a*n_elements(dust)+d
                    vmatrix[*,vindx]=spec*10.^(-0.4*dust[d]* $
                                               ext_ccm(lambda[0:nl-1],3.5))
                endfor
            endif
        endfor
    endif 
endfor

end
;------------------------------------------------------------------------------
