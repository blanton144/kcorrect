;+
; NAME:
;   lf_calc_vmax
; PURPOSE:
;   calculate vmax for an object, given the flux and redshift limits
; USAGE:
;   lf_calc_vmax,appm,
;   absm,coeffs,filtername,marea,mmin,mmax,sample_zmin, $
;     sample_zmax, [, band_shift=, im=, omega0=, omegal0=, vmax= ]
; INPUTS:
;   absm             [N] absolute magnitudes
;   coeffs           [NT, N] K-correction coeffs for each galaxy
;   filtername       name of filter to use
;   marea            [NA] area of mag limit region
;   mmin             [NA] minimum apparent mag for each mag limit region
;   mmax             [NA] maximum apparent mag for each mag limit region
;   sample_zmin      minimum redshift of sample
;   sample_zmax      maximum redshift of sample
; OPTIONAL INPUTS:
;   band_shift       shift in defn of absm 
;   im               [N] index of magnitude limit region (default 0)
;   omega0           omega_matter to use (default: 0.3)
;   omegal0          omega_lambda to use (default: 0.7)
; OUTPUTS:
;   zmin             [N] local zmin of each
;   zmax             [N] local zmax of each
;   vmax             [N] vmax of each integrated over all mag limit regions
; COMMENTS:
;   vmax is returned in h^{-3} Mpc^3 comoving
; DEPENDENCIES:
;   idlutils
;   kcorrect (in kcorrect product)
; REVISION HISTORY:
;   2002-11-17  written - Blanton
;   2008-09-04  added VNAME optional input - J. Moustakas
;-
pro lf_calc_vmax,appm,absm,coeffs,in_filtername,marea,mmin,mmax, $
                 sample_zmin,sample_zmax,im=im, vmax=vmax, omega0=in_omega0, $
                 omegal0=in_omegal0,zmin=zmin, zmax=zmax, $
                 band_shift=in_band_shift, q0=in_q0, q1=in_q1, qz0=in_qz0, $
                 actual_z=actual_z, actual_k=actual_k, $
                 magoffset=in_magoffset, absmagdep=absmagdep, $
                 ref_absmagdep=ref_absmagdep, rmatrix=rmatrix, zvals=zvals, $
                 smoothr=smoothr,vname=vname

; settings
pi=3.14159265358979D
dh=2.99792D+5/100.D

; defaults
ngals=n_elements(absm)
if(n_elements(in_omega0) eq 0) then in_omega0=0.3
if(n_elements(in_omegal0) eq 0) then in_omegal0=0.7
if(n_elements(in_band_shift) eq 0) then in_band_shift=0.
if(n_elements(in_q0) eq 0) then in_q0=0.
if(n_elements(in_q1) eq 0) then in_q1=0.
if(n_elements(in_qz0) eq 0) then in_qz0=0.
if(n_elements(in_magoffset) eq 0) then in_magoffset=0.
if(n_elements(im) eq 0) then im=lonarr(ngals)
omega0=in_omega0
omegal0=in_omegal0
filtername=in_filtername
band_shift=in_band_shift
q0=in_q0
q1=in_q1
qz0=in_qz0
magoffset=in_magoffset
nv=n_elements(coeffs)/n_elements(absm)
nk=1

; run kcorrect to get rmatrix and zvals
if(n_elements(rmatrix) eq 0 OR n_elements(zvals) eq 0) then $
  kcorrect,dummaggies,dummaggies_ivar,0.,dumk,band_shift=band_shift, $
  filterlist=[filtername],rmatrix=rmatrix,zvals=zvals,coeffs=coeffs[*,0],$
  vname=vname
nz=n_elements(rmatrix)/nv/nk

if(keyword_set(smoothr)) then begin
    old_rmatrix=rmatrix
    for i=0L, nv-1L do $
      rmatrix[*,i]=smooth(old_rmatrix[*,i], 4., /edge)
endif

; for each object...
vmax=fltarr(ngals)
zmin=fltarr(ngals)
zmax=fltarr(ngals)
        help, omega0
        help, omegal0
for i=0L,ngals-1L do begin
    if((i mod 1000) eq 0) then splog,' galaxy '+string(i)
    curr_absm=k_evolve(absm[i], actual_z[i], q0, q1, qz0)
    curr_coeffs=coeffs[*,i]
    if(keyword_set(smoothr)) then begin
        off=-2.5*alog10(interpol(rmatrix#abs(curr_coeffs), zvals, $
                                 actual_z[i]))+ $
          2.5*alog10(interpol(old_rmatrix#abs(curr_coeffs), zvals, $
                                 actual_z[i]))
        curr_absm=curr_absm-off
    endif
    
    for j=0L, n_elements(marea)-1L do begin

        curr_zmin=-1.
        curr_zmax=-1.
        soname=filepath('libkcorrect.'+kcorrect_so_ext(), $
                        root_dir=getenv('KCORRECT_DIR'), subdirectory='lib')
        curr_mmin=mmin[j]
        curr_mmax=mmax[j]
        retval=call_external(soname, 'idl_lf_calc_vmax', float(curr_absm), $
                             float(curr_coeffs),long(nv),float(zvals), $
                             long(nz),float(rmatrix),long(nk), $
                             float(sample_zmin),float(sample_zmax), $
                             float(curr_mmin),float(curr_mmax), $
                             float(q0),float(q1),float(qz0), $
                             float(band_shift),float(magoffset), $
                             float(omega0),float(omegal0), $
                             float(curr_zmin),float(curr_zmax))
        
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; calculate vmax; don't allow it to be less than zero
        if(im[i] eq j) then begin
            zmin[i]=curr_zmin
            zmax[i]=curr_zmax
            if(zmin[i] gt actual_z[i]) then begin
                if(abs(curr_mmin-appm[i]) gt 0.005 AND $
                   abs(zmin[i]-actual_z[i]) gt 1.e-4) then begin
                   message,'galaxy outside allowed putative selection limits', /info
;                  splog, i, actual_z[i], zmin[i], zmax[i], appm[i], curr_mmin, curr_mmax, curr_absm
                endif
            endif
            if(zmax[i] lt actual_z[i]) then begin
                if(abs(curr_mmax-appm[i]) gt 0.005 AND $
                  abs(zmax[i]-actual_z[i]) gt 1.e-4)then begin
                   message,'galaxy outside allowed putative selection limits', /info
;                  splog, i, actual_z[i], zmin[i], zmax[i], appm[i], curr_mmin, curr_mmax, curr_absm
;                  print, appm[i], absm[i] + lf_distmod(actual_z[i]) + actual_k[i]
                endif
            endif
        endif
        vmax[i]=vmax[i]+ $
          ((marea[j]*(lf_comvol(curr_zmax, omega0=omega0, omegal0=omegal0)- $
                      lf_comvol(curr_zmin, omega0=omega0, omegal0=omegal0))) $
           > 0.)/3.
    endfor

endfor

end
