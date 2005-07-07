;+
; NAME:
;   lf_select_eep
; PURPOSE:
;   calculate EEP selection function from luminosity function
; USAGE:
;   lf_select_eep,zz,absmag,kcorrect,mmin,mmax,sample_absmmin,sample_absmmax,
;       sample_zmin,sample_zmax, absmk, phi, sel [, q0=, $
;       omega0=, omegal0=, q1=, qz0=]
; INPUTS:
;   zz               [N] redshifts
;   absmag           [N] absolute magnitudes
;   kcorrect         [N] K-correction for each galaxy
;   mmin             [N] minimum apparent mag for each galaxy
;   mmax             [N] maximum apparent mag for each galaxy
;   sample_absmmin   observed absolute magnitude minimum of sample
;   sample_absmmax   observed absolute magnitude maximum of sample
;   absmk            [nbin] center of each bin
;   phi              [nbin] amplitude of each bin
; OPTIONAL INPUTS:
;   omega0           omega_matter to use (default: 0.3)
;   omegal0          omega_lambda to use (default: 0.7)
;   q0,q1            evolution (default 0.)
;   qz0              pivot for redshift evolution (default 0.)
; OUTPUTS:
;   sel              [N] selection function
; REVISION HISTORY:
;   2003-10-20  written - Blanton
;-
pro lf_select_eep,zz,absmag,kcorrect,mmin,mmax,sample_absmmin,sample_absmmax, $
                  absmk,phi,sel, q0=q0, q1=q1, qz0=qz0, $
                  omega0=omega0,omegal0=omegal0
                 
; set defaults
if(n_elements(omega0) eq 0) then omega0=0.3
if(n_elements(omegal0) eq 0) then omegal0=0.7
if(n_elements(q0) eq 0) then q0=0.
if(n_elements(q1) eq 0) then q1=0.
if(n_elements(qz0) eq 0) then qz0=0.

nbin=n_elements(phi)
ngals=n_elements(zz)
sel=fltarr(ngals)

; calculate absmag limits (without evolution, which is accounted for within)
dm=lf_distmod(zz,omega0=omega0,omegal0=omegal0)
dmK=dm+kcorrect
absmmin=k_evolve((mmin-dmK), zz, q0, q1, qz0)
absmmax=k_evolve((mmax-dmK), zz, q0, q1, qz0)
absmmin=(absmmin > sample_absmmin) < sample_absmmax
absmmax=(absmmax > sample_absmmin) < sample_absmmax
absmag_evol=k_evolve(absmag, zz, q0, q1, qz0)

soname=filepath('libkcorrect.'+kcorrect_so_ext(), $
                root_dir=getenv('KCORRECT_DIR'), subdirectory='lib')
retval=call_external(soname, 'idl_lf_select_eep', float(zz), $
                     float(absmag_evol), float(absmmin),float(absmmax), $
                     long(ngals), float(sample_absmmin),float(sample_absmmax),$
                     float(absmk),float(phi),float(sel),long(nbin))

end
