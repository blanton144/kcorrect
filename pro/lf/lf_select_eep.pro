;+
; NAME:
;   lf_select_eep
; PURPOSE:
;   calculate EEP selection function from luminosity function
; USAGE:
;   lf_select_eep,zz,absmag,kcorrect,mmin,mmax,sample_absmmin,sample_absmmax,
;       sample_zmin,sample_zmax, absmk, phi, sel [, qevolve=, $
;       omega0=, omegal0=, qz0=, absmagdep=, ref_absmagdep= ]
; INPUTS:
;   zz               [N] redshifts
;   absmag           [N] absolute magnitudes
;   kcorrect         [N] K-correction for each galaxy
;   mmin             [N] minimum apparent mag for each galaxy
;   mmax             [N] maximum apparent mag for each galaxy
;   sample_absmmin   observed absolute magnitude minimum of sample
;   sample_absmmax   observed absolute magnitude maximum of sample
; OPTIONAL INPUTS:
;   omega0           omega_matter to use (default: 0.3)
;   omegal0          omega_lambda to use (default: 0.7)
;   qevolve          evolution
;   qz0              pivot for redshift evolution
;   absmk            center of each bin
;   phi              amplitude of each bin
; REVISION HISTORY:
;   2003-10-20  written - Blanton
;-
pro lf_select_eep,zz,absmag,kcorrect,mmin,mmax,sample_absmmin,sample_absmmax, $
                  absmk,phi,sel, omega0=omega0,omegal0=omegal0
                 
; set defaults
if(n_elements(omega0) eq 0) then omega0=0.3
if(n_elements(omegal0) eq 0) then omegal0=0.7

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

soname=filepath('libkcorrect.'+idlutils_so_ext(), $
                root_dir=getenv('KCORRECT_DIR'), subdirectory='lib')
retval=call_external(soname, 'idl_lf_select_eep', float(zz), $
                     float(absmmin),float(absmmax), long(ngals), $
                     float(sample_absmmin),float(sample_absmmax), $
                     float(absmk),float(phi),float(sel),long(nbin))

end
