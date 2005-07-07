;+
; NAME:
;   lf_eep
; PURPOSE:
;   calculate EEP luminosity function 
; USAGE:
;   lf_eep,zz,absmag,kcorrect,mmin,mmax,sample_absmmin,sample_absmmax,
;       absmk, phi, phi_err, [ nbin=, q0=, q1=, qz0=, $
;       omega0=, omegal0=, /calc_err ]
; INPUTS:
;   zz               [N] redshifts
;   absmag           [N] absolute magnitudes (kcorrected, unevolved)
;   kcorrect         [N] K-correction for each galaxy
;   mmin             [N] minimum apparent mag for each galaxy
;   mmax             [N] maximum apparent mag for each galaxy
;   sample_absmmin   observed absolute magnitude minimum of sample
;   sample_absmmax   observed absolute magnitude maximum of sample
; OPTIONAL INPUTS:
;   q0, q1, qz0      evolution model (pass to k_evolve())
;   omega0           omega_matter to use (default: 0.3)
;   omegal0          omega_lambda to use (default: 0.7)
; KEYWORDS:
;   /calc_err        if set, calculate the error bars
; OUTPUTS:
;   absmk            center of each bin
;   phi              amplitude of each bin
;   phi_err          error in each bin
; REVISION HISTORY:
;   2003-10-20  written - Blanton
;-
pro lf_eep,zz,absmag,kcorrect,mmin,mmax,sample_absmmin,sample_absmmax, $
           absmk,phi,phi_err,nbin=nbin,q0=q0,q1=q1, qz0=qz0, $
           omega0=omega0, omegal0=omegal0, calc_err=calc_err, weight=weight

; set defaults
if(n_elements(nbin) eq 0) then nbin=40L
if(n_elements(omega0) eq 0) then omega0=0.3
if(n_elements(omegal0) eq 0) then omegal0=0.7
if(n_elements(q0) eq 0) then q0=0.0
if(n_elements(q1) eq 0) then q1=0.0
if(n_elements(qz0) eq 0) then qz0=0.1
if(n_elements(calc_err) eq 0) then calc_err=0L
if(n_elements(weight) eq 0) then weight=fltarr(n_elements(zz))+1.

; settings
absmk=fltarr(nbin+1L)
phi=fltarr(nbin)
phi_err=fltarr(nbin)
covar=fltarr(nbin,nbin)
ngals=n_elements(zz)

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
retval=call_external(soname, 'idl_lf_eep', float(zz), $
                     float(absmag_evol),float(absmmin),float(absmmax), $
                     long(ngals),float(sample_absmmin),float(sample_absmmax), $
                     float(absmk),float(phi),float(phi_err), $
                     float(covar),long(nbin),long(calc_err), $
                     float(weight))

end
