;+
; NAME:
;   k_make_v3_0_templates
; PURPOSE:
;   Code to make the templates for what I intend to be v3_0 of
;   kcorrect. 
; CALLING SEQUENCE:
;   k_make_v3_0_templates
; INPUTS:
; OPTIONAL INPUTS:
; KEYWORDS:
; OUTPUTS:
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
;   k_mkspec_pegase
; REVISION HISTORY:
;   18-Jan-2003  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_make_v3_0_templates

if(n_elements(nages) eq 0) then nages=12
if(n_elements(sigscale) eq 0) then sigscale=0.1

; settings
pi=!DPI

; make sfh pars
attime=0.D
minage=1.5D+6
maxage=2.6D+10
sfhpars1={sfhstr, age:0.d, agesigma:0.d}
sfhpars=replicate(sfhpars1,nages)
sfhpars.age=exp((alog(minage)+(alog(maxage)-alog(minage))* $
                 ((dindgen(nages)+0.5)/double(nages))))
sfhpars[0].agesigma=max([sigscale*(sfhpars[1].age-sfhpars[0].age),1.d+6])
for i=1, nages-1 do begin
    sfhpars[i].agesigma=max([sigscale*(sfhpars[i].age-sfhpars[i-1].age),1.d+6])
endfor

; make dust
dust1={dusty_str, geometry:'', dust:'', structure:'', tauv:0.}
dust=replicate(dust1,3)
dust.geometry=['dusty','dusty','dusty']
dust.dust=['MW','MW','MW']
dust.structure=['c','c','c']
dust.tauv=[0.,3.5,8.]

; make metallicity
metallicity=['02','004','001']

; make galaxy templates 
k_mkspec_pegase, galvmatrix, lambda, metallicity, dust, $
  'gaussian', sfhpars, minage=minage, maxage=maxage, attime=attime, $
  lmin=lmin, lmax=lmax, nl=nl,nolines=nolines

; make a qso template
k_mkspec_qso, qsovmatrix, lambda

; put together vmatrices
ngal=n_elements(galvmatrix)/(n_elements(lambda)-1L)
nqso=n_elements(qsovmatrix)/(n_elements(lambda)-1L)
ntotal=ngal+nqso
vmatrix=dblarr(n_elements(lambda)-1L,ntotal)
vmatrix[*,0:ngal-1]=galvmatrix
vmatrix[*,ngal:ntotal-1L]=qsovmatrix

; output the appropriate files
k_write_ascii_table,vmatrix,'vmatrix.test.dat'
k_write_ascii_table,lambda,'lambda.test.dat'

end
