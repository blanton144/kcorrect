;+
; NAME:
;   k_make_templates
; PURPOSE:
;   Code to make kcorrect templates 
; CALLING SEQUENCE:
;   k_make_templates [, nages=, sigscale=, minage=, maxage= ]
; INPUTS:
; OPTIONAL INPUTS:
;   nages - number of ages
;   sigscale - width of distribution
;   [min|max]age - age limits
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
pro k_make_templates, nages=nages, sigscale=sigscale, maxage=maxage, $
                      minage=minage, name=name

if(n_elements(name) eq 0) then name='test'
if(n_elements(nages) eq 0) then nages=12
if(n_elements(sigscale) eq 0) then sigscale=0.1
if(n_elements(minage) eq 0) then minage=1.5e+6
if(n_elements(maxage) eq 0) then maxage=2.6e+10

; make sfh pars
attime=0.D
sfhpars1={sfhstr, age:0.d, agesigma:0.d}
sfhpars=replicate(sfhpars1,nages)
sfhpars.age=exp((alog(minage)+(alog(maxage)-alog(minage))* $
                 ((findgen(nages)+0.5)/float(nages))))
sfhpars[0].agesigma=max([sigscale*(sfhpars[1].age-sfhpars[0].age),1.e+6])
for i=1, nages-1 do begin
    sfhpars[i].agesigma=max([sigscale*(sfhpars[i].age-sfhpars[i-1].age),1.e+6])
endfor

; make dust
dust1={dusty_str, geometry:'', dust:'', structure:'', tauv:0.}
dust=replicate(dust1,3)
dust.geometry=['dusty','dusty','dusty']
dust.dust=['MW','MW','MW']
dust.structure=['c','c','c']
dust.tauv=[0.,3.5,8.]

; make metallicity
metallicity=['03','02','008','004','001']

; make galaxy templates 
k_mkspec_pegase, galvmatrix, lambda, metallicity, dust, $
  sfhpars, minage=minage, maxage=maxage, attime=attime, $
  lmin=lmin, lmax=lmax, nl=nl,nolines=nolines

; put together vmatrices
ngal=n_elements(galvmatrix)/(n_elements(lambda)-1L)
vmatrix=fltarr(n_elements(lambda)-1L,ngal)
vmatrix[*,0:ngal-1]=galvmatrix

; output the appropriate files
k_write_ascii_table,vmatrix,'vmatrix.'+name+'.dat'
k_write_ascii_table,lambda,'lambda.'+name+'.dat'

end
