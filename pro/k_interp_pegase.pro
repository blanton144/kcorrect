;+
; NAME:
;   k_interp_pegase
; PURPOSE:
;   Interpolate a set of PEGASE spectra to a specific time
; CALLING SEQUENCE:
;   k_interp_pegase, peg, time [,nl=, lmin=, lmax=]
; INPUTS:
;   peg - structure returned by read_peg
;   time - desired time (in Myrs)
; OPTIONAL INPUTS:
;   nl, lmin, lmax - structure of output spectrum
; KEYWORDS:
; OUTPUTS:
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   25-Jul-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function k_interp_pegase,peg,time,nl=nl,lmin=lmin,lmax=lmax

; find bounds
ii=lindgen(n_elements(peg)-1L)
if(time lt peg[0].arr1[0]) then $
  itime=0 $
else if (time ge peg[n_elements(peg)-2L].arr1[0]) then $
  itime=n_elements(peg)-2L $
else $
  itime=(where(time ge peg[ii].arr1[0] and time lt peg[ii+1].arr1[0]))[0]
klog,'itime= '+string(itime)

; make spectra
k_spec_pegase,peg[itime+0],spec0,lambda,nl=nl,lmin=lmin,lmax=lmax
k_spec_pegase,peg[itime+1],spec1,lambda,nl=nl,lmin=lmin,lmax=lmax

; interpolate
stime=(time-peg[itime].arr1[0])/(peg[itime+1].arr1[0]-peg[itime].arr1[0])
spec=spec0+stime*(spec1-spec0)

return,spec

end
;------------------------------------------------------------------------------
