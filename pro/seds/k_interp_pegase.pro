;+
; NAME:
;   k_interp_pegase
; PURPOSE:
;   interpolate a set of PEGASE spectra to a specific time
; CALLING SEQUENCE:
;   k_interp_pegase, peg, time [,nl=, lmin=, lmax=, /nolines, $
;      /nocontinuum]
; INPUTS:
;   peg - structure returned by read_peg
;   time - desired time (in Myrs)
; OPTIONAL INPUTS:
;   nl, lmin, lmax - structure of output spectrum
; KEYWORDS:
;   /nolines - leave out the lines
;   /nocontinuum - leave out the continuum
; REVISION HISTORY:
;   25-Jul-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function k_interp_pegase,peg,time,nl=nl,lmin=lmin,lmax=lmax,nolines=nolines, $
                         nocontinuum=nocontinuum

; find bounds
ii=lindgen(n_elements(peg)-1L)
mintime=peg[0].arr1[0]
maxtime=peg[n_elements(peg)-1L].arr1[0]- $
  1.e-6*(peg[n_elements(peg)-1L].arr1[0]-peg[n_elements(peg)-2L].arr1[0])
use_time=(time > mintime) < maxtime
itime=(where(use_time ge peg[ii].arr1[0] and use_time lt peg[ii+1].arr1[0]))[0]
itime=((itime > 0) < (n_elements(peg)-2L))
klog,'itime= '+string(itime)

; make spectra
k_spec_pegase,peg[itime+0],spec0,lambda,nl=nl,lmin=lmin,lmax=lmax, $
  nolines=nolines,nocontinuum=nocontinuum
k_spec_pegase,peg[itime+1],spec1,lambda,nl=nl,lmin=lmin,lmax=lmax, $
  nolines=nolines,nocontinuum=nocontinuum

; interpolate
stime=(use_time-peg[itime].arr1[0])/(peg[itime+1].arr1[0]-peg[itime].arr1[0])
spec=spec0+stime*(spec1-spec0)

return,spec

end
;------------------------------------------------------------------------------
