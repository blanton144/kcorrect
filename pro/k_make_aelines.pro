;+
; NAME:
;   k_make_aelines
; PURPOSE:
;   Make the 3-template A+E+lines
; CALLING SEQUENCE:
; INPUTS:
; OPTIONAL INPUTS:
; KEYWORDS:
; OUTPUTS:
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   24-jan-2002  WRitten by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_make_aelines, vmatrix, lambda,nl=nl,lmin=lmin,lmax=lmax

if(NOT keyword_set(pegasepath)) then $
  pegasepath=getenv('DATA')+'/specmodels/PEGASE.2'
if(NOT keyword_set(nl)) then nl=5000L
if(NOT keyword_set(lmin)) then lmin=2500.
if(NOT keyword_set(lmax)) then lmax=12000.

vmatrix=dblarr(nl,2)
lambda=lmin+(lmax-lmin)*dindgen(nl+1)/double(nl)
dl=lambda[1]-lambda[0]

pegfile=pegasepath+'/mrb_spectra.0.008.dat'  
read_peg,pegfile,peg=peg
vmatrix[*,0]=k_interp_pegase(peg,12000.,nl=nl,lmin=lmin,lmax=lmax)
vmatrix[*,1]=k_interp_pegase(peg,100.,nl=nl,lmin=lmin,lmax=lmax)

dust1={dusty_str, geometry:'', dust:'', structure:'', tauv:0.}
dust=replicate(dust1,1)
dust.geometry=['dusty']
dust.dust=['MW']
dust.structure=['c']
dust.tauv=[1.5]
vmatrix[*,0]=vmatrix[*,0]*exp(-witt_ext(dust[0],dust[0].tauv,lambda[0:nl-1]))

nv=n_elements(vmatrix)/(n_elements(lambda)-1L)
for i=0,nv-1 do $
  vmatrix[*,i]=vmatrix[*,i]/total(vmatrix[*,i],/double)

k_write_ascii_table,vmatrix,'vmatrix.aelines.dat'
k_write_ascii_table,lambda,'lambda.aelines.dat'

end
;------------------------------------------------------------------------------
