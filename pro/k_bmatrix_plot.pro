;+
; NAME:
;   k_espec_plot
;
; PURPOSE:
;   Make plot of eigenspectra for K-correct paper
;
; CALLING SEQUENCE:
;   k_espec_plot,version,[vpath=]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL INPUT/OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   Change to using main sample only, and cutting on M_r; cannot
;   determine much from color-selected samples...
;
; PROCEDURES CALLED:
;   k_load_ascii_table
;
; REVISION HISTORY:
;   23-Jan-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_bmatrix_plot,version,vpath=vpath,lamlim=lamlim

if(n_elements(vpath) eq 0) then $
  vpath=getenv('KCORRECT_DIR')+'/data/etemplates'
if(n_elements(lamlim) eq 0) then lamlim=[2000.,10000.]

k_load_ascii_table,bmatrix,vpath+'/bmatrix.'+version+'.dat'
k_load_ascii_table,lambda,vpath+'/lambda.'+version+'.dat'
nl=n_elements(lambda)-1l
nb=long(n_elements(bmatrix))/nl

lam=0.5*(lambda[0l:nl-1l]+lambda[1l:nl])
!p.multi=[nb,1,nb]
for i=0, nb-1l do begin
    print,i
    plot,lam,bmatrix[*,i],xst=1,yst=1,xra=lamlim
endfor
    
end
;------------------------------------------------------------------------------
