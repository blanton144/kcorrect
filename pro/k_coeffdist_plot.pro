;+
; NAME:
;   k_coeffdist_plot
;
; PURPOSE:
;   Make plot of distribution of coefficients for K-correct paper
;
; CALLING SEQUENCE:
;   k_coeffdist_plot,version,[vpath=]
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
pro k_coeffdist_plot,version,vpath=vpath,basecoeff=basecoeff

if(n_elements(vpath) eq 0) then $
  vpath=getenv('KCORRECT_DIR')+'/data/etemplates'
if(n_elements(lamlim) eq 0) then lamlim=[2000.,10000.]
if(n_elements(basecoeff) eq 0) then basecoeff=2

k_load_ascii_table,coeff,vpath+'/coeff.'+version+'.dat'
k_load_ascii_table,ematrix,vpath+'/ematrix.'+version+'.dat'
k_load_ascii_table,bmatrix,vpath+'/bmatrix.'+version+'.dat'
k_load_ascii_table,lambda,vpath+'/lambda.'+version+'.dat'
nt=(size(ematrix))[2]
nl=n_elements(lambda)-1l

espec=bmatrix#ematrix
lam=0.5*(lambda[0l:nl-1l]+lambda[1l:nl])
!p.multi=[2*(nt-2),nt-2,2]
indx=(lindgen(nt-1)+1)[where(lindgen(nt-1)+1 ne basecoeff)]
for i=0, n_elements(indx)-1 do begin
    xrat=coeff[basecoeff,*]/coeff[0,*]
    yrat=coeff[indx[i],*]/coeff[0,*]
    xsig=djsig(xrat)
    ysig=djsig(yrat)
    xmean=mean(xrat)
    ymean=mean(yrat)
    plot,xrat,yrat,psym=3,xst=1,yst=1,xra=xmean+[-2.*xsig,4.*xsig], $
      yra=ymean+3.*[-xsig,xsig]
endfor

!p.multi=[2,1,2]
plot,xrat,yrat,psym=3,xst=1,yst=1,xra=xmean+[-2.*xsig,4.*xsig], $
  yra=ymean+3.*[-xsig,xsig]
    
end
;------------------------------------------------------------------------------
