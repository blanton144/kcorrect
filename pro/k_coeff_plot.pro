;+
; NAME:
;   k_coeff_plot
;
; PURPOSE:
;   Make plot of coefficients for K-correct paper
;
; CALLING SEQUENCE:
;   k_coeff_plot,version,[vpath=]
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
pro k_coeff_plot,version,vpath=vpath, nlum=nlum, omega0=omega0, omegal0=omegal0, lumlimits=lumlimits

if(n_elements(vpath) eq 0) then $
  vpath=getenv('KCORRECT_DIR')+'/data/etemplates'
if(n_elements(nlum) eq 0) then nlum=1l
if(n_elements(omega0) eq 0) then omega0=0.3d
if(n_elements(omegal0) eq 0) then omegal0=0.7d
if(n_elements(lumlimits) eq 0) then $ 
  lumlimits=[[-5.2,-4.8]]
help,nlum

k_load_ascii_table,coeff,vpath+'/coeff.'+version+'.dat'
nt=(size(coeff))[1]
k_load_ascii_table,z,vpath+'/z.'+version+'.dat'
;k_model_fluxes,coeff,replicate(0.1,n_elements(z)),model_fluxes, $
  ;version=version
;lum=lumdis(z,omega0,omegal0)^2*model_fluxes[2,*]
lum=lumdis(z,omega0,omegal0)^2*coeff[0,*]

!p.multi=[nlum*(nt-1l),nlum,nt-1l]
for i=0, nt-2l do begin
    for j=0, nlum-1l do begin
        indx=where(alog10(lum) gt lumlimits[0,j] and $
                   alog10(lum) lt lumlimits[1,j],count)
        rat=coeff[i+1l,indx]/coeff[0,indx]
        sig=djsig(rat)
        meanr=mean(rat)
        if(count gt 0) then $
          plot,z[indx],rat,psym=3,xra=[0.,0.5],yra=[meanr-15.*sig,meanr+15.*sig]
        isortz=indx[sort(z[indx])]
        sortz=z[isortz]
        sortc=transpose(coeff[i+1l,isortz]/coeff[0,isortz])
        result=regress(sortz,sortc,yfit=fitc)
        oplot,sortz,fitc
    endfor
endfor
    
end
;------------------------------------------------------------------------------
