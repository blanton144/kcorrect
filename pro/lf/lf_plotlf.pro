;+
; NAME:
;   lf_plotlf
; PURPOSE:
;   plot a luminosity function
; INPUTS:
;   nphi             number of "bins"
;   sample_absmmin   minimum abs mag for the sample
;   sample_absmmax   maximum abs mag for the sample
;   sigabsmag        width of gaussian
;   phi              amplitude of bin 
; OPTIONAL INPUTS:
;   subsample        how much to subsample for plot
; KEYWORDS:
; OUTPUTS:
; OPTIONAL OUTPUTS:
; BUGS:
; DEPENDENCIES:
; REVISION HISTORY:
;   2002-5-22  written - Blanton
;-
pro lf_plotlf,nphi,absmk,sample_absmmin,sample_absmmax,sigabsmag,phi, $
              data_absmag,subsample=subsample

if(NOT keyword_set(subsample)) then subsample=10L
pi=3.14159265358979

absmagvals=sample_absmmin+(dindgen(subsample*nphi)+0.5)/ $
  double(subsample*nphi)*(sample_absmmax-sample_absmmin)
phivals=dblarr(subsample*nphi)

factor=1.d/(sqrt(2*pi)*sigabsmag)
invsig=0.5d/(sigabsmag^2)
for i=0L, nphi-1L do begin
    phivals=phivals+phi[i]*factor*exp(-invsig*(absmk[i]-absmagvals)^2)
endfor

range=max(alog10(phivals))-min(alog10(phivals))
;plothist,data_absmag,bin=0.4*sigabsmag
;phivals=phivals*double(n_elements(data_absmag))*(0.4*sigabsmag)
plot,absmagvals,alog10(phivals),xst=1,yst=1,yra=[min(alog10(phivals))-0.03*range, $
                                                 max(alog10(phivals))+0.03*range]
modvals=exp(-0.5*((absmagvals+21.5)/0.4)^2)
modvals=lf_schechter(absmagvals,1.,-20.5,-1.)
modvals=total(phivals)*modvals/(total(modvals))
oplot,absmagvals,alog10(modvals),color=255

for i=0L, nphi-1L do begin
    phivals=phi[i]*factor*exp(-invsig*(absmk[i]-absmagvals)^2)
    oplot,absmagvals,alog10(phivals),color=190
endfor

end
