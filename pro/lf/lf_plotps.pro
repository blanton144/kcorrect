;+
; NAME:
;   lf_plotps
; PURPOSE:
;   plot a luminosity function
; INPUTS:
;   absmk, phi, sigabsmag - description of npgauss l.f.
;   sample_absmmin - minimum abs mag for the sample
;   sample_absmmax - maximum abs mag for the sample
; OPTIONAL INPUTS:
;   subsample - sample the function at this sampling (default 10)
;   filename - output file name (default 'npgauss_lf.ps')
;   schecter - structure with .PHISTAR, .MSTAR, .ALPHA to show for
;              comparison (default none)
; OPTIONAL KEYWORDS:
;   /noclose - do not close file
; OPTIONAL OUTPUTS:
;   pold, xold, yold - saved !P, !X, !Y before any changes
;   phivals, amvals - l.f. values at the output abs. mags. 
; REVISION HISTORY:
;   2002-5-22  written - Blanton
;-
pro lf_plotps,absmk,phi,sigabsmag,sample_absmmin,sample_absmmax, $
              subsample=subsample, filename=filename, phivals=phivals, $
              amvals=amvals, schechter=schechter, noclose=noclose, $
              pold=pold, xold=xold, yold=yold

; defaults
if(NOT keyword_set(subsample)) then subsample=10L
if(NOT keyword_set(filename)) then filename='npgauss_lf.ps'

; settings
pi=3.14159265358979D+0
nphi=n_elements(phi)

; set amvals and phivals
lf_vals,absmk,phi,sigabsmag,sample_absmmin,sample_absmmax,amvals,phivals, $
  subsample=subsample,factor=factor,invsig=invsig

; if defined plot schechter
if(keyword_set(schechter)) then begin
    splog,'trying schechter'
    if(schechter.phistar eq 0.) then begin
        schvals=0.4*alog(10.)*1.* $
          10.^(-0.4*(amvals-schechter.mstar)*(schechter.alpha+1.))* $
          exp(-10.^(-0.4*(amvals-schechter.mstar)))
        total_schvals=total(schvals,/double)*(amvals[1]-amvals[0])
        schvals=schvals/total_schvals
    endif else begin
        schvals=0.4*alog(10.)*schechter.phistar* $
          10.^(-0.4*(amvals-schechter.mstar)*(schechter.alpha+1.))* $
          exp(-10.^(-0.4*(amvals-schechter.mstar)))
    endelse 
endif

; setup postscript file
pold=!P
xold=!X
yold=!Y
!P.FONT= -1
set_plot, "PS"
!P.BACKGROUND= djs_icolor('white')
!P.COLOR= djs_icolor('black')
xsize= 7.5 & ysize= 7.5
device, file=filename,/inches,xsize=xsize,ysize=ysize, $
  xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color
!P.THICK= 2.0
!P.CHARTHICK= !P.THICK & !X.THICK= !P.THICK & !Y.THICK= !P.THICK
!P.CHARSIZE= 1.0
if(NOT keyword_set(axis_char_scale)) then axis_char_scale= 1.75
tiny= 1.d-4
!P.PSYM= 0
!P.LINESTYLE= 0
!P.TITLE= ''
!X.STYLE= 1
!X.CHARSIZE= axis_char_scale
!X.MARGIN= [1,1]*0.5
!X.OMARGIN= [5,5]*axis_char_scale
!X.RANGE= 0
!X.TICKS= 0
!Y.STYLE= 1
!Y.CHARSIZE= !X.CHARSIZE
!Y.MARGIN= 0.6*!X.MARGIN
!Y.OMARGIN= 0.6*!X.OMARGIN
!Y.RANGE= 0
!Y.TICKS= !X.TICKS
!P.MULTI= [1,1,1]
xyouts, 0,0,'!3'
colorname= ['red','green','blue','magenta','cyan','dark yellow', $
            'purple','light green','orange','navy','light magenta', $
            'yellow green']
ncolor= n_elements(colorname)
theta= 2.0D *double(!PI)*dindgen(31)/30.0D
x= cos(theta)
y= sin(theta)

minphi=max([min(alog10(phivals)),-8.])
range=max(alog10(phivals))-minphi
plot,amvals,alog10(phivals),xst=1,yst=1, $
  yra=[minphi-0.03*range,max(alog10(phivals))+0.03*range], $
  thick=8

for i=0L, nphi-1L do begin
    phivals_single=phi[i]*factor*exp(-invsig*(absmk[i]-amvals)^2)
    oplot,amvals,alog10(phivals_single),thick=1
endfor

if(keyword_set(schvals)) then $
  djs_oplot,amvals,alog10(schvals),color='red',thick=3

if(keyword_set(noclose)) then return

device,/close

!P=pold
!X=xold
!Y=yold


end
