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
pro k_coeffdist_plot,version,vpath=vpath,basecoeff=basecoeff,subsample=subsample,nsig=nsig,psfile=psfile

if(n_elements(vpath) eq 0) then $
  vpath=getenv('KCORRECT_DIR')+'/data/etemplates'
if(n_elements(lamlim) eq 0) then lamlim=[2000.,10000.]
if(n_elements(nsig) eq 0) then nsig=3.

k_load_ascii_table,coeff,vpath+'/coeff.'+version+'.dat'
k_load_ascii_table,ematrix,vpath+'/ematrix.'+version+'.dat'
k_load_ascii_table,bmatrix,vpath+'/bmatrix.'+version+'.dat'
k_load_ascii_table,lambda,vpath+'/lambda.'+version+'.dat'
nt=long((size(ematrix))[2])
if(n_elements(basecoeff) eq 0) then basecoeff=nt-1
ngalaxy=long(n_elements(coeff))/nt
nl=n_elements(lambda)-1l

coeffindx=(lindgen(nt-1)+1)[where(lindgen(nt-1)+1 ne basecoeff)]
coeffmean=dblarr(nt-2)
for i=0, nt-3 do $
  coeffmean[i]=djs_avsigclip(coeff[coeffindx[i],*]/coeff[0,*])
speccoeffs=dblarr(nt,3)
speccoeffs[0,*]=1.
speccoeffs[basecoeff,*]=[-0.20,-0.06,0.10]
for i=0, nt-3 do $
  speccoeffs[coeffindx[i],*]=coeffmean[i]
nspecs=(size(speccoeffs))[2]

if(n_elements(subsample) eq 0) then subsample=1l
indx=lindgen(ngalaxy/long(subsample))*long(subsample)
coeff=coeff[*,indx]

if(keyword_set(psfile)) then begin
    set_plot, "PS"
    xsize= 7.5 & ysize= 7.5
    device, file=psfile,/inches,xsize=xsize,ysize=ysize, $
      xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color,/encapsulated
    !P.FONT= -1 & !P.BACKGROUND= 255 & !P.COLOR= 0
endif else begin
    set_plot,'x'
endelse
!P.THICK= 2.0
!P.CHARTHICK= !P.THICK & !X.THICK= !P.THICK & !Y.THICK= !P.THICK
!P.CHARSIZE= 1.2
axis_char_scale= 1.0
tiny= 1.d-4
!P.PSYM= 0
!P.TITLE= ''
!X.STYLE= 1
!X.CHARSIZE= axis_char_scale
!X.MARGIN= [1,1]*0.5*axis_char_scale
!X.OMARGIN= [6,6]*axis_char_scale
!X.RANGE= 0
!Y.STYLE= 1
!Y.CHARSIZE= !X.CHARSIZE
!Y.MARGIN= 0.6*!X.MARGIN
!Y.OMARGIN= 0.6*!X.OMARGIN
!Y.RANGE= 0
xyouts, 0,0,'!3'


; Make a vector of 16 points, A[i] = 2pi/16:
A = FINDGEN(17) * (!PI*2/16.)
; Define the symbol to be a unit circle with 16 points, 
; and set the filled flag:
USERSYM, COS(A), SIN(A), /FILL

; make useful vectors for plotting
colorname= ['red','green','blue','magenta','cyan','dark yellow', $
            'purple','light green','orange','navy','light magenta', $
            'yellow green']
ncolor= n_elements(colorname)

!p.multi=[2*(nt-2),nt-2,2]
espec=bmatrix#ematrix
lam=0.5*(lambda[0l:nl-1l]+lambda[1l:nl])
for i=0, n_elements(coeffindx)-1 do begin
    klog,basecoeff,coeffindx[i]
    xrat=coeff[basecoeff,*]/coeff[0,*]
    yrat=coeff[coeffindx[i],*]/coeff[0,*]
    xsig=djsig(xrat)
    ysig=djsig(yrat)
    xmean=mean(xrat)
    ymean=mean(yrat)
    !X.CHARSIZE= tiny
    !Y.CHARSIZE= tiny
    !X.RANGE= xmean+nsig*[-0.7*xsig,1.3*xsig]
    !Y.RANGE= ymean+nsig*[-xsig,xsig]
    djs_plot,xrat,yrat,psym=3,xst=1,yst=1
    axis,!X.RANGE[0],!Y.RANGE[1],xaxis=1,xcharsize=axis_char_scale, $
      xtitle='a!d'+strtrim(string(basecoeff),2)+'!n/a!d0'
    if(i eq 0) then $
      axis,!X.RANGE[0],!Y.RANGE[0],yaxis=0,ycharsize=axis_char_scale, $
      ytitle='a!d'+strtrim(string(coeffindx[i]),2)+'!n/c!d0'
    if(i eq n_elements(coeffindx)-1) then $
      axis,!X.RANGE[1],!Y.RANGE[0],yaxis=1,ycharsize=axis_char_scale, $
      ytitle='a!d'+strtrim(string(coeffindx[i]),2)+'!n/c!d0'
    for j=0,nspecs-1 do begin
        djs_oplot,[speccoeffs[basecoeff,j]],[speccoeffs[coeffindx[i],j]],psym=8, $
          color=colorname[j],/fill
    endfor
endfor

!X.CHARSIZE= axis_char_scale
!Y.CHARSIZE= axis_char_scale
!p.multi=[1,1,2]
spec=espec#speccoeffs
out=spec[*,0]/max(spec[*,0])
djs_plot,lam,out,xst=1,yst=1,yra=[-0.1,1.1],xra=[2000.,11000.],/nodata, $
  xtitle='Wavelength (Angstroms)', ytitle='f!d!4k!3', $
  ycharsize=axis_char_scale, xcharsize=axis_char_scale
xout=lam[n_elements(lam)/9]
for i = 0, nt-3 do begin
    yout=0.95-double(i)*0.1
    xyouts,xout,yout,'a!d'+strtrim(string(coeffindx[i]),2)+'!n/a!d0!n='+ $
      strtrim(string(speccoeffs[coeffindx[i],0],format='(f5.2)'),2)
endfor
scalespec=[1.,2.,1.5]
for j=0, nspecs-1 do begin
    out=scalespec[j]*spec[*,j]/max(spec[*,0])* $
       interpolate(spec[*,0], $
									 double(nl)*(3000.-lambda[0])/(lambda[nl]-lambda[0]))/ $
       interpolate(spec[*,j], $
									 double(nl)*(3000.-lambda[0])/(lambda[nl]-lambda[0]))
    djs_oplot,lam,out,color=colorname[j],thick=6
    xout=lam[n_elements(lam)/2]
    yout=out[n_elements(lam)/2]+0.02
    xyouts,xout,yout,'a!d'+strtrim(string(basecoeff),2)+'!n/a!d0!n='+ $
      strtrim(string(speccoeffs[basecoeff,j],format='(f5.2)'),2)
endfor

if(keyword_set(psfile)) then begin
    device, /close
    set_plot,'x'
endif
    
end
;------------------------------------------------------------------------------
