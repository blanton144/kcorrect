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
pro k_coeff_plot,version,vpath=vpath, nlum=nlum, omega0=omega0, omegal0=omegal0, lumlimits=lumlimits, psfile=psfile, nsig=nsig

if(n_elements(vpath) eq 0) then $
  vpath=getenv('KCORRECT_DIR')+'/data/etemplates'
if(n_elements(nlum) eq 0) then nlum=1l
if(n_elements(omega0) eq 0) then omega0=0.3d
if(n_elements(omegal0) eq 0) then omegal0=0.7d
if(n_elements(nsig) eq 0) then nsig=3.
if(n_elements(lumlimits) eq 0) then $ 
  lumlimits=[[-22.5,-22.]]
help,nlum

k_load_ascii_table,coeff,vpath+'/coeff.'+version+'.dat'
nt=(size(coeff))[1]
k_load_ascii_table,z,vpath+'/z.'+version+'.dat'
k_reconstruct_maggies,coeff,replicate(0.2,n_elements(z)),reconstruct_maggies, $
  version=version
dm=2.5*alog10((2.99792e+8*lumdis(z,omega0,omegal0))^2)
lum=-2.5*alog10(reconstruct_maggies[2,*])-dm

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
axis_char_scale= 2.6
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

!p.multi=[nlum*(nt-1l),nlum,nt-1l]
for i=0, nt-2l do begin
    for j=0, nlum-1l do begin
        indx=where((lum) gt lumlimits[0,j] and $
                   (lum) lt lumlimits[1,j],count)
        rat=coeff[i+1l,indx]/coeff[0,indx]
        sig=djsig(rat)
        meanr=djs_avsigclip(rat)
        !X.RANGE=[0.,0.5]
        !Y.RANGE=[meanr-nsig*sig,meanr+nsig*sig]
        !X.TITLE=''
        if(i eq nt-2) then begin
            !X.CHARSIZE= axis_char_scale
            !X.TITLE='Redshift z'
        endif else begin
            !X.CHARSIZE= tiny
        endelse
        djs_plot,z[indx],rat,psym=3
        isortz=indx[sort(z[indx])]
        sortz=z[isortz]
        sortc=transpose(coeff[i+1l,isortz]/coeff[0,isortz])
        result=regress(sortz,sortc,yfit=fitc)
        djs_oplot,[0.,0.5],[meanr,meanr],thick=6
        xyouts, $
          !X.RANGE[0]+0.05*(!X.RANGE[1]-!X.RANGE[0]), $
          !Y.RANGE[1]-0.2*(!Y.RANGE[1]-!Y.RANGE[0]), $
          'a!d'+strtrim(string(i+1l),2)+'!n/a!d0', $
          charsize=0.85*axis_char_scale
        if(i eq 0) then $
          xyouts, $
          !X.RANGE[0]+0.05*(!X.RANGE[1]-!X.RANGE[0]), $
          !Y.RANGE[1]+0.06*(!Y.RANGE[1]-!Y.RANGE[0]), $
          strtrim(string(lumlimits[0],format='(f8.2)'),2) + $
          '<M!d!s!e0.2!r r!n<'+ $
          strtrim(string(lumlimits[1],format='(f8.2)'),2), $
          charsize=1.*axis_char_scale 
    endfor
endfor

if(keyword_set(psfile)) then begin
    device, /close
    set_plot,'x'
endif
    
end
;------------------------------------------------------------------------------
