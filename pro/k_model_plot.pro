;+
; NAME:
;   k_model_plot
;
; PURPOSE:
;   Make plot of model versus actual fluxes
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
pro k_model_plot,savfile,version=version,vpath=vpath,psfile=psfile,nsig=nsig, subsample=subsample, sdssfix=sdssfix, vconstraint=vconstraint, addgrgap=addgrgap

if(NOT keyword_set(version)) then version='default'
if(NOT keyword_set(vpath)) then vpath=getenv('KCORRECT_DIR')+'/data/etemplates'
if(NOT keyword_set(nsig)) then nsig=4.d

restore,savfile
if(n_elements(sp) gt 0) then galaxy_z=sp.z
if(n_elements(z) gt 0) then galaxy_z=z
ngalaxy=long(n_elements(galaxy_z))
nk=long(n_elements(galaxy_maggies))/ngalaxy

;if(n_elements(sp) gt 0) then begin
    ;indx=where(sp.petrocountserr[0,*] lt 0.1 and $
               ;sp.petrocountserr[1,*] lt 0.1 and $
               ;sp.petrocountserr[2,*] lt 0.1 and $
               ;sp.petrocountserr[3,*] lt 0.1 and $
               ;sp.petrocountserr[4,*] lt 0.1)
    ;galaxy_z=galaxy_z[indx]
    ;galaxy_maggies=galaxy_maggies[*,indx]
    ;coeff=coeff[*,indx]
;endif

if(n_elements(subsample) eq 0) then subsample=1l
indx=lindgen(ngalaxy/long(subsample))*long(subsample)
galaxy_z=galaxy_z[indx]
galaxy_maggies=galaxy_maggies[*,indx]
coeff=coeff[*,indx]

kcorrect,galaxy_maggies,galaxy_invvar,galaxy_z,galaxy_reconstruct_maggies, $
  coeff=coeff,version=version, vpath=vpath, /maggies, /invvar, $
  addgrgap=addgrgap, vconstraint=vconstraint, sdssfix=sdssfix

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
axis_char_scale= 2.0
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

!p.multi=[nk,1,nk]
bands=['u','g','r','i','z']
siglim=djsig(galaxy_reconstruct_maggies[2,*]/galaxy_maggies[2,*],sigrej=3)
for k=0l, nk-1 do begin
    !X.CHARSIZE = tiny
    !Y.CHARSIZE = 1.2*axis_char_scale
    !X.TITLE = ''
    !Y.TITLE = 'f!d'+bands[k]+',r!n/f!d'+bands[k]+',o!n'
    if (k eq nk-1) then !X.CHARSIZE = 1.2*axis_char_scale
    if (k eq nk-1) then !X.TITLE = 'Redshift z'
    !X.RANGE=[0.,0.5]
    !Y.RANGE=1.+nsig*[-siglim,siglim]
    sig=djsig(galaxy_reconstruct_maggies[k,*]/galaxy_maggies[k,*],sigrej=3)
    djs_plot,galaxy_z,galaxy_reconstruct_maggies[k,*]/galaxy_maggies[k,*], $
      psym=3,xst=1,yst=1
    xyouts,!X.RANGE[1]-0.18*(!X.RANGE[1]-!X.RANGE[0]), $
      !Y.RANGE[0]+0.08*(!Y.RANGE[1]-!Y.RANGE[0]), $
      '!4r!3='+strtrim(string(sig,format='(f8.2)'),2), $
      charsize=1.*axis_char_scale,charthick=5
endfor

if(keyword_set(psfile)) then begin
    device, /close
    set_plot,'x'
endif
    
end
;------------------------------------------------------------------------------
