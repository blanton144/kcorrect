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
pro k_speck_plot,savfile,version=version,vpath=vpath,psfile=psfile,nsig=nsig, subsample=subsample,to_z=to_z,zrange=zrange,usefiber=usefiber,lumlim=lumlim

if(NOT keyword_set(version)) then version='default'
if(NOT keyword_set(vpath)) then vpath=getenv('KCORRECT_DIR')+'/data/etemplates'
if(NOT keyword_set(nsig)) then nsig=1.7d
if(NOT keyword_set(to_z)) then to_z=0.25
if(NOT keyword_set(zrange)) then zrange=[0.,0.2]

if(keyword_set(usefiber)) then begin
    restore,'main.default.sav'
    spfull=sp
endif

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


if(keyword_set(usefiber)) then begin
    spherematch,sp.ra,sp.dec,spfull.ra,spfull.dec,1.d/3600.d, $
      match1,match2,distance12
    isort=sort(match1)
    tmpmags=spfull[match2[isort[indx]]].fibercounts[*]
    tmpmagserr=spfull[match2[isort[indx]]].fibercountserr[*]
    maglimit=dblarr(5)+22.5
    errlimit=dblarr(5)+0.8
    k_fix_mags,galaxy_z,tmpmags,tmpmagserr,maglimit,errlimit,8l
    galaxy_maggies=dblarr(5,n_elements(indx))
    galaxy_invvar=dblarr(5,n_elements(indx))
    errband=[0.05,0.02,0.02,0.02,0.03]
    for k = 0,4 do begin
        galaxy_maggies[k,*]= $
          10.d^(-0.4d*(tmpmags[k,*] $
                       -spfull[match2[isort[indx]]].reddening[k]))
        galaxy_invvar[k,*]=galaxy_maggies[k,*]*0.4d*alog(10.d)* $
          sqrt(tmpmagserr[k,*]^2+errband[k]^2)
        galaxy_invvar[k,*]=1.d/(galaxy_invvar[k,*]^2)
    endfor
endif else begin
    galaxy_maggies=galaxy_maggies[*,indx]
    galaxy_invvar=galaxy_invvar[*,indx]
    coeff=coeff[*,indx]
endelse

mags=mags[*,indx]
mags0=mags0[*,indx]
help,indx

if(n_elements(lumlim) gt 0) then begin
    k_fit_coeffs,galaxy_maggies,galaxy_invvar,galaxy_z,coeff,version=version, $
      vpath=vpath
    galaxy_z_k=replicate(to_z,n_elements(galaxy_z))
    k_reconstruct_maggies,coeff,galaxy_z_k,recmaggies, $
      version=version, vpath=vpath

    omega0=0.3
    omegal0=0.7
    dm=2.5*alog10((2.99792e+8*lumdis(galaxy_z,omega0,omegal0))^2)
    lum=-2.5*alog10(recmaggies[2,*])-dm
    indx=where(lum gt lumlim[0] and lum lt lumlim[1])
    
    galaxy_maggies=galaxy_maggies[*,indx]
    galaxy_invvar=galaxy_invvar[*,indx]
    coeff=coeff[*,indx]
    galaxy_z=galaxy_z[indx]
    mags=mags[*,indx]
    mags0=mags0[*,indx]
endif

k_fit_coeffs,galaxy_maggies,galaxy_invvar,galaxy_z,coeff,version=version, $
  vpath=vpath
galaxy_z_k=replicate(to_z,n_elements(galaxy_z))
k_reconstruct_maggies,coeff,galaxy_z_k,galaxy_reconstruct_maggies, $
  version=version, vpath=vpath
k_reconstruct_maggies,coeff,galaxy_z,galaxy_reconstruct_maggies0, $
  version=version, vpath=vpath

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

!p.multi=[nk-2,1,nk-2]
bands=['u','g','r','i','z']
k=1
nzindx=where(galaxy_reconstruct_maggies0[k,*] gt 0. and $
             galaxy_reconstruct_maggies[k,*] gt 0. and $
             mags0[k,*] eq mags0[k,*] and abs(mags0[k,*]) lt 1.d+30 and $
             mags[k,*] eq mags[k,*] and abs(mags[k,*]) lt 1.d+30)
kk=2.5*alog10(galaxy_reconstruct_maggies[k,nzindx]/ $
              galaxy_reconstruct_maggies0[k,nzindx])
kkspec=mags[k,nzindx]-mags0[k,nzindx]
kk=kkspec-kk
sig=djsig(kk,sigrej=5)
for k=1l, nk-2 do begin
    !X.CHARSIZE = tiny
    !Y.CHARSIZE = 1.35*axis_char_scale
    !X.TITLE = ''
    !Y.TITLE = 'K!d!s!e'+strtrim(string(to_z,format='(d4.2)'),2)+'!r  '+ $
      bands[k]+'!n(z)'
    if (k eq nk-2) then !X.CHARSIZE = 1.2*axis_char_scale
    if (k eq nk-2) then !X.TITLE = 'Redshift z'
    !X.RANGE=zrange
    nzindx=where(galaxy_reconstruct_maggies0[k,*] gt 0. and $
                 galaxy_reconstruct_maggies[k,*] gt 0. and $
                 mags0[k,*] eq mags0[k,*] and abs(mags0[k,*]) lt 1.d+30 and $
                 mags[k,*] eq mags[k,*] and abs(mags[k,*]) lt 1.d+30)
    kk=2.5*alog10(galaxy_reconstruct_maggies[k,nzindx]/ $
                   galaxy_reconstruct_maggies0[k,nzindx])
    kkspec=mags[k,nzindx]-mags0[k,nzindx]
    kk=kkspec-kk
    !Y.RANGE=[0.,0.]+nsig*sig*[-1.,1.]
    djs_plot,galaxy_z[nzindx],kk, psym=3,xst=1,yst=1
    ;xyouts,!X.RANGE[1]-0.18*(!X.RANGE[1]-!X.RANGE[0]), $
      ;!Y.RANGE[0]+0.08*(!Y.RANGE[1]-!Y.RANGE[0]), $
      ;'!4r!3='+strtrim(string(sig,format='(f8.2)'),2), $
      ;charsize=1.*axis_char_scale,charthick=5
endfor

if(keyword_set(psfile)) then begin
    device, /close
    set_plot,'x'
endif
    
end
;------------------------------------------------------------------------------
