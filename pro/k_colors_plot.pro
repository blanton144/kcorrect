;+
; NAME:
;   k_colors_plot
;
; PURPOSE:
;   Make plot of colors
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
;   06-Feb-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_colors_plot,savfile,version=version,vpath=vpath,psfile=psfile, $
                  to_z=to_z,omega0=omega0, omegal0=omegal0,lumlim=lumlim, $
                  zrange=zrange,nsig=nsig, zsplit=zsplit, $
                  nprimtargetmask=nprimtargetmask, colorlimits=colorlimits, $
                  subsample=subsample, addgrgap=addgrgap, sdssfix=sdssfix, $
                  vconstraint=vconstraint

if(n_elements(vpath) eq 0) then $
  vpath=getenv('KCORRECT_DIR')+'/data/etemplates'
if(n_elements(version) eq 0) then version='default'
if(n_elements(nsig) eq 0) then nsig=3.d
if(n_elements(omega0) eq 0) then omega0=0.3d
if(n_elements(omegal0) eq 0) then omegal0=0.7d
if(n_elements(lumlim) eq 0) then lumlim=[-21.3,-21.0]
if(n_elements(zrange) eq 0) then zrange=[0.,0.5]
if(n_elements(zsplit) eq 0) then zsplit=[0.05,0.1,0.15]

restore,savfile
if(n_elements(sp) gt 0) then galaxy_z=sp.z
if(n_elements(z) gt 0) then galaxy_z=z
ngalaxy=long(n_elements(galaxy_z))
nk=long(n_elements(galaxy_maggies))/ngalaxy

if(n_elements(nprimtargetmask) gt 0) then begin
    help,sp
    indx=where((sp.primtarget and nprimtargetmask) eq 0,count)
    if(count eq 0) then return
    sp=sp[indx]
    galaxy_z=galaxy_z[indx]
    galaxy_maggies=galaxy_maggies[*,indx]
    coeff=coeff[*,indx]
    help,indx
    ngalaxy=n_elements(indx)
endif

if(n_elements(subsample) eq 0) then subsample=1l
indx=lindgen(ngalaxy/long(subsample))*long(subsample)
galaxy_z=galaxy_z[indx]
galaxy_maggies=galaxy_maggies[*,indx]
galaxy_invvar=galaxy_invvar[*,indx]
coeff=coeff[*,indx]

if(keyword_set(to_z)) then begin
    to_z_str=strtrim(string(to_z,format='(f3.1)'),2)
    to_z=galaxy_z-galaxy_z+to_z
endif else begin
    to_z=galaxy_z
    to_z_str='z'
endelse

kcorrect,galaxy_maggies,galaxy_invvar,galaxy_z,recmaggies,kcorrectz=to_z, $
  version=version,vpath=vpath,/maggies,/invvar,addgrgap=addgrgap, $
  vconstraint=vconstraint, sdssfix=sdssfix
kcorrect,galaxy_maggies,galaxy_invvar,galaxy_z,recmaggies0, $
  version=version,vpath=vpath,/maggies,/invvar,addgrgap=addgrgap, $
  vconstraint=vconstraint, sdssfix=sdssfix
kcorrect_val=recmaggies/recmaggies0
recmaggies=galaxy_maggies*kcorrect_val

dm=2.5*alog10((2.99792e+8*lumdis(galaxy_z,omega0,omegal0))^2)
lum=-2.5*alog10(recmaggies[2,*])-dm
indx=where(lum gt lumlim[0] and lum lt lumlim[1] and $
           galaxy_z lt zrange[1] and galaxy_z gt zrange[0] ,count)
help,indx
if(count eq 0) then begin
    klog,'No objects.'
    return
endif

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

!p.multi=[2*(nk-1),2,nk-1]
bands=['u','g','r','i','z']
!Y.MARGIN=0.3*axis_char_scale*[1,5]
for k=0l, nk-2 do begin
    useindx=indx[where(recmaggies[k,indx] gt 0. and $
                       recmaggies[k+1,indx] gt 0.)]
    color=transpose(-2.5*alog10(recmaggies[k,useindx]/recmaggies[k+1,useindx]))
    mn=(djs_avsigclip(color,sigrej=5))[0]
    sig=djsig(color,sigrej=5)
    !X.CHARSIZE = 1.*axis_char_scale
    !Y.CHARSIZE = 1.*axis_char_scale
    !X.TITLE = ''
    !Y.TITLE = '!u'+to_z_str+'!n('+ $
      bands[k]+'-'+bands[k+1]+')'
    if (k eq nk-2) then !X.CHARSIZE = 1.*axis_char_scale
    if (k eq nk-2) then !X.TITLE = 'Redshift z'
    !X.RANGE=zrange
    !Y.RANGE=[mn,mn]+nsig*[-sig,sig]
    if(n_elements(colorlimits) gt 0) then begin
        !Y.RANGE=colorlimits[*,k]
    endif
    djs_plot,galaxy_z[useindx], color, psym=3,xst=1,yst=1
;    xyouts,!X.RANGE[1]-0.18*(!X.RANGE[1]-!X.RANGE[0]), $
;      !Y.RANGE[0]+0.08*(!Y.RANGE[1]-!Y.RANGE[0]), $
;      '!4r!3='+strtrim(string(sig,format='(f8.2)'),2), $
;      charsize=1.*axis_char_scale,charthick=5
    hindx=where(color gt mn-nsig*sig and color lt mn+nsig*sig and $
                galaxy_z[useindx] gt zsplit[0] and $
                galaxy_z[useindx] lt zsplit[1])
    hmax=colorlimits[1,k]
    hmin=colorlimits[0,k]
    hbins=30
    colorvals=hmin+(hmax-hmin)*(dindgen(hbins)+0.5)/double(hbins)
    colorhistlo=histogram(color[hindx],nbins=hbins,max=hmax,min=hmin)
    hindx=where(color gt mn-nsig*sig and color lt mn+nsig*sig and $
                galaxy_z[useindx] gt zsplit[1] and $
                galaxy_z[useindx] lt zsplit[2])
    colorhisthi=histogram(color[hindx],nbins=hbins,max=hmax,min=hmin)
    colorhistloerr=sqrt(colorhistlo)
    zindx=where(colorhistloerr eq 0,count)
    if(count gt 0) then colorhistloerr[zindx]=1.d
    colorhisthierr=sqrt(colorhisthi)
    zindx=where(colorhisthierr eq 0,count)
    if(count gt 0) then colorhisthierr[zindx]=1.d
    colorhistallerr=sqrt(colorhistloerr^2+colorhisthierr^2)
    numer=total(colorhisthi*colorhistlo/(colorhistallerr^2),/double)
    denom=total(colorhisthi*colorhisthi/(colorhistallerr^2),/double)
    scale=numer/denom
    colorhistlo=colorhistlo/scale
    maxval=max([colorhisthi,colorhistlo])
    colorhistlo=colorhistlo/maxval
    colorhisthi=colorhisthi/maxval
    !X.RANGE=!Y.RANGE
    !Y.RANGE=[0.,1.05*max([colorhistlo,colorhisthi])]
    !Y.CHARSIZE=tiny
    !X.CHARSIZE=1.*axis_char_scale
    if (k eq nk-2) then !X.TITLE = 'Color'
    djs_plot,colorvals,colorhistlo,psym=10,xst=1,yst=1
    djs_oplot,colorvals,colorhisthi,psym=10,linestyle=1
    axis,!X.RANGE[1],!Y.RANGE[0],yaxis=1,ycharsize=axis_char_scale, $
      ytitle='N!dgal!n'
    if(k eq 1) then begin
        xyouts,!X.RANGE[0]+(!X.RANGE[1]-!X.RANGE[0])*0.03, $
          !Y.RANGE[1]-(!Y.RANGE[1]-!Y.RANGE[0])*0.15, $
          'Solid:  '+strtrim(string(zsplit[0],format='(d5.2)'),2) + $
          ' z < '+strtrim(string(zsplit[1],format='(d5.2)'),2)
        xyouts,!X.RANGE[0]+(!X.RANGE[1]-!X.RANGE[0])*0.03, $
          !Y.RANGE[1]-(!Y.RANGE[1]-!Y.RANGE[0])*0.30, $
          'Dotted: '+strtrim(string(zsplit[1],format='(d5.2)'),2) + $
          ' z < '+strtrim(string(zsplit[2],format='(d5.2)'),2)
    endif
endfor

if(keyword_set(psfile)) then begin
    device, /close
    set_plot,'x'
endif
    
end
;------------------------------------------------------------------------------
