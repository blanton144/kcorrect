;+
; NAME:
;   lf_zdist_plot
; PURPOSE:
;   plot distribution of distance in bins of absmag, along with model
; CALLING_SEQUENCe:
;   lf_zdist_plot,zz_data,am_data,zz_model,am_model, filename= [ ,fadjust=, $
;                 amlimits=, sample_zmin=, sample_zmax=, nbins=, $
;                 band=, /ylog ]
; INPUTS:
;   zz_data - redshifts of data
;   am_data - absmags  of data
;   zz_model - redshifts of model
;   am_model - absmags  of model
;   filename - output postscript file name
; OPTIONAL INPUTS:
;   fadjust - adjust for sampling fraction (default 1.)
;   amlimits - [2] limiting abs, mag. limits of plot
;   sample_zmin, sample_zmax - z limits to use
;   nbins - number of redshift bins 
;   band - 0-4 corresponds to ugriz
; KEYWORDS:
;   /ylog - plot on log scale
; REVISION HISTORY:
;   2002-5-22  written - Blanton
;-
pro lf_zdist_plot,zz_data,am_data,zz_model,am_model,fadjust=fadjust, $
                  amlimits=amlimits,sample_absmmin=sample_absmmin, $
                  sample_absmmax=sample_absmmax,nbins=nbins, $
                  filename=filename,sample_zmin=sample_zmin, $
                  sample_zmax=sample_zmax,band=band,ylog=ylog

; defaults
if(NOT keyword_set(fadjust)) then fadjust=1.
if(NOT keyword_set(amlimits)) then amlimits=[min(am_model),max(am_model)]
if(NOT keyword_set(nbins)) then nbins=50
if(NOT keyword_set(filename)) then filename='lf_zdist.ps'
if(n_elements(sample_absmmin) eq 0) then sample_absmmin=min(am_model)
if(n_elements(sample_absmmax) eq 0) then sample_absmmax=max(am_model)
if(n_elements(sample_zmin) eq 0) then sample_zmin=min(zz_model)
if(n_elements(sample_zmax) eq 0) then sample_zmax=max(zz_model)
dh=2.99792e+5/100.
pi=3.14159265358979D+0
bandstr=['u','g','r','i','z']

; setup postscript file
pold=!P
xold=!X
yold=!Y
!P.FONT= -1
set_plot, "PS"
!P.BACKGROUND= djs_icolor('white')
!P.COLOR= djs_icolor('black')
xsize= 5.5 & ysize= 5.5
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
!P.MULTI= [n_elements(amlimits)-1,1,n_elements(amlimits)-1]
xyouts, 0,0,'!6'
colorname= ['red','green','blue','magenta','cyan','dark yellow', $
            'purple','light green','orange','navy','light magenta', $
            'yellow green']
ncolor= n_elements(colorname)
theta= 2.0D *double(!PI)*dindgen(31)/30.0D
x= cos(theta)
y= sin(theta)

for i=0, n_elements(amlimits)-2 do begin
    plot_ammax=amlimits[i+1]
    plot_ammin=amlimits[i]
    indx_data=where(am_data lt plot_ammax and am_data gt plot_ammin)
    indx_model=where(am_model lt plot_ammax and am_model gt plot_ammin)
    
    vals=sample_zmin+(dindgen(nbins)+0.5)*(sample_zmax-sample_zmin)/ $
      double(nbins)
    
    hist_data=double(histogram(zz_data[indx_data],nbins=nbins, $
                               min=sample_zmin,max=sample_zmax))
    indx=where(hist_data eq 0,count)
    if(count gt 0) then hist_data[indx]=0.1

    hist_model=fadjust*double(histogram(zz_model[indx_model],nbins=nbins, $
                                        min=sample_zmin,max=sample_zmax))
    indx=where(hist_model eq 0,count)
    if(count gt 0) then hist_model[indx]=0.1

    maxhist=max([hist_model,hist_data])
    minhist=min([hist_model,hist_data])
    if(NOT keyword_set(ylog)) then $
      !Y.range=[-0.1,(maxhist*1.35)] $
    else $
      !Y.range=[minhist*.95,(maxhist*1.05)]
    !X.charsize=tiny
    if(i eq n_elements(amlimits)-2) then !X.charsize=axis_char_scale
    plot,vals,(hist_model),psym=10,thick=2, $
      ytitle='!6!8N!6', xtitle='!6Redshift !8z!6', ylog=ylog
    oplot,vals,(hist_data),psym=10,thick=6
    
    !x.title='!6!8M!i!60.1!d!8'+bandstr[band]+'!n!6 - 5 log!d10!n!8h!6'
    outstr=strtrim(string(plot_ammin,format='(f12.2)'),2)+ $
      ' !8< M!i!60.1!d!8'+bandstr[band]+'!n <!6 '+ $
      strtrim(string(plot_ammax,format='(f12.2)'),2)
    xyouts,sample_zmin+(sample_zmax-sample_zmin)*0.05, $
      (0.1+(maxhist)*1.25)*0.87,outstr
endfor

!P.CHARSIZE= 0.92
outstr='!6Thick Histogram: Data ('+strtrim(string(n_elements(zz_data)),2)+ $
  ' galaxies)'
xyouts,sample_zmin+(sample_zmax-sample_zmin)*0.35, $
  (0.1+(maxhist))*0.82,outstr
outstr='!6Thin Histogram: Model ('+ $
  strtrim(string(long(n_elements(zz_model)*fadjust)),2)+' galaxies)'
xyouts,sample_zmin+(sample_zmax-sample_zmin)*0.35, $
  (0.1+(maxhist))*0.58,outstr
device,/close

!P=pold
!X=xold
!Y=yold


end
