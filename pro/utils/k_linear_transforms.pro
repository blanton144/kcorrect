;+
; NAME:
;   k_linear_transforms
; PURPOSE:
;   calculate transformations between filter systems
; CALLING SEQUENCE:
;   k_linear_tranforms
; COMMENTS:
;   Creates plots in k_linear_transforms.ps
;   Prints lineareqs.tex format linear equations to screen for input into
; REVISION HISTORY:
;   10-Feb-2004  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro klt_plot, color, bdiff, colorstr, bdiffstr, slope, offset, slopestr, $
              offsetstr, toband=toband, fromband=fromband, unit=unit, $
              fromcolor=fromcolor

meancolor=mean(color)
sigcolor=djsig(color,sigrej=6)

iuse=where(abs(color-meancolor) lt 5.*sigcolor)
color=color[iuse]
bdiff=bdiff[iuse]

coeffs=linfit(color-meancolor,bdiff)
slope=coeffs[1]
offset=coeffs[0]

xline=minmax(color)
yline=offset+slope*(xline-meancolor)

djs_plot, color, bdiff, xtitle=colorstr, ytitle=bdiffstr, psym=3
djs_oplot, xline, yline, color='red'
xst=!X.CRANGE[0]+0.05*(!X.CRANGE[1]-!X.CRANGE[0])
yst=!Y.CRANGE[0]+0.05*(!Y.CRANGE[1]-!Y.CRANGE[0])
offsetstr=strtrim(string(offset,format='(f40.4)'),2)
slopestr=strtrim(string(slope,format='(f40.4)'),2)
meancolorstr=strtrim(string(meancolor,format='(f40.4)'),2)
djs_xyouts, xst, yst, 'y = '+offsetstr+' + ('+slopestr+') (x - ('+ $
  meancolorstr+')'

slopestr=strtrim(string(abs(slope),format='(f40.4)'),2)
offsetstr=strtrim(string(abs(offset),format='(f40.4)'),2)
meancolorstr=strtrim(string(abs(meancolor),format='(f40.4)'),2)
sigcolorstr=strtrim(string((sigcolor),format='(f40.2)'),2)
if(slope lt 0.) then $
  slopestr='- '+slopestr $
else $
  slopestr='+ '+slopestr 
if(offset lt 0.) then $
  offsetstr='- '+offsetstr $
else $
  offsetstr='+ '+offsetstr 
if(meancolor gt 0.) then $
  meancolorstr='- '+meancolorstr $
else $
  meancolorstr='+ '+meancolorstr 
printf, unit, '$'+toband+' = '+fromband+' '+offsetstr+' '+slopestr+ $
  ' \left[ ('+fromcolor+') '+meancolorstr+' \right] '+ $
  '$ & $\sigma\left['+fromcolor+'\right] = '+sigcolorstr+'$ \cr' 

end
;
pro k_linear_transforms

sample='dr4'

post=mrdfits(vagc_name('post_catalog', sample=sample, letter='bsafe', $
                       post='1'),1)
ii=where(post.absm[0]-post.absm[1] lt 1.7)
post=post[ii]
post=post[shuffle_indx(n_elements(post), num_sub=10000)]
cat=mrdfits(vagc_name('object_sdss_imaging'),1,row=post.object_position)
kc01=sdss_kcorrect(post.z, calibobj=cat, band_shift=0.1, absmag=absmag01, $
                   mtol=mtol, mass=mass, intsfh=intsfh, $
                   b1000=b1000,b300=b300, coeffs=coeffs01)
kc=sdss_kcorrect(post.z, calibobj=cat, band_shift=0., absmag=absmag, $
                 mtol=mtol, mass=mass, intsfh=intsfh, $
                 b1000=b1000,b300=b300, coeffs=coeffs)
kcb=sdss2bessell(post.z, calibobj=cat, band_shift=0., absmag=absmagb, $
                 mtol=mtol, mass=mass, intsfh=intsfh, $
                 b1000=b1000,b300=b300)


bands=['u', $
       'g', $
       'r', $
       'i', $
       'z', $
       '\band{0.1}{u}', $
       '\band{0.1}{g}', $
       '\band{0.1}{r}', $
       '\band{0.1}{i}', $
       '\band{0.1}{z}', $
       'U', $
       'B', $
       'V', $
       'R', $
       'I']

filterlist=['sdss_u0.par', $
            'sdss_g0.par', $
            'sdss_r0.par', $
            'sdss_i0.par', $
            'sdss_z0.par', $
            'bessell_U.par', $
            'bessell_B.par', $
            'bessell_V.par', $
            'bessell_R.par', $
            'bessell_I.par']
leff00=k_lambda_eff(filterlist=filterlist)
leff01=k_lambda_eff(filterlist=filterlist,band_shift=0.1)
leff=[leff00[0:4], leff01[0:4], leff00[5:9]]

mags=fltarr(15, n_elements(post))
mags[0:4,*]=absmag
mags[5:9,*]=absmag01
mags[10:14,*]=absmagb

links=[[0, 5,6], $
       [1, 6,7], $
       [2, 7,8], $
       [3, 8,9], $
       [4, 8,9], $
       [0, 10,11], $
       [1, 11,12], $
       [2, 12,13], $
       [3, 13,14], $
       [4, 13,14], $
       [5, 10,11], $
       [6, 10,11], $
       [6, 11,12], $
       [7, 11,12], $
       [7, 12,13], $
       [8, 13,14], $
       [9, 13,14], $
       [5, 0,1], $
       [6, 0,1], $
       [6, 1,2], $
       [7, 1,2], $
       [8, 2,3], $
       [9, 3,4], $
       [10, 0,1], $
       [11, 0,1], $
       [11, 1,2], $
       [12, 1,2], $
       [13, 2,3], $
       [14, 3,4], $
       [10, 5,6], $
       [11, 6,7], $
       [12, 6,7], $
       [13, 7,8], $
       [14, 8,9]]

neq=(size(links,/dim))[1]

k_print, filename=getenv('KCORRECT_DIR')+ $
  '/docs/paper/k_linear_transform.ps'
openw, unit, getenv('KCORRECT_DIR')+'/docs/paper/lineareqs.tex', /get_lun
for i=0L, neq-1L do begin
    iorig=links[0,i]
    l0=links[1,i]
    l1=links[2,i]
    mag0=mags[l0,*]
    mag1=mags[l1,*]
    mag=mags[iorig,*]
    colorstr=bands[l0]+'-'+bands[l1]
    bdiffstr=bands[iorig]+'-'+bands[l0]
    klt_plot, mag0-mag1, mag-mag0, colorstr, bdiffstr, unit=unit, $
      fromband=bands[l0], toband=bands[iorig], $
      fromcolor=bands[l0]+'-'+bands[l1]
endfor
free_lun,unit
k_end_print


end
;------------------------------------------------------------------------------

