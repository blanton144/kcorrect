;+
; NAME:
;   k_linear_transforms
; PURPOSE:
;   calculate transformations between filter systems
; CALLING SEQUENCE:
;   k_linear_tranforms
; REVISION HISTORY:
;   10-Feb-2004  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro klt_plot, color, bdiff, colorstr, bdiffstr, slope, offset, slopestr, $
              offsetstr, toband=toband, fromband=fromband, fromcolor=fromcolor

chist=histogram(color,binsize=0.03,omin=cmin,omax=cmax)
cvals=cmin+(cmax-cmin)*(findgen(n_elements(chist))+0.5)/ $
  float(n_elements(chist))
weight=interpolate(1./chist, (color-cmin)/(cmax-cmin)*float(n_elements(chist)))
weight=weight/max(weight)
iuse=where(randomu(seed, n_elements(weight)) lt weight)

usecolor=reform(color[iuse]+0.01*randomn(seed,n_elements(iuse)), $
                n_elements(iuse))
usebdiff=reform(bdiff[iuse]+0.01*randomn(seed,n_elements(iuse)), $
                n_elements(iuse))
coeffs=linfit(usecolor,usebdiff)
slope=coeffs[1]
offset=coeffs[0]

xline=minmax(usecolor)
yline=offset+slope*xline

djs_plot, usecolor, usebdiff, xtitle=colorstr, ytitle=bdiffstr, psym=3
djs_oplot, xline, yline, color='red'
xst=!X.CRANGE[0]+0.05*(!X.CRANGE[1]-!X.CRANGE[0])
yst=!Y.CRANGE[0]+0.05*(!Y.CRANGE[1]-!Y.CRANGE[0])
offsetstr=strtrim(string(offset,format='(f40.4)'),2)
slopestr=strtrim(string(slope,format='(f40.4)'),2)
djs_xyouts, xst, yst, 'y = '+offsetstr+' + ('+slopestr+') x'

slopestr=strtrim(string(abs(slope),format='(f40.4)'),2)
offsetstr=strtrim(string(abs(offset),format='(f40.4)'),2)
if(slope lt 0.) then $
  slopestr='- '+slopestr $
else $
  slopestr='+ '+slopestr 
if(offset lt 0.) then $
  offsetstr='- '+offsetstr $
else $
  offsetstr='+ '+offsetstr 
print, toband+' &=& '+fromband+' '+offsetstr+' '+slopestr+ $
  ' \left[ '+fromcolor+' \right] \cr' 

end
;
pro k_linear_transforms

; get galaxy sample
kcorrect=hogg_mrdfits(getenv('VAGC_REDUX')+ $
                      '/kcorrect/kcorrect.none.petro.z0.10.fits',1, $
                      nrowchunk=100000)
ii=where(kcorrect.z gt 0.01 and kcorrect.z lt 0.3 and $
         kcorrect.coeffs[0] gt 0. and $
         kcorrect.coeffs[1] gt 0. and $
         kcorrect.coeffs[2] gt 0. and $
         kcorrect.abmaggies_ivar[7] gt 0., $
         ngals)
kcorrect=kcorrect[ii]

k_load_vmatrix, vmatrix, lambda, vfile=vfile, lfile=lfile,vpath=vpath

sfilters=['sdss_u0.par', 'sdss_g0.par', 'sdss_r0.par', 'sdss_i0.par', $
          'sdss_z0.par', 'twomass_J.par', 'twomass_H.par', 'twomass_Ks.par']
k_reconstruct_maggies, kcorrect.coeffs, fltarr(ngals), smaggies01, $
  filterlist=sfilters, lambda=lambda, vmatrix=vmatrix, band_shift=0.1
k_reconstruct_maggies, kcorrect.coeffs, fltarr(ngals), smaggies00, $
  filterlist=sfilters, lambda=lambda, vmatrix=vmatrix, band_shift=0.0

jfilters=['bessell_U.par', 'bessell_B.par', 'bessell_V.par', 'bessell_R.par', $
          'bessell_I.par']
k_reconstruct_maggies, kcorrect.coeffs, fltarr(ngals), jmaggies01, $
  filterlist=jfilters, lambda=lambda, vmatrix=vmatrix, band_shift=0.1
k_reconstruct_maggies, kcorrect.coeffs, fltarr(ngals), jmaggies00, $
  filterlist=jfilters, lambda=lambda, vmatrix=vmatrix, band_shift=0.0

jab2vega01=-k_vega2ab(filterlist=jfilters,/kurucz, band_shift=0.1)
jab2vega00=-k_vega2ab(filterlist=jfilters,/kurucz, band_shift=0.0)

sab2vega01=-k_vega2ab(filterlist=sfilters,/kurucz, band_shift=0.1)
sab2vega00=-k_vega2ab(filterlist=sfilters,/kurucz, band_shift=0.0)

johnson00=22.5-2.5*alog10(jmaggies00)
sdss00=22.5-2.5*alog10(smaggies00)
johnson01=22.5-2.5*alog10(jmaggies01)
sdss01=22.5-2.5*alog10(smaggies01)

for i=5, 7 do $
  sdss00[i,*]=sdss00[i,*]+sab2vega00[i]
for i=5, 7 do $
  sdss01[i,*]=sdss01[i,*]+sab2vega01[i]
for i=0, 4 do $
  johnson00[i,*]=johnson00[i,*]+jab2vega00[i]
for i=0, 4 do $
  johnson01[i,*]=johnson01[i,*]+jab2vega01[i]

u01='!8!u0.1!nu!6'
g01='!8!u0.1!ng!6'
r01='!8!u0.1!nr!6'
i01='!8!u0.1!ni!6'
z01='!8!u0.1!nz!6'
j01='!8!u0.1!nJ!6'
h01='!8!u0.1!nH!6'
k01='!8!u0.1!nK_s!6'

u00='!8!u0.0!nu!6'
g00='!8!u0.0!ng!6'
r00='!8!u0.0!nr!6'
i00='!8!u0.0!ni!6'
z00='!8!u0.0!nz!6'
j00='!8!u0.0!nJ!6'
h00='!8!u0.0!nH!6'
k00='!8!u0.0!nK_s!6'

bu00='!8!u0.0!nU!6'
bb00='!8!u0.0!nB!6'
bv00='!8!u0.0!nV!6'
br00='!8!u0.0!nR_c!6'
bi00='!8!u0.0!nI_c!6'

umg01='!8!u0.1!n(u-g)!6'
gmr01='!8!u0.1!n(g-r)!6'
rmi01='!8!u0.1!n(r-i)!6'
imz01='!8!u0.1!n(i-z)!6'
zmj01='!8!u0.1!n(z-J)!6'
jmh01='!8!u0.1!n(J-H)!6'
hmk01='!8!u0.1!n(H-K_s)!6'

umg00='!8!u0.0!n(u-g)!6'
gmr00='!8!u0.0!n(g-r)!6'
rmi00='!8!u0.0!n(r-i)!6'
imz00='!8!u0.0!n(i-z)!6'
zmj00='!8!u0.0!n(z-J)!6'
jmh00='!8!u0.0!n(J-H)!6'
hmk00='!8!u0.0!n(H-K_s)!6'

set_print, filename='k_linear_transform.ps'

; transformations from SDSS/TWOMASS at 0.1 to SDSS/TWOMASS at 0.0
klt_plot, sdss01[0,*]-sdss01[1,*],  sdss00[0,*]-sdss01[0,*], $
  umg01, u00+'-'+u01 , slope, offset, slopestr, offsetstr, $
  toband='\band{0.0}{u}', fromband='\band{0.1}{u}', $
  fromcolor='\band{0.1}{(u-g)}'
klt_plot, sdss01[1,*]-sdss01[2,*],  sdss00[1,*]-sdss01[1,*], $
  gmr01, g00+'-'+g01, slope, offset, $
  toband='\band{0.0}{g}', fromband='\band{0.1}{g}', $
  fromcolor='\band{0.1}{(g-r)}'
klt_plot, sdss01[2,*]-sdss01[3,*],  sdss00[2,*]-sdss01[2,*], $
  rmi01, r00+'-'+r01, slope, offset, $
  toband='\band{0.0}{g}', fromband='\band{0.1}{r}', $
  fromcolor='\band{0.1}{(r-i)}'
klt_plot, sdss01[3,*]-sdss01[4,*],  sdss00[3,*]-sdss01[3,*], $
  imz01, i00+'-'+i01, slope, offset, $
  toband='\band{0.0}{i}', fromband='\band{0.1}{i}', $
  fromcolor='\band{0.1}{(i-z)}'
klt_plot, sdss01[3,*]-sdss01[4,*],  sdss00[4,*]-sdss01[4,*], $
  imz01, z00+'-'+z01, slope, offset, $
  toband='\band{0.0}{z}', fromband='\band{0.1}{z}', $
  fromcolor='\band{0.1}{(i-z)}'
klt_plot, sdss01[4,*]-sdss01[5,*],  sdss00[4,*]-sdss01[4,*], $
  zmj01, z00+'-'+z01, slope, offset, $
  toband='\band{0.0}{z}', fromband='\band{0.1}{z}', $
  fromcolor='\band{0.1}{(z-J)}'
klt_plot, sdss01[5,*]-sdss01[6,*],  sdss00[5,*]-sdss01[5,*], $
  jmh01, j00+'-'+j01, slope, offset, $
  toband='\band{0.0}{J}', fromband='\band{0.1}{J}', $
  fromcolor='\band{0.1}{(J-H)}'
klt_plot, sdss01[6,*]-sdss01[7,*],  sdss00[6,*]-sdss01[6,*], $
  hmk01, h00+'-'+h01, slope, offset, $
  toband='\band{0.0}{H}', fromband='\band{0.1}{H}', $
  fromcolor='\band{0.1}{(H-K_s)}'
klt_plot, sdss01[6,*]-sdss01[7,*],  sdss00[7,*]-sdss01[7,*], $
  hmk01, k00+'-'+k01, slope, offset, $
  toband='\band{0.0}{K_s}', fromband='\band{0.1}{K_s}', $
  fromcolor='\band{0.1}{(H-K_s)}'

; transformations from SDSS/TWOMASS at 0.0 to SDSS/TWOMASS at 0.1
klt_plot, sdss00[0,*]-sdss00[1,*],  sdss01[0,*]-sdss00[0,*], $
  umg00, u01+'-'+u00 , slope, offset, $
  toband='\band{0.0}{u}', fromband='\band{0.1}{u}', $
  fromcolor='\band{0.1}{(u-g)}'
klt_plot, sdss00[0,*]-sdss00[1,*],  sdss01[1,*]-sdss00[1,*], $
  umg00, g01+'-'+g00, slope, offset, $
  toband='\band{0.0}{g}', fromband='\band{0.1}{g}', $
  fromcolor='\band{0.1}{(u-g)}'
klt_plot, sdss00[1,*]-sdss00[2,*],  sdss01[1,*]-sdss00[1,*], $
  gmr00, g01+'-'+g00, slope, offset, $
  toband='\band{0.0}{g}', fromband='\band{0.1}{g}', $
  fromcolor='\band{0.1}{(g-r)}'
klt_plot, sdss00[1,*]-sdss00[2,*],  sdss01[2,*]-sdss00[2,*], $
  gmr00, r01+'-'+r00, slope, offset, $
  toband='\band{0.0}{r}', fromband='\band{0.1}{r}', $
  fromcolor='\band{0.1}{(g-r)}'
klt_plot, sdss00[2,*]-sdss00[3,*],  sdss01[3,*]-sdss00[3,*], $
  rmi00, i01+'-'+i00, slope, offset, $
  toband='\band{0.0}{i}', fromband='\band{0.1}{i}', $
  fromcolor='\band{0.1}{(r-i)}'
klt_plot, sdss00[3,*]-sdss00[4,*],  sdss01[4,*]-sdss00[4,*], $
  imz00, z01+'-'+z00, slope, offset, $
  toband='\band{0.0}{z}', fromband='\band{0.1}{z}', $
  fromcolor='\band{0.1}{(i-z)}'
klt_plot, sdss00[4,*]-sdss00[5,*],  sdss01[5,*]-sdss00[5,*], $
  zmj00, j01+'-'+j00, slope, offset, $
  toband='\band{0.0}{J}', fromband='\band{0.1}{J}', $
  fromcolor='\band{0.1}{(z-J)}'
klt_plot, sdss00[5,*]-sdss00[6,*],  sdss01[5,*]-sdss00[5,*], $
  jmh00, j01+'-'+j00, slope, offset, $
  toband='\band{0.0}{J}', fromband='\band{0.1}{J}', $
  fromcolor='\band{0.1}{(J-H)}'
klt_plot, sdss00[5,*]-sdss00[6,*],  sdss01[6,*]-sdss00[6,*], $
  jmh00, h01+'-'+h00, slope, offset, $
  toband='\band{0.0}{H}', fromband='\band{0.1}{H}', $
  fromcolor='\band{0.1}{(J-H)}'
klt_plot, sdss00[6,*]-sdss00[7,*],  sdss01[7,*]-sdss00[7,*], $
  hmk00, k01+'-'+k00, slope, offset, $
  toband='\band{0.0}{K_s}', fromband='\band{0.1}{K_s}', $
  fromcolor='\band{0.1}{(H-K_s)}'

; transformations from SDSS at 0.1 to Johnson at 0.
klt_plot, sdss01[0,*]-sdss01[1,*],  johnson00[0,*]-sdss01[0,*], $
  umg01, bu00+'-'+u01 , slope, offset, $
  toband='\band{0.0}{U}', fromband='\band{0.1}{u}', $
  fromcolor='\band{0.1}{(u-g)}'
klt_plot, sdss01[1,*]-sdss01[2,*],  johnson00[1,*]-sdss01[1,*], $
  gmr01, bb00+'-'+g01 , slope, offset, $
  toband='\band{0.0}{B}', fromband='\band{0.1}{g}', $
  fromcolor='\band{0.1}{(g-r)}'
klt_plot, sdss01[1,*]-sdss01[2,*],  johnson00[2,*]-sdss01[2,*], $
  gmr01, bv00+'-'+r01 , slope, offset, $
  toband='\band{0.0}{V}', fromband='\band{0.1}{r}', $
  fromcolor='\band{0.1}{(g-r)}'
klt_plot, sdss01[2,*]-sdss01[3,*],  johnson00[3,*]-sdss01[3,*], $
  rmi01, br00+'-'+i01 , slope, offset, $
  toband='\band{0.0}{R_c}', fromband='\band{0.1}{i}', $
  fromcolor='\band{0.1}{(r-i)}'
klt_plot, sdss01[3,*]-sdss01[4,*],  johnson00[4,*]-sdss01[4,*], $
  imz01, bi00+'-'+z01 , slope, offset, $
  toband='\band{0.0}{I_c}', fromband='\band{0.1}{z}', $
  fromcolor='\band{0.1}{(i-z)}'

; transformations from SDSS at 0. to Johnson at 0.
klt_plot, sdss00[0,*]-sdss00[1,*],  johnson00[0,*]-sdss00[0,*], $
  umg00, bu00+'-'+u00 , slope, offset, $
  toband='\band{0.0}{U}', fromband='\band{0.0}{u}', $
  fromcolor='\band{0.0}{(u-g)}'
klt_plot, sdss00[1,*]-sdss00[2,*],  johnson00[1,*]-sdss00[1,*], $
  gmr00, bb00+'-'+g00 , slope, offset, $
  toband='\band{0.0}{B}', fromband='\band{0.0}{g}', $
  fromcolor='\band{0.0}{(g-r)}'
klt_plot, sdss00[1,*]-sdss00[2,*],  johnson00[2,*]-sdss00[2,*], $
  gmr00, bv00+'-'+r00 , slope, offset, $
  toband='\band{0.0}{V}', fromband='\band{0.0}{r}', $
  fromcolor='\band{0.0}{(g-r)}'
klt_plot, sdss00[2,*]-sdss00[3,*],  johnson00[3,*]-sdss00[3,*], $
  rmi00, br00+'-'+i00 , slope, offset, $
  toband='\band{0.0}{R_c}', fromband='\band{0.0}{i}', $
  fromcolor='\band{0.0}{(r-i)}'
klt_plot, sdss00[3,*]-sdss00[4,*],  johnson00[4,*]-sdss00[4,*], $
  imz00, bi00+'-'+z00 , slope, offset, $
  toband='\band{0.0}{I_c}', fromband='\band{0.0}{z}', $
  fromcolor='\band{0.0}{(i-z)}'

end_print


end
;------------------------------------------------------------------------------

