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
pro klt_plot, color, bdiff, colorstr, bdiffstr, slope, offset


djs_plot, sdss01[0,*]-sdss01[1,*],  sdss00[0,*]-sdss01[0,*], $
  xtitle=umg01, ytitle=u00+'-'+u01, psym=3
end
;
pro k_linear_transforms

; get galaxy sample
kcorrect=hogg_mrdfits(getenv('VAGC_REDUX')+ $
                      '/kcorrect/kcorrect.none.petro.z0.10.fits',1, $
                      nrowchunk=100000)
ii=where(kcorrect.z gt 0.01 and kcorrect.z lt 0.3 and $
         kcorrect.abmaggies_ivar[7] gt 0., ngals)
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
jab2vega=-k_vega2ab(filterlist=jfilters,/hayes)

johnson00=22.5-2.5*alog10(jmaggies00)
sdss00=22.5-2.5*alog10(smaggies00)
johnson01=22.5-2.5*alog10(jmaggies01)
sdss01=22.5-2.5*alog10(smaggies01)

for i=0, 4 do $
  johnson00[i,*]=johnson00[i,*]+jab2vega[i]
for i=0, 4 do $
  johnson01[i,*]=johnson01[i,*]+jab2vega[i]

u01='!8!u0.1!nu!6'
g01='!8!u0.1!ng!6'
r01='!8!u0.1!nr!6'
i01='!8!u0.1!ni!6'
z01='!8!u0.1!nz!6'
j01='!8!u0.1!nj!6'
h01='!8!u0.1!nh!6'
k01='!8!u0.1!nK_s!6'

u00='!8!u0.0!nu!6'
g00='!8!u0.0!ng!6'
r00='!8!u0.0!nr!6'
i00='!8!u0.0!ni!6'
z00='!8!u0.0!nz!6'
j00='!8!u0.0!nH!6'
h00='!8!u0.0!nJ!6'
k00='!8!u0.0!nK_s!6'

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
djs_plot, sdss01[0,*]-sdss01[1,*],  sdss00[0,*]-sdss01[0,*], $
  xtitle=umg01, ytitle=u00+'-'+u01, psym=3
djs_plot, sdss01[1,*]-sdss01[2,*],  sdss00[1,*]-sdss01[1,*], $
  xtitle=gmr01, ytitle=g00+'-'+g01, psym=3
djs_plot, sdss01[2,*]-sdss01[3,*],  sdss00[2,*]-sdss01[2,*], $
  xtitle=rmi01, ytitle=r00+'-'+r01, psym=3
djs_plot, sdss01[3,*]-sdss01[4,*],  sdss00[3,*]-sdss01[3,*], $
  xtitle=imz01, ytitle=i00+'-'+i01, psym=3
djs_plot, sdss01[3,*]-sdss01[4,*],  sdss00[4,*]-sdss01[4,*], $
  xtitle=imz01, ytitle=z00+'-'+z01, psym=3
djs_plot, sdss01[4,*]-sdss01[5,*],  sdss00[4,*]-sdss01[4,*], $
  xtitle=zmj01, ytitle=z00+'-'+z01, psym=3
djs_plot, sdss01[5,*]-sdss01[6,*],  sdss00[5,*]-sdss01[5,*], $
  xtitle=jmh01, ytitle=j00+'-'+j01, psym=3
djs_plot, sdss01[6,*]-sdss01[7,*],  sdss00[6,*]-sdss01[6,*], $
  xtitle=hmk01, ytitle=h00+'-'+h01, psym=3
djs_plot, sdss01[6,*]-sdss01[7,*],  sdss00[7,*]-sdss01[7,*], $
  xtitle=hmk01, ytitle=k00+'-'+k01, psym=3

; transformations from SDSS/TWOMASS at 0.0 to SDSS/TWOMASS at 0.1

; transformations from SDSS at 0.1 to Johnson at 0.

; transformations from SDSS at 0. to Johnson at 0.

; transformations from JOHNSON at 0.1 to Johnson at 0.

end_print
stop


end
;------------------------------------------------------------------------------

