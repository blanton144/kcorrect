
restore,'cnoc2_sdssmatch.sav'

a_indx=where(a.z lt 0.8 and $
             a.UBVRI[0] lt 29. and $
             a.UBVRI[1] lt 29. and $
             a.UBVRI[2] lt 29. and $
             a.UBVRI[3] lt 29. and $
             a.UBVRI[4] lt 29. and $
             obj.matchdist*3600. lt 3.)
a=a[a_indx]
childobj=childobj[a_indx]
shiftband=10.^(-0.4*([-0.042,0.036,0.015,0.013,-0.002]))
for i=0,4 do $
  childobj.modelflux[i]=childobj.modelflux[i]*shiftband[i]
for i=0,4 do $
  childobj.modelflux_ivar[i]=childobj.modelflux_ivar[i]/shiftband[i]^2

bessell_filters=['bessell_U.par', $
                 'bessell_B.par', $
                 'bessell_V.par', $
                 'bessell_R.par', $
                 'bessell_I.par']
vega2ab_UBVRI=k_vega2ab(filterlist=bessell_filters,/kurucz)
for i=0, 4 do $
  a.UBVRI[i]=a.UBVRI[i]+vega2ab_UBVRI[i]

set_print,filename='check_cnoc2.ps'
!X.MARGIN=5.
!Y.MARGIN=3.

kcorrect,childobj.modelflux,childobj.modelflux_ivar,a.z,kcorrect, $
  coeffs=coeffs, vmatrix=vmatrix, lambda=lambda, filterlist=sdss_filters
k_reconstruct_maggies,coeffs,a.z,ugriz_maggies,vmatrix=vmatrix,lambda=lambda, $
  filterlist=sdss_filters
k_reconstruct_maggies,coeffs,a.z,UBVRI_maggies,vmatrix=vmatrix,lambda=lambda, $
  filterlist=bessell_filters
ugriz_mags=22.5-2.5*alog10(ugriz_maggies)
!P.MULTI=[0,2,3]
for i=0, 4 do begin & $
  diff=ugriz_mags[i,*]-(22.5-2.5*alog10(childobj.modelflux[i])) & $
  plot,22.5-2.5*alog10(childobj.modelflux[i]), diff, $
  xra=[16.,25.],yra=[-0.7,0.7], psym=3, $
  xtitle='SDSS '+sdss_filters[i]+' (AB)', $
  ytitle='SDSS(model)-SDSS' & $
  oplot,[0,50],replicate(djs_avsigclip(diff[where(diff eq diff and diff lt 1.e+30)]),2), $
  color=djs_icolor('red') & $
  endfor
UBVRI_mags=22.5-2.5*alog10(UBVRI_maggies)
!P.MULTI=[0,2,3]
for i=0, 4 do begin & $
  diff=UBVRI_mags[i,*]-a.UBVRI[i] & $
  plot,a.UBVRI[i],diff,xra=[16.,25.], $
  yra=[-0.7,0.7],psym=3, $
  xtitle='CNOC2 '+bessell_filters[i]+' (AB)', $
  ytitle='SDSS-CNOC2' & $
  oplot,[0,50],replicate(djs_avsigclip(diff[where(diff eq diff and diff lt 1.e+30)]),2), $
  color=djs_icolor('red') & $
  endfor
!P.MULTI=[0,2,2]
for i=0, 3 do begin & $
  diff=(UBVRI_mags[i,*]-UBVRI_mags[i+1,*])-(a.UBVRI[i]-a.UBVRI[i+1]) & $
  plot,a.UBVRI[i],diff,xra=[16.,25.], $
  yra=[-0.7,0.7],psym=3, $
  xtitle='CNOC2 '+bessell_filters[i]+' (AB)', $
  ytitle='SDSS-CNOC2' & $
  oplot,[0,50],replicate(djs_avsigclip(diff[where(diff eq diff and diff lt 1.e+30)]),2), $
  color=djs_icolor('red') & $
  endfor
print,djs_avsigclip(diff[where(diff eq diff and diff lt 1.e+30)])

kcorrect,a.UBVRI,a.UBVRI_ivar,a.z,kcorrect, /magnitude, $
  coeffs=coeffs, vmatrix=vmatrix, lambda=lambda, $
  filterlist=bessell_filters
k_reconstruct_maggies,coeffs,a.z,ugriz_maggies,vmatrix=vmatrix,lambda=lambda, $
  filterlist=sdss_filters
k_reconstruct_maggies,coeffs,a.z,UBVRI_maggies,vmatrix=vmatrix,lambda=lambda, $
  filterlist=bessell_filters
UBVRI_mags=-2.5*alog10(UBVRI_maggies)
!P.MULTI=[0,2,3]
for i=0, 4 do begin & $
  diff=UBVRI_mags[i,*]-a.UBVRI[i] & $
  plot,a.UBVRI[i], diff, xra=[16.,25.], $
  yra=[-0.7,0.7],psym=3, $
  xtitle='CNOC2 '+bessell_filters[i]+' (AB)', $
  ytitle='CNOC2(model)-CNOC2' & $
  oplot,[0,50],replicate(djs_avsigclip(diff[where(diff eq diff and diff lt 1.e+30)]),2), $
  color=djs_icolor('red') & $
  endfor
ugriz_mags=-2.5*alog10(ugriz_maggies)
!P.MULTI=[0,2,3]
for i=0, 4 do begin & $
  diff=ugriz_mags[i,*]-(22.5-2.5*alog10(childobj.modelflux[i])) & $
  plot,22.5-2.5*alog10(childobj.modelflux[i]), diff, $
  xra=[16.,25.],yra=[-0.7,0.7], psym=3, $
  xtitle='SDSS '+sdss_filters[i]+' (AB)', $
  ytitle='CNOC2-SDSS' & $
  oplot,[0,50],replicate(djs_avsigclip(diff[where(diff eq diff and diff lt 1.e+30)]),2), $
  color=djs_icolor('red') & $
  endfor
!P.MULTI=[0,2,2]
for i=0, 3 do begin & $
  diff=(ugriz_mags[i,*]-ugriz_mags[i+1,*])- $
  ((22.5-2.5*alog10(childobj.modelflux[i]))- $
   (22.5-2.5*alog10(childobj.modelflux[i+1])) ) & $
  plot,22.5-2.5*alog10(childobj.modelflux[i]), diff, $
  xra=[16.,25.],yra=[-0.7,0.7], psym=3, $
  xtitle='SDSS '+sdss_filters[i]+' (AB)', $
  ytitle='CNOC2-SDSS' & $
  oplot,[0,50],replicate(djs_avsigclip(diff[where(diff eq diff and diff lt 1.e+30)]),2), $
  color=djs_icolor('red') & $
  endfor

end_print
