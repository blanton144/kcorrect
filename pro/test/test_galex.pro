pro test_galex

; get photometry
galex=mrdfits('/global/data/galex/ero/raw/galex-ero-photom-resolve.fits.gz',1)
abmaggies=galex.flux
abmaggies_ivar=galex.flux_ivar

; drop FUV
abmaggies_ivar[0,*]=0.

; ab correct
abmaggies[2:6,*]=sdssflux2ab(abmaggies[2:6,*])
abmaggies_ivar[2:6,*]=sdssflux2ab(abmaggies_ivar[2:6,*],/ivar)

; extinction correct
splog,'WARNING: not extinction correcting because i am too stupid'

; add errors to prevent overfitting sdss
abmaggies_ivar[2,*]=1./(1./abmaggies_ivar[2,*]+(0.05*abmaggies[2,*])^2)
abmaggies_ivar[2,*]=1./(1./abmaggies_ivar[2,*]+(0.05*abmaggies[2,*])^2)
abmaggies_ivar[3:6,*]=1./(1./abmaggies_ivar[3:6,*]+(0.02*abmaggies[3:6,*])^2)

kcorrect, abmaggies, abmaggies_ivar, photo[*,0].z, kcorrect, coeffs=coeffs, $
  filterlist=['galex_FUV.par','galex_NUV.par','sdss_u0.par', $
              'sdss_g0.par','sdss_r0.par','sdss_i0.par','sdss_z0.par'], $
  rmatrix=rmatrix, zvals=zvals
k_reconstruct_maggies, coeffs, photo[*,0].z, rmaggies, rmatrix=rmatrix, $
  zvals=zvals

band=['F','N','u','g','r','i','z']
set_print,filename='~/galex-default.ps'
for i=0, 6 do $
  djs_plot, photo[*,0].z, abmaggies[i,*]/rmaggies[i,*], psym=4, yra=[-0.5,4.5], $
  xtitle='z', ytitle='f_'+band[i]+'/f_{m,'+band[i]+'}'
end_print

undefine,rmatrix
undefine,zvals

k_load_vmatrix,vmatrix,lambda, $
  vpath='/global/data/sdss/kcorrect/data/v3_2_0', vfile='vmatrix.v3_2.dat', $
  lfile='lambda.v3_2.dat'
nv=n_elements(vmatrix)/(n_elements(lambda)-1L)
vmatrix=vmatrix/1.e+38
iuse=lindgen(n_elements(photo[*,0]))
coeffs=k_fit_nonneg(abmaggies[*,iuse], abmaggies_ivar[*,iuse],vmatrix, $
                    lambda,redshift=photo[iuse,0].z, $
                    filterlist=['galex_FUV.par','galex_NUV.par', $
                                'sdss_u0.par','sdss_g0.par','sdss_r0.par', $
                                'sdss_i0.par','sdss_z0.par'], $
                    chi2=chi2,rmatrix=rmatrix,zvals=zvals,maxiter=10000, $
                    /verbose)
k_reconstruct_maggies, coeffs, photo[iuse,0].z, rmaggies, rmatrix=rmatrix, $
  zvals=zvals

set_print,filename='~/galex-full.ps'
for i=0, 6 do $
  djs_plot, photo[*,0].z, abmaggies[i,*]/rmaggies[i,*], psym=4, yra=[-0.5,4.5], $
  xtitle='z', ytitle='f_'+band[i]+'/f_{m,'+band[i]+'}'
end_print

galfit1={abmaggies:fltarr(7), $
         abmaggies_ivar:fltarr(7), $
         z:0., $
         coeffs:fltarr(n_elements(coeffs)/n_elements(photo[*,0])), $
         rmaggies:fltarr(7) }
galfit=replicate(galfit1,n_elements(photo[*,0]))
galfit.abmaggies=abmaggies
galfit.abmaggies_ivar=abmaggies_ivar
galfit.z=photo[*,0].z
galfit.coeffs=coeffs
galfit.rmaggies=rmaggies
mwrfits,galfit,'~/galex_fit.fits',/create
stop

end
