pro test_galex

; get photometry
galex=mrdfits('/global/data/galex/ero/raw/galex-ero-photom-resolve.fits.gz',1)
abmaggies=galex.flux*10.^(0.4*galex.extinction)
abmaggies_ivar=galex.flux_ivar*10.^(-0.8*galex.extinction)

; ab correct
abmaggies[2:6,*]=sdssflux2ab(abmaggies[2:6,*])
abmaggies_ivar[2:6,*]=sdssflux2ab(abmaggies_ivar[2:6,*],/ivar)

; add errors to prevent overfitting 
abmaggies_ivar[0,*]=1./(1./abmaggies_ivar[0,*]+(0.1*abmaggies[0,*])^2)
abmaggies_ivar[1,*]=1./(1./abmaggies_ivar[1,*]+(0.1*abmaggies[1,*])^2)
abmaggies_ivar[2,*]=1./(1./abmaggies_ivar[2,*]+(0.05*abmaggies[2,*])^2)
abmaggies_ivar[3:6,*]=1./(1./abmaggies_ivar[3:6,*]+(0.02*abmaggies[3:6,*])^2)

kcorrect, abmaggies, abmaggies_ivar, galex.z, kcorrect, coeffs=coeffs, $
  filterlist=['galex_FUV.par','galex_NUV.par','sdss_u0.par', $
              'sdss_g0.par','sdss_r0.par','sdss_i0.par','sdss_z0.par'], $
  rmatrix=rmatrix, zvals=zvals, chi2=chi2
k_reconstruct_maggies, coeffs, galex.z, rmaggies, rmatrix=rmatrix, $
  zvals=zvals

band=['F','N','u','g','r','i','z']
set_print,filename='~/galex-default.ps'
for i=0, 6 do $
  djs_plot, galex.z, abmaggies[i,*]/rmaggies[i,*], psym=4, yra=[-0.5,4.5], $
  xtitle='z', ytitle='f_'+band[i]+'/f_{m,'+band[i]+'}'
end_print

undefine,rmatrix
undefine,zvals

k_load_vmatrix,vmatrix,lambda, $
  vpath='/global/data/sdss/kcorrect/data/v3_3', vfile='vmatrix.test9.dat', $
  lfile='lambda.test9.dat'
nv=n_elements(vmatrix)/(n_elements(lambda)-1L)
vmatrix=vmatrix/1.e+38
iuse=lindgen(n_elements(galex))
coeffs=k_fit_nonneg(abmaggies[*,iuse], abmaggies_ivar[*,iuse],vmatrix, $
                    lambda,redshift=galex[iuse].z, $
                    filterlist=['galex_FUV.par','galex_NUV.par', $
                                'sdss_u0.par','sdss_g0.par','sdss_r0.par', $
                                'sdss_i0.par','sdss_z0.par'], $
                    chi2=chi2,rmatrix=rmatrix,zvals=zvals,maxiter=10000, $
                    /verbose)
k_reconstruct_maggies, coeffs, galex[iuse].z, rmaggies, rmatrix=rmatrix, $
  zvals=zvals

band=['F','N','u','g','r','i','z']
set_print,filename='~/galex-full9.ps'
for i=0, 6 do $
  djs_plot, galex[iuse].z, abmaggies[i,iuse]/rmaggies[i,iuse], psym=4, $
  yra=[-0.5,4.5], xtitle='z', ytitle='f_'+band[i]+'/f_{m,'+band[i]+'}'
end_print

galfit1={abmaggies:fltarr(7), $
         abmaggies_ivar:fltarr(7), $
         z:0., $
         coeffs:fltarr(n_elements(coeffs)/n_elements(galex)), $
         rmaggies:fltarr(7) }
galfit=replicate(galfit1,n_elements(iuse))
galfit.abmaggies=abmaggies[*,iuse]
galfit.abmaggies_ivar=abmaggies_ivar[*,iuse]
galfit.z=galex[iuse].z
galfit.coeffs=coeffs
galfit.rmaggies=rmaggies
mwrfits,galfit,'~/galex_fit.fits',/create
stop

end
