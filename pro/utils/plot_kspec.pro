pro plot_kspec, coeffs, vmatrix=vmatrix, lambda=lambda, vfile=vfile, $
                vpath=vpath, lfile=lfile, _EXTRA=extra_for_splot

if(n_elements(vmatrix) eq 0 or n_elements(lambda) eq 0) then $
  k_load_vmatrix, vmatrix, lambda, vfile=vfile, vpath=vpath, lfile=lfile

spec=vmatrix#coeffs
nl=n_elements(lambda)-1L
lout=0.5*(lambda[0:nl-1L]+lambda[1:nl])
splot, lout, spec, _EXTRA=extra_for_splot

end
