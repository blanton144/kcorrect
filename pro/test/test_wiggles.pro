k_load_vmatrix, vmatrix, lambda
vmatrix=vmatrix[*,0]
vmatrix=smooth(vmatrix,20)
coeffs=replicate(1.,1000)
redshift=findgen(1000)/1000.
k_reconstruct_maggies, coeffs, redshift, rec, vmatrix=vmatrix, lambda=lambda, $
  filterlist='sdss_i0.par'
