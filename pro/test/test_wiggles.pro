k_load_vmatrix, vmatrix, lambda
vmatrix=vmatrix[*,0]
vmatrix=smooth(vmatrix,20)
coeffs=replicate(1.,1000)
redshift=findgen(1000)/1000.
k_reconstruct_maggies, coeffs, redshift, rec, vmatrix=vmatrix, lambda=lambda, $
  filterlist='sdss_i0.par'
set_print,filename='cnoc2_stack.ps'
plot,cnoc2.z,umg2,psym=4,xra=[0.,0.7],yra=[-0.5,3.],xtitle='z',ytitle='u-g'
plot,cnoc2.z,umg,psym=4,xra=[0.,0.7] ,yra=[-0.5,3.],xtitle='z',ytitle='u-g (stacked)'
plot,cnoc2.z,gmr2,psym=4,xra=[0.,0.7],yra=[-0.1,2.2],xtitle='z',ytitle='g-r'
plot,cnoc2.z,gmr,psym=4,xra=[0.,0.7] ,yra=[-0.1,2.2],xtitle='z',ytitle='g-r (stacked)'
plot,cnoc2.z,rmi2,psym=4,xra=[0.,0.7],yra=[-0.3,1.1],xtitle='z',ytitle='r-i'
plot,cnoc2.z,rmi,psym=4,xra=[0.,0.7] ,yra=[-0.3,1.1],xtitle='z',ytitle='r-i (stacked)'
plot,cnoc2.z,imz2,psym=4,xra=[0.,0.7],yra=[-0.7,1.1],xtitle='z',ytitle='i-z'
plot,cnoc2.z,imz,psym=4,xra=[0.,0.7] ,yra=[-0.7,1.1],xtitle='z',ytitle='i-z (stacked)'
end_print
