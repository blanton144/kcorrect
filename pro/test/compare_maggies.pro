pro compare_maggies,coeffs,redshift,maggies,maggies_ivar, $
                    rmatrix=rmatrix,zvals=zvals,recmaggies=recmaggies, $
                    filename=filename,vmatrix=vmatrix,lambda=lambda, $
                    filterlist=filterlist,filterpath=filterpath

k_reconstruct_maggies,coeffs,redshift,recmaggies, $
  rmatrix=rmatrix,zvals=zvals,vmatrix=vmatrix,lambda=lambda, $
  filterlist=filterlist,filterpath=filterpath

set_plot,'ps'
device,filename=filename,xsize=5.5,ysize=5.5,/inches
hogg_plot_defaults, axis_char_scale=axis_char_scale,xold=xold,yold=yold, $
  pold=pold,default_font=default_font
!P.MULTI=[0,1,5]
for i=0,4 do $
  plot,redshift,-2.5*alog10(recmaggies[i,*]/maggies[i,*]),yra=[-0.39,-0.39], $
  psym=3
!P.MULTI=[0,1,4]
for i=0,3 do $
  plot,redshift,-2.5*alog10(recmaggies[i,*]/maggies[i,*])+ $
  2.5*alog10(recmaggies[i+1,*]/maggies[i+1,*]),yra=[-0.39,0.39],psym=3
device,/close
set_plot,'x'

end
