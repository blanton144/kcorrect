pro compare_maggies,coeffs,redshift,maggies,maggies_ivar, $
                    rmatrix=rmatrix,zvals=zvals,recmaggies=recmaggies, $
                    filename=filename,vmatrix=vmatrix,lambda=lambda, $
                    filterlist=filterlist,filterpath=filterpath, $
                    psym=psym,xrange=xrange

if(NOT keyword_set(psym)) then psym=3

k_reconstruct_maggies,coeffs,redshift,recmaggies, $
  rmatrix=rmatrix,zvals=zvals,vmatrix=vmatrix,lambda=lambda, $
  filterlist=filterlist,filterpath=filterpath

nk=n_elements(filterlist)

set_plot,'ps'
device,filename=filename,xsize=5.5,ysize=5.5,/inches
hogg_plot_defaults, axis_char_scale=axis_char_scale,xold=xold,yold=yold, $
  pold=pold,default_font=default_font
!P.MULTI=[0,1,nk]
for i=0,nk-1 do $
  plot,redshift,-2.5*alog10(recmaggies[i,*]/maggies[i,*]),yra=[-0.99,0.99], $
  psym=psym,xra=xrange
!P.MULTI=[0,1,nk-1]
for i=0,nk-2 do $
  plot,redshift,-2.5*alog10(recmaggies[i,*]/maggies[i,*])+ $
  2.5*alog10(recmaggies[i+1,*]/maggies[i+1,*]),yra=[-0.99,0.99],psym=psym, $
  xra=xrange
device,/close
set_plot,'x'

end
