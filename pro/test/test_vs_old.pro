pro test_vs_old

restore,'testdata*.sav'
band_shift=0.1

create_old,maggies[0:4,0:n_elements(gals)-1L], $
  maggies_ivar[0:4,0:n_elements(gals)-1L], $
  gals.redshift,kcorrect_old, band_shift=band_shift

kcorrect, maggies[0:4,0:n_elements(gals)-1L], $
  maggies_ivar[0:4,0:n_elements(gals)-1L], $
  gals.redshift,kcorrect, band_shift=band_shift

set_print,filename='test_vs_old.ps'

yrangelist=[ $
             [-0.8,3.], $
             [-0.5,2.], $
             [-0.3,1.], $
             [-0.3,0.5], $
             [-0.3,0.4] ]

for i=0, 4 do begin & $
  yrange=yrangelist[*,i] & $
  !P.MULTI=[0,1,2] & $
  plot,gals.redshift,2.5*alog10(kcorrect[i,*]),psym=3, $
  yra=yrange,xra=[0.,0.5] & $
  oplot,gals.redshift,2.5*alog10(kcorrect_old[i,*]),psym=3, $
  color=djs_icolor('red') & $
  plot,gals.redshift, $
  2.5*alog10(kcorrect[i,*])-2.5*alog10(kcorrect_old[i,*]) ,psym=3, $
  yra=[-0.3,0.3],xra=[0.,0.5] & $
  endfor

end_print

end
