; make demos for web page
pro mk_demo

plate=455
fiberids=139+[105,300,104,114]

k_print,filename='demo.ps',axis_char_scale=axis_char_scale
!P.MULTI=[0,2,2]
xcharsize=[0.0001,0.0001,0.7*axis_char_scale,0.7*axis_char_scale]
ycharsize=[0.7*axis_char_scale,0.0001,0.7*axis_char_scale,0.0001]
for ifiberid=0, n_elements(fiberids)-1L do begin
    istr=strtrim(string(ifiberid),2)
    fiberid=fiberids[ifiberid]
    readspec,plate,fiberid,zans=zans,flux=flux,wave=wave
    flux=smooth(flux,20)
    kcorrect,zans.counts_spectro[1:3],zans.counts_spectro[1:3]*0.1,zans.z, $
      kcorrect,coeffs=coeffs,/stddev,vmatrix=vmatrix,lambda=lambda, $
      filterlist=['sdss_g0.par','sdss_r0.par','sdss_i0.par']
    k_load_filters,['sdss_g0.par','sdss_r0.par','sdss_i0.par'], filter_n, $
      filter_lambda, filter_pass
    spec=vmatrix#coeffs*1.e+17
    flimits=[-0.05,1.15]*max(flux)
    wlimits=[0.92*min(wave),0.95*max(wave)]
    !X.CHARSIZE=xcharsize[ifiberid]
    !Y.CHARSIZE=ycharsize[ifiberid]
    djs_plot,wave,flux,thick=2,xra=wlimits,yra=flimits, $
      xtitle='\lambda (\AA)', ytitle='f(\lambda) (10^{-17} ergs/cm^2/s/\AA)'
    djs_oplot,lambda*(1.+zans.z),spec/(1.+zans.z),color=djs_icolor('red')
    colors=['blue','green','red']
    band=['g','r','i']
    for i=0L, n_elements(colors)-1L do $
      djs_oplot,filter_lambda[0:filter_n[i]-1,i], $
      filter_pass[0:filter_n[i]-1,i]*max(flux)*0.8, color=colors[i], thick=5
    for i=0L, n_elements(colors)-1L do $
      xyouts,filter_lambda[filter_n[i]/2L,i], $
      filter_pass[filter_n[i]/2L,i]*max(flux)*0.4, $
      band[i],color=djs_icolor(colors[i]),charsize=2.,align=0.5
    xyouts,!X.CRANGE[0]+0.05*(!X.CRANGE[1]-!X.CRANGE[0]), $
      !Y.CRANGE[0]+0.9*(!Y.CRANGE[1]-!Y.CRANGE[0]), $
      'z='+strtrim(string(zans.z,format='(f40.2)'),2), $
      charsize=1.7
endfor
k_end_print
spawn,'psgif -r 90 demo.ps > demo.gif'

end
