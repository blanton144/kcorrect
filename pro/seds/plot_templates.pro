;+
; NAME:
;   plot_templates
; PURPOSE:
;   plot the results of k_make_templates
; CALLING SEQUENCE:
;   plot_templates
; INPUTS:
; OPTIONAL INPUTS:
; KEYWORDS:
; OUTPUTS:
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   28-Nov-2003  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro plot_templates

restore,'data_make_templates.sav'

attime=0.
ntsteps=2000L
time=exp((alog(1.e+5)+(alog(20.e+9)-alog(1.e+5))* $
          (findgen(ntsteps))/float(ntsteps-1L)))
set_print,filename='k_sfh.ps'
!P.MULTI=[2,1,2]
!Y.MARGIN=[5.7]
!X.OMARGIN=[13.7]
djs_plot,[0],[0],/nodata,xra=[0.0005,20.],yra=[0.,1.1],/xlog, $
  xtitle='!6time (Gyrs)', ytitle='!6SFR'
for a = 0, n_elements(sfhpars)-1L do begin
    contrib=1.e-0*exp(-0.5*((sfhpars[a].age-attime-time)/ $
                            sfhpars[a].agesigma)^2)
    if(keyword_set(maxage)) then begin
        indx=where(time gt maxage-attime,count)
        if(count gt 0) then contrib[indx]=0.d
    endif
    if(keyword_set(minage)) then begin
        indx=where(time lt minage-attime,count)
        if(count gt 0) then contrib[indx]=0.d
    endif

    djs_oplot,time/1.e+9,contrib
endfor

djs_plot,[0],[0],/nodata,xra=[2000.,25000.], yra=[0.009,1040.5], $
  ytitle='!8\lambda f(\lambda)!6', $
  xtitle='!8\lambda!6 (Angstroms)',/ylog, /xlog
for a = 1, n_elements(sfhpars)-1L,2 do begin
    vindx=a+n_elements(sfhpars)*n_elements(dust)*2
    djs_oplot,lambda, $
      lambda*vmatrix[*,vindx]/(lambda[4000]*vmatrix[4000,vindx])
endfor



end_print

end
