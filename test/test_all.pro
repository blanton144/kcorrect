pro test_all

; read in the test sample
testsp=mrdfits('testsp.fits',1)

; use model mags, fix errors
errband=[0.05,0.02,0.02,0.02,0.03]
mags=testsp.counts_model
magserr=mags-mags
for i = 0, 4 do magserr[i]=sqrt(testsp.counts_model[i]^2+errband[i]^2)

; K-correct the objects
kcorrect,mags,magserr,testsp.z,rec_mags,kcorrectz=testsp.z,coeff=coeffs_real, $
  /sdssfix,/vconstraint,/returnmag
kcorrect,mags,magserr,testsp.z,kcorrects,kcorrectz=0.1,coeff=coeffs_real, $
  /sdssfix,/vconstraint

; Get photo-z's
kphotoz,mags,magserr,photoz,coeff=coeffs_photoz,/sdssfix

set_plot, "PS"
xsize= 7.5 & ysize= 7.5
!P.FONT= -1 & !P.BACKGROUND= 255 & !P.COLOR= 0
!P.THICK= 2.0
!P.CHARTHICK= !P.THICK & !X.THICK= !P.THICK & !Y.THICK= !P.THICK
!P.CHARSIZE= 1.2
axis_char_scale= 1.0
tiny= 1.d-4
!P.PSYM= 0
!P.TITLE= ''
!X.STYLE= 1
!X.CHARSIZE= axis_char_scale
!X.MARGIN= [1,1]*0.5*axis_char_scale
!X.OMARGIN= [6,6]*axis_char_scale
!X.RANGE= 0
!Y.STYLE= 1
!Y.CHARSIZE= !X.CHARSIZE
!Y.MARGIN= 0.6*!X.MARGIN
!Y.OMARGIN= 0.6*!X.OMARGIN
!Y.RANGE= 0
xyouts, 0,0,'!3'

device, file='test_all1.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color,/encapsulated
!p.multi=[5,1,5]
for i = 0, 4 do begin
    mn=djs_avsigclip(rec_mags[i,*]-mags[i,*],sigrej=5)
    sig=djsig(rec_mags[i,*]-mags[i,*],sigrej=5)
    djs_plot,testsp.z,rec_mags[i,*]-mags[i,*],xra=[0.,0.5], $
      yra=[mn,mn]+3.*[-sig,sig],xst=1,yst=1, psym=3
endfor
device,/close

device, file='test_all2.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color,/encapsulated
!p.multi=[5,1,5]
for i = 0, 4 do begin
    mn=djs_avsigclip(kcorrects[i,*],sigrej=5)
    sig=djsig(kcorrects[i,*],sigrej=5)
    djs_plot,testsp.z,kcorrects[i,*],xra=[0.,0.5], $
      yra=[mn,mn]+3.*[-sig,sig],xst=1,yst=1, psym=3
endfor
device,/close

stop

end
