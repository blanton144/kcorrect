pro k_print,filename=filename,axis_char_scale=axis_char_scale, $
            tiny=tiny,pold=pold,xold=xold,yold=yold,colorname=colorname, $
            xsize=xsize,ysize=ysize

if(NOT keyword_set(axis_char_scale)) then axis_char_scale= 1.75
if(NOT keyword_set(tiny)) then tiny=1.d-4

; setup postscript file
pold=!P
xold=!X
yold=!Y
!P.FONT= -1
set_plot, "PS"
!P.BACKGROUND= djs_icolor('white')
!P.COLOR= djs_icolor('black')
if(NOT keyword_set(xsize)) then xsize= 7.5 
if(NOT keyword_set(ysize)) then ysize= 7.5 
device, file=filename,/inches,xsize=xsize,ysize=ysize, $
  xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color, $
  bits_per_pixel=64
!P.THICK= 2.0
!P.CHARTHICK= !P.THICK & !X.THICK= !P.THICK & !Y.THICK= !P.THICK
!P.CHARSIZE= 1.0
!P.PSYM= 0
!P.LINESTYLE= 0
!P.TITLE= ''
!X.STYLE= 1
!X.CHARSIZE= axis_char_scale
!X.MARGIN= [1,1]*0.5
!X.OMARGIN= [7,7]*axis_char_scale
!X.RANGE= 0
!X.TICKS= 0
!Y.STYLE= 1
!Y.CHARSIZE= !X.CHARSIZE
!Y.MARGIN= 0.6*!X.MARGIN
!Y.OMARGIN= 0.6*!X.OMARGIN
!Y.RANGE= 0
!Y.TICKS= !X.TICKS
!P.MULTI= [1,1,1]
xyouts, 0,0,'!6'
colorname= ['red','green','blue','magenta','cyan','dark yellow', $
            'purple','light green','orange','navy','light magenta', $
            'yellow green']
ncolor= n_elements(colorname)
loadct,0

end
