pro test_qld,covar,x0

nconstraints=2L
neconstraints=0L
nconstraints_max=3L
nvar=2L
nvar_max=2L

npix=200
xlimits=[[-3.,3.],[-3.,3.]]
xx=dblarr(2,npix)
ximage=dblarr(2,npix,npix)

xx[0,*]=xlimits[0,0]+(xlimits[1,0]-xlimits[0,0])*(dindgen(npix)+0.5)/ $
  double(npix)
ximage[0,*,*]=xx[0,*]##replicate(1.,npix)
xx[1,*]=xlimits[0,1]+(xlimits[1,1]-xlimits[0,1])*(dindgen(npix)+0.5)/ $
  double(npix)
ximage[1,*,*]=transpose(xx[1,*]##replicate(1.,npix))

image=dblarr(npix,npix)
for i=0, npix-1 do $
  for j=0, npix-1 do $
  image[i,j]=0.5*ximage[*,i,j]#invert(covar)#ximage[*,i,j]- $
  x0#invert(covar)#ximage[*,i,j]+0.5*x0#invert(covar)#x0
image=transpose(image)
contour,image,xx[0,*],xx[1,*],levels=[0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.]

cmatrix=invert(covar)
dmatrix=-x0#invert(covar)
amatrix=dblarr(nconstraints_max,nvar)
amatrix[0,*]=[1.,0.]
amatrix[1,*]=[0.,1.]
bmatrix=dblarr(nconstraints_max)
xl=replicate(-1.d+30,nvar)
xu=replicate(1.d+30,nvar)
x=dblarr(nvar)
lagrange=dblarr(nconstraints+nvar*2)
ifail=-1L
iprint=-1L

; Set source object name
soname=filepath('libkcorrect.so', $
                root_dir=getenv('KCORRECT_DIR'), subdirectory='lib')

retval=call_external(soname, 'idl_k_qld', long(nconstraints), $
                     long(neconstraints), long(nconstraints_max), $
                     long(nvar), long(nvar_max), double(cmatrix), $
                     double(dmatrix), double(amatrix), double(bmatrix), $
                     double(xl), double(xu), x, lagrange, ifail, $
                     iprint)

splog,x
stop

end
