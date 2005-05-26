pro swire_plot_bruce

hogg_usersym, 10, /fill

;; read in the template basics
mmatrix=mrdfits('k_nmf_mmatrix.fits',0,hdr)
lambda=mrdfits('k_nmf_mmatrix.fits',1)
dust=mrdfits('k_nmf_mmatrix.fits',2)
met=mrdfits('k_nmf_mmatrix.fits',3)
age=mrdfits('k_nmf_mmatrix.fits',4)
rawspec=mrdfits('k_nmf_rawspec.fits',0)
filterlist=string(mrdfits('k_nmf_mmatrix.fits',5))
early=mrdfits('k_nmf_early.fits',0,earlyhdr)
zf=mrdfits('k_nmf_mmatrix.fits',6)
nspec=long(sxpar(hdr, 'NSPEC'))
back=long(sxpar(hdr, 'BACK'))
nextra=long(sxpar(hdr, 'NEXTRA'))
nzf=long(sxpar(hdr, 'NZ'))
nfilter=long(sxpar(hdr, 'NFILTER'))
ndusts=long(sxpar(hdr, 'NDUST'))
nmets=long(sxpar(hdr, 'NMET'))
nages=long(sxpar(hdr, 'NAGE'))
filternames=['F', 'N', 'u', 'g', 'r', 'i', 'z', 'J', 'H', 'K_s', 'B', 'R', $
             'I', 'J', "H", 'K_s', 'B', 'V', 'i', 'z', '[3.6]', '[4.5]', $
             '[5.8]', '[8.0]', '[24]' ]
metallicities=[0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05]

;; read in the data
datastr=mrdfits('k_nmf_spdata.fits',1)
vals=mrdfits('k_nmf_spdata.fits',2)
ivar=mrdfits('k_nmf_spdata.fits',3)
xx=mrdfits('k_nmf_spdata.fits',4)
data=create_struct(datastr, $
                   'val', fltarr(n_elements(vals)), $
                   'x', fltarr(n_elements(vals)))
data.val=vals
data.x=xx
data_ivar=create_struct(datastr, $
                        'val', fltarr(n_elements(vals)), $
                        'x', fltarr(n_elements(vals)))
data_ivar.val=ivar
data_ivar.x=xx
zhelio=mrdfits('k_nmf_spdata.fits',7)
iz=long(floor((nzf-1.)*(zhelio-zf[0])/(zf[nzf-1]-zf[0])+0.5))

;; read in the results
templates=mrdfits('k_nmf_soln.fits')
coeffs=mrdfits('k_nmf_soln.fits',1)
nt=n_elements(coeffs)/n_elements(zhelio)
model=data
early=data
mcoeffs=templates#coeffs
mmeval, model, transpose(mmatrix), mcoeffs

;; make vmatrix and lambda
outvmatrix=mmatrix[0:nspec-1L,*]#templates
absrc=3.631*2.99792*1.e-2/(lambda[0:nspec-1L])^2
for i=0L, nt-1L do $
  outvmatrix[*,i]=outvmatrix[*,i]*absrc
outlambda=fltarr(nspec+1L)
dlg10l=alog10(lambda[1])-alog10(lambda[0])
outlambda[0:nspec-1L]= 10.^(alog10(lambda[0:nspec-1L])-0.5*dlg10l)
outlambda[nspec]= 10.^(alog10(lambda[nspec-1L])+0.5*dlg10l)

k_print, filename='swire_red.ps'

tspec=outvmatrix#templates
djs_plot, lambda[0:nspec-1], 1.e+31*tspec[*,i], /xlog, /ylog, $
  xra=[2001., 500000.], xtitle='\lambda (Angstroms)', $
  ytitle='f_\lambda'

ngals=n_elements(data.rowstart)
col=fltarr(4L,ngals)
ifilters=[20,21,22,23,24]
colerr=fltarr(nfilter-1L,ngals)
mcol=fltarr(nfilter-1L,ngals)

for i=0, 3 do begin
    f0start=nspec+ifilters[i]*nzf
    f0end=nspec+(ifilters[i]+1L)*nzf-1L
    f1start=nspec+(ifilters[i]+1L)*nzf
    f1end=nspec+(ifilters[i]+2L)*nzf-1L
    for j=0L, ngals-1L do begin
        currx=data.rowstart[j]+lindgen(data.nxrow[j])
        i0=where(data.x[currx] ge f0start AND data.x[currx] le f0end, n0)
        i1=where(data.x[currx] ge f1start AND data.x[currx] le f1end, n1)
        if(n1 gt 0 and n0 gt 0) then begin
            mcol[i,j]=-2.5*alog10(model.val[currx[i0[0]]]/ $
                                        model.val[currx[i1[0]]])
            col[i,j]=-2.5*alog10(data.val[currx[i0[0]]]/ $
                                       data.val[currx[i1[0]]])
            colerr[i,j]= $
              2.5/alog(10.)*sqrt(1./data_ivar.val[currx[i0[0]]]/ $
                                 data.val[currx[i0[0]]]^2+ $
                                 1./data_ivar.val[currx[i1[0]]]/ $
                                 data.val[currx[i1[0]]]^2)
        endif
    endfor
endfor
    
!P.MULTI=[0,1,4]
names=['[3.6]', '[4.5]', '[5.8]', '[8.0]', '[24]']
yra=[[-0.69, 0.29], $
     [-0.49, 0.69], $
     [-0.99, 2.49], $
     [-0.39, 2.79]]
for i=0, 3 do begin
    xcharsize=0.0001
    ycharsize=2.3
    if(i eq 3) then xcharsize=ycharsize
    igood=where(colerr[i,*] ne 0., ngood)
    djs_plot, zhelio[igood], col[i, igood], psym=8, symsize=0.5, $
      xcharsize=xcharsize, xtitle='z', $
      ycharsize=ycharsize, ytitle=names[i]+'-'+names[i+1], color='red', $
      xra=[0.01, 0.39], yra=yra[*,i]
    djs_oplot, zhelio[igood], mcol[i, igood], psym=8, symsize=0.5
endfor

!P.MULTI=[0,1,1]
igood=where(colerr[0,*] ne 0. and $
            colerr[1,*] ne 0., ngood)
djs_plot, col[0,igood], col[1,igood], psym=8, symsize=0.5, color='red', $
  ycharsize=ycharsize, ytitle=names[1]+'-'+names[1+1], $
  xcharsize=xcharsize, xtitle=names[0]+'-'+names[0+1], $
  xra=yra[*,0], yra=yra[*,1]
djs_oplot, mcol[0,igood], mcol[1,igood], psym=8, symsize=0.5

  
    
k_end_print

end
