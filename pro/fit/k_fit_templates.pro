;+
; NAME:
;   k_fit_templates
; PURPOSE:
;   run the template fitting code, with priming
; CALLING SEQUENCE:
;   k_fit_templates [, nt= ]
; OPTIONAL INPUTS:
;   nt - number of templates to fit for (default 6)
; COMMENTS:
;   Requires k_nmf_mmatrix and k_nmf_spdata to have been run. 
;   Runs k_run_nmf multiple times:
;     1. runs short set of iterations to get basic fit 
;     2. classifies galaxies based on optical (model) color into nt sets
;     3. runs k_run_nmf nt times, once on each set, to get starting points
;     4. then runs k_run_nmf over and over again (up to 100 times) with
;        1000 iterations each time
; REVISION HISTORY:
;   09-Apr-2005  Michael Blanton (NYU)
;-
;------------------------------------------------------------------------------
pro k_fit_templates, nt=nt

if(NOT keyword_set(nprime)) then nprime=2000
if(NOT keyword_set(subprime)) then subprime=1.

;; 1. initial fits
k_run_nmf, nt=nt, niter=100, /reset
k_run_nmf, nt=nt, niter=nprime

;; 2. split sample

;;     a. read in data
hdr=headfits('k_nmf_mmatrix.fits')
nspec=long(sxpar(hdr,'NSPEC'))
nzf=long(sxpar(hdr,'NZ'))
mmatrix=mrdfits('k_nmf_mmatrix.fits',0,hdr)
zf=mrdfits('k_nmf_mmatrix.fits',6)
datastr=mrdfits('k_nmf_spdata.fits',1)
vals=mrdfits('k_nmf_spdata.fits',2)
ivar=mrdfits('k_nmf_spdata.fits',3)
xx=mrdfits('k_nmf_spdata.fits',4)
zhelio=mrdfits('k_nmf_spdata.fits',7)
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

;;     b. read in the results
filterlist=mrdfits('k_nmf_mmatrix.fits',5)
templates=mrdfits('k_nmf_soln.fits')
coeffs=mrdfits('k_nmf_soln.fits',1)
nt=(size(templates, /dim))[1]

;;     c. get u-g for model
model={val:fltarr(2L*data.ny), $
       x:fltarr(2L*data.ny), $
       nx:data.nx, $
       ny:data.ny, $
       rowstart:lonarr(data.ny), $
       nxrow:lonarr(data.ny)}
iu=where(strtrim(string(filterlist),2) eq 'sdss_u0.par')
ir=where(strtrim(string(filterlist),2) eq 'sdss_r0.par')
xu=nspec+iu[0]*nzf+0L
xr=nspec+ir[0]*nzf+0L
for i=0L, model.ny-1L do begin
    model.rowstart[i]=i*2L
    model.nxrow[i]=2L
    model.x[model.rowstart[i]]=xu
    model.x[model.rowstart[i]+1L]=xr
endfor
mcoeffs=templates#coeffs
mmeval, model, transpose(mmatrix), mcoeffs
umag=-2.5*alog10(model.val[model.rowstart])
rmag=-2.5*alog10(model.val[model.rowstart+1L])
umr=umag-rmag

;;       c. classify by umr
isort=sort(umr)
ipos=lonarr(n_elements(umr))
ipos[isort]=lindgen(n_elements(umr))
iclass=long(float(nt)*float(ipos)/float(n_elements(umr)))

save, filename='prime.sav'

;; 3. run nmf fitting for each separately
for i=0L, nt-1L do begin
    spawn, 'mkdir -p prime'+strtrim(string(i),2)
    cd, 'prime'+strtrim(string(i),2)
    spawn, 'cp ../k_nmf_mmatrix.fits .'
    spawn, 'cp ../k_nmf_mmatrix.fits .'

    hdr=headfits('../k_nmf_spdata.fits')
    datastrfull=mrdfits('../k_nmf_spdata.fits',1)
    valsfull=mrdfits('../k_nmf_spdata.fits',2)
    ivarfull=mrdfits('../k_nmf_spdata.fits',3)
    xxfull=mrdfits('../k_nmf_spdata.fits',4)
    lambda=mrdfits('../k_nmf_spdata.fits',5)
    zdistfull=mrdfits('../k_nmf_spdata.fits',6)
    zheliofull=mrdfits('../k_nmf_spdata.fits',7)
    
    iuse=where(iclass eq i, nuse)
    nusenew=long(nuse*subprime)
    iuse=iuse[shuffle_indx(nuse,num_sub=nusenew)]
    nuse=n_elements(iuse)
    datastr={nx:datastrfull.nx, $
             ny:nuse, $
             rowstart:lonarr(nuse), $
             nxrow:lonarr(nuse)}
    vals=fltarr(long(total(datastrfull.nxrow[iuse])))
    ivar=fltarr(long(total(datastrfull.nxrow[iuse])))
    xx=lonarr(long(total(datastrfull.nxrow[iuse])))
    ix=0L
    for j=0L, nuse-1L do begin
        vals[ix:ix+datastrfull.nxrow[iuse[j]]-1]= $
          valsfull[datastrfull.rowstart[iuse[j]]: $
                   datastrfull.rowstart[iuse[j]]+datastrfull.nxrow[iuse[j]]-1]
        ivar[ix:ix+datastrfull.nxrow[iuse[j]]-1]= $
          ivarfull[datastrfull.rowstart[iuse[j]]: $
                   datastrfull.rowstart[iuse[j]]+datastrfull.nxrow[iuse[j]]-1]
        xx[ix:ix+datastrfull.nxrow[iuse[j]]-1]= $
          xxfull[datastrfull.rowstart[iuse[j]]: $
                 datastrfull.rowstart[iuse[j]]+datastrfull.nxrow[iuse[j]]-1]
        datastr.rowstart[j]=ix
        datastr.nxrow[j]=datastrfull.nxrow[iuse[j]]
        ix=ix+datastrfull.nxrow[iuse[j]]
    endfor
    zdist=zdistfull[iuse]
    zhelio=zheliofull[iuse]
    outfile='k_nmf_spdata.fits'
    mwrfits, 0, outfile, hdr, /create
    mwrfits, datastr, outfile
    mwrfits, vals, outfile
    mwrfits, ivar, outfile
    mwrfits, xx, outfile
    mwrfits, lambda, outfile
    mwrfits, zdist, outfile
    mwrfits, zhelio, outfile

    k_run_nmf, nt=1, niter=nprime, /reset, /qa
    cd, '../'
endfor

alltemplates=mrdfits('prime0/k_nmf_soln.fits')
for i=1L, nt-1L do $
  alltemplates=[alltemplates, mrdfits('prime'+strtrim(string(i),2)+ $
                                      '/k_nmf_soln.fits')]
alltemplates=reform(alltemplates, n_elements(alltemplates)/nt, nt)
mwrfits, alltemplates, 'k_nmf_soln.fits', /create

k_run_nmf, nt=nt, niter=nprime, /qa

for i=0L, 99L do $
  k_run_nmf, nt=nt, niter=nstep, /qa

end
