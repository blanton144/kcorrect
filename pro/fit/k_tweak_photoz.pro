;+
; NAME:
;   k_tweak_photoz
; PURPOSE:
;   Read in data and use it to tweak the photo-z templates 
; CALLING SEQUENCE:
;   k_tweak_photoz, vname=vname
; COMMENTS:
;   currently experimental
; REVISION HISTORY:
;   01-May-2003  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_tweak_photoz, vname=vname

datastr=mrdfits('k_nmf_spdata.fits',1)
vals=mrdfits('k_nmf_spdata.fits',2)
ivar=mrdfits('k_nmf_spdata.fits',3)
xx=mrdfits('k_nmf_spdata.fits',4)
hdr=headfits('k_nmf_mmatrix.fits')
nzf=long(sxpar(hdr, 'NZ'))
nspec=long(sxpar(hdr, 'NSPEC'))
zf=mrdfits('k_nmf_mmatrix.fits',6)
filterlist=strtrim(string(mrdfits('k_nmf_mmatrix.fits',5)),2)
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

maggies=fltarr(n_elements(filterlist), n_elements(zhelio))
ivar=fltarr(n_elements(filterlist), n_elements(zhelio))

for ifilter=0L, n_elements(filterlist)-1L do begin
    fstart=nspec+ifilter*nzf
    fend=nspec+(ifilter+1L)*nzf-1L
    for j=0L, n_elements(zhelio)-1L do begin
        currx=data.rowstart[j]+lindgen(data.nxrow[j])
        ii=where(data.x[currx] ge fstart AND data.x[currx] le fend, nii)
        if(nii gt 0) then begin
            maggies[ifilter, j]= data.val[currx[ii[0]]]
            ivar[ifilter, j]= data_ivar.val[currx[ii[0]]]
        endif
    endfor 
endfor

coeffs=mrdfits('k_nmf_soln.fits',1)
k_load_vmatrix, vmatrix, lambda, vname='test', vpath='.'
k_tweak_templates, maggies, ivar, zhelio, coeffs, vmatrix, lambda, $
  vmatrix_tweaked=vmatrix_tweaked, maggies_factor=maggies_factor, $
  filterlist=filterlist

k_write_ascii_table, vmatrix_tweaked, 'vmatrix.photoz.dat'
k_write_ascii_table, lambda, 'lambda.photoz.dat'

stop

end
;------------------------------------------------------------------------------
