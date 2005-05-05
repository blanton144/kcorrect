;+
; NAME:
;   k_plot_bc_colors
; PURPOSE:
;   plot colors vs time for various models
; CALLING SEQUENCE:
;   k_plot_bc_colors
; REVISION HISTORY:
;   02-May-2005  MRB, NYU
;-
;------------------------------------------------------------------------------
pro k_plot_bc_colors

salpeter=0
isolib='Padova1994'
if(NOT keyword_set(nagesmax)) then nagesmax=40
norm_lmin=800. ;; limits when testing whether two BC03 models are similar
norm_lmax=23000.

bc03= k_im_read_bc03(salpeter=salpeter, isolib=isolib)
imaxage=max(where(bc03.age lt 14.e+9))
iwave=where(bc03.wave gt norm_lmin and bc03.wave lt norm_lmax, nwave)
nages=3000L + nagesmax
tol=5.
while(nages gt nagesmax) do begin
    iuse=imaxage 
    i=iuse
    while i ge 1 do begin
        j=iuse[0]-1L
        scalefact=total(bc03.flux[iwave,i]*bc03.flux[iwave,j])/ $
          total(bc03.flux[iwave,i]*bc03.flux[iwave,i])
        diff=total(((bc03.flux[iwave,i]*scalefact- $
                     bc03.flux[iwave,j])/ $
                    bc03.flux[iwave,j])^2,/double) 
        while(diff lt tol AND j gt 0L) do begin
            j=j-1L
            scalefact=total(bc03.flux[iwave,i]*bc03.flux[iwave,j])/ $
              total(bc03.flux[iwave,i]*bc03.flux[iwave,i])
            diff=total(((bc03.flux[iwave,i]*scalefact- $
                         bc03.flux[iwave,j])/ $
                        bc03.flux[iwave,j])^2,/double) 
        endwhile
        if(j ge 0 and diff ge tol) then begin
            iuse=[j,iuse]
        endif
        i=j
    endwhile
    tol=tol/0.85
    ages=bc03.age[iuse]
    nages=n_elements(ages)
    help,tol, nages
endwhile
nuse=n_elements(iuse)

;;     b. now make the full grid for the desired ages
tmp_bc03= k_im_read_bc03(age=1., isolib=isolib, salpeter=salpeter)
nl=n_elements(tmp_bc03.flux)
lambda=k_lambda_to_edges(tmp_bc03.wave)
nmets=6   ;; metallicities to use
mets=[0,1,2,3,4,5]
umg=fltarr(nuse, nmets)
gmr=fltarr(nuse, nmets)
for im= 0L, nmets-1L do begin
    flux= (k_im_read_bc03(met=mets[im], isolib=isolib, salpeter=salpeter)).flux[*,iuse]
    for it=0, nuse-1L do begin
        maggies=k_project_filters(lambda, flux[*,it], band_shift=0.1, $
                                  filterlist=['sdss_u0.par', 'sdss_g0.par', $
                                              'sdss_r0.par'])
        umg[it,im]=-2.5*alog10(maggies[0]/maggies[1])
        gmr[it,im]=-2.5*alog10(maggies[1]/maggies[2])
    endfor
endfor

postcat=hogg_mrdfits(vagc_name('post_catalog', sample='sample15', $
                               letter='bsafe', post='1'), 1, nrow=28800)
ii=where(postcat.z gt 0.07 and postcat.z lt 0.13)
postcat=postcat[ii]
postcat=postcat[shuffle_indx(n_elements(postcat), num_sub=10000)]
im=mrdfits(vagc_name('object_sdss_imaging'),1,row=postcat.object_position)
kc=sdss_kcorrect(postcat.z, calibobj=im, absmag=mabsm, band_shift=0.1, $
                 flux='model')
kc=sdss_kcorrect(postcat.z, calibobj=im, absmag=pabsm, band_shift=0.1, $
                 flux='petro')

k_print, filename='bc_umg.ps'
djs_plot,[0],[0], /nodata, /xlog, xra=[500000., 1.4e+10], $
  yra=[-0.5, 2.8]
for im= 0L, nmets-1L do $
  djs_oplot,ages, umg[*,im], /xlog
k_end_print

k_print, filename='bc_umg_gmr.ps'
djs_plot,[0],[0], /nodata, xra=[-0.5, 2.8], yra=[-0.5,1.5]
djs_oplot, pabsm[0,*]-pabsm[1,*], $
  pabsm[1,*]-pabsm[2,*], psym=3, color='red'
djs_oplot, mabsm[0,*]-mabsm[1,*], $
  mabsm[1,*]-mabsm[2,*], psym=3, color='green'
for im= 0L, nmets-1L do $
  djs_oplot,umg[*,im], gmr[*,im]

k_end_print

stop

end
;------------------------------------------------------------------------------
