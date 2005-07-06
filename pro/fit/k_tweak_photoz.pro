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
function k_tweak_photoz, vname=vname

sample='drtwo14'
postcat=hogg_mrdfits(vagc_name('post_catalog', sample=sample, letter='safe', $
                               post='1'), 1, nrow=28800)
postcat=postcat[where((postcat.letter_mask and 32) eq 0)]
indx_photo=shuffle_indx(n_elements(postcat), num_sub=10000, seed=seed)
postcat=postcat[indx_photo]
im=mrdfits(vagc_name('object_sdss_imaging'),1,row=postcat.object_position)
kc=sdss_kcorrect(postcat.z, calibobj=im, omaggies=maggies, oivar=ivar, $
                 coeffs=coeffs, vname=vname, vmatrix=vmatrix, lambda=lambda)

k_load_vmatrix, vmatrix, lambda, vname=vname
k_tweak_templates,maggies,ivar, postcat.z, coeffs, vmatrix, lambda, $
  vmatrix_tweaked=vmatrix_tweaked, maggies_factor=maggies_factor

end
;------------------------------------------------------------------------------
