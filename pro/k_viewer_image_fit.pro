;+
; NAME:
;   k_viewer_image_fit
; PURPOSE:
;   Fit components to viewer image
; CALLING SEQUENCE:
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
;   25-jan-2002  WRitten by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_viewer_image_fit, filename, red, blue

; get id
fullid=long64((stregex(filename,'image(.*)\.fit',/subexpr,/extract))[1])
id=fullid mod 10000
field=(fullid/10000) mod 10000
camcol=(fullid/10000/10000) mod 100  
rerun=(fullid/10000/10000/100) mod 10000
run=(fullid/10000/10000/100/10000)    

; template information
k_load_ascii_table,vmatrix,'vmatrix.aelines.dat'
k_load_ascii_table,lambda,'lambda.aelines.dat'
nv=n_elements(vmatrix)/(n_elements(lambda)-1L)

; read in images
dummy=mrdfits(filename,0,hdr)
softbias=float(sxpar(hdr,'SOFTBIAS'))
im_g=float(mrdfits(filename,1))-softbias
im_r=float(mrdfits(filename,2))-softbias
im_i=float(mrdfits(filename,3))-softbias
nx=(size(im_g,/dimensions))[0]
ny=(size(im_g,/dimensions))[1]

; calibrate images
tsfieldfile=rdss_name('tsField',run,rerun,camcol,field)
tsfield=mrdfits(tsfieldfile+'*',1)
maggies=fltarr(3,n_elements(im_g))
maggies_err=fltarr(3,n_elements(im_g))+1. ; assume g r and i same S/N
maggies[0,*]=im_g/55.247* $
  10.^(0.4*(tsfield.aa[1]+tsfield.airmass[1]*tsfield.kk[1]))
maggies[1,*]=im_r/55.247* $
  10.^(0.4*(tsfield.aa[2]+tsfield.airmass[2]*tsfield.kk[2]))
maggies[2,*]=im_i/55.247* $
  10.^(0.4*(tsfield.aa[3]+tsfield.airmass[3]*tsfield.kk[3]))
for i=0,2 do maggies[i,*]=maggies[i,*]-median(maggies[i,*])

; now fit
k_fit_coeffs,maggies,1./maggies_err,replicate(0.03,n_elements(im_g)), $
  coeffs,bmatrix=vmatrix,lambda=lambda, $
  filterlist=['sdss_g0','sdss_r0','sdss_i0']+'.dat',ematrix=identity(nv)

red=reform(coeffs[0,*],nx,ny)
blue=reform(coeffs[1,*],nx,ny)

end
;------------------------------------------------------------------------------
