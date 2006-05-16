;+
; NAME:
;   k_imf
; COMMENTS:
;   return IMF as a grid of M bins with relative # in each bin
; REVISION HISTORY:
;   3-May-2006  Michael Blanton (NYU)
;-
;------------------------------------------------------------------------------
pro k_imf, name=name, mgrid=mgrid, imf=imf

num=100
mgrid=10.^(alog10(1.e-1)+findgen(num)*(alog10(125.)-alog10(1.e-1))/ $
           float(num-1L))

;; From Table 2, left column, Chabrier (2003) PASP, 115, 763
if(name eq 'chabrier03') then begin
    mnorm=1.
    al=0.158  ;; disk, +-0.05
    mc=0.079  ;; disk, +-0.02
    sigma=0.69 ;; disk +-0.05
    ap=4.4e-2
    xx=1.3 ;; +- 1.3
    imf=fltarr(num)
    ii=where(mgrid lt mnorm, nii)
    if(nii gt 0) then $
      imf[ii]=al*exp(-(alog10(mgrid[ii])-alog10(mc))^2/(2.*sigma^2))/mgrid[ii]
    ii=where(mgrid ge mnorm, nii)
    if(nii gt 0) then $
      imf[ii]=ap*(mgrid[ii])^(-1.-xx)
endif

;; Salpeter (1955)
if(name eq 'salpeter55') then begin
    xx=1.35
    ap=1.
    imf=fltarr(num)
    imf=ap*(mgrid)^(-1.-xx)
endif

;; Diet Salpeter v1
if(name eq 'diet1') then begin
    xx=1.35
    ap=1.
    imf=fltarr(num)
    ii=where(mgrid gt 0.35, nii)
    if(nii gt 0) then $
      imf[ii]=ap*(mgrid[ii])^(-1.-xx)
    ii=where(mgrid le 0.35, nii)
    if(nii gt 0) then $
      imf[ii]=ap*(mgrid[ii])^(-1.)*(0.35^(-xx))
endif

;; normalize 
norm=total(mgrid*imf)
imf=imf/norm
print,norm

end
