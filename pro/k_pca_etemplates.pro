;+
; NAME:
;   k_pca_etemplates
;
; PURPOSE:
;   PCA eigentemplates for SED fitting
;
; CALLING SEQUENCE:
;   k_pca_etemplates, coeff, ematrix, bflux
;
; INPUTS:
;   coeff  - coefficients to eigenspectrum
;   ematrix  - eigenspectra
;   bflux    - flux in each SED template
;   
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL INPUT/OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   05-Jan-2002  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_pca_etemplates, coeff, ematrix, bflux, nclip=nclip, niter=niter

; Need at least 2 parameters
if (N_params() LT 3) then begin
    klog, 'Syntax - k_pca_etemplates, coeff, ematrix, bflux'
    return 
endif

if (NOT keyword_set(nclip)) then nclip=5.d
if (NOT keyword_set(niter)) then niter=5l

nb=long(n_elements(bflux))
nt=long(n_elements(ematrix))/nb
ngalaxy=long(n_elements(coeff))/nt

if(nt le 1) then begin
    klog,'must have > 1 template'
    return
endif

usesample=lindgen(ngalaxy)
ascale=dblarr(nt-1,ngalaxy)
amean=dblarr(nt-1,ngalaxy)
for j = 1,nt-1 do begin
    ascale[j-1,*]=coeff[j,*]/coeff[0l,*]
endfor
for clip=0l, niter-1l do begin
    for j = 1,nt-1 do begin
        amean[j-1,usesample]=total(ascale[j-1,usesample],/double) $
          /double(n_elements(usesample))
    endfor
    covar=dblarr(nt-1l,nt-1l)
		for i=0, nt-2l do $	
		    for j=0, nt-2l do $	
				   covar[i,j]=total((ascale[i,usesample]-amean[i])* $
														(ascale[j,usesample]-amean[j]),/double)
    covar=covar/n_elements(usesample)
    i=lindgen(nt-1l)
    cliplimit=nclip^2*(total(covar(i,i),/double)/double(nt-1l))
    klog,cliplimit
    usesample=where(total((ascale[*,usesample]-amean[*,usesample])^2,1, $
                           /double) lt cliplimit,count)
    if(count eq 0) then begin
        klog,'all objects clipped! resetting.'
        usesample=lindgen(ngalaxy)
    endif
endfor 

; evals returned, evecs in rows of covar
i=lindgen(nt-1l)
trired,covar,eval, tmp, /double
triql,eval,tmp,covar

; rotate
ascale=transpose(transpose(ascale)#covar)
for i = 1l, nt-1l do begin
    coeff[i,*]=coeff[0,*]*ascale[i-1,*]
endfor
for i = 1l, nb-1l do begin
    ematrix[i,1:nt-1]=covar#transpose(ematrix[i,1:nt-1])
endfor

end 
