;+
; NAME:
;   k_load_templates
;
; PURPOSE:
;   Load the template information from a list of files, and 
;   interpolate them all onto the same wavelength scale
;
; CALLING SEQUENCE:
;   k_load_templates, templatelist, vmatrix, lambda
;
; INPUTS:
;   templatelist      - list of files with template information
;
; OPTIONAL INPUTS:
;   nl   - number of elements in output templates
;   lmin,lmax - lambda limits in output templates
;   
;
; OUTPUTS:
;   vmatrix   - template value [nl+1,nt]
;   lambda    - template wavelength [nl+1]
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
;   04-Jan-2002  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_load_templates, templatelist, vmatrix, lambda, lmin=lmin, lmax=lmax, nl=nl

; Need at least 3 parameters
if (N_params() LT 3) then begin
    klog, 'Syntax - k_load_templates, templatelist, vmatrix, lambda, [nl=, $'
    klog, '    lmin=, lmax=]'

    return
endif

; Set defaults
if (NOT keyword_set(nl)) then nl=500l
if (NOT keyword_set(lmin)) then lmin=1000.d
if (NOT keyword_set(lmax)) then lmax=12000.d
lambda=lmin+(lmax-lmin)*dindgen(nl+1l)/double(nl)
lint=0.5*(lambda[0l:nl-1l]+lambda[1l:nl])

; Read in vmatrix
vmatrix=dblarr(nl,n_elements(templatelist))
for i = 0l, n_elements(templatelist)-1l do begin
    spawn,'cat '+templatelist[i]+' | wc -l', nltmp
    vtmp=dblarr(2,nltmp[0])
    openr,unit,templatelist[i],/get_lun
    readf,unit,vtmp
    close,unit
    free_lun,unit
    indx=lindgen(nltmp[0]-1)
    for j = 0l, nl-1l do begin
        llo=where(vtmp[0l,indx] lt lint[j] and vtmp[0l,indx+1] ge lint[j])
        lhi=llo+1l
        sl=lint[j]-vtmp[0l,llo]
        vmatrix[j,i]=vtmp[1l,llo]+(vtmp[1l,lhi]-vtmp[1l,llo])*sl/ $
          (vtmp[0l,lhi]-vtmp[0l,llo])
    endfor
endfor 

end
