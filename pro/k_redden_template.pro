;+
; NAME:
;   k_redden_template
;
; PURPOSE:
;   Redden a template
;
; CALLING SEQUENCE:
;
; INPUTS:
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
;   k_load_ascii_table
;   k_write_ascii_table
;
; REVISION HISTORY:
;   02-Apr-2002  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_redden_template,intemplatefiles,outtemplatefiles,redden,lmin=lmin,lmax=lmax,nl=nl

k_load_templates,intemplatefiles,vmatrix,lambda,lmin=lmin,lmax=lmax,nl=nl
nl=n_elements(vmatrix)/n_elements(intemplatefiles)
lval=0.5*(lambda[0:nl-1]+lambda[1:nl])
vmatrix=reform(vmatrix,nl,n_elements(intemplatefiles))
outvmatrix=dblarr(nl,n_elements(intemplatefiles))
for i=0L, n_elements(intemplatefiles)-1L do begin
    outvmatrix[*,i]=10^(-0.4*ext_ccm(lval,redden))* $
      vmatrix[*,i]
endfor 

outarr=dblarr(2,nl)
for i=0L, n_elements(outtemplatefiles)-1L do begin
    outarr[0,*]=lval
    outarr[1,*]=outvmatrix[*,i]
    openw,unit,outtemplatefiles[i],/get_lun
    printf,unit,outarr
    close,unit
    free_lun,unit
endfor

end
