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
pro k_redden_template,intemplatefiles,outtemplatefiles,redden

k_load_templates,intemplatefiles,vmatrix,lambda
nl=n_elements(vmatrix)/n_elements(intemplatefiles)
vmatrix=reform(vmatrix,nl,n_elements(intemplatefiles))
outvmatrix=vmatrix-vmatrix
for i=0L, n_elements(intemplatefiles)-1L do begin
    outvmatrix[*,i]=10^(-0.4*ext_ccm(lambda[0:nl-1],redden))* $
      vmatrix[*,i]
endfor 

for i=0L, n_elements(outtemplatefiles)-1L do begin
    
endfor

end
