;+
; NAME:
;   k_load_filters
;
; PURPOSE:
;   Load the filter information from a list of files
;
; CALLING SEQUENCE:
;   k_load_filters, filterlist, filter_n, filter_lambda, filter_pass
;
; INPUTS:
;   filterlist      - list of files with filter information
;   filter_n        - number of elements in each filter specification
;   filter_lambda   - wavelength for each element
;   filter_pass     - value for each element
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
;   04-Jan-2002  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_load_filters, filterlist, filter_n, filter_lambda, filter_pass

; Need at least 4 parameters
if (N_params() LT 4) then begin
    klog, 'Syntax - k_load_filters, filterlist, filter_n, filter_lambda, filter_pass'
    return
endif

; Read in sizes and pick maximum
filter_n=lonarr(n_elements(filterlist))
for i = 0l, n_elements(filterlist)-1l do begin
    tmp=1l
    openr,unit,filterlist[i],/get_lun
    readf,unit,tmp
    filter_n[i]=tmp
    close,unit
    free_lun,unit
endfor 
maxn=max(filter_n)

; Read in each filter
filter_lambda=dblarr(maxn,n_elements(filterlist))
filter_pass=dblarr(maxn,n_elements(filterlist))
for i = 0l, n_elements(filterlist)-1l do begin
    openr,unit,filterlist[i],/get_lun
    tmp=1l
    readf,unit,tmp
    filter_n[i]=tmp
    tmp=dblarr(2,filter_n[i]) 
    readf,unit,tmp
    filter_lambda[0l:filter_n[i]-1l,i]=tmp[0,*]
    filter_pass[0l:filter_n[i]-1l,i]=tmp[1,*]
    close,unit
    free_lun,unit
endfor 
end
