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
filter_ncolumns=lonarr(n_elements(filterlist))
filter_default=lonarr(n_elements(filterlist))
for i = 0l, n_elements(filterlist)-1l do begin
    readcol,filterlist[i],nlines,ncolumns,defaultcol, $
      format='I,I,I',comment='#',/silent
    filter_default[i]=defaultcol[0]
    filter_n[i]=nlines[0]
    filter_ncolumns[i]=ncolumns[0]
endfor 
maxn=max(filter_n)

; Read in each filter
filter_lambda=dblarr(maxn,n_elements(filterlist))
filter_pass=dblarr(maxn,n_elements(filterlist))
for i = 0l, n_elements(filterlist)-1l do begin
    format='D'
    for j=1, filter_default[i]-2 do format=format+',X'
    format=format+',D'
    readcol,filterlist[i],tmp_lambda,tmp_pass, $
      format=format,comment='#',/silent
    if(n_elements(tmp_lambda) eq filter_n[i]+1) then begin 
        filter_lambda[0l:filter_n[i]-1l,i]=tmp_lambda[1L:filter_n[i]]
        filter_pass[0l:filter_n[i]-1l,i]=tmp_pass[1L:filter_n[i]]
    endif else begin
        filter_lambda[0l:filter_n[i]-1l,i]=tmp_lambda
        filter_pass[0l:filter_n[i]-1l,i]=tmp_pass
    endelse 
endfor 

end
