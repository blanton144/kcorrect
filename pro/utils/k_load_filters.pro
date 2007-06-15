;+
; NAME:
;   k_load_filters
; PURPOSE:
;   Load the filter information from a list of files
; CALLING SEQUENCE:
;   k_load_filters, filterlist, filter_nlambda, filter_lambda, filter_pass [, $
;       filterpath=filterpath
; INPUTS:
;   filterlist      - list of files with filter information
; OPTIONAL INPUTS:
;   filterpath      - path to use for filters if the filterlist files 
;                     do not exist (if this has more than one element
;                     they are checked in order)
; OUTPUTS:
;   filter_nlambda  - number of elements in each filter specification
;   filter_lambda   - wavelength for each element
;   filter_pass     - value for each element
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
;   Filters should be in the Yanny parameter format with columns:
;          lambda
;          pass
;   This code caches the filters it reads in. 
; PROCEDURES CALLED:
;   yanny_readone (idlutils)
; REVISION HISTORY:
;   04-Jan-2002  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function k_load_filters_path,filtername,filterpath=filterpath
curr_path='.'
ipath=-1
if(NOT file_test(filtername)) then begin
    for i=0L, n_elements(filterpath)-1L do begin
        if(ipath eq -1 AND $
           file_test(filterpath[i]+'/'+filtername)) then $
          ipath=i
    endfor 
    if(ipath eq -1) then $
      curr_path=getenv('KCORRECT_DIR')+'/data/filters' $
    else $
      curr_path=filterpath[ipath]
endif
return, curr_path
end
;
pro k_load_filters, filterlist, filter_nlambda, filter_lambda, filter_pass, $
                    filterpath=filterpath

common k_load_filters_com, filters

; Need at least 4 parameters
if (N_params() LT 4) then begin
    print, 'Syntax - k_load_filters, filterlist, filter_nlambda, $'

    print, '           filter_lambda, filter_pass [, filterpath= ]'
    return
endif

; Read in sizes and pick maximum
filter_nlambda=lonarr(n_elements(filterlist))
for i = 0l, n_elements(filterlist)-1l do begin
    for j=0L, n_elements(filters)-1L do begin
        if(filters[j].file eq filterlist[i]) then $
          filter_str=*filters[j].str
    endfor
    if(n_tags(filter_str) eq 0) then begin
        curr_path=k_load_filters_path(filterlist[i],filterpath=filterpath)
        if(NOT file_test(curr_path+'/'+filterlist[i])) then $
          message,'filter '+curr_path+'/'+filterlist[i]+' does not exist'
        filter_str=yanny_readone(curr_path+'/'+filterlist[i], hdr=hdr, $ 
                                 /anonymous) 
                                 
        filters_new={file:filterlist[i], $
                     str:ptr_new(filter_str)}
        if(n_elements(filters) gt 0) then $
          filters=[filters,filters_new] $
        else $
          filters=filters_new
    endif
    filter_nlambda[i]=n_elements(filter_str)
    filter_str=0
endfor 
maxn=max(filter_nlambda)

; Read in each filter
filter_lambda=fltarr(maxn,n_elements(filterlist))
filter_pass=fltarr(maxn,n_elements(filterlist))
for i = 0l, n_elements(filterlist)-1l do begin
    for j=0L, n_elements(filters)-1L do begin
        if(filters[j].file eq filterlist[i]) then $
          filter_str=*filters[j].str
    endfor
    filter_lambda[0L:filter_nlambda[i]-1L,i]=filter_str.lambda
    filter_pass[0L:filter_nlambda[i]-1L,i]=filter_str.pass
endfor 

end
