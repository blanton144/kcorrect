;+
; NAME:
;   k_load_ascii_table
;
; PURPOSE:
;   Read an ascii file in my standard format, which is:
; 
;  <ndim> <size_{0}> ... <size_{ndim-1}>
;  <entry_0>
;  <entry_1>
;  ...
;  <entry_n>
;  ...
;  <entry_{size_0*size_1*..*size_{ndim-1}-1>
; 
;  where the table element [k,j,i] (for ndim==3) would be the 
;  entry element n, where n=i*size_2*size1+j*size2+k
;
; CALLING SEQUENCE:
;   k_load_ascii_table,table,filename
;
; INPUTS:
;   filename  - file containing table
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;   table   - table to read
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
;   17-Jan-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_load_ascii_table,table,filename 

; Need at least 2 parameters
if (N_params() LT 2) then begin
    klog, 'Syntax - k_load_ascii_table, table, filename'
    return
endif

openr,unit,filename,/get_lun
ndim=0L
readf,unit,ndim
close,unit

openr,unit,filename
tmpsizes=lonarr(ndim)
sizes=lonarr(ndim)
readf,unit,dim,tmpsizes
for i=0, ndim-1 do $
   sizes[i]=tmpsizes[ndim-i-1]
table=make_array(/double,dimension=sizes)
readf,unit,table
close,unit
free_lun,unit

end
;------------------------------------------------------------------------------

