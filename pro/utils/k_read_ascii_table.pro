;+
; NAME:
;   k_read_ascii_table
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
; CALLING SEQUENCE:
;   k_read_ascii_table,table,filename
; INPUTS:
;   filename  - file containing table
; OPTIONAL INPUTS:
; OUTPUTS:
;   table   - table to read
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
;   Only reads floating point tables.
; EXAMPLES:
;    Here is an example of reading a [500,10] array from a file:
; 
;    IDL> $wc bmatrix.default.dat   
;       5001    5003  130037 bmatrix.default.dat
;    IDL> $head bmatrix.default.dat
;               2          10         500
;       4.1921895303265806e-04
;       4.1921245277927205e-04
;       4.1920432703966440e-04
;       4.1919374271866136e-04
;       4.1918060198276090e-04
;       4.1916390305939648e-04
;       4.1914261546946777e-04
;       4.1911632959028446e-04
;       4.1908366349717002e-04
;    IDL> k_read_ascii_table,bmatrix,'bmatrix.default.dat'
;    % Compiled module: K_READ_ASCII_TABLE.
;    IDL> help,bmatrix
;    BMATRIX         DOUBLE    = Array[500, 10]
;
; REVISION HISTORY:
;   17-Jan-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_read_ascii_table,table,filename 

; Need at least 2 parameters
if (N_params() LT 2) then begin
    print, 'Syntax - k_read_ascii_table, table, filename'
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
table=make_array(/float,dimension=sizes)
readf,unit,table
close,unit
free_lun,unit

end
;------------------------------------------------------------------------------

