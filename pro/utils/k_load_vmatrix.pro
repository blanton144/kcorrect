;+
; NAME:
;   k_load_vmatrix
; PURPOSE:
;   Load the template information 
; CALLING SEQUENCE:
;   k_load_vmatrix, vmatrix, lambda [, vfile=, vpath=, lfile= ]
; INPUTS:
; OPTIONAL INPUTS:
;   vfile      - ascii format file with vmatrix in it [default
;                vmatrix.default.dat]
;   lfile - ascii format file with wavelengths in it (in Angstroms)
;   vpath      - path to use for vfile if it does not exist 
;                (if this has more than one element they are checked in order)
;                [default $KCORRECT_DIR/data/templates]
; OUTPUTS:
;   vmatrix   - [nl,nv] template
;   lambda    - [nl+1] pixel edges of template
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   04-Jun-2003  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function k_load_vmatrix_path,vfile,vpath=vpath
curr_path='.'
ipath=-1
if(NOT file_test(vfile)) then begin
    for i=0L, n_elements(vpath)-1L do begin
        if(ipath eq -1 AND $
           file_test(vpath[i]+'/'+vfile)) then $
          ipath=i
    endfor 
    if(ipath eq -1) then $
      curr_path=getenv('KCORRECT_DIR')+'/data/templates'
endif
return, curr_path
end
;
pro k_load_vmatrix, vmatrix, lambda, lfile=lfile, vfile=vfile, $
                    vpath=vpath

; Need at least 4 parameters
if (N_params() LT 2) then begin
    print, 'Syntax - k_load_vmatrix, vmatrix, lambda [, vfile=, vpath=, lfile= ]'
    return
endif

if(NOT keyword_set(vfile)) then $
  vfile='vmatrix.default.dat'
if(NOT keyword_set(lfile)) then $
  lfile='lambda.default.dat'

curr_path=k_load_vmatrix_path(vfile,vpath=vpath)
k_read_ascii_table,vmatrix,curr_path+'/'+vfile
curr_path=k_load_vmatrix_path(lfile,vpath=vpath)
k_read_ascii_table,lambda,curr_path+'/'+lfile

end
