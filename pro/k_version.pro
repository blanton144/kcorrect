;+
; NAME:
;   k_version
;
; PURPOSE:
;   Return version of K-correction code
;
; CALLING SEQUENCE:
;   k_version
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
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   05-Mar-2002  Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_version

spawn,'cat '+getenv('KCORRECT_DIR')+'/VERSION',version
splog,version[0]

end
