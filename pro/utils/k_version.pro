;+
; NAME:
;   k_version
; PURPOSE:
;   Return the version name for the product kcorrect
; CALLING SEQUENCE:
;   vers = k_version()
; OUTPUTS:
;   vers       - Version name for the product kcorrect
; COMMENTS:
;   If this version is not tagged by CVS, then we return 'NOCVS:TOPLEVEL'
;   where TOPLEVEL is the last directory in the environment variable
;   $KCORRECT_DIR.  For example, if you are using a version of the code
;   in the directory '/u/schlegel/kcorrect/v0_0', then this returns
;   'NOCVS:v0_0'.
; REVISION HISTORY:
;   01-Dec-1999  Written by D. Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function k_version

   ; The following expression in dollar signs is expanded by CVS
   ; and replaced by the tag name for this version.
   name = '$Name$'

   words = str_sep(strcompress(name), ' ')

   if (words[0] EQ '$Name:' AND N_elements(words) EQ 3) then begin
      vers = words[1]
   endif else begin
      dirname = getenv('KCORRECT_DIR')
      if (dirname NE '') then begin
         words = str_sep(dirname,'/')
         nword = N_elements(words)
         vers = 'NOCVS:' + words[nword-1]
      endif else begin
         vers = 'NOCVS:Unknown'
      endelse
   endelse

   return, vers
end
;------------------------------------------------------------------------------
