; $Id$
;
; Copyright (c) 1992-1998, CreaSo Creative Software Systems GmbH,
;	and Research Systems, Inc.  All rights reserved.
;	Unauthorized reproduction prohibited.
;+
; NAME:
;    STR_SEP
;
; PURPOSE:
;    This routine cuts a string into pieces which are separated by the 
;    separator string.
; CATEGORY:
;    String processing.
; CALLING SEQUENCE:
;    arr = STR_SEP(str, separator)
;
; INPUTS:
;    str - The string to be separated.
;    separator - The separator.
;
; KEYWORDS:
;    ESC = escape character.  Only valid if separator is a single character.
;		Characters following the escape character are treated
;		literally and not interpreted as separators.
;		For example, if the separator is a comma,
;		and the escape character is a backslash, the character
;		sequence 'a\,b' is a single field containing the characters
;		'a,b'.
;    REMOVE_ALL = if set, remove all blanks from fields.
;    TRIM = if set, remove only leading and trailing blanks from fields.
;
; OUTPUT: 
;    An array of strings as function value.
;
; COMMON BLOCKS:
;    None
;
; SIDE EFFECTS:
;    No known side effects.
;
; RESTRICTIONS:
;    None.
;
; EXAMPLE:
;    array = STR_SEP ("ulib.usca.test", ".")
;
; MODIFICATION HISTORY:
;	July 1992, AH,	CreaSo		Created.
;	December, 1994, DMS, RSI	Added TRIM and REMOVE_ALL.
;	March, 2002, MRB, add to kcorrect to remove idlutils dependency
;-
function K_STR_SEP, str, separator, REMOVE_ALL = remove_all, TRIM = trim, ESC=esc


ON_ERROR, 2
if n_params() ne 2 then message,'Wrong number of arguments.'

spos = 0L
if n_elements(esc) gt 0 then begin		;Check for escape character?
  if strpos(str, esc) lt 0 then goto, no_esc	;None in string, use fast case
  besc = (byte(esc))[0]
  bsep = (byte(separator))[0]
  new = bytarr(strlen(str)+1)
  new[0] = byte(str)
  j = 0L
  for i=0L, n_elements(new)-2 do begin
    if new[i] eq besc then begin
	new[j] = new[i+1]
	i = i + 1
    endif else if new[i] eq bsep then new[j] = 1b $   ;Change seps to 1b char
    else new[j] = new[i]
    j = j + 1
    endfor
  new = string(new[0:j-1])
  w = where(byte(new) eq 1b, count)  ;where seps are...
  arr = strarr(count+1)
  for i=0L, count-1 do begin
	arr[i] = strmid(new, spos, w[i]-spos)
	spos = w[i] + 1
	endfor
  arr[count] = strmid(new, spos, strlen(str))  ;Last element
  goto, done
  endif			;esc

no_esc:
if strlen(separator) le 1 then begin	;Single character separator?
    w = where(byte(str) eq (byte(separator))[0], count)  ;where seps are...
    arr = strarr(count+1)
    for i=0, count-1 do begin
	arr[i] = strmid(str, spos, w[i]-spos)
	spos = w[i] + 1
	endfor
    arr[count] = strmid(str, spos, strlen(str))  ;Last element
endif else begin		;Multi character separator....
    n = 0L		   ; Determine number of seperators in string.
    repeat begin
	pos = strpos (str, separator, spos)
	spos = pos + strlen(separator)
	n = n+1
    endrep until pos eq -1

    arr = strarr(n)	   ; Create result array
    spos = 0L
    for i=0L, n-1 do begin   ; Separate substrings
      pos = strpos (str, separator, spos)
      if pos ge 0 then arr[i] = strmid (str, spos, pos-spos) $
      else arr[i] = strmid(str, spos, strlen(str))
      spos = pos+strlen(separator)
   endfor
endelse

done:
if keyword_set(trim) then arr = strtrim(arr,2) $
else if keyword_set(remove_all) then arr = strcompress(arr, /REMOVE_ALL)
return, arr
end
