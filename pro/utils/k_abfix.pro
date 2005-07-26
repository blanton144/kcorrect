;+
; NAME:
;   k_abfix
; PURPOSE:
;   Take SDSS pipeline maggies and return AB maggies
; CALLING SEQUENCE:
;   k_abfix, maggies, maggies_ivar [ aboff=]
; INPUTS/OUTPUTS:
;   maggies - input maggies (changed on output)
;   maggies_ivar - input inverse variances (changed on output)
; OPTIONAL KEYWORDS:
;   /ab02 - use 2002 version of AB corrections (overrides aboff keyword)
; OPTIONAL INPUTS:
;   aboff - [5] AB offsets (defaults to list below)
; COMMENTS:
;   Uses the AB conversions posted by D. Eisenstein (sdss-calib/???)
;     u(AB,2.5m) = u(database, 2.5m) - 0.036
;     g(AB,2.5m) = g(database, 2.5m) + 0.012
;     r(AB,2.5m) = r(database, 2.5m) + 0.010
;     i(AB,2.5m) = i(database, 2.5m) + 0.028
;     z(AB,2.5m) = z(database, 2.5m) + 0.040
;   You can specify your own with the "aboff" input.
;
;   Note that older ("2002") versions (v3_2 and previous) used a different set
;   of corrections: 
;     u(AB,2.5m) = u(database, 2.5m) - 0.042   (NOT WHAT WE DO HERE)
;     g(AB,2.5m) = g(database, 2.5m) + 0.036   (NOT WHAT WE DO HERE)
;     r(AB,2.5m) = r(database, 2.5m) + 0.015   (NOT WHAT WE DO HERE)
;     i(AB,2.5m) = i(database, 2.5m) + 0.013   (NOT WHAT WE DO HERE)
;     z(AB,2.5m) = z(database, 2.5m) - 0.002   (NOT WHAT WE DO HERE)
;   You can use these if you set the /ab02 flag
; REVISION HISTORY:
;   05-Apr-2005  Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_abfix, maggies, maggies_ivar, aboff=aboff, ab02=ab02

if(keyword_set(ab02)) then aboff=[-0.042, 0.036, 0.015, 0.013, -0.002]
if(NOT keyword_set(aboff)) then aboff=[-0.036, 0.012, 0.010, 0.028, 0.040]

for i=0L, n_elements(aboff)-1L do begin
    maggies[i,*]=maggies[i,*]*10.^(-0.4*aboff[i])
    maggies_ivar[i,*]=maggies_ivar[i,*]*10.^(0.8*aboff[i])
endfor

end
;------------------------------------------------------------------------------
