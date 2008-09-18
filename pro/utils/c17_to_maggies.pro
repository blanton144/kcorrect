;+
; NAME:
;   c17_to_maggies
; PURPOSE:
;   convert c17 table3 input to Galactic-extinction corrected AB maggies 
; CALLING SEQUENCE:
;   c17_to_maggies, c17, maggies, ivar
; INPUTS:
;    c17 - structure with COMBO17 data:
;            UF_S
;            BF_S
;            VF_D
;            RF_S
;            IF_D
;            W420F_E 
;            W462F_E 
;            W485F_D 
;            W518F_E 
;            W571F_S 
;            W604F_E 
;            W646F_D 
;            W696F_E 
;            W753F_E 
;            W815F_S 
;            W856F_D 
;            W914F_E 
;         plus all of the above with E_ prefixed, plus:
;            APD_RMAG
; OUTPUTS:
;   maggies - [5, N] output in AB maggies in UBVRI
;   ivar - [5, N] inverse variance of maggies
; COMMENTS:
;   It ALWAYS applies a minimum error of 0.01 mag in all bands
; REVISION HISTORY:
;   30-June-2005 D. Schiminovich
;   22-July-2008 Guangtun Zhu, NYU, Initiate ivar with 1E-32
;   17-Sep-2008 Guangtun Zhu, NYU, Initiate ivar with 0.
;-
;------------------------------------------------------------------------------
pro c17_to_maggies, c17, maggies, ivar

; based on c17 paper we first convert to AB maggies

;w420F_E       ; E
;w462F_E       ; E
;w485F_D       ; D
;w518F_E       ; E
;w571F_S       ; D/E/S  ; pick S
;w604F_E       ; E
;w646F_D       ; D
;w696F_E       ; E
;w753F_E       ; E
;w815F_S       ; E/G/S ; pick S
;w856F_D       ; D
;w914F_E       ; D/E  ; pick E


;-0.19   1.571
;-0.18   1.412
;-0.06   1.207
;-0.06   1.125
;0.04    0.932
;0.10    0.832
;0.22    0.703
;0.27    0.621
;0.36    0.525
;0.45    0.442
;0.56    0.386
;0.50    0.380  

abconv=[0.77,-0.13,-0.02,0.19,0.49,-0.19,-0.18,-0.06,-0.06,0.04, $
        0.1,0.22,0.27,0.36,0.45,0.56,0.50]
flux=1.0E8*[0.737,1.371,1.055,0.725,0.412,1.571,1.412,1.207,1.125, $
            0.932,0.832,0.703,0.621,0.525,0.442,0.386,0.380]
tags=['UF_S','BF_S','VF_D','RF_S','IF_D','W420F_E','W462F_E','W485F_D', $
      'W518F_E','W571F_S','W604F_E','W646F_D','W696F_E','W753F_E', $
      'W815F_S','W856F_D','W914F_E']

apd_rmag=dblarr(n_elements(c17))
apd_rmag=c17.apd_rmag
ii=where(apd_rmag ne apd_rmag, nii)
if(nii gt 0) then apd_rmag[ii]=0.

n=n_elements(c17)
maggies=dblarr(17,n)
ivar=dblarr(17,n)
for i=0L,16L do begin
    if tag_exist(c17,tags[i],index=ptr) then begin
        igood=where(c17.(ptr) eq c17.(ptr), ngood)
        if(ngood gt 0) then $
          maggies[i,igood]= $
          (c17[igood].(ptr)/flux[i])*10.0^((-.4)*(abconv[i]+apd_rmag[igood]))
    endif
    if tag_exist(c17,'E_'+tags[i],index=ptr) then begin
        igood=where(c17.(ptr) eq c17.(ptr) AND $
                    maggies[i,*] ne 0., ngood)
        if(ngood gt 0) then $
          ivar[i,igood]= $
          1.0/((c17[igood].(ptr)/flux[i])* $
               10.0^((-.4)*(abconv[i]+apd_rmag[igood])))^2
    endif
endfor

k_minerror, maggies, ivar, replicate(0.1, 17)

end
