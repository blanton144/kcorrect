;+
;NAME:
;  read_peg.pro
;PURPOSE:
;  reads the output document of the pegase program
;CALLING SEQUENCE:
;  read_peg, peg_file
;INTPUTS:
;  peg_file      - output file of the pegase program
;OUTPUTS:
;  peg           - an array of structures (one for each timestep)
;REVISION HISTORY:
;  2002-Aug-28 written by Quintero
;_
;------------------------------------------------------------
pro k_read_peg, peg_file, peg = peg

if n_params() ne 1 then begin
    print,'Syntax - read_peg, peg_file [, peg=]'
    return
endif

openr, unit, peg_file, /get_lun
print, unit

a = ' '
readf, unit, a
WHILE a NE '************************************************************' $
  DO BEGIN 
    readf, unit, a
;    print, a
ENDWHILE

point = lonarr(3)
readf, unit, point 

lam = lonarr(point[1])
readf, unit, lam

lines = lonarr(point[2])
readf, unit, lines

arr1 = fltarr(10)
arr2 = fltarr(9)
cont_lum_arr = fltarr(point[1])
line_lum_arr = fltarr(point[2])

peg = peg_struct(point[0], point[1], point[2])

FOR i=0, point[0]-1 DO BEGIN
    readf, unit, arr1 
    readf, unit, arr2
    readf, unit, cont_lum_arr
    readf, unit, line_lum_arr
    
    peg.Ntimesteps = point[0]
    peg.Ncont = point[1]
    peg.Nlines = point[2]
    peg[i].Arr1 = arr1
    peg[i].Arr2 = arr2
    peg.Cont = lam
    peg.Lines = lines
    peg[i].Contlum = cont_lum_arr
    peg[i].Linelum = line_lum_arr
    
;   if NOT keyword_set(arr1vec) then arr1vec = arr1 $
;    else arr1vec = [ [arr1vec], [arr1] ]
;   if NOT keyword_set(arr2vec) then arr2vec = arr2 $
;    else arr2vec = [ [arr2vec], [arr2] ]
;   if NOT keyword_set(cont_lum_vec) then cont_lum_vec = cont_lum_arr $
;    else cont_lum_vec = [ [cont_lum_vec], [cont_lum_arr] ]
;   if NOT keyword_set(line_lum_vec) then line_lum_vec = line_lum_arr $
;    else line_lum_vec = [ [line_lum_vec], [line_lum_arr] ]
    
;   index = where(lam GT 3400 AND lam LT 7000)
    
;   djs_plot, lam[index], cont_lum_arr[index]
;   djs_oplot, lines, line_cont_arr, psym = 1, color = 'red'
ENDFOR

end
