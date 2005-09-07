;+
; NAME:
;   read_mappings
; PURPOSE:
;   Read a MAPPINGS III file
; CALLING SEQUENCE:
;   model= read_mappings(filename [, /vac, /nolya] )
; INPUTS:
;   filename - name of file
;   /vac - convert to vacuum wavelengths
;   /nolya - omit lyman alpha
; OUTPUTS:
;   model - structure containing descriptions of model:
;     .RUN - name of run
;     .LAMBDA[N] - wavelengths (angstroms)
;     .FLUX[N] - fluxs (erg s^{-1} ?)
;     .NAMES[N] - names of lines
; COMMENTS:
;   I don't know the meaning of the model parameters.
; REVISION HISTORY:
;   06-May-2005  Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function read_mappings, filename, vac=vac, nolya=nolya

model={run:''}

openr, unit, filename, /get_lun

line=''
readf, unit, line
readf, unit, line
readf, unit, line
readf, unit, line
words=strsplit(line,/extr)
model.run=strtrim(strjoin(words[2:n_elements(words)-1L], ' '),2)

while(strtrim(words[0],2) ne 'Lambda(A)') do begin
    readf, unit, line
    words=strsplit(line,/extr)
endwhile
readf,unit,line

while(NOT eof(unit)) do begin
    readf,unit,line
    words=strsplit(line,/extr)
    lambda1= float(words[0])
    flux1= float(words[2])
    names1= strjoin(words[3:4])
    if(lambda1 lt 1210. OR lambda1 gt 1220. OR keyword_set(nolya) eq 0) then $
      begin
        if(NOT keyword_set(lambda)) then begin
            lambda=lambda1
            flux=flux1
            names=names1
        endif else begin
            lambda=[lambda,lambda1]
            flux=[flux,flux1]
            names=[names,names1]
        endelse
    endif
endwhile

free_lun, unit

if(keyword_set(vac)) then $
  airtovac, lambda

model=create_struct('lambda', lambda, $
                    'flux', flux, $
                    'names', names, $
                    model)
return,model

end
;------------------------------------------------------------------------------

