;+
; NAME:
;   read_draine
; PURPOSE:
;   Read one of Bruce Draine's dust spectrum files
; CALLING SEQUENCE:
;   model= read_draine(filename)
; INPUTS:
;   filename - name of file
; OUTPUTS:
;   model - structure containing descriptions of dust model:
;     .LAMBDA[N] - wavelengths (angstroms)
;     .FLUX[N] - ergs cm^{-2} s^{-1} A^{-1} at 10pc
;     .GRAIN_MODEL
;     .A01
;     .SIGMA_1
;     .B_C1
;     .A02
;     .SIGMA_2
;     .B_C2
;     .UMIN
;     .UMAX
;     .BETA
;     .AVGU
;     .RADFIELD
; COMMENTS:
;   I don't know the meaning of the model parameters.
; REVISION HISTORY:
;   06-May-2005  Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function read_draine, filename

model={grain_model:0L, $
       a01:0., $
       sigma_1:0., $
       b_c1:0., $
       a02:0., $
       sigma_2:0., $
       b_c2:0., $
       umin:0., $
       umax:0., $
       beta:0., $
       avgu:0., $
       radfield:''}

openr, unit, filename, /get_lun

line=''
readf, unit, line
words=strsplit(line,/extr)
model.grain_model=long(words[0])

readf, unit, line
words=strsplit(line,/extr)
model.a01=float(words[0])
model.sigma_1=float(words[1])
model.b_c1=float(words[2])

readf, unit, line
words=strsplit(line,/extr)
model.a02=float(words[0])
model.sigma_2=float(words[1])
model.b_c2=float(words[2])

readf, unit, line
words=strsplit(line,/extr)
model.umin=float(words[0])
model.umax=float(words[1])
model.beta=float(words[2])

readf, unit, line
words=strsplit(line,/extr)
model.avgu=float(words[0])

readf, unit, line
words=strsplit(line,/extr)
model.radfield=strtrim(words[0],2)

while(strtrim(words[0],2) ne 'lambda') do begin
    readf, unit, line
    words=strsplit(line,/extr)
endwhile
readf,unit,line

while(NOT eof(unit)) do begin
    readf,unit,line
    words=strsplit(line,/extr)
    lambda1= float(words[0])
    flux1= float(words[1])
    if(NOT keyword_set(lambda)) then begin
        lambda=lambda1
        flux=flux1
    endif else begin
        lambda=[lambda,lambda1]
        flux=[flux,flux1]
    endelse
endwhile

;; microns to angstroms
lambda=lambda*1.e+4

;; to per A
flux=flux/lambda

free_lun, unit

model=create_struct('lambda', lambda, $
                    'flux', flux, $
                    model)
return,model

end
;------------------------------------------------------------------------------

