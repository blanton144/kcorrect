;+
; NAME:
;       K_IM_READ_BC03()
;
; PURPOSE:
;       Read the Bruzual & Charlot (2003) binary format population
;       synthesis models into a convenient data structure.
;
; CALLING SEQUENCE:
;       bc03 = k_im_read_bc03(isedfile=,isedpath,metallicity=,$
;          age=,minwave=,maxwave=,bc03_extras=,/salpeter,/lr,$
;          /silent)
;
; INPUTS:
;       None required.  By default this routine reads the high
;       resolution (hr), solar metallicity (m62), Chabrier SSP
;       models. 
;
; OPTIONAL INPUTS:
;       isedfile - read this binary SED rather than the default SSP
;                  models 
;       isedpath - full data path (with trailing slash) to ISEDFILE
;       metallicity - SSP metallicity
;          0 - Z=0.0001 (m22)
;          1 - Z=0.0004 (m32)
;          2 - Z=0.004  (m42)
;          3 - Z=0.008  (m52)
;          4 - Z=0.02   (m62) (Solar, default)
;          5 - Z=0.05   (m72)
;       age - return the SSP(s) corresponding to this scalar or vector
;             age(s) [Gyr]
;       minwave - crop the SSP spectra to this minimum wavelength (Angstrom)
;       maxwave - crop the SSP spectra to this maximum wavelength (Angstrom)
;       isolib - Isochrone library to use (default 'Padova1994')
;
; KEYWORD PARAMETERS:
;       salpeter - read the Salpeter IMF models (default is to read
;                  the Chabrier models)
;       lr       - read the low resolution models (default is to read
;                  the high resolution models)
;       silent   - do not print any messages to STDOUT
;       vac      - translate wavelengths to vacuum
;
; OUTPUTS:
;       bc03 - data structure with the following fields:
;          age  - vector of SSP ages [NAGE] (yr)
;          wave - wavelength array [NPIX] (Angstrom)
;          flux - SSP models [NPIX,NAGE] (L_sun/M_sun/A)
;
; OPTIONAL OUTPUTS:
;       bc03_extras - data structure containing the extra parameters
;                     associated with each SSP (see the BC03
;                     documentation) 
;
; COMMENTS:
;       N.B.  The environment variable ${bc03_dir} must be defined in
;       your '.cshrc' indicating the *root* directory of the BC03
;       models.  For example, in my '.cshrc' I have
;
;              setenv bc03_dir ${HOME}/synthesis/bc03
;
;       Note that this routine works with the Padova (1994) binary
;       isochrones by default. Use ISOLIB input to use others. Using
;       the ISEDFILE and ISEDPATH optional inputs you can read in an
;       arbitrary BC03 SED in binary format.
;
; BUGS:
;       It does not appear that this routine works with BC03 models
;       *other* than the instantaneous burst models.  For example, if
;       you use 'csp_galexev.f' to generate SFH-convolved models then
;       those outputted models cannot be read with IM_READ_BC03().  I
;       think the reason is that the time steps and spacing are
;       modified, which changes the format of the binary file.
;
;       If you solve this please let me know!
;
; INTERNAL SUPPORT ROUTINES:
;       READ_EXTRAS(), GET_ELEMENT
;
; PROCEDURES USED:
;       READFAST, STRUCT_TRIMTAGS(), STRUCT_ADDTAGS(), MATCH
;
; EXAMPLES:
;       [1] Read the high-resolution, Salpeter IMF, Solar metallicity
;           (Z=0.02) SSP models:
;
;          IDL> bc03 = k_im_read_bc03()
;          IDL> help, bc03, /str
;
;       [2] Read the low-resolution, Salpeter IMF, LMC metallicity
;           (Z=0.004) SSP models:
;
;          IDL> bc03 = k_im_read_bc03(metallicity=2,/lr,/salpeter)
;
;       [3] Read the high-resolution, Chabrier IMF, twice-solar
;           metallicity (Z=0.05) SSP models and plot the 10 Gyr model: 
;
;          IDL> bc03 = k_im_read_bc03(metallicity=5)
;          IDL> indx = where(bc03.age eq 1E10)
;          IDL> plot, bc03.wave, bc03.flux[*,indx], /xlog, /ylog
;
;       [4] Read the high-resolution, Chabrier IMF, solar metallicity
;           (Z=0.02) SSP models and extras and plot D4000 versus
;           H-delta_A:  
;
;          IDL> bc03 = k_im_read_bc03(bc03_extras=bce)
;          IDL> plot, bce.d_4000_, bce.h_delta_a, xsty=3, ysty=3
;
;       [5] Read a non-SSP SED in my temp subdirectory:
;
;          IDL> bc03 = k_im_read_bc03(isedfile='mysed.ised',isedpath='temp/')
;
;       [6] Retrieve the 12 Gyr Salpeter IMF, high-resolution SSP
;           between 3000 and 10000 Angstroms.
;
;          IDL> bc03 = k_im_read_bc03(/salpeter,age=12.0,minwave=3000,maxwave=1E4)
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 October 29, U of A, based in some part on
;          IDL code originally written by C. Papovich
;       jm03nov12uofa - added ISEDFILE and ISEDPATH optional inputs
;                       and various bug fixes
;       jm04mar01uofa - added AGE, MINWAVE, and MAXWAVE optional
;                       inputs; changed wave, flux, and age to double
;                       precision 
;-

function read_extras, extrafile, extrapath=extrapath, silent=silent

if file_test(extrapath+extrafile) eq 0L then begin
    print, 'Extras file '+extrapath+extrafile+' not found.'
    return, -1L
endif

if not keyword_set(silent) then $
  print, 'Reading extras file '+extrapath+extrafile+'.'

; read the first 50 lines of the file to figure out the column names
; and also to determine on which row the data begin 

temp = strarr(50L)
tmpstr = ''
openr, lun1, extrapath+extrafile, /get_lun
for i = 0L, 49L do begin
    readf, lun1, tmpstr
    temp[i] = tmpstr
endfor
free_lun, lun1

headlines = where(strmatch(temp,'*#*') eq 1B,nheadlines)
if (nheadlines eq 0L) then $
  message, 'There is a problem with '+extrapath+extrafile+'.'

agerow = where(strmatch(temp,'*log-age*',/fold) eq 1B,nagerow)
if (nagerow ne 1L) then $
  message, 'There is a problem with '+extrapath+extrafile+'.'

; read the column names and convert them to valid structure tag names

cols = strupcase(strcompress(strsplit(strmid(temp[agerow],1),' ',/extract), $
                             /remove))
cols[0] = 'LOGAGE'
ncols = n_elements(cols)
for j = 0L, ncols-1L do cols[j] = idl_validname(cols[j],/convert_all)

; read the data

readfast, extrapath+extrafile, data, header, skipline=nheadlines, $
  /double, nlines=nage, ncols=ncheck
if (ncheck ne ncols) then $
  message, 'There is a problem reading '+extrapath+extrafile+'.'

; initialize the output data structure and fill it

extras = create_struct(cols[0],0.0D)
for k = 1L, ncols-1L do extras = create_struct(extras,cols[k],0.0D)
extras = replicate(extras,nage)

for iage = 0L, ncols-1L do extras.(iage) = reform(data[iage,*])

return, extras
end

pro get_element, x, value, position
; jm01jan28uofa
; a generalization of GETELEMENT_VECTOR, this routine will also accept
; an array of positions, and return an array of indices

position = long(value-value)
for i = 0L, n_elements(value)-1L do begin
    array_value = min((abs(x-value[i])),temp)
    position[i] = temp
endfor

return
end

function k_im_read_bc03, isedfile=isedfile, isedpath=isedpath, $
                         metallicity=metallicity, age=age, minwave=minwave, $
                         maxwave=maxwave, salpeter=salpeter, lr=lr, $
                         silent=silent, isolib=isolib, vac=vac, $
                         bc03_extras=bc03_extras

if n_elements(metallicity) eq 0L then metallicity = 4L

if(n_elements(isolib) eq 0L) then $
  isolib='Padova1994'
if n_elements(isedpath) eq 0L then begin
    isedpath = getenv('bc03_dir')+'/models/'+isolib+'/'
    if keyword_set(salpeter) then $
      isedpath = isedpath+'salpeter/' else $
      isedpath = isedpath+'chabrier/'
endif

mprefix=''
if(isolib eq 'Padova2000') then $
  mprefix='1'

case long(metallicity) of
    0L: begin
        Zstr = 'm'+mprefix+'22'
        Zinfo = 'Z=0.0001'
    end
    1L: begin
        Zstr = 'm'+mprefix+'32'
        Zinfo = 'Z=0.0004'
    end
    2L: begin
        Zstr = 'm'+mprefix+'42'
        Zinfo = 'Z=0.004'
    end
    3L: begin
        Zstr = 'm'+mprefix+'52'
        Zinfo = 'Z=0.008'
    end
    4L: begin
        Zstr = 'm'+mprefix+'62'
        Zinfo = 'Z=0.02'
    end
    5L: begin
        Zstr = 'm'+mprefix+'72'
        Zinfo = 'Z=0.05'
    end
    else: begin
        Zstr = 'm'+mprefix+'62'
        Zinfo = 'Z=0.02'
    end
endcase

if keyword_set(salpeter) then begin
    imfstr = 'salp'
    imfinfo = 'Salpeter IMF' 
endif else begin
    imfstr = 'chab'
    imfinfo = 'Chabrier IMF'
endelse

if keyword_set(lr) then begin   ; low resolution

    npix = 1221L        
    resstr = 'lr'
    resinfo = 'Low Resolution'

    if keyword_set(salpeter) then begin
        ifs = 306L           ; starting index of the wavelength vector
        ifs2 = 2807L            ; starting index of the first spectrum
    endif else begin
        ifs = 300L
        ifs2 = 2801L
    endelse

endif else begin                ; high resolution

    npix = 6900L 
    resstr = 'hr'
    resinfo = 'High Resolution'

    if keyword_set(salpeter) then begin
        ifs = 306L   
        ifs2 = 7209L 
    endif else begin
        ifs = 300L
        ifs2 = 7203L
    endelse
    
endelse

; initialize some indices

offa = 2L           ; offset from beginning of file to first age index
offs = 56L                      ; space between spectra
nage = 221L                     ; number of age bins (and models)

; initialize the output data structure

bc03 = {age: dblarr(nage), wave: dblarr(npix), flux: dblarr(npix,nage)}

; read the binary file

if n_elements(isedfile) eq 0L then $
  isedfile = 'bc2003_'+resstr+'_'+Zstr+'_'+imfstr+'_ssp.ised' else $
  Zinfo = 'Z=Unknown'

if file_test(isedpath+isedfile) then begin
    if not keyword_set(silent) then begin
        print, 'Reading SSP file '+isedpath+isedfile+':'
        print, imfinfo+', '+resinfo+', '+Zinfo
    endif
    tempbin = read_binary(isedpath+isedfile,data_type=4, endian='little')
endif else begin
    print, 'SSP file '+isedpath+isedfile+' not found.'
    return, -1L
endelse

bc03.age = tempbin[offa:offa+nage-1L] ; age vector
bc03.wave = tempbin[ifs:ifs+npix-1L] ; wavelength vector

for i = 0L, nage-1L do begin

    i1 = ifs2 + i*(npix+offs)
    i2 = ifs2 + i*(npix+offs) + npix-1L
;      print, i1, i2, (i2-i1)+1

    bc03.flux[*,i] = tempbin[i1:i2]

endfor

ninputage = n_elements(age)
if (ninputage ne 0L) then begin

    get_element, bc03.age/1E9, age, indx
    if not keyword_set(silent) then begin
        splog, 'Retrieving the following SSP(s):'
        for j = 0L, ninputage-1L do print, '   '+ $
          string(bc03.age[indx[j]]/1E9,format='(F12.5)')+' Gyr'
    endif

    bc03 = {age: bc03.age[indx], wave: bc03.wave, flux: bc03.flux[*,indx]}
    
endif

nminwave = n_elements(minwave)
nmaxwave = n_elements(maxwave)

if (nminwave ne 0L) or (nmaxwave ne 0L) then begin

    if (nminwave eq 0L) then minwave = min(bc03.wave)
    if (nmaxwave eq 0L) then maxwave = max(bc03.wave)

    get_element, bc03.wave, [minwave,maxwave], ww
    if not keyword_set(silent) then begin
        splog, 'Reducing wavelength range to ['+ $
          string(minwave,format='(I0)')+','+$
          string(maxwave,format='(I0)')+'] Angstrom.'
    endif

    bc03 = {age: bc03.age, $
            wave: bc03.wave[(ww[0]-1L)>0L:(ww[1]+1L)<(npix-1L)], $
            flux: bc03.flux[(ww[0]-1L)>0L:(ww[1]+1L)<(npix-1L),*]}
    
endif

if arg_present(bc03_extras) then begin

; read the extra parameters for this SSP

    base = 'bc2003_'+resstr+'_'+Zstr+'_'+imfstr+'_ssp.'
    colorfiles = base+string(lindgen(5)+1,format='(I0)')+'color'
    ABmagfile = base+'1ABmag'
    
    if keyword_set(lr) then begin

        extrafiles = [colorfiles,ABmagfile]

    endif else begin

        indxfiles6 = base+'6lsindx_'+'sed' ; ['ffn','sed','sed_lick_system']
        indxfiles7 = base+'7lsindx_'+'sed' ; ['ffn','sed','sed_lick_system']

        extrafiles = [colorfiles,indxfiles6,indxfiles7,ABmagfile]

    endelse
    
    nextra = n_elements(extrafiles)

    for k = 0L, nextra-1L do begin

        extras1 = read_extras(extrafiles[k],extrapath=isedpath,silent=silent)
        if (size(extras1,/type) eq 8L) then begin
            
            if k eq 0L then bc03_extras = extras1 else begin

; remove repeated structure tag names before appending the data             
                
                oldtags = tag_names(bc03_extras)
                newtags = tag_names(extras1)

                match, oldtags, newtags, oldindx, newindx, count=count
                if count ne 0L then $
                  extras1 = struct_trimtags(extras1,except=newtags[newindx])
                
                bc03_extras = struct_addtags(bc03_extras,extras1)

            endelse

        endif
        
    endfor

endif 

if(keyword_set(vac)) then begin
    wave=bc03.wave
    airtovac, wave
    bc03.wave=wave
endif

return, bc03
end
