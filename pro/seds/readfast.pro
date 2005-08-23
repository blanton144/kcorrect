;+
; NAME:
;	READFAST
;
; PURPOSE:
;	Routine to read in column-formated data very quickly.
;
; INPUTS:
;	filename - name of the file to be read
;
; OPTIONAL INPUTS:
;	skipline - number of lines to skip at the top of the file
;	
; KEYWORD PARAMETERS:
;	double - optionally read the data in as double precision
;
; OUTPUTS:
;	data - floating point 2D array containing the columns and rows
;              in the data file
;
; OPTIONAL OUTPUTS:
;	header - string array containing the lines skipped at the
;                beginning of the file
;	ncols  - number of columns read
;	nlines - number of lines read
;
; COMMON BLOCKS:
;	None.
;
; COMMENTS: 
;	This program will only read floating (or double) precision
;	data, and the file must be column-formated.  It will read in
;	hundreds of thousands of lines in a matter of seconds after
;	determining how many lines and columns there are in the file.
;	Note that the program will not skip blank lines nor comment
;	lines unless explicitly told to do so with the skipline
;	option.  The idea is speed over flexibility. 
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 May 24, UC Berkeley
;	jm01jan16uofa - added nlines keyword to optionally return the
;                       number of lines in the file; added error
;                       checking 
;	jm01jan25uofa - routine figures out how many columns there are
;                       in the file; added double keyword; documented,
;                       cleaned up, and optimized
;-

pro readfast, filename, data, header, double=double, skipline=skipline, $
  ncols=ncols, nlines=nlines

;	t1 = systime(1)

; some error checking

	on_error, 2 ; return to main program
        if n_params() eq 0L then begin 
           print, 'Syntax: readfast, filename, data, header, skipline=skipline, ncols=ncols, nlines=nlines'
           return
        endif
        if not keyword_set(filename) then message, 'Please specify a file to read.'
        
        file = file_search(filename,count=nfiles) ; check that the file exists
        if nfiles eq 1L then if file[0] eq '' then message, 'File not found.'
        if nfiles gt 1L then message, 'Multiple files found.'
           
; how many lines are there?

        openr, lun, filename, /get_lun     ; open the file

        nlines = 0L
        spawn, ['wc -l '+filename], wcount 
        reads, wcount, nlines              ; file length

        if keyword_set(skipline) then begin

            nlines = nlines-skipline
            header = strarr(skipline)
            readf, lun, header	           ; read the header

        endif

; read the next line to figure out how many columns there are

        position = (fstat(lun)).cur_ptr ; where is the file pointer?

        dummy = ' '
        readf, lun, dummy
        ncols = n_elements(strsplit(dummy,' ',/extract)) ; parse by blank spaces

        point_lun, lun, position ; rewind
        
        if keyword_set(double) then data = dblarr(ncols,nlines) else data = fltarr(ncols,nlines)
        readf, lun, data         ; read the data

        free_lun, lun

;       print, systime(1)-t1

return
end

