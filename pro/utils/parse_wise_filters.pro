;+
; NAME:
;   parse_wise_filters
; PURPOSE:
;   Read the WISE filters in raw form and output kcorrect format
; CALLING SEQUENCE:
;   parse_wise_filters
; COMMENTS:
;   Reads $KCORRECT_DIR/data/wise/RSR-W[1-4].txt and outputs:
;     wise_w1.par
;     wise_w2.par
;     wise_w3.par
;     wise_w4.par
; REVISION HISTORY:
;   14-Apr-2011 MRB NYU
;-
;------------------------------------------------------------------------------
pro parse_wise_filters


kf0= {KFILTER_WISE, lambda:0.D, pass:0., unc:0.}

for i=1L, 4L do begin
    istr=strtrim(string(i),2)
    infile= getenv('KCORRECT_DIR')+'/data/filters/wise/RSR-W'+istr+'.txt'
    outfile= getenv('KCORRECT_DIR')+'/data/filters/wise_w'+istr+'.par'
    readcol, infile, lmu, pass, unc, format='(f,f,f)'
    kf= replicate(kf0, n_elements(lmu))
    kf.lambda= lmu*10000.
    kf.pass= pass
    kf.unc= unc
    hdr = [$
            '# Units:',$
            '#  "lambda" is in Angstroms',$
            '#  "pass" is the contribution to the detector signal per photon',$
            '#',$
            '# Bandpass Name(s): W'+istr, $
            '#',$
            '# Instrument: WISE', $
            '#',$
            '# Determined by: WISE team',$
            '#',$
            '# Date of determination: sometime before 2009-Aug-04',$
            '#',$
            '# Meaning of/Reason for default column: Only column available',$
            '#',$
            '# Notes:',$
            '#  Obtained from the web site', $
            '#   http://www.astro.ucla.edu/~wright/WISE/passbands.html', $
            '#',$
            '# Mike Blanton, NYU 2011-Apr-11']
    yanny_write, outfile, ptr_new(kf), hdr=hdr
endfor

end
