;+
; NAME:
;   k_prepare_fit
; PURPOSE:
;   create input files for fit from mmatrix
; CALLING SEQUENCE:
;   k_prepare_fit [, mmatrix=mmatrix, filterlist=filterlist, /nospec ]
; OPTIONAL INPUTS:
;   mmatrix - name of input mmatrix file (default: 'mmatrix.fits') 
;   filterlist - list of filters to use (default to all listed in infile)
; OPTIONAL KEYWORDS:
;   /nospec - don't use the spectra
; COMMENTS:
; BUGS:
;   needs to include emission line predictions   
; REVISION HISTORY:
;   29-Jul-2004  Michael Blanton (NYU)
;-
;------------------------------------------------------------------------------
pro k_prepare_fit, infile=infile, filterlist=filterlist, nospec=nospec, $
                   

if(NOT keyword_set(infile)) then infile='mmatrix.fits'

grid=mrdfits(infile, 1)
avloglam=mrdfits(infile, 2)
tauv=mrdfits(infile, 3)
met=mrdfits(infile, 4)
age=mrdfits(infile, 5)
filterlist=mrdfits(infile, 6)
zf=mrdfits(infile, 7)

mwrfits, outgrid, outfile, /create
mwrfits, avloglam,outfile 
mwrfits, tauv,outfile 
mwrfits, met, outfile
mwrfits, age, outfile
mwrfits, filterlist, outfile
mwrfits, zf, outfile


end
