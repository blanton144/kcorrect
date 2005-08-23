;+
; NAME:
;   datasweep_lrg_photoz
; PURPOSE:
;   sweep the data sweeps and do LRG photoz's with everything 
; CALLING SEQUENCE:
;   datasweep_lrg_photoz [, /noclobber]
; OPTIONAL KEYWORDS:
;   /noclobber - don't overwrite 
; REVISION HISTORY:
;   07-Apr-2005  Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro datasweep_lrg_photoz, noclobber=noclobber

cd, getenv('PHOTO_SWEEP')

flist = file_search('calibObj0*-gal.fits.gz',/test_regular,count=count)
splog, 'Number of files found :', count

for i=0L, n_elements(flist)-1L do begin
    rname=strmid(flist[i],9,8)
    outname='lrg-photoz-'+rname+'-gal.fits'

    obj = mrdfits(flist[i],1)
    sdss_recalibrate, obj
    
    photoz=sdss_kphotoz(calibobj='obj', /lrg, chi2=chi2, rmaggies=rmaggies, $
                       omaggies=omaggies)

    outstr1={'photoz', 0., $
             'chi2', 0., $
             'omaggies':fltarr(5), $
             'rmaggies':fltarr(5) }
    photozstr=replicate(outstr1, n_elements(obj))
    photozstr.photoz=photoz
    photozstr.chi2=chi2
    photozstr.rmaggies=rmaggies
    photozstr.omaggies=omaggies
    
    hdr=['']
    sxaddpar, hdr, 'K_VERS', k_version(), 'kcorrect version'
    mwrfits, photozstr, outname, hdr, /create
    
endfor

end
