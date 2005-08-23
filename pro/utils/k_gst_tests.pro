;+
; NAME:
;   k_gst_tests
; PURPOSE:
;   runs tests on GALEX+SDSS+2MASS test data
; CALLING SEQUENCE:
;   k_gst_tests
; DATA DEPENDENCIES:
;   $KCORRECT_DIR/data/test/gst_tests.fits (builds if not there)
; REVISION HISTORY:
;   2005-04-07 MRB, NYU
;-
;------------------------------------------------------------------------------
pro k_gst_tests, vname=vname

gstfile=getenv('KCORRECT_DIR')+'/data/test/gst_tests.fits'
if(NOT file_test(gstfile)) then begin
    calibobj=hogg_mrdfits(getenv('VAGC_REDUX')+'/object_sdss_imaging.fits',1, $
                          nrow=28800)
    galex=hogg_mrdfits(getenv('VAGC_REDUX')+'/object_galex.fits', 1, $
                       nrow=28800, columns=['galex_tag'])
    twomass=hogg_mrdfits(getenv('VAGC_REDUX')+'/object_twomass.fits', 1, $
                         columns=['twomass_tag', 'ra', 'decl', $
                                 'j_m_ext', 'j_msig_ext', 'j_flg_ext', $
                                 'h_m_ext', 'h_msig_ext', 'h_flg_ext', $
                                 'k_m_ext', 'k_msig_ext', 'k_flg_ext'], $
                         nrow=28800)
    sp=hogg_mrdfits(vagc_name('object_sdss_spectro'),1,columns='z', nrow=28800)
    
    ii=where(sp.z gt 0.01 and sp.z lt 0.6 and twomass.twomass_tag ge 0 and $
             galex.galex_tag ge 0, nii)
    sp=sp[ii]
    calibobj=calibobj[ii]
    galex=galex[ii]
    twomass=twomass[ii]

    if(nii gt 50000) then begin
        jj=shuffle_indx(n_elements(twomass), num_sub=50000)
        sp=sp[jj]
        calibobj=calibobj[jj]
        galex=galex[jj]
        twomass=twomass[jj]
    endif

    galex=mrdfits(getenv('VAGC_REDUX')+'/object_galex.fits', 1, $
                  row=ii[jj])
    
    cat1=create_struct(calibobj[0], $
                       sp[0], $
                       struct_trimtags(twomass[0], select_tags='*', $
                                       except_tags='ra'), $
                       struct_trimtags(galex[0], select_tags='*', $
                                       except_tags='id'))
    cat=replicate(cat1, n_elements(calibobj))
    struct_assign, sp, cat
    struct_assign, calibobj, cat, /nozero
    struct_assign, twomass, cat, /nozero
    struct_assign, galex, cat, /nozero
    mwrfits, cat, gstfile, /create
endif else begin
    cat=mrdfits(gstfile,1)
endelse

kc=gst_kcorrect(cat.z, calibobj=cat, twomass=cat, galex=cat, $
                band_shift=0., vname=vname)

k_print, filename='kcorrect.ps', $
  pold=pold, xold=xold, yold=yold, $
  axis_char_scale=2.4

!X.MARGIN=[0,2]
!Y.MARGIN=[0,0]
!X.OMARGIN=0
!Y.OMARGIN=0
!P.MULTI=[0,2,5]

for i=0, 9 do begin  
  djs_plot, cat.z, kc[i,*], psym=3, xra=[0.009, 0.599]
  if (i mod 2) eq 1 then $
    axis,!X.CRANGE[1],!Y.CRANGE[0],yaxis=1, $
    ytitle=textoidl('y'),ycharsize=2.4 
endfor

k_end_print, pold=pold, xold=xold, yold=yold


end
