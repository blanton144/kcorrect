;+
; NAME:
;   k_gst_tests
; PURPOSE:
;   runs tests on GALEX+SDSS+2MASS test data
; CALLING SEQUENCE:
;   k_galex_tests
; REVISION HISTORY:
;   2005-04-07 MRB, NYU
;-
;------------------------------------------------------------------------------
pro k_gst_tests, vname=vname

if(NOT file_test('data_k_gst_tests.sav')) then begin
    calibobj=hogg_mrdfits(getenv('VAGC_REDUX')+'/object_sdss_imaging.fits',1, $
                          nrow=28800)
    galex=hogg_mrdfits(getenv('VAGC_REDUX')+'/object_galex.fits', 1, $
                       nrow=28800, columns=['e_bv', 'fuv_mag', 'fuv_magerr', $
                                            'nuv_mag', 'nuv_magerr'])
    twomass=hogg_mrdfits(getenv('VAGC_REDUX')+'/object_twomass.fits', 1, $
                         columns=['ra', 'decl', $
                                 'j_m_ext', 'j_msig_ext', 'j_flg_ext', $
                                 'h_m_ext', 'h_msig_ext', 'h_flg_ext', $
                                 'k_m_ext', 'k_msig_ext', 'k_flg_ext'], $
                         nrow=28800)
    sp=hogg_mrdfits(vagc_name('object_sdss_spectro'),1,columns='z', nrow=28800)
    
    ii=where(sp.z gt 0.01 and sp.z lt 0.6)
    sp=sp[ii]
    calibobj=calibobj[ii]
    galex=galex[ii]
    twomass=twomass[ii]
    
    save, filename='data_k_gst_tests.sav'
endif else begin
    restore, 'data_k_gst_tests.sav'
endelse

kc=gst_kcorrect(sp.z, calibobj=calibobj, twomass=twomass, galex=galex, $
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
  djs_plot, sp.z, kc[i,*], psym=3, xra=[0.009, 0.599]
  if (i mod 2) eq 1 then $
    axis,!X.CRANGE[1],!Y.CRANGE[0],yaxis=1, $
    ytitle=textoidl('y'),ycharsize=2.4 
endfor

k_end_print, pold=pold, xold=xold, yold=yold


end
