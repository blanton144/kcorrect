;+
; NAME:
;   k_run_tests
; PURPOSE:
;   runs tests on test data
; CALLING SEQUENCE:
;   k_run_tests
; REVISION HISTORY:
;   2005-04-07 MRB, NYU
;-
;------------------------------------------------------------------------------
pro k_run_tests, vname=vname

k_sdss_tests, vname=vname
k_galex_tests, vname=vname
k_twomass_tests, vname=vname
k_gst_tests, vname=vname
k_deep_tests, vname=vname
k_goods_tests, vname=vname

end
