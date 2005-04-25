;+
; NAME:
;   k_sdss_tests
; PURPOSE:
;   runs tests on SDSS test data
; CALLING SEQUENCE:
;   k_sdss_tests
; REVISION HISTORY:
;   2005-04-07 MRB, NYU
;-
;------------------------------------------------------------------------------
pro k_sdss_tests, vname=vname

k_sdss_tests_lrg, vname=vname
k_sdss_tests_main, vname=vname
k_sdss_tests_generic, vname=vname

end
