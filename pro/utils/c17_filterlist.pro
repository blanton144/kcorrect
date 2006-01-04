;+
; NAME:
;   c17_filterlist
; PURPOSE:
;   return list of COMBO-17 filter names
; CALLING SEQUENCE:
;   filterlist=c17_filterlist()
; OUTPUTS:
;   filterlist - [17] list of filter names
; REVISION HISTORY:
;   03-Jan-2006  Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function c17_filterlist

filterlist='epsi_'+['U', $
                    'B', $
                    'V', $
                    'R', $
                    'I', $
                    '420m', $
                    '464m', $
                    '485m', $
                    '518m', $
                    '571m', $
                    '604m', $
                    '646m', $
                    '696m', $
                    '753m', $
                    '815m', $
                    '855m', $
                    '915m']+'.par'

return, filterlist

end
