;+
; NAME:
;   k_lambda_to_edges
; PURPOSE:
;   convert set of pixel centers to equivalent pixels edges
; CALLING SEQUENCE:
;   lambda_edges=k_lambda_to_edges(lambda_centers)
; INPUTS:
;   lambda_centers - pixel centers 
; OPTIONAL INPUTS:
; OUTPUTS:
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   30-Apr-2003  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function k_lambda_to_edges, lambda_centers

nlambda=n_elements(lambda_centers)
lambda_edges=fltarr(nlambda+1)
dlambda=fltarr(nlambda)
dlambda[0]=0.5*(lambda_centers[1]-lambda_centers[0])
dlambda[nlambda-1]=0.5*(lambda_centers[nlambda-1]-lambda_centers[nlambda-2])
dlambda[1:nlambda-2]=0.5*(lambda_centers[2:nlambda-1]- $
                          lambda_centers[0:nlambda-3])
lambda_edges[0:nlambda-1]=lambda_centers-0.5*dlambda
lambda_edges[1:nlambda]=lambda_centers+0.5*dlambda

return,lambda_edges

end
;------------------------------------------------------------------------------

