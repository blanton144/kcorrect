;+
; NAME:
;   k_lambda_to_edges
; PURPOSE:
;   convert set of pixel centers to equivalent pixels edges
; CALLING SEQUENCE:
;   lambda_edges=k_lambda_to_edges(lambda_centers)
; INPUTS:
;   lambda_centers - [N] pixel centers 
; OUPUTS:
;   lambda_edges - [N+1] pixel edges 
; REVISION HISTORY:
;   30-Apr-2003  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function k_lambda_to_edges, lambda_centers

nlambda=n_elements(lambda_centers)
lambda_edges=fltarr(nlambda+1)
lambda_edges[1:nlambda-1]= $
  0.5*(lambda_centers[0:nlambda-2]+lambda_centers[1:nlambda-1])
lambda_edges[0]=lambda_centers[0]-(lambda_edges[1]-lambda_centers[0])
lambda_edges[nlambda]=lambda_centers[nlambda-1L]+ $
  (lambda_centers[nlambda-1]-lambda_edges[nlambda-1])

return,lambda_edges

end
;------------------------------------------------------------------------------

