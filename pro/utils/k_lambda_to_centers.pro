;+
; NAME:
;   k_lambda_to_centers
; PURPOSE:
;   convert set of pixel edges to equivalent pixels centers
; CALLING SEQUENCE:
;   lambda_centers=k_lambda_to_centers(lambda_edges)
; INPUTS:
;   lambda_edges - [N+1] pixel edges 
; OUPUTS:
;   lambda_centers - [N] pixel centers 
; REVISION HISTORY:
;   05-Apr-2005  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function k_lambda_to_centers, lambda_edges

nlambda=n_elements(lambda_edges)-1L
lambda_centers=0.5*(lambda_edges[0:nlambda-1]+lambda_edges[1:nlambda])

return,lambda_centers

end
;------------------------------------------------------------------------------

