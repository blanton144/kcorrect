function lf_erf,x
val=1.-lf_erfc(x)
indx=where(x lt 0.,count)
if(count gt 0) then val[indx]=-val[indx]
return,val
end
