pro k_end_print,pold=pold,xold=xold,yold=yold

device,/close

if(keyword_set(pold)) then !P=pold
if(keyword_set(xold)) then !X=xold
if(keyword_set(yold)) then !Y=yold

set_plot,'x'

end
