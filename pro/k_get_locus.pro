pro k_get_locus,points,locus,nlocus,reinit=reinit

nsig=3.
maxiter=80
sigrej=3.
mindiff=0.0001

; Create and sort covariance
ndim=(size(points))[1]
mean=dblarr(ndim)
covar=dblarr(ndim,ndim)
for i=0, ndim-1 do $
  mean[i]=djs_mean(points[i,*])
for i=0, ndim-1 do $
  for j=0, ndim-1 do $
  covar[i,j]=djsig(points[i,*]-mean[i],sigrej=sigrej)* $
  djsig(points[j,*]-mean[j],sigrej=sigrej)
diag=lindgen(ndim)
var=covar[diag,diag]
dim=sort(var)

if (keyword_set(reinit)) then begin
    locus=dblarr(ndim,nlocus)
    locus[dim[ndim-1],*]=mean[dim[ndim-1]]-nsig*sqrt(var[dim[ndim-1]])+ $
      (2.*nsig*sqrt(var[dim[ndim-1]]))*(dindgen(nlocus)+0.5)/double(nlocus)
    for i=0,ndim-2 do begin
        locus[dim[i],*]=mean[dim[i]]
    endfor
endif 

niter=0
device,pseudo_color=8
!p.multi=[3,1,3]
erase
plot,points[0,*],points[1,*],psym=3
oplot,locus[0,*],locus[1,*],color=255,thick=2
oplot,locus[0,*],locus[1,*],color=255,psym=4
plot,points[0,*],points[2,*],psym=3
oplot,locus[0,*],locus[2,*],color=255,thick=2
oplot,locus[0,*],locus[2,*],color=255,psym=4
plot,points[1,*],points[2,*],psym=3
oplot,locus[1,*],locus[2,*],color=255,thick=2
oplot,locus[1,*],locus[2,*],color=255,psym=4
klog,locus
;stop
while (niter lt maxiter) do begin
    uselocus=lonarr(nlocus)+1
    for i=0, nlocus-1 do begin
        if(uselocus[i]) then begin
; find tangent vector
            im1=i-1
            if(im1 lt 0) then im1=0
            ip1=i+1
            if(ip1 gt nlocus-1) then ip1=nlocus-1
            tangent=locus[*,ip1]-locus[*,im1]
            tangent=tangent/sqrt(total(tangent^2,/double))
            
; find points within region of tangent vector
            projlocus=(tangent##transpose(locus[*,i]))[0]
            projlocusm1=(tangent##transpose(locus[*,im1]))[0]
            projlocusp1=(tangent##transpose(locus[*,ip1]))[0]
            projlocusmid=0.5*(projlocusm1+projlocusp1)
            if(projlocusp1-projlocusm1 lt mindiff) then begin
                projlocusp1=projlocusmid+0.5*mindiff
                projlocusm1=projlocusmid-0.5*mindiff
            endif
            projpoints=transpose(tangent#points)
            indx=where(projpoints gt projlocusm1 and $
                       projpoints lt projlocusp1,count)
            
; find their mean
            if(count gt 5) then begin
                for j=0, ndim-1 do begin
                    tmpmean=djs_mean(points[j,indx])
                    sig=djsig(points[j,indx]-tmpmean)
                    within=where(points[j,indx]-tmpmean lt 5.1*sig)
                    locus[j,i]=djs_mean(points[j,indx[within]])
                endfor
                locus[*,i]=locus[*,i]+tangent* $
                  (projlocus-(tangent##transpose(locus[*,i]))[0])
            endif else begin
                uselocus[i]=0
            endelse
        endif
    endfor 
    locus=locus[*,where(uselocus gt 0)]
    nlocus=n_elements(locus)/ndim
    klog,nlocus
        
    device,pseudo_color=8
    !p.multi=[3,1,3]
    erase
    plot,points[0,*],points[1,*],psym=3
    oplot,locus[0,*],locus[1,*],color=255,thick=2
    oplot,locus[0,*],locus[1,*],color=255,psym=4
    plot,points[0,*],points[2,*],psym=3
    oplot,locus[0,*],locus[2,*],color=255,thick=2
    oplot,locus[0,*],locus[2,*],color=255,psym=4
    plot,points[1,*],points[2,*],psym=3
    oplot,locus[1,*],locus[2,*],color=255,thick=2
    oplot,locus[1,*],locus[2,*],color=255,psym=4
    klog,locus
    ;stop
    
    niter=niter+1
end

end

