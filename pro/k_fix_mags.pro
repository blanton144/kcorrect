;+
; NAME:
;   k_fix_mags
;
; PURPOSE:
;   Fix a set of magnitudes up
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL INPUT/OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   07-Feb-2002  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_fix_mags,z,mags,magserr,maglimit,errlimit,zstep

mags=abs(mags)
magserr=abs(magserr)

; get average colors of good galaxies
goodindx=where(abs(magserr[0,*]) le errlimit[0] and $
               abs(magserr[1,*]) le errlimit[1] and $
               abs(magserr[2,*]) le errlimit[2] and $
               abs(magserr[3,*]) le errlimit[3] and $
               abs(magserr[4,*]) le errlimit[4] and $
               abs(mags[0,*]) le maglimit[0] and $
               abs(mags[1,*]) le maglimit[1] and $
               abs(mags[2,*]) le maglimit[2] and $
               abs(mags[3,*]) le maglimit[3] and $
               abs(mags[4,*]) le maglimit[4])
nk=n_elements(errlimit)
good=lonarr(n_elements(mags)/nk)
good[goodindx]=1l
nz=long((max(z)-min(z))/zstep)
zbounds=min(z)+(max(z)-min(z))*dindgen(nz+1l)/double(nz)
avgcolors=dblarr(n_elements(zbounds),nk-1l)
for i=0l, nz-1l do begin
    help,i
    indxz=where(z[goodindx] ge zbounds[i] and $
                z[goodindx] le zbounds[i+1],count)
    if(count gt 0) then begin 
        for k=0l, nk-2l do $
          avgcolors[i,k]= $
          djs_avsigclip(mags[k,goodindx[indxz]]-mags[k+1,goodindx[indxz]])
    endif else begin
        for k=0l, nk-2l do $
          avgcolors[i,k]=0.5
    endelse
endfor

; use average colors to fix bad galaxies
badindx=where(good ne 1l,count)
if(count gt 0) then begin
    help,count
    for i=0l, count-1l do begin
        if ((i mod 1000) eq 0) then help,i
        zindx=long(double(nz)*(z[i]-min(z))/(max(z)-min(z)))
        fixed=lonarr(nk)
        tmpindx=where(magserr[*,badindx[i]] le errlimit+1.e-6 and $
                      mags[*,badindx[i]] le maglimit+1.e-6 or $
                      fixed,ngood)
        if(ngood eq 0) then begin
            magserr[*,badindx[i]]=errlimit
            tmpindx=where(magserr[*,badindx[i]] le errlimit+1.e-6 and $
                          mags[*,badindx[i]] le maglimit+1.e-6 or $
                          fixed,ngood)
            if(ngood eq 0) then begin
                fixed=lonarr(nk)+1l
                tmpindx=where(magserr[*,badindx[i]] le errlimit+1.e-6 and $
                              mags[*,badindx[i]] le maglimit+1.e-6 or $
                              fixed,ngood)
            endif
        endif
        isbad=lonarr(nk)+1l
        isbad[tmpindx]=0l
        iter=0
        while (ngood lt nk) do begin
; go backwards through the bands (for ugriz case, bases mags on most
; stable bands
            for k=ngood-1l, 0l, -1l do begin
                if(tmpindx[k] gt 0) then begin
                    if(isbad[tmpindx[k]-1]) then begin
                        mags[tmpindx[k]-1l,badindx[i]]= $
                          mags[tmpindx[k],badindx[i]]+ $
                          avgcolors[zindx,tmpindx[k]-1]
                        magserr[tmpindx[k]-1l,badindx[i]]= $
                          errlimit[tmpindx[k]-1l]
                        fixed[tmpindx[k]-1]=1l
                    endif
                endif
                if(tmpindx[k] lt nk-1) then begin
                    if(isbad[tmpindx[k]+1]) then begin
                        mags[tmpindx[k]+1l,badindx[i]]= $
                          mags[tmpindx[k],badindx[i]]- $
                          avgcolors[zindx,tmpindx[k]]
                        magserr[tmpindx[k]+1l,badindx[i]]= $
                          errlimit[tmpindx[k]+1l]
                        fixed[tmpindx[k]+1]=1l
                    endif
                endif
            endfor
            tmpindx=where(magserr[*,badindx[i]] le errlimit+1.e-6 and $
                          mags[*,badindx[i]] le maglimit+1.e-6 or $
                          fixed,ngood)
            isbad=lonarr(nk)+1l
            isbad[tmpindx]=0l
            iter=iter+1
            if(iter gt 6) then stop
        end
    endfor
endif

end
;------------------------------------------------------------------------------
