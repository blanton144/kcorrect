;+
; NAME:
;   k_read_tbl
; PURPOSE:
;   Read 2MASS style table into IDL structure
; CALLING SEQUENCE:
;   tbl= k_read_tbl(filename)
; INPUTS:
;   filename - file name
; BUGS:
;   not well commented
; REVISION HISTORY:
;   Spring 2003, Written Malcolm Britton
;-
;------------------------------------------------------------------------------
FUNCTION type2var, type, val

IF(strmatch(type,'d*') EQ 1) THEN $
  RETURN,double(val) $
ELSE IF(strmatch(type,'i*') EQ 1) THEN $
  RETURN,long(val) $
else if(strmatch(type,'f*') eq 1) then $
  return,float(val) $
else if(strmatch(type,'r*') eq 1) then $
  return,float(val) $
else if(strmatch(type,'c*') eq 1) then $
  return,string(val) $
else $
  message, 'No registered type '+type

end
;
FUNCTION k_read_tbl,file,chunksize=chunksize

; defaults
if(NOT keyword_set(chunksize)) then chunksize=1000

openr,unit,file,/get_lun

line=''
counter=0
stopper=0

WHILE  (stopper EQ 0) DO BEGIN
  readf,unit,line
    IF (STRMID(line,0,1) EQ '\') THEN $
      counter=counter+1 $
      ELSE $
      IF (STRMID(line,0,1) EQ '|') THEN BEGIN
        names=strtrim(strsplit(line,'|',/extract),2) 
        limits=strsplit(line,'|',length=length) 
        ncols=n_elements(names) 
        counter=counter+1       
        readf,unit,line 
        types=strtrim(strsplit(line,'|',/extract),2) 
        counter=counter+1       
        readf,unit,line 
        dimensions=strtrim(strsplit(line,'|',/extract),2) 
        counter=counter+1       
        readf,unit,line 
        counter=counter+1 
      ENDIF ELSE BEGIN 
        stopper=1 
      ENDELSE
ENDWHILE



; now build structure
FOR i=0, n_elements(names)-1 DO BEGIN
    dummy=type2var(types[i],'')
    IF(n_tags(instr1) EQ 0) THEN $
      instr1=create_struct(names[i],dummy) $
    ELSE $
      instr1=create_struct(instr1,names[i],dummy) 
ENDFOR

; create table 
nelem=numlines(file)-counter
instr=replicate(instr1,nelem)

; now read in table
chunk=strarr(ncols,chunksize)
nchunks=long(ceil(double(nelem)/double(chunksize)))
FOR i=0L, nchunks-1L DO BEGIN
    startchunk=i*chunksize
    endchunk=(((i+1L)*chunksize)<nelem)-1L
    ninchunk=endchunk-startchunk+1L
    splog,'startchunk= '+string(startchunk)+'; endchunk= '+string(endchunk)
    FOR j=0L, ninchunk-1L DO BEGIN
      IF (i EQ 0) AND (j EQ 0) THEN $
        chunk[*,0]=strmid(line,limits-1,length+1) $
      ELSE $
        readf,unit,line
        chunk[*,j]=strmid(line,limits-1,length+1)
    ENDFOR
    FOR j=0L,ncols-1L DO BEGIN
        nullindx=where(strtrim(chunk[j,0L:ninchunk-1L],2) eq 'null',nullcount)
        if(nullcount gt 0) then chunk[j,nullindx]='0'
        instr[startchunk:endchunk].(j)= $
          transpose(type2var(types[j],chunk[j,0L:ninchunk-1L]))
    ENDFOR
ENDFOR
free_lun,unit

RETURN,instr

STOP

END
