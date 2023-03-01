function read_hdf_sd,filename,vars=varname,list=list

FileID=HDF_SD_START(filename,/READ)
HDF_SD_FILEINFO,FileID,nvar,nattr

nvarname = n_elements(varname)
if (nvarname EQ 0) then varname=''
a=create_struct('FileName',filename)

var=strarr(nvar)
typ=strarr(nvar)
ndim=intarr(nvar)
dim=strarr(nvar,10)
varlen = 0
typlen = 0
dimnum = 0
ndmax  = 0
index  = 0
ivar   = 0
while (index LT nvar) do begin
  sds_id=HDF_SD_SELECT(FileID,index)
  HDF_SD_GETINFO, sds_id, NAME=name, DIMS=dims, NDIMS=ndims, TYPE=type
  ind = where(name EQ varname, cnt)
  if (cnt GT 0 OR nvarname EQ 0) then begin
    if (ndims GT 1 OR dims[0] GT 0) then begin
      var[ivar] = name
      typ[ivar] = type
      ndim[ivar] = ndims
      nd=0
      flag=0
      for i=0,ndims-1 do begin
        dim[ivar,i] = dims[i]
        nd=nd+alog10(dim[ivar,i])+3
        if (dims[i] EQ 0) then flag=1
      endfor
      if (flag EQ 0) then begin
        if (nd GT ndmax) then ndmax=nd
        if (strlen(var[ivar]) GT varlen) then varlen=strlen(var[ivar])
        if (strlen(typ[ivar]) GT typlen) then typlen=strlen(typ[ivar])
        if (ndims GT dimnum) then dimnum=ndims
        HDF_SD_GETDATA, sds_id, data
        name = var[ivar]
        c=['-','.','/','(',')','+']
        for n=0,n_elements(c)-1 do begin
          pos=strpos(name,c[n])
          while (pos GE 0) do begin
            strput,name,' ',pos
            pos=strpos(name,c[n])
          endwhile
        endfor
        name = strcompress(name,/remove_all)
        fchar=strmid(name,0,1)
        if (fchar GE '0' AND fchar LE '9') then name='Var_'+name
  
        svar=string(var[ivar],ivar,format='(a,"_",i3.3)')
        a = create_struct(a, name, data)
        ivar = ivar + 1
      endif
    endif
  endif
  index = index + 1
endwhile
nvar = ivar

if (nvar GT 0) then begin
  s1=string(fix(alog10(nvar))+1,format='("i",i1)')
  s2='2x,a'
  s3=string(typlen+1,format='("a",i2.2)')
  ndmax=0
  if (keyword_set(list)) then openw,2,'var.list'
  for index=0,nvar-1 do begin
    fmt='("SD[",'+s1+',"]: ",'+s3+','+s2
    s4=strarr(ndim[index])
    for i=0,ndim[index]-1 do begin
      s4[i] = string(alog10(dim[index,i])+1,format='("i",i1)')
      fmt=fmt+',"[",'+s4[i]+',"]"'
    endfor
    fmt=fmt+')'
    if (keyword_set(list)) then $
      printf,2,index,typ[index],var[index],dim[index,0:ndim[index]-1],format=fmt
  endfor
endif else begin
  a = -1
endelse

if (keyword_set(list)) then close,2
HDF_SD_END, FileID
return,a
end

