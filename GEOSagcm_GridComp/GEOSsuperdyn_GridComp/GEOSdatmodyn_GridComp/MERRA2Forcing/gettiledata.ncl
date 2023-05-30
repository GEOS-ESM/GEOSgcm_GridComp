; This routine finds tiles nearest to the specified alat and alon varaibles
; and writes them out to file.

;Revision History
; ???       R.Cullather - Original author
; 21Aug2020 E.Sherman   - Added to Github

;**********************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

; alat =   74.0
; alon = -155.0
;  alat =  75.
;  alon =  -140.00
;  alat =  71.00
;  alon =  -157.00

;  alat = 45.95
;  alon = -90.27

  alat = 36.00
  alon = -97.50

  g = 9.80665
  rho = 6371.220
  pi = 4. * atan(1.)

  nc_infile1 = addfile("/discover/nobackup/projects/ice/rcullath/SCM/merra2_cube_latlon.nc", "r")
  cube_lat2d = (/ nc_infile1->lat2d /)
  cube_lon2d = (/ nc_infile1->lon2d /)

  cube_lat2dll = (/ nc_infile1->lat2d_ll /)
  cube_lon2dll = (/ nc_infile1->lon2d_ll /)
  cube_lat2dur = (/ nc_infile1->lat2d_ur /)
  cube_lon2dur = (/ nc_infile1->lon2d_ur /)


  alat_tmp = cube_lat2d
  alat_tmp = alat
  alon_tmp = cube_lon2d
  alon_tmp = alon

  gcdist = gc_latlon(alat_tmp,alon_tmp,cube_lat2d,cube_lon2d,2,4)

  gcdist_1d = ndtooned(gcdist)
  dsizes_gcdist = dimsizes(gcdist)
  indices = ind_resolve( minind(gcdist_1d), dsizes_gcdist)
; print(indices)
; print(cube_lon2d(indices(0,0), indices(0,1))+" "+cube_lat2d(indices(0,0), indices(0,1)))

  ig = indices(0,1)
  jg = indices(0,0)

  cig = sprinti( "%10.0i", ig + 1)
  cjg = sprinti( "%10.0i", jg + 1)
  system("cat /discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus-NL/Icarus-NL_MERRA-2/CF0180x6C_DE1440xPE0720/CF0180x6C_DE1440xPE0720-Pfafstetter.til | nl | grep '"+cig+cjg+"' > tile.tmp")


  fsub = readAsciiTable("tile.tmp", 13, "float", 0)
  isub = readAsciiTable("tile.tmp", 13, "integer", 0)
; printVarSummary(fsub)

  tile_latsub = fsub(:,4)
  tile_lonsub = fsub(:,3)
  alat_tmp2 = tile_latsub
  alat_tmp2 = alat
  alon_tmp2 = tile_lonsub
  alon_tmp2 = alon


  gcdist2 = gc_latlon(alat_tmp2,alon_tmp2,tile_latsub,tile_lonsub,2,4)
  subind = isub(minind(gcdist2),0)
; print(subind)


  alatround = round(alat, 3)
  alonround = round(alon, 3)
  alatroundm1 = alatround - 1
  alatroundp1 = alatround + 1
  if (alatroundp1.gt.90) then
    alatroundp1 = 90
  end if
  if (alatroundm1.lt.-90) then
    alatroundm1 = -90
  end if
  alonroundm1 = alonround - 1
  alonroundp1 = alonround + 1
  if(alonroundm1.lt.-180)
    alonroundm1 = alonroundm1 + 360
  end if
  if(alonroundp1.gt.180)
    alonroundp1 = alonroundp1 - 360
  end if
  if (alatroundp1.lt.0) then
    clatp1sign = "S"
  else
    clatp1sign = "N"
  end if
  if (alatroundm1.lt.0) then
    clatm1sign = "S"
  else
    clatm1sign = "N"
  end if

  if (alonroundm1.lt.0) then
    clonm1sign = "W"
  else
    clonm1sign = "E"
  end if
  if (alonroundp1.lt.0) then
    clonp1sign = "W"
  else
    clonp1sign = "E"
  end if


  print("  1 1 1")
  print("  2")
  print("XY1x1-C_"+abs(alatroundm1)+clatm1sign+"_"+abs(alonroundm1)+clonm1sign+"_"+abs(alatroundp1)+clatp1sign+"_"+abs(alonroundp1)+clonp1sign)
  print("     1")
  print("     1")
  print("XY1x1-C_"+abs(alatroundm1)+clatm1sign+"_"+abs(alonroundm1)+clonm1sign+"_"+abs(alatroundp1)+clatp1sign+"_"+abs(alonroundp1)+clonp1sign)
  print("     1")
  print("     1")
  system("cat /discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus-NL/Icarus-NL_MERRA-2/CF0180x6C_DE1440xPE0720/CF0180x6C_DE1440xPE0720-Pfafstetter.til | head -"+subind+" | tail -1")



  itype = toint( systemfunc("cat /discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus-NL/Icarus-NL_MERRA-2/CF0180x6C_DE1440xPE0720/CF0180x6C_DE1440xPE0720-Pfafstetter.til | head -"+subind+" | tail -1 | cut -c 8-10") )
  if(itype.eq.100) then
    print("well, it's land "+itype)
  else
    print("well, it's not land "+itype)
  end if

  npfaff = toint( systemfunc("cat /discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus-NL/Icarus-NL_MERRA-2/CF0180x6C_DE1440xPE0720/CF0180x6C_DE1440xPE0720-Pfafstetter.til | cut -c 8-10 | grep 100 | wc -l") )
  print("npfaff = "+npfaff)
  ipfaff = toint( systemfunc("cat /discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus-NL/Icarus-NL_MERRA-2/CF0180x6C_DE1440xPE0720/CF0180x6C_DE1440xPE0720-Pfafstetter.til | head -"+subind+" | cut -c 8-10 | grep 100 | wc -l") )
  print("ipfaff = "+ipfaff)


  do k = 0, 13
    green_head = fbinrecread("/discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus-NLv2/Icarus-NLv2_MERRA-2/CF0180x6C_DE1440xPE0720/green_clim_180x1080.data", k * 2, 14, "float")
    ftile = fbinrecread("/discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus-NLv2/Icarus-NLv2_MERRA-2/CF0180x6C_DE1440xPE0720/green_clim_180x1080.data", k * 2 + 1, npfaff, "float")
    green_head(12) = 1.
    green_head(13) = 1.
    fbinrecwrite("green.data", k * 2, green_head)
    fbinrecwrite("green.data", k * 2 + 1, ftile(ipfaff-1))

  end do

  nirdf_recnum = fbinnumrec("/discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus-NLv2/Icarus-NLv2_MERRA-2/CF0180x6C_DE1440xPE0720/nirdf_180x1080.dat")
  nirdf_kend = nirdf_recnum / 2

  do k = 0, nirdf_kend - 1
    nirdf_head = fbinrecread("/discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus-NLv2/Icarus-NLv2_MERRA-2/CF0180x6C_DE1440xPE0720/nirdf_180x1080.dat", k * 2, 14, "float")
    ftile = fbinrecread("/discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus-NLv2/Icarus-NLv2_MERRA-2/CF0180x6C_DE1440xPE0720/nirdf_180x1080.dat", k * 2 + 1, npfaff, "float")
    nirdf_head(12) = 1.
    nirdf_head(13) = 1.
    fbinrecwrite("nirdf.dat", k * 2, nirdf_head)
    fbinrecwrite("nirdf.dat", k * 2 + 1, ftile(ipfaff-1))
  end do


  lai_recnum = fbinnumrec("/discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus-NLv2/Icarus-NLv2_MERRA-2/CF0180x6C_DE1440xPE0720/lai_clim_180x1080.data")
  lai_kend = lai_recnum / 2
; print(lai_kend)
  do k = 0, lai_kend - 1
    lai_head = fbinrecread("/discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus-NLv2/Icarus-NLv2_MERRA-2/CF0180x6C_DE1440xPE0720/lai_clim_180x1080.data", k * 2, 14, "float")
    ftile = fbinrecread("/discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus-NLv2/Icarus-NLv2_MERRA-2/CF0180x6C_DE1440xPE0720/lai_clim_180x1080.data", k * 2 + 1, npfaff, "float")
    lai_head(12) = 1.
    lai_head(13) = 1.
    fbinrecwrite("lai.data", k * 2, lai_head)
    fbinrecwrite("lai.data", k * 2 + 1, ftile(ipfaff-1))

  end do

  visdf_recnum = fbinnumrec("/discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus-NLv2/Icarus-NLv2_MERRA-2/CF0180x6C_DE1440xPE0720/visdf_180x1080.dat")
  visdf_kend = visdf_recnum / 2
; print(visdf_kend)
  do k = 0, visdf_kend - 1
    visdf_head = fbinrecread("/discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus-NLv2/Icarus-NLv2_MERRA-2/CF0180x6C_DE1440xPE0720/visdf_180x1080.dat", k * 2, 14, "float")
    ftile = fbinrecread("/discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus-NLv2/Icarus-NLv2_MERRA-2/CF0180x6C_DE1440xPE0720/visdf_180x1080.dat", k * 2 + 1, npfaff, "float")
    visdf_head(12) = 1.
    visdf_head(13) = 1.
    fbinrecwrite("visdf.dat", k * 2, visdf_head)
    fbinrecwrite("visdf.dat", k * 2 + 1, ftile(ipfaff-1))
  end do


  ndvi_recnum = fbinnumrec("/discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus-NLv2/Icarus-NLv2_MERRA-2/CF0180x6C_DE1440xPE0720/ndvi_clim_180x1080.data")
  ndvi_kend = ndvi_recnum / 2
; print(ndvi_kend)
  do k = 0, ndvi_kend - 1
    ndvi_head = fbinrecread("/discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus-NLv2/Icarus-NLv2_MERRA-2/CF0180x6C_DE1440xPE0720/ndvi_clim_180x1080.data", k * 2, 14, "float")
    ftile = fbinrecread("/discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus-NLv2/Icarus-NLv2_MERRA-2/CF0180x6C_DE1440xPE0720/ndvi_clim_180x1080.data", k * 2 + 1, npfaff, "float")
    ndvi_head(12) = 1.
    ndvi_head(13) = 1.
    fbinrecwrite("ndvi.data", k * 2, ndvi_head)
    fbinrecwrite("ndvi.data", k * 2 + 1, ftile(ipfaff-1))
  end do




  ftile = fbinrecread("/discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus-NLv2/Icarus-NLv2_MERRA-2/CF0180x6C_DE1440xPE0720/vegdyn_180x1080.dat", 0, npfaff, "float")
  fbinrecwrite("vegdyn.data", 0, ftile(ipfaff-1))


  fgcm = fbinrecread("/discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus-NLv2/Icarus-NLv2_MERRA-2/CF0180x6C_DE1440xPE0720/topo_DYN_ave_180x1080.data", 0, (/1080,180/), "float")
  fbinrecwrite("topo_dynave.data", 0, fgcm(toint(jg),toint(ig)))

  fgcm = fbinrecread("/discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus-NLv2/Icarus-NLv2_MERRA-2/CF0180x6C_DE1440xPE0720/topo_GWD_var_180x1080.data", 0, (/1080,180/), "float")
  fbinrecwrite("topo_gwdvar.data", 0, fgcm(toint(jg),toint(ig)))

  fgcm = fbinrecread("/discover/nobackup/projects/gmao/bcs_shared/legacy_bcs/Icarus-NLv2/Icarus-NLv2_MERRA-2/CF0180x6C_DE1440xPE0720/topo_TRB_var_180x1080.data", 0, (/1080,180/), "float")
  fbinrecwrite("topo_trbvar.data", 0, fgcm(toint(jg),toint(ig)))

end



