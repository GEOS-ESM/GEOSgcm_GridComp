function read_mod10a1_hdf, filename

; code to read MODIS snow albedo files, like: 
; MOD10A1.A2020365.h34v08.061.2021006142942.hdf
; or
; MOD10A1.A2022002.h12v11.006.2022004045329.hdf

; to be used in the stitching process

; Created March 2022 Biljana Orescanin NASA@SSAI

vars2read= [                                       $
             'NDSI_Snow_Cover'                   , $
             'NDSI_Snow_Cover_Basic_QA'          , $
             'NDSI_Snow_Cover_Algorithm_Flags_QA', $
             'NDSI'                              , $
             'Snow_Albedo_Daily_Tile'            , $
             'orbit_pnt'                         , $
             'granule_pnt'                         $
               ]

; read the hdf file; no scale and offsets, so it's ok to use read_hdf_sd
; note: read_hdf_sd does not allow for choosing which variables to read,
;       so I read them all
d=read_hdf_sd(filename)

return, d

stop
end
