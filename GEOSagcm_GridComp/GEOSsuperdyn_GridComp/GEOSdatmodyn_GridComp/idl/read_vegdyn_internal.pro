pro read_vegdyn_internal,f=f,nt=nt,veg=veg


; for MERRA 144x91 nt =45147L
;integer,parameter :: FSAT=1  !  Saturated subtile
;integer,parameter :: FTRN=2  !  Transition subtile
;integer,parameter :: FWLT=3  !  Wilting subtile
;integer,parameter :: FSNW=4  !  Snowcover subtile
;integer,parameter :: NUM_SUBTILES=4

; Vegetation type as follows:
;                  1:  BROADLEAF EVERGREEN TREES
;                  2:  BROADLEAF DECIDUOUS TREES
;                  3:  NEEDLELEAF TREES
;                  4:  GROUND COVER
;                  5:  BROADLEAF SHRUBS
;                  6:  DWARF TREES (TUNDRA)
;                  7:  BARE SOIL
;                  8:  DESERT

close,1
openr,1,/f77,f

veg={ ntiles: nt }

    a=fltarr(nt)

       ;SHORT_NAME = 'LAI_PREV'                                ,&
       ;LONG_NAME  = 'leaf_area_index_prev'                    ,&
       ;UNITS      = '1'                                       ,&
       ;DIMS       = MAPL_;DimsTileOnly                         ,&
       ;VLOCATION  = MAPL_VLocationNone                        ,&
   ;readu,1,a
   ;lai0=a
   ;veg=create_struct( veg , 'lai_prev', a )


       ;SHORT_NAME = 'LAI_NEXT'                                ,&
       ;LONG_NAME  = 'leaf_area_index_next'	              ,&
       ;UNITS      = '1'                                       ,&
       ;DIMS       = MAPL_;DimsTileOnly                         ,&
       ;VLOCATION  = MAPL_VLocationNone                        ,&
   ;readu,1,a
   ;lai2=a
   ;veg=create_struct( veg , 'lai_next', a )

       ;SHORT_NAME = 'GRN_PREV'                                ,&
       ;LONG_NAME  = 'greeness_fraction_prev'		      ,&
       ;UNITS      = '1'                                       ,&
       ;DIMS       = MAPL_;DimsTileOnly                         ,&
       ;VLOCATION  = MAPL_VLocationNone                        ,&
   ;readu,1,a
   ;grn0=a
   ;veg=create_struct( veg , 'grn_prev', a )

       ;SHORT_NAME = 'GRN_NEXT'                                ,&
       ;LONG_NAME  = 'greeness_fraction_next'                  ,&
       ;UNITS      = '1'                                       ,&
       ;DIMS       = MAPL_;DimsTileOnly                         ,&
       ;VLOCATION  = MAPL_VLocationNone                        ,&
   ;readu,1,a
   ;grn2=a
   ;veg=create_struct( veg , 'grn_next', a )


       ;SHORT_NAME = 'ITY'                                     ,&
       ;LONG_NAME  = 'vegetation_type'			      ,&
       ;UNITS      = '1'                                       ,&
       ;DIMS       = MAPL_;DimsTileOnly                         ,&
       ;VLOCATION  = MAPL_VLocationNone                        ,&
    readu,1,a
    ity=a
    veg=create_struct( veg , 'ity', a )

return
end
