pro read_catch_internal,f=f,nt=nt,cas=cas

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


cas={ ntiles: nt }

close,1
openr,1,/f77,f

    a=fltarr(nt)
    aa=fltarr(nt,4)

    ;LONG_NAME          = 'topo_baseflow_param_1'     ,&
    ;UNITS              = '1'                         ,&
    ;SHORT_NAME         = 'BF1'                       ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    bf1=a
    cas=create_struct( cas , 'basefl1', a )
 
    ;LONG_NAME          = 'topo_baseflow_param_2'     ,&
    ;UNITS              = '1'                         ,&
    ;SHORT_NAME         = 'BF2'                       ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    bf2=a
    cas=create_struct( cas , 'basefl2', a )
 
    ;LONG_NAME          = 'topo_baseflow_param_3'     ,&
    ;UNITS              = '1'                         ,&
    ;SHORT_NAME         = 'BF3'                       ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    bf3=a
    cas=create_struct( cas , 'basefl3', a )

    ;LONG_NAME          = 'max_rootzone_water_content',&
    ;UNITS              = 'kg m-2 s-1'                ,&
    ;SHORT_NAME         = 'VGWMAX'                    ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    vgwmax=a
    cas=create_struct( cas , 'vgwmax', a )
 
    ;LONG_NAME          = 'moisture_threshold'        ,&
    ;UNITS              = 'kg m-2 s-1'                ,&
    ;SHORT_NAME         = 'CDCR1'                     ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    cdcr1=a
    cas=create_struct( cas , 'cdcr1', a )
 
    ;LONG_NAME          = 'max_water_content_unsat_zone',&
    ;UNITS              = 'kg m-2 s-1'                ,&
    ;SHORT_NAME         = 'CDCR2'                     ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    cdcr2=a
    cas=create_struct( cas , 'cdcr2', a )
 
    ;LONG_NAME          = 'saturated_matrix_potential',&
    ;UNITS              = 'm'                         ,&
    ;SHORT_NAME         = 'PSIS'                      ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    psis=a
    cas=create_struct( cas , 'psis', a )
 
    ;LONG_NAME          = 'clapp_hornberger_b'        ,&
    ;UNITS              = 'm'                         ,&
    ;SHORT_NAME         = 'BEE'                       ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    bee=a
    cas=create_struct( cas , 'bee', a )
 
    ;LONG_NAME          = 'soil_porosity'             ,&
    ;UNITS              = '1'                         ,&
    ;SHORT_NAME         = 'POROS'                     ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    poros=a
    cas=create_struct( cas , 'poros', a )

    ;LONG_NAME          = 'wetness_at_wilting_point'  ,&
    ;UNITS              = '1'                         ,&
    ;SHORT_NAME         = 'WPWET'                     ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    wpwet=a
    cas=create_struct( cas , 'wpwet', a )

    ;LONG_NAME          = 'sfc_sat_hydraulic_conduct' ,&
    ;UNITS              = 'm s-1'                     ,&
    ;SHORT_NAME         = 'COND'                      ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    cond=a
    cas=create_struct( cas , 'cond', a )
 
    ;LONG_NAME          = 'vertical_transmissivity'   ,&
    ;UNITS              = 'm-1'                       ,&
    ;SHORT_NAME         = 'GNU'                       ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    gnu=a
    cas=create_struct( cas , 'gnu', a )
 
    ;LONG_NAME          = 'wetness_param_1'           ,&
    ;UNITS              = '1'                         ,&
    ;SHORT_NAME         = 'ARS1'                      ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    ars1=a
    cas=create_struct( cas , 'ars1', a )

    ;LONG_NAME          = 'wetness_param_2'           ,&
    ;UNITS              = '1'                         ,&
    ;SHORT_NAME         = 'ARS2'                      ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    ars2=a
    cas=create_struct( cas , 'ars2', a )
 
    ;LONG_NAME          = 'wetness_param_3'           ,&
    ;UNITS              = '1'                         ,&
    ;SHORT_NAME         = 'ARS3'                      ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    ars3=a
    cas=create_struct( cas , 'ars3', a )

    ;LONG_NAME          = 'shape_param_1'             ,&
    ;UNITS              = '1'                         ,&
    ;SHORT_NAME         = 'ARA1'                      ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    ara1=a
    cas=create_struct( cas , 'ara1', a )

    ;LONG_NAME          = 'shape_param_2'             ,&
    ;UNITS              = '1'                         ,&
    ;SHORT_NAME         = 'ARA2'                      ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    ara2=a
    cas=create_struct( cas , 'ara2', a )
 
    ;LONG_NAME          = 'shape_param_3'             ,&
    ;UNITS              = '1'                         ,&
    ;SHORT_NAME         = 'ARA3'                      ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    ara3=a
    cas=create_struct( cas , 'ara3', a )
 
    ;LONG_NAME          = 'shape_param_4'             ,&
    ;UNITS              = '1'                         ,&
    ;SHORT_NAME         = 'ARA4'                      ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    ara4=a
    cas=create_struct( cas , 'ara4', a )

    ;LONG_NAME          = 'min_theta_1'               ,&
    ;UNITS              = '1'                         ,&
    ;SHORT_NAME         = 'ARW1'                      ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    arw1=a
    cas=create_struct( cas , 'arw1', a )
 
    ;LONG_NAME          = 'min_theta_2'               ,&
    ;UNITS              = '1'                         ,&
    ;SHORT_NAME         = 'ARW2'                      ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    arw2=a
    cas=create_struct( cas , 'arw2', a )

    ;LONG_NAME          = 'min_theta_3'               ,&
    ;UNITS              = '1'                         ,&
    ;SHORT_NAME         = 'ARW3'                      ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    arw3=a
    cas=create_struct( cas , 'arw3', a )

    ;LONG_NAME          = 'min_theta_4'               ,&
    ;UNITS              = '1'                         ,&
    ;SHORT_NAME         = 'ARW4'                      ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    arw4=a
    cas=create_struct( cas , 'arw4', a )

    ;LONG_NAME          = 'water_transfer_1'          ,&
    ;UNITS              = '1'                         ,&
    ;SHORT_NAME         = 'TSA1'                      ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    tsa1=a
    cas=create_struct( cas , 'tsa1', a )

    ;LONG_NAME          = 'water_transfer_2'          ,&
    ;UNITS              = '1'                         ,&
    ;SHORT_NAME         = 'TSA2'                      ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    tsa2=a
    cas=create_struct( cas , 'tsa2', a )

    ;LONG_NAME          = 'water_transfer_3'          ,&
    ;UNITS              = '1'                         ,&
    ;SHORT_NAME         = 'TSB1'                      ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    tsb1=a      
    cas=create_struct( cas , 'tsb1', a )
 
    ;LONG_NAME          = 'water_transfer_4'          ,&
    ;UNITS              = '1'                         ,&
    ;SHORT_NAME         = 'TSB2'                      ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    tsb2=a
    cas=create_struct( cas , 'tsb2', a )

    ;LONG_NAME          = 'soil_param_1'              ,&
    ;UNITS              = '1'                         ,&
    ;SHORT_NAME         = 'ATAU'                      ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    atau=a
    cas=create_struct( cas , 'atau', a )

    ;LONG_NAME          = 'soil_param_2'              ,&
    ;UNITS              = '1'                         ,&
    ;SHORT_NAME         = 'BTAU'                      ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    btau=a
    cas=create_struct( cas , 'btau', a )

    ;LONG_NAME          = 'vegetation_type'           ,&
    ;UNITS              = '1'                         ,&
    ;SHORT_NAME         = 'ITY'                       ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    ity=a
    cas=create_struct( cas , 'ity', a )

    ;LONG_NAME          = 'canopy_temperature'        ,&
    ;UNITS              = 'K'                         ,&
    ;SHORT_NAME         = 'TC'                        ,&
    ;DIMS               = MAPL_DimsTileTile           ,&
    ;NUM_SUBTILES       = NUM_SUBTILES                ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,aa
    tc=aa
    cas=create_struct( cas , 'tc', aa )
 
    ;LONG_NAME          = 'canopy_specific_humidity'  ,&
    ;UNITS              = '1'                         ,&
    ;SHORT_NAME         = 'QC'                        ,&
    ;DIMS               = MAPL_DimsTileTile           ,&
    ;NUM_SUBTILES       = NUM_SUBTILES                ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,aa
    qc=aa
    cas=create_struct( cas , 'qc', aa )

    ;LONG_NAME          = 'interception_reservoir_capac',&
    ;UNITS              = 'kg m-2'                    ,&
    ;SHORT_NAME         = 'CAPAC'                     ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    capac=a
    cas=create_struct( cas , 'capac', a )
 
    ;LONG_NAME          = 'catchment_deficit'         ,&
    ;UNITS              = 'kg m-2'                    ,&
    ;SHORT_NAME         = 'CATDEF'                    ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    catdef=a
    cas=create_struct( cas , 'catdef', a )
 
    ;LONG_NAME          = 'root_zone_excess'          ,&
    ;UNITS              = 'kg m-2'                    ,&
    ;SHORT_NAME         = 'RZEXC'                     ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    rzexc=a
    cas=create_struct( cas , 'rzexc', a )
 
    ;LONG_NAME          = 'surface_excess'            ,&
    ;UNITS              = 'kg m-2'                    ,&
    ;SHORT_NAME         = 'SRFEXC'                    ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    srfexc=a
    cas=create_struct( cas , 'srfexc', a )
 
    ;LONG_NAME          = 'soil_heat_content_layer_1' ,&
    ;UNITS              = 'J m-2'                     ,&
    ;SHORT_NAME         = 'GHTCNT1'                   ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    ghtcn1=a
    cas=create_struct( cas , 'ghtcn1', a )
 
    ;LONG_NAME          = 'soil_heat_content_layer_2' ,&
    ;UNITS              = 'J_m-2'                     ,&
    ;SHORT_NAME         = 'GHTCNT2'                   ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    ghtcn2=a
    cas=create_struct( cas , 'ghtcn2', a )
 
    ;LONG_NAME          = 'soil_heat_content_layer_3' ,&
    ;UNITS              = 'J m-2'                     ,&
    ;SHORT_NAME         = 'GHTCNT3'                   ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    ghtcn3=a
    cas=create_struct( cas , 'ghtcn3', a )

    ;LONG_NAME          = 'soil_heat_content_layer_4' ,&
    ;UNITS              = 'J m-2'                     ,&
    ;SHORT_NAME         = 'GHTCNT4'                   ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    ghtcn4=a
    cas=create_struct( cas , 'ghtcn4', a )
 
    ;LONG_NAME          = 'soil_heat_content_layer_5' ,&
    ;UNITS              = 'J m-2'                     ,&
    ;SHORT_NAME         = 'GHTCNT5'                   ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    ghtcn5=a
    cas=create_struct( cas , 'ghtcn5', a )

    ;LONG_NAME          = 'soil_heat_content_layer_6' ,&
    ;UNITS              = 'J m-2'                     ,&
    ;SHORT_NAME         = 'GHTCNT6'                   ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    ghtcn6=a
    cas=create_struct( cas , 'ghtcn6', a )
 
    ;LONG_NAME          = 'mean_catchment_temp_incl_snw',&
    ;UNITS              = 'K'                         ,&
    ;SHORT_NAME         = 'TSURF'                     ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    tsurf=a
    cas=create_struct( cas , 'tsurf', a )
 
    ;LONG_NAME          = 'snow_mass_layer_1'         ,&
    ;UNITS              = 'kg m-2'                    ,&
    ;SHORT_NAME         = 'WESNN1'                    ,&
    ;FRIENDLYTO         = trim(COMP_NAME)             ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    wsenn1=a
    cas=create_struct( cas , 'wsenn1', a )
 
    ;LONG_NAME          = 'snow_mass_layer_2'         ,&
    ;UNITS              = 'kg m-2'                    ,&
    ;SHORT_NAME         = 'WESNN2'                    ,&
    ;FRIENDLYTO         = trim(COMP_NAME)             ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    wsenn2=a
    cas=create_struct( cas , 'wsenn2', a )

    ;LONG_NAME          = 'snow_mass_layer_3'         ,&
    ;UNITS              = 'kg m-2'                    ,&
    ;SHORT_NAME         = 'WESNN3'                    ,&
    ;FRIENDLYTO         = trim(COMP_NAME)             ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    wsenn3=a
    cas=create_struct( cas , 'wsenn3', a )
 
    ;LONG_NAME          = 'heat_content_snow_layer_1' ,&
    ;UNITS              = 'J m-2'                     ,&
    ;SHORT_NAME         = 'HTSNNN1'                   ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    htsnnn1=a
    cas=create_struct( cas , 'htsnnn1', a )
 
    ;LONG_NAME          = 'heat_content_snow_layer_2' ,&
    ;UNITS              = 'J m-2'                     ,&
    ;SHORT_NAME         = 'HTSNNN2'                   ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    htsnnn2=a
    cas=create_struct( cas , 'htsnnn2', a )

    ;LONG_NAME          = 'heat_content_snow_layer_3' ,&
    ;UNITS              = 'J m-2'                     ,&
    ;SHORT_NAME         = 'HTSNNN3'                   ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    htsnnn3=a
    cas=create_struct( cas , 'htsnnn3', a )
 
    ;LONG_NAME          = 'snow_depth_layer_1'        ,&
    ;UNITS              = 'm'                         ,&
    ;SHORT_NAME         = 'SNDZN1'                    ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    sndzn1=a
    cas=create_struct( cas , 'sndzn1', a )
 
    ;LONG_NAME          = 'snow_depth_layer_2'        ,&
    ;UNITS              = 'm'                         ,&
    ;SHORT_NAME         = 'SNDZN2'                    ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    sndzn2=a
    cas=create_struct( cas , 'sndzn2', a )
 
    ;LONG_NAME          = 'snow_depth_layer_3'        ,&
    ;UNITS              = 'm'                         ,&
    ;SHORT_NAME         = 'SNDZN3'                    ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,a
    sndzn3=a
    cas=create_struct( cas , 'sndzn3', a )
 
    ;LONG_NAME          = 'surface_heat_exchange_coefficient',&
    ;UNITS              = 'kg m-2 s-1'                ,&
    ;SHORT_NAME         = 'CH'                        ,&
    ;DIMS               = MAPL_DimsTileTile           ,&
    ;NUM_SUBTILES       = NUM_SUBTILES                ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&       
readu,1,aa
    ch=aa
    cas=create_struct( cas , 'ch', aa )
 
    ;LONG_NAME          = 'surface_momentum_exchange_coefficient',&
    ;UNITS              = 'kg m-2 s-1'                ,&
    ;SHORT_NAME         = 'CM'                        ,&
    ;DIMS               = MAPL_DimsTileTile           ,&
    ;NUM_SUBTILES       = NUM_SUBTILES                ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,aa
    cm=aa
    cas=create_struct( cas , 'cm', aa )
 
    ;LONG_NAME          = 'surface_moisture_exchange_coffiecient',&
    ;UNITS              = 'kg m-2 s-1'                ,&
    ;SHORT_NAME         = 'CQ'                        ,&
    ;DIMS               = MAPL_DimsTileTile           ,&
    ;NUM_SUBTILES       = NUM_SUBTILES                ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,aa
    cq=aa
    cas=create_struct( cas , 'cq', aa )
 
    ;LONG_NAME          = 'subtile_fractions'         ,&
    ;UNITS              = '1'                         ,&
    ;SHORT_NAME         = 'FR'                        ,&
    ;DIMS               = MAPL_DimsTileTile           ,&
    ;NUM_SUBTILES       = NUM_SUBTILES                ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
readu,1,aa
    fr=aa
    cas=create_struct( cas , 'fr', aa )
 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; New Catch restart for Fortuna-2_0 does not have these next 4

    ;LONG_NAME          = 'diffuse_visible_factor_prev',&
    ;UNITS              = '1'                         ,&
    ;SHORT_NAME         = 'VISDF_PREV'                ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
;readu,1,a
;    visdf_prev=a
;    cas=create_struct( cas , 'visdf_prev', a )

    ;LONG_NAME          = 'diffuse_visible_factor_next',&
    ;UNITS              = '1'                         ,&
    ;SHORT_NAME         = 'VISDF_NEXT'                ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
;readu,1,a
;    visdf_next=a
;    cas=create_struct( cas , 'visdf_next', a )
 
    ;LONG_NAME          = 'diffuse_NIR_factor_prev'   ,&
    ;UNITS              = '1'                         ,&
    ;SHORT_NAME         = 'NIRDF_PREV'                ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
;readu,1,a
;    nirdf_prev=a
;    cas=create_struct( cas , 'nirdf_prev', a )

    ;LONG_NAME          = 'diffuse_NIR_factor_next'    ,&
    ;UNITS              = '1'                         ,&
    ;SHORT_NAME         = 'NIRDF_NEXT'                ,&
    ;DIMS               = MAPL_DimsTileOnly           ,&
    ;VLOCATION          = MAPL_VLocationNone          ,&
;readu,1,a
;    nirdf_next=a
;    cas=create_struct( cas , 'nirdf_next', a )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ;SHORT_NAME         = 'WW',                        &
    ;LONG_NAME          = 'vertical_velocity_scale_squared', &
    ;UNITS              = 'm+2 s-2',                   &
    ;DIMS               = MAPL_DimsTileTile,           &
    ;NUM_SUBTILES       = NUM_SUBTILES                ,&
    ;VLOCATION          = MAPL_VLocationNone,          &
readu,1,aa
    ww=aa
    cas=create_struct( cas , 'ww', aa )

close,1



return
end
