! Pointers to inputs
!-------------------

!used  call MAPL_GetPointer(IMPORT,ALW    , 'ALW'    ,    RC=STATUS); VERIFY_(STATUS)
!used  call MAPL_GetPointer(IMPORT,BLW    , 'BLW'    ,    RC=STATUS); VERIFY_(STATUS)
!used  call MAPL_GetPointer(IMPORT,LWDNSRF, 'LWDNSRF',    RC=STATUS); VERIFY_(STATUS)
!used  call MAPL_GetPointer(IMPORT,DRPAR  , 'DRPAR'  ,    RC=STATUS); VERIFY_(STATUS)
!used  call MAPL_GetPointer(IMPORT,DFPAR  , 'DFPAR'  ,    RC=STATUS); VERIFY_(STATUS)
!used  call MAPL_GetPointer(IMPORT,DRNIR  , 'DRNIR'  ,    RC=STATUS); VERIFY_(STATUS)
!used  call MAPL_GetPointer(IMPORT,DFNIR  , 'DFNIR'  ,    RC=STATUS); VERIFY_(STATUS)
!used  call MAPL_GetPointer(IMPORT,DRUVR  , 'DRUVR'  ,    RC=STATUS); VERIFY_(STATUS)
!used  call MAPL_GetPointer(IMPORT,DFUVR  , 'DFUVR'  ,    RC=STATUS); VERIFY_(STATUS)
!used  call MAPL_GetPointer(IMPORT,EVAP   , 'EVAP'   ,    RC=STATUS); VERIFY_(STATUS)
!used  call MAPL_GetPointer(IMPORT,SH     , 'SH'     ,    RC=STATUS); VERIFY_(STATUS)
!   call MAPL_GetPointer(IMPORT,TAUX   , 'TAUX'   ,    RC=STATUS); VERIFY_(STATUS)
!   call MAPL_GetPointer(IMPORT,TAUY   , 'TAUY'   ,    RC=STATUS); VERIFY_(STATUS)
!used  call MAPL_GetPointer(IMPORT,DEV    , 'DEVAP'  ,    RC=STATUS); VERIFY_(STATUS)
!used  call MAPL_GetPointer(IMPORT,DSH    , 'DSH'    ,    RC=STATUS); VERIFY_(STATUS)
!used  call MAPL_GetPointer(IMPORT,SNO    , 'SNO'    ,    RC=STATUS); VERIFY_(STATUS)
!used  call MAPL_GetPointer(IMPORT,PLS    , 'PLS'    ,    RC=STATUS); VERIFY_(STATUS)
!used  call MAPL_GetPointer(IMPORT,PCU    , 'PCU'    ,    RC=STATUS); VERIFY_(STATUS)
!used  call MAPL_GetPointer(IMPORT,PS     , 'PS'     ,    RC=STATUS); VERIFY_(STATUS)
!   call MAPL_GetPointer(IMPORT,UU     , 'UU'     ,    RC=STATUS); VERIFY_(STATUS)
   !call MAPL_GetPointer(IMPORT,TF     , 'TFREEZE',    RC=STATUS); VERIFY_(STATUS)

   ! TODO: revisit for dual_ocean
   !   call MAPL_GetPointer(IMPORT,FI     , 'FRACICE',    RC=STATUS); VERIFY_(STATUS)

!used  call MAPL_GetPointer(IMPORT,UW     , 'UW'     ,    RC=STATUS); VERIFY_(STATUS)
!used  call MAPL_GetPointer(IMPORT,VW     , 'VW'     ,    RC=STATUS); VERIFY_(STATUS)
!used  call MAPL_GetPointer(IMPORT,UI     , 'UI'     ,    RC=STATUS); VERIFY_(STATUS)
!used  call MAPL_GetPointer(IMPORT,VI     , 'VI'     ,    RC=STATUS); VERIFY_(STATUS)
!used  call MAPL_GetPointer(IMPORT,THATM  , 'THATM'  ,    RC=STATUS); VERIFY_(STATUS)
!used  call MAPL_GetPointer(IMPORT,QHATM  , 'QHATM'  ,    RC=STATUS); VERIFY_(STATUS)
!   call MAPL_GetPointer(IMPORT,UHATM  , 'UHATM'  ,    RC=STATUS); VERIFY_(STATUS)
!   call MAPL_GetPointer(IMPORT,VHATM  , 'VHATM'  ,    RC=STATUS); VERIFY_(STATUS)
!used  call MAPL_GetPointer(IMPORT,UUA    , 'UUA'    ,    RC=STATUS); VERIFY_(STATUS)
!used  call MAPL_GetPointer(IMPORT,VVA    , 'VVA'    ,    RC=STATUS); VERIFY_(STATUS)
!used  call MAPL_GetPointer(IMPORT,CTATM  , 'CTATM'  ,    RC=STATUS); VERIFY_(STATUS)
!used  call MAPL_GetPointer(IMPORT,CQATM  , 'CQATM'  ,    RC=STATUS); VERIFY_(STATUS)
!   call MAPL_GetPointer(IMPORT,CMATM  , 'CMATM'  ,    RC=STATUS); VERIFY_(STATUS)

!used  call MAPL_GetPointer(IMPORT,TAUXBOT, 'TAUXBOT',    RC=STATUS); VERIFY_(STATUS)
!used  call MAPL_GetPointer(IMPORT,TAUYBOT, 'TAUYBOT',    RC=STATUS); VERIFY_(STATUS)
!used  call MAPL_GetPointer(IMPORT,TW     , 'TS_FOUND',    RC=STATUS); VERIFY_(STATUS)
!used  call MAPL_GetPointer(IMPORT,SW     , 'SS_FOUND',    RC=STATUS); VERIFY_(STATUS)
!used  call MAPL_GetPointer(IMPORT,FRZMLT , 'FRZMLT' ,    RC=STATUS); VERIFY_(STATUS)

! Pointers to internals
!----------------------

!used   call MAPL_GetPointer(INTERNAL,TI     ,'TSKINI',    RC=STATUS); VERIFY_(STATUS) !2D
!   call MAPL_GetPointer(INTERNAL,HI     ,'HSKINI',    RC=STATUS); VERIFY_(STATUS)
!   call MAPL_GetPointer(INTERNAL,SI     ,'SSKINI',    RC=STATUS); VERIFY_(STATUS)
!used   call MAPL_GetPointer(INTERNAL,QS     , 'QS'   ,    RC=STATUS); VERIFY_(STATUS) !2D
!used   call MAPL_GetPointer(INTERNAL,CH     , 'CH'   ,    RC=STATUS); VERIFY_(STATUS) !2D
!used   call MAPL_GetPointer(INTERNAL,CQ     , 'CQ'   ,    RC=STATUS); VERIFY_(STATUS) !2D
!used   call MAPL_GetPointer(INTERNAL,CM     , 'CM'   ,    RC=STATUS); VERIFY_(STATUS) !2D

!MAPL_R8
!used   call MAPL_GetPointer(INTERNAL,FR8    , 'FR'   ,    RC=STATUS); VERIFY_(STATUS) !2D
!used   call MAPL_GetPointer(INTERNAL,VOLICE ,'VOLICE',    RC=STATUS); VERIFY_(STATUS) !2D
!used   call MAPL_GetPointer(INTERNAL,VOLSNO ,'VOLSNO',    RC=STATUS); VERIFY_(STATUS) !2D
!used   call MAPL_GetPointer(INTERNAL,VOLPOND,'VOLPOND',   RC=STATUS); VERIFY_(STATUS) !2D
!used   call MAPL_GetPointer(INTERNAL,APONDN, 'APONDN',    RC=STATUS); VERIFY_(STATUS) !2D
!used   call MAPL_GetPointer(INTERNAL,HPONDN, 'HPONDN',    RC=STATUS); VERIFY_(STATUS) !2D
!used   call MAPL_GetPointer(INTERNAL,ERGICE ,'ERGICE',    RC=STATUS); VERIFY_(STATUS) !3D
!used   call MAPL_GetPointer(INTERNAL,ERGSNO ,'ERGSNO',    RC=STATUS); VERIFY_(STATUS) !3D

!used   call MAPL_GetPointer(INTERNAL,TAUAGE ,'TAUAGE',    RC=STATUS); VERIFY_(STATUS) !2D
!used   call MAPL_GetPointer(INTERNAL,SLMASK ,'SLMASK',    RC=STATUS); VERIFY_(STATUS)

! Pointers to outputs
!--------------------

  ! call MAPL_GetPointer(EXPORT,EMISS  , 'EMIS' , alloc=.true., RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,ALBVF  , 'ALBVF', alloc=.true., RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,ALBVR  , 'ALBVR', alloc=.true., RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,ALBNF  , 'ALBNF', alloc=.true., RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,ALBNR  , 'ALBNR', alloc=.true., RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,QST    , 'QST'     ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,TST    , 'TST'     ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,DELTS  , 'DELTS'   ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,DELQS  , 'DELQS'   ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,TAUXI  , 'TAUXI'   ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,TAUYI  , 'TAUYI'   ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,PENUVR , 'PENUVR'  , RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,PENUVF , 'PENUVF'  , RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,PENPAR , 'PENPAR'  , RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,PENPAF , 'PENPAF'  , RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,EVAPOUT, 'EVAPOUT' ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,SUBLIM,  'SUBLIM'  ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,SHOUT  , 'SHOUT'   ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,SHICE  , 'SHICE'   ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,HLATN  , 'HLATN'   ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,HLATICE, 'HLATICE' ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,FSURFe , 'FSURF'   ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,FSURFICE,'FSURFICE',    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,HLWUP  , 'HLWUP'   ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,HLWUPe , 'HLWUPICE',    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,LWNDSRF, 'LWNDSRF' ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,SWNDSRF, 'SWNDSRF' ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,LWNDICE, 'LWNDICE' ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,SWNDICE, 'SWNDICE' ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,FRACI  , 'FRACI'   ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,FRACINEW,'FRACINEW',    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,LWDNSRFe,'LWDNSRF' ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,SWDNSRFe,'SWDNSRF' ,    RC=STATUS); VERIFY_(STATUS)


  ! call MAPL_GetPointer(EXPORT,FRAZIL , 'FRAZIL'  ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,CONGELO, 'CONGEL'  ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,SNOICEO, 'SNOICE'  ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,FRESH  , 'FRESH'   ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,FSALT  , 'FSALT'   ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,FHOCN  , 'FHOCN'   ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,PICE   , 'PICE'    ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,FSWTRUO, 'FSWTHRU' ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,FSWABSO, 'FSWABS'  ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,MELTL  , 'MELTL'   ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,MELTTL , 'MELTT'   ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,MELTBL , 'MELTB'   ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,MELTSL , 'MELTS'   ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,HICE   , 'HICE'    ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,HSNO   , 'HSNO'    ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,HICEUNT, 'HICEUNT' ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,SNOONICE,'SNOONICE',    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,ISTSFC , 'ISTSFC'  ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,IAGE   , 'IAGE'    ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,DAIDTT , 'DAIDTT'  ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,DVIDTT , 'DVIDTT'  ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,FBOTL  , 'FBOT'    ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,USTARI , 'USTARI'  ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,FCONDTOP,'FCONDTOP',    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,FCONDB,  'FCONDBOT',    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,NIERG,  'NEWICEERG',    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,SBLXOUT,'SUBLIMFLX',    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,SIALB,   'SIALB'   ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,GHTSKIN, 'GHTSKIN' ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,FRZMLTe, 'FRZMLT'  ,    RC=STATUS); VERIFY_(STATUS)

   ! category dimensional exports
  ! call MAPL_GetPointer(EXPORT,FCONDBOTN,  'FCONDBOTN' ,  RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,FCONDTOPN,  'FCONDTOPN' ,  RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,TINZ     ,  'TINZ'      ,  RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,SHICEN   ,  'SHICEN'    ,  RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,HLWUPN   ,  'HLWUPN'    ,  RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,LWNDSRFN ,  'LWNDSRFN'  ,  RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,FSURFN   ,  'FSURFN'    ,  RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,TSURFN   ,  'TSURFN'    ,  RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,FSWSFCN  ,  'FSWSFCN'   ,  RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,ALBINe   ,  'ALBIN'     ,  RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,ALBSNe   ,  'ALBSN'     ,  RC=STATUS); VERIFY_(STATUS)

   ! CMIP5 exports
  ! call MAPL_GetPointer(EXPORT,EVAP_C5,        'evap_CMIP5' ,     RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,PR_C5,          'pr_CMIP5'   ,     RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,PRSN_C5,        'prsn_CMIP5' ,     RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,GRFRAZIL_C5,    'grFrazil_CMIP5' , RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,GRCONGEL_C5,    'grCongel_CMIP5' , RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,GRLATERAL_C5,   'grLateral_CMIP5', RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,SNOTOICE_C5,    'snoToIce_CMIP5' , RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,SNOMELT_C5,     'snomelt_CMIP5' ,  RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,TMELT_C5,       'tmelt_CMIP5' ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,BMELT_C5,       'bmelt_CMIP5' ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,SFDSI_C5,       'sfdsi_CMIP5' ,    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,HFSIFRAZIL_C5,  'hfsifrazil_CMIP5',RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,IALB_C5,        'ialb_CMIP5',      RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,RSDSSI_C5,      'rsdssi_CMIP5',    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,RSUSSI_C5,      'rsussi_CMIP5',    RC=STATUS); VERIFY_(STATUS)
  ! call MAPL_GetPointer(EXPORT,FSITHERM_CMIP5,'fsitherm_CMIP5',   RC=STATUS); VERIFY_(STATUS)
