! Pointers to inputs
!-------------------

   call MAPL_GetPointer(IMPORT,ALW    , 'ALW'    ,    _RC)
   call MAPL_GetPointer(IMPORT,BLW    , 'BLW'    ,    _RC)
   call MAPL_GetPointer(IMPORT,LWDNSRF, 'LWDNSRF',    _RC)
   call MAPL_GetPointer(IMPORT,DRPAR  , 'DRPAR'  ,    _RC)
   call MAPL_GetPointer(IMPORT,DFPAR  , 'DFPAR'  ,    _RC)
   call MAPL_GetPointer(IMPORT,DRNIR  , 'DRNIR'  ,    _RC)
   call MAPL_GetPointer(IMPORT,DFNIR  , 'DFNIR'  ,    _RC)
   call MAPL_GetPointer(IMPORT,DRUVR  , 'DRUVR'  ,    _RC)
   call MAPL_GetPointer(IMPORT,DFUVR  , 'DFUVR'  ,    _RC)
   call MAPL_GetPointer(IMPORT,EVAP   , 'EVAP'   ,    _RC)
   call MAPL_GetPointer(IMPORT,SH     , 'SH'     ,    _RC)
   call MAPL_GetPointer(IMPORT,TAUX   , 'TAUX'   ,    _RC)
   call MAPL_GetPointer(IMPORT,TAUY   , 'TAUY'   ,    _RC)
   call MAPL_GetPointer(IMPORT,DEV    , 'DEVAP'  ,    _RC)
   call MAPL_GetPointer(IMPORT,DSH    , 'DSH'    ,    _RC)
   call MAPL_GetPointer(IMPORT,SNO    , 'SNO'    ,    _RC)
   call MAPL_GetPointer(IMPORT,PLS    , 'PLS'    ,    _RC)
   call MAPL_GetPointer(IMPORT,PCU    , 'PCU'    ,    _RC)
   call MAPL_GetPointer(IMPORT,PS     , 'PS'     ,    _RC)
   call MAPL_GetPointer(IMPORT,UU     , 'UU'     ,    _RC)
   !call MAPL_GetPointer(IMPORT,TF     , 'TFREEZE',    _RC)

   ! TODO: revisit for dual_ocean
   !   call MAPL_GetPointer(IMPORT,FI     , 'FRACICE',    _RC)

   call MAPL_GetPointer(IMPORT,UW     , 'UW'     ,    _RC)
   call MAPL_GetPointer(IMPORT,VW     , 'VW'     ,    _RC)
   call MAPL_GetPointer(IMPORT,UI     , 'UI'     ,    _RC)
   call MAPL_GetPointer(IMPORT,VI     , 'VI'     ,    _RC)
   call MAPL_GetPointer(IMPORT,THATM  , 'THATM'  ,    _RC)
   call MAPL_GetPointer(IMPORT,QHATM  , 'QHATM'  ,    _RC)
   call MAPL_GetPointer(IMPORT,UHATM  , 'UHATM'  ,    _RC)
   call MAPL_GetPointer(IMPORT,VHATM  , 'VHATM'  ,    _RC)
   call MAPL_GetPointer(IMPORT,UUA    , 'UUA'    ,    _RC)
   call MAPL_GetPointer(IMPORT,VVA    , 'VVA'    ,    _RC)
   call MAPL_GetPointer(IMPORT,CTATM  , 'CTATM'  ,    _RC)
   call MAPL_GetPointer(IMPORT,CQATM  , 'CQATM'  ,    _RC)
   call MAPL_GetPointer(IMPORT,CMATM  , 'CMATM'  ,    _RC)

   call MAPL_GetPointer(IMPORT,TAUXBOT, 'TAUXBOT',    _RC)
   call MAPL_GetPointer(IMPORT,TAUYBOT, 'TAUYBOT',    _RC)
   call MAPL_GetPointer(IMPORT,TW     , 'TS_FOUND',    _RC)
   call MAPL_GetPointer(IMPORT,SW     , 'SS_FOUND',    _RC)
   call MAPL_GetPointer(IMPORT,FRZMLT , 'FRZMLT' ,    _RC)

! Pointers to internals
!----------------------

   call MAPL_GetPointer(INTERNAL,TI     ,'TSKINI',    _RC)
   call MAPL_GetPointer(INTERNAL,HI     ,'HSKINI',    _RC)
   call MAPL_GetPointer(INTERNAL,SI     ,'SSKINI',    _RC)
   call MAPL_GetPointer(INTERNAL,QS     , 'QS'   ,    _RC)
   call MAPL_GetPointer(INTERNAL,CH     , 'CH'   ,    _RC)
   call MAPL_GetPointer(INTERNAL,CQ     , 'CQ'   ,    _RC)
   call MAPL_GetPointer(INTERNAL,CM     , 'CM'   ,    _RC)

   call MAPL_GetPointer(INTERNAL,FR8    , 'FR'   ,    _RC)
   call MAPL_GetPointer(INTERNAL,VOLICE ,'VOLICE',    _RC)
   call MAPL_GetPointer(INTERNAL,VOLSNO ,'VOLSNO',    _RC)
   call MAPL_GetPointer(INTERNAL,VOLPOND,'VOLPOND',   _RC)
   call MAPL_GetPointer(INTERNAL,APONDN, 'APONDN',    _RC)
   call MAPL_GetPointer(INTERNAL,HPONDN, 'HPONDN',    _RC)
   call MAPL_GetPointer(INTERNAL,ERGICE ,'ERGICE',    _RC)
   call MAPL_GetPointer(INTERNAL,ERGSNO ,'ERGSNO',    _RC)
   call MAPL_GetPointer(INTERNAL,TAUAGE ,'TAUAGE',    _RC)
   call MAPL_GetPointer(INTERNAL,SLMASK ,'SLMASK',    _RC)

! Pointers to outputs
!--------------------

   call MAPL_GetPointer(EXPORT,EMISS  , 'EMIS'    ,    _RC)
   call MAPL_GetPointer(EXPORT,ALBVF  , 'ALBVF'   ,    _RC)
   call MAPL_GetPointer(EXPORT,ALBVR  , 'ALBVR'   ,    _RC)
   call MAPL_GetPointer(EXPORT,ALBNF  , 'ALBNF'   ,    _RC)
   call MAPL_GetPointer(EXPORT,ALBNR  , 'ALBNR'   ,    _RC)
   call MAPL_GetPointer(EXPORT,QST    , 'QST'     ,    _RC)
   call MAPL_GetPointer(EXPORT,TST    , 'TST'     ,    _RC)
   call MAPL_GetPointer(EXPORT,DELTS  , 'DELTS'   ,    _RC)
   call MAPL_GetPointer(EXPORT,DELQS  , 'DELQS'   ,    _RC)
   call MAPL_GetPointer(EXPORT,TAUXI  , 'TAUXI'   ,    _RC)
   call MAPL_GetPointer(EXPORT,TAUYI  , 'TAUYI'   ,    _RC)
   call MAPL_GetPointer(EXPORT,PENUVR , 'PENUVR'  , _RC)
   call MAPL_GetPointer(EXPORT,PENUVF , 'PENUVF'  , _RC)
   call MAPL_GetPointer(EXPORT,PENPAR , 'PENPAR'  , _RC)
   call MAPL_GetPointer(EXPORT,PENPAF , 'PENPAF'  , _RC)
   call MAPL_GetPointer(EXPORT,EVAPOUT, 'EVAPOUT' ,    _RC)
   call MAPL_GetPointer(EXPORT,SUBLIM,  'SUBLIM'  ,    _RC)
   call MAPL_GetPointer(EXPORT,SHOUT  , 'SHOUT'   ,    _RC)
   call MAPL_GetPointer(EXPORT,SHICE  , 'SHICE'   ,    _RC)
   call MAPL_GetPointer(EXPORT,HLATN  , 'HLATN'   ,    _RC)
   call MAPL_GetPointer(EXPORT,HLATICE, 'HLATICE' ,    _RC)
   call MAPL_GetPointer(EXPORT,FSURFe , 'FSURF'   ,    _RC)
   call MAPL_GetPointer(EXPORT,FSURFICE,'FSURFICE',    _RC)
   call MAPL_GetPointer(EXPORT,HLWUP  , 'HLWUP'   ,    _RC)
   call MAPL_GetPointer(EXPORT,HLWUPe , 'HLWUPICE',    _RC)
   call MAPL_GetPointer(EXPORT,LWNDSRF, 'LWNDSRF' ,    _RC)
   call MAPL_GetPointer(EXPORT,SWNDSRF, 'SWNDSRF' ,    _RC)
   call MAPL_GetPointer(EXPORT,LWNDICE, 'LWNDICE' ,    _RC)
   call MAPL_GetPointer(EXPORT,SWNDICE, 'SWNDICE' ,    _RC)
   call MAPL_GetPointer(EXPORT,FRACI  , 'FRACI'   ,    _RC)
   call MAPL_GetPointer(EXPORT,FRACINEW,'FRACINEW',    _RC)
   call MAPL_GetPointer(EXPORT,LWDNSRFe,'LWDNSRF' ,    _RC)
   call MAPL_GetPointer(EXPORT,SWDNSRFe,'SWDNSRF' ,    _RC)


   call MAPL_GetPointer(EXPORT,FRAZIL , 'FRAZIL'  ,    _RC)
   call MAPL_GetPointer(EXPORT,CONGELO, 'CONGEL'  ,    _RC)
   call MAPL_GetPointer(EXPORT,SNOICEO, 'SNOICE'  ,    _RC)
   call MAPL_GetPointer(EXPORT,FRESH  , 'FRESH'   ,    _RC)
   call MAPL_GetPointer(EXPORT,FSALT  , 'FSALT'   ,    _RC)
   call MAPL_GetPointer(EXPORT,FHOCN  , 'FHOCN'   ,    _RC)
   call MAPL_GetPointer(EXPORT,PICE   , 'PICE'    ,    _RC)
   call MAPL_GetPointer(EXPORT,FSWTRUO, 'FSWTHRU' ,    _RC)
   call MAPL_GetPointer(EXPORT,FSWABSO, 'FSWABS'  ,    _RC)
   call MAPL_GetPointer(EXPORT,MELTL  , 'MELTL'   ,    _RC)
   call MAPL_GetPointer(EXPORT,MELTTL , 'MELTT'   ,    _RC)
   call MAPL_GetPointer(EXPORT,MELTBL , 'MELTB'   ,    _RC)
   call MAPL_GetPointer(EXPORT,MELTSL , 'MELTS'   ,    _RC)
   call MAPL_GetPointer(EXPORT,HICE   , 'HICE'    ,    _RC)
   call MAPL_GetPointer(EXPORT,HSNO   , 'HSNO'    ,    _RC)
   call MAPL_GetPointer(EXPORT,HICEUNT, 'HICEUNT' ,    _RC)
   call MAPL_GetPointer(EXPORT,SNOONICE,'SNOONICE',    _RC)
   call MAPL_GetPointer(EXPORT,TSKINICE, 'TSKINICE'  ,    _RC)
   call MAPL_GetPointer(EXPORT,IAGE   , 'IAGE'    ,    _RC)
   call MAPL_GetPointer(EXPORT,DAIDTT , 'DAIDTT'  ,    _RC)
   call MAPL_GetPointer(EXPORT,DVIDTT , 'DVIDTT'  ,    _RC)
   call MAPL_GetPointer(EXPORT,FBOTL  , 'FBOT'    ,    _RC)
   call MAPL_GetPointer(EXPORT,USTARI , 'USTARI'  ,    _RC)
   call MAPL_GetPointer(EXPORT,FCONDTOP,'FCONDTOP',    _RC)
   call MAPL_GetPointer(EXPORT,FCONDB,  'FCONDBOT',    _RC)
   call MAPL_GetPointer(EXPORT,NIERG,  'NEWICEERG',    _RC)
   call MAPL_GetPointer(EXPORT,SBLXOUT,'SUBLIMFLX',    _RC)
   call MAPL_GetPointer(EXPORT,SIALB,   'SIALB'   ,    _RC)
   call MAPL_GetPointer(EXPORT,GHTSKIN, 'GHTSKIN' ,    _RC)
   call MAPL_GetPointer(EXPORT,FRZMLTe, 'FRZMLT'  ,    _RC)

   ! category dimensional exports
   call MAPL_GetPointer(EXPORT,FCONDBOTN,  'FCONDBOTN' ,  _RC)
   call MAPL_GetPointer(EXPORT,FCONDTOPN,  'FCONDTOPN' ,  _RC)
   call MAPL_GetPointer(EXPORT,TINZ     ,  'TINZ'      ,  _RC)
   call MAPL_GetPointer(EXPORT,SHICEN   ,  'SHICEN'    ,  _RC)
   call MAPL_GetPointer(EXPORT,HLWUPN   ,  'HLWUPN'    ,  _RC)
   call MAPL_GetPointer(EXPORT,LWNDSRFN ,  'LWNDSRFN'  ,  _RC)
   call MAPL_GetPointer(EXPORT,FSURFN   ,  'FSURFN'    ,  _RC)
   call MAPL_GetPointer(EXPORT,TSURFN   ,  'TSURFN'    ,  _RC)
   call MAPL_GetPointer(EXPORT,FSWSFCN  ,  'FSWSFCN'   ,  _RC)
   call MAPL_GetPointer(EXPORT,ALBINe   ,  'ALBIN'     ,  _RC)
   call MAPL_GetPointer(EXPORT,ALBSNe   ,  'ALBSN'     ,  _RC)

   ! CMIP5 exports
   call MAPL_GetPointer(EXPORT,EVAP_C5,        'evap_CMIP5' ,     _RC)
   call MAPL_GetPointer(EXPORT,PR_C5,          'pr_CMIP5'   ,     _RC)
   call MAPL_GetPointer(EXPORT,PRSN_C5,        'prsn_CMIP5' ,     _RC)
   call MAPL_GetPointer(EXPORT,GRFRAZIL_C5,    'grFrazil_CMIP5' , _RC)
   call MAPL_GetPointer(EXPORT,GRCONGEL_C5,    'grCongel_CMIP5' , _RC)
   call MAPL_GetPointer(EXPORT,GRLATERAL_C5,   'grLateral_CMIP5', _RC)
   call MAPL_GetPointer(EXPORT,SNOTOICE_C5,    'snoToIce_CMIP5' , _RC)
   call MAPL_GetPointer(EXPORT,SNOMELT_C5,     'snomelt_CMIP5' ,  _RC)
   call MAPL_GetPointer(EXPORT,TMELT_C5,       'tmelt_CMIP5' ,    _RC)
   call MAPL_GetPointer(EXPORT,BMELT_C5,       'bmelt_CMIP5' ,    _RC)
   call MAPL_GetPointer(EXPORT,SFDSI_C5,       'sfdsi_CMIP5' ,    _RC)
   call MAPL_GetPointer(EXPORT,HFSIFRAZIL_C5,  'hfsifrazil_CMIP5',_RC)
   call MAPL_GetPointer(EXPORT,IALB_C5,        'ialb_CMIP5',      _RC)
   call MAPL_GetPointer(EXPORT,RSDSSI_C5,      'rsdssi_CMIP5',    _RC)
   call MAPL_GetPointer(EXPORT,RSUSSI_C5,      'rsussi_CMIP5',    _RC)
   call MAPL_GetPointer(EXPORT,FSITHERM_CMIP5,'fsitherm_CMIP5',   _RC)