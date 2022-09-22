! Pointers to inputs
!-------------------

   call MAPL_GetPointer(IMPORT,UU     , 'UU'     ,    _RC)
   call MAPL_GetPointer(IMPORT,UWINDLMTILE     , 'UWINDLMTILE'     ,    _RC)
   call MAPL_GetPointer(IMPORT,VWINDLMTILE     , 'VWINDLMTILE'     ,    _RC)
   call MAPL_GetPointer(IMPORT,UI     , 'UI'     ,    _RC)
   call MAPL_GetPointer(IMPORT,VI     , 'VI'     ,    _RC)
   call MAPL_GetPointer(IMPORT,DZ     , 'DZ'     ,    _RC)
   call MAPL_GetPointer(IMPORT,TA     , 'TA'     ,    _RC)
   call MAPL_GetPointer(IMPORT,QA     , 'QA'     ,    _RC)
   call MAPL_GetPointer(IMPORT,PS     , 'PS'     ,    _RC)
   call MAPL_GetPointer(IMPORT,PCU    , 'PCU'    ,    _RC)
   call MAPL_GetPointer(IMPORT,SW     , 'SS_FOUND' ,  _RC)
   ! the call below may be needed in dual-ocean mode 
   !   call MAPL_GetPointer(IMPORT,FI     , 'FRACICE',    _RC)


! Pointers to internals
!----------------------

   call MAPL_GetPointer(INTERNAL,TI   , 'TSKINI' ,    _RC)
   call MAPL_GetPointer(INTERNAL,QS   , 'QS'     ,    _RC)
   call MAPL_GetPointer(INTERNAL,CH   , 'CH'     ,    _RC)
   call MAPL_GetPointer(INTERNAL,CM   , 'CM'     ,    _RC)
   call MAPL_GetPointer(INTERNAL,CQ   , 'CQ'     ,    _RC)
   call MAPL_GetPointer(INTERNAL,Z0   , 'Z0'     ,    _RC)
   call MAPL_GetPointer(INTERNAL,WW   , 'WW'     ,    _RC)
   call MAPL_GetPointer(INTERNAL,FR   , 'FR'     ,    _RC) 
   call MAPL_GetPointer(INTERNAL,VOLICE ,'VOLICE',    _RC) 
   call MAPL_GetPointer(INTERNAL,VOLSNO ,'VOLSNO',    _RC) 
   call MAPL_GetPointer(INTERNAL,VOLPOND,'VOLPOND',   _RC) 
   call MAPL_GetPointer(INTERNAL,ERGICE ,'ERGICE',    _RC)
   call MAPL_GetPointer(INTERNAL,ERGSNO ,'ERGSNO',    _RC) 
   call MAPL_GetPointer(INTERNAL,TAUAGE ,'TAUAGE',    _RC) 
   call MAPL_GetPointer(INTERNAL,SLMASK ,'SLMASK',    _RC) 

! Pointers to outputs
!--------------------

   call MAPL_GetPointer(EXPORT,QH    , 'QH'      ,    _RC)
   call MAPL_GetPointer(EXPORT,TH    , 'TH'      ,    _RC)
   call MAPL_GetPointer(EXPORT,UH    , 'UH'      ,    _RC)
   call MAPL_GetPointer(EXPORT,VH    , 'VH'      ,    _RC)
   call MAPL_GetPointer(EXPORT,QST   , 'QST'     ,    _RC)
   call MAPL_GetPointer(EXPORT,TST   , 'TST'     ,    _RC)
   call MAPL_GetPointer(EXPORT,CHT   , 'CHT'     ,    _RC)
   call MAPL_GetPointer(EXPORT,CMT   , 'CMT'     ,    _RC)
   call MAPL_GetPointer(EXPORT,CQT   , 'CQT'     ,    _RC)
   call MAPL_GetPointer(EXPORT,CNT   , 'CNT'     ,    _RC)
   call MAPL_GetPointer(EXPORT,RIT   , 'RIT'     ,    _RC)
   call MAPL_GetPointer(EXPORT,RET   , 'RET'     ,    _RC)
   call MAPL_GetPointer(EXPORT,Z0O   , 'Z0'      ,    _RC)
   call MAPL_GetPointer(EXPORT,Z0H   , 'Z0H'     ,    _RC)
   call MAPL_GetPointer(EXPORT,MOT2M, 'MOT2M'   ,    _RC)
   call MAPL_GetPointer(EXPORT,MOQ2M, 'MOQ2M'   ,    _RC)
   call MAPL_GetPointer(EXPORT,MOU2M, 'MOU2M'  ,    _RC)
   call MAPL_GetPointer(EXPORT,MOV2M, 'MOV2M'  ,    _RC)
   call MAPL_GetPointer(EXPORT,MOT10M, 'MOT10M'   ,    _RC)
   call MAPL_GetPointer(EXPORT,MOQ10M, 'MOQ10M'   ,    _RC)
   call MAPL_GetPointer(EXPORT,MOU10M, 'MOU10M'  ,    _RC)
   call MAPL_GetPointer(EXPORT,MOV10M, 'MOV10M'  ,    _RC)
   call MAPL_GetPointer(EXPORT,MOU50M, 'MOU50M'  ,    _RC)
   call MAPL_GetPointer(EXPORT,MOV50M, 'MOV50M'  ,    _RC)
   call MAPL_GetPointer(EXPORT,GST   , 'GUST'    ,    _RC)
   call MAPL_GetPointer(EXPORT,VNT   , 'VENT'    ,    _RC)
   call MAPL_GetPointer(EXPORT,QSAT1 , 'QSAT1'   ,    _RC)
   call MAPL_GetPointer(EXPORT,QSAT2 , 'QSAT2'   ,    _RC)

  ! export to openwater
   call MAPL_GetPointer(EXPORT,TF    , 'TFREEZE' ,    _RC)
   call MAPL_GetPointer(EXPORT,FRACI , 'FRACI'   ,    _RC)
