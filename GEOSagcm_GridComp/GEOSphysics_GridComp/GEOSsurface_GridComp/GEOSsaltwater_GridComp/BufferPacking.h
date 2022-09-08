! Code section for reordering original array pointers and packing them into 
! 1D buffer that is then subsequently redistributed for load balancing.
! The actual pointers used for computation are then made to point back to the
! relevant porion of the 1D buffer after redistribution.
! ***CRITTICAL*** If fields are addedd or deleted to any of the ESMF states, or
! if any pointer name(s) is changed, this code section needs to modified accordingly.

! import packing and redistribution

   numUsedImp = 34  ! should match the number of imports used in this subroutine + 2 (for LATS and LONS)

!  Allocate the buffer that will hold all balanced variables. The
!   dimension of its 1D representation must ne NUMMAX---the larger of the
!   balanced and unbalanced runs.
!------------------------------------------------------------------------

      allocate(BUFIMP(NUMMAX*numUsedImp),stat=STATUS)
      VERIFY_(STATUS)
      BUFIMP = MAPL_UNDEF
      LN = 0
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'ALW', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      ALW =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'BLW', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      BLW =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'LWDNSRF', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      LWDNSRF =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'DRPAR', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      DRPAR =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'DFPAR', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      DFPAR =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'DRNIR', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      DRNIR =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'DFNIR', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      DFNIR =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'DRUVR', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      DRUVR =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'DFUVR', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      DFUVR =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'EVAP', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      EVAP =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'SH', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      SH =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'DEVAP', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      DEV =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'DSH', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      DSH =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'SNO', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      SNO =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'PLS', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      PLS =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'PCU', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      PCU =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'PS', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      PS =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'UW', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      UW =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'VW', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      VW =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'UI', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      UI =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'VI', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      VI =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'THATM', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      THATM =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'QHATM', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      QHATM =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'UUA', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      UUA =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'VVA', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      VVA =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'CTATM', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      CTATM =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'CQATM', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      CQATM =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'TAUXBOT', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      TAUXBOT =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'TAUYBOT', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      TAUYBOT =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'TS_FOUND', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      TW =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'SS_FOUND', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      SW =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'FRZMLT', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      FRZMLT =>  PTR1(1:NT)
      L1 = LN + 1
      call CICEReorder(BUFIMP(L1),LATS_ORIGINAL,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      LATS =>  PTR1(1:NT) 
      L1 = LN + 1
      call CICEReorder(BUFIMP(L1),LONS_ORIGINAL,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      LONS =>  PTR1(1:NT)

      call MAPL_BalanceWork(BUFIMP, NUMMAX, Direction=MAPL_Distribute, Handle=CICECOREBalanceHandle, __RC__)

! REAL4 internal  packing and redistribution
      numIntSlices =                             &
         NUM_ICE_CATEGORIES                  +   & ! TSKINI
         NUM_SUBTILES                        +   & ! QS
         NUM_SUBTILES                        +   & ! CH
         NUM_SUBTILES                        +   & ! CQ
         NUM_SUBTILES                        +   & ! CM
         NUM_ICE_CATEGORIES                  +   & ! TAUAGE
         1                                         ! SLMASK
         
      allocate(BUFINT(NUMMAX*numIntSlices),stat=STATUS)
      VERIFY_(STATUS)
      BUFINT = MAPL_UNDEF
      LN = 0
      L1 = LN + 1
      call MAPL_GetPointer(INTERNAL,PTR2,'TSKINI', __RC__) 
      call CICEReorder(BUFINT(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),PACKIT)
      LN = L1 + NUMMAX*size(PTR2,2) - 1
      PTR2(1:NUMMAX,1:size(PTR2,2)) => BUFINT(L1:LN)
      TI =>  PTR2(1:NT,:)
      L1 = LN + 1
      call MAPL_GetPointer(INTERNAL,PTR2,'QS', __RC__) 
      call CICEReorder(BUFINT(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),PACKIT)
      LN = L1 + NUMMAX*size(PTR2,2) - 1
      PTR2(1:NUMMAX,1:size(PTR2,2)) => BUFINT(L1:LN)
      QS =>  PTR2(1:NT,:)
      L1 = LN + 1
      call MAPL_GetPointer(INTERNAL,PTR2,'CH', __RC__) 
      call CICEReorder(BUFINT(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),PACKIT)
      LN = L1 + NUMMAX*size(PTR2,2) - 1
      PTR2(1:NUMMAX,1:size(PTR2,2)) => BUFINT(L1:LN)
      CH =>  PTR2(1:NT,:)
      L1 = LN + 1
      call MAPL_GetPointer(INTERNAL,PTR2,'CQ', __RC__) 
      call CICEReorder(BUFINT(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),PACKIT)
      LN = L1 + NUMMAX*size(PTR2,2) - 1
      PTR2(1:NUMMAX,1:size(PTR2,2)) => BUFINT(L1:LN)
      CQ =>  PTR2(1:NT,:)
      L1 = LN + 1
      call MAPL_GetPointer(INTERNAL,PTR2,'CM', __RC__) 
      call CICEReorder(BUFINT(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),PACKIT)
      LN = L1 + NUMMAX*size(PTR2,2) - 1
      PTR2(1:NUMMAX,1:size(PTR2,2)) => BUFINT(L1:LN)
      CM =>  PTR2(1:NT,:)
      L1 = LN + 1
      call MAPL_GetPointer(INTERNAL,PTR2,'TAUAGE', __RC__) 
      call CICEReorder(BUFINT(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),PACKIT)
      LN = L1 + NUMMAX*size(PTR2,2) - 1
      PTR2(1:NUMMAX,1:size(PTR2,2)) => BUFINT(L1:LN)
      TAUAGE =>  PTR2(1:NT,:)
      L1 = LN + 1
      call MAPL_GetPointer(INTERNAL,PTR1,'SLMASK', __RC__) 
      call CICEReorder(BUFINT(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFINT(L1:LN)
      SLMASK =>  PTR1(1:NT)

      call MAPL_BalanceWork(BUFINT, NUMMAX, Direction=MAPL_Distribute, Handle=CICECOREBalanceHandle, __RC__)

! REAL8 internal  packing and redistribution
      numIntSlices8 =                            &
         NUM_SUBTILES                        +   & ! FR
         NUM_ICE_CATEGORIES                  +   & ! VOLICE
         NUM_ICE_CATEGORIES                  +   & ! VOLSNO
         NUM_ICE_CATEGORIES                  +   & ! VOLPOND
         NUM_ICE_CATEGORIES                  +   & ! APONDN
         NUM_ICE_CATEGORIES                  +   & ! HPONDN
         NUM_ICE_LAYERS*NUM_ICE_CATEGORIES   +   & ! ERGICE
         NUM_SNOW_LAYERS*NUM_ICE_CATEGORIES        ! ERGSNO

      allocate(BUFINT8(NUMMAX*numIntSlices8),stat=STATUS)
      VERIFY_(STATUS)
      BUFINT8 = MAPL_UNDEF
      LN = 0
      L1 = LN + 1
      call MAPL_GetPointer(INTERNAL,PTR2R8,'FR', __RC__) 
      call CICEReorder8(BUFINT8(L1),PTR2R8,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2R8,2),PACKIT)
      LN = L1 + NUMMAX*size(PTR2R8,2) - 1
      PTR2R8(1:NUMMAX,1:size(PTR2R8,2)) => BUFINT8(L1:LN)
      FR8 =>  PTR2R8(1:NT,:)
      L1 = LN + 1
      call MAPL_GetPointer(INTERNAL,PTR2R8,'VOLICE', __RC__) 
      call CICEReorder8(BUFINT8(L1),PTR2R8,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2R8,2),PACKIT)
      LN = L1 + NUMMAX*size(PTR2R8,2) - 1
      PTR2R8(1:NUMMAX,1:size(PTR2R8,2)) => BUFINT8(L1:LN)
      VOLICE =>  PTR2R8(1:NT,:)
      L1 = LN + 1
      call MAPL_GetPointer(INTERNAL,PTR2R8,'VOLSNO', __RC__) 
      call CICEReorder8(BUFINT8(L1),PTR2R8,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2R8,2),PACKIT)
      LN = L1 + NUMMAX*size(PTR2R8,2) - 1
      PTR2R8(1:NUMMAX,1:size(PTR2R8,2)) => BUFINT8(L1:LN)
      VOLSNO =>  PTR2R8(1:NT,:)
      L1 = LN + 1
      call MAPL_GetPointer(INTERNAL,PTR2R8,'VOLPOND', __RC__) 
      call CICEReorder8(BUFINT8(L1),PTR2R8,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2R8,2),PACKIT)
      LN = L1 + NUMMAX*size(PTR2R8,2) - 1
      PTR2R8(1:NUMMAX,1:size(PTR2R8,2)) => BUFINT8(L1:LN)
      VOLPOND =>  PTR2R8(1:NT,:)
      L1 = LN + 1
      call MAPL_GetPointer(INTERNAL,PTR2R8,'APONDN', __RC__) 
      call CICEReorder8(BUFINT8(L1),PTR2R8,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2R8,2),PACKIT)
      LN = L1 + NUMMAX*size(PTR2R8,2) - 1
      PTR2R8(1:NUMMAX,1:size(PTR2R8,2)) => BUFINT8(L1:LN)
      APONDN =>  PTR2R8(1:NT,:)
      L1 = LN + 1
      call MAPL_GetPointer(INTERNAL,PTR2R8,'HPONDN', __RC__) 
      call CICEReorder8(BUFINT8(L1),PTR2R8,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2R8,2),PACKIT)
      LN = L1 + NUMMAX*size(PTR2R8,2) - 1
      PTR2R8(1:NUMMAX,1:size(PTR2R8,2)) => BUFINT8(L1:LN)
      HPONDN =>  PTR2R8(1:NT,:)
      L1 = LN + 1
      call MAPL_GetPointer(INTERNAL,PTR3R8,'ERGICE', __RC__) 
      call CICEReorder8(BUFINT8(L1),PTR3R8,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR3R8,2)*size(PTR3R8,3),PACKIT)
      LN = L1 + NUMMAX*size(PTR3R8,2)*size(PTR3R8,3) - 1
      PTR3R8(1:NUMMAX,1:size(PTR3R8,2),1:size(PTR3R8,3)) => BUFINT8(L1:LN)
      ERGICE =>  PTR3R8(1:NT,:,:)
      L1 = LN + 1
      call MAPL_GetPointer(INTERNAL,PTR3R8,'ERGSNO', __RC__) 
      call CICEReorder8(BUFINT8(L1),PTR3R8,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR3R8,2)*size(PTR3R8,3),PACKIT)
      LN = L1 + NUMMAX*size(PTR3R8,2)*size(PTR3R8,3) - 1
      PTR3R8(1:NUMMAX,1:size(PTR3R8,2),1:size(PTR3R8,3)) => BUFINT8(L1:LN)
      ERGSNO =>  PTR3R8(1:NT,:,:)

      call MAPL_BalanceWork(BUFINT8, NUMMAX, Direction=MAPL_Distribute, Handle=CICECOREBalanceHandle, __RC__)


! export packing and redistribution
! not all exports will necessarily be needed but memory will be allocated just in case
      numExpSlices =                           &
         79                                +   & ! 79 total single slice exports
         NUM_ICE_CATEGORIES                +   & ! FCONDBOTN
         NUM_ICE_CATEGORIES                +   & ! FCONDTOPN
         NUM_ICE_LAYERS*NUM_ICE_CATEGORIES +   & ! TINZ
         NUM_ICE_CATEGORIES                +   & ! SHICEN
         NUM_ICE_CATEGORIES                +   & ! HLWUPN
         NUM_ICE_CATEGORIES                +   & ! LWNDSRFN
         NUM_ICE_CATEGORIES                +   & ! FSURFN
         NUM_ICE_CATEGORIES                +   & ! TSURFN
         NUM_ICE_CATEGORIES                +   & ! FSWSFCN
         NUM_ICE_CATEGORIES                +   & ! ALBIN
         NUM_ICE_CATEGORIES                      ! ALBSN

      allocate(BUFEXP(NUMMAX*numExpSlices),stat=STATUS)
      VERIFY_(STATUS)

      LN = 0
      L1 = LN + 1
      call MAPL_GetPointer(EXPORT,PTR1,'EMIS', __RC__)
      LN = L1 + NUMMAX -1
      PTR1(1:NUMMAX) => BUFEXP(L1:LN)
      EMISS => PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(EXPORT,PTR1,'ALBVF', __RC__)
      LN = L1 + NUMMAX -1
      PTR1(1:NUMMAX) => BUFEXP(L1:LN)
      ALBVF => PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(EXPORT,PTR1,'ALBVR', __RC__)
      LN = L1 + NUMMAX -1
      PTR1(1:NUMMAX) => BUFEXP(L1:LN)
      ALBVR => PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(EXPORT,PTR1,'ALBNF', __RC__)
      LN = L1 + NUMMAX -1
      PTR1(1:NUMMAX) => BUFEXP(L1:LN)
      ALBNF => PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(EXPORT,PTR1,'ALBNR', __RC__)
      LN = L1 + NUMMAX -1
      PTR1(1:NUMMAX) => BUFEXP(L1:LN)
      ALBNR => PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(EXPORT,PTR1,'QST', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         QST => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'TST', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         TST => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'DELTS', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         DELTS => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'DELQS', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         DELQS => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'TAUXI', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         TAUXI => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'TAUYI', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         TAUYI => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'PENUVR', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         PENUVR => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'PENUVF', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         PENUVF => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'PENPAR', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         PENPAR => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'PENPAF', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         PENPAF => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'EVAPOUT', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         EVAPOUT => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SUBLIM', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         SUBLIM => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SHOUT', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         SHOUT => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SHICE', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         SHICE => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'HLATN', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         HLATN => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'HLATICE', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         HLATICE => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FSURF', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         FSURFe => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FSURFICE', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         FSURFICE => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'HLWUP', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         HLWUP => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'HLWUPICE', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         HLWUPe => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'LWNDSRF', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         LWNDSRF => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SWNDSRF', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         SWNDSRF => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'LWNDICE', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         LWNDICE => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SWNDICE', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         SWNDICE => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FRACINEW', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         FRACINEW => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'LWDNSRF', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         LWDNSRFe => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SWDNSRF', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         SWDNSRFe => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FRAZIL', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         FRAZIL => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'CONGEL', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         CONGELO => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SNOICE', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         SNOICEO => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FRESH', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         FRESH => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FSALT', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         FSALT => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FHOCN', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         FHOCN => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'PICE', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         PICE => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FSWTHRU', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         FSWTRUO => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FSWABS', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         FSWABSO => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MELTL', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         MELTL => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MELTT', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         MELTTL => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MELTB', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         MELTBL => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MELTS', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         MELTSL => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'HICE', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         HICE => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'HSNO', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         HSNO => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'HICEUNT', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         HICEUNT => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SNOONICE', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         SNOONICE => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'TSKINICE', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         TSKINICE => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'IAGE', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         IAGE => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'DAIDTT', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         DAIDTT => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'DVIDTT', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         DVIDTT => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FBOT', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         FBOTL => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'USTARI', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         USTARI => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FCONDTOP', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         FCONDTOP => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FCONDBOT', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         FCONDB => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'NEWICEERG', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         NIERG => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SUBLIMFLX', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         SBLXOUT => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SIALB', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         SIALB => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'GHTSKIN', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         GHTSKIN => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FRZMLT', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         FRZMLTe => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'evap_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         EVAP_C5 => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'pr_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         PR_C5 => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'prsn_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         PRSN_C5 => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'grFrazil_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         GRFRAZIL_C5 => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'grCongel_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         GRCONGEL_C5 => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'grLateral_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         GRLATERAL_C5 => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'snoToIce_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         SNOTOICE_C5 => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'snomelt_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         SNOMELT_C5 => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'tmelt_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         TMELT_C5 => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'bmelt_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         BMELT_C5 => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'sfdsi_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         SFDSI_C5 => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'hfsifrazil_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         HFSIFRAZIL_C5 => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'ialb_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         IALB_C5 => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'rsdssi_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         RSDSSI_C5 => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'rsussi_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         RSUSSI_C5 => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'fsitherm_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         FSITHERM_CMIP5 => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'FCONDBOTN', __RC__)
      if ( associated(PTR2) ) then
         LN = L1 + NUMMAX*size(PTR2,2) -1
         PTR2(1:NUMMAX,1:size(PTR2,2)) => BUFEXP(L1:LN)
         FCONDBOTN => PTR2(1:NT,:)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'FCONDTOPN', __RC__)
      if ( associated(PTR2) ) then
         LN = L1 + NUMMAX*size(PTR2,2) -1
         PTR2(1:NUMMAX,1:size(PTR2,2)) => BUFEXP(L1:LN)
         FCONDTOPN => PTR2(1:NT,:)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'TINZ', __RC__)
      if ( associated(PTR2) ) then
         LN = L1 + NUMMAX*size(PTR2,2) -1
         PTR2(1:NUMMAX,1:size(PTR2,2)) => BUFEXP(L1:LN)
         TINZ => PTR2(1:NT,:)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'SHICEN', __RC__)
      if ( associated(PTR2) ) then
         LN = L1 + NUMMAX*size(PTR2,2) -1
         PTR2(1:NUMMAX,1:size(PTR2,2)) => BUFEXP(L1:LN)
         SHICEN => PTR2(1:NT,:)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'HLWUPN', __RC__)
      if ( associated(PTR2) ) then
         LN = L1 + NUMMAX*size(PTR2,2) -1
         PTR2(1:NUMMAX,1:size(PTR2,2)) => BUFEXP(L1:LN)
         HLWUPN => PTR2(1:NT,:)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'LWNDSRFN', __RC__)
      if ( associated(PTR2) ) then
         LN = L1 + NUMMAX*size(PTR2,2) -1
         PTR2(1:NUMMAX,1:size(PTR2,2)) => BUFEXP(L1:LN)
         LWNDSRFN => PTR2(1:NT,:)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'FSURFN', __RC__)
      if ( associated(PTR2) ) then
         LN = L1 + NUMMAX*size(PTR2,2) -1
         PTR2(1:NUMMAX,1:size(PTR2,2)) => BUFEXP(L1:LN)
         FSURFN => PTR2(1:NT,:)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'TSURFN', __RC__)
      if ( associated(PTR2) ) then
         LN = L1 + NUMMAX*size(PTR2,2) -1
         PTR2(1:NUMMAX,1:size(PTR2,2)) => BUFEXP(L1:LN)
         TSURFN => PTR2(1:NT,:)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'FSWSFCN', __RC__)
      if ( associated(PTR2) ) then
         LN = L1 + NUMMAX*size(PTR2,2) -1
         PTR2(1:NUMMAX,1:size(PTR2,2)) => BUFEXP(L1:LN)
         FSWSFCN => PTR2(1:NT,:)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'ALBIN', __RC__)
      if ( associated(PTR2) ) then
         LN = L1 + NUMMAX*size(PTR2,2) -1
         PTR2(1:NUMMAX,1:size(PTR2,2)) => BUFEXP(L1:LN)
         ALBINe => PTR2(1:NT,:)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'ALBSN', __RC__)
      if ( associated(PTR2) ) then
         LN = L1 + NUMMAX*size(PTR2,2) -1
         PTR2(1:NUMMAX,1:size(PTR2,2)) => BUFEXP(L1:LN)
         ALBSNe => PTR2(1:NT,:)
      end if
