! Code section for reordering original array pointers and packing them into 
! 1D buffer that is then subsequently redistributed for load balancing.
! The actual pointers used for computation are then made to point back to the
! relevant porion of the 1D buffer after redistribution.
! ***CRITTICAL*** If fields are addedd or deleted to any of the ESMF states, or
! if any pointer name(s) is changed, this code section needs to modified accordingly.

! import packing and redistribution

   numUsedImp = 16  ! should match the number of imports used in this subroutine + 2 (for LATS and LONS)

!  Allocate the buffer that will hold all balanced variables. The
!   dimension of its 1D representation must ne NUMMAX---the larger of the
!   balanced and unbalanced runs.
!------------------------------------------------------------------------

      allocate(BUFIMP(NUMMAX*numUsedImp),stat=STATUS)
      VERIFY_(STATUS)
      BUFIMP = MAPL_UNDEF
      LN = 0
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'UU', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      UU =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'UWINDLMTILE', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      UWINDLMTILE =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'VWINDLMTILE', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      VWINDLMTILE =>  PTR1(1:NT)
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
      call MAPL_GetPointer(IMPORT,PTR1,'DZ', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      DZ =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'TA', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      TA =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'QA', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      QA =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'PS', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      PS =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'PCU', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      PCU =>  PTR1(1:NT)
      L1 = LN + 1
      call MAPL_GetPointer(IMPORT,PTR1,'SS_FOUND', __RC__) 
      call CICEReorder(BUFIMP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,PACKIT)
      LN = L1 + NUMMAX - 1
      PTR1(1:NUMMAX) => BUFIMP(L1:LN)
      SW =>  PTR1(1:NT)
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
         NUM_SUBTILES                        +   & ! CM
         NUM_SUBTILES                        +   & ! CQ
         NUM_SUBTILES                        +   & ! WW
         NUM_SUBTILES                        +   & ! Z0
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
      call MAPL_GetPointer(INTERNAL,PTR2,'CM', __RC__) 
      call CICEReorder(BUFINT(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),PACKIT)
      LN = L1 + NUMMAX*size(PTR2,2) - 1
      PTR2(1:NUMMAX,1:size(PTR2,2)) => BUFINT(L1:LN)
      CM =>  PTR2(1:NT,:)
      L1 = LN + 1
      call MAPL_GetPointer(INTERNAL,PTR2,'CQ', __RC__) 
      call CICEReorder(BUFINT(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),PACKIT)
      LN = L1 + NUMMAX*size(PTR2,2) - 1
      PTR2(1:NUMMAX,1:size(PTR2,2)) => BUFINT(L1:LN)
      CQ =>  PTR2(1:NT,:)
      L1 = LN + 1
      call MAPL_GetPointer(INTERNAL,PTR2,'WW', __RC__) 
      call CICEReorder(BUFINT(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),PACKIT)
      LN = L1 + NUMMAX*size(PTR2,2) - 1
      PTR2(1:NUMMAX,1:size(PTR2,2)) => BUFINT(L1:LN)
      WW =>  PTR2(1:NT,:)
      L1 = LN + 1
      call MAPL_GetPointer(INTERNAL,PTR2,'Z0', __RC__) 
      call CICEReorder(BUFINT(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),PACKIT)
      LN = L1 + NUMMAX*size(PTR2,2) - 1
      PTR2(1:NUMMAX,1:size(PTR2,2)) => BUFINT(L1:LN)
      Z0 =>  PTR2(1:NT,:)
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
      FR =>  PTR2R8(1:NT,:)
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
         28                                +   & ! 79 total single slice exports
         NUM_ICE_CATEGORIES                +   & ! QSAT1
         NUM_ICE_CATEGORIES                      ! QSAT2

      allocate(BUFEXP(NUMMAX*numExpSlices),stat=STATUS)
      VERIFY_(STATUS)

      LN = 0
      L1 = LN + 1
      call MAPL_GetPointer(EXPORT,PTR1,'QH', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         QH => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'TH', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         TH => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'UH', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         UH => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'VH', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         VH => PTR1(1:NT)
         L1 = LN + 1
      end if
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
      call MAPL_GetPointer(EXPORT,PTR1,'CHT', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         CHT => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'CMT', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         CMT => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'CQT', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         CQT => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'CNT', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         CNT => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'RIT', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         RIT => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'RET', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         RET => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'Z0', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         Z0O => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'Z0H', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         Z0H => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOT2M', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         MOT2M => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOQ2M', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         MOQ2M => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOU2M', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         MOU2M => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOV2M', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         MOV2M => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOT10M', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         MOT10M => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOQ10M', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         MOQ10M => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOU10M', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         MOU10M => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOV10M', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         MOV10M => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOU50M', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         MOU50M => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOV50M', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         MOV50M => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'GUST', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         GST => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'VENT', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         VNT => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'TFREEZE', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         TF => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FRACI', __RC__)
      if ( associated(PTR1) ) then
         LN = L1 + NUMMAX -1
         PTR1(1:NUMMAX) => BUFEXP(L1:LN)
         FRACI => PTR1(1:NT)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'QSAT1', __RC__)
      if ( associated(PTR2) ) then
         LN = L1 + NUMMAX*size(PTR2,2) -1
         PTR2(1:NUMMAX,1:size(PTR2,2)) => BUFEXP(L1:LN)
         QSAT1 => PTR2(1:NT,:)
         L1 = LN + 1
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'QSAT2', __RC__)
      if ( associated(PTR2) ) then
         LN = L1 + NUMMAX*size(PTR2,2) -1
         PTR2(1:NUMMAX,1:size(PTR2,2)) => BUFEXP(L1:LN)
         QSAT2 => PTR2(1:NT,:)
      end if
