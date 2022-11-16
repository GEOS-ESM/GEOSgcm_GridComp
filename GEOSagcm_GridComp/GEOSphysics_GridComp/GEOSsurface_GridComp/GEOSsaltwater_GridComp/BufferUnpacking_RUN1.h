! This section reverses what was done in BufferPacking by redistributing the 
! the updated 1D buffer back to the original distribution and then map relevant
! portions to the original pointers in the INTERNAL and EXPORT fields. 
! There is no need to do IMPORT on exit.
! ***CRITTICAL*** For any changes (pointer variable name change, field addition or 
! deletion) made in BufferPacking, corresponding changes need to be made here 
! as well.

! retrieve real4 internal
      call MAPL_BalanceWork(BUFINT, NUMMAX, Direction=MAPL_Retrieve, Handle=CICECOREBalanceHandle, _RC)
      L1 = 1
      call MAPL_GetPointer(INTERNAL,PTR2,'TSKINI', _RC) 
      call CICEReorder(BUFINT(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR2,2)
      call MAPL_GetPointer(INTERNAL,PTR2,'QS', _RC) 
      call CICEReorder(BUFINT(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR2,2)
      call MAPL_GetPointer(INTERNAL,PTR2,'CH', _RC) 
      call CICEReorder(BUFINT(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR2,2)
      call MAPL_GetPointer(INTERNAL,PTR2,'CM', _RC) 
      call CICEReorder(BUFINT(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR2,2)
      call MAPL_GetPointer(INTERNAL,PTR2,'CQ', _RC) 
      call CICEReorder(BUFINT(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR2,2)
      call MAPL_GetPointer(INTERNAL,PTR2,'WW', _RC) 
      call CICEReorder(BUFINT(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR2,2)
      call MAPL_GetPointer(INTERNAL,PTR2,'Z0', _RC) 
      call CICEReorder(BUFINT(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR2,2)
      call MAPL_GetPointer(INTERNAL,PTR2,'TAUAGE', _RC) 
      call CICEReorder(BUFINT(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR2,2)
      call MAPL_GetPointer(INTERNAL,PTR1,'SLMASK', _RC) 
      call CICEReorder(BUFINT(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)

! retrieve real8 internal
      call MAPL_BalanceWork(BUFINT8, NUMMAX, Direction=MAPL_Retrieve, Handle=CICECOREBalanceHandle, _RC)
      L1 = 1
      call MAPL_GetPointer(INTERNAL,PTR2R8,'FR', _RC) 
      call CICEReorder8(BUFINT8(L1),PTR2R8,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2R8,2),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR2R8,2)
      call MAPL_GetPointer(INTERNAL,PTR2R8,'VOLICE', _RC) 
      call CICEReorder8(BUFINT8(L1),PTR2R8,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2R8,2),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR2R8,2)
      call MAPL_GetPointer(INTERNAL,PTR2R8,'VOLSNO', _RC) 
      call CICEReorder8(BUFINT8(L1),PTR2R8,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2R8,2),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR2R8,2)
      call MAPL_GetPointer(INTERNAL,PTR2R8,'VOLPOND', _RC) 
      call CICEReorder8(BUFINT8(L1),PTR2R8,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2R8,2),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR2R8,2)
      call MAPL_GetPointer(INTERNAL,PTR3R8,'ERGICE', _RC) 
      call CICEReorder8(BUFINT8(L1),PTR3R8,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR3R8,2)*size(PTR3R8,3),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR3R8,2)*size(PTR3R8,3)
      call MAPL_GetPointer(INTERNAL,PTR3R8,'ERGSNO', _RC) 
      call CICEReorder8(BUFINT8(L1),PTR3R8,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR3R8,2)*size(PTR3R8,3),UNPACKIT)

! retrieve export
      call MAPL_BalanceWork(BUFEXP, NUMMAX, Direction=MAPL_Retrieve, Handle=CICECOREBalanceHandle, _RC)
      L1 = 1
      call MAPL_GetPointer(EXPORT,PTR1,'QH', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'TH', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'UH', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'VH', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'QST', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'TST', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'CHT', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'CMT', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'CQT', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'CNT', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'RIT', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'RET', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'Z0', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'Z0H', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOT2M', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOQ2M', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOU2M', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOV2M', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOT10M', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOQ10M', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOU10M', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOV10M', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOU50M', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOV50M', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'GUST', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'VENT', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'TFREEZE', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FRACI', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'QSAT1', _RC)
      if ( associated(PTR2) ) then
         call CICEReorder(BUFEXP(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
         L1 = L1 + NUMMAX*size(PTR2,2)
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'QSAT2', _RC)
      if ( associated(PTR2) ) then
         call CICEReorder(BUFEXP(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
      end if
