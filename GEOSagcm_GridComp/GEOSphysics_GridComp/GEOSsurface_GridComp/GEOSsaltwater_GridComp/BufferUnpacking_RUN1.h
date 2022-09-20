! This section reverses what was done in BufferPacking by redistributing the 
! the updated 1D buffer back to the original distribution and then map relevant
! portions to the original pointers in the INTERNAL and EXPORT fields. 
! There is no need to do IMPORT on exit.
! ***CRITTICAL*** For any changes (pointer variable name change, field addition or 
! deletion) made in BufferPacking, corresponding changes need to be made here 
! as well.

! retrieve real4 internal
      call MAPL_BalanceWork(BUFINT, NUMMAX, Direction=MAPL_Retrieve, Handle=CICECOREBalanceHandle, __RC__)
      L1 = 1
      call MAPL_GetPointer(INTERNAL,PTR2,'TSKINI', __RC__) 
      call CICEReorder(BUFINT(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR2,2)
      call MAPL_GetPointer(INTERNAL,PTR2,'QS', __RC__) 
      call CICEReorder(BUFINT(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR2,2)
      call MAPL_GetPointer(INTERNAL,PTR2,'CH', __RC__) 
      call CICEReorder(BUFINT(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR2,2)
      call MAPL_GetPointer(INTERNAL,PTR2,'CM', __RC__) 
      call CICEReorder(BUFINT(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR2,2)
      call MAPL_GetPointer(INTERNAL,PTR2,'CQ', __RC__) 
      call CICEReorder(BUFINT(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR2,2)
      call MAPL_GetPointer(INTERNAL,PTR2,'WW', __RC__) 
      call CICEReorder(BUFINT(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR2,2)
      call MAPL_GetPointer(INTERNAL,PTR2,'Z0', __RC__) 
      call CICEReorder(BUFINT(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR2,2)
      call MAPL_GetPointer(INTERNAL,PTR2,'TAUAGE', __RC__) 
      call CICEReorder(BUFINT(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR2,2)
      call MAPL_GetPointer(INTERNAL,PTR1,'SLMASK', __RC__) 
      call CICEReorder(BUFINT(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)

! retrieve real8 internal
      call MAPL_BalanceWork(BUFINT8, NUMMAX, Direction=MAPL_Retrieve, Handle=CICECOREBalanceHandle, __RC__)
      L1 = 1
      call MAPL_GetPointer(INTERNAL,PTR2R8,'FR', __RC__) 
      call CICEReorder8(BUFINT8(L1),PTR2R8,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2R8,2),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR2R8,2)
      call MAPL_GetPointer(INTERNAL,PTR2R8,'VOLICE', __RC__) 
      call CICEReorder8(BUFINT8(L1),PTR2R8,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2R8,2),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR2R8,2)
      call MAPL_GetPointer(INTERNAL,PTR2R8,'VOLSNO', __RC__) 
      call CICEReorder8(BUFINT8(L1),PTR2R8,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2R8,2),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR2R8,2)
      call MAPL_GetPointer(INTERNAL,PTR2R8,'VOLPOND', __RC__) 
      call CICEReorder8(BUFINT8(L1),PTR2R8,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2R8,2),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR2R8,2)
      call MAPL_GetPointer(INTERNAL,PTR3R8,'ERGICE', __RC__) 
      call CICEReorder8(BUFINT8(L1),PTR3R8,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR3R8,2)*size(PTR3R8,3),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR3R8,2)*size(PTR3R8,3)
      call MAPL_GetPointer(INTERNAL,PTR3R8,'ERGSNO', __RC__) 
      call CICEReorder8(BUFINT8(L1),PTR3R8,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR3R8,2)*size(PTR3R8,3),UNPACKIT)

! retrieve export
      call MAPL_BalanceWork(BUFEXP, NUMMAX, Direction=MAPL_Retrieve, Handle=CICECOREBalanceHandle, __RC__)
      L1 = 1
      call MAPL_GetPointer(EXPORT,PTR1,'QH', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'TH', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'UH', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'VH', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'QST', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'TST', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'CHT', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'CMT', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'CQT', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'CNT', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'RIT', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'RET', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'Z0', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'Z0H', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOT2M', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOQ2M', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOU2M', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOV2M', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOT10M', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOQ10M', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOU10M', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOV10M', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOU50M', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MOV50M', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'GUST', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'VENT', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'TFREEZE', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FRACI', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'QSAT1', __RC__)
      if ( associated(PTR2) ) then
         call CICEReorder(BUFEXP(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
         L1 = L1 + NUMMAX*size(PTR2,2)
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'QSAT2', __RC__)
      if ( associated(PTR2) ) then
         call CICEReorder(BUFEXP(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
      end if
