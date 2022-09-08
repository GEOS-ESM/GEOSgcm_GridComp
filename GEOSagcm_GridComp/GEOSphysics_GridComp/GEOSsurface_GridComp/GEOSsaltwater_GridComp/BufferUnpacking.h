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
      call MAPL_GetPointer(INTERNAL,PTR2,'CQ', __RC__) 
      call CICEReorder(BUFINT(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR2,2)
      call MAPL_GetPointer(INTERNAL,PTR2,'CM', __RC__) 
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
      call MAPL_GetPointer(INTERNAL,PTR2R8,'APONDN', __RC__) 
      call CICEReorder8(BUFINT8(L1),PTR2R8,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2R8,2),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR2R8,2)
      call MAPL_GetPointer(INTERNAL,PTR2R8,'HPONDN', __RC__) 
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
      call MAPL_GetPointer(EXPORT,PTR1,'EMIS', __RC__)
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      call MAPL_GetPointer(EXPORT,PTR1,'ALBVF', __RC__)
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      call MAPL_GetPointer(EXPORT,PTR1,'ALBVR', __RC__)
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      call MAPL_GetPointer(EXPORT,PTR1,'ALBNF', __RC__)
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      call MAPL_GetPointer(EXPORT,PTR1,'ALBNR', __RC__)
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
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
      call MAPL_GetPointer(EXPORT,PTR1,'DELTS', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'DELQS', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'TAUXI', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'TAUYI', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'PENUVR', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'PENUVF', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'PENPAR', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'PENPAF', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'EVAPOUT', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SUBLIM', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SHOUT', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SHICE', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'HLATN', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'HLATICE', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FSURF', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FSURFICE', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'HLWUP', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'HLWUPICE', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'LWNDSRF', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SWNDSRF', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'LWNDICE', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SWNDICE', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FRACINEW', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'LWDNSRF', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SWDNSRF', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FRAZIL', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'CONGEL', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SNOICE', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FRESH', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FSALT', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FHOCN', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'PICE', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FSWTHRU', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FSWABS', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MELTL', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MELTT', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MELTB', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MELTS', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'HICE', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'HSNO', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'HICEUNT', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SNOONICE', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'TSKINICE', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'IAGE', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'DAIDTT', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'DVIDTT', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FBOT', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'USTARI', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FCONDTOP', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FCONDBOT', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'NEWICEERG', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SUBLIMFLX', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SIALB', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'GHTSKIN', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FRZMLT', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'evap_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'pr_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'prsn_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'grFrazil_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'grCongel_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'grLateral_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'snoToIce_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'snomelt_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'tmelt_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'bmelt_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'sfdsi_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'hfsifrazil_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'ialb_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'rsdssi_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'rsussi_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'fsitherm_CMIP5', __RC__)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'FCONDBOTN', __RC__)
      if ( associated(PTR2) ) then
         call CICEReorder(BUFEXP(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
         L1 = L1 + NUMMAX*size(PTR2,2)
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'FCONDTOPN', __RC__)
      if ( associated(PTR2) ) then
         call CICEReorder(BUFEXP(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
         L1 = L1 + NUMMAX*size(PTR2,2)
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'TINZ', __RC__)
      if ( associated(PTR2) ) then
         call CICEReorder(BUFEXP(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
         L1 = L1 + NUMMAX*size(PTR2,2)
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'SHICEN', __RC__)
      if ( associated(PTR2) ) then
         call CICEReorder(BUFEXP(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
         L1 = L1 + NUMMAX*size(PTR2,2)
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'HLWUPN', __RC__)
      if ( associated(PTR2) ) then
         call CICEReorder(BUFEXP(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
         L1 = L1 + NUMMAX*size(PTR2,2)
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'LWNDSRFN', __RC__)
      if ( associated(PTR2) ) then
         call CICEReorder(BUFEXP(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
         L1 = L1 + NUMMAX*size(PTR2,2)
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'FSURFN', __RC__)
      if ( associated(PTR2) ) then
         call CICEReorder(BUFEXP(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
         L1 = L1 + NUMMAX*size(PTR2,2)
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'TSURFN', __RC__)
      if ( associated(PTR2) ) then
         call CICEReorder(BUFEXP(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
         L1 = L1 + NUMMAX*size(PTR2,2)
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'FSWSFCN', __RC__)
      if ( associated(PTR2) ) then
         call CICEReorder(BUFEXP(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
         L1 = L1 + NUMMAX*size(PTR2,2)
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'ALBIN', __RC__)
      if ( associated(PTR2) ) then
         call CICEReorder(BUFEXP(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
         L1 = L1 + NUMMAX*size(PTR2,2)
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'ALBSN', __RC__)
      if ( associated(PTR2) ) then
         call CICEReorder(BUFEXP(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
      end if
