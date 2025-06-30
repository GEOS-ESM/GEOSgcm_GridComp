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
      call MAPL_GetPointer(INTERNAL,PTR2,'CQ', _RC) 
      call CICEReorder(BUFINT(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR2,2)
      call MAPL_GetPointer(INTERNAL,PTR2,'CM', _RC) 
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
      call MAPL_GetPointer(INTERNAL,PTR2R8,'APONDN', _RC) 
      call CICEReorder8(BUFINT8(L1),PTR2R8,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2R8,2),UNPACKIT)
      L1 = L1 + NUMMAX*size(PTR2R8,2)
      call MAPL_GetPointer(INTERNAL,PTR2R8,'HPONDN', _RC) 
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
      call MAPL_GetPointer(EXPORT,PTR1,'EMIS', _RC)
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      call MAPL_GetPointer(EXPORT,PTR1,'ALBVF', _RC)
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      call MAPL_GetPointer(EXPORT,PTR1,'ALBVR', _RC)
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      call MAPL_GetPointer(EXPORT,PTR1,'ALBNF', _RC)
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      call MAPL_GetPointer(EXPORT,PTR1,'ALBNR', _RC)
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
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
      call MAPL_GetPointer(EXPORT,PTR1,'DELTS', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'DELQS', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'TAUXI', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'TAUYI', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'PENUVR', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'PENUVF', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'PENPAR', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'PENPAF', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'EVAPOUT', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SUBLIM', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SHOUT', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SHICE', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'HLATN', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'HLATICE', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FSURF', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FSURFICE', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'HLWUP', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'HLWUPICE', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'LWNDSRF', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SWNDSRF', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'LWNDICE', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SWNDICE', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FRACINEW', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'LWDNSRF', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SWDNSRF', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FRAZIL', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'CONGEL', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SNOICE', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FRESH', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FSALT', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FHOCN', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'PICE', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FSWTHRU', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FSWABS', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MELTL', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MELTT', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MELTB', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'MELTS', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'HICE', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'HSNO', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'HICEUNT', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SNOONICE', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'TSKINICE', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'IAGE', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'DAIDTT', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'DVIDTT', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FBOT', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'USTARI', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FCONDTOP', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FCONDBOT', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'NEWICEERG', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SUBLIMFLX', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'SIALB', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'GHTSKIN', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'FRZMLT', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'evap_CMIP5', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'pr_CMIP5', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'prsn_CMIP5', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'grFrazil_CMIP5', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'grCongel_CMIP5', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'grLateral_CMIP5', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'snoToIce_CMIP5', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'snomelt_CMIP5', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'tmelt_CMIP5', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'bmelt_CMIP5', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'sfdsi_CMIP5', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'hfsifrazil_CMIP5', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'ialb_CMIP5', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'rsdssi_CMIP5', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'rsussi_CMIP5', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR1,'fsitherm_CMIP5', _RC)
      if ( associated(PTR1) ) then
         call CICEReorder(BUFEXP(L1),PTR1,TILE_WITH_ICE,NUMMAX,HorzDims,1,UNPACKIT)
         L1 = L1 + NUMMAX
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'FCONDBOTN', _RC)
      if ( associated(PTR2) ) then
         call CICEReorder(BUFEXP(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
         L1 = L1 + NUMMAX*size(PTR2,2)
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'FCONDTOPN', _RC)
      if ( associated(PTR2) ) then
         call CICEReorder(BUFEXP(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
         L1 = L1 + NUMMAX*size(PTR2,2)
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'DELTAVOL1', _RC)
      if ( associated(PTR2) ) then
         call CICEReorder(BUFEXP(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
         L1 = L1 + NUMMAX*size(PTR2,2)
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'TINZ', _RC)
      if ( associated(PTR2) ) then
         call CICEReorder(BUFEXP(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
         L1 = L1 + NUMMAX*size(PTR2,2)
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'SHICEN', _RC)
      if ( associated(PTR2) ) then
         call CICEReorder(BUFEXP(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
         L1 = L1 + NUMMAX*size(PTR2,2)
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'HLWUPN', _RC)
      if ( associated(PTR2) ) then
         call CICEReorder(BUFEXP(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
         L1 = L1 + NUMMAX*size(PTR2,2)
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'LWNDSRFN', _RC)
      if ( associated(PTR2) ) then
         call CICEReorder(BUFEXP(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
         L1 = L1 + NUMMAX*size(PTR2,2)
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'FSURFN', _RC)
      if ( associated(PTR2) ) then
         call CICEReorder(BUFEXP(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
         L1 = L1 + NUMMAX*size(PTR2,2)
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'TSURFN', _RC)
      if ( associated(PTR2) ) then
         call CICEReorder(BUFEXP(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
         L1 = L1 + NUMMAX*size(PTR2,2)
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'FSWSFCN', _RC)
      if ( associated(PTR2) ) then
         call CICEReorder(BUFEXP(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
         L1 = L1 + NUMMAX*size(PTR2,2)
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'ALBIN', _RC)
      if ( associated(PTR2) ) then
         call CICEReorder(BUFEXP(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
         L1 = L1 + NUMMAX*size(PTR2,2)
      end if
      call MAPL_GetPointer(EXPORT,PTR2,'ALBSN', _RC)
      if ( associated(PTR2) ) then
         call CICEReorder(BUFEXP(L1),PTR2,TILE_WITH_ICE,NUMMAX,HorzDims,size(PTR2,2),UNPACKIT)
      end if
