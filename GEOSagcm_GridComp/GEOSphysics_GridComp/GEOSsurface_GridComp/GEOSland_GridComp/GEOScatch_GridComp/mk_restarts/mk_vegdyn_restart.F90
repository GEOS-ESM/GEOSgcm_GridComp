PROGRAM mk_vegdyn_internal
implicit none
real, allocatable :: dummy(:),ity0(:)
integer,allocatable  :: ity0_int(:) 
real :: filler, dum0
integer :: vvv
integer ::  bi,li
integer :: nland, nt, index, id, dum
character*256  outpath, sarithdir, dateline, restag, vegname, numland
character*256  landir

!---------------------------------------------------------------------------
    call GETENV ( 'LANDIR'  , landir   )
    call GETENV ( 'rslv'    , restag   )
    call GETENV ( 'dateline', dateline )
    call GETENV ( 'nland'   , numland  )
    read(numland,*)nland

outpath  = 'vegdyn_internal_restart.' // trim(restag) // '_' // trim(dateline)
!---------------------------------------------------------------------------

allocate(dummy   (nland))
allocate(ity0    (nland))
allocate(ity0_int(nland))

dummy(:)=-999.0
sarithdir = trim(landir) // '/' // trim(dateline) // '/FV_' // trim(restag) // '/'
vegname   = trim(sarithdir)//'mosaic_veg_typs_fracs'
write (*,*) 'Reading '//vegname

open(unit=21, file=trim(vegname),form='formatted')
DO nt=1,nland
!  read (21, *) index, id, ity0_int(nt), dum, dum0, dum0, dum0
   read (21, *) index, id, ity0_int(nt), dum, dum0, dum0  ! version 2 doesn't have frc3
   print *, ity0_int(nt)
ENDDO
ity0=ity0_int*1.0
close(21)


open(unit=30, file=trim(outpath),form='unformatted')
! write out dummy lai_prev, lai_next, grn_prev, grn_next
print *, '    VEGTYPES', minval(ity0), maxval(ity0)
write (30) dummy
write (30) dummy
write (30) dummy
write (30) dummy
write (30) ity0
close (30)
deallocate(ity0)
deallocate(ity0_int)
deallocate(dummy)
END

