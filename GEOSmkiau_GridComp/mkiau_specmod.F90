module mkiau_specmod
!$$$   module documentation block
!                .      .    .                                       .
! module:    specmod
!   prgmmr: treadon          org: np23                date: 2003-11-24
!
! abstract: module containing spectral related variables
!
! program history log:   
!   2003-11-24  treadon
!   2004-04-28  d. kokron, updated SGI's fft to use scsl
!   2004-05-18  kleist, documentation
!   2004-08-27  treadon - add/initialize variables/arrays needed by 
!                         splib routines for grid <---> spectral 
!                         transforms
!   2008-02-01  whitaker - modifications for use in ensemble kalman filter.
!   2010-10-27  whitaker - made thread safe (can now be called within OMP parallel regions).
!
! remarks: variable definitions below
!   def jcap         - spectral (assumed triangular) truncation
!   def nc           - (N+1)*(N+2); N=truncation
!   def ncd2         - [(N+1)*(N+2)]/2; N=truncation
!   def idrt         - integer grid identifier
!                      (idrt=4 for gaussian grid,
!                       idrt=0 for equally-spaced grid including poles,
!                       idrt=256 for equally-spaced grid excluding poles)
!   def imax         - integer even number of longitudes for transform
!   def jmax         - integer number of latitudes for transform
!   def jn           - integer skip number between n.h. latitudes from north
!   def js           - integer skip number between s.h. latitudes from south
!   def kw           - integer skip number between wave fields
!   def jb           - integer latitude index (from pole) to begin transform
!   def je           - integer latitude index (from pole) to end transform
!   def gaulats      - sin of latitudes on grid.
!   def gauwts       - gaussian weights.
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
! this version is thread safe (can be called within an OMP parallel region)
!
!$$$
  use MAPL_ConstantsMod, only: MAPL_I4, MAPL_R4
  implicit none

  integer(MAPL_I4) jcap
  integer(MAPL_I4) idrt,imax,jmax,jn,js,jb,je,ioffset,nc,ncd2,ijmax
  integer(MAPL_I4) :: iromb=0 ! only triangular truncation used
  real(MAPL_R4),allocatable,dimension(:):: eps,epstop,enn1,elonn1,eon,eontop
  real(MAPL_R4),allocatable,dimension(:):: clat,slat,wlat,gaulats,gauwts
  real(MAPL_R4),allocatable,dimension(:,:):: pln,plntop
  real(8),allocatable,dimension(:):: afft_save
  logical :: isinitialized=.false.

contains
  
  subroutine init_spec_vars(nlon,nlat,jcapin,idrtin)

!   Declare passed variables
    integer(MAPL_I4),intent(in):: nlat,nlon,jcapin,idrtin

!   Set constants
    jcap = jcapin
    idrt = idrtin
    imax = nlon
    jmax = nlat
    ijmax = imax*jmax
    ioffset=imax*(jmax-1)
    jn=imax
    js=-jn
    jb=1
    je=(jmax+1)/2
    nc=(jcap+1)*(jcap+2)
    ncd2=nc/2

!   Allocate arrays
    if (isinitialized) then
       call destroy_spec_vars()
    end if
    allocate( eps(ncd2) )
    allocate( epstop(jcap+1) )
    allocate( enn1(ncd2) )
    allocate( elonn1(ncd2) )
    allocate( eon(ncd2) )
    allocate( eontop(jcap+1) )
    allocate( afft_save(50000+4*imax) )
    allocate( gaulats(jmax) )
    allocate( gauwts(jmax) )
    allocate( clat(jb:je) )
    allocate( slat(jb:je) ) 
    allocate( wlat(jb:je) ) 
    allocate( pln(ncd2,jb:je) )
    allocate( plntop(jcap+1,jb:je) )

!   Initialize arrays used in transforms
    call sptranf0(iromb,jcap,idrt,imax,jmax,jb,je, &
       eps,epstop,enn1,elonn1,eon,eontop, &
       afft_save,clat,slat,wlat,pln,plntop)
    call splat(idrt,jmax,gaulats,gauwts)
    
    isinitialized = .true.

  end subroutine init_spec_vars

  subroutine destroy_spec_vars
    deallocate(eps,epstop,enn1,elonn1,eon,eontop,&
       clat,slat,wlat,pln,plntop,gaulats,gauwts,afft_save)
  end subroutine destroy_spec_vars

subroutine sptez_s(wave,grid,idir)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  sptez_s       perform a simple scalar spherical transform
!   prgmmr: iredell          org: np23                date: 1996-02-29
!
! absract: this subprogram performs a spherical transform
!           between spectral coefficients of a scalar quantity
!           and a field on a global cylindrical grid.
!           the wave-space is triangular.
!           the grid-space can be either an equally-spaced grid
!           (with or without pole points) or a gaussian grid.
!           the wave field is in sequential 'ibm order'.
!           the grid field is indexed east to west, then north to south.
!           for more flexibility and efficiency, call sptran.
!           subprogram can be called from a multiprocessing environment.
!
!           This routine differs from splib routine sptez in that
!              1) the calling list only contains the in/out arrays and 
!                 flag for the direction in which to transform
!              2) it calls a version of sptranf that does not invoke 
!                 initialization routines on each entry
!              3) some generality built into the splib version is
!                 removed in the code below
!
! program history log:
!   1996-02-29  iredell
!   2004-08-23  treadon - adapt splib routine sptez for gsi use
!   2008-02-01  whitaker - modifications for use in ensemble kalman filter.
!
!   input arguments:
!     wave     - real (2*mx) wave field if idir>0
!                where mx=(jcap+1)*(jcap+2)/2
!     grid     - real (imax,jmax) grid field (e->w,n->s) if idir<0
!     idir     - integer transform flag
!                (idir>0 for wave to grid, idir<0 for grid to wave)
!
!   output arguments:
!     wave     - real (2*mx) wave field if idir<0
!                where mx=(maxwv+1)*(maxwv+2)/2
!     grid     - real (imax,jmax) grid field (e->w,n->s) if idir>0
!
! subprograms called:
!   sptranf_s  -  perform a scalar spherical transform
!
! remarks: minimum grid dimensions for unaliased transforms to spectral:
!   dimension                    linear              quadratic
!   -----------------------      ---------           -------------
!   imax                         2*maxwv+2           3*maxwv/2*2+2
!   jmax (idrt=4)                1*maxwv+1           3*maxwv/2+1
!   jmax (idrt=0)                2*maxwv+3           3*maxwv/2*2+3
!   jmax (idrt=256)              2*maxwv+1           3*maxwv/2*2+1
!   -----------------------      ---------           -------------
!
! attributes:
!   language: fortran 77
!
!$$$
  implicit none

! Declare passed variables
  integer(MAPL_I4),intent(in):: idir
  real(MAPL_R4),dimension(nc),intent(inout):: wave
  real(MAPL_R4),dimension(ijmax),intent(inout):: grid

! Declare local variables
  integer(MAPL_I4) i

  if (.not. isinitialized) then
     print *,'call init_spec_vars first to initialize'
     stop
  end if
! Zero appropriate output array based on direction of transform
  if (idir < 0) then
     do i=1,nc
        wave(i)=0._MAPL_R4
     end do
  elseif (idir > 0) then
     do i=1,ijmax
        grid(i)=0._MAPL_R4
     end do
  endif

! Call spectral <--> grid transform
  call sptranf_s(wave,grid,grid,idir)

  return
end subroutine sptez_s

subroutine sptranf_s(wave,gridn,grids,idir)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  sptranf_s     perform a scalar spherical transform
!   prgmmr: iredell          org: np23                date: 1996-02-29
!
! abstract: this subprogram performs a spherical transform
!           between spectral coefficients of scalar quantities
!           and fields on a global cylindrical grid.
!           the wave-space is triangular.
!           the grid-space can be either an equally-spaced grid
!           (with or without pole points) or a gaussian grid.
!           the wave and grid fields may have general indexing,
!           but each wave field is in sequential 'ibm order',
!           i.e. with zonal wavenumber as the slower index.
!           transforms are done in latitude pairs for efficiency;
!           thus grid arrays for each hemisphere must be passed.
!           if so requested, just a subset of the latitude pairs
!           may be transformed in each invocation of the subprogram.
!           the transforms are all multiprocessed over latitude except
!           the transform from fourier to spectral is multiprocessed
!           over zonal wavenumber to ensure reproducibility.
!           transform several fields at a time to improve vectorization.
!           subprogram can be called from a multiprocessing environment.
!
!           This routine differs from splib routine sptranf in that
!           it does not call sptranf0 (an initialization routine).
!
! program history log:
!   1996-02-29  iredell
!   1998-12-15  iredell  generic fft used
!   2004-08-23  treadon - adapt splib routine sptranf for gsi use
!   2006-05-03  treadon - remove jc from specmod list since not used
!   2008-02-01  whitaker - modifications for use in ensemble kalman filter.
!
!   input arguments:
!     wave     - real (*) wave fields if idir>0
!     gridn    - real (*) n.h. grid fields (starting at jb) if idir<0
!     grids    - real (*) s.h. grid fields (starting at jb) if idir<0
!     idir     - integer transform flag
!                (idir>0 for wave to grid, idir<0 for grid to wave)
!
!   output arguments:
!     wave     - real (*) wave fields if idir<0
!     gridn    - real (*) n.h. grid fields (starting at jb) if idir>0
!     grids    - real (*) s.h. grid fields (starting at jb) if idir>0
!
! subprograms called:
!   sptranf1     sptranf spectral transform
!   
! remarks: 
!   This routine assumes that splib routine sptranf0 has been 
!   previously called.  sptranf0 initializes arrays needed in
!   the transforms.
!
!   minimum grid dimensions for unaliased transforms to spectral:
!   dimension                    linear              quadratic
!   -----------------------      ---------           -------------
!   imax                         2*maxwv+2           3*maxwv/2*2+2
!   jmax (idrt=4)                1*maxwv+1           3*maxwv/2+1
!   jmax (idrt=0)                2*maxwv+3           3*maxwv/2*2+3
!   jmax (idrt=256)              2*maxwv+1           3*maxwv/2*2+1
!   -----------------------      ---------           -------------
!
! attributes:
!   language: fortran 77
!
!$$$
  implicit none

! Declare passed variables
  integer(MAPL_I4),intent(in):: idir
  real(MAPL_R4),dimension(nc),intent(inout):: wave
  real(MAPL_R4),dimension(ijmax),intent(inout):: gridn
  real(MAPL_R4),dimension(ijmax),intent(inout):: grids

! Declare local variables
  integer(MAPL_I4) i,j,jj,ijn,ijs,mp
  real(MAPL_R4),dimension(2*(jcap+1)):: wtop
  real(MAPL_R4),dimension(imax,2):: g
  real(8),dimension(50000+4*imax):: afft

  if (.not. isinitialized) then
     print *,'call init_spec_vars first to initialize'
     stop
  end if

! this is needed for thread safety.
  afft = afft_save

! Initialize local variables
  mp=0

  do i=1,2*(jcap+1)
     wtop(i)=0._MAPL_R4
  end do

! Transform wave to grid
  if(idir > 0) then
     do j=jb,je
        call sptranf1(iromb,jcap,idrt,imax,jmax,j,j, &
             eps,epstop,enn1,elonn1,eon,eontop, &
             afft,clat(j),slat(j),wlat(j), &
             pln(1,j),plntop(1,j),mp, &
             wave,wtop,g,idir)
        do i=1,imax
           jj  = j-jb
           ijn = i + jj*jn
           ijs = i + jj*js + ioffset
           gridn(ijn)=g(i,1)
           grids(ijs)=g(i,2)
        enddo
     enddo

! Transform grid to wave
  else
     do j=jb,je
        if(wlat(j) > 0._MAPL_R4) then
           do i=1,imax
              jj  = j-jb
              ijn = i + jj*jn
              ijs = i + jj*js + ioffset
              g(i,1)=gridn(ijn)
              g(i,2)=grids(ijs)
              
           enddo
           call sptranf1(iromb,jcap,idrt,imax,jmax,j,j, &
                eps,epstop,enn1,elonn1,eon,eontop, &
                afft,clat(j),slat(j),wlat(j), &
                pln(1,j),plntop(1,j),mp, &
                wave,wtop,g,idir)
        endif
     enddo
  endif
end subroutine sptranf_s


subroutine sptezv_s(waved,wavez,gridu,gridv,idir)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  sptez_v       perform a simple vector spherical transform
!   prgmmr: iredell          org: np23                date: 1996-02-29
!
! abstract: this subprogram performs a spherical transform
!           between spectral coefficients of divergence and curl
!           and a vector field on a global cylindrical grid.
!           the wave-space is triangular.
!           the grid-space can be either an equally-spaced grid
!           (with or without pole points) or a gaussian grid.
!           the wave field is in sequential 'ibm order'.
!           the grid fiels is indexed east to west, then north to south.
!           for more flexibility and efficiency, call sptran.
!           subprogram can be called from a multiprocessing environment.
!
!           This routine differs from splib routine sptezv in that
!              1) the calling list only contains the in/out arrays and
!                 flag for the direction in which to transform
!              2) it calls a version of sptranfv that does not invoke
!                 initialization routines on each entry
!              3) some generality built into the splib version is
!                 removed in the code below
!
! program history log:
!   1996-02-29  iredell
!   2004-08-23  treadon - adapt splib routine sptezv for gsi use
!   2008-02-01  whitaker - modifications for use in ensemble kalman filter.
!
!   input arguments:
!     waved    - real (2*mx) wave divergence field if idir>0
!                where mx=(maxwv+1)*(maxwv+2)/2
!     wavez    - real (2*mx) wave vorticity field if idir>0
!                where mx=(maxwv+1)*(maxwv+2)/2
!     gridu    - real (imax,jmax) grid u-wind (e->w,n->s) if idir<0
!     gridv    - real (imax,jmax) grid v-wind (e->w,n->s) if idir<0
!     idir     - integer transform flag
!                (idir>0 for wave to grid, idir<0 for grid to wave)
!
!   output arguments:
!     waved    - real (2*mx) wave divergence field if idir<0
!                where mx=(maxwv+1)*(maxwv+2)/2
!     wavez    - real (2*mx) wave vorticity field if idir>0
!                where mx=(maxwv+1)*(maxwv+2)/2
!     gridu    - real (imax,jmax) grid u-wind (e->w,n->s) if idir>0
!     gridv    - real (imax,jmax) grid v-wind (e->w,n->s) if idir>0
!
! subprograms called:
!   sptranf_v  - perform a vector spherical transform
!
! remarks: minimum grid dimensions for unaliased transforms to spectral:
!   dimension                    linear              quadratic
!   -----------------------      ---------           -------------
!   imax                         2*maxwv+2           3*maxwv/2*2+2
!   jmax (idrt=4)                1*maxwv+1           3*maxwv/2+1
!   jmax (idrt=0)                2*maxwv+3           3*maxwv/2*2+3
!   jmax (idrt=256)              2*maxwv+1           3*maxwv/2*2+1
!   -----------------------      ---------           -------------
!
! attributes:
!   language: fortran 77
!
!$$$
  implicit none

! Declare passed variables
  integer(MAPL_I4),intent(in):: idir
  real(MAPL_R4),dimension(nc),intent(inout):: waved,wavez
  real(MAPL_R4),dimension(ijmax),intent(inout):: gridu,gridv

! Declare local variables
  integer(MAPL_I4) i

! Zero appropriate output array based on direction of transform
  if (idir < 0) then
     do i=1,nc
        waved(i)=0._MAPL_R4
        wavez(i)=0._MAPL_R4
     end do
  elseif (idir > 0) then
     do i=1,ijmax
        gridu(i)=0._MAPL_R4
        gridv(i)=0._MAPL_R4
     end do
  endif

! Call spectral <--> grid transform
  call sptranf_v(waved,wavez,gridu,gridu,gridv,gridv,idir)

end subroutine sptezv_s

subroutine sptranf_v(waved,wavez,gridun,gridus,gridvn,gridvs,idir)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  sptranf_v     perform a vecor spherical transform
!   prgmmr: iredell          org: np23                date: 1996-02-29
!
! abstract: this subprogram performs a spherical transform
!           between spectral coefficients of divergences and curls
!           and vector fields on a global cylindrical grid.
!           the wave-space is triangular.
!           the grid-space can be either an equally-spaced grid
!           (with or without pole points) or a gaussian grid.
!           the wave and grid fields may have general indexing,
!           but each wave field is in sequential 'ibm order',
!           i.e. with zonal wavenumber as the slower index.
!           transforms are done in latitude pairs for efficiency;
!           thus grid arrays for each hemisphere must be passed.
!           if so requested, just a subset of the latitude pairs
!           may be transformed in each invocation of the subprogram.
!           the transforms are all multiprocessed over latitude except
!           the transform from fourier to spectral is multiprocessed
!           over zonal wavenumber to ensure reproducibility.
!           transform several fields at a time to improve vectorization.
!           subprogram can be called from a multiprocessing environment.
!
!           This routine differs from splib routine sptranfv in that
!           it does not call sptranf0 (an initialization routine).
!
! program history log:
!   1996-02-29  iredell
!   1998-12-15  iredell  generic fft used
!   2004-08-23  treadon - adapt splib routine sptranfv for gsi use
!   2006-05-03  treadon - remove jc from specmod list since not used
!   2006-07-07  kleist - correct bug in indexing of j=1,2*ncd2 loop
!   2008-02-01  whitaker - modifications for use in ensemble kalman filter.
!
!   input arguments:
!     waved    - real (*) wave divergence fields if idir>0
!     wavez    - real (*) wave vorticity fields if idir>0
!     gridun   - real (*) n.h. grid u-winds (starting at jb) if idir<0
!     gridus   - real (*) s.h. grid u-winds (starting at jb) if idir<0
!     gridvn   - real (*) n.h. grid v-winds (starting at jb) if idir<0
!     gridvs   - real (*) s.h. grid v-winds (starting at jb) if idir<0
!     idir     - integer transform flag
!                (idir>0 for wave to grid, idir<0 for grid to wave)
!
!   output arguments:
!     waved    - real (*) wave divergence fields if idir<0
!                [waved=(d(gridu)/dlam+d(clat*gridv)/dphi)/(clat*rerth)]
!     wavez    - real (*) wave vorticity fields if idir<0
!                [wavez=(d(gridv)/dlam-d(clat*gridu)/dphi)/(clat*rerth)]
!     gridun   - real (*) n.h. grid u-winds (starting at jb) if idir>0
!     gridus   - real (*) s.h. grid u-winds (starting at jb) if idir>0
!     gridvn   - real (*) n.h. grid v-winds (starting at jb) if idir>0
!     gridvs   - real (*) s.h. grid v-winds (starting at jb) if idir>0
!
! subprograms called:
!   sptranf1     sptranf spectral transform
!   spdz2uv      compute winds from divergence and vorticity
!   spuv2dz      compute divergence and vorticity from winds
!
! remarks: 
!   This routine assumes that splib routine sptranf0 has been
!   previously called.  sptranf0 initializes arrays needed in
!   the transforms.
!
!   minimum grid dimensions for unaliased transforms to spectral:
!   dimension                    linear              quadratic
!   -----------------------      ---------           -------------
!   imax                         2*maxwv+2           3*maxwv/2*2+2
!   jmax (idrt=4)                1*maxwv+1           3*maxwv/2+1
!   jmax (idrt=0)                2*maxwv+3           3*maxwv/2*2+3
!   jmax (idrt=256)              2*maxwv+1           3*maxwv/2*2+1
!   -----------------------      ---------           -------------
!
! attributes:
!   language: fortran 77
!
!$$$
  implicit none

! Declare passed variables
  integer(MAPL_I4),intent(in):: idir
  real(MAPL_R4),dimension(nc):: waved,wavez
  real(MAPL_R4),dimension(ijmax):: gridun,gridus,gridvn,gridvs


! Declare local variables
  integer(MAPL_I4) i,j,jj,ijn,ijs
  integer(MAPL_I4),dimension(2):: mp
  real(MAPL_R4),dimension(ncd2*2,2):: w
  real(MAPL_R4),dimension(2*(jcap+1),2):: wtop
  real(MAPL_R4),dimension(imax,2,2):: g
  real(8),dimension(50000+4*imax):: afft

! Set parameters
  mp=1
! this is needed for thread safety.
  afft = afft_save

! Transform wave to grid
  if(idir > 0) then
     call spdz2uv(iromb,jcap,enn1,elonn1,eon,eontop, &
          waved,wavez, &
          w(1,1),w(1,2),wtop(1,1),wtop(1,2))
     do j=jb,je
        call sptranf1(iromb,jcap,idrt,imax,jmax,j,j, &
             eps,epstop,enn1,elonn1,eon,eontop, &
             afft,clat(j),slat(j),wlat(j), &
             pln(1,j),plntop(1,j),mp, &
             w(1,1),wtop(1,1),g(1,1,1),idir)
        call sptranf1(iromb,jcap,idrt,imax,jmax,j,j, &
             eps,epstop,enn1,elonn1,eon,eontop, &
             afft,clat(j),slat(j),wlat(j), &
             pln(1,j),plntop(1,j),mp, &
             w(1,2),wtop(1,2),g(1,1,2),idir)
        do i=1,imax
           jj   = j-jb
           ijn = i + jj*jn
           ijs = i + jj*js + ioffset
           gridun(ijn)=g(i,1,1)
           gridus(ijs)=g(i,2,1)
           gridvn(ijn)=g(i,1,2)
           gridvs(ijs)=g(i,2,2)
           
        enddo
     enddo

!  Transform grid to wave
  else
     w=0
     wtop=0
     do j=jb,je
        if(wlat(j) > 0._MAPL_R4) then
           do i=1,imax
              jj   = j-jb
              ijn = i + jj*jn
              ijs = i + jj*js + ioffset

              g(i,1,1)=gridun(ijn)/clat(j)**2
              g(i,2,1)=gridus(ijs)/clat(j)**2
              g(i,1,2)=gridvn(ijn)/clat(j)**2
              g(i,2,2)=gridvs(ijs)/clat(j)**2
           enddo
           call sptranf1(iromb,jcap,idrt,imax,jmax,j,j, &
                eps,epstop,enn1,elonn1,eon,eontop, &
                afft,clat(j),slat(j),wlat(j), &
                pln(1,j),plntop(1,j),mp, &
                w(1,1),wtop(1,1),g(1,1,1),idir)
           call sptranf1(iromb,jcap,idrt,imax,jmax,j,j, &
                eps,epstop,enn1,elonn1,eon,eontop, &
                afft,clat(j),slat(j),wlat(j), &
                pln(1,j),plntop(1,j),mp, &
                w(1,2),wtop(1,2),g(1,1,2),idir)
        endif
     enddo
     call spuv2dz(iromb,jcap,enn1,elonn1,eon,eontop, &
          w(1,1),w(1,2),wtop(1,1),wtop(1,2), &
          waved(1),wavez(1))
  endif

 end subroutine sptranf_v

end module mkiau_specmod
