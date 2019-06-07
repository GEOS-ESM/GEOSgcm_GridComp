#include "MAPL_Generic.h"

module fv_regridding_utils

   use ESMF 
   use fv_arrays_mod,     only: fv_atmos_type, fv_grid_type, fv_grid_bounds_type, FVPRC, REAL4, REAL8
   use fv_diagnostics_mod,only: prt_maxmin
   use fv_mp_mod,         only: is_master, ng
   use fv_mapz_mod,       only: mappm
   use mpp_mod,            only: mpp_error, FATAL, NOTE, mpp_broadcast,mpp_npes
   !use MAPL_MOD,          only: MAPL_PI_R8, MAPL_OMEGA, MAPL_GRAV, &
         !MAPL_KAPPA, MAPL_RGAS, MAPL_RVAP, &
         !MAPL_CP
   use MAPL_MOD

   implicit none

   private

   public remap_scalar
   public fv_rst
   public copy_fv_rst
   public simple_cs_grid_creator

   real(FVPRC), parameter :: PI           = MAPL_PI_R8
   real(FVPRC), parameter :: OMEGA        = MAPL_OMEGA
   real(FVPRC), parameter :: GRAV         = MAPL_GRAV
   real(FVPRC), parameter :: KAPPA        = MAPL_KAPPA
   real(FVPRC), parameter :: RDGAS        = MAPL_RGAS
   real(FVPRC), parameter :: RVGAS        = MAPL_RVAP
   real(FVPRC), parameter :: CP_AIR       = MAPL_CP
   real(FVPRC), parameter:: zvir = rvgas/rdgas - 1.

   type fv_var
      character(len=128)   :: name
      integer              :: nlev
      real(REAL8), pointer :: ptr2d(:,:) => null()
      real(REAL8), pointer :: ptr3d(:,:,:) => null()
   end type fv_var

   type fv_rst
      character(len=1024)   :: file_name
      logical               :: isBin
      logical               :: have_descriptor
      type(fv_var), pointer :: vars(:) => null()
   end type fv_rst
      

contains

 subroutine copy_fv_rst(in_rst,out_rst)
  type(fv_rst), pointer, intent(inout) :: in_rst(:)
  type(fv_rst), pointer, intent(inout) :: out_rst(:)
  
  integer :: ifile,ivar
  allocate(out_rst(size(in_rst)) )
  do ifile=1,size(in_rst)
     allocate( out_rst(ifile)%vars(size(in_rst(ifile)%vars) ) )
     out_rst(ifile)%isBin=in_rst(ifile)%isBin
     out_rst(ifile)%file_name=in_rst(ifile)%file_name
     out_rst(ifile)%have_descriptor=in_rst(ifile)%have_descriptor
     do ivar=1,size(in_rst(ifile)%vars)
        out_rst(ifile)%vars(ivar)%name=in_rst(ifile)%vars(ivar)%name
        out_rst(ifile)%vars(ivar)%nlev=in_rst(ifile)%vars(ivar)%nlev
     enddo
  enddo
   
 end subroutine copy_fv_rst

 subroutine remap_scalar(im, jm, km, npz, nq, ncnst, ak0, bk0, psc, gzc, ta, qa, Atm, in_fv_rst,out_fv_rst)
  type(fv_atmos_type), intent(inout) :: Atm
  integer, intent(in):: im, jm, km, npz, nq, ncnst
  real(FVPRC),    intent(in):: ak0(km+1), bk0(km+1)
  real(FVPRC),    intent(in), dimension(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je):: psc, gzc
  real(FVPRC),    intent(in), dimension(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je,km):: ta
  real(FVPRC),    intent(in), dimension(Atm%bd%is:Atm%bd%ie,Atm%bd%js:Atm%bd%je,km,ncnst):: qa
  type(fv_rst), pointer,   intent(inout) :: in_fv_rst(:)
  type(fv_rst), pointer,   intent(inout) :: out_fv_rst(:)
! local:
  real(FVPRC), dimension(Atm%bd%is:Atm%bd%ie,km):: tp
  real(FVPRC), dimension(Atm%bd%is:Atm%bd%ie,km+1):: pe0, pn0
  real(FVPRC), dimension(Atm%bd%is:Atm%bd%ie,npz):: qn1
  real(FVPRC), dimension(Atm%bd%is:Atm%bd%ie,npz+1):: pe1, pn1
  real(FVPRC), dimension(Atm%bd%is:Atm%bd%ie,npz+1):: q_edge_old,q_edge_new
  real(FVPRC) pt0(km), gz(km+1), pk0(km+1)
  real(FVPRC) qp( Atm%bd%is:Atm%bd%ie,km,ncnst)
  real(FVPRC) qp1(Atm%bd%is:Atm%bd%ie,km)
  real(FVPRC) pst, p1, p2, alpha, rdg
  integer i,j,k, iq
  integer  sphum,ifile,ivar
  integer :: is,  ie,  js,  je
  integer :: isd, ied, jsd, jed
  logical :: doVert

  is  = Atm%bd%is
  ie  = Atm%bd%ie
  js  = Atm%bd%js
  je  = Atm%bd%je
  isd = Atm%bd%isd
  ied = Atm%bd%ied
  jsd = Atm%bd%jsd
  jed = Atm%bd%jed
  sphum   = 1
  if ( sphum/=1 ) then
       call mpp_error(FATAL,'SPHUM must be 1st tracer')
  endif

  do j=js,je
     do i=is,ie

       do iq=1,ncnst
          do k=1,km
             qp(i,k,iq) = qa(i,j,k,iq)
          enddo
       enddo
       do k=1,km
          tp(i,k) = ta(i,j,k)*(1.+zvir*qp(i,k,sphum))
       enddo

! Tracers:

       do k=1,km+1
          pe0(i,k) = ak0(k) + bk0(k)*psc(i,j)
          pn0(i,k) = log(pe0(i,k))
            pk0(k) = pe0(i,k)**kappa
       enddo

! * Adjust interpolated ps to model terrain
       gz(km+1) = gzc(i,j)
       do k=km,1,-1
           gz(k) = gz(k+1) + rdgas*tp(i,k)*(pn0(i,k+1)-pn0(i,k))
       enddo
! Only lowest layer potential temp is needed
          pt0(km) = tp(i,km)/(pk0(km+1)-pk0(km))*(kappa*(pn0(i,km+1)-pn0(i,km)))
       if( Atm%phis(i,j)>gzc(i,j) ) then
           do k=km,1,-1
              if( Atm%phis(i,j) <  gz(k)  .and.    &
                  Atm%phis(i,j) >= gz(k+1) ) then
                  pst = pk0(k) + (pk0(k+1)-pk0(k))*(gz(k)-Atm%phis(i,j))/(gz(k)-gz(k+1))
                  go to 123
              endif
           enddo
       else
! Extrapolation into the ground
           pst = pk0(km+1) + (gzc(i,j)-Atm%phis(i,j))/(cp_air*pt0(km))
       endif

123    Atm%ps(i,j) = pst**(1./kappa)

     enddo   !i-loop

     do i=is,ie
        pe1(i,1) = Atm%ak(1)
        pn1(i,1) = log(pe1(i,1))
     enddo
     do k=2,npz+1
       do i=is,ie
          pe1(i,k) = Atm%ak(k) + Atm%bk(k)*Atm%ps(i,j)
          pn1(i,k) = log(pe1(i,k))
       enddo
     enddo

! * Compute delp
     do k=1,npz
        do i=is,ie
           Atm%delp(i,j,k) = pe1(i,k+1) - pe1(i,k)
        enddo
     enddo

!---------------
! map tracers
!----------------
     do iq=1,ncnst
        call mappm(km, pe0, qp(is,1,iq), npz, pe1,  qn1, is,ie, 0, 11, Atm%ptop)
        do k=1,npz
           do i=is,ie
              Atm%q(i,j,k,iq) = qn1(i,k)
           enddo
        enddo
     enddo

!---------------
! map extra 3d variables
!----------------
     
     do ifile=1,size(out_fv_rst)

        if (out_fv_rst(ifile)%have_descriptor) then
           do ivar=1,size(out_fv_rst(ifile)%vars)
              if (out_fv_rst(ifile)%vars(ivar)%nLev==npz) then
                 do k=1,in_fv_rst(ifile)%vars(ivar)%nLev
                    qp1(is:ie,k)=in_fv_rst(ifile)%vars(ivar)%ptr3d(is:ie,j,k)
                 enddo
                 call mappm(km, pe0, qp1, npz, pe1,  qn1, is,ie, 0, 11, Atm%ptop)
                 do k=1,npz
                    do i=is,ie
                       out_fv_rst(ifile)%vars(ivar)%ptr3d(i,j,k) = qn1(i,k)
                    enddo
                 enddo
              else if (out_fv_rst(ifile)%vars(ivar)%nLev==npz+1) then
                 do k=1,in_fv_rst(ifile)%vars(ivar)%nLev
                    q_edge_old(is:ie,k)=in_fv_rst(ifile)%vars(ivar)%ptr3d(is:ie,j,k)
                 enddo
                 call remap_edge(q_edge_old,q_edge_new,is,ie,km,npz,Atm%ak,Atm%bk)
                 do k=1,npz+1
                    do i=is,ie
                       out_fv_rst(ifile)%vars(ivar)%ptr3d(i,j,k) = q_edge_new(i,k)
                    enddo
                 enddo
              else
                 out_fv_rst(ifile)%vars(ivar)%ptr2d(is:ie,j)=in_fv_rst(ifile)%vars(ivar)%ptr2d(is:ie,j)
              end if
           enddo
        else
           do ivar=1,size(out_fv_rst(ifile)%vars)
              out_fv_rst(ifile)%vars(ivar)%ptr3d(is:ie,j,:)=in_fv_rst(ifile)%vars(ivar)%ptr3d(is:ie,j,:)
           end do
        end if
     enddo

!-------------------------------------------------------------
! map virtual temperature using geopotential conserving scheme.
!-------------------------------------------------------------
     call mappm(km, pn0, tp, npz, pn1, qn1, is,ie, 1, 9, Atm%ptop)
     do k=1,npz
        do i=is,ie
           Atm%pt(i,j,k) = qn1(i,k)/(1.+zvir*Atm%q(i,j,k,sphum))
        enddo
     enddo

  enddo

  call prt_maxmin('PS_model', Atm%ps, is, ie, js, je, ng, 1, 0.01_FVPRC)

  if (is_master()) write(*,*) 'done remap_scalar'

end subroutine remap_scalar

subroutine remap_edge(q1,q2,is,ie,km,kn,ak,bk)

!  q1,km - old levels
!  q2,kn - new levels
   integer, intent(in) :: is,ie,km,kn
   real(FVPRC),intent(in) :: ak(kn+1), bk(kn+1)
   real(FVPRC),intent(in) :: q1(is:ie,km+1)
   real(FVPRC),intent(out) :: q2(is:ie,kn+1)

   integer i,k
   do i=is,ie
      do k=1,kn+1
         q2(i,k)=ak(k)+bk(k)*q1(i,km+1)
      enddo
   enddo

end subroutine remap_edge

function simple_cs_grid_creator(IM_World,JM_World,NX,NY,rc) result(esmfgrid)

   integer,           intent(IN)    :: IM_WORLD, JM_WORLD
   integer,           intent(IN)    :: NX, NY
   integer, optional, intent(OUT)   :: rc
   type (ESMF_Grid)                 :: esmfgrid

   integer, allocatable             :: IMS(:), JMS(:)
   integer                          :: n
   integer                          :: STATUS
   character(len=ESMF_MAXSTR), parameter :: Iam="simple_cs_grid_creator"

   allocate( IMS(0:NX-1) )
   allocate( JMS(0:NY-1) )

   call MAPL_DecomposeDim ( IM_WORLD  , IMS             , NX  , symmetric=.true. )
   call MAPL_DecomposeDim ( JM_WORLD/6, JMS(0:NY/6 -1)  , NY/6, symmetric=.true. )
   do n=2,6
      JMS((n-1)*NY/6 : n*NY/6 -1) = JMS(0:NY/6 -1)
   enddo

   esmfgrid = ESMF_GridCreate(             &
        name="dummy",                  &
        countsPerDEDim1=IMS,           &
        countsPerDEDim2=JMS,           &
        indexFlag = ESMF_INDEX_USER,   &
        gridMemLBound = (/1,1/),       &
        gridEdgeLWidth = (/0,0/),      &
        gridEdgeUWidth = (/0,0/),      &
        coordDep1 = (/1,2/),           &
        coordDep2 = (/1,2/),           &
        rc=status)
   VERIFY_(STATUS)

   RETURN_(ESMF_SUCCESS)

end function simple_cs_grid_creator

end module fv_regridding_utils

