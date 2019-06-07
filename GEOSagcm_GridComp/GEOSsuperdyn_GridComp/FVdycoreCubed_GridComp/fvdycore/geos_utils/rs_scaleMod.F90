#include "MAPL_Generic.h"
module rs_scaleMod

   use fv_arrays_mod
   use MAPL_ConstantsMod, only: MAPL_PSDRY
   use MAPL_Mod
   use ESMF
   use pFIO_StringIntegerMapMod 
   use, intrinsic :: iso_fortran_env, only: REAL64, REAL32 
   ! bma added
   implicit none

! ************************************************************************
! ************************************************************************
! ****                                                                ****
! ****   Program to insert DRY MASS value of 983.05 mb into Restarts  ****
! ****   (to within 1e-10 Pa)                                         ****
! ****                                                                ****
! ************************************************************************
! ************************************************************************

   private

   public scale_drymass

contains


   subroutine scale_drymass(FVAtm,moist_tracers,rc)
      type(FV_atmos_type), intent(inout) :: FVAtm(:)
      type(StringIntegerMap), intent(inout) :: moist_tracers
      integer, optional,   intent(out  ) :: rc

      integer i,j,L,n,im,jm,lm
      integer iter

! restart variables
! -----------------
      real(REAL64),   allocatable ::   pe(:,:,:)
      real(REAL64),   allocatable ::   pk(:,:,:)
      real(REAL64),   allocatable ::   dp(:,:,:)
      real(REAL64),   allocatable ::  pke(:,:,:)

      real(REAL32),   allocatable ::   qv(:,:,:)
      real(REAL32),   allocatable :: qlls(:,:,:)
      real(REAL32),   allocatable :: qlcn(:,:,:)
      real(REAL32),   allocatable :: cfls(:,:,:)
      real(REAL32),   allocatable :: cfcn(:,:,:)
      real(REAL32),   allocatable :: qils(:,:,:)
      real(REAL32),   allocatable :: qicn(:,:,:)
      real(REAL32),   allocatable :: area(:,:)
      integer iq,iqlls,iqlcn,iqils,iqicn,icfls,icfcn

      real(REAL64), allocatable ::    gsum(:,:)
      real(REAL32), allocatable ::    qsum(:,:)
      real(REAL64), allocatable ::   psold(:,:)
      real(REAL64), allocatable ::   psnew(:,:)
      real(REAL32), allocatable :: pdryold(:,:)
      real(REAL32), allocatable :: pdrynew(:,:)

      real(REAL64), parameter   ::    pdry_ave = MAPL_PSDRY
      real(REAL64)              :: pdryold_ave
      real(REAL64)              :: pdrynew_ave
      real(REAL64)              :: pdrydif_ave
      real(REAL64), parameter   :: eps = epsilon(1.0d-10)
      real                :: kappa
      integer             :: ie,is,je,js,status
      character(len=ESMF_MAXSTR), parameter :: Iam="scale_drymass"

      kappa = 2.0/7.0

      is = FVAtm(1)%bd%isc
      ie = FVAtm(1)%bd%iec
      js = FVAtm(1)%bd%jsc
      je = FVAtm(1)%bd%jec
      im = ie-is+1
      jm = je-js+1
      lm = FVatm(1)%flagstruct%npz

! **********************************************************************


      allocate (   pk(im,jm,lm)   )
      allocate (   pe(im,jm,lm+1)   )
      allocate (  pke(im,jm,lm+1) )

      allocate ( dp(im,jm,lm) )
      do L=1,lm+1
         pe(:,:,L) = FVatm(1)%pe(is:ie,L,js:je)
      enddo
      pk = fvatm(1)%pkz(is:ie,js:je,:)

      do L=1,lm
         dp(:,:,L) = pe(:,:,L+1)-pe(:,:,L)
      enddo

! **********************************************************************
! ****                   Read moist internal Restart                ****
! **********************************************************************

      allocate (   qv(im,jm,lm) )
      allocate ( qlls(im,jm,lm) )
      allocate ( qlcn(im,jm,lm) )
      allocate ( cfls(im,jm,lm) )
      allocate ( cfcn(im,jm,lm) )
      allocate ( qils(im,jm,lm) )
      allocate ( qicn(im,jm,lm) )

      iq = moist_tracers%at("Q")
      iqlls = moist_tracers%at("QLLS")
      iqlcn = moist_tracers%at("QLCN")
      iqicn = moist_tracers%at("QICN")
      iqils = moist_tracers%at("QILS")
      icfls = moist_tracers%at("CLLS")
      icfcn = moist_tracers%at("CLCN")
      qv = fvatm(1)%q(is:ie,js:je,:,iq)
      qlls = fvatm(1)%q(is:ie,js:je,:,iqlls)
      qlcn = fvatm(1)%q(is:ie,js:je,:,iqlcn)
      cfls = fvatm(1)%q(is:ie,js:je,:,icfls)
      cfcn = fvatm(1)%q(is:ie,js:je,:,icfcn)
      qils = fvatm(1)%q(is:ie,js:je,:,iqils)
      qicn = fvatm(1)%q(is:ie,js:je,:,iqicn)


! **********************************************************************
! ****                 Compute/Import Grid-Cell Area                ****
! **********************************************************************

      allocate ( area(im,jm) )
      area = FVatm(1)%gridstruct%area(is:ie,js:je)

! **********************************************************************
! ****                       Constrain PDRY                         ****
! **********************************************************************

      allocate ( gsum(5,2)   )
      allocate ( qsum(im,jm) )

      allocate (   psold(im,jm) )
      allocate (   psnew(im,jm) )
      allocate ( pdryold(im,jm) )
      allocate ( pdrynew(im,jm) )

      pdrydif_ave = 1.0d0
      iter = 1

      do while ( dabs( pdrydif_ave ).gt.eps .and. iter.le.20 )

         do n=1,5
            qsum = 0.0_8
            do L=1,lm
               if( n.eq.1 ) qsum = qsum +   qv(:,:,L)*dp(:,:,L)
               if( n.eq.2 ) qsum = qsum + qlls(:,:,L)*dp(:,:,L)
               if( n.eq.3 ) qsum = qsum + qlcn(:,:,L)*dp(:,:,L)
               if( n.eq.4 ) qsum = qsum + qils(:,:,L)*dp(:,:,L)
               if( n.eq.5 ) qsum = qsum + qicn(:,:,L)*dp(:,:,L)
            enddo
            call AreaMean( qsum, area, gsum(n,1), rc=status )
         enddo

         qsum = 0.0_8
         do L=1,lm
            qsum = qsum + (  qv(:,:,L) + &
                  qlls(:,:,L) + &
                  qlcn(:,:,L) + &
                  qils(:,:,L) + &
                  qicn(:,:,L) ) * dp(:,:,L)
         enddo
         psold = pe(:,:,lm+1)
         pdryold = pe(:,:,lm+1) - qsum  ! Subtract Total Water Content

         call AreaMean( pdryold, area, pdryold_ave, rc=status )

         pdrynew = pdryold * ( pdry_ave/pdryold_ave )

         psnew = pdrynew + qsum

         do L=1,lm+1
            do j=1,jm
               do i=1,im
                  pe(i,j,L) = pe(i,j,L) + FVatm(1)%bk(L)*( psnew(i,j)-psold(i,j) )
               enddo
            enddo
         enddo

         pke = pe**kappa
         do L=1,lm
            dp(:,:,L) =             pe(:,:,L+1)-pe(:,:,L)
            pk(:,:,L) =            (pke(:,:,L+1)-pke(:,:,L)) &
                  / ( kappa*log(pe(:,:,L+1)/pe(:,:,L)) )
         enddo

! --------------------------------------

         do n=1,5
            qsum = 0.0_8
            do L=1,lm
               if( n.eq.1 ) qsum = qsum +   qv(:,:,L)*dp(:,:,L)
               if( n.eq.2 ) qsum = qsum + qlls(:,:,L)*dp(:,:,L)
               if( n.eq.3 ) qsum = qsum + qlcn(:,:,L)*dp(:,:,L)
               if( n.eq.4 ) qsum = qsum + qils(:,:,L)*dp(:,:,L)
               if( n.eq.5 ) qsum = qsum + qicn(:,:,L)*dp(:,:,L)
            enddo
            call AreaMean( qsum, area, gsum(n,2), rc=status )
         enddo

         qv =   qv * ( gsum(1,1)/gsum(1,2) )
         qlls = qlls * ( gsum(2,1)/gsum(2,2) )
         qlcn = qlcn * ( gsum(3,1)/gsum(3,2) )
         qils = qils * ( gsum(4,1)/gsum(4,2) )
         qicn = qicn * ( gsum(5,1)/gsum(5,2) )

         qsum = 0.0_8
         do L=1,lm
            qsum = qsum + (  qv(:,:,L) + &
                  qlls(:,:,L) + &
                  qlcn(:,:,L) + &
                  qils(:,:,L) + &
                  qicn(:,:,L) ) * dp(:,:,L)
         enddo
         pdrynew = pe(:,:,lm+1) - qsum  ! Subtract Total Water Content

         call AreaMean( pdrynew, area, pdrynew_ave, rc=status )

         pdrydif_ave = pdrynew_ave - pdryold_ave

         if (mapl_am_I_root()) write(6,1001) pdrynew_ave/100,pdryold_ave/100,pdrynew_ave/pdryold_ave,pdrydif_ave/100
1001     format(1x,'PSDRY_NEW: ',g21.14,'  PSDRY_OLD: ',g21.14,'  RATIO: ',g25.18,'  DIF: ',g21.14)

         iter = iter + 1
      enddo

      fvatm(1)%q(is:ie,js:je,:,iq) = qv
      fvatm(1)%q(is:ie,js:je,:,iqlls) = qlls
      fvatm(1)%q(is:ie,js:je,:,iqlcn) = qlcn
      fvatm(1)%q(is:ie,js:je,:,icfls) = cfls
      fvatm(1)%q(is:ie,js:je,:,icfcn) = cfcn
      fvatm(1)%q(is:ie,js:je,:,iqils) = qils
      fvatm(1)%q(is:ie,js:je,:,iqicn) = qicn

      do L=1,lm+1
         FVatm(1)%pe(is:ie,L,js:je) = pe(:,:,L)
      enddo
      fvatm(1)%pkz(is:ie,js:je,:) = pk

      deallocate(qv,qlls,qlcn,cfls,cfcn,qils,qicn,pk,pe,pke,area,dp,psold,psnew,pdryold,pdrynew)

      RETURN_(ESMF_SUCCESS)

   end subroutine scale_drymass

  subroutine AreaMean ( q, area, qave, rc )

    implicit none

    ! Arguments
    real(REAL64),            intent(  OUT) :: qave
    real(REAL32),              intent(IN   ) :: q(:,:)
    real(REAL32),              intent(IN   ) :: area(:,:)
    integer, optional, intent(  OUT) :: rc

    ! Log err vars
    integer :: status
    character(len=ESMF_MAXSTR), parameter :: Iam='AreaMean'

    ! Local vars
    real(REAL64)  :: qdum(2)
    real(REAL64)  :: qdumloc(2)
    integer :: im,jm

    integer :: i,j

    type(ESMF_VM) :: vm


    ! get VM (should get from the grid, but this is quicker)
    call ESMF_VmGetCurrent(vm, rc=status)
    VERIFY_(STATUS)

    im = size(area,1) ! local grid dim
    jm = size(area,2) ! local grid dim

    ! do calculation on every PE

    ! compute local sum
    qdumloc = 0.0_8
    do j=1,jm
       do i=1,im
          if (q(i,j) == MAPL_Undef) cycle ! exclude any undefs
          qdumloc(1) = qdumloc(1) + q(i,j)*area(i,j)
          qdumloc(2) = qdumloc(2) + area(i,j)
       enddo
    end do

    call MAPL_CommsAllReduceSum(vm, sendbuf=qdumloc, recvbuf=qdum, &
         cnt=2, RC=status)
    VERIFY_(STATUS)

    if (qdum(2) /= 0.0_8) then

       qave = qdum(1) / qdum(2)

       !ALT: convert the the result to single precision to get rid of
       !     numerical non-associativity in floating point numbers
       !     qave = real(qave, kind=4)
    else
       qave = MAPL_Undef
    end if

    RETURN_(ESMF_SUCCESS)
  end subroutine AreaMean

end module rs_scaleMod
