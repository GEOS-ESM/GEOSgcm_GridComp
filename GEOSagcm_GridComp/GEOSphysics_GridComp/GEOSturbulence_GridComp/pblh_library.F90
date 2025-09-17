module GEOSturbulence_PBLH_Library

  use ESMF
  use MAPL
  use MAPL_ConstantsMod, only: mapl_grav, mapl_cp, mapl_karman, mapl_rdry, mapl_pi, mapl_p00
  use MAPL_SatVaporMod, only: MAPL_EQsat

  implicit none
  private

  character(len=ESMF_MAXSTR)              :: IAm="GEOSturbulence_PBLH_Library"

  public :: find_bulk_ri_pblh
  public :: find_kh2_pblh
  public :: find_kh10p_pblh
  public :: find_ri_pblh
  public :: find_thv_pblh
  public :: find_rfrct_pblh
  public :: find_qv_pblh
  public :: return_surface_inversion_stats
  public :: return_trade_inversion_stats
  
contains

   !----------------------------------------------------------------------------
   ! Bulk-Richardson-Number PBLH (aka Transcom), described by Seidel et al 2012.
   ! Scans from surface to find level at which bulk Ri first exceeds a critical
   ! value specified by bulkri_crit parameter.
   !----------------------------------------------------------------------------
   subroutine find_bulk_ri_pblh(im,jm,lm,u,v,z,thv,zpbl_bulkri,kpbl_bulkri,bulkri_crit)

      integer,                      intent(in)  :: im, jm, lm
      real                        , intent(in)  :: bulkri_crit
      real   , dimension(im,jm,lm), intent(in)  :: u, v, z, thv
      real   , dimension(im,jm)   , intent(out) :: zpbl_bulkri
      real   , dimension(im,jm)   , intent(out) :: kpbl_bulkri

      ! Local vars
      integer             :: i, j, k
      real                :: uv2h
      real, dimension(lm) :: tcrib
     
      zpbl_bulkri = MAPL_UNDEF
      tcrib(lm) = 0.0
      do i = 1, im
         do j = 1, jm
            do k=lm-1,1,-1
               uv2h = max(u(i,j,k)**2+v(i,j,k)**2,1.0E-8)
               tcrib(k) = MAPL_GRAV * (thv(i,j,k)-thv(i,j,LM)) * z(i,j,k) / (thv(i,j,LM)*uv2h)
               if (tcrib(k) >= bulkri_crit) then
                  zpbl_bulkri(i,j) = z(i,j,k+1)+(bulkri_crit-tcrib(k+1))/(tcrib(k)-tcrib(k+1))*(z(i,j,k)-z(i,j,k+1))
                  kpbl_bulkri(i,j) = float(k)
                  exit
               end if
            end do
         end do
      end do
      where (zpbl_bulkri<0.)
         zpbl_bulkri = z(:,:,lm)
         kpbl_bulkri = float(lm)
      end where
    
   end subroutine find_bulk_ri_pblh
  

   !----------------------------------------------------------------------------
   ! Diffusivity-based PBLH with threshold of 2 m2/s2
   ! Scans from surface to find level at which KH drops below threshold.
   !----------------------------------------------------------------------------
   subroutine find_kh2_pblh(im,jm,lm,kpblmin,z,kh,zpbl_kh2,kpbl_kh2)

      integer,                        intent(in)  :: im, jm, lm, kpblmin
      real   , dimension(im,jm,lm)  , intent(in)  :: z  ! surface-relative heights
      real   , dimension(im,jm,0:lm), intent(in)  :: kh
      real   , dimension(im,jm)     , intent(out) :: zpbl_kh2
      real   , dimension(im,jm)     , intent(out) :: kpbl_kh2

      ! Local vars
      integer :: i, j, k

      zpbl_kh2 = MAPL_UNDEF

      do i = 1,im
         do j = 1,jm
            do k=lm,2,-1
               if ((kh(i,j,k-1) < 2.).and.(kh(i,j,k) >= 2.).and.(zpbl_kh2(i,j)==MAPL_UNDEF)) then
                  zpbl_kh2(i,j) = z(i,j,k)
                  kpbl_kh2(i,j) = float(k)
               end if
            end do
         end do
      end do

      where ( zpbl_kh2 .eq. MAPL_UNDEF )
         zpbl_kh2 = z(:,:,lm)
         kpbl_kh2 = float(lm)
      end where
      zpbl_kh2 = min(zpbl_kh2,z(:,:,kpblmin))

   end subroutine find_kh2_pblh
  

   !----------------------------------------------------------------------------
   ! Diffusivity-based PBLH with threshold of 10 percent of maximum
   ! Scans from surface to find level at which KH drops below threshold.
   !----------------------------------------------------------------------------
   subroutine find_kh10p_pblh(im,jm,lm,kpblmin,z,ze,kh,zpbl_kh10p,kpbl_kh10p)

      integer,                        intent(in)  :: im, jm, lm, kpblmin
      real   , dimension(im,jm,lm)  , intent(in)  :: z  ! surface-relative heights
      real   , dimension(im,jm,0:lm), intent(in)  :: ze  ! surface-relative edge heights
      real   , dimension(im,jm,0:lm), intent(in)  :: kh
      real   , dimension(im,jm)     , intent(out) :: zpbl_kh10p
      real   , dimension(im,jm)     , intent(out) :: kpbl_kh10p

      ! Local vars
      real, dimension(lm+1) :: temparray
      integer :: i, j, k, locmax
      real    :: minlval, maxkh

      zpbl_kh10p = MAPL_UNDEF

      do i = 1,im
         do j = 1,jm

            temparray(1:lm+1) = kh(I,J,0:lm)
            do k = lm,2,-1
               locmax = maxloc(temparray,1)
               minlval = max(0.001,0.0001*maxval(temparray))
               if(temparray(locmax-1)<minlval.and.temparray(locmax+1)<minlval) temparray(locmax) = minlval
            enddo
            maxkh = temparray(lm)
            do k = lm-1,2,-1
               if(temparray(k)>maxkh) maxkh = temparray(k)
               if(temparray(k-1)<minlval) exit
            end do
            do k=lm-1,2,-1
               if ( (temparray(k) < 0.1*maxkh) .and. (temparray(k+1) >= 0.1*maxkh)  &
                  .and. (zpbl_kh10p(i,j) == MAPL_UNDEF ) ) then
                  zpbl_kh10p(i,j) = ze(i,j,k)+ &
                  ((ze(i,j,k-1)-ze(i,j,k))/(temparray(k)-temparray(k+1))) * (0.1*maxkh-temparray(k+1))
                     kpbl_kh10p(i,j) = float(k)
               end if
            end do
            if (  zpbl_kh10p(i,j) .eq. MAPL_UNDEF .or. (maxkh.lt.1.)) then
               zpbl_kh10p(i,j) = z(i,j,lm)
               kpbl_kh10p(i,j) = float(lm)
            endif
         end do
      end do
      
      zpbl_kh10p = min(zpbl_kh10p,ze(:,:,kpblmin))

   end subroutine find_kh10p_pblh
  

   !----------------------------------------------------------------------------
   ! Local Richardson number PBLH with a threshold of ri_crit.
   ! Scans from surface to find level at which RI drops below threshold.
   !----------------------------------------------------------------------------
   subroutine find_ri_pblh(im,jm,lm,kpblmin,z,ri,zpbl_ri,ri_crit)

      integer,                        intent(in)  :: im, jm, lm, kpblmin
      real                          , intent(in)  :: ri_crit
      real   , dimension(im,jm,lm)  , intent(in)  :: z  ! surface-relative heights
      real   , dimension(im,jm,0:lm), intent(in)  :: ri ! local Richardson number
      real   , dimension(im,jm)     , intent(out) :: zpbl_ri

      ! local variables
      integer :: i, j, k

      zpbl_ri = MAPL_UNDEF 
      where (ri(:,:,lm-1)>ri_crit) zpbl_ri = z(:,:,lm)

      do i = 1, im
         do j = 1, jm
            do k=lm-1,1,-1
               if( (ri(i,j,k-1)>ri_crit) .and. (zpbl_ri(i,j) == MAPL_UNDEF) ) then
                  zpbl_ri(i,j) = z(i,j,k+1)+(ri_crit-ri(i,j,k))/(ri(i,j,k-1)-ri(i,j,k))*(z(i,j,k)-z(i,j,k+1))
               end if
            end do
         end do
      end do

      where ( zpbl_ri .eq. MAPL_UNDEF ) zpbl_ri = z(:,:,lm)
      zpbl_ri = min(zpbl_ri,z(:,:,kpblmin))
      where ( zpbl_ri < 0.0 ) zpbl_ri = z(:,:,lm)

   end subroutine find_ri_pblh

   
   !-----------------------------------------------------------------------------
   ! Potential temperature gradient PBLH.
   ! Scans from surface to find height of maximum virtual potential temperature
   ! gradient.
   !-----------------------------------------------------------------------------
   subroutine find_thv_pblh(im,jm,lm,kpblmin,z,thv,zpbl_thv)

      integer,                      intent(in)  :: im, jm, lm, kpblmin
      real   , dimension(im,jm,lm), intent(in)  :: z  ! surface-relative heights
      real   , dimension(im,jm,lm), intent(in)  :: thv ! virtual potential temperature
      real   , dimension(im,jm)   , intent(out) :: zpbl_thv

      ! local variables
      integer :: i, j, k
      real    :: dthvdz, maxdthvdz
      
      zpbl_thv = MAPL_UNDEF

      do i = 1, im
         do j = 1, jm

            maxdthvdz = 0

            do k=lm-1,1,-1
               if(z(i,j,k)<=z(i,j,kpblmin)) then
                  dthvdz = (thv(i,j,k+1)-thv(i,j,k))/(z(i,j,k+1)-z(i,j,k))
                  if(dthvdz>maxdthvdz) then
                     maxdthvdz = dthvdz
                     zpbl_thv(i,j) = 0.5*(z(i,j,k+1)+z(i,j,k))
                  end if
               end if
            end do

         end do
      end do      

   end subroutine find_thv_pblh


   !-----------------------------------------------------------------------------
   ! Specific humidity gradient PBLH.
   ! Scans from surface to find height of maximum gradient in specific humidity.
   !-----------------------------------------------------------------------------
   subroutine find_qv_pblh(im,jm,lm,kpblmin,z,qv,zpbl_qv)

      integer,                      intent(in)  :: im, jm, lm, kpblmin
      real   , dimension(im,jm,lm), intent(in)  :: z  ! surface-relative heights
      real   , dimension(im,jm,lm), intent(in)  :: qv ! specific humidity
      real   , dimension(im,jm)   , intent(out) :: zpbl_qv

      ! local variables
      integer :: i, j, k
      real    :: dqvdz, maxdqvdz

      zpbl_qv = MAPL_UNDEF

      do i = 1, im
         do j = 1, jm

            maxdqvdz = 0.
            do k=lm-1,1,-1
               if(z(i,j,k)<=z(i,j,kpblmin)) then
                  dqvdz = -1.*(qv(i,j,k+1)-qv(i,j,k))/(z(i,j,k+1)-z(i,j,k))
                  if(dqvdz>maxdqvdz) then
                     maxdqvdz = dqvdz
                     zpbl_qv(i,j) = 0.5*(z(i,j,k+1)+z(i,j,k))
                  end if
               end if
            end do

         end do
      end do

   end subroutine find_qv_pblh


   !=========================================================================
   !  ZPBL defined by minimum in vertical gradient of refractivity.
   !  As shown in Ao, et al, 2012: "Planetary boundary layer heights from
   !  GPS radio occultation refractivity and humidity profiles", Climate and
   !  Dynamics.  https://doi.org/10.1029/2012JD017598
   !=========================================================================
   subroutine find_rfrct_pblh(im,jm,lm,z,p,t,qv,zpbl_rfrct)

      integer,                      intent(in)  :: im, jm, lm
      real   , dimension(im,jm,lm), intent(in)  :: z, &  ! surface-relative heights
                                                   p, &
                                                   t, &
                                                   qv
      real   , dimension(im,jm)   , intent(out) :: zpbl_rfrct
      
      ! local variables
      integer :: i, j, k
      real    :: a1, a2
      real, dimension(im,jm,lm) :: wvp, dum3d, tmp3d

      a1 = 0.776    ! K/Pa
      a2 = 3.73e3   ! K2/Pa
      
      wvp = qv * p / (qv*(1.-0.622)+0.622)  ! water vapor partial pressure
      
      ! pressure gradient term
      dum3d(:,:,2:lm-1) = (p(:,:,1:lm-2)-p(:,:,3:lm)) / (z(:,:,1:lm-2)-z(:,:,3:lm))
      dum3d(:,:,1) = (p(:,:,1)-p(:,:,2)) / (z(:,:,1)-z(:,:,2))
      dum3d(:,:,lm) = (p(:,:,lm-1)-p(:,:,lm)) / (z(:,:,lm-1)-z(:,:,lm))
      tmp3d = a1 * dum3d / t

      ! add temperature gradient term
      dum3d(:,:,2:lm-1) = (t(:,:,1:lm-2)-t(:,:,3:lm)) / (z(:,:,1:lm-2)-z(:,:,3:lm))
      dum3d(:,:,1) = (t(:,:,1)-t(:,:,2)) / (z(:,:,1)-z(:,:,2))
      dum3d(:,:,lm) = (t(:,:,lm-1)-t(:,:,lm)) / (z(:,:,lm-1)-z(:,:,lm))
      tmp3d = tmp3d - (a1*p/t**2 + 2.*a2*wvp/t**3)*dum3d

      ! add vapor pressure gradient term
      dum3d(:,:,2:lm-1) = (wvp(:,:,1:lm-2)-wvp(:,:,3:lm)) / (z(:,:,1:lm-2)-z(:,:,3:lm))
      dum3d(:,:,1) = (wvp(:,:,1)-wvp(:,:,2)) / (z(:,:,1)-z(:,:,2))
      dum3d(:,:,lm) = (wvp(:,:,lm-1)-wvp(:,:,lm)) / (z(:,:,lm-1)-z(:,:,lm))
      tmp3d = tmp3d + (a2/t**2)*dum3d

      ! zpbl is height of minimum in refractivity (tmp3d)
      do i = 1,im
        do j = 1,jm
           k = minloc(tmp3d(i,j,:),dim=1,back=.true.)   ! return last index, if multiple
           zpbl_rfrct(i,j) = z(i,j,k)
        end do
      end do   

   end subroutine find_rfrct_pblh


   !-----------------------------------------------------------------------------
   ! Surface-based inversion (top) height and frequency
   ! Scans up from surface while absolute temperature increases with height.
   !-----------------------------------------------------------------------------
   subroutine return_surface_inversion_stats(im,jm,lm,z,t,sbitop,sbifrq)

      integer,                      intent(in)  :: im, jm, lm
      real   , dimension(im,jm,lm), intent(in)  :: z  ! surface-relative heights
      real   , dimension(im,jm,lm), intent(in)  :: t  ! temperature
      real   , dimension(im,jm)   , intent(out) :: sbitop  ! sbi top height
      real   , dimension(:,:)     , pointer     :: sbifrq  ! sbi frequency   
      
      ! local variables
      integer :: i, j, k

      sbitop = MAPL_UNDEF

      do i = 1, im
         do j = 1, jm
            if (t(i,j,lm-1).gt.t(i,j,lm)) then
               do k = lm-1,1,-1
                  if (t(i,j,k).gt.t(i,j,k+1)) then
                     sbitop(i,j) = z(i,j,k)
                  else
                     exit
                  end if
               end do
            end if
         end do
      end do
      if (associated(sbifrq)) then
         sbifrq = 0.
         where(sbitop.ne.MAPL_UNDEF) sbifrq=1.
      end if

   end subroutine return_surface_inversion_stats


   !-----------------------------------------------------------------------------
   ! Trade inversion height, magnitude and frequency
   ! Searches upward from 950mb to 600mb for either an absolute temperature
   ! increase of at least 0.5 K, or a segment of at least 25mb with no
   ! temperature decrease.
   !-----------------------------------------------------------------------------
   subroutine return_trade_inversion_stats(im,jm,lm,p,t,trinvbs,trinvdelt,trinvfrq)

      integer,                      intent(in)  :: im, jm, lm
      real   , dimension(im,jm,lm), intent(in)  :: p  ! pressure
      real   , dimension(im,jm,lm), intent(in)  :: t  ! temperature
      real   , dimension(im,jm)   , intent(out) :: trinvbs   ! inversion base height
      real   , dimension(:,:)     , pointer     :: trinvdelt ! inversion temperature jump
      real   , dimension(:,:)     , pointer     :: trinvfrq  ! inversion frequency
      
      ! local variables
      integer :: i, j, k, ktop
      real, dimension(im,jm) :: tmpdelt, tmpfrq
      
      trinvbs = MAPL_UNDEF
      tmpdelt = MAPL_UNDEF
      tmpfrq = 0.
      do i = 1,im
         do j = 1,jm
            k = lm

            do while (p(i,j,k).gt.95000.)
               k = k-1
            end do
            do k = k,1,-1    ! k is first level above 950mb
               if (p(i,j,k).lt.60000.) exit
               if (t(i,j,k-1).ge.t(i,j,k)) then ! if next level is warmer...
                  ktop = k                      ! k is index of minimum t so far
                  do while (t(i,j,ktop).ge.t(i,j,k)) ! find depth of warm layer                                                                                                         
                     ktop = ktop-1
                  end do
                  ktop = ktop+1   ! ktop is index of highest level inside warm layer
                    
                  if (  maxval(t(i,j,ktop:k))-t(i,j,k).ge.0.5 .or. &
                       (maxval(t(i,j,ktop:k))-t(i,j,k).gt.0.01 .and. p(i,j,k)-p(i,j,ktop)>2500.) ) then

                     ! only save if delta-t exceeds any previous inversion
                     if ( tmpfrq(i,j).eq.0. .or. &
                         (tmpfrq(i,j).ne.0. .and. maxval(t(i,j,ktop:k))-t(i,j,k).gt.tmpdelt(i,j)) ) then
                        trinvbs(i,j) = p(i,j,k)
                        tmpdelt(i,j) = maxval(t(i,j,ktop:k))-t(i,j,k)
                        tmpfrq(i,j)  = 1.
                     end if

                  end if
               end if ! next level warmer
                 
	    end do ! k vert loop
         end do ! j
      end do ! i
      if (associated(trinvfrq))  trinvfrq  = tmpfrq
      if (associated(trinvdelt)) trinvdelt = tmpdelt
      
   end subroutine return_trade_inversion_stats
   
end module GEOSturbulence_PBLH_Library

  
