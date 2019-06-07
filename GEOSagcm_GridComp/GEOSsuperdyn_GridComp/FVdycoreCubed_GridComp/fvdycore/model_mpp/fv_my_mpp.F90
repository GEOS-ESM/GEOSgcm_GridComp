module fv_my_mpp

use fv_mp_mod, only : domain
use mpp_domains_mod,  only: CGRID_NE, DGRID_NE, mpp_get_boundary,   &
                                mpp_update_domains

  public mpp_update_domains_dummy,mpp_get_boundary_dummy,mp_reduce_max_dummy
  public mpp_update_domains_dummy4,mpp_update_domains_CGRID,mpp_update_domains_DGRID
  public mpp_get_boundary_DGRID

  interface mpp_update_domains_dummy
     module procedure mpp_update_domains_dummy1
     module procedure mpp_update_domains_dummy2
     module procedure mpp_update_domains_dummy3
  end interface

  interface mpp_get_boundary_dummy
     module procedure mpp_get_boundary_dummy1
     module procedure mpp_get_boundary_dummy2
  end interface

  interface mp_reduce_max_dummy
     module procedure mp_reduce_max_dummy1
     module procedure mp_reduce_max_dummy2
  end interface

CONTAINS

 subroutine mp_reduce_max_dummy1(cmax)
 real :: cmax
! cmax=1.*cmax
 end subroutine mp_reduce_max_dummy1
 
 subroutine mp_reduce_max_dummy2(cmax,npz)
 real :: cmax(npz)
! cmax(:)=1.*cmax(:)
 end subroutine mp_reduce_max_dummy2

 subroutine mpp_update_domains_dummy1(field,is,ie,js,je,isd,ied,jsd,jed,km)
  integer :: is,ie,js,je,isd,ied,jsd,jed,km
  real :: field(isd:ied,jsd:jed,km)
  real :: field_temp

  integer :: i,j,k

      call mpp_update_domains(  field, domain, complete=.true.)

!      do k=1,km
!         field_temp=field(is,js,k)
!         field(isd:is-1,:,k)=field_temp
!         field(:,jsd:js-1,k)=field_temp

!         field_temp=field(ie,je,k)
!         field(ie+1:ied,:,k)=field_temp
!         field(:,je+1:jed,k)=field_temp
!      enddo

 end subroutine mpp_update_domains_dummy1

 subroutine mpp_update_domains_dummy2(field1,field2,&
                                      is1,ie1,js1,je1,isd1,ied1,jsd1,jed1,km1,&
                                      is2,ie2,js2,je2,isd2,ied2,jsd2,jed2,km2)!,my_flag)

  integer :: is1,ie1,js1,je1,isd1,ied1,jsd1,jed1,km1
  integer :: is2,ie2,js2,je2,isd2,ied2,jsd2,jed2,km2
  real :: field1(isd1:ied1,jsd1:jed1,km1)
  real :: field2(isd2:ied2,jsd2:jed2,km2)
  real :: field_temp

  integer :: i,j,k,my_flag

!  if (my_flag.eq.1) call mpp_update_domains(field2,field1, domain, gridtype=CGRID_NE, complete=.true.)
!  if (my_flag.eq.2) call mpp_update_domains(field1,field2, domain, gridtype=DGRID_NE)

      do k=1,km2
         field_temp=field2(is2,js2,k)
         field2(isd2:is2-1,:,k)=field_temp
         field2(:,jsd2:js2-1,k)=field_temp

         field_temp=field2(ie2,je2,k)
         field2(ie2+1:ied2-1,:,k)=field_temp
         field2(:,je2+1:jed2,k)=field_temp
     enddo

     do k=1,km1
         field_temp=field1(is1,js1,k)
         field1(isd1:is1-1,:,k)=field_temp
         field1(:,jsd1:js1-1,k)=field_temp

         field_temp=field1(ie1,je1,k)
         field1(ie1+1:ied1,:,k)=field_temp
         field1(:,je1+1:jed1-1,k)=field_temp
      enddo

 end subroutine mpp_update_domains_dummy2

 subroutine mpp_update_domains_CGRID(field1,field2,&
                                      is1,ie1,js1,je1,isd1,ied1,jsd1,jed1,km1,&
                                      is2,ie2,js2,je2,isd2,ied2,jsd2,jed2,km2)

  integer :: is1,ie1,js1,je1,isd1,ied1,jsd1,jed1,km1
  integer :: is2,ie2,js2,je2,isd2,ied2,jsd2,jed2,km2
  real :: field1(isd1:ied1,jsd1:jed1,km1)
  real :: field2(isd2:ied2,jsd2:jed2,km2)
  real :: field_temp

  integer :: i,j,k,my_flag

  call mpp_update_domains(field2,field1, domain, gridtype=CGRID_NE, complete=.true.)

 end subroutine mpp_update_domains_CGRID

 subroutine mpp_update_domains_DGRID(field1,field2,&
                                      is1,ie1,js1,je1,isd1,ied1,jsd1,jed1,km1,&
                                      is2,ie2,js2,je2,isd2,ied2,jsd2,jed2,km2)

  integer :: is1,ie1,js1,je1,isd1,ied1,jsd1,jed1,km1
  integer :: is2,ie2,js2,je2,isd2,ied2,jsd2,jed2,km2
  real :: field1(isd1:ied1,jsd1:jed1,km1)
  real :: field2(isd2:ied2,jsd2:jed2,km2)
  real :: field_temp

  integer :: i,j,k,my_flag

  call mpp_update_domains(field1,field2, domain, gridtype=DGRID_NE)

 end subroutine mpp_update_domains_DGRID

 subroutine mpp_get_boundary_DGRID(field1,field2,&
                                   is1,ie1,js1,je1,isd1,ied1,jsd1,jed1,km1,&
                                   is2,ie2,js2,je2,isd2,ied2,jsd2,jed2,km2,&
                                   npx,npy,npz,is,ie,js,je)

  integer :: is1,ie1,js1,je1,isd1,ied1,jsd1,jed1,km1
  integer :: is2,ie2,js2,je2,isd2,ied2,jsd2,jed2,km2
  real :: field1(isd1:ied1,jsd1:jed1,km1)
  real :: field2(isd2:ied2,jsd2:jed2,km2)
  real :: field_temp
  integer :: npx,npy,npz,is,ie,js,je
  real :: wbuffer(npy+2,npz)
  real :: sbuffer(npx+2,npz)

  integer :: i,j,k,my_flag

  wbuffer=0
  sbuffer=0
  call mpp_get_boundary(field1,field2,domain,wbuffery=wbuffer,&
                        ebuffery=field2(ie+1,js:je,1:npz),&
                                             sbufferx=sbuffer,&
                        nbufferx=field1(is:ie,je+1,1:npz),gridtype=DGRID_NE )

 end subroutine mpp_get_boundary_DGRID


 subroutine mpp_update_domains_dummy3(field1,is1,ie1,js1,je1,isd1,ied1,jsd1,jed1,km1,nq)
  integer :: is1,ie1,js1,je1,isd1,ied1,jsd1,jed1,km1,nq
  real :: field1(isd1:ied1,jsd1:jed1,km1,nq)
  real :: field_temp

  integer :: i,j,k,m


      do m=1,nq
      do k=1,km1
         field_temp=field1(is1,js1,k,m)
         field1(isd1:is1-1,:,k,m)=field_temp
         field1(:,jsd1:js1-1,k,m)=field_temp

         field_temp=field1(ie1,je1,k,m)
         field1(ie1+1:ied1,:,k,m)=field_temp
         field1(:,je1+1:jed1-1,k,m)=field_temp
      enddo
      enddo

 end subroutine mpp_update_domains_dummy3

 subroutine mpp_update_domains_dummy4(field1,is1,ie1,js1,je1,isd1,ied1,jsd1,jed1,km1,nq)
  integer :: is1,ie1,js1,je1,isd1,ied1,jsd1,jed1,km1,nq
  real :: field1(isd1:ied1,jsd1:jed1,km1,nq)
  real :: field_temp

  integer :: i,j,k,m


      do m=1,nq
      call mpp_update_domains(  field1(:,:,:,m), domain, complete=.true.)
!      do k=1,km1
!         field_temp=field1(is1,js1,k,m)
!         field1(isd1:is1-1,:,k,m)=field_temp
!         field1(:,jsd1:js1-1,k,m)=field_temp

!         field_temp=field1(ie1,je1,k,m)
!         field1(ie1+1:ied1,:,k,m)=field_temp
!         field1(:,je1+1:jed1-1,k,m)=field_temp
!      enddo
      enddo

 end subroutine mpp_update_domains_dummy4

 subroutine mpp_get_boundary_dummy1(field,is,ie,js,je,isd,ied,jsd,jed,km)
  integer :: is,ie,js,je,isd,ied,jsd,jed,km
  real :: field(isd:ied,jsd:jed,km)
  real :: field_temp

  integer :: i,j,k

      do k=1,km
         field_temp=field(is,js,k)
         field(isd:is-1,:,k)=field_temp
         field(:,jsd:js-1,k)=field_temp

         field_temp=field(ie,je,k)
         field(ie+1:ied,:,k)=field_temp
         field(:,je+1:jed,k)=field_temp
      enddo

 end subroutine mpp_get_boundary_dummy1

 subroutine mpp_get_boundary_dummy2(field1,field2,&
                                      is1,ie1,js1,je1,isd1,ied1,jsd1,jed1,km1,&
                                      is2,ie2,js2,je2,isd2,ied2,jsd2,jed2,km2)
  integer :: is1,ie1,js1,je1,isd1,ied1,jsd1,jed1,km1
  integer :: is2,ie2,js2,je2,isd2,ied2,jsd2,jed2,km2
  real :: field1(isd1:ied1,jsd1:jed1,km1)
  real :: field2(isd2:ied2,jsd2:jed2,km2)
  real :: field_temp

  integer :: i,j,k

      do k=1,km2
         field_temp=field2(is2,js2,k)
         field2(isd2:is2-1,:,k)=field_temp
         field2(:,jsd2:js2-1,k)=field_temp

         field_temp=field2(ie2,je2,k)
         field2(ie2+1:ied2-1,:,k)=field_temp
         field2(:,je2+1:jed2,k)=field_temp
      enddo

      do k=1,km1
         field_temp=field1(is1,js1,k)
         field1(isd1:is1-1,:,k)=field_temp
         field1(:,jsd1:js1-1,k)=field_temp

         field_temp=field1(ie1,je1,k)
         field1(ie1+1:ied1,:,k)=field_temp
         field1(:,je1+1:jed1-1,k)=field_temp
      enddo

 end subroutine mpp_get_boundary_dummy2


end module fv_my_mpp
