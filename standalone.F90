program vtrisolvesurf_standalone

   implicit none

   real, dimension(:,:,:), allocatable :: BKS, CKS, DKS, DKS_ref
   real, dimension(:,:,:), allocatable :: BKQ, CKQ, DKQ, DKQ_ref
   real, dimension(:,:,:), allocatable :: BKV, CKV, DKV, DKV_ref

   real, dimension(:,:,:), allocatable :: BKSS, CKSS, DKSS, DKSS_ref
   real, dimension(:,:,:), allocatable :: BKQQ, CKQQ, DKQQ, DKQQ_ref
   real, dimension(:,:,:), allocatable :: BKUU, CKUU, DKUU, DKUU_ref

   integer :: IM, JM, LM, fileID

   character*100 :: dirName, rank_str

   real :: t1, t2

   dirName = './c90-90-90-72'
   rank_str = '0'

   ! C12 Data Set
   if (dirName == './c12-6-6-72') then
      IM = 6
      JM = 6
      LM = 72
   else if (dirName == './c90-90-90-72') then
      IM = 90
      JM = 90
      LM = 72
   else if(dirName == './c180-180-180-72') then
      IM = 180
      JM = 180
      LM = 72
   endif

   allocate(BKS(IM, JM, LM))
   allocate(CKS(IM, JM, LM))
   allocate(DKS(IM, JM, LM))
   allocate(DKS_ref(IM, JM, LM))

   open(newunit=fileID, file=trim(dirName) // "/BKS_" // trim(rank_str) // ".out", &
        status='old', form="unformatted", action="read")
   read(fileID) BKS
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // "/CKS_" // trim(rank_str) // ".in", &
        status='old', form="unformatted", action="read")
   read(fileID) CKS
   close(fileID)
    
   open(newunit=fileID, file=trim(dirName) // "/DKS_" // trim(rank_str) // ".out", &
        status='old', form="unformatted", action="read")
   read(fileID) DKS_ref
   close(fileID)

!$acc data copyin(BKS, CKS) copyout(DKS)

   call cpu_time(t1)
   call VTRISOLVESURF(BKS,CKS,DKS)
   call cpu_time(t2)

!$acc end data

   print*,'SUM DIFF DKS = ', sum(DKS-DKS_ref)
   print*,'Runtime = ', t2-t1

   allocate(BKQ(IM, JM, LM))
   allocate(CKQ(IM, JM, LM))
   allocate(DKQ(IM, JM, LM))
   allocate(DKQ_ref(IM, JM, LM))

   open(newunit=fileID, file=trim(dirName) // "/BKQ_" // trim(rank_str) // ".out", &
        status='old', form="unformatted", action="read")
   read(fileID) BKQ
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // "/CKQ_" // trim(rank_str) // ".in", &
        status='old', form="unformatted", action="read")
   read(fileID) CKQ
   close(fileID)
    
   open(newunit=fileID, file=trim(dirName) // "/DKQ_" // trim(rank_str) // ".out", &
        status='old', form="unformatted", action="read")
   read(fileID) DKQ_ref
   close(fileID)

!$acc data copyin(BKQ, CKQ) copyout(DKQ)

   call cpu_time(t1)
   call VTRISOLVESURF(BKQ,CKQ,DKQ)
   call cpu_time(t2)

!$acc end data

   print*,'SUM DIFF DKQ = ', sum(DKQ-DKQ_ref)
   print*,'Runtime = ', t2-t1

   allocate(BKV(IM, JM, LM))
   allocate(CKV(IM, JM, LM))
   allocate(DKV(IM, JM, LM))
   allocate(DKV_ref(IM, JM, LM))

   open(newunit=fileID, file=trim(dirName) // "/BKV_" // trim(rank_str) // ".out", &
        status='old', form="unformatted", action="read")
   read(fileID) BKV
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // "/CKV_" // trim(rank_str) // ".in", &
        status='old', form="unformatted", action="read")
   read(fileID) CKV
   close(fileID)
    
   open(newunit=fileID, file=trim(dirName) // "/DKV_" // trim(rank_str) // ".out", &
        status='old', form="unformatted", action="read")
   read(fileID) DKV_ref
   close(fileID)

!$acc data copyin(BKV, CKV) copyout(DKV)

   call cpu_time(t1)
   call VTRISOLVESURF(BKV,CKV,DKV)
   call cpu_time(t2)

!$acc end data

   print*,'SUM DIFF DKV = ', sum(DKV-DKV_ref)
   print*,'Runtime = ', t2-t1

   allocate(BKSS(IM, JM, LM))
   allocate(CKSS(IM, JM, LM))
   allocate(DKSS(IM, JM, LM))
   allocate(DKSS_ref(IM, JM, LM))

   open(newunit=fileID, file=trim(dirName) // "/BKSS_" // trim(rank_str) // ".out", &
        status='old', form="unformatted", action="read")
   read(fileID) BKSS
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // "/CKSS_" // trim(rank_str) // ".in", &
        status='old', form="unformatted", action="read")
   read(fileID) CKSS
   close(fileID)
    
   open(newunit=fileID, file=trim(dirName) // "/DKSS_" // trim(rank_str) // ".out", &
        status='old', form="unformatted", action="read")
   read(fileID) DKSS_ref
   close(fileID)

!$acc data copyin(BKSS, CKSS) copyout(DKSS)

   call cpu_time(t1)
   call VTRISOLVESURF(BKSS,CKSS,DKSS)
   call cpu_time(t2)

!$acc end data

   print*,'SUM DIFF DKSS = ', sum(DKSS-DKSS_ref)
   print*,'Runtime = ', t2-t1

   allocate(BKQQ(IM, JM, LM))
   allocate(CKQQ(IM, JM, LM))
   allocate(DKQQ(IM, JM, LM))
   allocate(DKQQ_ref(IM, JM, LM))

   open(newunit=fileID, file=trim(dirName) // "/BKQQ_" // trim(rank_str) // ".out", &
        status='old', form="unformatted", action="read")
   read(fileID) BKQQ
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // "/CKQQ_" // trim(rank_str) // ".in", &
        status='old', form="unformatted", action="read")
   read(fileID) CKQQ
   close(fileID)
    
   open(newunit=fileID, file=trim(dirName) // "/DKQQ_" // trim(rank_str) // ".out", &
        status='old', form="unformatted", action="read")
   read(fileID) DKQQ_ref
   close(fileID)

!$acc data copyin(BKQQ, CKQQ) copyout(DKQQ)

   call cpu_time(t1)
   call VTRISOLVESURF(BKQQ,CKQQ,DKQQ)
   call cpu_time(t2)

!$acc end data

   print*,'SUM DIFF DKQQ = ', sum(DKQQ-DKQQ_ref)
   print*,'Runtime = ', t2-t1

   allocate(BKUU(IM, JM, LM))
   allocate(CKUU(IM, JM, LM))
   allocate(DKUU(IM, JM, LM))
   allocate(DKUU_ref(IM, JM, LM))

   open(newunit=fileID, file=trim(dirName) // "/BKUU_" // trim(rank_str) // ".out", &
        status='old', form="unformatted", action="read")
   read(fileID) BKUU
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // "/CKUU_" // trim(rank_str) // ".in", &
        status='old', form="unformatted", action="read")
   read(fileID) CKUU
   close(fileID)
    
   open(newunit=fileID, file=trim(dirName) // "/DKUU_" // trim(rank_str) // ".out", &
        status='old', form="unformatted", action="read")
   read(fileID) DKUU_ref
   close(fileID)

!$acc data copyin(BKUU, CKUU) copyout(DKUU)

   call cpu_time(t1)
   call VTRISOLVESURF(BKUU,CKUU,DKUU)
   call cpu_time(t2)

!$acc end data

   print*,'SUM DIFF DKUU = ', sum(DKUU-DKUU_ref)
   print*,'Runtime = ', t2-t1

   contains
! !IROUTINE:  VTRISOLVESURF -- Solves for sensitivity to surface value


! !INTERFACE:

   subroutine VTRISOLVESURF(B,C,Y)

    ! !ARGUMENTS:
    
      real,    dimension(:,:,:), intent(IN   ) :: B, C
      real,    dimension(:,:,:), intent(  OUT) :: Y
    
    ! !DESCRIPTION: Solves tridiagonal system that has been LU decomposed
    !   for the special case
    !   where the surface Y (YG) is 1 and the rest of the input Ys are 0.
    !   Everything else is as in {\tt VTRISOLVE}. This gives the sensitivity of the
    !   solution to a unit change in the surface values.
    
    !EOP
    
      integer :: IM, JM, LM, I, J, L

      IM = size(B,1)
      JM = size(B,2)
      LM = size(B,3)

!$acc data present(B, C, Y)

!$acc parallel loop gang vector collapse(2)
      do J = 1,JM
         do I = 1,IM
            Y(I,J,LM) = -C(I,J,LM) * B(I,J,LM)
         enddo
      enddo
!$acc end parallel loop

!$acc parallel
!$acc loop seq
      do L = LM-1,1,-1
!$acc loop gang vector collapse(2)
         do J = 1, JM
            do I = 1,IM
               Y(I,J,L) = -C(I,J,L) * Y(I,J,L+1) * B(I,J,L)
            enddo
         enddo
      enddo
!$acc end parallel

!$acc end data
    
   end subroutine VTRISOLVESURF
end program

! NASA Docket No. GSC-15,354-1, and identified as "GEOS-5 GCM Modeling Software”

! “Copyright © 2008 United States Government as represented by the Administrator
! of the National Aeronautics and Space Administration. All Rights Reserved.”

! Licensed under the Apache License, Version 2.0 (the "License"); you may not use
! this file except in compliance with the License. You may obtain a copy of the
! License at

! http://www.apache.org/licenses/LICENSE-2.0

! Unless required by applicable law or agreed to in writing, software distributed
! under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
! CONDITIONS OF ANY KIND, either express or implied. See the License for the
! specific language governing permissions and limitations under the License.
