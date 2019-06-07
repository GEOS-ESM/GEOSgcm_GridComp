#include "_gpudef.inc"  
#include "_memdef.inc"

  
      
      module rrlw_kg01

       !use parkind ,only : im => kind , rb => kind 


      use memory
      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 1
! band 1:  10-250 cm-1 (low - h2o; high - h2o)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!fracrefao: real    
!fracrefbo: real
! kao     : real     
! kbo     : real     
! kao_mn2 : real     
! kbo_mn2 : real     
! selfrefo: real     
! forrefo : real
!-----------------------------------------------------------------

      integer , parameter :: no1  = 16

      real  :: fracrefao(no1)  , fracrefbo(no1)
      real  :: kao(5,13,no1)
      real  :: kbo(5,13:59,no1)
      real  :: kao_mn2(19,no1) , kbo_mn2(19,no1)
      real  :: selfrefo(10,no1), forrefo(4,no1)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 1
! band 1:  10-250 cm-1 (low - h2o; high - h2o)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!fracrefa : real    
!fracrefb : real
! ka      : real     
! kb      : real     
! absa    : real
! absb    : real
! ka_mn2  : real     
! kb_mn2  : real     
! selfref : real     
! forref  : real
!-----------------------------------------------------------------

      integer , parameter :: ng1  = 10

      
      real  _cpusnp :: ka(5,13,ng1)   , absa(65,ng1)
      real  _cpusnp :: kb(5,13:59,ng1), absb(235,ng1)
      real  _cpus :: fracrefa(ng1)  , fracrefb(ng1)
      real  _cpus :: ka_mn2(19,ng1) , kb_mn2(19,ng1)
      real  _cpus :: selfref(10,ng1), forref(4,ng1)

      
      real  _gpudevanp :: kad(:,:,:), absad(:,:), absbd(:,:)
      real  _gpudevanp :: kbd(:,:,:)
      
      real  _gpudeva :: fracrefad(:)  , fracrefbd(:)
      real  _gpudeva :: ka_mn2d(:,:) , kb_mn2d(:,:)
      real  _gpudeva :: selfrefd(:,:), forrefd(:,:)

      equivalence (ka(1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

      contains

      subroutine copyToGPU1

     
        dbcop(fracrefa,fracrefad)
        dbcop(fracrefb,fracrefbd)
        dbcop(ka_mn2,ka_mn2d)
        dbcop(kb_mn2,kb_mn2d)
        dbcop(selfref,selfrefd)
        dbcop(forref,forrefd)
        dbcopnp(absa,absad, 65, ng1)
        dbcopnp(absb,absbd, 235, ng1)
     
      end subroutine 

      subroutine reg1

        dbreg(fracrefa)
        dbreg(fracrefb)
        dbreg(ka_mn2)
        dbreg(kb_mn2)
        dbreg(selfref)
        dbreg(forref)
        dbreg(absa)
        dbreg(absb)
        

      end subroutine 

      end module rrlw_kg01




