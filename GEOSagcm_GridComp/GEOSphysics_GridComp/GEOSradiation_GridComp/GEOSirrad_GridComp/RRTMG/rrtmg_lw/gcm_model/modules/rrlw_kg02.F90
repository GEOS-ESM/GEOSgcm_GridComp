#include "_gpudef.inc"
#include "_memdef.inc"
      module rrlw_kg02

       !use parkind ,only : im => kind , rb => kind 
      use memory
      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 2
! band 2:  250-500 cm-1 (low - h2o; high - h2o)
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
! selfrefo: real     
! forrefo : real
!-----------------------------------------------------------------

      integer , parameter :: no2  = 16
      real  _cpus :: kao(5,13,no2)
      real  _cpus :: kbo(5,13:59,no2)
      real  _cpus :: fracrefao(no2)   , fracrefbo(no2)
      real  _cpus :: selfrefo(10,no2) , forrefo(4,no2)

      real  _gpudeva :: fracrefaod(:)   , fracrefbod(:)
      real  _gpudeva :: selfrefod(:,:) , forrefod(:,:)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 2
! band 2:  250-500 cm-1 (low - h2o; high - h2o)
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
! selfref : real     
! forref  : real
!
! refparam: real
!-----------------------------------------------------------------

      integer , parameter :: ng2  = 12

      real  _cpus :: fracrefa(ng2)  , fracrefb(ng2)
      real  _cpusnp :: ka(5,13,ng2)   , absa(65,ng2)
      real  _cpusnp :: kb(5,13:59,ng2), absb(235,ng2)
      real  _cpus :: selfref(10,ng2), forref(4,ng2)

      real  _gpudeva :: fracrefad(:)  , fracrefbd(:)
      real  _gpudevanp :: absad(:,:)
      real  _gpudevanp :: absbd(:,:)
      real  _gpudeva :: selfrefd(:,:), forrefd(:,:)

      real  :: refparam(13)

      equivalence (ka(1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

      contains
      
      subroutine copyToGPU2

        dbcop(fracrefao,fracrefaod)
        dbcop(fracrefbo,fracrefbod)       
        dbcop(selfrefo, selfrefod)
        dbcop(forrefo, forrefod)

        dbcop(fracrefa,fracrefad)
        dbcop(fracrefb,fracrefbd)       
        dbcopnp(absa,absad, 65,ng2)
        dbcopnp(absb,absbd, 235,ng2)
        dbcop(selfref, selfrefd)
        dbcop(forref, forrefd)
        

       
        
      end subroutine 
        
      subroutine reg2
         ! 9
        dbreg(fracrefao)
        dbreg(fracrefbo)        
        dbreg(selfrefo)
        dbreg(forrefo)
         
        dbreg(fracrefa)
        dbreg(fracrefb)        
        dbreg(absa)        
        dbreg(absb)
        dbreg(selfref)
        dbreg(forref)

      end subroutine

      end module rrlw_kg02


