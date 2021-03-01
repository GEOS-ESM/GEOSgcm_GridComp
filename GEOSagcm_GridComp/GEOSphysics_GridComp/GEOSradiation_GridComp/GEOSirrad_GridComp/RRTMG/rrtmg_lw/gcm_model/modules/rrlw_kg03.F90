#include "_gpudef.inc"
#include "_memdef.inc" 
 module rrlw_kg03

       !use parkind ,only : im => kind , rb => kind 

      use memory
 implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 3
! band 3:  500-630 cm-1 (low - h2o,co2; high - h2o,co2)
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
! kao_mn2o: real     
! kbo_mn2o: real     
! selfrefo: real     
! forrefo : real
!-----------------------------------------------------------------

      integer , parameter :: no3  = 16

      real  _cpus :: fracrefao(no3,9) ,fracrefbo(no3,5)
      real  _cpus :: kao(9,5,13,no3)
      real  _cpus :: kbo(5,5,13:59,no3)
      real  _cpus :: kao_mn2o(9,19,no3), kbo_mn2o(5,19,no3)
      real  _cpus :: selfrefo(10,no3)
      real  _cpus :: forrefo(4,no3)

      real  _gpudeva :: fracrefaod(:,:) ,fracrefbod(:,:)
      !real  _gpudeva :: kaod(9,5,13,no3)
      !real  _gpudeva :: kbod(5,5,13:59,no3)
      real  _gpudeva :: kao_mn2od(:,:,:), kbo_mn2od(:,:,:)
      real  _gpudeva :: selfrefod(:,:)
      real  _gpudeva :: forrefod(:,:)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 3
! band 3:  500-630 cm-1 (low - h2o,co2; high - h2o,co2)
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
! ka_mn2o : real     
! kb_mn2o : real     
! selfref : real     
! forref  : real
!
! absa    : real
! absb    : real
!-----------------------------------------------------------------

      integer , parameter :: ng3  = 16

      real  _cpus :: fracrefa(ng3,9) ,fracrefb(ng3,5)
      real  _cpusnp :: ka(9,5,13,ng3)  ,absa(585,ng3)
      real  _cpusnp :: kb(5,5,13:59,ng3),absb(1175,ng3)
      real  _cpus :: ka_mn2o(9,19,ng3), kb_mn2o(5,19,ng3)
      real  _cpus :: selfref(10,ng3)
      real  _cpus :: forref(4,ng3)

      real  _gpudeva :: fracrefad(:,:) ,fracrefbd(:,:)
      real  _gpudevanp :: absad(:,:)
      real  _gpudevanp :: absbd(:,:)
      real  _gpudeva :: ka_mn2od(:,:,:), kb_mn2od(:,:,:)
      real  _gpudeva :: selfrefd(:,:)
      real  _gpudeva :: forrefd(:,:)

      equivalence (ka(1,1,1,1),absa(1,1)),(kb(1,1,13,1),absb(1,1))

      contains 

      subroutine copyToGPU3

      
       dbcop( fracrefao , fracrefaod )
       dbcop( fracrefbo , fracrefbod )
       dbcop( kao_mn2o , kao_mn2od )
       dbcop( kbo_mn2o , kbo_mn2od )
       dbcop( selfrefo , selfrefod )
       dbcop( forrefo , forrefod )

       dbcop( fracrefa , fracrefad )
       dbcop( fracrefb , fracrefbd )
   
       dbcopnp( absa , absad, 585,ng3 )
    
       dbcopnp( absb , absbd, 1175, ng3 )
       dbcop( ka_mn2o , ka_mn2od )
       dbcop( kb_mn2o , kb_mn2od )
       dbcop( selfref , selfrefd )
       dbcop( forref , forrefd )

      end subroutine 

      subroutine reg3
       !19
       dbreg( fracrefao )
       dbreg( fracrefbo )
     
       dbreg( kao_mn2o )
       dbreg( kbo_mn2o )
       dbreg( selfrefo )
       dbreg( forrefo )

       dbreg( fracrefa )
       dbreg( fracrefb )
      
       dbreg( absa )
     
       dbreg( absb )
       dbreg( ka_mn2o )
       dbreg( kb_mn2o )
       dbreg( selfref )
       dbreg( forref )


      end subroutine

      end module rrlw_kg03


