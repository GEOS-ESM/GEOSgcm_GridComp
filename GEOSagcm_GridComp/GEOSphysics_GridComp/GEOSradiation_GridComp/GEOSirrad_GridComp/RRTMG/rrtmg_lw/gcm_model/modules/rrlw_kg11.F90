#include "_gpudef.inc"
#include "_memdef.inc" 
 module rrlw_kg11

       !use parkind ,only : im => kind , rb => kind 

      use memory
 implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 11
! band 11:  1480-1800 cm-1 (low - h2o; high - h2o)
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
! kao_mo2 : real     
! kbo_mo2 : real     
! selfrefo: real     
! forrefo : real     
!-----------------------------------------------------------------

      integer , parameter :: no11 = 16

      real  _cpus, dimension(no11) :: fracrefao
      real  _cpus, dimension(no11) :: fracrefbo

      real  _cpus :: kao(5,13,no11)
      real  _cpus :: kbo(5,13:59,no11)
      real  _cpus :: kao_mo2(19,no11)
      real  _cpus :: kbo_mo2(19,no11)
      real  _cpus :: selfrefo(10,no11)
      real  _cpus :: forrefo(4,no11)

      real  _gpudeva , dimension(:) :: fracrefaod
      real  _gpudeva , dimension(:) :: fracrefbod

      real  _gpudev :: kaod(5,13,no11)
      real  _gpudev :: kbod(5,13:59,no11)
      real  _gpudeva :: kao_mo2d(:,:)
      real  _gpudeva :: kbo_mo2d(:,:)
      real  _gpudeva :: selfrefod(:,:)
      real  _gpudeva :: forrefod(:,:)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 11
! band 11:  1480-1800 cm-1 (low - h2o; high - h2o)
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
! ka_mo2  : real     
! kb_mo2  : real     
! selfref : real     
! forref  : real     
!
! absa    : real
! absb    : real
!-----------------------------------------------------------------

      integer , parameter :: ng11 = 8

      real  _cpus, dimension(ng11) :: fracrefa
      real  _cpus, dimension(ng11) :: fracrefb

      real  _cpusnp :: ka(5,13,ng11)   , absa(65,ng11)
      real  _cpusnp :: kb(5,13:59,ng11), absb(235,ng11)
      real  _cpus :: ka_mo2(19,ng11)
      real  _cpus :: kb_mo2(19,ng11)
      real  _cpus :: selfref(10,ng11)
      real  _cpus :: forref(4,ng11)

      real  _gpudeva , dimension(:) :: fracrefad
      real  _gpudeva , dimension(:) :: fracrefbd

      real  _gpudevanp ::   absad(:,:)
      real  _gpudevanp ::   absbd(:,:)
      real  _gpudeva :: ka_mo2d(:,:)
      real  _gpudeva :: kb_mo2d(:,:)
      real  _gpudeva :: selfrefd(:,:)
      real  _gpudeva :: forrefd(:,:)

      equivalence (ka(1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

      contains 

      subroutine copyToGPU11


     dbcop( fracrefa , fracrefad )
     dbcop( fracrefb , fracrefbd )

     
     dbcopnp( absa , absad, 65,  ng11 )     
     dbcopnp( absb , absbd, 235, ng11 )
     dbcop( ka_mo2 , ka_mo2d )
     dbcop( kb_mo2 , kb_mo2d )
     dbcop( selfref , selfrefd )
     dbcop( forref , forrefd )

      end subroutine 

     subroutine reg11


     dbreg( fracrefa )
     dbreg( fracrefb )

     !dbreg( ka ) 
     dbreg( absa )
     !dbreg( kb )
     dbreg( absb )
     dbreg( ka_mo2 )
     dbreg( kb_mo2 )
     dbreg( selfref )
     dbreg( forref )

      end subroutine 


      end module rrlw_kg11
