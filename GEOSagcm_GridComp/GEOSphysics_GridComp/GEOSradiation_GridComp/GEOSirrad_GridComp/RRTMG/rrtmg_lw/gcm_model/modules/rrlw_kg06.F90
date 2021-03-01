#include "_gpudef.inc"
#include "_memdef.inc" 
 module rrlw_kg06

       !use parkind ,only : im => kind , rb => kind 

      use memory

 implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 6
! band 6:  820-980 cm-1 (low - h2o; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!fracrefao: real    
! kao     : real     
! kao_mco2: real     
! selfrefo: real     
! forrefo : real     
!cfc11adjo: real
! cfc12o  : real
!-----------------------------------------------------------------

      integer , parameter :: no6  = 16
      integer , parameter :: ng6  = 8

      real  _cpusnp :: ka(5,13,ng6),absa(65,ng6)
      real  _cpus, dimension(no6) :: fracrefao
      real  _cpus :: kao(5,13,no6)
      real  _cpus :: kao_mco2(19,no6)
      real  _cpus :: selfrefo(10,no6)
      real  _cpus :: forrefo(4,no6)

      real  _cpus, dimension(no6) :: cfc11adjo
      real  _cpus, dimension(no6) :: cfc12o

      real  _gpudeva , dimension(:) :: fracrefaod
      real  _gpudeva :: kaod(:,:,:)
      real  _gpudeva :: kao_mco2d(:,:)
      real  _gpudeva :: selfrefod(:,:)
      real  _gpudeva :: forrefod(:,:)

      real  _gpudeva , dimension(:) :: cfc11adjod
      real  _gpudeva , dimension(:) :: cfc12od

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 6
! band 6:  820-980 cm-1 (low - h2o; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!fracrefa : real    
! ka      : real     
! ka_mco2 : real     
! selfref : real     
! forref  : real     
!cfc11adj : real
! cfc12   : real
!
! absa    : real
!-----------------------------------------------------------------

      

      real  _cpus, dimension(ng6) :: fracrefa
      
      real  _cpus :: ka_mco2(19,ng6)
      real  _cpus :: selfref(10,ng6)
      real  _cpus :: forref(4,ng6)

      real  _cpus, dimension(ng6) :: cfc11adj
      real  _cpus, dimension(ng6) :: cfc12

      real  _gpudeva , dimension(:) :: fracrefad
      real  _gpudevanp :: absad(:,:)
      real  _gpudeva :: ka_mco2d(:,:)
      real  _gpudeva :: selfrefd(:,:)
      real  _gpudeva :: forrefd(:,:)

      real  _gpudeva , dimension(:) :: cfc11adjd
      real  _gpudeva , dimension(:) :: cfc12d

      
      equivalence (ka(1,1,1),absa(1,1))

      contains 

      subroutine copyToGPU6

       dbcop( fracrefao , fracrefaod )    
       dbcop( kao , kaod )      
       dbcop( kao_mco2 , kao_mco2d )
       dbcop( selfrefo , selfrefod )
       dbcop( forrefo , forrefod )
       dbcop( cfc11adjo , cfc11adjod )
       dbcop( cfc12o , cfc12od )
      
       dbcop( fracrefa , fracrefad )
      
       dbcopnp( absa , absad, 65, ng6 )
       dbcop( ka_mco2 , ka_mco2d )
       dbcop( selfref , selfrefd )
       dbcop( forref , forrefd )
       dbcop( cfc11adj , cfc11adjd )
       dbcop( cfc12 , cfc12d )

      end subroutine 

      subroutine reg6
       !53
       dbreg( fracrefao )    
       dbreg( kao )      
       dbreg( kao_mco2 )
       dbreg( selfrefo )
       dbreg( forrefo )
       dbreg( cfc11adjo )
       dbreg( cfc12o )
      
       dbreg( fracrefa )
     
       dbreg( absa )
       dbreg( ka_mco2 )
       dbreg( selfref )
       dbreg( forref )
       dbreg( cfc11adj )
       dbreg( cfc12 )

      end subroutine 

      end module rrlw_kg06
