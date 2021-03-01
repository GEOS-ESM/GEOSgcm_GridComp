#include "_gpudef.inc"
#include "_memdef.inc" 
 module rrlw_kg08

       !use parkind ,only : im => kind , rb => kind 

      use memory
 implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 8
! band 8:  1080-1180 cm-1 (low (i.e.>~300mb) - h2o; high - o3)
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
! kao_mco2: real     
! kbo_mco2: real     
! kao_mn2o: real     
! kbo_mn2o: real     
! kao_mo3 : real     
! selfrefo: real     
! forrefo : real     
! cfc12o  : real     
!cfc22adjo: real     
!-----------------------------------------------------------------

      integer , parameter :: no8  = 16

      real  _cpus, dimension(no8) :: fracrefao
      real  _cpus, dimension(no8) :: fracrefbo
      real  _cpus, dimension(no8) :: cfc12o
      real  _cpus, dimension(no8) :: cfc22adjo

      real  _cpus :: kao(5,13,no8)
      real  _cpus :: kao_mco2(19,no8)
      real  _cpus :: kao_mn2o(19,no8)
      real  _cpus :: kao_mo3(19,no8)
      real  _cpus :: kbo(5,13:59,no8)
      real  _cpus :: kbo_mco2(19,no8)
      real  _cpus :: kbo_mn2o(19,no8)
      real  _cpus :: selfrefo(10,no8)
      real  _cpus :: forrefo(4,no8)

      real  _gpudeva , dimension(:) :: fracrefaod
      real  _gpudeva , dimension(:) :: fracrefbod
      real  _gpudeva , dimension(:) :: cfc12od
      real  _gpudeva , dimension(:) :: cfc22adjod

      real  _gpudev :: kaod(5,13,no8)
      real  _gpudeva :: kao_mco2d(:,:)
      real  _gpudeva :: kao_mn2od(:,:)
      real  _gpudeva :: kao_mo3d(:,:)
      real  _gpudev :: kbod(5,13:59,no8)
      real  _gpudeva :: kbo_mco2d(:,:)
      real  _gpudeva :: kbo_mn2od(:,:)
      real  _gpudeva :: selfrefod(:,:)
      real  _gpudeva :: forrefod(:,:)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 8
! band 8:  1080-1180 cm-1 (low (i.e.>~300mb) - h2o; high - o3)
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
! ka_mco2 : real     
! kb_mco2 : real     
! ka_mn2o : real     
! kb_mn2o : real     
! ka_mo3  : real     
! selfref : real     
! forref  : real     
! cfc12   : real     
! cfc22adj: real     
!
! absa    : real
! absb    : real
!-----------------------------------------------------------------

      integer , parameter :: ng8  = 8

      real  _cpus, dimension(ng8) :: fracrefa
      real  _cpus, dimension(ng8) :: fracrefb
      real  _cpus, dimension(ng8) :: cfc12
      real  _cpus, dimension(ng8) :: cfc22adj

      real  _cpusnp :: ka(5,13,ng8)    ,absa(65,ng8)
      real  _cpusnp :: kb(5,13:59,ng8) ,absb(235,ng8)
      real  _cpus :: ka_mco2(19,ng8)
      real  _cpus :: ka_mn2o(19,ng8)
      real  _cpus :: ka_mo3(19,ng8)
      real  _cpus :: kb_mco2(19,ng8)
      real  _cpus :: kb_mn2o(19,ng8)
      real  _cpus :: selfref(10,ng8)
      real  _cpus :: forref(4,ng8)

      real  _gpudeva  , dimension(:) :: fracrefad
      real  _gpudeva  , dimension(:) :: fracrefbd
      real  _gpudeva  , dimension(:) :: cfc12d
      real  _gpudeva  , dimension(:) :: cfc22adjd

      real  _gpudevanp  ::  absad(:,:)
      real  _gpudevanp  ::  absbd(:,:)
      real  _gpudeva  :: ka_mco2d(:,:)
      real  _gpudeva  :: ka_mn2od(:,:)
      real  _gpudeva  :: ka_mo3d(:,:)
      real  _gpudeva  :: kb_mco2d(:,:)
      real  _gpudeva  :: kb_mn2od(:,:)
      real  _gpudeva  :: selfrefd(:,:)
      real  _gpudeva  :: forrefd(:,:)

      equivalence (ka(1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

       contains 

      subroutine copyToGPU8

      kaod = kao
      kbod = kbo

     dbcop( fracrefao , fracrefaod )
     dbcop( fracrefbo , fracrefbod )
     dbcop( cfc12o , cfc12od )
     dbcop( cfc22adjo , cfc22adjod )
   
     dbcop( kao_mco2 , kao_mco2d )
     dbcop( kao_mn2o , kao_mn2od )
     dbcop( kao_mo3 , kao_mo3d )
     
     dbcop( kbo_mco2 , kbo_mco2d )
     dbcop( kbo_mn2o , kbo_mn2od )
     dbcop( selfrefo , selfrefod )
     dbcop( forrefo , forrefod )


     dbcop( fracrefa , fracrefad )
     dbcop( fracrefb , fracrefbd )
     dbcop( cfc12 , cfc12d )
     dbcop( cfc22adj , cfc22adjd )
    
     dbcopnp( absa , absad, 65,ng8 )
     dbcopnp( absb , absbd, 235,ng8 )
     dbcop( ka_mco2 , ka_mco2d )
     dbcop( ka_mn2o , ka_mn2od )
     dbcop( ka_mo3 , ka_mo3d )
     dbcop( kb_mco2 , kb_mco2d )
     dbcop( kb_mn2o , kb_mn2od )
     dbcop( selfref , selfrefd )
     dbcop( forref , forrefd )


      end subroutine 

     subroutine reg8
     
     dbreg( fracrefao )
     dbreg( fracrefbo )
     dbreg( cfc12o )
     dbreg( cfc22adjo )

     dbreg( kao_mco2 )
     dbreg( kao_mn2o )
     dbreg( kao_mo3 )

     dbreg( kbo_mco2 )
     dbreg( kbo_mn2o )
     dbreg( selfrefo )
     dbreg( forrefo )

     dbreg( fracrefa )
     dbreg( fracrefb )
     dbreg( cfc12 )
     dbreg( cfc22adj )
     dbreg( absa )
     dbreg( absb )
     dbreg( ka_mco2 )
     dbreg( ka_mn2o )
     dbreg( ka_mo3 )
     dbreg( kb_mco2 )
     dbreg( kb_mn2o )
     dbreg( selfref )
     dbreg( forref )

      end subroutine 

      end module rrlw_kg08

