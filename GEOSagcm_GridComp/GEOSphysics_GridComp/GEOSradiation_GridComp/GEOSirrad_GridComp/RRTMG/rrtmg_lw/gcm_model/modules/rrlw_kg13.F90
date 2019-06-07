#include "_gpudef.inc"
#include "_memdef.inc" 
 module rrlw_kg13

       !use parkind ,only : im => kind , rb => kind 

      use memory
 implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 13
! band 13:  2080-2250 cm-1 (low - h2o,n2o; high - nothing)
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
! kao_mco : real     
! kbo_mo3 : real     
! selfrefo: real     
! forrefo : real     
!-----------------------------------------------------------------

      integer , parameter :: no13 = 16

      real  _cpus, dimension(no13) :: fracrefbo

      real  _cpus :: fracrefao(no13,9)
      real  _cpus :: kao(9,5,13,no13)
      real  _cpus :: kao_mco2(9,19,no13)
      real  _cpus :: kao_mco(9,19,no13)
      real  _cpus :: kbo_mo3(19,no13)
      real  _cpus :: selfrefo(10,no13)
      real  _cpus :: forrefo(4,no13)

      real  _gpudeva , dimension(:) :: fracrefbod

      real  _gpudeva  :: fracrefaod(:,:)
      real  _gpudev  :: kaod(9,5,13,no13)
      real  _gpudeva  :: kao_mco2d(:,:,:)
      real  _gpudeva  :: kao_mcod(:,:,:)
      real  _gpudeva  :: kbo_mo3d(:,:)
      real  _gpudeva  :: selfrefod(:,:)
      real  _gpudeva  :: forrefod(:,:)


!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 13
! band 13:  2080-2250 cm-1 (low - h2o,n2o; high - nothing)
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
! ka_mco  : real     
! kb_mo3  : real     
! selfref : real     
! forref  : real     
!
! absa    : real
!-----------------------------------------------------------------

      integer , parameter :: ng13 = 4

      real  _cpus, dimension(ng13) :: fracrefb

      real  _cpus :: fracrefa(ng13,9)
      real  _cpusnp :: ka(9,5,13,ng13) ,absa(585,ng13)
      real  _cpus :: ka_mco2(9,19,ng13)
      real  _cpus :: ka_mco(9,19,ng13)
      real  _cpus :: kb_mo3(19,ng13)
      real  _cpus :: selfref(10,ng13)
      real  _cpus :: forref(4,ng13)

      real  _gpudeva , dimension(:) :: fracrefbd

      real  _gpudeva :: fracrefad(:,:)
      real  _gpudevanp ::  absad(:,:)
      real  _gpudeva :: ka_mco2d(:,:,:)
      real  _gpudeva :: ka_mcod(:,:,:)
      real  _gpudeva :: kb_mo3d(:,:)
      real  _gpudeva :: selfrefd(:,:)
      real  _gpudeva :: forrefd(:,:)

      equivalence (ka(1,1,1,1),absa(1,1))

      contains
      
      subroutine copyToGPU13
      kaod = kao
     dbcop( fracrefbo , fracrefbod )
     dbcop( fracrefao , fracrefaod )
    
     dbcop( kao_mco2 , kao_mco2d )
     dbcop( kao_mco , kao_mcod )
     dbcop( kbo_mo3 , kbo_mo3d )
     dbcop( selfrefo , selfrefod )
     dbcop( forrefo , forrefod )

     dbcop( fracrefb , fracrefbd )
     dbcop( fracrefa , fracrefad )

     dbcopnp( absa , absad , 585,ng13)
     dbcop( ka_mco2 , ka_mco2d )
     dbcop( ka_mco , ka_mcod )
     dbcop( kb_mo3 , kb_mo3d )
     dbcop( selfref , selfrefd )
     dbcop( forref , forrefd )




      end subroutine

            
      subroutine reg13

     dbreg( fracrefbo )
     dbreg( fracrefao )
     !dbreg( kao ) 
     dbreg( kao_mco2 )
     dbreg( kao_mco )
     dbreg( kbo_mo3 )
     dbreg( selfrefo )
     dbreg( forrefo )

     dbreg( fracrefb )
     dbreg( fracrefa )
     !dbreg( ka ) 
     dbreg( absa )
     dbreg( ka_mco2 )
     dbreg( ka_mco )
     dbreg( kb_mo3 )
     dbreg( selfref )
     dbreg( forref )




      end subroutine
      end module rrlw_kg13
