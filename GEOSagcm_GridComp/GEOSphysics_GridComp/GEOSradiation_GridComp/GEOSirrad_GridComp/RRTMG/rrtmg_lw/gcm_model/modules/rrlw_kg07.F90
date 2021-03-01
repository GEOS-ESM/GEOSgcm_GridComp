#include "_gpudef.inc"
#include "_memdef.inc" 
 module rrlw_kg07

       !use parkind ,only : im => kind , rb => kind 

      use memory
 implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 7
! band 7:  980-1080 cm-1 (low - h2o,o3; high - o3)
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
! selfrefo: real     
! forrefo : real     
!-----------------------------------------------------------------

      integer , parameter :: no7  = 16
      integer , parameter :: ng7  = 12
        real  _gpudev :: kaod(9,5,13,no7)
      real  _gpudev :: kbod(5,13:59,no7)
      real  _cpusnp :: ka(9,5,13,ng7) ,kb(5,13:59,ng7),absa(585,ng7)
      real  _cpusnp :: absb(235,ng7)

      real  _cpus, dimension(no7) :: fracrefbo
      real  _cpus :: fracrefao(no7,9)
      real  _cpus :: kao(9,5,13,no7)
      real  _cpus :: kbo(5,13:59,no7)
      real  _cpus :: kao_mco2(9,19,no7)
      real  _cpus :: kbo_mco2(19,no7)
      real  _cpus :: selfrefo(10,no7)
      real  _cpus :: forrefo(4,no7)

      real  _gpudeva , dimension(:) :: fracrefbod
      real  _gpudeva :: fracrefaod(:,:)
    
      real  _gpudeva :: kao_mco2d(:,:,:)
      real  _gpudeva :: kbo_mco2d(:,:)
      real  _gpudeva :: selfrefod(:,:)
      real  _gpudeva :: forrefod(:,:)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 7
! band 7:  980-1080 cm-1 (low - h2o,o3; high - o3)
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
! selfref : real     
! forref  : real     
!
! absa    : real
!-----------------------------------------------------------------

  

      real  _cpus, dimension(ng7) :: fracrefb
      real  _cpus :: fracrefa(ng7,9)
      
      real  _cpus :: ka_mco2(9,19,ng7)
      real  _cpus :: kb_mco2(19,ng7)
      real  _cpus :: selfref(10,ng7)
      real  _cpus :: forref(4,ng7)

      real  _gpudeva , dimension(:) :: fracrefbd
      real  _gpudeva :: fracrefad(:,:)
      real  _gpudevanp ::  absad(:,:)
      real  _gpudevanp ::  absbd(:,:)
      real  _gpudeva :: ka_mco2d(:,:,:)
      real  _gpudeva :: kb_mco2d(:,:)
      real  _gpudeva :: selfrefd(:,:)
      real  _gpudeva :: forrefd(:,:)
      equivalence (ka(1,1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

      contains

      subroutine copyToGPU7

       dbcop( fracrefb , fracrefbd )    
       dbcop( fracrefa , fracrefad )

        
       dbcopnp( absa , absad, 585,ng7 )
  
       dbcopnp( absb , absbd, 235,ng7 )
       dbcop( ka_mco2 , ka_mco2d )
       dbcop( kb_mco2 , kb_mco2d )
       dbcop( selfref , selfrefd )
       dbcop( forref , forrefd )

       dbcop( fracrefbo , fracrefbod )    
       dbcop( fracrefao , fracrefaod )
     
       dbcop( kao_mco2 , kao_mco2d )
       dbcop( kbo_mco2 , kbo_mco2d )
       dbcop( selfrefo , selfrefod )
       dbcop( forrefo , forrefod )


      end subroutine 

       subroutine reg7
       !67
       dbreg( fracrefb )    
       dbreg( fracrefa )

       !dbreg( ka )      
       dbreg( absa )
       !dbreg( kb )
       dbreg( absb )
       dbreg( ka_mco2 )
       dbreg( kb_mco2 )
       dbreg( selfref )
       dbreg( forref )

       dbreg( fracrefbo )    
       dbreg( fracrefao )
       !dbreg( kao )      
       !dbreg( kbo )
       !dbreg( absbo )
       dbreg( kao_mco2 )
       dbreg( kbo_mco2 )
       dbreg( selfrefo )
       dbreg( forrefo )


      end subroutine 
    

      end module rrlw_kg07
