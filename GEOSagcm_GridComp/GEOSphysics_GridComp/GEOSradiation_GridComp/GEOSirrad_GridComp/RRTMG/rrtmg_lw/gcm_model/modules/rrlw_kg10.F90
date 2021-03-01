#include "_gpudef.inc"
#include "_memdef.inc" 
 module rrlw_kg10

       !use parkind ,only : im => kind , rb => kind 

      use memory
 implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 10
! band 10:  1390-1480 cm-1 (low - h2o; high - h2o)
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

      integer , parameter :: no10 = 16

      real  _cpus, dimension(no10) :: fracrefao
      real  _cpus, dimension(no10) :: fracrefbo

      real  _cpus :: kao(5,13,no10)
      real  _cpus :: kbo(5,13:59,no10)
      real  _cpus :: selfrefo(10,no10)
      real  _cpus :: forrefo(4,no10)

      real  _gpudeva , dimension(:) :: fracrefaod
      real  _gpudeva , dimension(:) :: fracrefbod

      real  _gpudev :: kaod(5,13,no10)
      real  _gpudev :: kbod(5,13:59,no10)
      real  _gpudeva :: selfrefod(:,:)
      real  _gpudeva :: forrefod(:,:)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 10
! band 10:  1390-1480 cm-1 (low - h2o; high - h2o)
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
! selfref : real     
! forref  : real     
!
! absa    : real
! absb    : real
!-----------------------------------------------------------------

      integer , parameter :: ng10 = 6

      real  _cpus , dimension(ng10) :: fracrefa
      real  _cpus , dimension(ng10) :: fracrefb

      real  _cpusnp :: ka(5,13,ng10)   , absa(65,ng10)
      real  _cpusnp :: kb(5,13:59,ng10), absb(235,ng10)
      real  _cpus :: selfref(10,ng10)
      real  _cpus :: forref(4,ng10)

      real  _gpudeva , dimension(:) :: fracrefad
      real  _gpudeva , dimension(:) :: fracrefbd

      real  _gpudevanp ::   absad(:,:)
      real  _gpudevanp ::   absbd(:,:)
      real  _gpudeva :: selfrefd(:,:)
      real  _gpudeva :: forrefd(:,:)
      
      equivalence (ka(1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

      contains 

      subroutine copyToGPU10

      kaod = kao
      kbod = kbo

     dbcop( fracrefao , fracrefaod )
     dbcop( fracrefbo , fracrefbod )

     !dbcop( kao , kaod ) 
     !dbcop( kbo , kbod )
     dbcop( selfrefo , selfrefod )
     dbcop( forrefo , forrefod )

     dbcop( fracrefa , fracrefad )
     dbcop( fracrefb , fracrefbd )

     !dbcop( ka , kad ) 
     !dbcop( kb , kbd )
     dbcopnp( absa , absad, 65,ng10 )
     dbcopnp( absb , absbd, 235,ng10)
     dbcop( selfref , selfrefd )
     dbcop( forref , forrefd )

      end subroutine 

    subroutine reg10

     dbreg( fracrefao )
     dbreg( fracrefbo )

     !dbreg( kao ) 
     !dbreg( kbo )
     dbreg( selfrefo )
     dbreg( forrefo )

     dbreg( fracrefa )
     dbreg( fracrefb )

     !dbreg( ka ) 
     !dbreg( kb )
     dbreg( absa )
     dbreg( absb )
     dbreg( selfref )
     dbreg( forref )

      end subroutine 

      end module rrlw_kg10
