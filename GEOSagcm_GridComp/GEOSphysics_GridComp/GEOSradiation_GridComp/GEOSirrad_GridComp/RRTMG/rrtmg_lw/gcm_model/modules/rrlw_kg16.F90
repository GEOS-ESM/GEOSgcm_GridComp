#include "_gpudef.inc"
#include "_memdef.inc" 
 module rrlw_kg16

       !use parkind ,only : im => kind , rb => kind 

      use memory
 implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 16
! band 16:  2600-3000 cm-1 (low - h2o,ch4; high - nothing)
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
! kbo     : real     
! selfrefo: real     
! forrefo : real     
!-----------------------------------------------------------------

      integer , parameter :: no16 = 16

      real  _cpus, dimension(no16) :: fracrefbo

      real  _cpus :: fracrefao(no16,9)
      real  _cpus :: kao(9,5,13,no16)
      real  _cpus :: kbo(5,13:59,no16)
      real  _cpus :: selfrefo(10,no16)
      real  _cpus :: forrefo(4,no16)
      
      real  _gpudeva , dimension(:) :: fracrefbod
      real  _gpudeva :: fracrefaod(:,:)
      real  _gpudev :: kaod(9,5,13,no16)
      real  _gpudev :: kbod(5,13:59,no16)
      real  _gpudeva :: selfrefod(:,:)
      real  _gpudeva :: forrefod(:,:)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 16
! band 16:  2600-3000 cm-1 (low - h2o,ch4; high - nothing)
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
! kb      : real     
! selfref : real     
! forref  : real     
!
! absa    : real
! absb    : real
!-----------------------------------------------------------------

      integer , parameter :: ng16 = 2

      real  _cpus, dimension(ng16) :: fracrefb

      real  _cpus :: fracrefa(ng16,9)
      real  _cpusnp :: ka(9,5,13,ng16) ,absa(585,ng16)
      real  _cpusnp :: kb(5,13:59,ng16), absb(235,ng16)
      real  _cpus :: selfref(10,ng16)
      real  _cpus :: forref(4,ng16)

      real  _gpudeva , dimension(:) :: fracrefbd

      real  _gpudeva :: fracrefad(:,:)
      real  _gpudevanp ::  absad(:,:)
      real  _gpudevanp ::   absbd(:,:)
      real  _gpudeva :: selfrefd(:,:)
      real  _gpudeva :: forrefd(:,:)

      equivalence (ka(1,1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

            contains 

      subroutine copyToGPU16

      kaod = kao
      kbod = kbo
     dbcop( fracrefao , fracrefaod )
     

     !dbcop( kao , kaod ) 
     !dbcop( kbo , kbod )
     dbcop( selfrefo , selfrefod )
     dbcop( forrefo , forrefod )

     dbcop( fracrefa , fracrefad )
     dbcop( fracrefb , fracrefbd )

     !dbcop( ka , kad ) 
     !dbcop( kb , kbd )
     dbcopnp( absa , absad , 585, ng16)
     dbcopnp( absb , absbd , 235, ng16)
     dbcop( selfref , selfrefd )
     dbcop( forref , forrefd )

      end subroutine 

     subroutine reg16
 
     dbreg( fracrefao )
     

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

      end module rrlw_kg16

