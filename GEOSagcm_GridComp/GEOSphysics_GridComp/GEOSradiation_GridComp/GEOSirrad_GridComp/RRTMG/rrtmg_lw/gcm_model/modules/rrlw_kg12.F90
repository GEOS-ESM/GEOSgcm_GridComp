#include "_gpudef.inc"
#include "_memdef.inc" 
 module rrlw_kg12

       !use parkind ,only : im => kind , rb => kind 

      use memory
 implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 12
! band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
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
! selfrefo: real     
! forrefo : real     
!-----------------------------------------------------------------

      integer , parameter :: no12 = 16

      real  _cpus :: fracrefao(no12,9)
      real  _cpus :: kao(9,5,13,no12)
      real  _cpus :: selfrefo(10,no12)
      real  _cpus :: forrefo(4,no12)

      real  _gpudeva :: fracrefaod(:,:)
      real  _gpudev :: kaod(9,5,13,no12) 
      real  _gpudeva :: selfrefod(:,:)
      real  _gpudeva :: forrefod(:,:)


!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 12
! band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
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
! selfref : real     
! forref  : real     
!
! absa    : real
!-----------------------------------------------------------------

      integer , parameter :: ng12 = 8

      real  _cpus :: fracrefa(ng12,9)
      real  _cpusnp :: ka(9,5,13,ng12) ,absa(585,ng12)
      real  _cpus :: selfref(10,ng12)
      real  _cpus :: forref(4,ng12)

      real  _gpudeva :: fracrefad(:,:)
      real  _gpudevanp ::  absad(:,:)
      real  _gpudeva :: selfrefd(:,:)
      real  _gpudeva :: forrefd(:,:)

      equivalence (ka(1,1,1,1),absa(1,1))

      contains 

      subroutine copyToGPU12
      kao = kaod
     dbcop( fracrefao , fracrefaod )
     !dbcop( kao , kaod ) 
     dbcop( selfrefo , selfrefod )
     dbcop( forrefo , forrefod )



     dbcop( fracrefa , fracrefad )
     !dbcop( ka , kad ) 
     dbcopnp( absa , absad, 585,ng12 )
     
     dbcop( selfref , selfrefd )
     dbcop( forref , forrefd )

      end subroutine

     subroutine reg12

     dbreg( fracrefao )
     !dbreg( kao ) 
     dbreg( selfrefo )
     dbreg( forrefo )



     dbreg( fracrefa )
     !dbreg( ka ) 
     dbreg( absa )
     
     dbreg( selfref )
     dbreg( forref )

      end subroutine

      end module rrlw_kg12
