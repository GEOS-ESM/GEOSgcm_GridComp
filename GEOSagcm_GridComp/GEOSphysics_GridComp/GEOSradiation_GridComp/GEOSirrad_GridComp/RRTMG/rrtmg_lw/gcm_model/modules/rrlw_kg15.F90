#include "_gpudef.inc"
#include "_memdef.inc" 
 module rrlw_kg15

       !use parkind ,only : im => kind , rb => kind 

      use memory
 implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 15
! band 15:  2380-2600 cm-1 (low - n2o,co2; high - nothing)
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
! kao_mn2 : real     
! selfrefo: real     
! forrefo : real     
!-----------------------------------------------------------------

      integer , parameter :: no15 = 16

      real  _cpus :: fracrefao(no15,9)
      real  _cpus :: kao(9,5,13,no15)
      real  _cpus :: kao_mn2(9,19,no15)
      real  _cpus :: selfrefo(10,no15)
      real  _cpus :: forrefo(4,no15)

      real  _gpudeva :: fracrefaod(:,:)
      real  _gpudev :: kaod(9,5,13,no15)
      real  _gpudeva :: kao_mn2d(:,:,:)
      real  _gpudeva :: selfrefod(:,:)
      real  _gpudeva :: forrefod(:,:)


!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 15
! band 15:  2380-2600 cm-1 (low - n2o,co2; high - nothing)
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
! ka_mn2  : real     
! selfref : real     
! forref  : real     
!
! absa    : real
!-----------------------------------------------------------------

      integer , parameter :: ng15 = 2

      real  _cpus :: fracrefa(ng15,9)
      real  _cpusnp :: ka(9,5,13,ng15) ,absa(585,ng15)
      real  _cpus :: ka_mn2(9,19,ng15)
      real  _cpus :: selfref(10,ng15)
      real  _cpus :: forref(4,ng15)

      real  _gpudeva :: fracrefad(:,:)
      real  _gpudevanp ::  absad(:,:)
      real  _gpudeva :: ka_mn2d(:,:,:)
      real  _gpudeva :: selfrefd(:,:)
      real  _gpudeva :: forrefd(:,:)


      equivalence (ka(1,1,1,1),absa(1,1))

      contains 

      subroutine copyToGPU15
      kaod = kao
     dbcop( fracrefao , fracrefaod )
     !dbcop( kao , kaod ) 
     dbcop( kao_mn2 , kao_mn2d )
     dbcop( selfrefo , selfrefod )
     dbcop( forrefo , forrefod )

     dbcop( fracrefa , fracrefad )
     !dbcop( ka , kad ) 
     dbcopnp( absa , absad, 585, ng15 )
     dbcop( ka_mn2 , ka_mn2d )
     dbcop( selfref , selfrefd )
     dbcop( forref , forrefd )

      end subroutine 

     subroutine reg15

     dbreg( fracrefao )
     !dbreg( kao ) 
     dbreg( kao_mn2 )
     dbreg( selfrefo )
     dbreg( forrefo )

     dbreg( fracrefa )
     !dbreg( ka ) 
     dbreg( absa )
     dbreg( ka_mn2 )
     dbreg( selfref )
     dbreg( forref )

      end subroutine 

      end module rrlw_kg15
