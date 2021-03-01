#include "_gpudef.inc"
#include "_memdef.inc" 
 module rrlw_kg09

       !use parkind ,only : im => kind , rb => kind 

      use memory
 implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 9
! band 9:  1180-1390 cm-1 (low - h2o,ch4; high - ch4)
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

      integer , parameter :: no9  = 16

      real  _cpus, dimension(no9) :: fracrefbo

      real  _cpus :: fracrefao(no9,9)
      real  _cpus :: kao(9,5,13,no9)
      real  _cpus :: kbo(5,13:59,no9)
      real  _cpus :: kao_mn2o(9,19,no9)
      real  _cpus :: kbo_mn2o(19,no9)
      real  _cpus :: selfrefo(10,no9)
      real  _cpus :: forrefo(4,no9)

      real  _gpudeva , dimension(:) :: fracrefbod

      real  _gpudeva :: fracrefaod(:,:)
      real  _gpudev :: kaod(9,5,13,no9)
      real  _gpudev :: kbod(5,13:59,no9)
      real  _gpudeva :: kao_mn2od(:,:,:)
      real  _gpudeva :: kbo_mn2od(:,:)
      real  _gpudeva :: selfrefod(:,:)
      real  _gpudeva :: forrefod(:,:)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 9
! band 9:  1180-1390 cm-1 (low - h2o,ch4; high - ch4)
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

      integer , parameter :: ng9  = 12

      real  _cpus, dimension(ng9) :: fracrefb
      real  _cpus :: fracrefa(ng9,9)
      real  _cpusnp :: ka(9,5,13,ng9) ,absa(585,ng9)
      real  _cpusnp :: kb(5,13:59,ng9) ,absb(235,ng9)
      real  _cpus :: ka_mn2o(9,19,ng9)
      real  _cpus :: kb_mn2o(19,ng9)
      real  _cpus :: selfref(10,ng9)
      real  _cpus :: forref(4,ng9)

      real  _gpudeva , dimension(:) :: fracrefbd
      real  _gpudeva :: fracrefad(:,:)
      real  _gpudevanp ::  absad(:,:)
      real  _gpudevanp ::  absbd(:,:)
      real  _gpudeva :: ka_mn2od(:,:,:)
      real  _gpudeva :: kb_mn2od(:,:)
      real  _gpudeva :: selfrefd(:,:)
      real  _gpudeva :: forrefd(:,:)

      equivalence (ka(1,1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

         contains 

      subroutine copyToGPU9
      kaod = kao
      kbod = kbo

     dbcop( fracrefao , fracrefaod )
     dbcop( fracrefbo , fracrefbod )


     dbcopnp( absa , absad , 585,ng9  )

     dbcopnp( absb , absbd, 235,ng9 )
     dbcop( kao_mn2o , kao_mn2od )
     dbcop( kbo_mn2o , kbo_mn2od )
     dbcop( selfref , selfrefd )
     dbcop( forref , forrefd )


     dbcop( fracrefa , fracrefad )
     dbcop( fracrefb , fracrefbd )

     dbcop( ka_mn2o , ka_mn2od )
     dbcop( kb_mn2o , kb_mn2od )
     dbcop( selfrefo , selfrefod )
     dbcop( forrefo , forrefod )


      end subroutine 

         subroutine reg9

    !105
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
      end module rrlw_kg09
