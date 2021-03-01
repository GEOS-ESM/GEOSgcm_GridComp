#include "_gpudef.inc"
#include "_memdef.inc" 
 module rrlw_kg04

       !use parkind ,only : im => kind , rb => kind 

      use memory
 implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 4
! band 4:  630-700 cm-1 (low - h2o,co2; high - o3,co2)
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
        integer , parameter :: ng4  = 14
      integer , parameter :: no4  = 16

      real  _cpus :: kao(9,5,13,no4)
      real  _cpus :: kbo(5,5,13:59,no4)
      real  _cpusnp :: ka(9,5,13,ng4)   ,absa(585,ng4)
      real  _cpusnp :: kb(5,5,13:59,ng4),absb(1175,ng4)

      real  _cpus :: fracrefao(no4,9)  ,fracrefbo(no4,5)
 
      real  _cpus :: selfrefo(10,no4)  ,forrefo(4,no4)

      real  _gpudeva :: fracrefaod(:,:)  ,fracrefbod(:,:)
      !real  _gpudev :: kaod(9,5,13,no4)
      !real  _gpudev :: kbod(5,5,13:59,no4)
      real  _gpudeva :: selfrefod(:,:)  ,forrefod(:,:)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 4
! band 4:  630-700 cm-1 (low - h2o,co2; high - o3,co2)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
! absa    : real
! absb    : real
!fracrefa : real    
!fracrefb : real
! ka      : real     
! kb      : real     
! selfref : real     
! forref  : real     
!-----------------------------------------------------------------

      

      real  _cpus :: fracrefa(ng4,9)  ,fracrefb(ng4,5)
      
      real  _cpus :: selfref(10,ng4)  ,forref(4,ng4)

      real  _gpudeva :: fracrefad(:,:)  ,fracrefbd(:,:)
      real  _gpudevanp ::  absad(:,:)
      real  _gpudevanp ::  absbd(:,:)
      real  _gpudeva :: selfrefd(:,:)  ,forrefd(:,:)

      equivalence (ka(1,1,1,1),absa(1,1)),(kb(1,1,13,1),absb(1,1))

      contains 

      subroutine copyToGPU4

       dbcop( fracrefa , fracrefad )
       dbcop( fracrefb , fracrefbd )
      
       dbcopnp( absa , absad, 585, ng4 )
    
       dbcopnp( absb , absbd , 1175, ng4)
       dbcop( selfref , selfrefd )
       dbcop( forref , forrefd )

      end subroutine 

      subroutine reg4
       !33
       dbreg( fracrefa )
       dbreg( fracrefb )
    
       dbreg( absa )

       dbreg( absb )
       dbreg( selfref )
       dbreg( forref )

      end subroutine 

      end module rrlw_kg04
