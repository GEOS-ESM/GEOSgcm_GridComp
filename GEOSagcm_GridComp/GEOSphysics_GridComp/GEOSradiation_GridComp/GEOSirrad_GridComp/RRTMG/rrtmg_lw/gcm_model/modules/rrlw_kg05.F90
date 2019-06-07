#include "_gpudef.inc"
#include "_memdef.inc" 
 module rrlw_kg05

       !use parkind ,only : im => kind , rb => kind 

      use memory
 implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 5
! band 5:  700-820 cm-1 (low - h2o,co2; high - o3,co2)
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
! kao_mo3 : real     
! selfrefo: real     
! forrefo : real     
! ccl4o   : real
!-----------------------------------------------------------------

      integer , parameter :: no5  = 16
      integer , parameter :: ng5  = 16
      real  _cpusnp :: ka(9,5,13,ng5),kb(5,5,13:59,ng5)  
      real  _cpus :: kao(9,5,13,no5)
      real  _cpus :: kbo(5,5,13:59,no5)

      real  _cpus :: fracrefao(no5,9) ,fracrefbo(no5,5) 
      real  _cpusnp :: absa(585,ng5)
 
      real  _cpus :: kao_mo3(9,19,no5)
      real  _cpus :: selfrefo(10,no5)
      real  _cpus :: forrefo(4,no5)
      real  _cpus :: ccl4o(no5)


      real  _gpudeva :: fracrefaod(:,:) ,fracrefbod(:,:)
      real  _gpudev :: kaod(9,5,13,no5)
      real  _gpudev :: kbod(5,5,13:59,no5)
      real  _gpudeva :: kao_mo3d(:,:,:)
      real  _gpudeva :: selfrefod(:,:)
      real  _gpudeva :: forrefod(:,:)
      real  _gpudeva :: ccl4od(:)
!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 5
! band 5:  700-820 cm-1 (low - h2o,co2; high - o3,co2)
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
! ka_mo3  : real     
! selfref : real     
! forref  : real     
! ccl4    : real
!
! absa    : real
! absb    : real
!-----------------------------------------------------------------


 
      real  _cpusnp :: absb(1175,ng5)

      real  _cpus :: fracrefa(ng5,9) ,fracrefb(ng5,5)
      
      real  _cpus :: ka_mo3(9,19,ng5)
      real  _cpus :: selfref(10,ng5)
      real  _cpus :: forref(4,ng5)
      real  _cpus :: ccl4(ng5)

      real  _gpudeva :: fracrefad(:,:) ,fracrefbd(:,:)
      real  _gpudevanp ::  absad(:,:)
      real  _gpudevanp ::  absbd(:,:)
      real  _gpudeva :: ka_mo3d(:,:,:)
      real  _gpudeva :: selfrefd(:,:)
      real  _gpudeva :: forrefd(:,:)
      real  _gpudeva :: ccl4d(:)
      
      equivalence (ka(1,1,1,1),absa(1,1)),(kb(1,1,13,1),absb(1,1))

      contains 

      subroutine copyToGPU5

       dbcop( fracrefao , fracrefaod )
       dbcop( fracrefbo , fracrefbod )
    
       dbcop( kao_mo3 , kao_mo3d )
       dbcop( selfrefo , selfrefod )
       dbcop( forrefo , forrefod )
       dbcop( ccl4o , ccl4od )

       dbcop( fracrefa , fracrefad )
       dbcop( fracrefb , fracrefbd )

       dbcopnp( absa , absad, 585, ng5 )
       dbcopnp( absb , absbd, 1175, ng5 )

       dbcop( ka_mo3 , ka_mo3d )
       dbcop( selfref , selfrefd )
       dbcop( forref , forrefd )
       dbcop( ccl4 , ccl4d )

      end subroutine 

      subroutine reg5
    
       dbreg( fracrefao )
       dbreg( fracrefbo )
     
       dbreg( kao_mo3 )
       dbreg( selfrefo )
       dbreg( forrefo )
       dbreg( ccl4o )

       dbreg( fracrefa )
       dbreg( fracrefb )
      
       dbreg( absa )
     
       dbreg( absb )
       dbreg( ka_mo3 )
       dbreg( selfref )
       dbreg( forref )
       dbreg( ccl4 )

      end subroutine 

      end module rrlw_kg05

