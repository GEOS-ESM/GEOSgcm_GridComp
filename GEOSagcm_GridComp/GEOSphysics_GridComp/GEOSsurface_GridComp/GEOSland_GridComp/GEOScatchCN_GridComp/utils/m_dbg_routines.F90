MODULE DBG_ROUTINES
  USE catch_constants, ONLY:                  &
       N_SNOW         => CATCH_N_SNOW,        &
       N_GT           => CATCH_N_GT     
  implicit none
  private 
  public :: c_inputs,c_params,c_updates,read_catch_inputs,read_catch_params,read_catch_updates

  type :: c_inputs

      real :: PCU       
      real :: PLS       
      real :: SNO       
      real :: UUU       
      real :: EVSBTS    
      real :: DEVSBTS   
      real :: TILEZERO  
      real :: SHSBTS    
      real :: DSHSBTS   
      real :: EVSBTT    
      real :: DEVSBTT   
      real :: SHSBTT    
      real :: DSHSBTT   
      real :: EVSBTW    
      real :: DEVSBTW   
      real :: SHSBTW    
      real :: DSHSBTW   
      real :: EVSBTSN   
      real :: DEVSBTSN  
      real :: SHSBTSN   
      real :: DSHSBTSN  
      real :: TA        
      real :: QA        
      real :: RAS       
      real :: RAT       
      real :: RAW       
      real :: RASN      
      real :: ZTH       
      real :: DRPAR     
      real :: DFPAR     
      real :: SWNETFREE 
      real :: SWNETSNOW 
      real :: LWDNSRF   
      real :: PS        
      real :: LAI0      
      real :: GRN0      
      real :: Z2CH      
      real :: SQSCAT    
      real :: RSL1      
      real :: RSL2      
      real :: RDC      
      real :: QSATS    
      real :: DQSS     
      real :: ALWX     
      real :: BLWX     
      real :: QSATT    
      real :: DQST     
      real :: QSATW    
      real :: DQSW     
      real :: QSATSN   
      real :: DQSSN    

   end type c_inputs

   type :: c_params

      integer :: VEG     
      real :: BF1     
      real :: BF2     
      real :: BF3     
      real :: VGWMAX  
      real :: CDCR1   
      real :: CDCR2   
      real :: PSIS    
      real :: BEE     
      real :: POROS   
      real :: WPWET   
      real :: COND    
      real :: GNU     
      real :: ARS1    
      real :: ARS2    
      real :: ARS3    
      real :: ARA1    
      real :: ARA2    
      real :: ARA3    
      real :: ARA4    
      real :: ARW1    
      real :: ARW2    
      real :: ARW3    
      real :: ARW4    
      real :: TSA1    
      real :: TSA2    
      real :: TSB1    
      real :: TSB2    
      real :: ATAU    
      real :: BTAU     

   end type c_params

   type :: c_updates

      real :: TC1
      real :: TC2
      real :: TC4
      real :: QC1
      real :: QC2
      real :: QC4
      real :: CAPAC
      real :: CATDEF
      real :: RZEXC
      real :: SRFEXC
      real :: GHTCNT(N_gt)
      real :: TSURF
      real :: WESNN (N_snow)
      real :: HTSNNN(N_snow)
      real :: SNDZN (N_snow)

   end type c_updates

   contains
!
! -------------------------------------------------------------------
!
   subroutine read_catch_params (unit,ntiles,catchin)
     
     implicit none
     
     integer, intent (in) :: ntiles,unit
     integer :: n
     type(c_params), dimension (ntiles), intent(inout) :: catchin
        read (unit) catchin%VEG
        read (unit) catchin%BF1     
        read (unit) catchin%BF2     
        read (unit) catchin%BF3     
        read (unit) catchin%VGWMAX  
        read (unit) catchin%CDCR1   
        read (unit) catchin%CDCR2   
        read (unit) catchin%PSIS    
        read (unit) catchin%BEE     
        read (unit) catchin%POROS   
        read (unit) catchin%WPWET   
        read (unit) catchin%COND    
        read (unit) catchin%GNU     
        read (unit) catchin%ARS1    
        read (unit) catchin%ARS2    
        read (unit) catchin%ARS3    
        read (unit) catchin%ARA1    
        read (unit) catchin%ARA2    
        read (unit) catchin%ARA3    
        read (unit) catchin%ARA4    
        read (unit) catchin%ARW1    
        read (unit) catchin%ARW2    
        read (unit) catchin%ARW3    
        read (unit) catchin%ARW4    
        read (unit) catchin%TSA1    
        read (unit) catchin%TSA2    
        read (unit) catchin%TSB1    
        read (unit) catchin%TSB2    
        read (unit) catchin%ATAU    
        read (unit) catchin%BTAU    

      end subroutine read_catch_params

!
! ------------------------------------------------------
!
!
! -------------------------------------------------------------------
!
subroutine read_catch_updates (unit,ntiles,catchin)

implicit none

integer, intent (in) :: ntiles,unit
type(c_updates), dimension (ntiles), intent(inout) :: catchin

        read (unit) catchin%TC1
        read (unit) catchin%TC2
        read (unit) catchin%TC4
        read (unit) catchin%QC1
        read (unit) catchin%QC2
        read (unit) catchin%QC4
        read (unit) catchin%CAPAC
        read (unit) catchin%CATDEF
        read (unit) catchin%RZEXC
        read (unit) catchin%SRFEXC
        read (unit) catchin%GHTCNT(1)
        read (unit) catchin%GHTCNT(2)
        read (unit) catchin%GHTCNT(3)
        read (unit) catchin%GHTCNT(4)
        read (unit) catchin%GHTCNT(5)
        read (unit) catchin%GHTCNT(6)
        read (unit) catchin%TSURF
        read (unit) catchin%WESNN(1)
        read (unit) catchin%WESNN(2)
        read (unit) catchin%WESNN(3)
        read (unit) catchin%HTSNNN(1)
        read (unit) catchin%HTSNNN(2)
        read (unit) catchin%HTSNNN(3)
        read (unit) catchin%SNDZN(1)
        read (unit) catchin%SNDZN(2)
        read (unit) catchin%SNDZN(3)

      END subroutine read_catch_updates

!
! ------------------------------------------------------
!
      subroutine read_catch_inputs (unit,ntiles,catchin)

        implicit none
        
        integer, intent (in) :: ntiles,unit
        type(c_inputs), dimension (ntiles), intent(inout) :: catchin
        
        read (unit) catchin%PCU      
        read (unit) catchin%PLS      
        read (unit) catchin%SNO      
        read (unit) catchin%UUU      
        read (unit) catchin%EVSBTS   
        read (unit) catchin%DEVSBTS  
        read (unit) catchin%TILEZERO 
        read (unit) catchin%SHSBTS   
        read (unit) catchin%TILEZERO 
        read (unit) catchin%DSHSBTS  
        read (unit) catchin%EVSBTT   
        read (unit) catchin%DEVSBTT  
        read (unit) catchin%TILEZERO 
        read (unit) catchin%SHSBTT   
        read (unit) catchin%TILEZERO 
        read (unit) catchin%DSHSBTT  
        read (unit) catchin%EVSBTW   
        read (unit) catchin%DEVSBTW  
        read (unit) catchin%TILEZERO 
        read (unit) catchin%SHSBTW   
        read (unit) catchin%TILEZERO 
        read (unit) catchin%DSHSBTW  
        read (unit) catchin%EVSBTSN  
        read (unit) catchin%DEVSBTSN 
        read (unit) catchin%TILEZERO 
        read (unit) catchin%SHSBTSN  
        read (unit) catchin%TILEZERO 
        read (unit) catchin%DSHSBTSN 
        read (unit) catchin%TA       
        read (unit) catchin%QA       
        read (unit) catchin%RAS      
        read (unit) catchin%RAT      
        read (unit) catchin%RAW      
        read (unit) catchin%RASN     
        read (unit) catchin%ZTH      
        read (unit) catchin%DRPAR    
        read (unit) catchin%DFPAR    
        read (unit) catchin%SWNETFREE
        read (unit) catchin%SWNETSNOW
        read (unit) catchin%LWDNSRF  
        read (unit) catchin%PS       
        read (unit) catchin%LAI0     
        read (unit) catchin%GRN0     
        read (unit) catchin%Z2CH     
        read (unit) catchin%SQSCAT   
        read (unit) catchin%RSL1     
        read (unit) catchin%RSL2     
        read (unit) catchin%RDC      
        read (unit) catchin%QSATS    
        read (unit) catchin%DQSS     
        read (unit) catchin%ALWX     
        read (unit) catchin%BLWX     
        read (unit) catchin%QSATT    
        read (unit) catchin%DQST     
        read (unit) catchin%ALWX     
        read (unit) catchin%BLWX     
        read (unit) catchin%QSATW    
        read (unit) catchin%DQSW     
        read (unit) catchin%ALWX     
        read (unit) catchin%BLWX     
        read (unit) catchin%QSATSN   
        read (unit) catchin%DQSSN    
        read (unit) catchin%ALWX     
        read (unit) catchin%BLWX     

      end subroutine read_catch_inputs

    END MODULE DBG_ROUTINES
