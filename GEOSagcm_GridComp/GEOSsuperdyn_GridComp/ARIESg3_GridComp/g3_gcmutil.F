      subroutine getpl (p,pl,grid,im,jm)
C***********************************************************************
C
C  PURPOSE
C     Compute pressure at model levels
C
C     p ........ Dynamics Surface Pressure State Variable
C     pl ....... Pressure at model levels
C     grid ..... Dynamics Grid Structure
C     im ....... Dimension in x
C     jm ....... Dimension in y
C
C***********************************************************************
C*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
C***********************************************************************

c Declare Modules Used
c --------------------
      use g3_dynamics_state_module

      implicit none
c Declare Input Variables
c -----------------------
      type ( dynamics_grid_type ) grid

      integer  i,j,L
      integer  im,jm,lm
      real(kind=8)   p(im,jm)
      real(kind=8)  pl(im,jm,grid%lm)
      real(kind=8) ple(im,jm,grid%lm+1)

      lm = grid%lm

      call getple (p,ple,grid,im,jm)
      do L = 1,lm
      do j = 1,jm
      do i = 1,im
      pl(i,j,L) = ( ple(i,j,L)+ple(i,j,L+1) )*0.5
      enddo
      enddo
      enddo

      return
      end

      subroutine getpke (p,pke,grid,im,jm)
C***********************************************************************
C
C  PURPOSE
C     Compute pressure at model levels
C
C     p ........ Dynamics Surface Pressure State Variable
C     pke ...... Pressure**kappa at model edges
C     grid ..... Dynamics Grid Structure
C     im ....... Dimension in x
C     jm ....... Dimension in y
C
C***********************************************************************
C*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
C***********************************************************************

c Declare Modules Used
c --------------------
      use g3_dynamics_state_module

      implicit none
c Declare Input Variables
c -----------------------
      type ( dynamics_grid_type ) grid

      integer  i,j,L
      integer  im,jm,lm
      real(kind=8)   p(im,jm)
      real(kind=8) ple(im,jm,grid%lm+1)
      real(kind=8) pke(im,jm,grid%lm+1)
      real(kind=8) getcon,kappa

         lm = grid%lm
      kappa = getcon('KAPPA')

      call getple (p,ple,grid,im,jm)
      do L = 1,lm+1
      do j = 1,jm
      do i = 1,im
      pke(i,j,L) = ple(i,j,L)**kappa
      enddo
      enddo
      enddo

      return
      end

      subroutine getple (p,ple,grid,im,jm)
C***********************************************************************
C
C  PURPOSE
C     Compute pressure at model edge-levels
C
C     p ........ Dynamics Surface Pressure State Variable
C     ple ...... Pressure at model edge-levels
C     grid ..... Dynamics Grid Structure
C     im ....... Dimension in x
C     jm ....... Dimension in y
C
C***********************************************************************
C*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
C***********************************************************************

c Declare Modules Used
c --------------------
      use g3_dynamics_state_module

      implicit none
c Declare Input Variables
c -----------------------
      type ( dynamics_grid_type ) grid

      integer  i,j,L
      integer  im,jm,lm
      real(kind=8)   p(im,jm)
      real(kind=8) ple(im,jm,grid%lm+1)

      lm = grid%lm

      do L = 1,lm+1
      do j = 1,jm
      do i = 1,im
!     ple(i,j,L) = grid%alf(L) + grid%bet(L)*p(i,j)
      ple(i,j,L) = grid%sige(L)*p(i,j) + grid%ptop
      enddo
      enddo
      enddo

      return
      end

      subroutine getpk (p,pke,pk,grid,im,jm)
C***********************************************************************
C
C  PURPOSE
C     Compute pressure**kappa at model levels (Phillips formulation)
C
C     p ........ Dynamics Surface Pressure State Variable
C     pke ...... p**kappa at model edge-levels
C     pk ....... Pressure at model levels
C     grid ..... Dynamics Grid Structure
C     im ....... Dimension in x
C     jm ....... Dimension in y
C
C***********************************************************************
C*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
C***********************************************************************

c Declare Modules Used
c --------------------
      use g3_dynamics_state_module

      implicit none
c Declare Input Variables
c -----------------------
      type ( dynamics_grid_type ) grid

      integer   i,j,L
      integer   im,jm,lm
      real(kind=8)    p(im,jm)
      real(kind=8)   pk(im,jm,grid%lm)
      real(kind=8)  pke(im,jm,grid%lm+1)
      real(kind=8)  ple(im,jm,grid%lm+1)
      real(kind=8) pinv(im,jm)

      real(kind=8)    kappa,getcon,fact
      kappa = getcon('KAPPA')
      fact  = 1.0/( 1.0+kappa )

      lm = grid%lm

c Get ple at model edges
c ----------------------
      call getple (p,ple,grid,im,jm)

c Compute pk at model mid-levels
c ------------------------------
      do L = 1,lm
      do j = 1,jm
      do i = 1,im
      pk(i,j,L) = ( ple(i,j,L+1)*pke(i,j,L+1)-ple(i,j,L)*pke(i,j,L) )*fact
     .          / ( ple(i,j,L+1)             -ple(i,j,L)            )
      enddo
      enddo
      enddo

      return
      end

      subroutine harm (F,FNAME,L,IM)                                            
                                                                                
      PARAMETER ( FOUR  = 4.0 )                                                 
      PARAMETER ( ONE   = 1.0 )                                                 
      PARAMETER ( AHALF = 0.5 )                                                 
      PARAMETER ( ZERO  = 0.0 )                                                 
      PARAMETER ( P001  = 0.001 )                                               
      integer   IFAX (100)                                                       
      REAL(kind=8)      TRIGS(3*(IM+1))                                                
                                                                                
      real(kind=8)   F(IM), FZ(IM+2), C(2*(IM+2))                                    
      CHARACTER*4 FNAME                                                         
                                                                                
      PI   = FOUR * ATAN(ONE)                                                   
      PID2 = PI   * AHALF                                                       
      CALL FFTFAX (IM,IFAX,TRIGS)                                               
      WRITE(6,100) FNAME, L                                                     
                                                                                
      DO 10 I=1,IM                                                             
      FZ(I) = F(I)                                                              
   10 CONTINUE                                                                  
                                                                                
      CALL RFFTMLT (FZ,C,TRIGS,IFAX,1 ,IM+2,IM,1,-1)                           
                                                                                
      DO 20 I=1,IM+2                                                           
      C(I) = FZ(I)                                                              
   20 CONTINUE                                                                  
                                                                                
      WRITE(6,200) C(1)                                                         
      WRITE(6,300)                                                              
      KM  = IM/2                                                                
      KM1 = IM/2-1                                                              
      DO 30 K=1,KM1                                                             
      SIZE = FLOAT(IM)/FLOAT(K)                                                 
      K2   =  2*K+1                                                             
      AK   =  2*C(K2)                                                           
      BK   = -2*C(K2+1)                                                         
      AMPL = SQRT( AK**2 + BK**2 )                                              
                                                                                
      IF(AK.EQ.ZERO) TEMP = PID2 - ( ONE - SIGN(ONE,BK) )*PID2                  
      IF(AK.NE.ZERO)                                                            
     .TEMP = ATAN(BK/AK)                                                        
     .     + SIGN(ONE,BK) * (ONE-SIGN(ONE,AK)) * PID2                           
                                                                                
                            PHAZ = -TEMP                                        
      IF( AMPL.LT.P001) PHAZ =  ZERO                                            
      WRITE(6,400) K,SIZE,AK,BK,AMPL,PHAZ                                       
   30 CONTINUE                                                                  
                                                                                
      AK   =  C(IM+1)                                                         
      BK   =  ZERO                                                              
      PHAZ = -(ONE-SIGN(ONE,AK)) * PID2                                         
      AMPL = ABS(AK)                                                            
                                                                                
      IF( AMPL.LT.P001) PHAZ =  ZERO                                            
      SIZE = FLOAT(IM)/FLOAT(KM)                                                
      WRITE(6,400) KM,SIZE,AK,BK,AMPL,PHAZ                                      
                                                                                
  100 FORMAT(/,/,5X,'HARMONIC ANALYSIS FOR ',A4,' FIELD,  LEVEL = ',I2)
  200 FORMAT(/1X,'MEAN VALUE:',2X,F13.4)                                        
  300 FORMAT(/1X,'WAVENUMBER',3X,'WAVELENGTH',                                  
     $ 6X,' COS ',7X,' SIN ',7X,' AMPL ',5X,' PHASE ')                          
  400 FORMAT(4X,I3,6X,F8.1,6X,F9.3,3X,F9.3,3X,F9.3,3X,F9.3)                     
      RETURN                                                                    
      END                                                                       
