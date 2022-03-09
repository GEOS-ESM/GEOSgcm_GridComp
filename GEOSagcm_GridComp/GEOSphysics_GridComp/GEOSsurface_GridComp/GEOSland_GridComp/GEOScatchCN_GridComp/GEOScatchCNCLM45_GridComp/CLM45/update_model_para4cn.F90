MODULE update_model_para4cn
  
  implicit none

  private
  INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE :: LocalTileID
  INTEGER, PUBLIC :: curr_year,curr_month,curr_day,curr_dofyr,curr_hour,curr_min,curr_sec 

  SAVE curr_year,curr_month,curr_day,curr_dofyr,curr_hour,curr_min,curr_sec, LocalTileID

  public :: upd_curr_date_time, upd_tileid

  contains

   ! ---------------------------------------

    subroutine upd_tileid (tileid)

      implicit none
      integer    :: NT
      integer, intent (in) :: tileid (:)

      NT = size (tileid)
      allocate (LocalTileID(1:NT))
      LocalTileID = tileid

    end subroutine upd_tileid

    ! ---------------------------------------
    
    subroutine upd_curr_date_time( year,month,day,dofyr,hour,min,sec )
      
      ! Return the current date_time.
      
      implicit none
      integer, intent(in) :: year,month,day,dofyr,hour,min,sec  
      
      curr_year  = year
      curr_month = month
      curr_day   = day
      curr_dofyr = dofyr
      curr_hour  = hour
      curr_min   = min
      curr_sec   = sec 
      
    end subroutine upd_curr_date_time
    
  end MODULE update_model_para4cn
