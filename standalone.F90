program test_moist_subroutines

    use test_process_library_subroutines

    implicit none

    character*100 :: dirName, rank_str
    integer       :: IM, JM, LM, fileID, num_args

    num_args = command_argument_count()
    
    if(num_args.ne.2) then
        print*, 'Missing arguments : <executable> <data directory> <trim(rank_str)>'
        call exit(1)
    else
        call get_command_argument(1, dirName)
        call get_command_argument(2, rank_str)
    endif

    ! if (dirName(1:10) == './c24_data') then
    !     IM = 24
    !     JM = 24
    !     LM = 72
    ! elseif (dirName(1:10) == './c90_data') then
    !     IM = 90
    !     JM = 90
    !     LM = 72

    ! elseif (dirName(1:11) == './c180_data') then
        IM = 180
        JM = 180
        LM = 72
    ! endif

    ! if(trim(dirName(12:21)) == 'fillq2zero' .or. trim(dirName(13:22)) == 'fillq2zero') then
        call test_FILLQ2ZERO(IM, JM, LM, dirName, rank_str)
    ! endif
    
end program

! NASA Docket No. GSC-15,354-1, and identified as "GEOS-5 GCM Modeling Software”
  
! “Copyright © 2008 United States Government as represented by the Administrator
! of the National Aeronautics and Space Administration. All Rights Reserved.”
  
! Licensed under the Apache License, Version 2.0 (the "License"); you may not use
! this file except in compliance with the License. You may obtain a copy of the
! License at
  
! http://www.apache.org/licenses/LICENSE-2.0
  
! Unless required by applicable law or agreed to in writing, software distributed
! under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
! CONDITIONS OF ANY KIND, either express or implied. See the License for the
! specific language governing permissions and limitations under the License.