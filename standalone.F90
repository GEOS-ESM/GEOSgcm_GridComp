program test_moist_subroutines

    use test_gf_run

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
       
    ! if (trim(dirName(12:38)) == 'gfdl_cloud_microphys_driver' .or. trim(dirName(13:39)) == 'gfdl_cloud_microphys_driver') then
        call test_gf_run(IM, JM, LM, dirName, rank_str)
    ! endif

end program