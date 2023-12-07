program test_moist_subroutines

    use GEOS_GF_InterfaceMod

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

    IM = 180
    JM = 180
    LM = 72

    ! Note : GF_Run parameters have been edited from original
    call GF_Run(IM, JM, LM, dirName, rank_str)

end program