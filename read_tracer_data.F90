module tracer_module
    use GEOSmoist_Process_Library

    implicit none

    private

    public read_tracers

    contains

subroutine read_tracers(IM, JM, LM, dirName, rank_str)
    integer, intent(in) :: IM, JM, LM
    character *100, intent(in) :: dirName, rank_str

    ! Note : This value is based on C180 run
    integer :: tracer_size = 40
    integer :: ii, fileID
    character * 2 :: tracer_str

    allocate(CNV_Tracers(tracer_size))

    do ii = 1, tracer_size
        if(ii < 10) then
            write(tracer_str,'(I1)') ii
        else
            write(tracer_str,'(I2)') ii
        endif
        allocate(CNV_Tracers(ii)%Q(IM, JM, LM))
        open(newunit=fileID, file=trim(dirName) // '/CNV_Tracers_' // trim(rank_str) // '_' // trim(tracer_str) //'_Q.in', status='old', form='unformatted', action='read')
        read(fileID) CNV_Tracers(ii)%Q
        close(fileID)
    enddo

end subroutine
end module