      module RASPARAMS

      implicit none

! Define RASPARAM
! ---------------
      type RASPARAM_TYPE
           real               :: CUFRICFAC             ! 1
           real               :: SHR_LAMBDA_FAC        ! 2

           real               :: QC_CRIT_CN            ! 4
           real               :: RASAL1                ! 5
           real               :: RASAL2                ! 6
           real               :: RASNCL                ! 7
           real               :: LAMBDA_FAC            ! 8
           real               :: LAMBMX_FAC            ! 9
           real               :: MIN_DIAMETER          ! 10
           real               :: CUFRICLAMBDA          ! 11
           real               :: RDTLEXPON             ! 12
           real               :: STRAPPING             ! 13
           real               :: SDQV2                 ! 14
           real               :: SDQV3                 ! 15
           real               :: SDQVT1                ! 16
           real               :: ACRITFAC              ! 17
           real               :: HMINTRIGGER           ! 18
           real               :: LLDISAGGXP            ! 19
           real               :: PBLFRAC               ! 20
           real               :: RASAUTORAMPB          ! 21
           real               :: AUTOC_CN_ZDEP         ! 22
           real               :: MAXDALLOWED_S         ! 23
           real               :: MAXDALLOWED_D         ! 24
           real               :: MAXDALLOWED_E         ! 25
           real               :: RAS_RHMIN             ! 26
           real               :: RAS_RHFULL            ! 27
           real               :: CLDMICRO              ! 28
           real               :: FDROP_DUST            ! 29
           real               :: FDROP_SOOT            ! 30
           real               :: FDROP_SEASALT            ! 30
           real               :: RASAL_SLOPE
      endtype RASPARAM_TYPE

  end module RASPARAMS
