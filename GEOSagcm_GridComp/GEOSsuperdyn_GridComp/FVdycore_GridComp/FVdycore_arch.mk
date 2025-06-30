#
# System dependent FLAGS for FVdycore.
#


  ifeq ($(ESMA_FC), gfortran)
        USER_FFLAGS = -DNO_R16 -fcray-pointer
  endif

  ifeq ($(ESMA_FC), ftn)
        USER_FFLAGS = -DNO_R16
  endif

  ifeq ($(ESMA_FC), pgfortran)
        USER_FFLAGS = -DNO_R16
  endif

