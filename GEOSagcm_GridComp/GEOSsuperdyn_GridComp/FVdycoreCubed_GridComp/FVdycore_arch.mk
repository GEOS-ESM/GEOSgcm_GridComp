#
# System dependent FLAGS for FVdycore.
#


  ifeq ($(ESMA_FC), gfortran)

        USER_FFLAGS = -fcray-pointer

  endif

  ifeq ($(ESMA_FC), pgfortran)
        USER_FFLAGS = -DNO_QUAD_PRECISION
  endif
