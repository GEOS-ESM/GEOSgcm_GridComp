#
# System dependent FLAGS for FVdycore.
#


ifeq ($(ARCH),Linux)

# Intel Fortran Compiler (ifort or mpiifort)
# ------------------------------------------
  ifeq ($(subst mpi,,$(FC)), ifort)

       ifeq ($(IFORT_MAJOR), 8)
          USER_FFLAGS = -mp -stack_temps -fno-alias -ftz -auto
       else
       ifeq ($(IFORT_MAJOR), 9)
          USER_FFLAGS = -fp-model precise
       else
       ifeq ($(IFORT_MAJOR),10)
          USER_FFLAGS = -fno-inline-functions -assume protect_parens,minus0 -prec-div -prec-sqrt -no-ftz
       endif
       endif
       endif

  endif

  ifeq ($(ESMA_FC), gfortran)

        USER_FFLAGS = -DNO_R16 -fcray-pointer

  endif

  ifeq ($(ESMA_FC), ftn)
        USER_FFLAGS = -DNO_R16
  endif

  ifeq ($(ESMA_FC), pgfortran)
        USER_FFLAGS = -DNO_R16
  endif
endif
