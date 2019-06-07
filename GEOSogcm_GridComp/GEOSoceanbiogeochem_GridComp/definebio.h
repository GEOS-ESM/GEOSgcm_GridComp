#define NNUT_DEFINED      4
#define NCHL_DEFINED      6
#define NZOO_DEFINED      1 
#define NDET_DEFINED 	  3 
#define NCAR_DEFINED 	  5 
      integer, parameter   :: nnut=NNUT_DEFINED
      integer, parameter   :: nchl=NCHL_DEFINED
      integer, parameter   :: nzoo=NZOO_DEFINED
      integer, parameter   :: ndet=NDET_DEFINED
      integer, parameter   :: ncar=NCAR_DEFINED
      integer, parameter   :: ntyp=nnut+nchl+nzoo+ndet+ncar
      integer, parameter   :: npe=nnut+nchl
      integer, parameter   :: nds=nnut+nchl+nzoo+1
      integer, parameter   :: nde=nds+(ndet-1)
      integer, parameter   :: ncs=ntyp-(ncar-1)

