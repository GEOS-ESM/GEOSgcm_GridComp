CBOP
C !ROUTINE: SEAICE_LAYERS.h

C !DESCRIPTION: \bv
C     *==========================================================*
C     | SEAICE_LAYERS.h
C     | o header file for sea ice multi-layer variables
C     |   for refined thermodynamics with vertical discretisation
C     |   of seaice and snow cover, accounting for heat content
C     |   and brine pockets contribution to seaice enthalpy
C     *==========================================================*
C     | Note: for now, only used when coupled to to GEOS AGCM
C     *==========================================================*
C \ev
CEOP

#ifdef HACK_FOR_GMAO_CPL
C     SIqIce       :: Seaice enthalpy for each ice layer and each category [J/m^2]
C     SIqSnow      :: Snow  enthalpy for each snow layer and each category [J/m^2]
C     SImeltPd     :: Melt Pond volume for each category [m]
C     SIiceAge     :: Seaice Age       for each category [s]
C     SIskinS      :: seaice skin salinity [psu]
C     SIskinH      :: seaice skin-layer depth [m]
      _RL SIqIce  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nIceLayers,
     &             nITD,nSx,nSy)
      _RL SIqSnow (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSnowLayers,
     &             nITD,nSx,nSy)
      _RL SImeltPd(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nITD,nSx,nSy)
      _RL SIiceAge(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nITD,nSx,nSy)
      _RL SIskinS (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL SIskinH (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /SEAICE_LAYERS/
     &        SIqIce,   SIqSnow,
     &        SImeltPd, SIiceAge,
     &        SIskinS,  SIskinH

C     SIwindTauX   :: wind stress over seaice, X-component at U point
C     SIwindTauY   :: wind stress over seaice, Y-component at V point
      _RL SIwindTauX(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL SIwindTauY(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /SEAICE_DYN_FORCING/
     &        SIwindTauX, SIwindTauY

C     SIadv_Area    :: advection increment of Seaice fraction  [-]
C     SIadv_Heff    :: advection increment of Seaice thickness [m]
C     SIadv_Hsnow   :: advection increment of snow thickness   [m]
C     SIadv_tIces   :: advection increment of ice surface temperature
C     SIadv_qIce    :: advection increment of Seaice enthalpy [J/m^2]
C     SIadv_qSnow   :: advection increment of Snow  enthalpy  [J/m^2]
C     SIadv_meltPd  :: advection increment of Melt Pond volume [m]
C     SIadv_iceAge  :: advection increment of Seaice Age       [s]
C     SIadv_skinS   :: advection increment of seaice skin salinity [psu]
C     SIadv_skinH   :: advection increment of seaice skin-layer depth [m]
      COMMON /SEAICE_ADV_INCREMENT/
     &        SIadv_Area, SIadv_Heff, SIadv_Hsnow,
     &        SIadv_tIces, SIadv_qIce, SIadv_qSnow,
     &        SIadv_meltPd, SIadv_iceAge, SIadv_skinS, SIadv_skinH
      _RL SIadv_Area  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nITD,nSx,nSy)
      _RL SIadv_Heff  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nITD,nSx,nSy)
      _RL SIadv_Hsnow (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nITD,nSx,nSy)
      _RL SIadv_tIces (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nITD,nSx,nSy)
      _RL SIadv_qIce  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nIceLayers,
     &                 nITD,nSx,nSy)
      _RL SIadv_qSnow (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSnowLayers,
     &                 nITD,nSx,nSy)
      _RL SIadv_meltPd(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nITD,nSx,nSy)
      _RL SIadv_iceAge(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nITD,nSx,nSy)
      _RL SIadv_skinS (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL SIadv_skinH (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

C     oceWeight    :: grid-cell ocean fraction from GEOS [0-1]
      _RL oceWeight(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /DIAGS_GMAO_CPL/
     &        oceWeight

C     SI_FRZMLT    :: available heat (W/m^2) to freeze (>0) or melt (<0) sea ice
C                     so that surface level ocean reaches freezing temperature
      _RL SI_FRZMLT(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /SEAICE_FRZMLT/
     &        SI_FRZMLT

#endif /* HACK_FOR_GMAO_CPL */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
