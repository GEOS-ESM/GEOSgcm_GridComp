
      diagName  = 'CPLoWGHT'
      diagTitle = 'grid-cell Ocean fraction from GEOS'
      diagUnits = '1               '
      diagCode  = 'SM      M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I    diagName, diagCode, diagUnits, diagTitle, 0, myThid )

#ifdef SEAICE_ITD

      diagName  = 'SItIces '
      diagTitle = 'Surface Temperature over Seaice for each category'
      diagUnits = 'K               '
      diagCode  = 'SM      MX      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      CALL DIAGNOSTICS_SETKLEV( diagName, nITD, myThid )

      diagName  = 'SIqIce  '
      diagTitle = 'SEAICE enthalpy for each layer and category'
      diagUnits = 'J/m^2           '
      diagCode  = 'SM      MX      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      numArea = nIceLayers*nITD
      CALL DIAGNOSTICS_SETKLEV( diagName, numArea, myThid )

      diagName  = 'SIqSnow '
      diagTitle = 'Snow enthalpy for each layer and category'
      diagUnits = 'J/m^2           '
      diagCode  = 'SM      MX      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      numArea = nSnowLayers*nITD
      CALL DIAGNOSTICS_SETKLEV( diagName, numArea, myThid )

      diagName  = 'SImeltPd'
      diagTitle = 'Melt Pond volume for each category'
      diagUnits = 'm               '
      diagCode  = 'SM      MX      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      CALL DIAGNOSTICS_SETKLEV( diagName, nITD, myThid )

      diagName  = 'SIiceAge'
      diagTitle = 'Seaice Age for each category'
      diagUnits = 's               '
      diagCode  = 'SM      MX      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      CALL DIAGNOSTICS_SETKLEV( diagName, nITD, myThid )

C-    Advection Increment:
      diagName  = 'SI_dArea'
      diagTitle = 'Seaice fraction Advection Increment per cat.'
      diagUnits = '1               '
      diagCode  = 'SM      MX      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      CALL DIAGNOSTICS_SETKLEV( diagName, nITD, myThid )

      diagName  = 'SI_dHeff'
      diagTitle = 'Seaice thickness Advection Increment per cat.'
      diagUnits = 'm               '
      diagCode  = 'SM      MX      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      CALL DIAGNOSTICS_SETKLEV( diagName, nITD, myThid )

      diagName  = 'SI_dHsnw'
      diagTitle = 'Snow thickness Advection Increment per cat.'
      diagUnits = 'm               '
      diagCode  = 'SM      MX      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      CALL DIAGNOSTICS_SETKLEV( diagName, nITD, myThid )

      diagName  = 'SI_dTIce'
      diagTitle = 'Seaice surf. temp. Advection Increment per cat.'
      diagUnits = 'K               '
      diagCode  = 'SM      MX      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      CALL DIAGNOSTICS_SETKLEV( diagName, nITD, myThid )

      diagName  = 'SI_dQIce'
      diagTitle = 'Seaice enthalpy Advect. Increment per layer and cat.'
      diagUnits = 'J/m^2           '
      diagCode  = 'SM      MX      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      numArea = nIceLayers*nITD
      CALL DIAGNOSTICS_SETKLEV( diagName, numArea, myThid )

      diagName  = 'SI_dQSnw'
      diagTitle = 'Snow enthalpy Advect. Increment per layer and cat.'
      diagUnits = 'J/m^2           '
      diagCode  = 'SM      MX      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      numArea = nSnowLayers*nITD
      CALL DIAGNOSTICS_SETKLEV( diagName, numArea, myThid )

      diagName  = 'SI_dMPnd'
      diagTitle = 'Melt Pond volume Advection Increment per cat.'
      diagUnits = 'm               '
      diagCode  = 'SM      MX      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      CALL DIAGNOSTICS_SETKLEV( diagName, nITD, myThid )

      diagName  = 'SI_dIcAg'
      diagTitle = 'Seaice Age Advection Increment per cat.'
      diagUnits = 's               '
      diagCode  = 'SM      MX      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      CALL DIAGNOSTICS_SETKLEV( diagName, nITD, myThid )

#endif /* SEAICE_ITD */
