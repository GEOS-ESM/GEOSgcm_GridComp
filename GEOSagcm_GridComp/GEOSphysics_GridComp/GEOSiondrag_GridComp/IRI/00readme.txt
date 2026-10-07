		International Reference Ionosphere Software 		----------------------------------------------------------

================================================================================

The IRI is a joined project of the Committee on Space Research (COSPAR) and the 
International Union of Radio Science (URSI). IRI is an empirical model specifying
monthly averages of electron density, ion composition, electron temperature, and 
ion temperature in the altitude range from 50 km to 1500 km.It also provides a 
few additional parameters: Total Electron Content (TEC) between a user specified 
lower and upper limit, equatorial vertical ion drift, occurrance probabilities for
spread-F and for the existance of an F1-layer, and auroral boundaries. 

This directory includes the FORTRAN program, coefficients, and indices files for 
the latest version of the International Reference Ionosphere model. This
version includes several options for different parts and parameters. A logical
array JF(50) is used to set these options. The IRI-recommended set of options
can be found in the COMMENT section at the beginning of IRISUB.FOR. IRITEST.FOR 
sets these options as the default.

The compilation/link command in Fortran 77 is:
f77 -o iri iritest.for irisub.for irifun.for iritec.for iridreg.for igrf.for 
 cira.for iriflip.for rocdrift.for
 
Directory Contents:
-----------------------------------------------------------------------------------

00_iri.zip	Compressed file that includes all files from this directory 
				 
irisub.for      This file contains the main subroutine iri_sub. It computes height  
                profiles of IRI output parameters (Ne, Te, Ti, Ni, vi) for specified 
                date and location. Also included is the subroutine iri_web that
                computes output parameters versus latitude, longitude (geog. or 
                geom.), year, day of year, hour (LT or UT), and month. An example  
                of how to call and use iri_web is shown in iritest.for.                

irifun.for      This file contains the subroutines and functions that are 
                required for running IRI.

iriflip.for     Subroutines for the FLIP-related new model for the bottomside
                ion composition of Richards et al., 2010 (doi:10.1029/2009RS004332)
                
iridreg.for     Subroutines for the D region models of Friedrich et al., 2018
		(doi:10.1029/2018JA025437)

iritec.for      This file includes the subroutines for computing the ionospheric 
                electron content between user-specified lower and upper limits. 

rocdrift.for    Equatorial vertical ion drift model of Fejer et al., 2008 
		(doi:10.1029/2007JA012801)

cira.for        This file includes the subroutines and functions for computing 
                the neutral temperature and densities from the COSPAR International 
		Reference Atmosphere (CIRA = NRLMSIS-00) of Picone et al., 2002 
		(doi:10.1029/2002JA009430) 

igrf.for        This file includes the subroutines for the International
                Geomagnetic Reference Field (IGRF) of Alken et al., 2021 
		(doi:10.1186/s40623-020-01288-x)

dgrf%%%%.dat    Definitive IGRF coefficients for the years 1945 to 2020 in steps
                of 5 years (%%%%=1945, 1950, etc.)(ASCII).
igrf%%%%.dat    Prelimenary IGRF coefficients for most recent IGRF year (ASCII).
igrf%%%%s.dat   IGRF coefficients for extrapolating 5 years into the future (ASCII).
 
MCSAT%%.dat     Monthly coefficient files for the Shubin (2015) 
		(doi:10.1016/j.asr.2015.05.029) COSMIC-based hmF2 model
                %%=month+10
                
iritest.for     Test program indicating how to use the IRI subroutines and what are the
		recommended defaults for the switches jf(50). 
                Compilation of iritest.for requires irisub.for, irifun.for, iritec.for, 
		iridreg.for, iriflip.for, igrf.for, rocdrift.for and cira.for.

irirtam.for     IRI subroutines that read the IRTAM coefficients and calculate NmF2, 
                hmF2, and B0 for the Real-Time IRI. The process of assimilating ionsonde 
                and Digisonde data into a Real-Time IRI is described in: Galkin et al.,
                Radio Sci., 47, RS0L07, doi:10.1029/2011RS004952, 2012.

irirtam-test.for  Test program for iritam.for

CCIR%%.dat	Monthly coefficient files for the CCIR foF2 and M(3000)F2 models
                (%% = month+10)

URSI%%.dat	Monthly coefficient files for the URSI foF2 model (%% = month+10)
                
IN ADDITION THE FOLLOWING INDICES FILES ARE REQUIRED THAT ARE NOT INCLUDED IN THIS 
DIRECTORY. THESE FILES ARE AVAILABLE ON THE IRI HOMEPAGE IRIMODEL.ORG:

ig_rz.dat       This file(s) contains the solar and ionospheric indices (IG12, Rz12) 
                for the time period from Jan 1958 onward in ASCII format. The file 
                is read by subroutine tcon that is included in irifun.for. 

apf107.dat      This file provides the 3-hour ap magnetic index and F10.7 daily
                81-day and 365-day index from 1960 onward (ASCII).

These indices files will be updated at close to quarterly intervals. Daily updates 
of these two files are available from the ECHAIM website (David Themens) at 
https://chain-new.chain-project.net/echaim_downloads/apf107.dat and
https://chain-new.chain-project.net/echaim_downloads/ig_rz.dat

The zipped file for IRI-2016 and earlier versions did not include the coefficients 
files CCIR%%.DAT and URSI%%.DAT. These files were made available at 
http://irimodel.org/COMMON_FILES/. Starting with IRI-2020 the files are no included
in the zipped file 00_iri.zip.

-----------------------------------------------------------------------------------
-----------------------------------------------------------------------------------

NOTE: Please consult the 'listing of changes' in the comments section at the top 
of each one of these programs for recent corrections and improvements.

More information about the IRI project can be found at irimodel.org including
references of articles that describe the model and access to sites that allow online 
computation and plotting of IRI parameters. 

The IRI output parameters are described at irimodel.org/IRI-output-arrays.docs
The available options are described at irimodel.org/IRI-Switches-options.docs
Answers to Frequently Asked Questions are available at irimodel.org/docs/IRI_FAQ.pdf

----------------------------- dbilitza@gmu.edu -------------------------------------
