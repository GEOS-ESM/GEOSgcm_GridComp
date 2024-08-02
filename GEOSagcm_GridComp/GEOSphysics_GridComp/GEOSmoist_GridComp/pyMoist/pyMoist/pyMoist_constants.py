#Script containing all constants used for pyMoist

from ndsl.dsl.typing import Float

#Generic Global Constants
taufrz  =  450.0


# Following taken from MathConstants.F90
PI = Float(3.14159265358979323846)
deg_2_rad = PI / 180
rad_to_deg = 180 / PI

# Following taken from PhysicalConstants.F90
# Universal Constants
CODATA_2018_CONSTANTS = True # set for now, needs to be dynamic
if CODATA_2018_CONSTANTS == True:
    stfbol = Float(5.670374419E-8) # W/(m^2 K^4)
    avogad = Float(6.02214076E26) # 1/kmol
    runiv = Float(8314.462618) # J/(Kmole K)
else:
    stfbol = Float(5.6734E-8) # W/(m^2 K^4)
    avogad = Float(6.023E26) # 1/kmol
    runiv = Float(8314.47) # J/(Kmole K)

# Earth Constants
psdry = Float(98305) # Pa
seconds_per_sidereal_day = Float(86164.0) # s
gravity = Float(9.80665) # m^2/s
radius = Float(6371.0E3) # m
omega = 2 * PI / seconds_per_sidereal_day # 1/s
eath_eccentricity = Float(8.181919084262200E-2) # --
earth_semimajor_axis = Float(6378137) # m
km_per_deg = (1 / (radius/1000)) * rad_to_deg
deg_per_km = (radius / 1000) * deg_2_rad

# Physical Properties
h2omw =  Float(18.015) # kg/Kmole
o3mw = Float(47.9982) # kg/Kmole
latent_heat_vaporization = Float(2.4665E6) # J/kg @15C @1atm
latent_heat_fusion = Float(3.3370E5) # J/kg @1atm
latent_heat_sublimation = latent_heat_vaporization+latent_heat_fusion # J/kg

# Earth Specific Chemistry and Thermodynamic Constants
airmw = Float(28.965) # kg/Kmole
rdry = runiv/airmw # J/(kg K)
cpdry = 3.5 * rdry # J/(kg K)
cvdry = cpdry-rdry # J/(kg K)
rvap = runiv/h2omw # J/(kg K)
cpvap = 4 * rvap # J/(kg K)
cvvap = cpvap-rvap # J/(kg K)
kappa = rdry / cpdry # (2.0/7.0) # translator note no clue what 2/7 references
epsilon = h2omw / airmw # --
deltap = cpvap / cpdry # --
deltav = cvvap / cvdry # --
gammad = cpdry / cvdry # --
reference_pressure = Float(100000) # Pa
capice = Float(2000) # J/(K kg)
capwtr = Float(4218) # J/(K kg)
rhowtr = Float(1000) # kg/m^3
nuair = Float(1.533E-5) # m^2/s (@ 18C)
t_ice = Float(273.16) # K
surface_pressure = Float(98470) # Pa
karman = Float(0.40) # --
usmin = Float(1.00) # m/s
rho_seawater = Float(1026.0) # sea water density [kg/m^3]
rho_seaice = Float(917.0) # sea ice density [kg/m^3]
rho_snow = Float(330.0) # snow density [kg/m^3]
celsius_to_kelvin = Float(273.15) # K


# ice_fraction constants
# In anvil/convective clouds
aT_ICE_ALL = Float(252.16)
aT_ICE_MAX = Float(268.16)
aICEFRPWR  = Float(2.0)
# Over snow/ice SRF_TYPE = 2
iT_ICE_ALL = Float(236.16)
iT_ICE_MAX = Float(261.16)
iICEFRPWR  = Float(6.0)
# Over Land     SRF_TYPE = 1
lT_ICE_ALL = Float(239.16)
lT_ICE_MAX = Float(261.16)
lICEFRPWR  = Float(2.0)
# Over Oceans   SRF_TYPE = 0
oT_ICE_ALL = Float(238.16)
oT_ICE_MAX = Float(263.16)
oICEFRPWR  = Float(4.0)
# Jason method constants (translator note: I have no clue who jason is)
# In anvil/convective clouds
JaT_ICE_ALL = Float(245.16)
JaT_ICE_MAX = Float(261.16)
JaICEFRPWR  = Float(2.0)
# Over snow/ice
JiT_ICE_ALL = Float(t_ice-40.0)
JiT_ICE_MAX = Float(t_ice)
JiICEFRPWR  = Float(4.0)
# Other miscellaneous parameters
k_cond = Float(2.4E-2) # 1/(J*s*K)
diffu =  Float(2.2E-5) # m^2/s
taufrz  = Float(450)
dqcmax  = Float(1E-4)