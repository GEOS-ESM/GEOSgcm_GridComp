from ndsl.dsl.typing import Float, Int

CP = Float(1004)  # J K-1 kg-1
XLV = Float(2.5e6)  # J kg-1
T_0 = Float(273.16)  # K
ccnclean = Float(250)

smaller_qv = Float(1.0e-16)  # kg/kg

MAXENS1 = Int(1)  # ensemble one on cap_max
MAXENS2 = Int(1)  # ensemble two on precip efficiency
MAXENS3 = Int(16)  # ensemble three done in cup_forcing_ens16 for G3d

deep = Int(1)
shallow = Int(2)
mid = Int(3)
