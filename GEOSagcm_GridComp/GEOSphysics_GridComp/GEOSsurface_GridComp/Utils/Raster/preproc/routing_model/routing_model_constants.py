# routing_model_constants.py

# Number of total river network nodes/segments
nc = 291284
# Number of reservoir nodes
nres = 7250
# Number of lake nodes
nlake = 3917
# Number of lat in 1m grid
nlat1m = 10800
# Number of lon in 1m grid
nlon1m = 21600
# M parameter in routing model
RRM_M = 0.45
# m parameter in routing model
RRM_mm = 0.35
# Liquid water density
RHO =1000. 
#Coefficient for hydropower calculation
fac_elec_a   = 0.30   
#Exponent for hydropower calculation
fac_elec_b   = 2.00   
#Coefficient for irrigation calculation (arid areas)
fac_irr_a    = 0.01   
#Scaling factor for irrigation (arid areas)
fac_irr_b    = 3.00   
#Coefficient for water supply calculation
fac_sup_a    = 0.03   
#Exponent for water supply calculation
fac_sup_b    = 2.00   
#Coefficient for other reservoir types
fac_other_a  = 0.20   
#Exponent for other reservoir types
fac_other_b  = 2.00  
#Coefficient for small lakes
fac_a_slake  = 0.003 
#Exponent for small lakes
fac_b_slake  = 0.40  
#Coefficient for large lakes
fac_a_llake  = 0.01  
#Exponent for large lakes
fac_b_llake  = 0.60   


