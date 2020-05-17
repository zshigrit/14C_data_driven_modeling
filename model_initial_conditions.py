#%% initial conditions
Aabs = 10**-12
scalar_input = 0

C_agg0 = 100.
C_agg = 458.42 * scalar_input+1 # gC/m2
D14C_agg = -50.
C_agg_rad = ((D14C_agg/1000+1)*Aabs)*C_agg 

C_pom0 = 200.
C_pom = 27.82 * scalar_input+1
D14C_pom = 75.
C_pom_rad = ((D14C_pom/1000+1)*Aabs)*C_pom 

C_mic0 = 100.
C_mic = 34.46* scalar_input+1
D14C_mic = 10.
C_mic_rad = ((D14C_mic/1000+1)*Aabs)*C_mic

C_maom0 = 100.
C_maom = 2088.40* scalar_input+1
D14C_maom = -400.
C_maom_rad = ((D14C_maom/1000+1)*Aabs)*C_maom

C_lmwc0 = 100.
C_lmwc = 60.34* scalar_input+1
D14C_lmwc = 100.
C_lmwc_rad = ((D14C_lmwc/1000+1)*Aabs)*C_lmwc