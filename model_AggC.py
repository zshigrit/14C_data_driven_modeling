# model_AggC
# from model_parameters import *
import model_scalars as m_scal 

def millennial(C_maom,C_pom,C_agg,day,Vma=0.07,Kma=200,Amax=500,Vpa=0.002,Kpa=50,Kagg=0.0002):
    St, Sw, CUE = m_scal.millennial(day)
# inputs
    Fma = Vma * C_maom /(Kma + C_maom) * (1 - C_agg/Amax) * St * Sw
    Fpa = Vpa * (C_pom / (Kpa + C_pom)) * (1 - (C_agg/Amax)) * St * Sw 
# outputs
    Fagg = Kagg * C_agg * St * Sw 

    C_agg_local = C_agg + Fma + Fpa - Fagg
    return C_agg_local

# in Xiaofeng Xu github (https://github.com/email-clm/Millennial/blob/master/main.F90)
# f_PO_SO_agg = Vpom_agg * POM / (kpom_agg + POM) * (1. - SOILAGG / AGGmax) * t_scalar * w_scalar
# f_MI_SO_agg = Vmin_agg * MINERAL / (kmin_agg + MINERAL) * (1. - SOILAGG / AGGmax) !* t_scalar * w_scalar
# Fagg = Kagg * C_agg * St * Sw

# from model_parameters_xu import *
import model_scalars as m_scal 

def millennial_xu(C_maom,C_pom,C_agg,day,Vma=0.07,Kma=2000,Amax=500,Vpa=0.002,Kpa=50,Kagg=0.0002):
    St, Sw, CUE = m_scal.millennial_xu(day)
# inputs
    Fma = Vma * C_maom / (Kma + C_maom) * (1. - C_agg / Amax)
    Fpa = Vpa * C_pom / (Kpa + C_pom) * (1. - C_agg / Amax) * St * Sw 
    # Fpa = Vpa * (C_pom / (Kpa + C_pom)) * (1 - (C_agg/Amax)) * St * Sw 
# outputs
    Fagg = Kagg * C_agg * St * Sw 

    C_agg_local = C_agg + Fma + Fpa - Fagg
    return C_agg_local

#
#
# radiocarbon processes
def millennial_rad(C_maom,C_pom,C_agg,C_maom_rad,C_pom_rad,C_agg_rad,day,Vma=0.07,Kma=200,Amax=500,Vpa=0.002,Kpa=50,Kagg=0.0002,lambdax = 8267, Aabs=10**-12):
    
    St, Sw, CUE = m_scal.millennial(day)
    Asn_maom = C_maom_rad/C_maom
    Asn_pom = C_pom_rad/C_pom
    Asn_agg = C_agg_rad/C_agg

# inputs
    Fma = Vma * C_maom /(Kma + C_maom) * (1 - C_agg/Amax) * St * Sw
    Fpa = Vpa * (C_pom / (Kpa + C_pom)) * (1 - (C_agg/Amax)) * St * Sw 

    Fma_rad = Vma * C_maom /(Kma + C_maom) * (1 - C_agg/Amax) * St * Sw * Asn_maom
    Fpa_rad = Vpa * (C_pom / (Kpa + C_pom)) * (1 - (C_agg/Amax)) * St * Sw * Asn_pom
    # Fma_rad = Vma * C_maom_rad /(Kma + C_maom_rad) * (1 - C_agg_rad/Amax) * St * Sw
    # Fpa_rad = Vpa * (C_pom_rad / (Kpa + C_pom_rad)) * (1 - (C_agg_rad/Amax)) * St * Sw
# outputs
    Fagg = Kagg * C_agg * St * Sw
    Fagg_rad = Kagg * C_agg * St * Sw * Asn_agg 
    # Fagg_rad = Kagg * C_agg_rad * St * Sw   

    C_agg_local = C_agg + Fma + Fpa - Fagg
    C_agg_rad_local = C_agg_rad + Fma_rad + Fpa_rad - Fagg_rad - (1/(lambdax*365))*C_agg_rad

    # change unit to D14C per mil
    D14C_agg_local = ((C_agg_rad_local/C_agg_local)/Aabs-1)*1000

    return C_agg_local,C_agg_rad_local,D14C_agg_local