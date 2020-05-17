# model_POM
# from model_parameters import *
import model_scalars as m_scal
from model_input import Fi,rad_atm

def millennial (C_pom,C_agg,C_mic,day,Kagg=0.0002,pa=1./3,Vpa=0.002,Kpa=50,Amax=500,Vpl=10,Kpl=150,Kpe=12,pi=2./3):
    St, Sw, CUE = m_scal.millennial(day)
# inputs 
    Fin_l = Fi[day-(day//365)*365]
    # temporary: start
    # if day ==0:
    #     Fin_l = 172.
    # else:
    #     Fin_l = 0
    ########### end temporary
    Fagg = Kagg * C_agg * St * Sw
    Fap = Fagg * pa
# outputs
    Fpa = Vpa * (C_pom / (Kpa + C_pom)) * (1 - (C_agg/Amax)) * St * Sw
    Fpl = Vpl * (C_pom/(Kpl+C_pom)) * (C_mic / (Kpe + C_mic)) * St * Sw

    C_pom_local = C_pom + pi * Fin_l + Fap - Fpa - Fpl 
    return C_pom_local,Fin_l

# in Xiaofeng Xu github (https://github.com/email-clm/Millennial/blob/master/main.F90)
# f_PO_SO_agg = Vpom_agg * POM / (kpom_agg + POM) * (1. - SOILAGG / AGGmax) * t_scalar * w_scalar
# f_PO_LM_dep = Vpom_lmc * POM / (POM + kpom) * t_scalar * w_scalar !* (1. - MB / (MB + k_POMes))
def millennial_xu (C_pom,C_agg,C_mic,day,Kagg=0.0002,pa=1.5/3,Vpa=0.002,Kpa=50,Amax=500,Vpl=10,Kpl=150,pi=2./3):
    St, Sw, CUE = m_scal.millennial_xu(day)
# inputs 
    Fin_l = Fi[day-(day//365)*365]
    # temporary: start
    # if day ==0:
    #     Fin_l = 172.
    # else:
    #     Fin_l = 0
    ########### end temporary
    Fagg = Kagg * C_agg * St * Sw
    Fap = Fagg * pa
# outputs
    #Fpa = Vpa * (C_pom / (Kpa + C_pom)) * (1 - (C_agg/Amax)) * St * Sw
    Fpa = Vpa * C_pom / (Kpa + C_pom) * (1. - C_agg / Amax) * St * Sw 
    #Fpl = Vpl * (C_pom/(Kpl+C_pom)) * (C_mic / (Kpe + C_mic)) * St * Sw
    Fpl = Vpl * C_pom / (C_pom + Kpl) * St * Sw # !* (1. - MB / (MB + k_POMes)) 

    C_pom_local = C_pom + pi * Fin_l + Fap - Fpa - Fpl 
    return C_pom_local,Fin_l  


# radiocarbon processes
def millennial_rad (C_pom,C_agg,C_mic,C_pom_rad,C_agg_rad,C_mic_rad,day,years,Kagg=0.0002,pa=1./3,\
    Vpa=0.002,Kpa=50,Amax=500,Vpl=10,Kpl=150,Kpe=12,pi=2./3,lambdax = 8267, Aabs=10**-12):
    St, Sw, CUE = m_scal.millennial(day)
    Asn_pom = C_pom_rad/C_pom 
    Asn_agg = C_agg_rad/C_agg
    # Asn_mic = C_mic_rad/C_mic

# inputs 
    Fin_l = Fi[day-(day//365)*365]
    Fin_l_rad = Fin_l * ((rad_atm[-1*years+day//365]/1000+1)*Aabs)

    # temporary: start
    # if day ==0:
    #     Fin_l = 172.
    # else:
    #     Fin_l = 0
    ########### end temporary
    Fagg = Kagg * C_agg * St * Sw
    Fap = Fagg * pa

    Fagg_rad = Kagg * C_agg * St * Sw *Asn_agg
    Fap_rad = Fagg_rad * pa
    # Fagg_rad = Kagg * C_agg_rad * St * Sw
    # Fap_rad = Fagg_rad * pa
# outputs
    Fpa = Vpa * (C_pom / (Kpa + C_pom)) * (1 - (C_agg/Amax)) * St * Sw
    Fpl = Vpl * (C_pom/(Kpl+C_pom)) * (C_mic / (Kpe + C_mic)) * St * Sw

    Fpa_rad = Vpa * (C_pom / (Kpa + C_pom)) * (1 - (C_agg/Amax)) * St * Sw * Asn_pom 
    Fpl_rad = Vpl * (C_pom/(Kpl+C_pom)) * (C_mic / (Kpe + C_mic)) * St * Sw * Asn_pom 

    # Fpa_rad = Vpa * (C_pom_rad / (Kpa + C_pom_rad)) * (1 - (C_agg_rad/Amax)) * St * Sw
    # Fpl_rad = Vpl * (C_pom_rad/(Kpl+C_pom_rad)) * (C_mic_rad / (Kpe + C_mic_rad)) * St * Sw

    C_pom_local = C_pom + pi * Fin_l + Fap - Fpa - Fpl 
    C_pom_rad_local = C_pom_rad + pi * Fin_l_rad + Fap_rad - Fpa_rad - Fpl_rad- (1/(lambdax*365))*C_pom_rad
    
    # change unit to D14C per mil
    D14C_pom_local = ((C_pom_rad_local/C_pom_local)/Aabs-1)*1000
    return C_pom_local,Fin_l,C_pom_rad_local,D14C_pom_local