#model_LMWC
# from model_parameters import *
import math
import model_scalars as m_scal
from model_input import Fi,rad_atm

def millennial(C_pom,C_mic,C_maom,C_lmwc,day,Vpl=10,Kpl=150,Kpe=12,Vml=0.01,M_Lmin=10,Kml=25,\
    Kl=0.0015,Klm=0.25,c1=0.4833,c2=2.3282,clay_per=40,BulkD=1150,Vlm=0.35,Klb=7.2,pi=2./3.):
    St, Sw, CUE = m_scal.millennial(day)
    Qmax1 = 10**(c1 * math.log(clay_per,10) + c2)
    Qmax = BulkD * Qmax1 * 10**(-3) # convert to gC/m2 (flagged)
# inputs
    Fin_l = Fi[day-(day//365)*365]
    # temporary: start
    #if day ==0:
     #   Fin_l = 172.
    #else:
     #   Fin_l = 0
    ########### end temporary
    Fpl = Vpl * (C_pom/(Kpl+C_pom)) * (C_mic / (Kpe + C_mic)) * St * Sw
    Fml = Vml * (C_maom - M_Lmin) / (Kml + C_maom - M_Lmin) * St * Sw
# outputs
    Fl = Kl * C_lmwc * St * Sw # leaching
    Flm = C_lmwc * (Klm * Qmax * C_lmwc/(1+(Klm * C_lmwc)) - C_maom ) / Qmax * St *Sw 
    Flb = Vlm * C_lmwc * C_mic / (C_mic + Klb) * St * Sw 

    C_lmwc_local = C_lmwc + (1-pi) * Fin_l + Fpl + Fml - Fl - Flm - Flb
    return C_lmwc_local 

# in Xu github (https://github.com/email-clm/Millennial/blob/master/main.F90)
# litter input: forc_npp / 3
# f_PO_LM_dep = Vpom_lmc * POM / (POM + kpom) * t_scalar * w_scalar !* (1. - MB / (MB + k_POMes)) 
# f_MI_LM_des = Vm_l * (MINERAL - M_Lmin) / (km_l + MINERAL - M_Lmin) * t_scalar * w_scalar
# f_LM_leaching = LMWC * k_leaching * t_scalar
# f_LM_MI_sor = (temp / Qmax + 0.0015) * LMWC / 50. * t_scalar * w_scalar
# temp = (klmc_min * Qmax * LMWC ) / (2. + klmc_min * LMWC) - MINERAL
#f_LM_MB_uptake = LMWC * klmc * t_scalar * w_scalar * MB / (MB + kes) * LMWC / (20. + LMWC)

# from model_parameters_xu import *
#import model_scalars as m_scal
# from model_input import Fi

def millennial_xu(C_pom,C_mic,C_maom,C_lmwc,day,Vpl=10,Kpl=150,Vml=0.01,M_Lmin=10,Kml=25,\
    Kl=0.0015,Klmc_min=0.25,clay_per=40,Vlm=0.35,Klb=7.2,pi=2./3.):
    St, Sw, CUE = m_scal.millennial_xu(day)
    Qmax = 10.0 ** (0.297 * math.log(clay_per,10) + 2.355 + 0.50) #!* 1.25
# inputs
    Fin_l = Fi[day-(day//365)*365]
    # temporary: start
    #if day ==0:
     #   Fin_l = 172.
    #else:
     #   Fin_l = 0
    ########### end temporary
    # Fpl = Vpl * (C_pom/(Kpl+C_pom)) * (C_mic / (Kpe + C_mic)) * St * Sw
    Fpl = Vpl * C_pom / (C_pom + Kpl) * St * Sw # !* (1. - MB / (MB + k_POMes)) 
    Fml = Vml * (C_maom - M_Lmin) / (Kml + C_maom - M_Lmin) * St * Sw
# outputs
    # Fl = Kl * C_lmwc * St * Sw # leaching
    Fl = C_lmwc * Kl * St
    #Flm = C_lmwc * (Klm * Qmax * C_lmwc/(1+(Klm * C_lmwc)) - C_maom ) / Qmax * St *Sw 
    temp_Flm = (Klmc_min * Qmax * C_lmwc ) / (2. + Klmc_min * C_lmwc) - C_maom
    Flm = (temp_Flm / Qmax + 0.0015) * C_lmwc / 50. * St * Sw
    # Flb = Vlm * C_lmwc * C_mic / (C_mic + Klb) * St * Sw 
    Flb = C_lmwc * Vlm * St * Sw * C_mic / (C_mic + Klb) * C_lmwc / (20. + C_lmwc) 

    C_lmwc_local = C_lmwc + (1-pi) * Fin_l + Fpl + Fml - Fl - Flm - Flb
    return C_lmwc_local 

# radiocarbon processes
def millennial_rad(C_pom,C_mic,C_maom,C_lmwc,C_pom_rad,C_mic_rad,C_maom_rad,C_lmwc_rad,years,\
    day,Vpl=10,Kpl=150,Kpe=12,Vml=0.01,M_Lmin=10,Kml=25,lambdax = 8267, Aabs=10**-12,\
    Kl=0.0015,Klm=0.25,c1=0.4833,c2=2.3282,clay_per=40,BulkD=1150,Vlm=0.35,Klb=7.2,pi=2./3.):
    St, Sw, CUE = m_scal.millennial(day)
    Qmax1 = 10**(c1 * math.log(clay_per,10) + c2)
    Qmax = BulkD * Qmax1 * 10**(-3) # convert to gC/m2 (flagged)
    #
    Asn_pom = C_pom_rad/C_pom
    # Asn_mic = C_mic_rad/C_mic
    Asn_maom = C_maom_rad/C_maom
    Asn_lmwc = C_lmwc_rad/C_lmwc
    
# inputs
    Fin_l = Fi[day-(day//365)*365]
    Fin_l_rad = Fin_l * ((rad_atm[-1*years+day//365]/1000+1)*Aabs)
    rad_atm_record = rad_atm[-1*years+day//365]
    # temporary: start
    #if day ==0:
     #   Fin_l = 172.
    #else:
     #   Fin_l = 0
    ########### end temporary
    Fpl = Vpl * (C_pom/(Kpl+C_pom)) * (C_mic / (Kpe + C_mic)) * St * Sw
    Fml = Vml * (C_maom - M_Lmin) / (Kml + C_maom - M_Lmin) * St * Sw

    Fpl_rad = Vpl * (C_pom/(Kpl+C_pom)) * (C_mic / (Kpe + C_mic)) * St * Sw * Asn_pom
    Fml_rad = Vml * (C_maom - M_Lmin) / (Kml + C_maom - M_Lmin) * St * Sw * Asn_maom
    # Fpl_rad = Vpl * (C_pom_rad/(Kpl+C_pom_rad)) * (C_mic_rad / (Kpe + C_mic_rad)) * St * Sw
    # Fml_rad = Vml * (C_maom_rad - M_Lmin) / (Kml + C_maom_rad - M_Lmin) * St * Sw
# outputs
    Fl = Kl * C_lmwc * St * Sw # leaching
    Flm = C_lmwc * (Klm * Qmax * C_lmwc/(1+(Klm * C_lmwc)) - C_maom ) / Qmax * St *Sw 
    Flb = Vlm * C_lmwc * C_mic / (C_mic + Klb) * St * Sw

    Fl_rad = Kl * C_lmwc * St * Sw * Asn_lmwc # leaching
    Flm_rad = C_lmwc * (Klm * Qmax * C_lmwc/(1+(Klm * C_lmwc)) - C_maom ) / Qmax * St *Sw * Asn_lmwc 
    Flb_rad = Vlm * C_lmwc * C_mic / (C_mic + Klb) * St * Sw * Asn_lmwc
    # Fl_rad = Kl * C_lmwc_rad * St * Sw # leaching
    # Flm_rad = C_lmwc_rad * (Klm * Qmax * C_lmwc_rad/(1+(Klm * C_lmwc_rad)) - C_maom_rad ) / Qmax * St *Sw 
    # Flb_rad = Vlm * C_lmwc_rad * C_mic_rad / (C_mic_rad + Klb) * St * Sw 

    C_lmwc_local = C_lmwc + (1-pi) * Fin_l + Fpl + Fml - Fl - Flm - Flb
    C_lmwc_rad_local = C_lmwc_rad + (1-pi) * Fin_l_rad + Fpl_rad + Fml_rad - Fl_rad - Flm_rad - Flb_rad - (1/(lambdax*365))*C_lmwc_rad

    # change unit to D14C per mil
    D14C_lmwc_local = ((C_lmwc_rad_local/C_lmwc_local)/Aabs-1)*1000
    return C_lmwc_local,C_lmwc_rad_local,D14C_lmwc_local,rad_atm_record
