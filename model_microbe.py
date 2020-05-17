#model_microbe
#def millennial(C_lmwc,C_mic,St,Sw):
# from model_parameters import *
#from model_input import *
#from model_scalars import *
# import model_parameters as m_par
import model_scalars as m_scal 
# test_list = [Vlm,Klb,Kmic,Kmm]
#print(test_list)
# print('test')
# def millennial(C_lmwc,C_mic,day):
def millennial(C_lmwc,C_mic,day,Vlm=0.35,Klb=7.2,Kmic=0.036,Kmm=0.025):
    St, Sw, CUE = m_scal.millennial(day)
# input
    Flb = Vlm * C_lmwc * C_mic / (C_mic + Klb) * St * Sw 
# output
    Fbr = Kmic * C_mic * St * Sw # microbial mortality
    Fbm = Kmm * C_mic * St* Sw # adsorption of mineral surface
    #
    #CUE = 0.4
    C_mic_local = C_mic + Flb * CUE - Fbm - Fbr
    F_res = Fbr+Flb*(1-CUE)
    return C_mic_local,F_res  
    #return Flb, Fbr, Fbm

# in Xiaofeng Xu github (https://github.com/email-clm/Millennial/blob/master/main.F90)
#f_LM_MB_uptake = LMWC * klmc * t_scalar * w_scalar * MB / (MB + kes) * LMWC / (20. + LMWC)
#f_MB_MI_sor = MB * kmic * 0.15 * t_scalar_mb * w_scalar  !* (MB / 200) * (MB / 200)
#f_MB_atm = temp2 + MB * kmic * t_scalar_mb * w_scalar
#temp2 = f_LM_MB_uptake * (1. - (CUEref + CUET * (forc_st - Taeref)))

#model_microbe from Xu github website

# from model_parameters_xu import *
# import model_parameters as m_par
import model_scalars as m_scal 


def millennial_xu(C_lmwc,C_mic,day,Vlm=0.35,Klb=7.2,Kmic=0.025):
    St, Sw, CUE = m_scal.millennial_xu(day)
# input
    # Flb = Vlm * C_lmwc * C_mic / (C_mic + Klb) * St * Sw
    Flb = C_lmwc * Vlm * St * Sw * C_mic / (C_mic + Klb) * C_lmwc / (20. + C_lmwc) 
# output
    # Fbr = Kmic * C_mic * St * Sw # microbial mortality
    Fbr= C_mic * Kmic * St * Sw
    #Fbm = Kmm * C_mic * St* Sw # adsorption of mineral surface
    Fbm = C_mic * Kmic * 0.15 * St * Sw  
    #
    #CUE = 0.4
    C_mic_local = C_mic + Flb * CUE - Fbm - Fbr

    F_res = Fbr+Flb*(1-CUE)
    return C_mic_local,F_res   
    #return Flb, Fbr, Fbm
#
#
# raidocarbon processes
def millennial_rad(C_lmwc,C_mic,C_lmwc_rad,C_mic_rad,day,Vlm=0.35,Klb=7.2,Kmic=0.036,Kmm=0.025,lambdax = 8267, Aabs=10**-12):
    
    St, Sw, CUE = m_scal.millennial(day)
   
    Asn_lmwc = C_lmwc_rad/C_lmwc
    Asn_mic = C_mic_rad/C_mic

# input
    Flb = Vlm * C_lmwc * C_mic / (C_mic + Klb) * St * Sw 
    Flb_rad = Vlm * C_lmwc * C_mic / (C_mic + Klb) * St * Sw * Asn_lmwc 
    # Flb_rad = Vlm * C_lmwc_rad * C_mic_rad / (C_mic_rad + Klb) * St * Sw 
# output
    Fbr = Kmic * C_mic * St * Sw # microbial mortality
    Fbm = Kmm * C_mic * St* Sw # adsorption of mineral surface

    Fbr_rad = Kmic * C_mic * St * Sw * Asn_mic # microbial mortality
    Fbm_rad = Kmm * C_mic * St* Sw * Asn_mic # adsorption of mineral surface
    # Fbr_rad = Kmic * C_mic_rad * St * Sw # microbial mortality
    # Fbm_rad = Kmm * C_mic_rad * St* Sw  # adsorption of mineral surface
    #
    C_mic_local = C_mic + Flb * CUE - Fbm - Fbr
    F_res = Fbr + Flb * (1-CUE)
    
    # radiocarbon
    C_mic_rad_local = C_mic_rad + Flb_rad * CUE - Fbm_rad - Fbr_rad - (1/(lambdax*365)) * C_mic_rad
    F_res_rad = Fbr_rad+Flb_rad*(1-CUE)
    #
    # change unit to D14C per mil
    D14C_mic_local = ((C_mic_rad_local/C_mic_local)/Aabs-1)*1000
    D14C_res_local = ((F_res_rad/F_res)/Aabs-1)*1000
    return C_mic_local,F_res,C_mic_rad_local,F_res_rad,D14C_mic_local,D14C_res_local   
    #return Flb, Fbr, Fbm