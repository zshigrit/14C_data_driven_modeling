#model_MAOM
# from model_parameters import *
import model_scalars as m_scal 
import math

def millennial(C_lmwc,C_maom,C_mic,C_agg,day,c1=0.4833,clay_per=40,c2=2.3282,BulkD=1150,Klm=0.25,Kmm=0.025,\
    Kagg=0.0002,pa=1./3.,Vma=0.07,Kma=200,Amax=500,Vml=0.01,M_Lmin=10,Kml=25):
    St, Sw, CUE = m_scal.millennial(day)
    Qmax1 = 10**(c1 * math.log(clay_per,10) + c2)
    Qmax = BulkD * Qmax1 * 10**(-3) # convert to gC/m2 (flagged)
#input
    Flm = C_lmwc * (Klm * Qmax * C_lmwc/(1+(Klm * C_lmwc)) - C_maom ) / Qmax * St *Sw # adsorption of LMWC
    Fbm = Kmm * C_mic * St * Sw # adsorption by mineral surface
    Fagg = Kagg * C_agg * St * Sw # break of aggregated carbon
    Fam = Fagg * (1-pa) # allocation of Fagg (flag: pa = 0.5 in Xu and 1/3 in Abramoff)
# output
    Fma = Vma * C_maom /(Kma + C_maom) * (1 - C_agg/Amax) * St * Sw # flag: Kma is 2000 in Xu (200 in Abramoff)
    Fml = Vml * (C_maom - M_Lmin) / (Kml + C_maom - M_Lmin) * St * Sw

    C_maom_local = C_maom + Flm + Fbm + Fam - Fma - Fml
    return C_maom_local 

# in Xiaofeng Xu github (https://github.com/email-clm/Millennial/blob/master/main.F90)
# f_LM_MI_sor = (temp / Qmax + 0.0015) * LMWC / 50. * t_scalar * w_scalar
# temp = (klmc_min * Qmax * LMWC ) / (2. + klmc_min * LMWC) - MINERAL
# f_SO_MI_break = f_SO_break * 1.5 / 3. (Fam)
# f_SO_break = SOILAGG * kagg * t_scalar * w_scalar
# f_MI_SO_agg = Vmin_agg * MINERAL / (kmin_agg + MINERAL) * (1. - SOILAGG / AGGmax)
# f_MI_LM_des = Vm_l * (MINERAL - M_Lmin) / (km_l + MINERAL - M_Lmin) * t_scalar * w_scalar

#model_MAOM_Xu
# from model_parameters_xu import *
# import model_scalars as m_scal 

def millennial_xu(C_lmwc,C_maom,C_mic,C_agg,day,clay_per=40,Klmc_min=0.25,Kmic=0.036,\
    Kagg=0.0002,pa=1.5/3.,Vma=0.07,Kma=2000,Amax=500,Vml=0.01,M_Lmin=10,Kml=25):
    St, Sw, CUE = m_scal.millennial_xu(day)
    Qmax = 10.0 ** (0.297 * math.log(clay_per,10) + 2.355 + 0.50) #!* 1.25
#input
    # Flm = C_lmwc * (Klm * Qmax * C_lmwc/(1+(Klm * C_lmwc)) - C_maom ) / Qmax * St *Sw # adsorption of LMWC
    temp_Flm = (Klmc_min * Qmax * C_lmwc ) / (2. + Klmc_min * C_lmwc) - C_maom
    Flm = (temp_Flm / Qmax + 0.0015) * C_lmwc / 50. * St * Sw
    # Fbm = Kmm * C_mic * St * Sw # adsorption by mineral surface
    Fbm = C_mic * Kmic * 0.15 * St * Sw  
    # Fagg = Kagg * C_agg * St * Sw # break of aggregated carbon
    Fagg = C_agg * Kagg * St * Sw
    Fam = Fagg * (1-pa) # allocation of Fagg
# output
    # Fma = Vma * C_maom /(Kma + C_maom) * (1 - C_agg/Amax) * St * Sw
    Fma = Vma * C_maom / (Kma + C_maom) * (1. - C_agg / Amax)
    Fml = Vml * (C_maom - M_Lmin) / (Kml + C_maom - M_Lmin) * St * Sw

    C_maom_local = C_maom + Flm + Fbm + Fam - Fma - Fml
    return C_maom_local 

#
# radiocarbon processes
def millennial_rad(C_lmwc,C_maom,C_mic,C_agg,C_lmwc_rad,C_maom_rad,C_mic_rad,C_agg_rad,day,c1=0.4833,\
    clay_per=40,c2=2.3282,BulkD=1150,Klm=0.25,Kmm=0.025,lambdax = 8267, Aabs=10**-12,\
    Kagg=0.0002,pa=1./3.,Vma=0.07,Kma=200,Amax=500,Vml=0.01,M_Lmin=10,Kml=25):
    St, Sw, CUE = m_scal.millennial(day)
    Qmax1 = 10**(c1 * math.log(clay_per,10) + c2)
    Qmax = BulkD * Qmax1 * 10**(-3) # convert to gC/m2 (flagged)
    #       
    Asn_lmwc = C_lmwc_rad/C_lmwc
    Asn_maom = C_maom_rad/C_maom 
    Asn_mic = C_mic_rad/C_mic
    Asn_agg = C_agg_rad/C_agg
#input
    Flm = C_lmwc * (Klm * Qmax * C_lmwc/(1+(Klm * C_lmwc)) - C_maom ) / Qmax * St *Sw # adsorption of LMWC
    Fbm = Kmm * C_mic * St * Sw # adsorption by mineral surface
    Fagg = Kagg * C_agg * St * Sw # break of aggregated carbon
    Fam = Fagg * (1-pa) # allocation of Fagg (flag: pa = 0.5 in Xu and 1/3 in Abramoff) 

    Flm_rad = C_lmwc * (Klm * Qmax * C_lmwc/(1+(Klm * C_lmwc)) - C_maom ) / Qmax * St * Sw * Asn_lmwc  # adsorption of LMWC
    Fbm_rad = Kmm * C_mic * St * Sw * Asn_mic # adsorption by mineral surface
    Fagg_rad = Kagg * C_agg * St * Sw * Asn_agg # break of aggregated carbon
    Fam_rad = Fagg_rad * (1-pa) # allocation of Fagg (flag: pa = 0.5 in Xu and 1/3 in Abramoff)
    # Flm_rad = C_lmwc_rad * (Klm * Qmax * C_lmwc_rad/(1+(Klm * C_lmwc_rad)) - C_maom_rad) / Qmax * St *Sw # adsorption of LMWC
    # Fbm_rad = Kmm * C_mic_rad * St * Sw # adsorption by mineral surface
    # Fagg_rad = Kagg * C_agg_rad * St * Sw # break of aggregated carbon
    # Fam_rad = Fagg_rad * (1-pa) # allocation of Fagg (flag: pa = 0.5 in Xu and 1/3 in Abramoff)

# output
    Fma = Vma * C_maom /(Kma + C_maom) * (1 - C_agg/Amax) * St * Sw # flag: Kma is 2000 in Xu (200 in Abramoff)
    Fml = Vml * (C_maom - M_Lmin) / (Kml + C_maom - M_Lmin) * St * Sw 

    Fma_rad = Vma * C_maom /(Kma + C_maom) * (1 - C_agg/Amax) * St * Sw * Asn_maom # flag: Kma is 2000 in Xu (200 in Abramoff)
    Fml_rad = Vml * (C_maom - M_Lmin) / (Kml + C_maom - M_Lmin) * St * Sw * Asn_maom

    # Fma_rad = Vma * C_maom_rad /(Kma + C_maom_rad) * (1 - C_agg_rad/Amax) * St * Sw # flag: Kma is 2000 in Xu (200 in Abramoff)
    # Fml_rad = Vml * (C_maom_rad - M_Lmin) / (Kml + C_maom_rad - M_Lmin) * St * Sw

    C_maom_local = C_maom + Flm + Fbm + Fam - Fma - Fml

    C_maom_rad_local = C_maom_rad + Flm_rad + Fbm_rad + Fam_rad - Fma_rad - Fml_rad - (1/(lambdax*365))*C_maom_rad

    # change unit to D14C per mil
    D14C_maom_local = ((C_maom_rad_local/C_maom_local)/Aabs-1)*1000
    return C_maom_local,C_maom_rad_local,D14C_maom_local  