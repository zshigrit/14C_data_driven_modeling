#%% import libraries
#import os 
#os.chdir('/Users/shanghuasun/Documents/my_projects/14C_modeling/model_lib')
#cd model_lib
import numpy as np
import math
#
#from model_parameters import *
from model_initial_conditions import *
#from model_input import *

import model_microbe as m_mic
import model_MAOM as m_maom
import model_AggC as m_agg 
import model_LMWC as m_lmwc 
import model_POM as m_pom 
#%% pre-assign carbon pools
C_pom_record = []
C_agg_record = []
C_lmwc_record = []
C_maom_record = []
C_mic_record = []
Fin_record = []
F_res_record = []

D14C_pom_record = []
D14C_lmwc_record = []
D14C_maom_record = []
D14C_mic_record = []
D14C_agg_record = []

D14C_res_record = []
F_res_rad_record = []
rad_atm_record = []
#%% model processes 
years = 100
for day in np.arange(years*365):
    #day = days.copy()
    # state variables 
# Abramoff
    # C_pom,Fin = m_pom.millennial(C_pom,C_agg,C_mic,day)
    # C_lmwc = m_lmwc.millennial(C_pom,C_mic,C_maom,C_lmwc,day)
    # C_agg = m_agg.millennial(C_maom,C_pom,C_agg,day)
    # C_maom = m_maom.millennial(C_lmwc,C_maom,C_mic,C_agg,day)
    # C_mic,F_res = m_mic.millennial(C_lmwc,C_mic,day)
# Xiaofeng Xu
    # C_pom,Fin = m_pom.millennial_xu(C_pom,C_agg,C_mic,day)
    # C_lmwc = m_lmwc.millennial_xu(C_pom,C_mic,C_maom,C_lmwc,day)
    # C_agg = m_agg.millennial_xu(C_maom,C_pom,C_agg,day)
    # C_maom = m_maom.millennial_xu(C_lmwc,C_maom,C_mic,C_agg,day)
    # C_mic,F_res = m_mic.millennial_xu(C_lmwc,C_mic,day)

# # radiocarbon
    C_pom,Fin,C_pom_rad,D14C_pom = m_pom.millennial_rad(C_pom,C_agg,C_mic,C_pom_rad,C_agg_rad,C_mic_rad,day,years)
    C_lmwc,C_lmwc_rad,D14C_lmwc,rad_atm = m_lmwc.millennial_rad(C_pom,C_mic,C_maom,C_lmwc,C_pom_rad,C_mic_rad,C_maom_rad,C_lmwc_rad,years,day)
    C_agg,C_agg_rad,D14C_agg = m_agg.millennial_rad(C_maom,C_pom,C_agg,C_maom_rad,C_pom_rad,C_agg_rad,day)
    C_maom,C_maom_rad,D14C_maom = m_maom.millennial_rad(C_lmwc,C_maom,C_mic,C_agg,C_lmwc_rad,C_maom_rad,C_mic_rad,C_agg_rad,day)
    C_mic,F_res,C_mic_rad,F_res_rad,D14C_mic,D14C_res  = m_mic.millennial_rad(C_lmwc,C_mic,C_lmwc_rad,C_mic_rad,day)
    #

# all records
    C_pom_record = np.append(C_pom_record,C_pom) 
    C_agg_record = np.append(C_agg_record,C_agg)
    C_lmwc_record = np.append(C_lmwc_record,C_lmwc)
    C_maom_record = np.append(C_maom_record,C_maom)  
    C_mic_record = np.append(C_mic_record,C_mic)

    Fin_record = np.append(Fin_record,Fin)
    F_res_record = np.append(F_res_record,F_res)

# radiocarbon
    D14C_pom_record = np.append(D14C_pom_record,D14C_pom)
    D14C_lmwc_record = np.append(D14C_lmwc_record,D14C_lmwc)
    D14C_agg_record = np.append(D14C_agg_record,D14C_agg)
    D14C_maom_record = np.append(D14C_maom_record,D14C_maom)
    D14C_mic_record = np.append(D14C_mic_record,D14C_mic)

    D14C_res_record = np.append(D14C_res_record,D14C_res)
    F_res_rad_record = np.append(F_res_rad_record,F_res_rad)
    rad_atm_record = np.append(rad_atm_record,rad_atm)
    




C_total = C_pom_record + C_agg_record + C_lmwc_record\
    + C_maom_record + C_mic_record

#C_pom_record_annual = C_pom_record[0:365]

tot_stock_res = C_total[365*years-1]+sum(F_res_record)
tot_input = sum (Fin_record)
#%% plotting figures note: horizontally and in year unit
import numpy as np
import matplotlib.pyplot as plt
import matplotlib 
#from model_main import years,C_pom_record,C_agg_record,\
#    C_lmwc_record,C_maom_record,C_mic_record,C_total
matplotlib.rcParams.update({'font.size': 18})

time_step = 365
C_pom_record_tem = np.empty(int(years*365/time_step))
C_mic_record_tem = np.empty(int(years*365/time_step))
C_maom_record_tem = np.empty(int(years*365/time_step))
C_agg_record_tem = np.empty(int(years*365/time_step))
C_lmwc_record_tem = np.empty(int(years*365/time_step))
C_total_tem = np.empty(int(years*365/time_step))
Fin_record_tem = np.empty(int(years*365/time_step))
F_res_record_tem = np.empty(int(years*365/time_step)) 

for i in np.arange(int(years*365/time_step)):
    C_pom_record_tem[i] = np.mean(C_pom_record[i*time_step:(i+1)*time_step])
    C_mic_record_tem[i] = np.mean(C_mic_record[i*time_step:(i+1)*time_step])
    C_maom_record_tem[i] = np.mean(C_maom_record[i*time_step:(i+1)*time_step])
    C_agg_record_tem[i] = np.mean(C_agg_record[i*time_step:(i+1)*time_step])
    C_lmwc_record_tem[i] = np.mean(C_lmwc_record[i*time_step:(i+1)*time_step])
    C_total_tem[i] = np.mean(C_total[i*time_step:(i+1)*time_step])
    Fin_record_tem[i] = np.sum(Fin_record[i*time_step:(i+1)*time_step])
    F_res_record_tem[i] = np.sum(F_res_record[i*time_step:(i+1)*time_step])
# %%    
fig, axs = plt.subplots(2,3,figsize = (20,12))

x = np.arange(0,years,time_step/365)

axs[0,0].plot(x,C_pom_record_tem)
axs[0,0].set_ylabel('C_pom (gC/m2)',fontsize = 18)
axs[0,0].set_xlabel('Time (years)')

axs[0,1].plot(x,C_agg_record_tem)
axs[0,1].set_ylabel('C_agg (gC/m2)',fontsize = 18)
axs[0,1].set_xlabel('Time (years)')

axs[0,2].plot(x,C_lmwc_record_tem)
axs[0,2].set_ylabel('C_lmwc (gC/m2)',fontsize = 18)
axs[0,2].set_xlabel('Time (years)')

axs[1,0].plot(x,C_maom_record_tem)
axs[1,0].set_ylabel('C_maom (gC/m2)',fontsize = 18)
axs[1,0].set_xlabel('Time (years)')

axs[1,1].plot(x,C_mic_record_tem)
axs[1,1].set_ylabel('C_mic (gC/m2)',fontsize = 18)
axs[1,1].set_xlabel('Time (years)')

axs[1,2].plot(x,C_total_tem)
axs[1,2].set_ylabel('C_total (gC/m2)',fontsize = 18)
axs[1,2].set_xlabel('Time (years)')
#plt.show()
#st = plt.suptitle('initial condition AGG'+ '%d' %C_agg0)
#st.set_y(0.75)
#st.set_x(0.55)
plt.tight_layout()
plt.savefig('soil_carbon_horizontal_annual.png')
#%% input. output and soil total carbon
fig, axs = plt.subplots(1,3,figsize = (20,6))
axs[0].plot(x,C_total_tem)
axs[0].set_ylabel('C_total (gC/m2)',fontsize = 18)
axs[0].set_xlabel('Time (years)')

axs[1].plot(x,Fin_record_tem)
axs[1].set_ylabel('Carbon input (gC/m2)',fontsize = 18)
axs[1].set_xlabel('Time (years)')

axs[2].plot(x,F_res_record_tem)
axs[2].set_ylabel('Respiration (gC/m2)',fontsize = 18)
axs[2].set_xlabel('Time (years)')

plt.tight_layout()
plt.savefig('carbon_input_and_output_horizontal.png')
# %% radiocarbon
#
# radiocarbon
matplotlib.rcParams.update({'font.size': 18})
fig, axs = plt.subplots(2,3,figsize = (20,12))
x = np.arange(0,years,1/365)

axs[0,0].plot(x,D14C_pom_record)
axs[0,0].set_ylabel('D14C_pom (per mil)',fontsize = 18)
axs[0,0].set_xlabel('Time (years)')

axs[0,1].plot(x,D14C_agg_record)
axs[0,1].set_ylabel('D14C_agg (per mil)',fontsize = 18)
axs[0,1].set_xlabel('Time (years)')

axs[0,2].plot(x,D14C_lmwc_record)
axs[0,2].set_ylabel('D14C_lmwc (per mil)',fontsize = 18)
axs[0,2].set_xlabel('Time (years)')

axs[1,0].plot(x,D14C_maom_record)
axs[1,0].set_ylabel('D14C_maom (per mil)',fontsize = 18)
axs[1,0].set_xlabel('Time (years)')

axs[1,1].plot(x,D14C_mic_record)
axs[1,1].set_ylabel('D14C_mic (per mil)',fontsize = 18)
axs[1,1].set_xlabel('Time (years)')

axs[1,2].plot(x,rad_atm_record)
axs[1,2].set_ylabel('rad_atm (per mil)',fontsize = 18)
axs[1,2].set_xlabel('Time (years)')

# axs[1,2].plot(x,D14C_res_record)
# axs[1,2].set_ylabel('D14C_res (per mil)',fontsize = 18)
# axs[1,2].set_xlabel('Time (years)')
#plt.show()
#st = plt.suptitle('initial condition AGG'+ '%d' %C_agg0)
#st.set_y(0.75)
#st.set_x(0.55)
plt.tight_layout()
plt.savefig('radiocarbon_horizontal_with_atm_rad.png')

#%% atmospheric radiocarbon
plt.figure(figsize=(8,6))
x = np.arange(0,years,1/365)
plt.plot(x,rad_atm_record)
plt.ylabel('atom_radiocarbon (per mil)',fontsize = 18)
plt.xlabel('Time (years)')
plt.savefig('atm_radiocarbon.png')
# plt.show()

# %%
