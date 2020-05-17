#%% carbon input into the system
#Fi = 1720. # gC/m2/yr
#Fi_day = Fi/365.# gC/m2/day suject to change (flagged)
#%% read in from Xu website
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
Fii = np.loadtxt('/Users/shanghuasun/Documents/my_projects/14C_modeling/Millennial-master/simulation/globalaverage.txt')
Fi = Fii[:,2]# input 3 is the carbon input in the Abramoff paper
forc_st = Fii[:,0]
forc_sw = Fii [:,1]

# load atmospheric radiocarbon 10000BP to 2010
rad_atm_meta = np.loadtxt('/Users/shanghuasun/Documents/my_projects/14C_modeling/model_lib/mergedatm14C.txt',skiprows=1,delimiter=',')
rad_atm = rad_atm_meta[:,1]
#%matplotlib inline
#fig,axes = plt.subplots(1,2,figsize = (10,12))
# plt.figure()
# plt.subplot(121)
# plt.scatter(np.arange(365),forc_sw)
# plt.subplot(122)
# plt.scatter(np.arange(365),forc_st)
# plt.show(block = True) 



# %%
