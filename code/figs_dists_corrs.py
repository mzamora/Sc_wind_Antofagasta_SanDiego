#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 13:29:33 2024

@author: monica
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
mat = scipy.io.loadmat('IndicesAntofaga.mat')
mat2 = scipy.io.loadmat('IndicesSd.mat')

#Antofagasta
an_vars=['tinicialruptura','tfinalruptura','deltaruptura','turb','kruptura', #4
      'tinicio','rglbprom','rglbmax','rglbmin','stdrglb','promkrup','tempromrup', #11
      'tempminrup','tempmaxrup','stdtemp','hrelpromrup','hrelmaxrup','hrelminrup', #17
      'stdhrel','prespromrup','presmaxrup','presminrup','stdpres','tpotpromrup', #23
      'tpotmaxrup','tpotminrup','stdtpot','tdewpromrup','tdewminrup','tdewmaxrup', #29
      'stdtdew','n_clouds','LCL_srf','z_cloudbase','z_inv_base','z_inv_top','qT_BL', #36
      'qT_jump','qT_3km','thetaL_BL','thetaL_jump','thetaL_3km','BLwnd_dir_avg', #42
      'BLwnd_spd_avg','dq_decoupling','dtheta_decoupling','promedionb','mediananb', #47
      'desvnb','maxnb','minnb','promvelv','promvelu','prommag','maxmagnitudviento', #54
      'minmagnitudviento','inicioviento','maginicioviento','Vamanecer','promagruptura', #59
      'maxmagruptura','minmagruptura0','maginicioruptura','magfirup'] #63

an_vars2=['$t_i$','$t_f$','$\Delta t$','LT','$k_i$', #4
         '$t_0$','$\overline{GHI}_{if}$','$GHI_{max,if}$','$GHI_{min,if}$','$SD_{GHI,if}$','$\overline{k}_{if}$', #10
         '$\overline{T}_{if}$','$T_{min,if}$','$T_{max,if}$','$SD_{T,if}$','$\overline{RH}_{if}$', #15
         '$RH_{max,if}$','$RH_{min,if}$','$SD_{RH,if}$',
         '$\overline{p}_{if}$','$p_{min,if}$','$p_{max,if}$','$SD_{p,if}$',
         '$\overline{\\theta}_{if}$','$\\theta_{min,if}$','$\\theta_{max,if}$','$SD_{\\theta,if}$',
         '$\overline{T_d}_{if}$','${T}_{d,min,if}$','${T}_{d,max,if}$','$SD_{T_d,if}$',
         '$n_{clouds}$','LCL','$z_b$', #34
         '$z_{ib}$','$z_{it}$','$q$','$\Delta q$','$q_{3km}$','$\\theta$', #40
         '$\Delta \\theta$','$\\theta_{3km}$','$\\alpha_{ABL}$','$w_{ABL}$', #44
         '$\Delta_{dec} q$','$\Delta_{dec} \\theta$',
         '$\overline{z}_{b,ceil}$','$Med. z_{b,ceil}$', '$SD_{z_{b,ceil}}$', 
         '$z_{b,ceil,max}$','$z_{b,ceil,min}$',
         '$\overline{v}$','$\overline{u}$','$\overline{w}$', #21
         '$w_{max}$','$w_{min}$','$t_{sb}$','$w_{sb}$', #25
         '$w_0$','$\overline{w}_{if}$','$w_{max,if}$','$w_{min,if}$', #29
         '$w_{i}$','$w_{f}$',] #46

an_data=mat['resultados']

sd_vars=['tinicialruptura','tfinalruptura','deltaruptura','turb','kruptura', #4
         'tinicio','rglbprom','rglbmax','rglbmin','stdrglb','promkrup', #10
         'tempromrup','tempminrup','tempmaxrup','stdtemp','hrelpromrup', #15
         'hrelmaxrup','hrelminrup','stdhrel','promvelv','promvelu','prommag', #21
         'maxmagnitudviento','minmagnitudviento','inicioviento','maginicioviento', #25
         'Vamanecer','promagruptura','maxmagruptura','minmagruptura', #29
         'maginicioruptura','magfirup','n_clouds','LCL_srf','z_cloudbase', #34
         'z_inv_base','z_inv_top','qT_BL','qT_jump','qT_3km','thetaL_BL', #40
         'thetaL_jump','thetaL_3km','BLwnd_dir_avg','BLwnd_spd_avg', #44
         'dq_decoupling','dtheta_decoupling'] #46

sd_vars2=['$t_i$','$t_f$','$\Delta t$','LT','$k_i$', #4
         '$t_0$','$\overline{GHI}_{if}$','$GHI_{max,if}$','$GHI_{min,if}$','$SD_{GHI,if}$','$\overline{k}_{if}$', #10
         '$\overline{T}_{if}$','$T_{min,if}$','$T_{max,if}$','$SD_{T,if}$','$\overline{RH}_{if}$', #15
         '$RH_{max,if}$','$RH_{min,if}$','$SD_{RH,if}$','$\overline{v}$','$\overline{u}$','$\overline{w}$', #21
         '$w_{max}$','$w_{min}$','$t_{sb}$','$w_{sb}$', #25
         '$w_0$','$\overline{w}_{if}$','$w_{max,if}$','$w_{min,if}$', #29
         '$w_{i}$','$w_{f}$','$n_{clouds}$','LCL','$z_b$', #34
         '$z_{ib}$','$z_{it}$','$q$','$\Delta q$','$q_{3km}$','$\\theta$', #40
         '$\Delta \\theta$','$\\theta_{3km}$','$\\alpha_{ABL}$','$w_{ABL}$', #44
         '$\Delta_{dec} q$','$\Delta_{dec} \\theta$'] #46

sd_data=mat2['RSd']

#filtering out events shorter than 5 min
an_data=an_data[:,an_data[2]>5]
sd_data=sd_data[:,sd_data[2]>5]

# to utc
sd_data[24,:]=sd_data[24,:]-8*60

# ceilometer: agl to asl
an_data[46,:]=an_data[46,:]+135
an_data[47,:]=an_data[47,:]+135
an_data[49,:]=an_data[49,:]+135
an_data[50,:]=an_data[50,:]+135

#%% plots
# met conditions
fig,ax=plt.subplots(2,3,figsize=(10,5))
ax[1][0].hist(an_data[33,:],fc='b',alpha=0.8,histtype='step',density=True)
ax[1][0].hist(sd_data[34,:],fc='b',alpha=0.8,histtype='step',density=True)
ax[1][0].set_xlabel('Cloud base height (m)')
ax[0][0].hist(an_data[34,:],fc='b',alpha=0.8,histtype='step',density=True)
ax[0][0].hist(sd_data[35,:],fc='b',alpha=0.8,histtype='step',density=True)
ax[0][0].set_xlabel('Inversion height (m)')
ax[0][1].hist(an_data[36,:],fc='b',alpha=0.8,histtype='step',density=True)
ax[0][1].hist(sd_data[37,:],fc='b',alpha=0.8,histtype='step',density=True)
ax[0][1].set_xlabel('ABL mean $q$ (g/kg)')
ax[0][2].hist(an_data[37,:],fc='b',alpha=0.8,histtype='step',density=True)
ax[0][2].hist(sd_data[38,:],fc='b',alpha=0.8,histtype='step',density=True)
ax[0][2].set_xlabel('Inversion jump $\Delta q$ (g/kg)')
ax[1][1].hist(an_data[39,:],fc='b',alpha=0.8,histtype='step',density=True)
ax[1][1].hist(sd_data[40,:],fc='b',alpha=0.8,histtype='step',density=True)
ax[1][1].set_xlabel('ABL mean $\\theta$ (K)')
ax[1][2].hist(an_data[40,:],fc='b',alpha=0.8,histtype='step',density=True)
ax[1][2].hist(sd_data[41,:],fc='b',alpha=0.8,histtype='step',density=True)
ax[1][2].set_xlabel('Inversion jump $\Delta \\theta$ (K)')
for i1 in range(2): 
    for i2 in range(3): 
        ax[i1][i2].set_ylabel('Frequency')
        #if i2>0: ax[i1][i2].set_yticks([])
ax[1][0].legend(['Antofagasta','San Diego'])
plt.tight_layout()
plt.savefig('Fig3_methists.png',dpi=300)

#%% Cross relations
fig,ax=plt.subplots(2,2,figsize=(8,5))
ax[0][0].plot(an_data[40,:],an_data[34,:],'.',alpha=0.2)
ax[0][0].plot(sd_data[41,:],sd_data[35,:],'.',alpha=0.2)
ax[0][0].set_xlabel('Inversion jump $\Delta \\theta$ (K)')
ax[0][0].set_ylabel('Inversion height (m)')
ax[0][1].plot(an_data[36,:],an_data[34,:],'.',alpha=0.2)
ax[0][1].plot(sd_data[37,:],sd_data[35,:],'.',alpha=0.2)
ax[0][1].set_xlabel('ABL mean $q$ (g/kg)')
ax[0][1].set_ylabel('Inversion height (m)')
ax[1][0].plot(an_data[39,:],an_data[34,:],'.',alpha=0.2)
ax[1][0].plot(sd_data[40,:],sd_data[35,:],'.',alpha=0.2)
ax[1][0].set_xlabel('ABL mean $\\theta$ (K)')
ax[1][0].set_ylabel('Inversion height (m)')
ax[1][1].plot(an_data[36,:],an_data[39,:],'.',alpha=0.2)
ax[1][1].plot(sd_data[37,:],sd_data[40,:],'.',alpha=0.2)
ax[1][1].set_xlabel('ABL mean $q$ (g/kg)')
ax[1][1].set_ylabel('ABL mean $\\theta$ (K)')
# ax[1][1].plot(an_data[36,:],an_data[33,:],'.',alpha=0.2)
# ax[1][1].plot(sd_data[37,:],sd_data[34,:],'.',alpha=0.2)
# ax[1][1].set_xlabel('ABL mean $q_t$ (g/kg)')
# ax[1][1].set_ylabel('Cloud base height (m)')
# ax[1][2].plot(an_data[34,:]-an_data[33,:],an_data[34,:],'.',alpha=0.2)
# ax[1][2].plot(sd_data[35,:]-sd_data[34,:],sd_data[35,:],'.',alpha=0.2)
# ax[1][2].set_xlabel('Cloud thickness (m)')
# ax[1][2].set_ylabel('Inversion height (m)')
# ax[1][2].set_xlim([-10,1000])
ax[0][0].legend(['Antofagasta','San Diego'])
plt.tight_layout()
plt.savefig('Fig4_metcross.png',dpi=300)

#%% Create dataframes and save them
import pandas as pd
sd=pd.DataFrame(sd_data.T,columns=sd_vars)
an=pd.DataFrame(an_data.T,columns=an_vars)

sd.tinicialruptura=sd.tinicialruptura-sd.tinicio
sd.tfinalruptura=sd.tfinalruptura-sd.tinicio
sd.inicioviento=sd.inicioviento-sd.tinicio

an.tinicialruptura=an.tinicialruptura-an.tinicio
an.tfinalruptura=an.tfinalruptura-an.tinicio
an.inicioviento=an.inicioviento-an.tinicio

sd.to_csv('sd_final_data.csv',sep=',')
an.to_csv('an_final_data.csv',sep=',')

#%%
import seaborn as sb

sd2=sd.set_axis(sd_vars2,axis=1)
sd2['h']=sd.z_inv_base-sd.z_cloudbase
sdcorr = sd2.corr()
an2=an.set_axis(an_vars2,axis=1)
an2['h']=an.z_inv_base-an.z_cloudbase
an2['$h_{ceil}$']=an.z_inv_base-an.promedionb
ancorr = an2.corr()

fig,ax=plt.subplots(figsize=(10,8))  
sb.heatmap(sdcorr, cmap="RdBu", annot=False, xticklabels=1, yticklabels=1, ax=ax)
plt.yticks(fontsize=9)
plt.xticks(fontsize=9)
plt.tight_layout()
plt.savefig('Fig10b_corr_sd.png', dpi=300)
plt.close()

fig,ax=plt.subplots(figsize=(10,8))  
sb.heatmap(ancorr, cmap="RdBu", annot=False, xticklabels=1, yticklabels=1, ax=ax)
plt.yticks(fontsize=8)
plt.xticks(fontsize=8)
plt.tight_layout()
plt.savefig('Fig10a_corr_an.png', dpi=300)
plt.close()

sdrelevant=np.logical_or(sdcorr>0.5,sdcorr<-0.5)*sdcorr
anrelevant=np.logical_or(ancorr>0.5,ancorr<-0.5)*ancorr

fig,ax=plt.subplots(1,2,figsize=(17,8))  
sb.set(font_scale=0.6)
sb.heatmap(anrelevant, cmap='RdBu', cbar=False, annot=False, xticklabels=1, yticklabels=1, ax=ax[0])
sb.heatmap(sdrelevant, cmap='RdBu', annot=False, xticklabels=1, yticklabels=1, ax=ax[1])
ax[0].set_title('Antofagasta')
ax[1].set_title('San Diego')
plt.tight_layout()
plt.savefig('Fig_rel_corrs.png', dpi=300)
plt.close()

sb.set(font_scale=1)
fig,ax=plt.subplots(figsize=(10,8))  
sb.heatmap(anrelevant, cmap="RdBu", annot=False, xticklabels=1, yticklabels=1, ax=ax)
plt.yticks(fontsize=9); plt.xticks(fontsize=9)
plt.tight_layout()
plt.savefig('Fig_relcorrs_an.png', dpi=300)
plt.close()

fig,ax=plt.subplots(figsize=(10,8))  
sb.heatmap(sdrelevant, cmap="RdBu", annot=False, xticklabels=1, yticklabels=1, ax=ax)
plt.yticks(fontsize=9); plt.xticks(fontsize=9)
plt.tight_layout()
plt.savefig('Fig_relcorrs_sd.png', dpi=300)
plt.close()

# to keep using matplotlib normally
sb.reset_defaults()

#%%

ip1=2; ip2=4;
fig,ax=plt.subplots(1,1,figsize=(6,4))
ax.plot(an_data[ip1,:],an_data[ip2,:],'.',alpha=0.25,markersize=15)
ax.set_xlabel(an_vars[ip1])
ax.set_ylabel(an_vars[ip2])

#%%

ip1=1; ip2=4;
fig,ax=plt.subplots(1,1,figsize=(6,4))
ax.plot(an_data[ip1,:],an_data[ip2,:],'.',alpha=0.25,markersize=15)
ax.set_xlabel(an_vars[ip1])
ax.set_ylabel(an_vars[ip2])

an_data[0,:].mean

#%%

fig,ax=plt.subplots(1,2,figsize=(9,3))
ax[0].hist((an_data[56,:]-an_data[5,:])/60,histtype=u'step')
ax[0].hist((an_data[0,:]-an_data[5,:])/60,histtype=u'step')
ax[0].hist((an_data[1,:]-an_data[5,:])/60,histtype=u'step')
ax[0].legend(['$t_{sb}$','$t_i$','$t_f$'])
ax[0].set_xlabel('Time after sunrise (h)')
ax[0].set_ylabel('Number of days')
ax[0].set_title('Antofagasta')
ax[0].plot()

ax[1].hist((sd_data[24,:]-sd_data[5,:])/60,histtype=u'step')
ax[1].hist((sd_data[0,:]-sd_data[5,:])/60,histtype=u'step')
ax[1].hist((sd_data[1,:]-sd_data[5,:])/60,histtype=u'step')
ax[1].legend(['$t_{sb}$','$t_i$','$t_f$'])
ax[1].set_ylabel('Number of days')
ax[1].set_xlabel('Time after sunrise (h)')
ax[1].set_title('San Diego')
plt.tight_layout()

fig.savefig('Fig8_disstimes.png',dpi=300)


#%%

fig,ax=plt.subplots(1,2,figsize=(9,3))
fw=an.promagruptura>an.promagruptura.mean()
ax[0].plot(an.tinicialruptura[fw]/60,an.tfinalruptura[fw]/60,
           '*r',alpha=0.6)
ax[0].plot(an.tinicialruptura[~fw]/60,an.tfinalruptura[~fw]/60,
           '.b',alpha=0.6)
fw=sd.promagruptura>sd.promagruptura.mean()
ax[1].plot(sd.tinicialruptura[fw]/60,sd.tfinalruptura[fw]/60,
           '*r',alpha=0.6)
ax[1].plot(sd.tinicialruptura[~fw]/60,sd.tfinalruptura[~fw]/60,
           '.b',alpha=0.6)
ax[0].set_xlabel('$t_i$ (h after sunrise)')
ax[0].set_ylabel('$t_f$ (h after sunrise)')
ax[1].set_xlabel('$t_i$ (h after sunrise)')
ax[1].set_ylabel('$t_f$ (h after sunrise)')
ax[0].set_title('Antofagasta')
ax[1].set_title('San Diego')
ax[1].legend(['$\overline{w}_{if}$ above mean','$\overline{w}_{if}$ below mean'])
plt.tight_layout()
fig.savefig('Fig_fragtimes_w.png',dpi=300)

#%% dth dq decoupling
fig,ax=plt.subplots(1,2,figsize=(9,3))
ax[0].plot(an.z_inv_base,an.dtheta_decoupling,'.',alpha=0.6)
ax[1].plot(sd.z_inv_base,sd.dtheta_decoupling,'.',alpha=0.6)


#%% ceil vs radiosonde decoupling
fig,ax=plt.subplots(1,2,figsize=(9,3))
ax[0].plot(an.z_cloudbase,an.promedionb,'.',alpha=0.6)
ax[0].plot([0,2000],[0,2000],'-r')
ax[1].plot(an.z_inv_base-an.z_cloudbase,an.z_inv_base-an.promedionb,'.',alpha=0.6)
ax[1].plot([0,2000],[0,2000],'-r')
