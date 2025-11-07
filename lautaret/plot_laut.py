import numpy as np
# import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import csv
import pandas as pd

#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# Ouverture des fichiers output

data = pd.read_csv("lautaret/output_2017-2024_lautaret_halfhour_not-corr.txt", delimiter=',')
data_slf = pd.read_csv('data\output_SLF_5WJ_corr-dt.txt', delimiter=',')
data['time']=pd.to_datetime(data['time'])
data_slf['time']=pd.to_datetime(data_slf['time'])

obs = pd.read_csv("lautaret/2017-2024_Lautaret_halfhour_input_corr.csv", delimiter=',')
obs.rename(columns={"DateTime": "time"}, inplace=True)                             
obs['time']=pd.to_datetime(obs['time'])
start = "2021-09-01 06:00:00"
end = "2022-08-01 05:00:00"
obs_laut = obs[(obs['time']>=start) & (obs['time']<=end)]
obs_slf = pd.read_csv("data\input_SLF_5WJ.txt", delim_whitespace=True)
obs_slf['time']=pd.to_datetime(obs_slf[['year','month','day','hour']])
print(obs_slf)

#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# Définition des variables à tracer:
# hs - hauteur de neige (m)
# Tsnow - temperature des couches de neige
# Ts - temperature de surface
# albedo
# I - masse de glace dans la couche de neige
# W - masse d'eau liquide dans la couche de neige
# snow_depth - 14 derniers jours
# swe - 14 derniers jours

# #''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# # Plot Lautaret

# fig0, (ax0, ax1, ax2, ax3) = plt.subplots(nrows=4, ncols=1, layout='constrained', figsize=(8, 6))
# ax0.plot(obs['time'], obs['Snow_Depth'], color='red')
# ax0.plot(data['time'], data['hs'])
# ax1.plot(data['time'], data['I'])
# ax2.plot(data['time'], data['W'])
# ax3.scatter(obs['time'], obs['Albedo_Avg'], marker='.', label='obs', color='red')
# ax3.scatter(data['time'], data['albedo'], marker='.', label='simu')
# fig0.suptitle('Data Col du Lautaret - FSM/Obs', fontsize=14, va='top')

# for ax in [ax0, ax1, ax2, ax3]:
#     ax.xaxis.set_major_locator(mdates.MonthLocator(interval=6))
#     ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
#     ax.tick_params("x", rotation=30)

# ax3.set_ylim(0.5, 1)
# ax3.legend()

# ax0.set_ylabel('Snow depth (m)')
# ax1.set_ylabel('Ice content (kg/m^2)')
# ax2.set_ylabel('Water content (kg/m^2)')
# ax3.set_ylabel('Albedo')

#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# Plot SLF

fig1, (ax10, ax11, ax12, ax13) = plt.subplots(nrows=4, ncols=1, layout='constrained', figsize=(8, 6))
ax10.plot(data_slf['time'], data_slf['hs'])
# ax11.plot(data_slf['time'], data_slf['I'])
# ax12.plot(data_slf['time'], data_slf['W'])
# ax13.scatter(data_slf['time'], data_slf['albedo'], marker='.', label='simu')
fig1.suptitle('Data SLF - FSM', fontsize=14, va='top')

# for ax in [ax10, ax11, ax12, ax13]:
#     ax.xaxis.set_major_locator(mdates.MonthLocator(interval=6))
#     ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
#     ax.tick_params("x", rotation=30)

# ax13.set_ylim(0.5, 1)
# ax13.legend()

ax10.set_ylabel('Snow depth (m)')
# ax11.set_ylabel('Ice content (kg/m^2)')
# ax12.set_ylabel('Water content (kg/m^2)')
# ax13.set_ylabel('Albedo')

# #''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# # Comparaison obs Lautaret/SLF

# fig2, (ax22, ax23, ax24, ax25, ax26, ax27) = plt.subplots(nrows=6, ncols=1, layout='constrained', figsize=(8, 6)) # ax20, ax21, 
# # ax20.plot(obs_slf['time'], obs_slf['Sf'], label='slf')
# # ax20.plot(obs_laut['time'], obs_laut['Sf']*109, label='laut')
# # ax21.plot(obs_slf['time'], obs_slf['Rf'], label='slf')
# # ax21.plot(obs_laut['time'], obs_laut['Rf']*109, label='laut')
# ax22.plot(obs_laut['time'], obs_laut['AirTC_Avg'], label='laut')
# ax22.plot(obs_slf['time'], obs_slf['Ta'], label='slf')
# ax23.plot(obs_laut['time'], obs_laut['WindSpeed_Avg'], label='laut')
# ax23.plot(obs_slf['time'], obs_slf['Ua'], label='slf')
# ax24.plot(obs_laut['time'], obs_laut['Patm_Avg'], label='laut')
# ax24.plot(obs_slf['time'], obs_slf['Ps'], label='slf')
# ax25.plot(obs_laut['time'], obs_laut['HRair_Avg'], label='laut')
# ax25.plot(obs_slf['time'], obs_slf['RH'], label='slf')
# ax26.plot(obs_laut['time'], obs_laut['short_up_Avg'], label='laut')
# ax26.plot(obs_slf['time'], obs_slf['Sdif'], label='slf')
# ax27.plot(obs_laut['time'], obs_laut['long_up_cor_Avg'], label='laut')
# ax27.plot(obs_slf['time'], obs_slf['LW'], label='slf')

# fig2.suptitle('Comparaison SLF/Lautaret - Obs', fontsize=14, va='top')

# # for ax in [ax27]: # ax20, ax21, ax22, ax23
# #     ax.xaxis.set_major_locator(mdates.MonthLocator(interval=6))
# #     ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
# #     ax.tick_params("x", rotation=30)

# # ax23.set_ylim(0.5, 1)
# # ax23.legend()

# # ax20.set_ylabel('Sf')
# # ax21.set_ylabel('Rf')
# ax22.set_ylabel('Ta')
# ax23.set_ylabel('Ua')
# ax24.set_ylabel('Ps')
# ax25.set_ylabel('RH')
# ax26.set_ylabel('SW')
# ax27.set_ylabel('LW')

# fig3, (ax30, ax31) = plt.subplots(nrows=2, ncols=1, layout='constrained', figsize=(8, 6))
# ax30.plot(obs_slf['time'], obs_slf['Sf'], label='slf')
# ax30.plot(obs_laut['time'], obs_laut['Sf']*109, label='laut')
# ax31.plot(obs_slf['time'], obs_slf['Rf'], label='slf')
# ax31.plot(obs_laut['time'], obs_laut['Rf']*109, label='laut')

# fig3.suptitle('Comparaison SLF/Lautaret - Obs', fontsize=14, va='top')

# ax30.set_ylabel('Sf')
# ax31.set_ylabel('Rf')

#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

plt.tight_layout()
plt.show()