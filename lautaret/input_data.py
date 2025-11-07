# Col du Lautaret 
# 1 fichier par année hydro de 2012 à 2024
# Données toutes les 30min, du 17/10-00:00:00 au 27/11-15:00:00

import csv
import numpy as np
import pandas as pd
import glob
import os

# #*****************************************************************************************
# # Concaténation des fichiers annuels du Lautaret
# #*****************************************************************************************
# years = ['2017-2018', '2018-2019', '2019-2020', '2020-2021', '2021-2022', '2022-2023', '2023-2024']

# csv_files = [f for yr in years for f in glob.glob(os.path.join('lautaret/', yr+'_Lautaret_halfhour.csv'))]
# data = pd.concat((pd.read_csv(f) for f in csv_files), ignore_index=True)
# data['DateTime']=data['TIMESTAMP'].combine_first(data['Unnamed: 0'])
# data.drop(columns=['TIMESTAMP', 'Unnamed: 0'], inplace=True)

# # for an in ['2012-2013', '2013-2014', '2014-2015', '2015-2016', '2016-2017', '2017-2018', '2018-2019', 
# #            '2019-2020', '2020-2021', '2021-2022', '2022-2023', '2023-2024']
# #     with open('lautaret/'+an+'_Lautaret_halfhour.csv', newline='') as file:
# #         data = pd.read_csv(file)

# (row,col)=data.shape
# print(row,col)

# #*****************************************************************************************
# # Calcul des précipitations et gap filling par fonction d'interpolation
# # Changement d'unités
# #*****************************************************************************************
# # data.rename(columns={"Unnamed: 0": "DateTime"}, inplace=True)
# data['Sf']=np.zeros(row)
# data['Rf']=np.zeros(row)
# data['Sf24h']=np.zeros(row)
# data['Null']=np.zeros(row)

# data.interpolate(method='linear', inplace=True)  # interpolation linéaire des données manquantes
# data.fillna(method='bfill', inplace=True)        # gap fill propagation backward pour les données manquantes restantes

# is_rain=1.0
# init=0
# for r in range(row):
#     data.loc[r,'AirTC_Avg']+=273.15              # K
#     data.loc[r,'Patm_Avg']*=100                  # Pa
#     data.loc[r,'Quantity_klok']*=0.001           # m
#     if data['is_rain_klok'][r]==is_rain:
#         data.loc[r,'Rf']=data.loc[r,'Quantity_klok']
#     else:
#         data.loc[r,'Sf']=data.loc[r,'Quantity_klok']

# init=sum(data['Sf'].iloc[0:47])
# for r in range(48,row):
#     data.loc[r,'Sf24h']=sum(data['Sf'].iloc[r-47:r])

# print(data.loc[:,['DateTime','Sf','Sf24h']])
# print(init)

# #*****************************************************************************************
# # Sélection des variables d'intéret
# #*****************************************************************************************
# data_input=data[['DateTime', 'AirTC_Avg', 'HRair_Avg', 'Patm_Avg', 'WindSpeed_Low_Avg', 
#                  'WindSpeed_Avg', 'short_up_Avg', 'short_dn_Avg', 'long_up_cor_Avg',
#                  'long_dn_cor_Avg', 'Rs_net_Avg', 'Rl_net_Avg', 'Albedo_Avg', 'Snow_Depth', 'Sf', 'Rf', 'Sf24h', 'Null']]

# data_input.to_csv("lautaret/2017-2024_Lautaret_halfhour_input.csv")

#*****************************************************************************************
# Mise à zéros des valeurs négatives de short_up_Avg
#*****************************************************************************************
with open('lautaret/2017-2024_Lautaret_halfhour_input.csv', newline='') as file:
    full_data = pd.read_csv(file)

(row,col)=full_data.shape

full_data = full_data.iloc[:,1:]

for r in range(row):
    if full_data.loc[r, 'short_up_Avg']<0:
        full_data.loc[r, 'short_up_Avg']=0
    if pd.to_datetime(full_data.loc[r, 'DateTime']).month>4 and pd.to_datetime(full_data.loc[r, 'DateTime']).month<10:
        full_data.loc[r, 'Snow_Depth']=0
    
print(full_data)

full_data.to_csv("lautaret/2017-2024_Lautaret_halfhour_input_corr.csv")





