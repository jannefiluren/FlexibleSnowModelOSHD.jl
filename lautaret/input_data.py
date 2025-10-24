# Col du Lautaret 
# 1 fichier par année hydro de 2012 à 2024
# Données toutes les 30min, du 17/10-00:00:00 au 27/11-15:00:00

import csv
import numpy as np
import pandas as pd

with open('lautaret/2223_Lautaret_halfhour.csv', newline='') as file:
    data = pd.read_csv(file)

(row,col)=data.shape
print(row,col)

data.rename(columns={"Unnamed: 0": "DateTime"}, inplace=True)
data['Sf']=np.zeros(row)
data['Rf']=np.zeros(row)
data['Sf24h']=np.zeros(row)
data['Null']=np.zeros(row)

data.interpolate(method='linear', inplace=True)  # interpolation linéaire des données manquantes
data.fillna(method='bfill', inplace=True)        # gap fill propagation backward pour les données manquantes restantes

is_rain=1
init=0
for r in range(row):
    data.loc[r,'AirTC_Avg']+=273.15              # K
    data.loc[r,'Patm_Avg']*=100                  # Pa
    data.loc[r,'Quantity_klok']*=0.001           # m
    if data['is_rain_klok'][r]==is_rain:
        data.loc[r,'Rf']=data.loc[r,'Quantity_klok']
    else:
        data.loc[r,'Sf']=data.loc[r,'Quantity_klok']

init=sum(data['Sf'].iloc[0:47])
for r in range(48,row):
    data.loc[r,'Sf24h']=sum(data['Sf'].iloc[r-47:r])

print(data.loc[:,['DateTime','Sf','Sf24h']])
print(init)

data_input=data[['DateTime', 'AirTC_Avg', 'HRair_Avg', 'Patm_Avg', 'WindSpeed_Low_Avg', 
                 'WindSpeed_Avg', 'short_up_Avg', 'short_dn_Avg', 'long_up_cor_Avg',
                 'long_dn_cor_Avg', 'Rs_net_Avg', 'Rl_net_Avg', 'Sf', 'Rf', 'Sf24h', 'Null']]

data_input.to_csv("lautaret/2223_Lautaret_halfhour_input.csv")
