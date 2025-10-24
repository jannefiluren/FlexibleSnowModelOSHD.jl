import numpy as np
# import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import csv
import pandas as pd

#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# Ouverture des fichiers output

data = pd.read_csv("lautaret/output_lautaret_halfhour.txt", delimiter=',')
data['time']=pd.to_datetime(data['time'])
print(data)

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

fig, ax = plt.subplots()
ax.plot(data['time'], data['Ts'])

ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))

# plt.locator_params(axis='x', nbins=10)
plt.xticks(rotation=45)  
plt.tight_layout()
plt.show()