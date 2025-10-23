# Col du Lautaret 
# 1 fichier par année hydro de 2012 à 2024
# Données toutes les 30min, du 17/10-00:00:00 au 27/11-15:00:00

import csv
import numpy as np

with open('2223_Lautaret_halfhour.csv', newline='') as datacsv:
    file = csv.reader(datacsv)
    data = np.array(list(file))

print(data[0])

i=0
while data[0][i] != 'Snow_Depth':
    i+=1
print(i)
print(data[1:,i])

j=0
while data[0][j] != 'Quantity_klok':
    j+=1
print(j)
print(data[1:,j])

row = np.shape(data)[0]
col = np.shape(data)[1]

data_input=np.append(data, np.zeros((row,3)), axis=1)
data_input[0, 49:]=['Sf', 'Rf', 'Null']
data_input[0,0]='Date'

snow_depth=data_input[1,i]
for r in range(2,row):
    if data_input[r,i]!=snow_depth:
        snow_depth=data_input[r,i]
        data_input[r, 49:-1]=[data_input[r,j],'0.0']
    else:
        data_input[r, 49:-1]=['0.0', data_input[r,j]]

print(data_input)
print(np.shape(data_input))

with open("2223_Lautaret_halfhour_input.csv", mode="w", newline="") as csvfile:
    file = csv.writer(csvfile)
    file.writerows(data_input)
