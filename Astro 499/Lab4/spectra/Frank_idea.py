import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import glob
import re
from astropy.io import ascii
import csv
import pandas as pd
from scipy import stats
from astropy.table import Table
import pickle



    #To find the linear fit for the different spectrometers

    #Vis_Calibration
PATH4='\\Users\\Jonathan\\Documents\\ASTR499\\Lab4\\graphs\\'
PATH1='\\Users\\Jonathan\\Documents\\ASTR499\\Lab4\\Vis_calibration\\'
os.chdir(PATH1)
filename='mercury4.fit'
hdulist=fits.open(filename)
scidata = hdulist[0].data
tot=np.zeros(scidata.shape[1])

for row in scidata:
    tot+=row
    
y=np.arange(1530)
y=y[::-1]

    #This block makes a list of maxima
L=np.argsort(-tot)
i=0
bad_idx=[]
    #To remove multiple maxima for the same peak. 10 is the critical number, since 5 didn't work
for i in range(len(L)-1):
     if abs(L[i]-L[i+1])<10:
         bad_idx.append(i+1)
L=np.delete(L,bad_idx)

    #Get 3 highest peaks
peaks=y[L[:3]]
    

theo_Vis=[5460, 5769,5790]
theo_Vis=theo_Vis[::-1]


plt.plot(y,tot)
#plt.show()

    #Find the linear fit
fit_Vis=np.polyfit(peaks,theo_Vis,1)
    

y_cal_vis=y*fit_Vis[0]+fit_Vis[1]
    

    #IR-Calibration

PATH2='\\Users\\Jonathan\\Documents\\ASTR499\\Lab4\\IR_calibration\\'
os.chdir(PATH2)
filename='neon6.fit'
hdulist=fits.open(filename)
scidata = hdulist[0].data
tot=np.zeros(scidata.shape[1])



for row in scidata:
    tot+=row
y=np.arange(1530)
y=y[::-1]

L=np.argsort(-tot)

i=0
for i in range(len(L)-1):
    if abs(L[i]-L[i+1])<10:
        bad_idx.append(i+1)
L=np.delete(L,bad_idx)
peaks=y[L[:3]]
    

theo_IR=[5852, 6402,6678]
theo_IR=theo_IR[::-1]

peaks=np.array(peaks)

fit_IR=np.polyfit(peaks,theo_IR,1)
    

plt.plot(y,tot)
#plt.show()

y_cal_ir=y*fit_IR[0]+fit_IR[1]

    


    #HST


PATH1='\\Users\\Jonathan\\Documents\\ASTR499\\Lab4\\Vis_calibration\\'
os.chdir(PATH1)
filename='sirius_stis_002.fits'


table = Table.read(filename)
df = table.to_pandas()
wavelength=[]
wavelength_ir=[]
wavelength_vis=[]
flux_ir=[]
flux_vis=[]

flux=df['FLUX']
wavelength=df['WAVELENGTH']

i=0
for i in range(len(flux)):
    if min(y_cal_vis) < wavelength[i] < max(y_cal_vis):
        wavelength_vis.append(wavelength[i])
        flux_vis.append(flux[i])
    elif min(y_cal_ir) < wavelength[i] < max(y_cal_ir):
        wavelength_ir.append(wavelength[i])
        flux_ir.append(flux[i])

    #For some reason the slope was in the wrong direction. It would be
    #discontinuous if you didn't flip the slope.        
flux_vis=flux_vis[::-1]





    #Our IR
PATH2='\\Users\\Jonathan\\Documents\\ASTR499\\Lab4\\spectra\\'
os.chdir(PATH2)
filename='SiriusIR2.fit'
hdulist=fits.open(filename)
scidata = hdulist[0].data
tot_IR=np.zeros(scidata.shape[1])

for row in scidata:
    tot_IR+=row
    
y=np.arange(1530)
y_cal_IR=y*fit_IR[0]+fit_IR[1]

    #Our Vis
PATH2='\\Users\\Jonathan\\Documents\\ASTR499\\Lab4\\spectra\\'
os.chdir(PATH2)
filename='Sirius2.fit'
hdulist=fits.open(filename)
scidata = hdulist[0].data
tot_Vis=np.zeros(scidata.shape[1])

for row in scidata:
    tot_Vis+=row
    
y=np.arange(1530)
y_cal_Vis=y*fit_Vis[0]+fit_Vis[1]

#Flux Calibrations

i=0
j=0
good_lines_theo_vis=[]
good_lines_exp_vis=[]

for i in range(len(flux_vis)):
    for j in range(len(y_cal_Vis)): 
        if abs(wavelength_vis[i]-y_cal_Vis[j])<0.2453:
            good_lines_theo_vis.append(flux_vis[i])
            good_lines_exp_vis.append(tot_Vis[i])



#IR

i=0
j=0
good_lines_theo_ir=[]
good_lines_exp_ir=[]
good_wavelengths=[]

for i in range(len(flux_ir)):
    for j in range(len(y_cal_IR)): 
        if abs(wavelength_ir[i]-y_cal_IR[j])<3.025:
            good_lines_theo_ir.append(flux_ir[i])
            good_lines_exp_ir.append(tot_IR[i])
            good_wavelengths.append(y_cal_IR[i])
            
good_lines_exp_vis=np.array(good_lines_exp_vis)
good_lines_exp_ir=np.array(good_lines_exp_ir)

sensitivity_IR=good_lines_exp_ir/good_lines_theo_ir
sensitivity_Vis=good_lines_exp_vis/good_lines_theo_vis

print(len(sensitivity_Vis))
print(len(sensitivity_IR))
print(len(good_wavelengths))


plt.savefig(PATH4+ 'Sirius_comparison_graph.pdf')
plt.show()

