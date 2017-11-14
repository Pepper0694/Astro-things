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
from scipy.optimize import curve_fit
import time
from scipy.signal import savgol_filter


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


#plt.plot(y,tot)
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
    

#plt.plot(y,tot)
#plt.show()

y_cal_ir=y*fit_IR[0]+fit_IR[1]

PATH3='\\Users\\Jonathan\\Documents\\ASTR499\\Lab4\\objs_to_run\\'
os.chdir(PATH3)
for filename in glob.iglob('*.fit'):
    
    tot=np.zeros(scidata.shape[1])
    
    for row in scidata:
        tot+=row
        
    y=np.arange(1530)
    
    y=np.array(y)
    #A check to see if the spectrum is in the visible or IR. The dual purpose is to name the graphs with obj later.
    obj_name=os.path.splitext(filename)[0]
    obj=re.sub(r'\W+', '', obj_name)

    if 'IR' in obj:
        y_spec_cal=y_cal_ir
        #tot_IR2_cal=(tot_IR2*flux_cal_ir[1]+flux_cal_ir[0])/1e10
        #tot_IR2_cal=np.array(tot_IR2_cal)[::-1]
        
    #I specified else instead of elif since there were a couple of objects without Vis or IR in the filenames.
    else:
        y_spec_cal=y_cal_vis
        #tot_cal=tot*flux_cal_vis[0]+flux_cal_vis[1]
    #Code to find local maxima checking in 25 angstrom intervals
        
    tot_new=savgol_filter(tot, len(tot)-1, len(tot)/2)
    plt.plot(y_spec_cal, tot)
    #plt.plot(y_spec_cal, tot_new, 'r')
    plt.title('Spectrum of'+ ' '+ obj)
    plt.xlabel('Wavelength ($\AA$)')
    plt.ylabel('Total Count')
    #plt.savefig(PATH4+ obj_name  +'_test_graph.png')
    plt.show()
