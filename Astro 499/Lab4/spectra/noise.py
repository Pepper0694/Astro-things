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

start_time = time.time()

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
    print(filename)

    hdulist=fits.open(filename)
    scidata = hdulist[0].data
    
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
    tot_new=savgol_filter(tot, len(tot)-1, 75)
    #plt.plot(y_spec_cal, tot_new, 'r')
    x=50
    spacing=50

    peaks=[]
    count=[]
    idx=[]
    i=0
    spread=[]
    ugh=0
    blah=0
    stds=[]
    max_stds=[]
    for i in range(int(len(tot)/spacing)-1):
        peak_wavelength=y_spec_cal[np.argmax(tot[x-spacing:x])+x]
        
        peaks.append(y_spec_cal[np.argmax(tot[x-spacing:x])+x])
        count=tot[np.argmax(tot[x-spacing:x])+x]
        #y_spec_cal=list(y_spec_cal)
        #for wavelength in y_spec_cal[x-spacing:x]:
        #    print(len(y_spec_cal[x-spacing:x]))
        #    if wavelength!=peak_wavelength:
        #        idx=y_spec_cal.index(wavelength)
        #        tot[idx]=0
        #    else:
        #        print('You are cool'+str(ugh))
        #        ugh+=1
        #y_spec_cal=np.array(y_spec_cal)

        for a in y_spec_cal:
            if peak_wavelength-spacing <a< peak_wavelength+spacing:
                spread.append(a)
        y_spec_cal=list(y_spec_cal)
        spread_idx=[]
        for wavelength in spread:
            spread_idx.append(y_spec_cal.index(wavelength))
        print(spread)
        max_stds=[]    
        for idx in spread_idx:
            stds.append((tot[idx]-np.mean(tot[spread_idx]))/np.std(tot[spread_idx]))
        #stds=list(stds)    
        #max_stds=np.max(stds)   
        #max_stds=[]
        max_stds.append(np.max(stds))
        max_stds_ind=np.max(stds)
      #  print(max_stds_ind)
        #print(np.max(max_stds))
        #print(np.min(max_stds))
        #print(np.mean(max_stds))
        #print(max_stds)
        #print(len(max_stds))
        not_noise=[]
        not_noise_wavelength=[]
        for idx in spread_idx:
           # print(max_stds_ind)
            if tot[idx]>np.mean(tot[spread_idx])+ (2.00)*np.std(tot[spread_idx]):   # old value (max_stds_ind)
                not_noise.append(tot[idx])
                not_noise_wavelength.append(y_spec_cal[idx])
        
        
        y_spec_cal=np.array(y_spec_cal)             
        i+=1
        x+=spacing
        #print(max_stds)
        #print('beep'+str(ugh))
        ugh+=1
        spread=list(spread)
        del spread[:]
        spread_idx=list(spread_idx)
        del spread_idx[:]
        stds=list(stds)
        del stds[:]
    #flux_count=[]
    #tot=list(tot)
    #for a in tot:
        #print('boop'+str(blah))
        #if a <np.mean(tot)+np.mean(max_stds)*np.std(tot):
            #if a <(0.5)*np.max(tot)+(1.0)*np.mean(max_stds)*np.std(tot):
        #    tot=list(tot)
        #    tot[tot.index(a)]=0
               # flux_count.append(a)
                #peaks.append(y_spec_cal[tot.index(a)])
        #    if a!=0:
        #        flux_count.append(a)
        #        peaks.append(y_spec_cal[tot.index(a)])
            #blah+=1    
    tot=np.array(tot)
    #print(len(max_stds))
    #print(max_stds)
    #print(max(max_stds))
    #print(min(max_stds))
    max_stds=list(max_stds)
    del max_stds[:]
    #highest=tot[np.argmax(tot)]
    #tot=tot/highest
    #print(peaks)
    #print(noise_wavelength)
    print(not_noise_wavelength)
    plt.scatter(not_noise_wavelength, not_noise)
    #plt.plot(y_spec_cal, tot)
    plt.title('Noise Spectrum of'+ ' '+ obj)
    #if 'IR' in obj and obj!='Sirius2':
    #    plt.xlim(4000,16000)
    #else:
    #    plt.xlim(5400,6300)
    if len(not_noise)>0:
        plt.ylim(0,max(not_noise)+1000)
    plt.xlabel('Wavelength ($\AA$)')
    plt.ylabel('Total Count')
    plt.savefig(PATH4+ obj_name  +'_not_noise_graph.png')
    plt.show()
print("--- %s seconds ---" % (time.time() - start_time))

