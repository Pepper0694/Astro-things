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
from scipy.signal import savgol_filter

def redchisqg(ydata,ymod,deg=2,sd=None):
    # Chi-square statistic  
    if sd==None:
        chisq=np.sum((ydata-ymod)**2)
    else:
        chisq=np.sum( ((ydata-ymod)/sd)**2 )
             
      # Number of degrees of freedom assuming 2 free parameters  
    nu=ydata.size-1-deg
        
    return chisq/nu
def norm(array):
    array=np.array(array)
    return array/np.max(array)

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
    
x=np.arange(1530)
x=x[::-1]

    #This block makes a list of maxima
L=np.argsort(-tot)
i=0
bad_idx=[]
    #To remove multiple maxima for the same peak. 10 is the critical number, since 5 didn't work
for i in range(len(L)-1):
     if abs(L[i]-L[i+1])<10:
         bad_idx.append(i+1)
L=np.delete(L,bad_idx)

    #Get wavelengths of 3 highest peaks
peaks=x[L[:3]]
    

theo_Vis=[5460, 5769,5790]
theo_Vis=theo_Vis[::-1]


#plt.plot(x,tot)
#plt.show()

    #Find the linear fit
fit_Vis=np.polyfit(peaks,theo_Vis,1)
    

x_cal_vis=x*fit_Vis[0]+fit_Vis[1]
    

    #IR-Calibration

PATH2='\\Users\\Jonathan\\Documents\\ASTR499\\Lab4\\IR_calibration\\'
os.chdir(PATH2)
filename='neon6.fit'
hdulist=fits.open(filename)
scidata = hdulist[0].data
tot=np.zeros(scidata.shape[1])



for row in scidata:
    tot+=row
x=np.arange(1530)
x=x[::-1]

L=np.argsort(-tot)

i=0
for i in range(len(L)-1):
    if abs(L[i]-L[i+1])<10:
        bad_idx.append(i+1)
L=np.delete(L,bad_idx)
peaks=x[L[:3]]
    

theo_IR=[5852, 6402,6678]
theo_IR=theo_IR[::-1]

peaks=np.array(peaks)

fit_IR=np.polyfit(peaks,theo_IR,1)
    

#plt.plot(x,tot)
#plt.show()

x_cal_ir=x*fit_IR[0]+fit_IR[1]

    


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
print(min(wavelength))
print(max(wavelength))
plt.plot(wavelength,flux)
plt.title('Sirius Spectrum')
plt.xlim(0,12000)
plt.xlabel('Wavelength ($\AA$)')
plt.show()

i=0
for i in range(len(flux)):
    if min(x_cal_vis) < wavelength[i] < max(x_cal_vis):
        wavelength_vis.append(wavelength[i])
        flux_vis.append(flux[i])
    elif min(x_cal_ir) < wavelength[i] < max(x_cal_ir):
        wavelength_ir.append(wavelength[i])
        flux_ir.append(flux[i])

    #For some reason the slope was in the wrong direction. It would be
    #discontinuous if you didn't flip the slope.        
flux_vis=flux_vis[::-1]
print(len(flux_vis))
print(len(flux_ir))
filename='alpha_lyr_004.fits'


table = Table.read(filename)
df2 = table.to_pandas()
wavelength2=[]
wavelength_ir2=[]
wavelength_vis2=[]
flux_ir2=[]
flux_vis2=[]

flux2=df2['FLUX']
wavelength2=df2['WAVELENGTH']
plt.plot(wavelength2,flux2)
plt.title('Vega Spectrum')
plt.xlabel('Wavelength ($\AA$)')
plt.show()
i=0
for i in range(len(flux2)):
    if min(x_cal_vis) < wavelength2[i] < max(x_cal_vis):
        wavelength_vis2.append(wavelength2[i])
        flux_vis2.append(flux2[i])
    elif min(x_cal_ir) < wavelength2[i] < max(x_cal_ir):
        wavelength_ir2.append(wavelength2[i])
        flux_ir2.append(flux2[i])

    #For some reason the slope was in the wrong direction. It would be
    #discontinuous if you didn't flip the slope.        
flux_vis=flux_vis[::-1]
print(len(flux_vis2))
print(len(flux_ir2))

filename='sun_reference_stis_002.fits'


table = Table.read(filename)
df3 = table.to_pandas()
wavelength3=[]
wavelength_ir3=[]
wavelength_vis3=[]
flux_ir3=[]
flux_vis3=[]

flux3=df3['FLUX']
wavelength3=df3['WAVELENGTH']
print(min(wavelength3))
print(max(wavelength3))
plt.plot(wavelength3,flux3)
plt.title('Sun Spectrum')
plt.xlim(0,12000)
plt.xlabel('Wavelength ($\AA$)')
plt.show()

i=0
for i in range(len(flux3)):
    if min(x_cal_vis) < wavelength3[i] < max(x_cal_vis):
        wavelength_vis3.append(wavelength3[i])
        flux_vis3.append(flux3[i])
    elif min(x_cal_ir) < wavelength3[i] < max(x_cal_ir):
        wavelength_ir3.append(wavelength3[i])
        flux_ir3.append(flux3[i])

    #For some reason the slope was in the wrong direction. It would be
    #discontinuous if you didn't flip the slope.        
flux_vis=flux_vis[::-1]
print(len(flux_vis3))
print(len(flux_ir3))


    #Our IR
PATH2='\\Users\\Jonathan\\Documents\\ASTR499\\Lab4\\spectra\\'
os.chdir(PATH2)
filename='Sirius_IR2.fit'
hdulist=fits.open(filename)
scidata = hdulist[0].data
tot_IR=np.zeros(scidata.shape[1])

for row in scidata:
    tot_IR+=row
    
x=np.arange(1530)
x_cal_IR=x*fit_IR[0]+fit_IR[1]

    #Our Vis
PATH2='\\Users\\Jonathan\\Documents\\ASTR499\\Lab4\\spectra\\'
os.chdir(PATH2)
filename='Sirius2.fit'
hdulist=fits.open(filename)
scidata = hdulist[0].data
tot_Vis=np.zeros(scidata.shape[1])

for row in scidata:
    tot_Vis+=row
    
x=np.arange(1530)
x_cal_Vis=x*fit_Vis[0]+fit_Vis[1]

#Flux Calibrations

i=0
j=0
good_lines_theo_vis=[]
good_lines_exp_vis=[]
good_wavelengths_vis=[]
idx_vis=[]
for i in range(len(flux_vis)):
    for j in range(len(x_cal_Vis)): 
        if abs(wavelength_vis[i]-x_cal_Vis[j])<0.2453:
            good_lines_theo_vis.append(flux_vis[i])
            good_lines_exp_vis.append(tot_Vis[i])
            good_wavelengths_vis.append(x_cal_Vis[i])
            idx_vis.append(i)
########################################################################
#print(len(good_wavelengths_vis)); print(len(wavelength_vis)); print(len(x_cal_Vis))          
#Linear           
flux_cal_vis_linear=stats.linregress(good_lines_exp_vis, good_lines_theo_vis)
flux_cal_vis_linear=flux_cal_vis_linear[:2]

#Exponential
flux_cal_vis_exponential2=stats.linregress(good_lines_exp_vis, np.log(good_lines_theo_vis))
flux_cal_vis_exponential2=flux_cal_vis_exponential2[:2]
flux_cal_vis_exponential = []
for item in flux_cal_vis_exponential2:
    flux_cal_vis_exponential.append(item)
#Machine Learning



#IR

i=0
j=0
good_lines_theo_ir=[]
good_lines_exp_ir=[]
good_wavelengths_ir=[]
idx_ir=[]
for i in range(len(flux_ir)):
    for j in range(len(x_cal_IR)): 
        if abs(wavelength_ir[i]-x_cal_IR[j])<3.025:
            good_lines_theo_ir.append(flux_ir[i])
            good_lines_exp_ir.append(tot_IR[i])
            good_wavelengths_ir.append(x_cal_IR[i])
            idx_ir.append(i)
            
good_lines_exp_vis=np.array(good_lines_exp_vis)
good_lines_exp_ir=np.array(good_lines_exp_ir)

#print(len(good_lines_exp_ir))
#print(len(good_lines_exp_ir))
#print(len(good_wavelengths_ir))

#Frank's Idea
sensitivity_IR=good_lines_exp_ir/good_lines_theo_ir
sensitivity_Vis=good_lines_exp_vis/good_lines_theo_vis
#print(sensitivity_IR)
#print(sensitivity_Vis)

#Linear
flux_cal_ir_linear=stats.linregress(good_lines_exp_ir, good_lines_theo_ir)
flux_cal_ir_linear=flux_cal_ir_linear[:2]

            
#Exponential            
flux_cal_ir_exponential2=stats.linregress(good_lines_exp_ir, np.log(good_lines_theo_ir))
flux_cal_ir_exponential2=flux_cal_ir_exponential2[:2]
flux_cal_ir_exponential=[]
for item in flux_cal_ir_exponential2:
    flux_cal_ir_exponential.append(item)
   
#Machine Learning
flux_cal_ir_ml=pickle.load(open('values.p','rb'))



#Linear
tot_Vis_linear=tot_Vis*flux_cal_vis_linear[0]+flux_cal_vis_linear[1]
tot_IR_linear=tot_IR*flux_cal_ir_linear[0]+flux_cal_ir_linear[1]
tot_IR_linear=tot_IR_linear*-1
#Exponential
tot_IR2_cal_exponential=(good_lines_exp_ir*flux_cal_ir_exponential[1]+flux_cal_ir_exponential[0])/1e10
tot_IR2_cal_exponential=np.array(tot_IR2_cal_exponential)[::-1]
tot_IR2_cal_exponential=tot_IR2_cal_exponential*-1
tot_Vis2_cal_exponential=(good_lines_exp_vis*flux_cal_vis_exponential[1]+flux_cal_vis_exponential[0])/1e10
tot_Vis2_cal_exponential=np.array(tot_Vis2_cal_exponential)[::-1]



#Machine Learning
#tot_Vis_cal_ml=
tot_IR_cal_ml=good_lines_exp_ir*flux_cal_ir_ml[1]+flux_cal_ir_ml[0]
tot_IR_cal_ml=tot_IR_cal_ml[::-1]

#Plots
f, ax = plt.subplots(2,5)
######################################################################
ax[0,0].plot(wavelength_vis,flux_vis)
ax[0,0].set_title('Sirius Visible HST',fontsize=12)
ax[0,0].tick_params(labelsize=8)


ax[1,0].plot(wavelength_ir,flux_ir)
ax[1,0].set_title('Sirius IR HST',fontsize=12)
ax[1,0].tick_params(labelsize=8)


ax[0,1].plot(x_cal_Vis,tot_Vis)
ax[0,1].set_title('Sirius Visible Experiment',fontsize=12)
ax[0,1].tick_params(labelsize=8)

tot_new1=savgol_filter(tot_IR, len(tot_IR)-1, int(len(tot_IR)/75))
ax[1,1].plot(x_cal_IR, tot_new1, 'r')
ax[1,1].set_ylim(min(tot_IR)-10000, max(tot_IR)+10000)
ax[1,1].plot(x_cal_IR,tot_IR)
ax[1,1].set_title('Sirius IR Experiment',fontsize=12)
ax[1,1].tick_params(labelsize=8)
#print(len(tot_Vis_linear))

#print(len(tot_Vis_linear));print(len(flux_vis))
tot_Vis_linear2=tot_Vis_linear[idx_vis]
rcs_lin_vis=redchisqg(norm(tot_Vis_linear2), norm(good_lines_theo_vis))
rcs_lin_vis="{0:.3f}".format(rcs_lin_vis)
ax[0,2].plot(x_cal_Vis,tot_Vis_linear)
ax[0,2].set_title('Sirius Vis Linear Calibration',fontsize=12)
ax[0,2].tick_params(labelsize=8)
ax[0,2].text( 5800, -0.5e-8, '$\chi ^2$ = '+str(rcs_lin_vis))

tot_Vis2_cal_exponential2=tot_Vis2_cal_exponential[idx_vis]
rcs_exp_vis=redchisqg(norm(tot_Vis2_cal_exponential2), norm(good_lines_theo_vis))
rcs_exp_vis="{0:.3f}".format(rcs_exp_vis)
ax[0,3].plot(wavelength_vis,tot_Vis2_cal_exponential)
ax[0,3].set_title('Sirius Vis Exponential Calibration', fontsize=12)
ax[0,3].tick_params(labelsize=8)
ax[0,3].text( 5900, -0.000545, '$\chi ^2$ = '+ str(rcs_exp_vis))

#rcs_ml_vis=redchisqg(tot_Vis_cal_ml, flux_vis)
#ax[0,4].plot(good_wavelengths,tot_Vis_cal_ml)
#ax[0,4].set_title('Calibrated with machine learning Sirius Vis')
#ax[0,4].tick_params(labelsize=6)
#ax[0,4].text( , , str(rcs_ml_vis))

tot_IR_linear2=tot_IR_linear[idx_ir]
rcs_lin_ir=redchisqg(norm(tot_IR_linear2), norm(good_lines_theo_ir))
rcs_lin_ir="{0:.3f}".format(rcs_lin_ir)
ax[1,2].plot(x_cal_IR,tot_IR_linear)
ax[1,2].set_title('Sirius IR Linear Calibration',fontsize=12)
ax[1,2].tick_params(labelsize=8)
ax[1,2].text( 11000,-2e-9 ,'$\chi ^2$ = '+ str(rcs_lin_ir))

tot_IR2_cal_exponential2=tot_IR2_cal_exponential[idx_ir]
rcs_exp_ir=redchisqg(norm(tot_IR2_cal_exponential2), norm(good_lines_theo_ir))
rcs_exp_ir="{0:.3f}".format(rcs_exp_ir)
ax[1,3].plot(wavelength_ir,tot_IR2_cal_exponential)
ax[1,3].set_title('Sirius IR Exponential Calibration',fontsize=12)
ax[1,3].tick_params(labelsize=8)
ax[1,3].text(10000 , 0.000105,'$\chi ^2$ = '+ str(rcs_exp_ir))

#print(len(tot_IR_cal_ml)); print(len(flux_ir))

rcs_ml_ir=redchisqg(norm(tot_IR_cal_ml), norm(good_lines_theo_ir))
rcs_ml_ir="{0:.3f}".format(rcs_ml_ir)
ax[1,4].plot(good_wavelengths_ir,tot_IR_cal_ml)
ax[1,4].set_title('Calibrated with machine learning Sirius IR EXP',fontsize=12)
ax[1,4].tick_params(labelsize=8)
ax[1,4].text( 11000, 90,'$\chi ^2$ = '+ str(rcs_ml_ir))
plt.savefig(PATH4+ 'Sirius_comparison_graph.pdf')
plt.show()



#Glob allowed me to loop through the directory I created for the data
PATH3='\\Users\\Jonathan\\Documents\\ASTR499\\Lab4\\spectra\\'
os.chdir(PATH3)
for filename in glob.iglob('*.fit'):
    

    hdulist=fits.open(filename)
    scidata = hdulist[0].data
    
    tot=np.zeros(scidata.shape[1])
    
    for row in scidata:
        tot+=row
        
    x=np.arange(1530)
    
    x=np.array(x)

    #A check to see if the spectrum is in the visible or IR. The dual purpose is to name the graphs with obj later.
    obj_name=os.path.splitext(filename)[0]
    obj=re.sub(r'\W+', '', obj_name)

    
    if 'IR' in obj:   #Fix Sirius PPPPPRRROOOBBLLEEMMM
        x_spec_cal=x*fit_IR[0]+fit_IR[1]
        i=0
        j=0
        good_lines_theo_ir=[]
        good_lines_exp_ir=[]
        good_wavelengths_ir=[]
        
        
        frank_good_ir=[]
        for i in range(len(flux_ir)):           #len 1230
            for j in range(len(x_spec_cal)):    #len 1530
                if abs(wavelength_ir[i]-x_spec_cal[j])<3.025:
                    good_lines_theo_ir.append(flux_ir[i])
                    good_lines_exp_ir.append(tot[i])
                    good_wavelengths_ir.append(x_spec_cal[i])
                    frank_good_ir.append(sensitivity_IR[i])
        frank_good_ir=np.array(frank_good_ir)           
        good_lines_exp_vis=np.array(good_lines_exp_vis)
        good_lines_exp_ir=np.array(good_lines_exp_ir)
        #tot_IR2_cal=(tot*flux_cal_ir_exponential[1]+flux_cal_ir_exponential[0])/1e10
        #tot_IR2_cal=np.array(tot_IR2_cal)[::-1]
        tot_IR2_cal_ml=good_lines_exp_ir*flux_cal_ir_ml[1]+flux_cal_ir_ml[0]
        tot_IR_cal_ml=tot_IR2_cal_ml[::-1]

        #Frank
        good_lines_exp_ir=np.array(good_lines_exp_ir)    
        frank_applied_ir=good_lines_exp_ir/frank_good_ir
    #I specified else instead of elif since there were a couple of objects without Vis or IR in the filenames.
    else:
        x_spec_cal=x*fit_Vis[0]+fit_Vis[1]
        tot_cal=tot*flux_cal_vis_linear[0]+flux_cal_vis_linear[1]
        i=0
        j=0
        good_lines_theo_vis=[]
        good_lines_exp_vis=[]
        good_wavelengths_vis=[]
        frank_good_vis=[]
        for i in range(len(flux_vis)):
            for j in range(len(x_cal_Vis)): 
                if abs(wavelength_vis[i]-x_cal_Vis[j])<0.2453:
                    good_lines_theo_vis.append(flux_vis[i])
                    good_lines_exp_vis.append(tot_Vis[i])
                    good_wavelengths_vis.append(x_cal_Vis[i])
                    frank_good_vis.append(sensitivity_Vis[i])
       #Linear           
        flux_cal_vis_linear=stats.linregress(good_lines_exp_vis, good_lines_theo_vis)
        flux_cal_vis_linear=flux_cal_vis_linear[:2]

        #Exponential
        flux_cal_vis_exponential2=stats.linregress(good_lines_exp_vis, np.log(good_lines_theo_vis))
        flux_cal_vis_exponential2=flux_cal_vis_exponential2[:2]
        flux_cal_vis_exponential = []
        for item in flux_cal_vis_exponential2:
            flux_cal_vis_exponential.append(item)
        #Frank
        good_lines_exp_vis=np.array(good_lines_exp_vis)    
        frank_applied_vis=good_lines_exp_vis/frank_good_vis   
    #Savgol smooth
    tot_new=savgol_filter(tot_cal, len(tot)-1, int(len(tot)/75))
    
    
    #Code to find local maxima checking in 25 angstrom intervals
    w=25
    spacing=25

    peaks=[]
    count=[]
    i=0
    print(x)
    for i in range(int(len(tot)/spacing)-1):
   
       peaks.append(x_spec_cal[np.argmax(tot[w-spacing:w])+w])
       count.append(tot[np.argmax(tot[w-spacing:w])+w])
       i+=1
       w+=spacing
    print(peaks)
    #Sort by showing the peaks with the highest counts first.
    count=np.array(count)
    peaks=np.array(peaks)
    
    ind_sort=np.argsort(count)

    peaks=peaks[ind_sort]
    count=count[ind_sort]
    
    count=count[::-1]
    
    #x=x[::-1]
    if '_IR' in obj:
        f, ax = plt.subplots(1,3)
        if 'm42_IR' in obj:
            tot_cal=tot_cal[::-1]
        
        ax[0].plot(x_spec_cal, tot_cal)
        ax[0].set_title('Spectrum of'+ ' '+ obj)
        ax[0].set_xlabel('Wavelength ($\AA$)')
        ax[0].set_ylabel('Total Count')
        ax[0].tick_params(labelsize=10)
        tot_new=tot_new[::-1]
        ax[0].plot(x_spec_cal, tot_new, 'r')
        tot_IR_cal_ml=tot_IR_cal_ml[::-1]

        if 'vega_IR2' in obj:
            tot_IR_cal_ml=tot_IR_cal_ml[::-1]
        
            
        ax[1].scatter(good_wavelengths_ir,tot_IR_cal_ml, s=5)
        #ax[1].plot(x_spec_cal, tot_new, 'r')
        ax[1].set_ylim(0,max(tot_IR_cal_ml))
        ax[1].set_title('Calibrated with machine learning '+obj+' IR EXP',fontsize=12)
        ax[1].set_xlabel('Wavelength ($\AA$)')
        ax[1].set_ylabel('Total Count')
        ax[1].tick_params(labelsize=10)
        frank_applied_ir=good_lines_exp_ir*frank_applied_ir
        frank_applied_ir=frank_applied_ir[::-1]
        
        #print(min(good_wavelengths_ir)); print(max(good_wavelengths_ir))
        
        ax[2].plot(good_wavelengths_ir,frank_applied_ir)
        ax[2].set_title('Spectrum of'+ ' '+ obj+'with QE divided out')
        ax[2].set_xlabel('Wavelength ($\AA$)')
        ax[2].set_ylabel('Total Count/QE')
        ax[2].tick_params(labelsize=10)
    else:
        
        f, ax = plt.subplots(1,2)
        #print(min(good_wavelengths_vis)); print(max(good_wavelengths_vis))
        ax[0].plot(x_spec_cal, tot_cal)
        ax[0].plot(x_spec_cal, tot_new, 'r')
        ax[0].set_title('Spectrum of '+  obj)
        ax[0].set_xlabel('Wavelength ($\AA$)')
        ax[0].set_ylabel('Total Count')
        ax[0].tick_params(labelsize=10)
        #print(len(good_wavelengths_vis)); print(len(frank_applied_vis))
        frank_applied_vis=good_lines_exp_vis*frank_applied_vis
        frank_applied_vis=frank_applied_vis[::-1]
        ax[1].plot(good_wavelengths_vis,frank_applied_vis)
        ax[1].set_title('Spectrum of'+ ' '+ obj+'with QE divided out')
        ax[1].set_xlabel('Wavelength ($\AA$)')
        
        ax[1].set_ylabel('Total Count/QE')
        ax[1].tick_params(labelsize=10)

        
    #Creating an ascii file to read off the peaks easily
    fmt={'count': '%i', 'wavelength': '%0.2f'}
    output={'count': count, 'wavelength':peaks}
    ascii.write(output,PATH3+ obj_name +'_calibrated_ascii', names=['count','wavelength'], formats=fmt) 
    
    
    plt.savefig(PATH4+ obj_name  +'_calibrated_graph.png')
    plt.show()
