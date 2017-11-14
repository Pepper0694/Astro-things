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
from scipy.interpolate import UnivariateSpline
import pickle


def flux_cal():
    #HST

    PATH1='\\Users\\Jonathan\\Documents\\ASTR499\\Lab4\\Vis_calibration\\'
    os.chdir(PATH1)
    #filename='alpha_lyr_004.fits'
    #filename='sirius_stis_002.fits'
    filename='sun_reference_stis_002.fits'

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
        if 5460.857541 < wavelength[i] < 6211.125434:
            wavelength_vis.append(wavelength[i])
            flux_vis.append(flux[i])
        elif 5282.39351 <= wavelength[i] <= 14530.80324:
            wavelength_ir.append(wavelength[i])
            flux_ir.append(flux[i])
    print(len(flux_ir))
    #zeroes=np.zeros(1230-179)
    
    
    
    #flux_ir=np.append(flux_ir,zeroes)


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
    y_cal_IR=y*(-6.04866562)+1.45308032e+04

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
    y_cal_Vis=y*(-4.90691886e-01)+6.21112543e+03

    #Plots
    #f, ax = plt.subplots(2,2)

   # ax[0,0].plot(wavelength_vis,flux_vis)
    #ax[0,0].set_title('Sirius_Vis_HST')

#    ax[1,0].plot(wavelength_ir,flux_ir)
#    ax[1,0].set_title('Sirius_IR_HST')


#    ax[0,1].plot(y_cal_Vis,tot_Vis)
 #   ax[0,1].set_title('Sirius_Vis_EXP')

  #  ax[1,1].plot(y_cal_IR,tot_IR)
   # ax[1,1].set_title('Sirius_IR_EXP')

    #plt.savefig(PATH2+ 'Sirius_comparison_graph.png')
    #plt.show()

    #Optimization

    i=0
    j=0
    good_lines_theo_vis=[]
    good_lines_exp_vis=[]
    good_wavelengths_vis=[]
    for i in range(len(flux_vis)):
        for j in range(len(y_cal_Vis)): 
            if abs(wavelength_vis[i]-y_cal_Vis[j])<0.2453:
                good_lines_theo_vis.append(flux_vis[i])
                good_lines_exp_vis.append(tot_Vis[i])
                good_wavelengths_vis.append(y_cal_Vis[i])

    

          
    flux_cal_vis=stats.linregress(good_lines_exp_vis, good_lines_theo_vis)
    flux_cal_vis=flux_cal_vis[:2]

    i=0
    j=0
    good_lines_theo_ir=[]
    good_lines_exp_ir=[]
    good_wavelengths_ir=[]
    rekt=[]
    print(len(flux_ir))
    print(len(y_cal_IR))
    for i in range(len(flux_ir)):
        for j in range(len(y_cal_IR)): 
            if abs(wavelength_ir[i]-y_cal_IR[j])<3.025:
                good_lines_theo_ir.append(flux_ir[i])
                good_lines_exp_ir.append(tot_IR[i])
                good_wavelengths_ir.append(y_cal_IR[i])
                
    #flux_cal_ir=stats.linregress(good_lines_exp_ir, np.log(good_lines_theo_ir))
    #flux_cal_ir=flux_cal_ir[:2]
    #print(flux_cal_ir)
    
    
    #spl_ir=UnivariateSpline(good_lines_exp_ir, good_lines_theo_ir)
    #plt.plot(good_lines_exp_ir, spl_ir(good_lines_exp_ir))
    
    #good_idx=np.arange(0,len(good_lines_exp_ir), 25)
    #print(good_idx)

    #good_lines_theo_ir=np.array(good_lines_theo_ir)
    #good_lines_exp_ir=np.array(good_lines_exp_ir)
    
    #good_lines_theo_ir2=good_lines_theo_ir[good_idx]
    #good_lines_exp_ir2=good_lines_exp_ir[good_idx]
    #w=np.arange(len(good_lines_exp_ir2))
    #plt.plot(w, good_lines_exp_ir2,'ro')
    #plt.plot(w, good_lines_theo_ir2, 'go')
    #plt.show()
    print(len(good_wavelengths))
    print(len(good_lines_exp_ir))
    array= np.array([good_wavelengths_vis, good_lines_exp_vis, good_lines_theo_vis])
    with open('sunVis.pkl', 'wb') as fout:
        pickle.dump(array ,fout)
    return flux_cal_vis, flux_cal_ir
flux_cal()
