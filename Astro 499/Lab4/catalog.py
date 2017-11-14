import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import glob
import re
from astropy.io import ascii
import csv
import pandas as pd
PATH='\\Users\\Jonathan\\Documents\\ASTR499\\Lab4\\'

name=[]
min_mag=[]
max_mag=[]
period=[]
location=[]

header=['names', 'location', 'max mag', 'min mag', 'period']
df=pd.read_csv('gcvs_cat.dat', sep='|', header=None, usecols=[1,2,4,5,9], names=header)
print(df)
df=df.dropna(axis=0, inplace=False)
df=df.sort_values(by=['period'])

df.to_csv('catalogs.csv', header=header, index=False)




          
