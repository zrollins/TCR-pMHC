#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 11:20:46 2021

@author: zrollins
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gromacs
from gromacs.formats import XPM

#MART1, 100ns Occupancy
hb = gromacs.formats.XPM("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/GVA/90/hmap.xpm", reverse=True)
hmap = hb.array
hmap2 = hmap.astype(np.int)
hb_fraction = hb.array.mean(axis=0)


#High Occupancy
h, d, a= [],[],[]
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/GVA/90/h.log") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 3:
            d.append(str(cols[0]))
            h.append(str(cols[1]))
            a.append(str(cols[2]))

#Atoms
h1,d1,a1=[],[],[]
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/GVA/90/h.ndx") as f:
    for line in f:
        if line.startswith('['):
            continue
        searchphrase = '[ hbonds_chain_DE-chain_ABC ]'
        if searchphrase in line:
            continue
        cols = line.split()
        if len(cols) == 3:
            d1.append(float(cols[0]))
            h1.append(float(cols[1]))
            a1.append(float(cols[2]))
            
print(len(h))
print(h1)
print(hmap2.shape)
            
hlabels={'Donor':d,'Donor Atom':d1,'Acceptor':a,'Acceptor Atom':a1,'Occupancy (%)':100*hb_fraction}
df =  pd.DataFrame(hlabels)

#Labelling Functions
def label_donor (row):
    if row['Donor Atom'] <= 4343:
        return 'MHC_'+ row['Donor']
    if row['Donor Atom'] > 4343 and row['Donor Atom'] <= 5983:
        return 'β2M_' + row['Donor']
    if row['Donor Atom'] >5983 and row['Donor Atom'] <= 6095:
        return 'GVA_' + row['Donor']
    if row['Donor Atom'] >6095 and row['Donor Atom'] <= 9102:
        return 'TCRα_' + row['Donor']
    if row['Donor Atom'] >9102 and row['Donor Atom'] <= 12811:
        return 'TCRβ_'+ row['Donor']
    return 'False'

def label_acceptor (row):
    if row['Acceptor Atom'] <= 4343:
        return 'MHC_'+ row['Acceptor']
    if row['Acceptor Atom'] > 4343 and row['Acceptor Atom'] <= 5983:
        return 'β2M_' + row['Acceptor']
    if row['Acceptor Atom'] >5983 and row['Acceptor Atom'] <= 6095:
        return 'GVA_' + row['Acceptor']
    if row['Acceptor Atom'] >6095 and row['Acceptor Atom'] <= 9102:
        return 'TCRα_' + row['Acceptor']
    if row['Acceptor Atom'] >9102 and row['Acceptor Atom'] <= 12811:
        return 'TCRβ_'+ row['Acceptor']
    return 'False'

df['Donor:Acceptor'] = df.apply (lambda row: label_donor(row), axis=1) + ':' + df.apply (lambda row: label_acceptor(row), axis=1)

#Convert map to dataframe with labelling
hmaps= pd.DataFrame(data=hmap2.T, index=df['Donor:Acceptor'])
print(hmaps)
hmaps.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/GVA/90/hmap.csv', header=None)
