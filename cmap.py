#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 16:54:02 2021

@author: zrollins
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gromacs
from gromacs.formats import XPM

#MART1, 100ns Occupancy
hb = gromacs.formats.XPM("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/L1/90/contacts.xpm", reverse=True)
hmap = hb.array
hmap2 = hmap.astype(np.int)
hb_fraction = hb.array.mean(axis=0)


#High Occupancy
h, d, a= [],[],[]
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/L1/90/contacts.log") as f:
    for line in f:
        if line.startswith('#'):
            continue
        if line.startswith('@'):
            continue
        cols = line.split()

        if len(cols) == 2:
            d.append(str(cols[0]))
            a.append(str(cols[1]))

#Atoms
h1,d1,a1=[],[],[]
with open("/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/L1/90/contacts.ndx") as f:
    for line in f:
        if line.startswith('['):
            continue
        searchphrase = '[ contacts_chain_DE-chain_ABC ]'
        if searchphrase in line:
            continue
        
        cols = line.split()
        if len(cols) == 2:
            d1.append(float(cols[0]))
            a1.append(float(cols[1]))
            
print(len(a))
print(len(a1))
print(hmap2.shape)
#d1 = d1[1:]
#a1 = a1[1:]

print(len(a1))
#print(d1)

            
hlabels={'Donor':d,'Donor Atom':d1,'Acceptor':a,'Acceptor Atom':a1,'Occupancy (%)':100*hb_fraction}
df =  pd.DataFrame(hlabels)

#Labelling Functions
def label_donor (row):
    if row['Donor Atom'] <= 4343:
        return 'MHC_'+ row['Donor']
    if row['Donor Atom'] > 4343 and row['Donor Atom'] <= 5983:
        return 'β2M_' + row['Donor']
    if row['Donor Atom'] >5983 and row['Donor Atom'] <= 6116:
        return 'L1_' + row['Donor']
    if row['Donor Atom'] >6116 and row['Donor Atom'] <= 9123:
        return 'TCRα_' + row['Donor']
    if row['Donor Atom'] >9123 and row['Donor Atom'] <= 12832:
        return 'TCRβ_'+ row['Donor']
    return 'False'

def label_acceptor (row):
    if row['Acceptor Atom'] <= 4343:
        return 'MHC_'+ row['Acceptor']
    if row['Acceptor Atom'] > 4343 and row['Acceptor Atom'] <= 5983:
        return 'β2M_' + row['Acceptor']
    if row['Acceptor Atom'] >5983 and row['Acceptor Atom'] <= 6116:
        return 'L1_' + row['Acceptor']
    if row['Acceptor Atom'] >6116 and row['Acceptor Atom'] <= 9123:
        return 'TCRα_' + row['Acceptor']
    if row['Acceptor Atom'] >9123 and row['Acceptor Atom'] <= 12832:
        return 'TCRβ_'+ row['Acceptor']
    return 'False'

df['Donor:Acceptor'] = df.apply (lambda row: label_donor(row), axis=1) + ':' + df.apply (lambda row: label_acceptor(row), axis=1)

#Convert map to dataframe with labelling
hmaps= pd.DataFrame(data=hmap2.T, index=df['Donor:Acceptor'])
print(hmaps)
hmaps.to_csv('/Users/zrollins/Documents/Documents/DMF5_MART1/rndm/L1/90/cmap.csv', header=None)
