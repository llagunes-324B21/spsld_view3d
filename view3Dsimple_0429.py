#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  21 17:41:12 2021

@author: leolagunes
"""

# This script file will analyze the data from SpringSaLaD 
# focusing on the cluster data! 
# Date: 04/21/2021
# Name: Leo Lagunes

# In this simulation, I ran a basic simulation with fixed monomers (no synth/deg)

# This specific script is to analyze a simulation I ran on Zeus 
# I want to plot the spheres in 3D and compare with the pymol pdb outputs

# import needed packages
from __future__ import division # for averages
#from matplotlib import pyplot as plt # for ploting
import numpy as np
import csv # for looping through files
#import glob
#import os
#import sys
# import math
#import regex as re
# from operator import add
# import operator
#from functools import reduce
import matplotlib.pyplot as plt # to plot in 3D
import pandas as pd


# --------------------------------------
# --------------------------------------
# CODE FOR VIEWING SPSLD IN PYMOL
# --------------------------------------
# --------------------------------------


print('------------------------')

print('------------------------')
print('. . . . . . . . . . .')
print('Ready...')
print('. . . . . . . . . . .')
# ===========================================================================================
# GOAL: view SpringSaLad outputs in PyMol to check if intermediates are in correct assembly 
# ===========================================================================================

# OUTLINE: 
    # 1. load in SpSld output .txt file 
    # 2. Plot each sphere of each molecule in 3D  
# ================================================================================
# List of all functions 
# ================================================================================

# convert from .txt to .csv
def convertFile(runN):
    # converst the .txt file to a .csv to make it easier to read and loop through
    # returns the file name 
    # --- 
    # 1. load in SpSld output .txt file 
    read_file = pd.read_csv (r'alpha_1_trial_5.txt')
    read_file.to_csv (r'alpha_f_run_'+str(runN)+'.csv', index=None)
    fileName = 'alpha_f_run_'+str(runN)+'.csv'
    return fileName

def getMolID(line):
    # this function takes in the line and returns the molecule ID and the sphere ID  
    
    # - molecule ID
    molcID_pre = line[1]; molcID_sp = molcID_pre[0:5]; molcID_int=int(molcID_sp[-1])+1; molcID=str(molcID_int)
    # - sphere ID 
    sphrID = molcID_pre[-2::]
    return molcID, sphrID


def getCoords(line):
    # this function will return the coordinates rounded to 4 decimal places 
    allCoords_raw = line[4::]
    allCoords=[]
    for coord in allCoords_raw:
        coord_new = round(float(coord),4)
        allCoords.append(coord_new)
    return allCoords

def getIndiv(line):
    # this function will take in the raw line and return the molecule number, sphere number and coordinates 
    mol_sphID_i= getMolID(line); 
    molID_i=mol_sphID_i[0];sphrID_i = mol_sphID_i[1]
    coords_i = getCoords(line); #print('Coords: ', coords_i)
    
    return molID_i,sphrID_i,coords_i

def getMolecDictionaryCoords(runN):
    # this function retunrs a dictionary of molecules with xyz coords 
    fileName = convertFile(runN)
    raw_coords = open(fileName);raw_coords_t = csv.reader(raw_coords, delimiter = '\t')

    dictionary_Tps = {}; dictionary_molecs={}
    for line in raw_coords_t:
        # print('line: ', line)
        if line[0] == 'SceneNumber': # get time point from here 
            timePt = line[3]; # print('TimePt = ', timePt)
            # for each time point, create a dictionary that holds all the molecules, and for each molecule the (x,y,z) coords seperately in a list of lists 
            # -- save new tp
            dictionary_Tps[timePt]={};  
        elif line[0] == 'ID': # contains coords of each sphere   
            # x_list =[]; y_list =[]; z_list =[]     
            allCoordInf = getIndiv(line); 
            molID_i=allCoordInf[0];coords_i=allCoordInf[2]
            # print('line: ', line)
            # print('Molecule ID: ', molID_i)
            # print('Molc: ', molID_i, ' sphere: ', sphrID_i, ' at coords: ', coords_i)
            xCoord = coords_i[0]; yCoord = coords_i[1]; zCoord = coords_i[2]; 
            # --- add x,y,z coords 
            # -- create the key if not already there 
            if molID_i not in dictionary_molecs: # if not a key yet, create 
                # print('not a key yet')
                dictionary_molecs[molID_i]={'x': [], 'y':[], 'z':[]}
            else: # add the new coords 
                dictionary_molecs[molID_i]['x'].append(xCoord); dictionary_molecs[molID_i]['y'].append(yCoord)
                dictionary_molecs[molID_i]['z'].append(zCoord)
                # dictionary_molecs[molID_i]={'x': [], 'y':[], 'z':[]}
        
        
    return dictionary_molecs



def getColor(indexValue, molecType):
    # this function return the color for the sphere 
    molecColors_alpha=['lime', 'forestgreen']; molecColors_beta=['steelblue', 'deepskyblue']
    if molecType ==1: # alpha
        if (indexValue % 2) == 0: # if even first color 
            colorVal = molecColors_alpha[0]
        else:
            colorVal = molecColors_alpha[1]
    else: # beta
         if (indexValue % 2) == 0: # if even first color 
            colorVal = molecColors_beta[0]
         else:
            colorVal = molecColors_beta[1]

    return colorVal



runN = 0; 

allCoords = getMolecDictionaryCoords(runN); #print(allCoords)
  
# OKAY NOW UNPACK TO PLOT EACH MOLECULE (KEY) IN 3D


ax = plt.axes(projection='3d')
for key in allCoords:
    print('Molecule number: ', key)
    molecType=1 # an alpha for this current code (eventually either 1 or 2)
    indexValue=int(key)-1;colorVal = getColor(indexValue, molecType)
    coordsList = allCoords[key]; # print('coords: ', coordsList)
    xData = coordsList['x']; yData = coordsList['y']; zData = coordsList['z']; 
    ax.scatter3D(xData, yData, zData, c=colorVal);
    
# plot/figure formats
ax.set_xlabel('x coord (nm)'); ax.set_ylabel('y coord (nm)'); ax.set_zlabel('z coord (nm)')
# xData1=16.3093; yData1=-36.3312; zData1=75.1913; ax.scatter3D(xData1, yData1, zData1, c='black')
# ax.set_zlim(0,80)

# ADD ALL MOLECULES + TIME POINTS IN EACH NEW FIGURE 
  
# ==============================
# START OF SCRIPT 
# ==============================

# Fixed variables
runN = 0; 
tEnd = 1; nTimes = 4 # end point and number of time points to output 

#print('Making 3D plot for run ', runN); 

#print('Done with run ', runN)
# NOW FIGURE IT OUT FOR ALL THE TIME POINTS I CHAVE IN THE SPSLD FILE
#  ALSO CREATE A PROGRESS OUTPUT BAR PERCENTAGE 

# practice 3D plot: 

#ax = plt.axes(projection='3d')
# Data for three-dimensional scattered points
#zdata = 15 * np.random.random(100)
#xdata = np.sin(zdata) + 0.1 * np.random.randn(100)
#ydata = np.cos(zdata) + 0.1 * np.random.randn(100)
#ax.scatter3D(xdata, ydata, zdata, c=zdata, cmap='Greens');



print('---')
print('Pinche done bitch!')
# ==============================
# END OF SCRIPT 
# ==============================

# TRY THIS INSTEAD:
    # LOOP FOR EVERY MOLECULE (1-20) AND ADD IN THE COORDIDATES 
    # CREATE A TABLE THAT HAS ATOM # (SITE ID) , MOLECULE # (MOLECULE ID), AND COORDINATES FIRST 
    # THIS WAY I CAN JUST CALL THE ITEM OF THE LIST 









