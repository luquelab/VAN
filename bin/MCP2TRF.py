# -*- coding: utf-8 -*-
"""
 Author:
    Diana Y. Lee, Luque Lab, SDSU
    dlee@sdsu.edu

 Purpose:
    Base random forest model for predicting the capsid architecture (as measured by the 
    T-number) of a tailed phage from the MCP sequence

 Requires: 

    phage_functions.ipynb  :  Functions for calculating T based on genome size
    MCP2T_RF_state.db  :  Trained random forest model database

Input: 
    MCP phage data to analyze. Must include the following columns:
        'genome_length'
        'MCP_len'
        'Virus_ID'
        'MCP_Sequence'
        'IPC'

Output:
    results\phageResults.csv  :  results of the random forest prediction
    results\phageResults_db.db  :  database state after the prediction
"""

# basic imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
np.random.seed(42)
import csv
import os
import copy

# ML imports
from sklearn.metrics import mean_squared_error
from sklearn import metrics
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.base import BaseEstimator
from sklearn.metrics import precision_score, recall_score, f1_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import mean_squared_error
from sklearn.metrics import accuracy_score
from sklearn.model_selection import cross_val_score, GridSearchCV
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import cross_val_score
from sklearn import metrics

# biopython imports
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC

# phage function imports
from phageFunctions import tNearest
from phageFunctions import tNearestValid
from phageFunctions import tModel
from phageFunctions import tNum
from phageFunctions import tList
from phageFunctions import tDictAll

# custom function to count amino acids
# amino acids are hardcoded to avoid broken dependencies, since they do not change
def createFreq(acidSeq, normF=None):
    normF = normF or 0
    if (normF > 1):
        print("Valid tTypes are 0 (raw) or 1 (normalized by sample). Defaults to 0.")
        return
    AA = []
    aaList = np.asarray(['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'])
    aaLen=len(aaList)
    n = len(acidSeq)
    for i in range(n):
        for acid in (aaList):
            trans1=str(acidSeq[i])
            a = trans1.count(str(acid))
            AA.append(a)
    rFreq = np.asarray(AA).reshape((n, aaLen))
    if (normF == 0):
#        print("Success! Created an nx20 array, where n is the length of the list provided:",n)
#        print("Columns are frequency totals for each amino acid:",aaList)
        return rFreq
    if (normF == 1):
        nFreq = copy.copy(rFreq).astype(float)
        fff3 = copy.copy(rFreq).astype(float)
        nf = rFreq.shape[1]
        for i in range(nf):
            nFreq[:,i] = fff3[:,i]/fff3.sum(axis=1)
#        print("Success! Created an nx20 array, where n is the length of the list provided:",n)
#        print("Columns are frequency percentages for each amino acid:",aaList)
        return nFreq

## load kernel state
## rerun imports cell as well
import dill
dill.load_session('MCP2T_RF_state.db')

# set the error margin
errMar = 0.09

MCP_Type = input("Get in the VAN! Enter MCP input type (1: fasta file, 2: csv file): ")

if MCP_Type=="1":
    print("Sad Honk. Sorry, I'm currently in the shop.")

elif MCP_Type=="2":
    print("Vroom! Let's go!")
    print("Your .csv file will require four columns: Virus ID, MCP_Sequence, MCP_len, and IPC")
    MCP_File_Loc = input("Enter file location: ")
    
    assert os.path.exists(MCP_File_Loc), "Error: file does not exist at "+str(MCP_File_Loc)
    MCP_data_file = open(MCP_File_Loc,'r+')
    MCPData = pd.read_csv(MCP_data_file)
    MCP_data_file.close()

    # count number of records
    n=MCPData.shape[0]

    # create the whole dataset including Virus ID
    freq = createFreq(MCPData["MCP_Sequence"], 1)
    AAT = []
    for i in range(n):
        AAT.append(MCPData.iloc[i]["Virus_ID"])
        AAT.append(MCPData.iloc[i]["IPC"])
        AAT.append(MCPData.iloc[i]["MCP_len"])
        for j in range(20):
            AAT.append(freq[i][j])
    AAT = np.reshape(np.ravel(AAT), (n, 23))
    x_Phage = np.asarray(AAT)
    # grab the subset of the dataset with just the features
    x_actual = (x_Phage[0:n,1:23]).astype(float)

    # predict T-numbers (ignore the error; the trained forest is imported in line 104)
    y_Pred = rfBest_clf.predict(x_actual)  

    # create an output file
    y_PredTemp = []
    for i in range(n):
        y_PredTemp.append(x_Phage[i,0])
        y_PredTemp.append(tdict2rev[y_Pred[i]])
    y_PredTemp = pd.DataFrame(np.reshape(np.ravel(y_PredTemp),(n,2)))
    y_PredTemp = y_PredTemp.rename(columns={0: 'Virus_ID', 1: 'T_pred'})
    phageMCP2TResult = MCPData.merge(y_PredTemp, how='left', on='Virus_ID')
    phageMCP2TResult.to_csv(r'MCP2TResults.csv', index=False)
    print("Success! see MCP2TResults.csv")
  
    
else:
    print("Sad honk. Invalid input.")
