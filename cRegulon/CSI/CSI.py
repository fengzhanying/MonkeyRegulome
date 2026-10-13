#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 13 15:34:48 2023

@author: Zhanying
"""

import numpy as np

f = open("SampleNames.txt")
ct = f.readlines();f.close()
ct = [ct[i].strip('\n') for i in range(len(ct))]

f = open(f'../../PECA2/Results/{ct[0]}/TGName.txt')
TG = f.readlines();f.close()
g = open('./Clean0/TGName.txt','w')
[g.write(TG[i]) for i in range(len(TG))]
g.close()
f = open(f'../../PECA2/Results/{ct[0]}/TFName.txt')
TF = f.readlines();f.close()
g = open('./Clean0/TFName.txt','w')
[g.write(TF[i]) for i in range(len(TF))]
g.close()

for c in range(len(ct)):
    TRS = np.loadtxt(f'../../PECA2/Results/{ct[c]}/TFTG_regulationScore.txt')
    TRS = np.log2(TRS+1)
    
    TRS1 = TRS.copy()
    for i in range(TRS1.shape[0]):
        if np.max(TRS1[i]) > 0:
            TRS1[i] = (TRS1[i]-np.mean(TRS1[i]))/np.std(TRS1[i])
    TRS2 = TRS.copy();TRS2 = TRS2.T
    for i in range(TRS2.shape[0]):
        if np.max(TRS2[i]) > 0:
            TRS2[i] = (TRS2[i]-np.mean(TRS2[i]))/np.std(TRS2[i])
    TRS2 = TRS2.T
    TRS3 = (TRS1+TRS2)/2;TRS3[TRS3<0] = 0
    g = open(f'./Clean0/{ct[c]}_NTRS.txt','w')
    for i in range(TRS3.shape[0]):
        s = []
        for j in range(TRS3.shape[1]):
            s.append(str(TRS3[i][j]))
        g.write('\t'.join(s)+'\n')
    g.close()
    
    TRS = np.loadtxt(f'./Clean0/{ct[c]}_NTRS.txt')
    PCC = np.corrcoef(TRS);PCC[np.isnan(PCC)] = 0.0
    
    CSI = np.zeros((TRS.shape[0],TRS.shape[0]))
    for i in range(TRS.shape[0]):
        for j in range(TRS.shape[0]):
            n = 0;cutoff = PCC[i][j] - 0.05
            for k in range(TRS.shape[0]):
                if PCC[i][k] < cutoff and PCC[j][k] < cutoff:
                    n += 1
            CSI[i][j] = n/TRS.shape[0]
    
    g = open(f'./Clean0/{ct[c]}_CSI.txt','w')
    for i in range(CSI.shape[0]):
        st = ''
        for j in range(CSI.shape[1]):
            st += (str(CSI[i][j])+'\t')
        g.write(st.strip('\t')+'\n')
    g.close()


f = open('./Clean0/TFName.txt')
TF = f.readlines();f.close()
TF = [TF[i].strip('\n') for i in range(len(TF))]
f = open('./Clean0/TGName.txt')
TG = f.readlines();f.close()
TG = [TG[i].strip('\n') for i in range(len(TG))]

TFSTD = [];TGSTD = []
for c in range(len(ct)):
    TRS = np.loadtxt(f'../../PECA2/Results/{ct[c]}/TFTG_regulationScore.txt')
    TFSTD.append(np.std(TRS,1))
    TGSTD.append(np.std(TRS,0))
TFSTD = np.array(TFSTD);TFSTD = np.max(TFSTD,0)
g = open('./Clean1/TFName.txt','w')
IndelTF = []
for i in range(len(TF)):
    if TFSTD[i] > 0:
        g.write(TF[i]+'\n')
        IndelTF.append(i)
g.close()
TGSTD = np.array(TGSTD);TGSTD = np.max(TGSTD,0)
g = open('./Clean1/TGName.txt','w')
IndelTG = []
for i in range(len(TG)):
    if TGSTD[i] > 0:
        g.write(TG[i]+'\n')
        IndelTG.append(i)
g.close()

for c in range(len(ct)):
    TRS = np.loadtxt(f'../../PECA2/Results/{ct[c]}/TFTG_regulationScore.txt')
    TRS = TRS[IndelTF,][:,IndelTG];TRS = np.log2(TRS+1)
    TRS1 = TRS.copy()
    for i in range(TRS1.shape[0]):
        if np.max(TRS1[i]) > 0:
            TRS1[i] = (TRS1[i]-np.mean(TRS1[i]))/np.std(TRS1[i])
    TRS2 = TRS.copy();TRS2 = TRS2.T
    for i in range(TRS2.shape[0]):
        if np.max(TRS2[i]) > 0:
            TRS2[i] = (TRS2[i]-np.mean(TRS2[i]))/np.std(TRS2[i])
    TRS2 = TRS2.T
    TRS3 = (TRS1+TRS2)/2;TRS3[TRS3<0] = 0
    g = open(f'./Clean1/{ct[c]}_NTRS.txt','w')
    for i in range(TRS3.shape[0]):
        s = []
        for j in range(TRS3.shape[1]):
            s.append(str(TRS3[i][j]))
        g.write('\t'.join(s)+'\n')
    g.close()
    
    TRS = np.loadtxt(f'./Clean1/{ct[c]}_NTRS.txt')
    
    PCC = np.corrcoef(TRS);PCC[np.isnan(PCC)] = 0.0
    
    CSI = np.zeros((TRS.shape[0],TRS.shape[0]))
    for i in range(TRS.shape[0]):
        for j in range(TRS.shape[0]):
            n = 0;cutoff = PCC[i][j] - 0.05
            for k in range(TRS.shape[0]):
                if PCC[i][k] < cutoff and PCC[j][k] < cutoff:
                    n += 1
            CSI[i][j] = n/TRS.shape[0]

    g = open(f'./Clean1/{ct[c]}_CSI.txt','w')
    for i in range(CSI.shape[0]):
        st = ''
        for j in range(CSI.shape[1]):
            st += (str(CSI[i][j])+'\t')
        g.write(st.strip('\t')+'\n')
    g.close()

f = open('./Clean1/TFName.txt')
TF = f.readlines();f.close()
TF = [TF[i].strip('\n') for i in range(len(TF))]
STD = []
for c in range(len(ct)):
    CSI = np.loadtxt(f'./Clean1/{ct[c]}_CSI.txt')
    STD.append(np.std(CSI,1))
STD = np.array(STD);STD = np.max(STD,0)
g = open('./Clean2/TFName.txt','w')
for i in range(len(TF)):
    if STD[i] > 0:
        g.write(TF[i]+'\n')
g.close()
f = open('./Clean1/TGName.txt')
a = f.readlines();f.close()
g = open('./Clean2/TGName.txt','w')
for i in range(len(a)):
    g.write(a[i])
g.close()

f = open('./Clean0/TFName.txt')
TF = f.readlines();f.close()
TF = [TF[i].strip('\n') for i in range(len(TF))]
f = open('./Clean0/TGName.txt')
TG = f.readlines();f.close()
TG = [TG[i].strip('\n') for i in range(len(TG))]

f = open('./Clean2/TFName.txt')
TF1 = f.readlines();f.close()
IndelTF = []
for i in range(len(TF1)):
    IndelTF.append(TF.index(TF1[i].strip('\n')))
f = open('./Clean2/TGName.txt')
TG1 = f.readlines();f.close()
IndelTG = []
for i in range(len(TG1)):
    IndelTG.append(TG.index(TG1[i].strip('\n')))
    
for c in range(len(ct)):
    TRS = np.loadtxt(f'../../PECA2/Results/{ct[c]}/TFTG_regulationScore.txt')
    TRS = TRS[IndelTF,][:,IndelTG];TRS = np.log2(TRS+1)
    TRS1 = TRS.copy()
    for i in range(TRS1.shape[0]):
        if np.max(TRS1[i]) > 0:
            TRS1[i] = (TRS1[i]-np.mean(TRS1[i]))/np.std(TRS1[i])
    TRS2 = TRS.copy();TRS2 = TRS2.T
    for i in range(TRS2.shape[0]):
        if np.max(TRS2[i]) > 0:
            TRS2[i] = (TRS2[i]-np.mean(TRS2[i]))/np.std(TRS2[i])
    TRS2 = TRS2.T
    TRS3 = (TRS1+TRS2)/2;TRS3[TRS3<0] = 0
    g = open(f'./Clean2/{ct[c]}_NTRS.txt','w')
    for i in range(TRS3.shape[0]):
        s = []
        for j in range(TRS3.shape[1]):
            s.append(str(TRS3[i][j]))
        g.write('\t'.join(s)+'\n')
    g.close()
    
    TRS = np.loadtxt(f'./Clean2/{ct[c]}_NTRS.txt')
    
    PCC = np.corrcoef(TRS);PCC[np.isnan(PCC)] = 0.0
    
    CSI = np.zeros((TRS.shape[0],TRS.shape[0]))
    for i in range(TRS.shape[0]):
        for j in range(TRS.shape[0]):
            n = 0;cutoff = PCC[i][j] - 0.05
            for k in range(TRS.shape[0]):
                if PCC[i][k] < cutoff and PCC[j][k] < cutoff:
                    n += 1
            CSI[i][j] = n/TRS.shape[0]

    g = open(f'./Clean2/{ct[c]}_CSI.txt','w')
    for i in range(CSI.shape[0]):
        st = ''
        for j in range(CSI.shape[1]):
            st += (str(CSI[i][j])+'\t')
        g.write(st.strip('\t')+'\n')
    g.close()
