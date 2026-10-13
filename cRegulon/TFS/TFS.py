#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 13 17:29:03 2023

@author: zaine
"""
import numpy as np

f = open('../CSI/Clean2/TFName.txt')
TF = f.readlines();f.close()
TF = [TF[i].strip('\n') for i in range(len(TF))]

eps = 1e-4
def FC(e1,e2):
    if e2 > 0:
        fc = e1/e2-0.5
    else:
        fc = e1/eps-0.5
    if fc > 1:
        fc = 1
    if fc < 0:
        fc = 0
    return fc

f = open('../CSI/SampleNames.txt')
ct = f.readlines();f.close()
ct = [ct[i].strip('\n') for i in range(len(ct))]

E = []
for c in range(len(ct)):
    f = open('../../PECA2/Input/'+ct[c]+'.txt')
    a = f.readlines();f.close()
    exp = [[],[]]
    for i in range(len(a)):
        a[i] = a[i].strip('\n').split('\t')
        exp[0].append(a[i][0])
        exp[1].append(a[i][1])
    tmp = []
    for i in range(len(TF)):
        indel = exp[0].index(TF[i])
        tmp.append(exp[1][indel])
    E.append(tmp)

E = np.array(E).astype('float').T
g = open('TFES.txt','w')
g.write('TF\t'+'\t'.join(ct)+'\n')
for i in range(len(TF)):
    g.write(TF[i])
    for j in range(len(ct)):
        xx = E[i][j];yy = np.mean(np.delete(E[i],j))
        fc = FC(xx,yy)
        g.write('\t'+str(fc))
    g.write('\n')
