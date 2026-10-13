#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue 17th Oct 10:05:38 2023

@author: Zhanying
"""

import numpy as np
from sklearn.decomposition import NMF
from numpy import linalg as LA
import os
import warnings
from sklearn.exceptions import ConvergenceWarning
warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
import argparse

parser = argparse.ArgumentParser(description="cRegulon")
parser.add_argument('--k', type=int, required=True, help='models')
args = parser.parse_args()

f = open('../CSI/SampleNames.txt')
ct = f.readlines();f.close()
ct = [ct[i].strip('\n') for i in range(len(ct))]
I = len(ct)

C = []
for c in range(I):
    C1 = np.loadtxt("../CSIS/"+ct[c]+'_CSIS.txt')
    C.append(C1);

f = open('../CSI/Clean2/TFName.txt')
TF = f.readlines();f.close()
TF = [TF[i].strip('\n') for i in range(len(TF))]

def Norm(X):
    X = X/np.sqrt(np.sum(X**2,0).reshape(1,X.shape[1]))
    X = X/(np.sum(X,1).reshape(X.shape[0],1)+eps)
    return X
def NormL(L):
    L = L/np.sum(L)
    return L
eps = 0.001;
def Loss(C,x0,l0,mu1):
    t1 = 0;t2 = 0
    for c in range(I):
        t1 += pow(LA.norm(C[c] - x0.dot(l0[c]).dot(x0.T), ord='fro'), 2)
        t2 += mu1*np.sum(C[c] * x0.dot(l0[c]).dot(x0.T))
    return t1 - t2
def RunTFModule(CC,II,KK,Out):
    C = CC.copy();I = II;K = KK;OutFold = Out
    C0 = np.zeros(C[0].shape)
    for c in range(I):
        C0 += C[c]
    C0 = C0/I
    err1 = np.zeros(50)
    for i in range(50):
        model = NMF(n_components=K, init='random', random_state=i, solver='cd', max_iter=50)
        X0 = model.fit_transform(C0);X1 = model.components_
        err1[i] = LA.norm(C0 - np.dot(X0, X1), ord='fro')
    model = NMF(n_components=K, init='random', random_state=np.argmin(err1), solver='cd', max_iter=1000)
    X0 = model.fit_transform(C0);X1 = model.components_
    np.savetxt(OutFold+'X0.txt',X0,delimiter='\t')
    np.savetxt(OutFold+'X1.txt',X1,delimiter='\t')
    L0 = []
    for c in range(I):
        Li = np.diag([1/K for i in range(K)])
        L0.append(Li)
    L1 = 0;L2 = 0
    for c in range(I):
        L1 += pow(LA.norm(C[c] - X0.dot(L0[c]).dot(X0.T), ord='fro'), 2)
        L2 += np.sum(C[c] * (X0.dot(L0[c]).dot(X0.T)))
    mu1 = L2/L1

    err = 1000;
    X = X0.copy();L = L0.copy();loss1 = Loss(C,X,L,mu1);it = 0
    while err > 1e-4:
        XTX = X.T.dot(X)
        XX1 = np.zeros(X.shape);XX2 = np.zeros(X.shape)
        for c in range(I):
            XX1 += C[c].dot(X).dot(L[c])
            XX2 += X.dot(L[c]).dot(XTX).dot(L[c])
        XNext = X * ((1+0.5*mu1)*XX1)/(eps+XX2)
        XNext = Norm(XNext)
        LNext = []
        for c in range(I):
            LNext.append(L[c]*((1+mu1)*X.T.dot(C[c]).dot(X))/(eps+XTX.dot(L[c]).dot(XTX)))
            LNext[c] = NormL(LNext[c])
                
        loss2 = Loss(C,XNext,LNext,mu1)
        err = np.abs(loss2-loss1)
        loss1 = loss2;it += 1
        X = XNext.copy();L = LNext.copy()
    print("The final loss is:",loss1,it)

    g = open(OutFold+'X.txt','w')
    g.write('TF\t'+'\t'.join(["C"+str(c+1) for c in range(K)])+'\n')
    for i in range(len(TF)):
        g.write(TF[i])
        for j in range(K):
            g.write('\t'+str(X[i][j]))
        g.write('\n')
    g.close()

    for i in range(X.shape[1]):
        D1 = []
        for j in range(X.shape[0]):
            if X[j,i] > 0.05:
                D1.append(X[j,i])
        me = np.mean(D1);sd = np.std(D1)
        for j in range(X.shape[0]):
            X[j,i] = (X[j,i]-me)/sd

    g = open(OutFold+'XN.txt','w')
    g.write('TF\t'+'\t'.join(["C"+str(c+1) for c in range(K)])+'\n')
    for i in range(len(TF)):
        g.write(TF[i])
        for j in range(K):
            g.write('\t'+str(X[i][j]))
        g.write('\n')
    g.close()
    LL = []
    for c in range(I):
        LL.append(np.diag(L[c]))
    LL = np.array(LL)
    g = open(OutFold+'A.txt','w')
    g.write('CellType\t'+'\t'.join(["C"+str(c+1) for c in range(K)])+'\n')
    for c in range(I):
        g.write(ct[c].strip('\n'))
        for j in range(K):
            g.write('\t'+str(LL[c][j]))
        g.write('\n')
    g.close()

K=args.k
OutFold = "./TFModule"+str(K)+"/";os.makedirs(OutFold,exist_ok=True)
print("Finding",K,"TF modules in fetal mouse")
RunTFModule(C,I,K,OutFold)
