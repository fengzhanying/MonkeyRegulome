import numpy as np

f = open('../CSI/SampleNames.txt')
ct = f.readlines();f.close()
ct = [ct[i].strip('\n') for i in range(len(ct))]

f = open('../CSI/Clean2/TFName.txt')
TF = f.readlines();f.close()
TF = [TF[i].strip('\n') for i in range(len(TF))]

f = open('../TFS/TFES.txt')
GS = f.readlines();f.close()
ct1 = GS[0];del GS[0]
ct1 = ct1.strip('\n').split('\t');del ct1[0]
gene = []
for i in range(len(GS)):
    GS[i] = GS[i].strip('\n').split('\t')
    gene.append(GS[i][0]);del GS[i][0]
    for j in range(len(ct1)):
        GS[i][j] = float(GS[i][j])
Indel = []
for i in range(len(TF)):
    Indel.append(gene.index(TF[i]))
TFES = []
for i in range(len(TF)):
    TFES.append(GS[Indel[i]])
TFES = np.array(TFES).T
for c in range(len(ct)):
    indel = ct1.index(ct[c]);
    TFES1 = np.dot(TFES[indel].reshape((TFES[indel].shape[0],1)),TFES[indel].reshape((1,TFES[indel].shape[0])))
    CSI = np.loadtxt('../CSI/Clean2/'+ct[c]+'_CSI.txt')
    
    CSIS = (CSI*TFES1)**1/3
    g = open(ct[c]+'_CSIS.txt','w')
    for i in range(CSIS.shape[0]):
        st = ''
        for j in range(CSIS.shape[1]):
            st += (str(CSIS[i][j])+'\t')
        g.write(st.strip('\t')+'\n')
    g.close()
