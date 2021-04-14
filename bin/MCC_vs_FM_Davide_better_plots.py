#!/usr/bin/env python3.6
#

# coding: utf-8
# # Compare MCC to Fowlkes-Mallows index
# import packages
# In[ ]:

import numpy as np
# import matplotlib.pyplot as plt

import math

import matplotlib
from matplotlib import pyplot as plt

from random import seed
from random import randint
# seed(11)

# Generate confusion matrices by varying TPR, TNR and Fowlkes-Mallows index
# In[ ]:
stepsize = 0.01 # 0.01
vals = np.arange(0, 1 + stepsize, stepsize)
N = len(vals) - 1
M = len(vals)

print("Number of values considered", N, "\n")

# In[ ]:
TP = []
TN = []
FN = []
FP = []
count = 1

for prev in vals:
    for TPR in vals:
        for TNR in vals:
            TP.append(prev*TPR)
            TN.append((1-prev) * TNR)
            FN.append(prev*(1-TPR))
            FP.append((1-prev)*(1-TNR))
            count = count + 1
            
TP = np.array(TP)
TN = np.array(TN)
FN = np.array(FN)
FP = np.array(FP)


# Calculate metrics
# In[ ]:
MCC_upper = TP * TN - FP * FN
MCC_lower = np.sqrt((TP +FP) * (TP + FN) * (TN + FP) * (TN +FN))
# MCC = MCC_upper / MCC_lower

MCC=np.zeros_like(MCC_upper)
for i in range(len(MCC_upper)):
    if(MCC_lower[i]!=0):
        MCC[i]=MCC_upper[i]/MCC_lower[i]


normMCC = (MCC + 1) / 2

# prevalence threshold
TPR = TP / (TP + FN)
TNR = TN / (TN + FP)
PPV = TP / (TP + FP)
# PT = (np.sqrt(TPR * (-TNR +1)) + TNR - 1)  / (TPR + TNR - 1)

FMI = np.sqrt(TPR * PPV)

#PT_upper = (np.sqrt(TPR * (-TNR +1)) + TNR - 1)  
#PT_lower = (TPR + TNR - 1)
#PT=np.zeros_like(PT_upper)

#for i in range(len(PT_upper)):
    #if(PT_lower[i]!=0):
        #PT[i]=PT_upper[i]/PT_lower[i]


#complPT = (1 - PT)
 
# plotting

#lengthComplPT = len(complPT)
lengthNormMCC = len(normMCC)

lengthFMI = len(FMI)

print("lengthFMI = ", lengthFMI, "\n")
print("lengthNormMCC = ", lengthNormMCC, "\n")


mySize = 0.2
myColor = "black"

# In[ ]:

matplotlib.rc('axes', edgecolor='white')


plt.rcParams['axes.axisbelow'] = True

plt.rcParams['axes.facecolor'] = 'whitesmoke'
plt.grid(color='white')
plt.scatter(normMCC, FMI, s=mySize, c=myColor)
plt.xlabel('normalized Matthews correlation coefficient (normMCC)')
plt.ylabel('Fowlkes-Mallows index')



# plt.show()

value = randint(1, 1000)

fileName = '../results/normMCC_vs_FMI_python_rand'+str(value)+'.png' 

print("saved file: ", fileName, "\n")
plt.savefig(fileName)

