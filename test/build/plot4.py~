import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
import sys
fig, ax = plt.subplots(1,1, figsize=(6,4.5))
fig.suptitle("Smooth coefficients", fontsize=15)

f = open("solutionOK40.txt", "r")
l = f.readlines()
x=[]
p=[]
s=[]
for line in l:
    data = line.split(" ")
    x.append(float(data[0]))
    p.append(float(data[1]))
ax.plot(x,p,lw=2, label="s = 1.8",ls="-")
ax.set_ylabel("Value", fontsize=12)
ax.set_xlabel("$x$-axis", fontsize=12)

f2 = open("solutionOK25.txt", "r")
l2 = f2.readlines()
x=[]
p=[]
s=[]
for line in l2:
    data = line.split(" ")
    x.append(float(data[0]))
    p.append(float(data[1]))
ax.plot(x,p,lw=2, label="s = 1.5",ls="--")


f3 = open("solutionOK10.txt", "r")
l3 = f3.readlines()
x=[]
p=[]
s=[]
for line in l3:
    data = line.split(" ")
    x.append(float(data[0]))
    p.append(float(data[1]))
ax.plot(x,p,lw=2, label = "s = 1.2", ls = "-.")

ax.legend()
plt.savefig("/Users/raylwin/Dropbox/fractional/solExp3.eps", format="eps")
#plt.show()
 
