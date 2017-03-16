import matplotlib.pyplot as plt
import numpy as np
from math import sqrt

E = open("Energies.dat","r")
times = [[],[],[],[]]
next(E)
for line in E:
    cols = line.split()
    vals = cols[0].split(",")
    vals +=cols[1:]
    n = int(vals[0])
    if n == 2:
        times[0].append((float(vals[4])))        
    elif n==6:
        times[1].append((float(vals[4])))
    elif n== 12:
        times[2].append((float(vals[4])))
    elif n == 20:
        times[3].append((float(vals[4])))


plt.figure()
plt.hold(1)
for n in range(1,5):
    particles = n*(n+1)
    
    plt.plot(range(n,len(times[n-1])+n),times[n-1],"-x",label="%i particles" %(particles))
    

plt.title("Time used for $\omega = 1$",fontsize=20)
plt.semilogy()
plt.legend(loc = "lower right")
plt.xlabel("shells",fontsize = 14)
plt.ylabel("time [s]",fontsize=14)
#fig.tight_layout()
plt.show()
