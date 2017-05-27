import numpy as np
import matplotlib.pyplot as plt
import sys

nargs = len(sys.argv)
if nargs != 3:
    print "Needs 2 arguments, readfile, #Particles"
    sys.exit(1)
args = map(str,sys.argv)
samples = np.loadtxt(args[1])   #read

#samples = np.random.normal(3,0.001,int(1e6))

N = len(samples)
nSteps = int(np.log10(N))


sizes = sorted([(2**i)*(5**j) for j in range(nSteps+1) for i in range(nSteps+1)])[:-2]

varlist = []

for size in sizes:
    nBlocks = N/size
    blocks = np.zeros((nBlocks))
    for i in range(nBlocks):
        blocks[i] = (samples[i*size:(i+1)*size]).mean()
    varlist.append(np.sqrt(blocks.var()/(nBlocks-1)))

varlist = np.array(varlist)
sizes = np.array(sizes)
#    print newblocks.var()

# print blocks.var()             


fig,ax = plt.subplots()

ax.plot(sizes,varlist,'-x')
# ax.errorbar(sizes,varlist,yerr = varlist*(1/(2*(N/sizes))**0.5),fmt='x-') 
# ax.set_xticks(sizes)
# ax.xlabel('xlabel')
plt.xlabel('Blocksize',fontsize = 15)
plt.ylabel('Error', fontsize = 15)
plt.title('Blocking with '+args[2] +' particles and %g points' %(N))

plt.show()
