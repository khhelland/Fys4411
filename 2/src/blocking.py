import numpy as np
import matplotlib.pyplot as plt
    

blocks = np.loadtxt("blocking.dat")   #read

# blocks = np.arange(1e6)
nBlocks = len(blocks)
nSteps = int(np.log10(nBlocks))

# blocks = blocks[:2**nSteps]

# mean = blocks.mean()
# mean2 = mean*mean
step = 2

print nBlocks

sizes = sorted([(2**i)*(5**j) for j in range(nSteps+1) for i in range(nSteps+1)])

varlist = []

for size in sizes:
    newblocks = np.zeros((nBlocks/size))
    for i in range(size):
        newblocks += blocks[i::size]
    newblocks /= size
    varlist.append(newblocks.var())
    print newblocks.var()



# print blocks.var()             

# while(True):
#     try:
#         newblocks = np.zeros((len(blocks)/step))
#         for i in range(step):
#             newblocks += blocks[i::step]
#         blocks = newblocks/step
#         blockvar = blocks.var()
#         varlist.append(blockvar)
#         print blockvar
#     except ValueError:
#         break

plt.figure()
plt.plot(sizes,varlist)
plt.show()
