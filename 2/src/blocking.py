import numpy as np
import matplotlib.pyplot as plt

# samples = np.loadtxt("blocking.dat")   #read

samples = np.random.normal(3,0.001,int(1e6))

nBlocks = len(blocks)
nSteps = int(np.log10(nBlocks))

# blocks = blocks[:2**nSteps]

# mean = blocks.mean()
# mean2 = mean*mean
step = 2

print nBlocks, blocks.mean(), (blocks*blocks).mean(), np.sqrt(blocks.var())

sizes = sorted([(2**i)*(5**j) for j in range(nSteps) for i in range(nSteps)])
print sizes
varlist = []

for size in sizes:
    newblocks = np.zeros((nBlocks/size))
    newblocks = blocks((
    newblocks /= size
    varlist.append(newblocks.var())
#    print newblocks.var()

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
