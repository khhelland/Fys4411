import numpy as np
import matplotlib.pyplot as plt


blocks = np.random.randn(2**20)   #read

nBlocks = len(blocks)
nSteps = int(np.log2(len(blocks)))
# print nSteps
blocks = blocks[:2**nSteps]

mean = blocks.mean()
mean2 = mean*mean

varlist = []

for s in xrange(nSteps):
    blocks = 0.5*(blocks[::2] + blocks[1::2])
    blockvar = (blocks*blocks).mean() - mean2
    varlist.append(blockvar)


plt.figure()
plt.plot(xrange(nSteps),varlist)
plt.show()
    
