import numpy as np
import matplotlib.pyplot as plt
import sys

def readData(filename):
    
    infile = open("%s" %filename, 'r')
    energies = []
    
    for line in infile:
        energies.append(float(line))
    
    infile.close()
    return np.asarray(energies)
    
def blocking(energies, nBlocks, blockSize):
    
    meansOfBlocks = np.zeros(nBlocks)
    
    for i in range(nBlocks):
        energiesOfBlock = energies[i*blockSize:(i+1)*blockSize]
        meansOfBlocks[i] = sum(energiesOfBlock)/blockSize
    
    mean = sum(meansOfBlocks)/nBlocks
    mean2 = sum(meansOfBlocks**2)/nBlocks
    variance = mean2 - mean**2
    
    return mean, variance
    
if __name__ == "__main__":
    N = int(sys.argv[1])
    energies = readData(sys.argv[3])
    
    deltaBlockSize = 100
    minBlockSize = 10
    maxBlockSize = int(sys.argv[2])
    numberOfSizes = (maxBlockSize-minBlockSize)/deltaBlockSize + 1
    largestBlockSize = minBlockSize + (numberOfSizes-1)*deltaBlockSize #9910
    
    #blockSizes = np.zeros(numberOfSizes)
    blockSizes = np.linspace(minBlockSize, largestBlockSize, numberOfSizes).astype(int)
    blockAmounts = len(energies)/blockSizes#np.zeros(numberOfSizes)
    means = np.zeros(numberOfSizes)
    variances = np.zeros(numberOfSizes)
    
    for i in range(numberOfSizes):
        #blockSize = minBlockSize + i*deltaBlockSize
        #blockAmount = len(energies)/blockSize
        #mean, variance = blocking(energies, blockAmount, blockSize)
        mean, variance = blocking(energies, blockAmounts[i], blockSizes[i])
        means[i] = mean
        variances[i] = variance
    
    standardDeviation = np.sqrt(abs(variances)/(blockAmounts-1.))

    font = {'family' : 'serif',
            'size'   : 15}

    plt.rc('font',**font)

    fig, ax = plt.subplots()

    ax.plot(blockSizes, standardDeviation,'b')

    ax.set_ylabel(r"Standard Deviation $\sigma$")
    ax.set_xlabel(r'Block size')

    ax.grid('on')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.savefig('/home/alexanfl/uio/masterthesis/doc/part2/figures/blocking-example.eps')
    plt.show()
