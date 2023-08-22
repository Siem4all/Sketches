import math
import mmh3
import numpy as np
import os
import pickle
import PclFileParser

from Morris import MorrisCntr
from CEDAR import CEDARCntr
from RealCounter import RealCntr
from printf import printf

def remove_existing_files():
    # Removes the specified files if they exist.
    if os.path.exists('../res/pcl_files/RdRmse.pcl'):
        os.remove('../res/pcl_files/RdRmse.pcl')
    if os.path.exists('../res/RdRMSE.res'):
        os.remove('../res/RdRMSE.res')

class CountMinSketch:
    def __init__(self, width, depth, num_flows, mode, cntrSize, cntrMaxVal):
        """
        Parameters:
        - width (int): the number of counters per row.
        - depth (int): the number of rows or hash functions.
        - num_flows (int): the total number of flows to be estimated.
        - mode : It is the counter types that i use.
       Attributes:
        - counter_type: a counter object that holds the counter values for the sketch.
        - pair: a tuple of row and column indices used to compute the corresponding counter index.
       Notes:
        - The pair attribute is used to compute the index of each counter. Instead of adding row with column which gives us wrong mapping, i have
        paired the depth number with the hash value and i have inserted them to the pair list to use the index of the pair as a counter index.
         """
        self.mode, self.counter_type, self.cntrSize, self.cntrMaxVal=mode, mode, cntrSize, cntrMaxVal
        self.width, self.depth, self.num_flows = width, depth, num_flows
        self.numCntrs=self.width*self.depth
        if self.mode=='Morris':
            self.counter_type = MorrisCntr(cntrSize=self.cntrSize, numCntrs=self.numCntrs, a=None, cntrMaxVal=self.cntrMaxVal, verbose=[]) # Initialize the Morris counter
        elif self.mode=='realCounter':
             self.counter_type=RealCntr(cntrSize=self.cntrSize, numCntrs=self.numCntrs)
        elif self.mode=='CEDAR':
             self.counter_type = CEDARCntr(cntrSize=self.cntrSize, numCntrs=self.numCntrs, cntrMaxVal=self.cntrMaxVal)  # Initialize the CEDAR counter
        else:
            print(f'Sorry, the mode {self.mode} that you requested is not supported')

        self.pair = [(i, j) for i in range(depth) for j in range(width)]

    def incNQueryFlow(self, flow):
        """
        When a flow arrives, it is hashed using the hash functions, and the corresponding counters are incremented.
        At the end,  the minimum value of the corresponding counters is turned as the estimate.
        """
        cntrValAfterInc = [0]*self.depth
        for mappedCntr in range(self.depth):
            counterIndex               = self.pair.index((mappedCntr, mmh3.hash(str(flow), seed=mappedCntr) % self.width))
            cntrValAfterInc[mappedCntr]=self.counter_type.incCntr(cntrIdx=counterIndex, factor=1)
        return min(cntrValAfterInc)
    def incFlow(self, flow):
        # increment the mapped counters values
        for seed in range(self.depth):
                counterIndex = self.pair.index((seed,mmh3.hash(str(flow), seed) % self.width))
                self.counter_type.incCntr(cntrIdx=counterIndex, factor=1)
    def queryFlow(self, flow):
        # Query the minimum Morris, CEDAR, real counter value for the given flow by hashing and finding the minimum value among the appropriate counters
        minNum = math.inf
        for seed in range(self.depth):
            counterIndex = self.pair.index((seed,mmh3.hash(str(flow), seed) % self.width))
            minNum = min(self.counter_type.queryCntr(counterIndex), minNum)
        return minNum
    def calculateNormalizedRMSE(self):
        """
        This simulation calculates the Read Root Mean Square Error (RMSE) and Normalized RMSE for different counters types.
        It stores the real frequency of a flow at the index of a flow and it queries the estimated value of the flow from the array counter after
        incrementing it using the methodology of count min sketch.
        """
        numOfIncrements = 1000
        realCntr        = np.zeros(self.num_flows)
        sumOfAllErors   = 0
        # Increment the counters randomly and update the real values
        for _ in range(numOfIncrements): #$$$ Great. As you can realize, the code now is much simpler, faster and more readable!
            # Choose a random flow to increment
            flow = np.random.randint(self.num_flows)
            realCntr[flow] += 1  # increment the real counter value of the flow upon its arrival
            sumOfAllErors  +=((realCntr[flow] - self.incNQueryFlow(flow))/realCntr[flow])**2
        # Compute the Root Mean Square Error (RMSE) and Normalized RMSE over all flows using the sumOfAllErors
        RMSE            = math.sqrt(sumOfAllErors/numOfIncrements)
        Normalized_RMSE =  RMSE/numOfIncrements
        resFile = open (f'../res/RdRMSE.res', 'a+')
        printf(resFile, '\nmode={}, numOfIncrements={}, numCntrs={}, sumOfAllErors={}, RMSE={}, Normalized_RMSE={}\n'.format(self.mode,
               numOfIncrements, self.numCntrs, sumOfAllErors, RMSE, Normalized_RMSE))
        # Write the results to a file and return them as a dictionary
        simulation_results = {
            'mode':self.mode,
            'width': self.width,
            'Normalized_RMSE': Normalized_RMSE,
        }
        # Append the simulation results to the appropriate file
        with open(f'../res/pcl_files/RdRmse.pcl', 'ab') as f:
            pickle.dump(simulation_results, f)

def main(counter_types):
    """
    It iterate the first for loop based on the list of counter types and tranfer the counter type as a mode to the CountMinSketch class
    to calculate Normalized_RMSE using Morris, real counter or CEDAR architecture.
    """
    remove_existing_files()
    num_flows     =20
    depth         = 2
    for counter_type in counter_types:
        for width in range(2, num_flows//2, 2):
            cmc = CountMinSketch(width=width, depth=depth, num_flows=num_flows, mode=counter_type, cntrSize=6, cntrMaxVal=1000)
            cmc.calculateNormalizedRMSE()

if __name__ == '__main__':
    counter_types = ['CEDAR', 'Morris','realCounter']
    main(counter_types)
    parser = PclFileParser.PclFileParser()
    parser.rdPcl()
    parser.NRMSEVsWidthPlot(modes=counter_types)
