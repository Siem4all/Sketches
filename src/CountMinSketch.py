import math
import mmh3
import numpy as np
import os
import pickle

from Morris import MorrisCntr
from CEDAR import CEDARCntr
from RealCounter import RealCntr

from PclFileParser import plot_counters
from printf import printf

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
        self.mode=mode
        self.counter_type=mode
        self.width, self.depth, self.num_flows = width, depth, num_flows
        self.numCntrs=self.width*self.depth
        self.cntrSize=cntrSize
        self.cntrMaxVal=cntrMaxVal
        if self.mode=='Morris':
            self.counter_type = MorrisCntr(cntrSize=self.cntrSize, numCntrs=self.numCntrs, a=None, cntrMaxVal=self.cntrMaxVal, verbose=[]) # Initialize the Morris counter
        elif self.mode=='realCounter':
             self.counter_type=RealCntr(cntrSize=self.cntrSize, numCntrs=self.numCntrs)
        elif self.mode=='CEDAR':
             self.counter_type = CEDARCntr(cntrSize=self.cntrSize, numCntrs=self.numCntrs, cntrMaxVal=self.cntrMaxVal)  # Initialize the CEDAR counter
        else:
            print(f'Sorry, the mode {self.mode} that you requested is not supported')

        self.pair = [(i, j) for i in range(depth) for j in range(width)]
    def queryFlow(self, flow):
        # Query the minimum Morris, CEDAR, real counter value for the given flow by hashing and finding the minimum value among the appropriate counters
        minNum = math.inf
        for seed in range(self.depth):
            counterIndex = self.pair.index((seed,mmh3.hash(str(flow), seed) % self.width))
            minNum = min(self.counter_type.queryCntr(counterIndex), minNum)
        return minNum

    def calculateWrRMSEForDiffModesAndCntrs(self):
        """
        This simulation calculates the write relative errors in order to determine the Root Mean Square Error (RMSE) and Normalized RMSE for different types of counters.
        To calculate the write relative errors for each increment, the real and estimated counter values of each increment need to be determined.
        In the Count Min Sketch, every flow increments multiple counters in the array.
        For example, if the depth (number of hash functions) is 2, a single flow will increment exactly two counters.
        Therefore, it is necessary to calculate the write relative error for each increment of the counters.
        This means that for a single flow, the write relative error needs to be calculated for each of the hash functions.
        """
        numOfIncrements = 1000
        realCntrVal =[0] * self.numCntrs # this holds the real counter value of each counter
        relativeErrors=[]
        # Increment the counters randomly and update the real values
        for incNum in range(numOfIncrements):
            # Choose a random flow to increment
            flow = np.random.randint(self.num_flows)
            for seed in range(self.depth):
                cntrIdx = self.pair.index((seed, mmh3.hash(str(flow), seed) % self.width))
                realCntrVal[cntrIdx]+=1  # this stores the real value of a counter at the specified conter index
                estimatedCntrVal=self.counter_type.incCntr(cntrIdx=cntrIdx, factor=1)
                relativeErrors.append(((realCntrVal[cntrIdx] - estimatedCntrVal)/realCntrVal[cntrIdx])**2)
        # Compute the Root Mean Square Error (RMSE) and Normalized RMSE over all flows using the error_value list
        RMSE = math.sqrt(sum(relativeErrors)/len(relativeErrors))
        Normalized_RMSE =  RMSE/len(relativeErrors)
        # Write the results to a file and return them as a dictionary
        simulation_results = {
            'mode':self.mode,
            'width': self.width,
            'Normalized_RMSE': Normalized_RMSE
         }
        # Append the simulation results to the appropriate file
        with open(f'../res/pcl_files/WrRmse.pcl', 'ab') as f:
            pickle.dump(simulation_results, f)

    def calculateRdRMSEForDiffModesAndCntrs(self):
        """
        This simulation calculates the Read Root Mean Square Error (RMSE) and Normalized RMSE for different counters types.
        It stores the real frequency of the flow at the index of the flow and it queries the estimated value of the flow from the array counter after
        incrementing it using the methodology of count min sketch.
        """
        numOfIncrements = 1000
        realFreq = np.zeros(self.num_flows)
        relativeError = [0]* numOfIncrements
        # Increment the counters randomly and update the real values

        for incNum in range(numOfIncrements):
            # Choose a random flow to increment
            flow = np.random.randint(self.num_flows)
            realFreq[flow] += 1
            for seed in range(self.depth):
                counterIndex = self.pair.index((seed,mmh3.hash(str(flow), seed) % self.width))
                self.counter_type.incCntr(cntrIdx=counterIndex, factor=1)
            relativeError[incNum] =abs(realFreq[flow] - self.queryFlow(flow))/realFreq[flow]
        # Compute the Root Mean Square Error (RMSE) and Normalized RMSE over all flows using the error_value list
        RMSE = math.sqrt(sum([relativeErrors[incNum] for incNum in range(numOfIncrements)])/numOfIncrements)
        Normalized_RMSE =  RMSE/numOfIncrements
        # Write the results to a file and return them as a dictionary
        simulation_results = {
            'mode':self.mode,
            'width': self.width,
            'Normalized_RMSE': Normalized_RMSE,
        }
        # Append the simulation results to the appropriate file
        with open(f'../res/pcl_files/RdRmse.pcl', 'ab') as f:
            pickle.dump(simulation_results, f)

def main():
    """
    It iterate the first for loop based on the list of counter types and tranfer the counter type as a mode to the CountMinSketchWithCU class
    to calculate Normalized_RMSE using Morris or real counter.
    """
    num_flows=20
    depth = 2  # default value
    counter_types = ['Morris','realCounter']
    for counter_type in counter_types:
        for width in range(2, num_flows//2, 2):
            cmc = CountMinSketch(width=width, depth=depth, num_flows=num_flows, mode=counter_type, cntrSize=6, cntrMaxVal=1000)
            cmc.calculateWrRMSEForDiffModesAndCntrs()
    plot_counters()


if __name__ == '__main__':
    main()