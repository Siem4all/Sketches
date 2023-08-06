import math
import mmh3
import numpy as np
import os
import pickle

from Morris import MorrisCounter
from CEDAR import CEDARCounter
from PclFileParser import plot_counters
from printf import printf


def remove_existing_files():
    # Removes the specified files if they exist.
    if os.path.exists('../res/pcl_files/Morris_increments.pcl'):
        os.remove('../res/pcl_files/Morris_increments.pcl')
    if os.path.exists('../res/pcl_files/CEDAR_increments.pcl'):
        os.remove('../res/pcl_files/CEDAR_increments.pcl')

    if os.path.exists('../res/Normalized_Morris_error.res'):
        os.remove('../res/Normalized_Morris_error.res')
    if os.path.exists('../res/Normalized_CEDAR_error.res'):
        os.remove('../res/Normalized_CEDAR_error.res')


class CountMinSketch:
    def __init__(self, width, depth, num_flows, mode):
        """
        SketchStatistics: A Count-Min Sketch data structure for estimating flow frequencies.
        Parameters:
        - width (int): the number of counters per row.
        - depth (int): the number of rows or hash functions.
        - num_flows (int): the total number of flows to be estimated.
        - mode : It is either of the counter types that i use.
       Attributes:
        - counters: a CntrMaster object that holds the counter values for the sketch.
        - pair: a tuple of row and column indices used to compute the corresponding counter index.
       Notes:
        - The pair attribute is used to compute the index of each counter. Instead of adding row with column which gives us wrong mapping, i have
        paired the depth number with the hash value and i have inserted them to the pair list to use the index of the pair as a counter index.
         """
        self.mode=mode
        self.width, self.depth, self.num_flows = width, depth, num_flows  # depth is equals to number of hash functions
        if self.mode=='Morris':  
            self.estimate_counter = MorrisCounter(cntrSize=4, numCntrs=depth * width, a=10, cntrMaxVal=1000, verbose=[]) # Initialize the Morris counter
        else:
            self.flow_counter = CEDARCounter(cntrSize=4, numCntrs=depth * width, cntrMaxVal=1000)  # Initialize the CEDAR counter
        self.pair = [f"{i},{j}" for i in range(depth) for j in range(width)]

    def incFlow(self, flow):
        """
        Increment either the Morris or CEDAR counters for the given flow by hashing and updating the appropriate counter values
        """
        for seed in range(self.depth):
            counterIndex = self.pair.index(f'{seed},{mmh3.hash(str(flow), seed) % self.width}')
            if self.mode=='Morris':
                self.estimate_counter.incCntr(cntrIdx=counterIndex, factor=1, verbose=[], mult=False)  
            else:
                self.flow_counter.cntrIncrement(flowIdx=counterIndex, factor=1)  

    def queryFlow(self, flow):
        # Query the minimum Morris or CEDAR counter value for the given flow by hashing and finding the minimum value among the appropriate counters
        minNum = math.inf
        for seed in range(self.depth):
            counterIndex = self.pair.index(f'{seed},{mmh3.hash(str(flow), seed) % self.width}')
            if self.mode=='Morris':
                estimate = self.estimate_counter.queryCntr(counterIndex)
                minNum = min(estimate['val'], minNum)
            else:
                estimate = self.flow_counter.queryCntr(counterIndex)
                minNum = min(estimate, minNum)
        return minNum

    def runSimulationAndCalculateNormalizedRMSE(self):
        """"
        This simulation calculates the Root Mean Square Error (RMSE) and Normalized RMSE for different counters types.
        i have updates the counter numOfIncrements times by generating random flow from 0 to self.num_flows.
        It stores the real value at the index of the flow and it queries the estimated value of the flow from the array counter after
        incrementing it using the methodology of count min sketch. 
        """
        numOfIncrements = 1000
        realValCntr = np.zeros(self.num_flows)
        errorValue = [0]* numOfIncrements 
        # Increment the counters randomly and update the real values
        for i in range(numOfIncrements):
            # Choose a random flow to increment
            flow = np.random.randint(self.num_flows)
            # Increment the corresponding counter in the Count-Min Sketch
            self.incFlow(flow)
            # Update the corresponding real value counter
            realValCntr[flow] += 1
            # Compute the error between the estimated and real frequencies for the flow and update the error_value at that index
            errorValue[i]= abs(realValCntr[flow] - self.queryFlow(flow)) / realValCntr[flow]
        # Compute the Root Mean Square Error (RMSE) and Normalized RMSE over all flows using the error_value list 
        RMSE = math.sqrt(sum(error ** 2 for error in errorValue) / numOfIncrements)
        Normalized_RMSE = RMSE / numOfIncrements
        if self.mode=='Morris':
            # Write the results to a file and return them as a dictionary
            resFile = open(f'../res/Normalized_Morris_error.res', 'a+')
            printf(resFile, '\nnumOfIncrements={},  numCntrs={}, RMSE={}, Normalized_RMSE={}\n'.format(numOfIncrements,
                                                                                                       self.width * self.depth,
                                                                                                       RMSE,
                                                                                                       Normalized_RMSE))
            simulation_results = {
                'width': self.width,
                'Normalized_RMSE': Normalized_RMSE,
            }
            # Append the simulation results to the appropriate file
            with open(f'../res/pcl_files/Morris_increments.pcl', 'ab') as f:
                pickle.dump(simulation_results, f)
        else:
            # Write the results to a file and return them as a dictionary
            resFile = open(f'../res/Normalized_CEDAR_error.res', 'a+')
            printf(resFile, '\nnumOfIncrements={},  numCntrs={}, RMSE={}, Normalized_RMSE={}\n'.format(numOfIncrements,
                                                                                                       self.width * self.depth,
                                                                                                       RMSE,
                                                                                                       Normalized_RMSE))
            simulation_results = {
                'width': self.width,
                'Normalized_RMSE': Normalized_RMSE,
            }
            # Append the simulation results to the appropriate file
            with open(f'../res/pcl_files/CEDAR_increments.pcl', 'ab') as f:
                pickle.dump(simulation_results, f)
        # Overall, this function simulates a Count-Min Sketch and calculates the Normalized RMSE over a set of flows.
        # It provides a way to evaluate the accuracy of the Count-Min Sketch and compare different increments.
        # The Normalized RMSE is a useful metric for measuring the accuracy of the Count-Min Sketch in estimating the frequencies.


def main(num_flows=100):
    remove_existing_files()
    depth = 2  # default value
    counter_types = ['Morris', 'CEDAR']
    for counter_type in counter_types:
        for width in range(2, num_flows//2, 3):
            cmc = CountMinSketch(width=width, depth=depth, num_flows=num_flows, mode=counter_type)
            cmc.runSimulationAndCalculateNormalizedRMSE()
    plot_counters()


if __name__ == '__main__':
    main()
