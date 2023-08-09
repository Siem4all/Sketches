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

def remove_existing_files():
    # Removes the specified files if they exist.
    if os.path.exists('../res/pcl_files/Morris_increments.pcl'):
        os.remove('../res/pcl_files/Morris_increments.pcl')
    if os.path.exists('../res/pcl_files/CEDAR_increments.pcl'):
        os.remove('../res/pcl_files/CEDAR_increments.pcl')

    if os.path.exists('../res/pcl_files/realCounter_increments.pcl'):
        os.remove('../res/pcl_files/realCounter_increments.pcl')


class CountMinSketch:
    def __init__(self, width, depth, num_flows, mode):
        """
        CountMinSketch: A Count-Min Sketch data structure for estimating flow frequencies.
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
        self.width, self.depth, self.num_flows = width, depth, num_flows  # depth is equals to number of hash functions
        self.numCntrs=self.width*self.depth
        if self.mode=='Morris':
            self.counter_type = MorrisCntr(cntrSize=6, numCntrs=self.numCntrs, a=None, cntrMaxVal=1000, verbose=[]) # Initialize the Morris counter
        elif self.mode=='realCounter':
             self.counter_type=RealCntr(cntrSize=6, numCntrs=self.numCntrs)
        elif self.mode=='CEDAR':
             self.counter_type = CEDARCntr(cntrSize=6, numCntrs=self.numCntrs, cntrMaxVal=1000)  # Initialize the CEDAR counter
        else:
            print(f'Sorry, the mode {self.mode} that you requested is not supported')

        self.pair = [f"{i},{j}" for i in range(depth) for j in range(width)]
    def queryFlow(self, flow):
        # Query the minimum Morris, CEDAR, real counter value for the given flow by hashing and finding the minimum value among the appropriate counters
        minNum = math.inf
        for seed in range(self.depth):
            counterIndex = self.pair.index(f'{seed},{mmh3.hash(str(flow), seed) % self.width}')
            if self.mode in ['Morris', 'CEDAR', 'realCounter']:
                estimate = self.counter_type.queryCntr(counterIndex)
                minNum = min(estimate['val'] if (self.mode=='Morris') else estimate, minNum)
            else:
                print(f'Sorry, the mode {self.mode} that you requested is not supported')
        return minNum

    def calculateWriteRMSEForDifferentCounters(self):
        """"
        This simulation calculates write relative errors to find the Root Mean Square Error (RMSE) and Normalized RMSE for different counters types.
        It stores the real value at the index of the flow plus row number or depth pf the array counter and it calculates relative in error for every increment
        incrementing it using the methodology of count min sketch.
        """
        numOfIncrements = 1000
        realValCntr = np.zeros(self.num_flows*self.depth)  # Count the real counter for each increment of the flow in the array counter
        writeErrorVar = []  # At each increment each flow has their own values at each array counters and this holds calculated write error variance of the flow in each array counter
        cntrValue=[]  # it hold the value that the counter can represent
        relativeErrors=[]
        # Increment the counters randomly and update the real values
        for i in range(numOfIncrements):
            # Choose a random flow to increment
            flow = np.random.randint(self.num_flows)
            for seed in range(self.depth):
                realCntrIdx=sum([flow+seed])
                realValCntr[realCntrIdx]+=1
                counterIndex = self.pair.index(f'{seed},{mmh3.hash(str(flow), seed) % self.width}')
                if self.mode in ['Morris', 'CEDAR', 'realCounter']:
                    cntrAfterInc=self.counter_type.incCntr(cntrIdx=counterIndex, factor=1)
                    print(cntrAfterInc, self.mode, counterIndex, i)
                    writeErrorVar.append((realValCntr[realCntrIdx] - cntrAfterInc) ** 2)
                    cntrValue.append(cntrAfterInc)
                else:
                    print(f'Sorry, the mode {self.mode} that you requested is not supported')
        # to calculate write relative error we have to calculate standard variance using the formula stdVar=sqrt(sum((Xbar-Xi)^2)/N)=sqrt(var)
        relativeError = [math.sqrt(writeErrorVar[i]) / cntrValue[i] for i in range(len(cntrValue))]
        # Compute the Root Mean Square Error (RMSE) and Normalized RMSE over all flows using the error_value list
        RMSE = math.sqrt(sum(error ** 2 for error in relativeError) / numOfIncrements)
        Normalized_RMSE = RMSE / numOfIncrements
        if self.mode=='Morris':
            # Write the results to a file and return them as a dictionary
            simulation_results = {
                'width': self.width,
                'Normalized_RMSE': Normalized_RMSE,
            }
            # Append the simulation results to the appropriate file
            with open(f'../res/pcl_files/Morris_increments.pcl', 'ab') as f:
                pickle.dump(simulation_results, f)
        elif self.mode=='realCounter':
            simulation_results = {
                'width': self.width,
                'Normalized_RMSE': Normalized_RMSE,
            }
            # Append the simulation results to the appropriate file
            with open(f'../res/pcl_files/realCounter_increments.pcl', 'ab') as f:
                pickle.dump(simulation_results, f)
        else:
            simulation_results = {
                'width': self.width,
                'Normalized_RMSE': Normalized_RMSE,
            }
            # Append the simulation results to the appropriate file
            with open(f'../res/pcl_files/CEDAR_increments.pcl', 'ab') as f:
                pickle.dump(simulation_results, f)

    def calculateReadRMSEForDifferentCounters(self):
        """
        This simulation calculates the Read Root Mean Square Error (RMSE) and Normalized RMSE for different counters types.
        It stores the real value at the index of the flow and it queries the estimated value of the flow from the array counter after
        incrementing it using the methodology of count min sketch.
        """
        numOfIncrements = 1000
        realValCntr = np.zeros(self.num_flows)
        relativeError = [0]* numOfIncrements
        # Increment the counters randomly and update the real values

        for i in range(numOfIncrements):
            # Choose a random flow to increment
            flow = np.random.randint(self.num_flows)
            realValCntr[flow] += 1
            for seed in range(self.depth):
                counterIndex = self.pair.index(f'{seed},{mmh3.hash(str(flow), seed) % self.width}')
                if self.mode in ['Morris', 'CEDAR', 'realCounter']:
                    self.counter_type.incCntr(cntrIdx=counterIndex, factor=1)
                else:
                    print(f'Sorry, the mode {self.mode} that you requested is not supported')
            relativeError[i] =abs(realValCntr[flow] - self.queryFlow(flow)) / realValCntr[flow]
        # Compute the Root Mean Square Error (RMSE) and Normalized RMSE over all flows using the error_value list
        RMSE = math.sqrt(sum(error ** 2 for error in relativeError) / numOfIncrements)
        Normalized_RMSE = RMSE / numOfIncrements
        if self.mode=='Morris':
            # Write the results to a file and return them as a dictionary
            simulation_results = {
                'width': self.width,
                'Normalized_RMSE': Normalized_RMSE,
            }
            # Append the simulation results to the appropriate file
            with open(f'../res/pcl_files/Morris_increments.pcl', 'ab') as f:
                pickle.dump(simulation_results, f)
        elif self.mode=='realCounter':
            simulation_results = {
                'width': self.width,
                'Normalized_RMSE': Normalized_RMSE,
            }
            # Append the simulation results to the appropriate file
            with open(f'../res/pcl_files/realCounter_increments.pcl', 'ab') as f:
                pickle.dump(simulation_results, f)
        else:
            print(f'Sorry, the mode {self.mode} is not found in the counter types list')

def main():
    """
    This starts by removing the old files if they exists using remove_existing_files method. It iterate the first for loop
    based on the list of counter types and tranfer the counter type as a mode to the CountMinSketch class to calculate Normalized_RMSE 
    using Morris or real counter.
    """
    num_flows=20
    remove_existing_files()
    depth = 2  # default value
    counter_types = ['Morris', 'realCounter']
    for counter_type in counter_types:
        for width in range(2, num_flows//2):
            cmc = CountMinSketch(width=width, depth=depth, num_flows=num_flows, mode=counter_type)
            cmc.calculateReadRMSEForDifferentCounters()
    plot_counters()


if __name__ == '__main__':
    main()
