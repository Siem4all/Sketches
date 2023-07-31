import math, mmh3, os, pickle, numpy as np
from Morris import CntrMaster
from PclFileParser import plot_counters
from printf import printf


def remove_existing_files():
    # Removes the specified files if they exist.
    if os.path.exists('../res/pcl_files/first_increments.pcl'):
        os.remove('../res/pcl_files/first_increments.pcl')

    if os.path.exists('../res/Normalized_error.res'):
        os.remove('../res/Normalized_error.res')


class SketchStatistics:
    def __init__(self, width, depth, num_flows, numOfIncrements):
        """
        SketchStatistics: A Count-Min Sketch data structure for estimating flow frequencies.
        Parameters:
        - width (int): the number of counters per row.
        - depth (int): the number of rows or hash functions.
        - num_flows (int): the total number of flows to be estimated.
        - cntrMaxVal (int): the maximum counter value.
       Attributes:
        - counters: a CntrMaster object that holds the counter values for the sketch.
        - row_column_pair: a tuple of row and column indices used to compute the corresponding counter index.
       Notes:
        - The row_column_pair attribute is used to compute the index of each counter. Instead of adding row with column which gives us wrong mapping, i have
        paired the depth number with the hash value and i have inserted them to the pair list to use the index of the pair as a counter index.
         """
        self.width, self.depth, self.num_flows, self.numOfIncrements = width, depth, num_flows, numOfIncrements  # depth is equals to number of hash functions
        self.counters = CntrMaster(cntrSize=4, numCntrs=depth * width, a=10, cntrMaxVal=1000, verbose=[])
        self.row_column_pair = [f"{i},{j}" for i in range(depth) for j in range(width)]

    def incFlow(self, flow):
        # Increment the counters for the given flow by hashing and updating the appropriate counter values
        for seed in range(self.depth):
            index_in_1D_array = self.row_column_pair.index(f'{seed},{mmh3.hash(str(flow), seed) % self.width}')
            self.counters.incCntr(cntrIdx=index_in_1D_array, factor=1, verbose=[], mult=False)

    def queryFlow(self, flow):
        # Query the minimum counter value for the given flow by hashing and finding the minimum value among the appropriate counters
        min_num = math.inf
        for seed in range(self.depth):
            index_in_1D_array = self.row_column_pair.index(f'{seed},{mmh3.hash(str(flow), seed) % self.width}')
            estimate = self.counters.queryCntr(index_in_1D_array)
            min_num = min(estimate['val'], min_num)
        return min_num

    def runSimulationAndCalculateNormalizedRMSE(self):
        """"
        This simulation calculates the Root Mean Square Error (RMSE) and Normalized RMSE for different number of increments.
        i have updates the counter self.numOfIncrements times by generating random flow from 0 to self.num_flows.
        It stores the real value at the index of the flow and it queries the estimated value of the flow from the array counter after
        incrementing it using the methodology of count min sketch.
        """
        realValCntr = np.zeros(self.num_flows)
        error_value = []
        # Increment the counters randomly and update the real values
        for i in range(self.numOfIncrements):
            # Choose a random flow to increment
            flow = np.random.randint(self.num_flows)
            # Increment the corresponding counter in the Count-Min Sketch
            self.incFlow(flow)
            # Update the corresponding real value counter
            realValCntr[flow] += 1
            # Compute the error between the estimated and real frequencies for the flow and append it to the error_value list
            error_value.append(abs(realValCntr[flow] - self.queryFlow(flow)) / realValCntr[flow])
        # Compute the Root Mean Square Error (RMSE) and Normalized RMSE over all flows using the error_value list
        RMSE = math.sqrt(sum(error ** 2 for error in error_value) / self.numOfIncrements)
        Normalized_RMSE = RMSE / self.numOfIncrements
        # Write the results to a file and return them as a dictionary
        resFile = open(f'../res/Normalized_error.res', 'a+')
        printf(resFile, '\nnumOfIncrements={},  numCntrs={}, RMSE={}, Normalized_RMSE={}\n'.format(self.numOfIncrements,
                                                                                               self.width * self.depth,
                                                                                               RMSE,
                                                                                               Normalized_RMSE))
        simulation_results = {
            'counters': self.width*self.depth,
            'Normalized_RMSE': Normalized_RMSE,
            'numOfIncrements': self.numOfIncrements
        }
        # Append the simulation results to the appropriate file
        with open(f'../res/pcl_files/first_increments.pcl', 'ab') as f:
            pickle.dump(simulation_results, f)


        # Overall, this function simulates a Count-Min Sketch and calculates the Normalized RMSE over a set of flows.
        # It provides a way to evaluate the accuracy of the Count-Min Sketch and compare different increments.
        # The Normalized RMSE is a useful metric for measuring the accuracy of the Count-Min Sketch in estimating the frequencies.

def main():
    remove_existing_files()
    num_flows=500  # default value
    width = 5  # default value
    depth = 2   # default value
    maxIncrement=10000
    for num in range(100, maxIncrement, 500):
        cmc = SketchStatistics(width=width, depth=depth, num_flows=num_flows, numOfIncrements=num)
        cmc.runSimulationAndCalculateNormalizedRMSE()
    plot_counters()


if __name__ == '__main__':
    main()
