import numpy as np


class CEDARCntr(object):

    def __init__(self, cntrSize, numCntrs, cntrMaxVal):
        """
        This counting mechanism separates the counter estimation into two arrays, namely flow and shared estimator array. Here, i have calculated
        the values estimator array using delta and the array difference between estimators. at the beginning, i initialize the first and the last
         estimator values (i.e. A[0] and A[L-1] to 0 and cntrMaxVal respectively). As my primary interest is to apply the CEDAR counter
         to Count Min Sketch, i need to make any flow to point to estimated array number of hash function times. The flow array increments its value (pointer)
         at the arrival of the same flow but the value of the shared estimator array is fixed so that i need to apply the count min sketch to the
         flow array and it the end, i will use the min value of the flow array as a pointer to the estimator.
        """
        self.cntrSize = cntrSize
        self.estimator_size = 2 ** cntrSize  # The estimator size is 2**q, where q is counter size
        self.numCntrs = numCntrs  # number of counters in the flow array, which is width*depth
        self.cntrZeroVec = '0' * self.cntrSize  # Initialize all the counter to 0, Eg. '0000'
        self.cntrs = [self.cntrZeroVec for i in range(self.numCntrs)]  # repeat the zero counter number of counter times
        self.cntrMaxVec = '1' * self.cntrSize  # the max counter vector is '1111' or 2**cntrSize(4)=15
        self.calcCntrMaxVal = (1 << self.cntrSize) - 1  # the max counter value which is 2**cntrSize
        self.estimate_array = np.zeros(self.estimator_size)  # Initialize the estimator array to zero
        self.estimate_array[self.estimator_size - 1] = cntrMaxVal
        self.array_diff = np.zeros(self.estimator_size - 1)  # array of differences between estimators
        self.delta = 0.01
        self.array_diff[0] = 1  # set first difference
        for l in range(1, self.estimator_size - 1):
            self.array_diff[l] = 1 + 2 * (self.delta ** 2) * sum(
                [self.array_diff[j] for j in range(l)])  # set remaining differences
        for i in range(1, self.estimator_size - 1):
            self.estimate_array[i] = self.estimate_array[i - 1] + self.array_diff[i - 1]  # D[i]=A[i+1]-A[i]

    def incCntr(self, cntrIdx, factor=1, mult=False, verbose=[]):
        """

        This converts the counter binary value to integer and check if that value can increment or reaches its max value. If it not reaches max
        value, i add 1 to the target value and save it as binary.
        """
        targetVal = int(self.cntrs[cntrIdx], 2)
        if targetVal < self.calcCntrMaxVal:
            self.cntrs[cntrIdx] = str(bin(targetVal + factor)[2:].zfill(self.cntrSize))

        return self.estimate_array[int(self.cntrs[cntrIdx], 2)]

    def queryCntr(self, cntrIdx):
        """

        Here i used the variable flowIdx to get the binary number from cntr list, i converted it to number in order to use it as
         a pointer to the estimator array to get the value of the estimator array.
        """
        return self.estimate_array[int(self.cntrs[cntrIdx], 2)]


