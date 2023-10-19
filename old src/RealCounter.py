import numpy as np


class RealCntr(object):

    def __init__(self, cntrSize, numCntrs):
        """
        first i have initialized  all the counters which has cntrSize bits to zero. Eg '000000' if cntrSize is 6.
        As count min sketch is dimensional array with rows equals to number of depth and columns equals to number of width,
         i have converted it to row and column pair list.
        """
        self.cntrSize    = cntrSize
        self.numCntrs    = numCntrs  # number of counters in the flow array, which is width*depth
        self.cntrZeroVec = '0' * self.cntrSize  # Initialize all the counter to 0, Eg. '0000'
        self.counters    = [self.cntrZeroVec for i in range(self.numCntrs)]  # repeat the zero counter number of counter times
        self.cntrMaxVec  = '1' * self.cntrSize # the max counter vector is '1111' or 2**cntrSize(4)=15
        self.cntrMaxVal  = 1 << self.cntrSize - 1  # the max counter value which is 2**cntrSize

    def incCntr(self, cntrIdx=0, factor=1, mult=False, verbose=[]):
        """

        This converts the counter binary value to integer and check if that value can increment or reaches its max value. If it not reaches max
        value, it added 1 to the target value and save it as binary.
        """
        targetVal = int(self.counters[cntrIdx], 2)
        if targetVal < self.cntrMaxVal:
            self.counters[cntrIdx] = str(bin(targetVal + factor)[2:].zfill(self.cntrSize))
        else:
            self.counters[cntrIdx] = self.cntrMaxVec
        return int(self.counters[cntrIdx], 2)

    def queryCntr(self, cntrIdx):
        """

        Here i used the variable flowIdx to get the binary number from counters list, and converted it to number
        """
        return int(self.counters[cntrIdx], 2)



