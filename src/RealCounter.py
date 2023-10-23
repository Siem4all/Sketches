import numpy as np


class RealCntr(object):

    def __init__(self, cntrSize, numCntrs):
        """
        first i have initialized  all the counters which has cntrSize bits to zero. Eg '000000' if cntrSize is 6.
        As count min sketch is dimensional array with rows equals to number of depth and columns equals to number of width,
         i have converted it to row and column pair list.
        """
        self.numCntrs    = numCntrs  # number of counters in the flow array, which is width*depth
        self.cntrs    = [0 for i in range(self.numCntrs)]  # repeat the zero counter number of counter times

    def rstAllCntrs(self):
        """
        """
        self.cntrs = [0 for i in range(self.numCntrs)]
    def rstCntr (self, cntrIdx=0):
        """
        """
        self.cntrs[cntrIdx] = 0

    def incCntr(self, cntrIdx=0, factor=1, mult=False, verbose=[]):
        """

        This converts the counter binary value to integer and check if that value can increment or reaches its max value. If it not reaches max
        value, it added 1 to the target value and save it as binary.
        """
        self.cntrs[cntrIdx] += factor
        return self.cntrs[cntrIdx]

    def queryCntr(self, cntrIdx):
        """

        Here i used the variable flowIdx to get the binary number from counters list, and converted it to number
        """
        return self.cntrs[cntrIdx]



