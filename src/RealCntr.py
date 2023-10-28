import numpy as np


class CntrMaster(object):

    cntr2num = lambda self, cntrVec: int(cntrVec, 2)

    def __init__(self, cntrSize, numCntrs):
        """
        first i have initialized  all the counters which has cntrSize bits to zero. Eg '000000' if cntrSize is 6.
        As count min sketch is dimensional array with rows equals to number of depth and columns equals to number of width,
         i have converted it to row and column pair list.
        """
        self.cntrSize    = cntrSize
        self.numCntrs    = numCntrs  # number of counters in the flow array, which is width*depth
        self.cntrZeroVec = '0' * self.cntrSize  # Initialize all the counter to 0, Eg. '0000'
        self.cntrs       = [self.cntrZeroVec for i in range(self.numCntrs)]  # repeat the zero counter number of counter times

    def rstAllCntrs(self):
        """
        """
        self.cntrs = [self.cntrZeroVec for i in range(self.numCntrs)]

    def incCntr(self, cntrIdx=0, factor=1, mult=False, verbose=[]):
        """
        This converts the counter binary value to integer and increment the value by factor. Return the binary and integer value
        as dictionary.
        """
        self.cntrs[cntrIdx] = str(bin(int(self.cntrs[cntrIdx], 2) + factor)[2:])
        return {'cntrVec': self.cntrs[cntrIdx], 'val': self.cntr2num(self.cntrs[cntrIdx])}

    def queryCntr(self, cntrIdx):
        """

        Here i used the variable flowIdx to get the binary number from counters list, and converted it to number
        """
        return {'cntrVec': self.cntrs[cntrIdx], 'val': self.cntr2num(self.cntrs[cntrIdx])}





