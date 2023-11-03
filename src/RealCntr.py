
class CntrMaster(object):
    """
    This class manages a set of counters.
    """

    def __init__(self, numCntrs):
        """
        Constructor method that initializes the CntrMaster object with a specified number of counters.
        """
        self.numCntrs = numCntrs  # number of counters in the flow array, which is width*depth
        self.cntrs = [0 for i in range(self.numCntrs)]  # repeat the zero counter number of counter times

    def rstAllCntrs(self):
        """
        Resets all the counter values to zero.
        """
        self.cntrs = [0 for i in range(self.numCntrs)]

    def incCntr(self, cntrIdx=0, factor=1, mult=False, verbose=[]):
        """
        Increments or multiplies a specific counter value based on the provided parameters.
        Returns the updated counter value as a dictionary of integer and binary.
        """
        self.cntrs[cntrIdx] = (self.cntrs[cntrIdx] * factor) if mult else (self.cntrs[cntrIdx] + factor)
        return {'cntrVec' : bin(self.cntrs[cntrIdx])[2:], 'val' : self.cntrs[cntrIdx]}

    def queryCntr(self, cntrIdx):
        """
        Retrieves the value of a specific counter from the counters list.
        Returns the counter value as a dictionary of integer and binary.
        """
        return {'cntrVec' : bin(self.cntrs[cntrIdx])[2:], 'val' : self.cntrs[cntrIdx]}

