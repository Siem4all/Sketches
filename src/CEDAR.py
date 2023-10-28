import random
import numpy as np

from printf import printf
import settings
# import commonFuncs 

# The 'delta' parameter determines CEDAR's accuracy.
# Given a counter size and maximum value to count, the function findMinDeltaByMaxVal finds the minimal delta. using binary search.
# To prevent overflows, the search range should be limited. 
# This is done using the list of dicts below. 
deltaSearchRanges = [
                 {'cntrSize' : 5,    'deltaLo' : 0.0001,    'deltaHi' : 0.3},
                 {'cntrSize' : 6,    'deltaLo' : 0.0001,    'deltaHi' : 0.3},
                 {'cntrSize' : 7,    'deltaLo' : 0.0001,    'deltaHi' : 0.3},
                 {'cntrSize' : 8,    'deltaLo' : 0.0001,    'deltaHi' : 0.3},
                 {'cntrSize' : 9,    'deltaLo' : 0.0001,    'deltaHi' : 0.2},
                 {'cntrSize' : 10,   'deltaLo' : 0.0001,    'deltaHi' : 0.2},
                 {'cntrSize' : 11,   'deltaLo' : 0.0001,    'deltaHi' : 0.15},
                 {'cntrSize' : 12,   'deltaLo' : 0.0001,    'deltaHi' : 0.13},
                 {'cntrSize' : 13,   'deltaLo' : 0.0001,    'deltaHi' : 0.1},
                 {'cntrSize' : 14,   'deltaLo' : 0.00001,   'deltaHi' : 0.1},
                 {'cntrSize' : 15,   'deltaLo' : 0.00001,   'deltaHi' : 0.08},
                 {'cntrSize' : 16,   'deltaLo' : 0.00001,   'deltaHi' : 0.07},
                 ]

class CntrMaster(object):
    """
    Generate, check and parse counters
    """
    # Generates a string that details the counter's settings (param vals).
    genSettingsStr = lambda self : 'Cedar_n{}_d{:.6f}'.format(self.cntrSize, self.delta)
    
    # This is the CEDAR formula to calculate the diff given the delta and the sum of the previous diffs
    calc_diff = lambda self, sum_of_prev_diffs: (1 + 2 * self.delta ** 2 * sum_of_prev_diffs) / (1 - self.delta ** 2)

    # print the details of the counter in a convenient way
    printCntrLine = lambda self, cntrSize, delta, numCntrs, mantVal, cntrVal: print('cntrSize={}, delta={}' .format(cntrSize, delta))

    cntr2num = lambda self, i: self.sharedEstimators[i]
    
    calcDiff = lambda self, estimator : (1 + 2*self.delta^2 * estimator) / (1 - self.delta^2)

    def __init__(self, cntrSize=8, delta=None, numCntrs=1, verbose=[], cntrMaxVal=None):
        """
        Initialize an array of cntrSize counters. The cntrs are initialized to 0.
        Inputs:
        cntrSize  - num of bits in each counter.
        Delta - the max relative error. 
        cntrMaxVal - requested max value to be reached by the counter. When Delta is not given, the initiator uses this value,
                     and calculates (using binary search) the minimum delta that allows reaching this maximum value.
        numCntrs - number of counters in the array.
        """
        self.cntrSize      = cntrSize
        self.numCntrs      = numCntrs
        self.numEstimators = 2**self.cntrSize
        self.verbose       = verbose
        self.cntrs = [0 for i in range(self.numCntrs)]
        if (delta==None):
            if (cntrMaxVal==None):
                print ('error: the input arguments should include either delta or cntrMaxVal')
                exit ()
            self.cntrMaxVal = cntrMaxVal
            if (settings.VERBOSE_DETAILS in self.verbose):
                self.detailFile = open ('../log/CEDAR_details.log', 'w')
            self.findMinDeltaByMaxVal(targetMaxVal=self.cntrMaxVal)
            if (settings.VERBOSE_DETAILS in self.verbose):
                printf (self.detailFile, 'cntrSize={}, cntrMaxVal={}, found delta={}\n' .format (self.cntrSize, self.cntrMaxVal, self.delta))
                for i in range(len(self.sharedEstimators)):
                    printf (self.detailFile, 'sharedEstimator[{}]={:.4f}\n' .format(i, self.sharedEstimators[i]))
                
        else:
            self.delta         = delta
            self.calcDiffsNSharedEstimators ()
        
    def calcDiffsNSharedEstimators (self):
        self.sharedEstimators = np.zeros (self.numEstimators)
        self.diffs             = np.zeros (self.numEstimators-1) 
        for i in range (1, 2**self.cntrSize):
            self.diffs[i-1] = self.calc_diff(self.sharedEstimators[i-1])
            self.sharedEstimators[i] = self.sharedEstimators[i-1] + self.diffs[i-1] 
        self.cntrMaxVal = self.sharedEstimators[-1]

    
    def rstCntr (self, cntrIdx=0):
        """
        """
        self.cntrs[cntrIdx] = 0
    def rstAllCntrs(self):
        """
        """
        self.cntrs = [0 for i in range(self.numCntrs)]

    def incCntr(self, cntrIdx=0, factor=1, mult=False, verbose=[]):
        """
        Increase a counter by a given factor.
        Input:
        cntrIdx - index of the cntr to increment, in the array of cntrs.
        mult - if true, multiply the counter by factor. Else, increase the counter by factor.
        factor - the additive/multiplicative coefficient.
        verbose - determines which data will be written to the screen.
        Output:
        cntrDict: a dictionary representing the modified counter where: 
            - cntrDict['cntrVec'] is the counter's binary representation; cntrDict['val'] is its value.
        Operation:
        Define cntrVal as the current counter's value. 
        Then, targetValue = cntrVal*factor (if mult==True), and targetValue = cntrVal + factor (otherwise).  
        If targetValue > maximum cntr's value, return the a cntr representing the max possible value. 
        If targetValue < 0, return a cntr representing 0.
        If targetValue can be represented correctly by the counter, return the exact representation.
        Else, use probabilistic cntr's modification.
        
        If verbose==settings.VERBOSE_DETAILS, the function will print to stdout:
        - the target value (the cntr's current value + factor)
          - cntrDict['cntrVec'] - the binary counter.
          - cntrDict['val']  - the counter's value.
        """
        settings.checkCntrIdx(cntrIdx=cntrIdx, numCntrs=self.numCntrs, cntrType='CEDAR')
        for i in range(factor):
            # The probability to increment is calculated  according to the diff
            if (self.cntrs[cntrIdx] == self.numEstimators-1): # reached the largest estimator --> cannot further inc
                if (settings.VERBOSE_NOTE in self.verbose):
                    print ('note: tried to inc cntr {} above the maximal estimator value of {}' .format (cntrIdx, self.sharedEstimators[-1]))
                break
            probOfFurtherInc = 1/self.diffs[self.cntrs[cntrIdx]]
            if random.random() < probOfFurtherInc:
                if (settings.VERBOSE_DETAILS in verbose): 
                    print ('oldVal={:.0f}, incedVal={:.0f}, probOfFurtherInc={:.6f}'
                            .format (self.sharedEstimators[self.cntrs[cntrIdx]], self.sharedEstimators[self.cntrs[cntrIdx]+1], probOfFurtherInc))
                self.cntrs[cntrIdx] += 1

        return {'cntrVec': np.binary_repr(self.cntrs[cntrIdx], self.cntrSize), 'val': self.sharedEstimators[self.cntrs[cntrIdx]]}

    def queryCntr(self, cntrIdx=0) -> dict:
        """
        Query a cntr.
        Input:
        cntrIdx - the counter's index.
        Output:
        cntrDic: a dictionary, where:
            - cntrDict['cntrVec'] is the counter's binary representation; cntrDict['val'] is its value.
        """
        settings.checkCntrIdx(cntrIdx=cntrIdx, numCntrs=self.numCntrs, cntrType='CEDAR')
        return {'cntrVec': np.binary_repr(self.cntrs[cntrIdx], self.cntrSize), 'val': self.sharedEstimators[self.cntrs[cntrIdx]]}

    def findMinDeltaByMaxVal (self, targetMaxVal):
        """
        Given a target maximum countable value, return the minimal 'delta' parameter that reaches this value, 
        for the current counter's size.
        delta value determines the expected error: a higher delta implies a higher estimated error.
        The min necessary delta is found through a binary search.
        Inputs:   
        * deltaLo - initial lower val for the binary search
        * deltaHi - initial higher val for the binary search
        * resolution = minimum difference (deltaHi-deltaLo); when reached - break the binary search.
        """

        deltaSearchRange = [item for item in deltaSearchRanges if item['cntrSize']==self.cntrSize]
        if len(deltaSearchRange)==0:
            print ('Sorry, but the requested cntrSize {self.cntrSize} is currently not supported by CEDAR')
            return
        deltaLo, deltaHi = deltaSearchRange[0]['deltaLo'], deltaSearchRange[0]['deltaHi']
        resolution = deltaLo

        # check first the extreme cases
        self.delta = deltaHi
        self.calcDiffsNSharedEstimators ()
        if (self.cntrMaxVal < targetMaxVal):
            print ('cannot reach maxVal={} even with highest delta, deltaHi={}. Skipping binary search' .format (targetMaxVal, deltaHi))
            return

        while (True):
            if (deltaHi - deltaLo < resolution): # converged. Still, need to check whether this delta is high enough.
                self.calcDiffsNSharedEstimators ()
                if (self.cntrMaxVal >= targetMaxVal): # can reach maxVal with this delta --> Good
                    return
                # now we know that cannot reach targetMaxVal with the current delta
                self.delta += resolution
                self.calcDiffsNSharedEstimators ()
                if (self.cntrMaxVal < targetMaxVal): 
                    print ('problem at binary search')
                    exit ()
                return
                
            self.delta = (deltaLo + deltaHi)/2
            if (settings.VERBOSE_DETAILS in self.verbose):
                printf (self.detailFile, 'delta={}\n' .format (self.delta))
            self.calcDiffsNSharedEstimators ()
            if (self.cntrMaxVal==targetMaxVal): # found exact match 
                break
            if (self.cntrMaxVal < targetMaxVal): # can't reach maxVal with this delta --> need larger delta value
                deltaLo = self.delta
            else: # maxVal > targetMaxVal --> reached the maximum value - try to decrease delta, to find a tighter value.
                deltaHi = self.delta
        settings.error (self.delta) #$$$
        return self.delta             


def printAllVals(cntrSize=8, delta=None, cntrMaxVal=None, verbose=[]):
    """
    Loop over all the binary combinations of the given counter size.
    For each combination, print to file the respective counter, and its value.
    The prints are sorted in an increasing order of values.
    """
    listOfVals = []
    myCntrMaster = CntrMaster(cntrSize=cntrSize, delta=delta, cntrMaxVal=cntrMaxVal, numCntrs=1)
    for num in range(2 ** cntrSize):
        val = myCntrMaster.cntr2num(num)
        listOfVals.append ({'cntrVec' : np.binary_repr(num, cntrSize), 'val' : val})


    if settings.VERBOSE_RES in verbose:
        outputFile = open('../res/{}.res'.format(myCntrMaster.genSettingsStr()), 'w')
        for item in listOfVals:
            printf(outputFile, '{}={:.1f}\n'.format(item['cntrVec'], item['val']))

