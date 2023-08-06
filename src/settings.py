import math, random, os
import numpy as np, scipy.stats as st, pandas as pd

from printf import printf
# import commonFuncs 

VERBOSE_COUT_CNTRLINE = 1
VERBOSE_DEBUG = 2
VERBOSE_RES = 3
VERBOSE_PCL = 4
VERBOSE_DETAILS = 5
VERBOSE_NOTE = 6
VERBOSE_LOG  = 7
VERBOSE_DETAILED_LOG = 8

Confs = [{'cntrSize' : 8,  'cntrMaxVal' :  1488888,  'hyperSize' : 2, 'hyperMaxSize' : 2, 'expSize' : 3, 'tetraSize' : 1,'tetraMaxSize' : 1},
         {'cntrSize' : 9,  'cntrMaxVal' :  2994160,  'hyperSize' : 2, 'hyperMaxSize' : 3, 'expSize' : 4},
         {'cntrSize' : 10, 'cntrMaxVal' :  6004704,  'hyperSize' : 2, 'hyperMaxSize' : 3, 'expSize' : 4},
         {'cntrSize' : 11, 'cntrMaxVal' :  12025792, 'hyperSize' : 2, 'hyperMaxSize' : 3, 'expSize' : 4},
         {'cntrSize' : 12, 'cntrMaxVal' :  24067968, 'hyperSize' : 2, 'hyperMaxSize' : 3, 'expSize' : 4},
         {'cntrSize' : 13, 'cntrMaxVal' :  48152320, 'hyperSize' : 2, 'hyperMaxSize' : 3, 'expSize' : 4},
         {'cntrSize' : 14, 'cntrMaxVal' :  96321024, 'hyperSize' : 2, 'hyperMaxSize' : 3, 'expSize' : 4},
         {'cntrSize' : 15, 'cntrMaxVal' : 192658432, 'hyperSize' : 2, 'hyperMaxSize' : 3, 'expSize' : 4},
         {'cntrSize' : 16, 'cntrMaxVal' : 385333248, 'hyperSize' : 2, 'hyperMaxSize' : 3, 'expSize' : 4}]

# Calculate the confidence interval of an array of values ar, given its avg. Based on 
# https://stackoverflow.com/questions/15033511/compute-a-confidence-interval-from-sample-data
confInterval = lambda ar, avg, conf_lvl=0.99 : st.t.interval (conf_lvl, len(ar)-1, loc=avg, scale=st.sem(ar)) if np.std(ar)>0 else [avg, avg]
   
def idxOfLeftmostZero (ar, maxIdx):
    """
    if the index of the leftmost 0 in the array >= maxIdx, return maxIdx.
    else, return the index of the leftmost 0 in the array.
    """ 
    if (ar == '1' * len(ar)): 
        return maxIdx
    return min (ar.index('0'), maxIdx)
    
def checkCntrIdx (cntrIdx, numCntrs, cntrType):
    """
    Check if the given cntr index is feasible.
    If not - print error msg and exit.
    """
    if (cntrIdx < 0 or cntrIdx>(numCntrs-1)):
        print ('error in {}: wrong cntrIdx. Please select cntrIdx between 0 and {}' .format (cntrType, numCntrs-1))
        exit ()
    
def sortcntrMaxVals ():
    """
    Read the file '../res/cntrMaxVals.txt". Sort it in an increasing fashion of the max cntr vals.
    Print the results to '../res/maxC
    """
    input_file      = open ('../res/cntrMaxVals.txt', 'r')
    lines           = (line.rstrip() for line in input_file) # "lines" contains all lines in input file
    lines           = (line for line in lines if line)       # Discard blank lines
    
    list_of_dicts = []
    for line in lines:

        # Discard lines with comments / verbose data
        if (line.split ("//")[0] == ""):
            continue
        
        splitted_line = line.split ()
        mode = splitted_line[0]
        list_of_dict = [item['mode'] for item in list_of_dicts if item['mode']==mode]
        if (mode not in [item['mode'] for item in list_of_dicts if item['mode']==mode]):
            list_of_dicts.append ({'mode' : mode, 'cntrSize' : int(mode.split('_n')[1].split('_')[0]), 'maxVal' : float(splitted_line[1].split('=')[1])})
    list_of_dicts =  sorted (list_of_dicts, key = lambda item : (item['cntrSize'], item['maxVal']))
    output_file   = open ('../res/cntrMaxValsSorted.txt', 'w')
    for item in list_of_dicts:
        val = item['maxVal']
        if (val < 10**8):
            printf (output_file, '{}\t{:.0f}\n' .format (item['mode'], item['maxVal']))
        else:
            printf (output_file, '{}\t{}\n' .format (item['mode'], item['maxVal']))            

def error (str):
    """
    Print an error msg and exit.
    """
    print (str)
    exit  ()
    