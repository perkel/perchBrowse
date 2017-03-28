import numpy as np

def findEvents(data, thr, lowBins, highBins):
# Takes a time series that is sometimes below and sometimes above a threshold
# to detect an event or "bout", the signal needs to be above threshold for highBins
# surrounded by silence for at least lowBins

    n = len(data)
    boutIndex = 0
    boutStart = []
    boutEnd = []
    
    i = 0 #i is index to the current data point
    dtf = np.array(data) > thr # make array of true/false values
    # first find all candidate silent periods
    while i < n:
        
        if max(data([i:i+lowBins-1])) < thr: # found a silent period sufficiently long
            # mark beginning; look for end of silent period
            boutEnd.append(i)
            i += lowBins
            
            while j < nChunks: # go until we find a threshold crossing
                if data[j] > thr: # candidate bout start
                    if min(data[j:j+highBins]) > thr: # all above threshold
    #                    boutStart.append(j)
    #                else: # false alarm
    #                    print "false alarm"
    #            j += 1
    #    
    #            # when we find it, mark the end. Now above threshold
    #            # if above threshold for long enough, start a new event.
    #
    #            # m
                        
        
        i+= 1
        
    print "boutEnd ", boutEnd
