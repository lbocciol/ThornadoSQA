import numpy as np
import sys
import time

def FrequencyAnalysis( Time, fMatrixOsc ):
    
    timer_tot = time.time()
    print( Time )
    
    Tolerance = 5.0e-2 
   
    nTime, nM, nE, nF, nF = np.shape(fMatrixOsc)
    nComplex = 2

    nSample = nTime
    s = np.linspace( 0.0, max(Time), nSample)
    
    f = np.zeros((nM,nE,nF,nF,nComplex,nSample))

    timer_interp = time.time()
    #do interpolations -------------------------------

    f[:,:,:,:,0,:] = np.array([[[[np.interp(s, Time, np.real(fMatrixOsc[:,iM,iE,iF1,iF2])) for iF1 in range(nF)] for iF2 in range(nF)] for iE in range(nE) ] for iM in range(nM)])
        
    f[:,:,:,:,1,:] = np.array([[[[np.interp(s, Time, np.imag(fMatrixOsc[:,iM,iE,iF1,iF2])) for iF1 in range(nF)] for iF2 in range(nF)] for iE in range(nE) ] for iM in range(nM)])
        
    timer_interp = time.time() - timer_interp
    # -------------------------------------------------

    #sampling interval
    T = s[1] - s[0]
    freq = np.linspace(0, 1 / T, nSample)

    period = 1 / (freq)

    timer_fft = time.time()
    
    ffte  = np.zeros((nSample,nM,nE))
    fftmu = np.zeros((nSample,nM,nE))
    for iM in range(nM):
        for iE in range(nE):

            ffte [:,iM,iE] = np.fft.fft( f[iM,iE,0,0,0,:] )
            fftmu[:,iM,iE] = np.fft.fft( f[iM,iE,1,1,0,:] )
            
            # Now perform normalization
            ffte [:,iM,iE] = np.abs(ffte) [:,iM,iE]/ np.sum(np.abs(ffte) [:,iM,iE] )
            fftmu[:,iM,iE] = np.abs(fftmu)[:,iM,iE]/ np.sum(np.abs(fftmu)[:,iM,iE] )

    timer_fft = time.time() - timer_fft

    #Calculate Periods:
    MaxInde  = [np.argmax(ffte [1:nTime//2,iM,iE]) + 1 for iE in range(nE) for iM in range(nM)]
    MaxIndmu = [np.argmax(fftmu[1:nTime//2,iM,iE]) + 1 for iE in range(nE) for iM in range(nM)]

    MaxInd = MaxInde + MaxIndmu

    print(MaxInd)
    timer_tot = time.time() - timer_tot
    print( 'tot interp fft',timer_tot, timer_interp, timer_fft)
    # Get most frequent index from MaxInd

    return period[ max(set(MaxInd), key = MaxInd.count) ]
