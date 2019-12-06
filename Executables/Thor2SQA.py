from Thor2SQA import thornadosqainterfacemodule as T2SQA
from Thor2SQA import programstartendmodule as Prog
from FrequencyUtils import FrequencyAnalysis
from MatrixUtils import UnitarityCheck
from Timers import TimerOsc, TimerFreq, TimerIMEX
import Timers
import numpy as np
#import sys

UnitarityTolerance = 1.0e-4

(Kilometer, Second, MeV, SpeedOfLightCGS) = Prog.initializeunits()

Centimeter = Kilometer * 1e-5
MilliSecond = Second * 1e-3

nNodes = 2
nSpecies = 4
nF = 2
nM = 2

nX    = np.array([ 128, 1, 1 ])
xL    = np.array([ 0.0 * Kilometer, 0.0, 0.0 ])
xR    = np.array([ 500.0 * Kilometer, np.pi, 2.0 * np.pi  ])
bcX   = np.array([ 32, 0, 0 ])
swX   = np.array([  1, 1, 1 ])
ZoomX = np.array([ 1.0, 1.0, 1.0 ])

nE    = 70
nE_G  = nE * nNodes
eL    = np.array( 2.5e0 ) * MeV
eR    = np.array( 2.5e2 ) * MeV
bcE   = 0
swE   = 0
ZoomE = np.array( 1.0 )

TimeSlice = np.array( 0.3 )

t      = np.array( 0.0 )
dt_wrt = np.array( 2.0e-2 ) * MilliSecond
t_SQA  = np.array( 3.0e-9 ) * MilliSecond
t_end  = np.array( 3.0e0 )  * MilliSecond

iCycleD = 1

WriteGF, WriteFF, WriteRF, WriteOp = False, False, True, False

# Do Initialization Stuff
Prog.programstart( nNodes, nSpecies, nF, \
              nX, xL, xR, bcX, swX, ZoomX, \
              nE, eL, eR, bcE, swE, ZoomE, \
              TimeSlice, t )

nDOFX = nNodes
R_Shock = T2SQA.getshockradius()

R   = np.zeros((nDOFX, nX[0], nX[1], nX[2])) #Km
dR  = np.zeros((nX[0], nX[1], nX[2])) #Km
Ye  = np.zeros((nDOFX, nX[0], nX[1], nX[2]))
Rho = np.zeros((nDOFX, nX[0], nX[1], nX[2])) #g/cm^3

for iNodeX in range(nDOFX):
    for iX1 in range(nX[0]):  
        for iX2 in range(nX[1]):
            for iX3 in range(nX[2]):

                R[iNodeX,iX1,iX2,iX3], \
                dR[iX1,iX2,iX3], \
                Rho[iNodeX,iX1,iX2,iX3], \
                Ye[iNodeX,iX1,iX2,iX3] = T2SQA.get_profile(iNodeX+1,iX1+1,iX2+1,iX3+1)

# Do Evolution
iCycle   = 0
t_wrt    = dt_wrt
wrt      = False
StartSQA = False
StopSQA  = False
tmax_freq = np.array( 1.0e-9, dtype=np.float64 ) # Set time (in seconds) safely above period 

while t < t_end:

    iCycle += 1

    dt = 0.5 * np.min( dR ) / ( 2.0 * (nNodes - 1) + 1.0 )

    if t + dt > t_SQA:
      StartSQA = True

    if t + dt > t_wrt :
      dt    = t_wrt - t
      t_wrt = t_wrt + dt_wrt
      wrt   = True

    if iCycle % iCycleD == 0:

      print('%3s' %'', 'Cycle =','%8.8i' %iCycle, \
            ' t =',  '%1.6e' %(t / Second),  '(s)', \
            ' dt =', '%1.6e' %(dt / Second), '(s)')

    if t + dt > t_end:
      dt = t_end - t

# SQA PART -------------------
    if StartSQA and not StopSQA:

      print('Starting SQA oscillations')
      StopSQA = True

      DoOscillations = False
      InitializeFirstZone = True
      
      T2SQA.resetsmatrix()
        
      STotal = np.array([[ np.identity(nF, dtype=np.complex128) for iN_E in range(nE_G) ] for iM in range(nM)])
      iX = 0
      for iX1 in range(nX[0]):
        for iX2 in range(nX[1]):
          for iX3 in range(nX[2]):
            for iNodeX in range(nDOFX):
                      
              # --- OSCILLATIONS DRIVER --------------------------------------------------
              #os.system('rm ../Output/Relaxation_SQA_f*')
              iX += 1
              
              if R[iNodeX,iX1,iX2,iX3] > 50.0e0 * Kilometer and not DoOscillations :
                DoOscillations = True
          
              if  R[iNodeX,iX1,iX2,iX3] > R_Shock:                
                DoOscillations = False
                  
              if DoOscillations:
               
                if InitializeFirstZone:                
                  T2SQA.initializefmatrixosc( iNodeX+1,iX1+1,iX2+1,iX3+1, STotal, nM, nE_G, nF )
                  InitializeFirstZone = False
        
                else:
                  # iX1-1 because you initialize fMatrixOsc by getting the 
                  # spectrum coming from the previous zone
                  T2SQA.initializefmatrixosc( iNodeX+1,iX1,iX2+1,iX3+1, STotal, nM, nE_G, nF )
        
                dr = dR[iX1,iX2,iX3] / nNodes / Centimeter 
                
                tZone = dr / SpeedOfLightCGS 

                T2SQA.initializeoscinterface( Rho[iNodeX,iX1,iX2,iX3], Ye[iNodeX,iX1,iX2,iX3])
                
                # CALL OSCILLATIONS TO FIND PERIOD
                FinishOsc = False
                tOsc = np.array( 0.0 )
                dt_osc = np.array( 1.0e-13 )
                tmax_freq = dt_osc * 1.e3
                f, t_freq = [], []
                while not FinishOsc:
                  
                  tmax_freq = dt_osc * 1.e3
                  if tOsc + dt_osc >= tmax_freq:    
                    dt_osc = tmax_freq - tOsc
                    FinishOsc = True

                  Timers.start()
                  ftemp, Stemp = T2SQA.oscillationsinterface( tOsc, dt_osc, nM, nE_G, nF )
                  TimerOsc = Timers.stop( TimerOsc )  

                  t_freq.append( float(tOsc) )
                  f.append(ftemp)
            
                print( 'tmax',tmax_freq )
                Timers.start()
                Period = FrequencyAnalysis( np.array(t_freq) , np.array(f) )
                TimerFreq = Timers.stop( TimerFreq )

                T2SQA.finalizeoscinterface()
                
                # -------------------------------------
                # RESTART OSCILLATIONS FOR ~100 PERIODS
                T2SQA.initializefmatrixosc( iNodeX+1,iX1,iX2+1,iX3+1, STotal, nM, nE_G, nF )

                T2SQA.initializeoscinterface( Rho[iNodeX,iX1,iX2,iX3], Ye[iNodeX,iX1,iX2,iX3])

                # Run for 100 periods
                Period = Period * 100

                nPeriods = tZone // Period
            
                print('PERIOD', Period, dt_osc, nPeriods)
                STotalThisZone = np.array([[ np.identity(nF, dtype=np.complex128) for iN_E in range(nE_G) ] for iM in range(nM)])
                FinishOsc = False
                tOsc = np.array( 0.0 )
                dt_osc = np.array( 1.0e-13 )
                while not FinishOsc:

                  if tOsc + dt_osc >= Period:

                    dt_osc = tmax_freq - tOsc
                    FinishOsc = True
                  
                  Timers.start()
                  ftemp, Stemp = T2SQA.oscillationsinterface( tOsc, dt_osc, nM, nE_G, nF )
                  TimerOsc = Timers.stop( TimerOsc )
                    
                  STotalThisZone = np.array([[ np.matmul( STotalThisZone[m,iN_E,:,:], Stemp[m,iN_E,:,:]) \
                                      for iN_E in range(nE_G) ] for m in range(nM) ] )
                
                  for m in range(nM):
                    for iN_E in range(nE_G):
                      UnitarityCheck( STotalThisZone[m,iN_E,:,:], UnitarityTolerance, message='STotalThisZone1' )
                
                STotalThisZone = np.array( [[np.linalg.matrix_power( STotalThisZone[iM,iN_E,:,:], np.int(nPeriods) ) \
                                      for iN_E in range(nE_G)] for iM in range(nM)] )

                UnitarityCheck( STotalThisZone[m,iN_E,:,:], UnitarityTolerance, message='STotalThisZone2' )
                
                # CALL OSCILLATIONS FOR THE REMAINDER OF THE ZONE
                tRemainder = tZone % Period
                FinishOsc = False
                tOsc = np.array( 0.0 )
                dt_osc = np.array( 1.0e-13 )
                while not FinishOsc:

                  if tOsc + dt_osc >= tRemainder:

                    dt_osc = tmax_freq - tOsc
                    FinishOsc = True

                  Timers.start()
                  ftemp, Stemp = T2SQA.oscillationsinterface( tOsc, dt_osc, nM, nE_G, nF )
                  TimerOsc = Timers.stop( TimerOsc )

                  STotalThisZone = np.array([[ np.matmul( STotalThisZone[m,iN_E,:,:], Stemp[m,iN_E,:,:]) \
                                      for iN_E in range(nE_G) ] for m in range(nM) ] )

                T2SQA.finalizeoscinterface()

                STotal = np.array([[ np.matmul( STotal[m,iN_E,:,:], STotalThisZone[m,iN_E,:,:]) \
                                      for iN_E in range(nE_G) ] for m in range(nM) ] )
               
                for m in range(nM):
                    for iN_R in range(nE_G):
                        UnitarityCheck( STotal[m,iN_E,:,:], UnitarityTolerance, message='STotal' )

                
                #T2SQA.calculateopacitiesosc( dr, iX, iNodeX+1, iX1+1, iX2+1, iX3+1, np.array( [[0.9999 for iN_E in range(nE_G)] for m in range(nM)]), nM, nE_G )
                T2SQA.calculateopacitiesosc( dr, iX, iNodeX+1, iX1+1, iX2+1, iX3+1, np.abs(STotal[:,:,0,0])**2, nM, nE_G )
        
                print( 'Done with Zone', iNodeX+1, iX1+1, iX2+1, iX3+1 )
                    
#---------------------------------------------------------
    Timers.start()
    T2SQA.integrationdriver( dt )
    TimerIMEX = Timers.stop( TimerIMEX )

    t = t + dt
    
    if wrt:
        Prog.dumpfields( t, WriteGF, WriteFF, WriteRF, WriteOp )
        wrt = False

# Finalize Stuff
Prog.programend( t )

print('%3s' %'', 'Python Timers')
print('%3s' %'', 'TimerOsc', '%1.6e' %TimerOsc  )
print('%3s' %'', 'TimerFreq','%1.6e' %TimerFreq )
print('%3s' %'', 'TimerIMEX','%1.6e' %TimerIMEX )
