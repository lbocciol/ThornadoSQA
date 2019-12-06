import time

TimerOsc, TimerFreq, TimerIMEX = 0., 0., 0.
TimerTemp = 0.

def start():

    global TimerTemp
    TimerTemp = time.time()

def stop( Timer ):

    global TimerTemp
    Timer += time.time() - TimerTemp
    
    TimerTemp = 0.
    
    return Timer
