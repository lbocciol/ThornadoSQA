import numpy as np
import sys

def UnitarityCheck( Matrix, tolerance, message = 'None'):

    nF, = np.shape( Matrix[0] )

    det = np.linalg.det(Matrix)
    
    if np.abs( np.abs( det ) - 1. ) > tolerance:
        print( message )
        print( 'Determinant exceeds tolerance = ' + str(tolerance) + '; Det =', det )
        sys.exit()

    SumOfColumns = np.array( [ np.sum( np.abs(Matrix[iF1,:])**2 ) for iF1 in range(nF) ] )

    for iF,Sum in zip(range(nF),SumOfColumns):
        if np.abs( Sum - 1.0 ) > tolerance:
            print( message )
            print( 'Sum of row', iF,'is too big:', Sum)
            sys.exit()
