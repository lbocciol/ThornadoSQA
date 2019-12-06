MODULE InterpolationModule

  USE KindModule, ONLY: &
    DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: LinearInterpolation1D
  PUBLIC :: LinearInterpolation2D
  PUBLIC :: QuadraticInterpolation1D

CONTAINS

  SUBROUTINE LinearInterpolation1D(X_Old, Y_Old, nX, &
                                   X_New, Y_new) 

    REAL(DP), INTENT(OUT)   :: Y_New
    REAL(DP), INTENT(IN)    :: X_New
    REAL(DP), INTENT(IN)    :: X_Old(nX)
    REAL(DP), INTENT(IN)    :: Y_Old(nX)
    INTEGER, INTENT(IN)     :: nX
    
    INTEGER :: i
     
    IF (X_New .le. X_Old(1)) THEN

      !Extrapolate down
      Y_New = Y_Old(1) - (X_Old(1) - X_New) * (Y_Old(2) - Y_Old(1)) &
                                            / (X_Old(2) - X_Old(1)) 

      RETURN

    ELSE IF (X_New .ge. X_Old(nX)) THEN
      
      !Exptrapolate up
      Y_New = Y_Old(nX) + (X_New - X_Old(nX)) * (Y_Old(nX) - Y_Old(nX-1)) &
                                              / (X_OLD(nX) - X_Old(nX-1))
      RETURN
    
    ELSE
      DO i = 1,nX
        
        IF(X_Old(i) .ge. X_New) THEN
            Y_New = Y_Old(i-1) + (X_New    - X_Old(i-1)) &
                               * (Y_Old(i) - Y_Old(i-1)) &
                               / (X_Old(i) - X_Old(i-1))
            RETURN
        
        END IF

      END DO
    END IF                   

  END SUBROUTINE LinearInterpolation1D

  SUBROUTINE LinearInterpolation2D( X_Old, Y_Old, Z_Old, & 
                            nX, nY, X_New, Y_new, Z_New )

    REAL(DP), INTENT(OUT)   :: Z_New
    REAL(DP), INTENT(IN)    :: X_New, Y_New
    REAL(DP), INTENT(IN)    :: X_Old(nX), Y_Old(nY)
    REAL(DP), INTENT(IN)    :: Z_Old(nX,nY)
    INTEGER, INTENT(IN)     :: nX, nY

    REAL(DP) :: Z_NewY1, Z_NewY2
    INTEGER :: iX, iY
    
    IF ( X_New <= X_Old(1) ) THEN
       
      IF ( Y_New <= Y_Old(1) ) THEN
  
        ! Extrapolate down along X
        Z_NewY1 = Z_Old(1,1) - ( X_Old(1)   - X_New     ) &
                             * ( Z_Old(2,1) - Z_Old(1,1)) &
                             / ( X_Old(2)   - X_Old(1)  )
  
        Z_NewY2 = Z_Old(1,2) - ( X_Old(1)   - X_New     ) &
                             * ( Z_Old(2,2) - Z_Old(1,2)) &
                             / ( X_Old(2)   - X_Old(1)  ) 
        
        ! Extrapolate down along Y
        Z_New =  Z_NewY1 - ( Y_Old(1) - Y_New   ) &
                         * ( Z_NewY2  - Z_NewY1 ) &
                         / ( Y_Old(2) - Y_Old(1) )
  
      ELSE IF ( Y_New >= Y_Old(nY) ) THEN
  
        ! Extrapolate down along X
        Z_NewY1 = Z_Old(1,nY-1) - ( X_Old(1)      - X_New         ) &
                                * ( Z_Old(2,nY-1) - Z_Old(1,nY-1) ) &
                                / ( X_Old(2)      - X_Old(1)      )
  
        Z_NewY2 = Z_Old(1,nY) - ( X_Old(1)    - X_New        ) &
                              * ( Z_Old(2,nY) - Z_Old(1,nY) ) &
                              / ( X_Old(2)    - X_Old(1)     )
  
        ! Extrapolate up along Y
        Z_New =  Z_NewY2 + ( Y_New     - Y_Old(nY)   ) &
                         * ( Z_NewY2   - Z_NewY1     ) &
                         / ( Y_Old(nY) - Y_Old(nY-1) )
  
      ELSE 
  
        DO iY = 1,nY
          
          IF ( Y_Old(iY) >= Y_New ) THEN
  
            ! Extrapolate down along X
            Z_NewY1 = Z_Old(1,iY-1) - ( X_Old(1)      - X_New         ) &
                                    * ( Z_Old(2,iY-1) - Z_Old(1,iY-1) ) &
                                    / ( X_Old(2)      - X_Old(1)      )
    
            Z_NewY2 = Z_Old(1,iY) - ( X_Old(1)    - X_New        ) &
                                  * ( Z_Old(2,iY) - Z_Old(1,iY) ) &
                                  / ( X_Old(2)    - X_Old(1)     )
        
            ! Interpolate Y 
            Z_New =  Z_NewY1 + (Y_New     - Y_Old(iY-1)) &
                             * (Z_NewY2   - Z_NewY1    ) &
                             / (Y_Old(iY) - Y_Old(iY-1))
  
            RETURN
  
          END IF
      
        END DO
  
      END IF

    ELSE IF ( X_New >= X_Old(nX) ) THEN

      IF ( Y_New <= Y_Old(1) ) THEN

        ! Extrapolate up along X
        Z_NewY1 = Z_Old(nX,1) + ( X_new       - X_Old(nX)    ) &
                              * ( Z_Old(nX,1) - Z_Old(nX-1,1)) &
                              / ( X_Old(nX)   - X_Old(nX-1)  )

        Z_NewY2 = Z_Old(nX,2) + ( X_new       - X_Old(nX)    ) &
                              * ( Z_Old(nX,2) - Z_Old(nX-1,2)) &
                              / ( X_Old(nX)   - X_Old(nX-1)  )

        ! Extrapolate down along Y
        Z_New =  Z_NewY1 - ( Y_Old(1) - Y_New   ) &
                         * ( Z_NewY2  - Z_NewY1 ) &
                         / ( Y_Old(2) - Y_Old(1) )

      ELSE IF ( Y_New >= Y_Old(nY) ) THEN

        ! Extrapolate up along X
        Z_NewY1 = Z_Old(nX,nY-1) + ( X_new          - X_Old(nX)       ) &
                                 * ( Z_Old(nX,nY-1) - Z_Old(nX-1,nY-1)) &
                                 / ( X_Old(nX)      - X_Old(nX-1)     )

        Z_NewY2 = Z_Old(nX,nY) + ( X_new        - X_Old(nX)     ) &
                               * ( Z_Old(nX,nY) - Z_Old(nX-1,nY)) &
                               / ( X_Old(nX)    - X_Old(nX-1)   )

        ! Extrapolate up along Y
        Z_New =  Z_NewY2 + ( Y_New     - Y_Old(nY)   ) &
                         * ( Z_NewY2   - Z_NewY1     ) &
                         / ( Y_Old(nY) - Y_Old(nY-1) )


      ELSE 

        DO iY = 1,nY

          IF ( Y_Old(iY) >= Y_New ) THEN

            ! Extrapolate up along X
            Z_NewY1 = Z_Old(nX,iY-1) + ( X_new          - X_Old(nX)       ) &
                                     * ( Z_Old(nX,iY-1) - Z_Old(nX-1,iY-1)) &
                                     / ( X_Old(nX)      - X_Old(nX-1)     )

            Z_NewY2 = Z_Old(nX,iY) + ( X_new        - X_Old(nX)     ) &
                                   * ( Z_Old(nX,iY) - Z_Old(nX-1,iY)) &
                                   / ( X_Old(nX)    - X_Old(nX-1)   )
            
            ! Interpolate along Y
            Z_New =  Z_NewY1 + (Y_New     - Y_Old(iY-1)) &
                             * (Z_NewY2   - Z_NewY1    ) &
                             / (Y_Old(iY) - Y_Old(iY-1))

            RETURN

          END IF

        END DO

      END IF

    ELSE IF ( Y_New <= Y_Old(1) ) THEN
        
      DO iX = 1,nX

        IF ( X_Old(iX) >= X_New ) THEN

          ! Interpolate along x
          Z_NewY1 = Z_Old(iX-1,1) + (X_New       - X_Old(iX-1)     ) &
                                  * (Z_Old(iX,1) - Z_Old(iX-1,1)) &
                                  / (X_Old(iX)   - X_Old(iX-1)     )
  
          Z_NewY2 = Z_Old(iX-1,2) + (X_New       - X_Old(iX-1)   ) &
                                  * (Z_Old(iX,2) - Z_Old(iX-1,2)) &
                                  / (X_Old(iX)   - X_Old(iX-1)   )

          ! Extrapolate down along y
          Z_New =  Z_NewY1 - (Y_Old(1) - Y_New   ) &
                           * (Z_NewY2  - Z_NewY1 ) &
                           / (Y_Old(2) - Y_Old(1))
  
  
          RETURN
  
        END IF

      END DO
    
    ELSE IF ( Y_New >= Y_Old(nY) ) THEN

      DO iX = 1,nX

        IF ( X_Old(iX) >= X_New ) THEN

          ! Interpolate along x
          Z_NewY1 = Z_Old(iX-1,nY-1) + (X_New          - X_Old(iX-1)     ) &       
                                     * (Z_Old(iX,nY-1) - Z_Old(iX-1,nY-1)) &          
                                     / (X_Old(iX)      - X_Old(iX-1)     )      
 
          Z_NewY2 = Z_Old(iX-1,nY) + (X_New        - X_Old(iX-1)   ) &    
                                   * (Z_Old(iX,nY) - Z_Old(iX-1,nY)) &      
                                   / (X_Old(iX)    - X_Old(iX-1)   )   
 
          ! Extrapolate up along y                     
          Z_New =  Z_NewY2 + (Y_New     - Y_Old(nY)  ) &     
                           * (Z_NewY2   - Z_NewY1    ) &  
                           / (Y_Old(nY) - Y_Old(nY-1))    
 
 
          RETURN

        END IF

      END DO

    ELSE
    
      DO iX = 1,nX
        DO iY = 1,nY
          
          IF ( X_Old(iX) >= X_New .AND. Y_Old(iY) >= Y_New ) THEN
            
            ! Interpolate along X
            Z_NewY1 = Z_Old(iX-1,iY-1) + (X_New          - X_Old(iX-1)     ) &
                                       * (Z_Old(iX,iY-1) - Z_Old(iX-1,iY-1)) &
                                       / (X_Old(iX)      - X_Old(iX-1)     )
            
            Z_NewY2 = Z_Old(iX-1,iY)   + (X_New        - X_Old(iX-1)   ) &
                                       * (Z_Old(iX,iY) - Z_Old(iX-1,iY)) &
                                       / (X_Old(iX)    - X_Old(iX-1)   )
    
            ! Interpolate along Y
            Z_New =  Z_NewY1 + (Y_New     - Y_Old(iY-1)) &
                             * (Z_NewY2   - Z_NewY1    ) &
                             / (Y_Old(iY) - Y_Old(iY-1))
                                   
            
            RETURN

          END IF

        END DO
      END DO

    END IF

  END SUBROUTINE LinearInterpolation2D


  SUBROUTINE QuadraticInterpolation1D(X_Old, Y_Old, nX, &
                                   X_New, Y_new)

    REAL(DP), INTENT(OUT)   :: Y_New
    REAL(DP), INTENT(IN)    :: X_New
    REAL(DP), INTENT(IN)    :: X_Old(nX)
    REAL(DP), INTENT(IN)    :: Y_Old(nX)
    INTEGER, INTENT(IN)     :: nX

    INTEGER  :: i
    REAL(DP) :: L2_0, L2_1, L2_2

    IF (X_New <= X_Old(1)) THEN

      !Extrapolate down
      Y_New = Y_Old(1) - (X_Old(1) - X_New) * (Y_Old(2) - Y_Old(1)) &
                                            / (X_Old(2) - X_Old(1))

      RETURN

    ELSE IF ( X_New <= X_Old(2) ) THEN

      Y_New = Y_Old(1) + (X_New    - X_Old(1)) &
                       * (Y_Old(2) - Y_Old(1)) &
                       / (X_Old(2) - X_Old(1))
    
    ELSE IF (X_New >= X_Old(nX)) THEN

      !Exptrapolate up
      Y_New = Y_Old(nX) + (X_New - X_Old(nX)) * (Y_Old(nX) - Y_Old(nX-1)) &
                                              / (X_OLD(nX) - X_Old(nX-1))
      RETURN

    ELSE IF ( X_New >= X_Old(nX-1) ) THEN

      Y_New = Y_Old(nX-1) + (X_New     - X_Old(nX-1)) &
                          * (Y_Old(nX) - Y_Old(nX-1)) &
                          / (X_Old(nX) - X_Old(nX-1))

    ELSE
  
      DO i = 1,nX
    
        IF( X_Old(i) >= X_New ) THEN
        
          L2_0 = ( (X_New      - X_Old(i)) * (X_New      - X_Old(i+1)) ) / &
                 ( (X_Old(i-1) - X_Old(i)) * (X_Old(i-1) - X_Old(i+1)) )
      
          L2_1 = ( (X_New    - X_Old(i-1)) * (X_New    - X_Old(i+1)) ) / &
                 ( (X_Old(i) - X_Old(i-1)) * (X_Old(i) - X_Old(i+1)) )
          
          L2_2 = ( (X_New      - X_Old(i-1)) * (X_New      - X_Old(i)) ) / &
                 ( (X_Old(i+1) - X_Old(i-1)) * (X_Old(i+1) - X_Old(i)) )  
          
          Y_New = Y_Old(i-1)*L2_0 + Y_Old(i)*L2_1 + Y_Old(i+1)*L2_2
          
          RETURN

        END IF

      END DO

    END IF

  END SUBROUTINE QuadraticInterpolation1D



END MODULE InterpolationModule
