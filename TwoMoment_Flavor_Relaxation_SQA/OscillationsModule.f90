MODULE OscillationsModule
  
  USE KindModule, ONLY: &
    DP, One, Zero, Five, &
    Two, Pi, TwoPi
  USE UnitsModule, ONLY: &
    MeV, Gram, Centimeter, &
    Second
  USE ProgramHeaderModule, ONLY: &
    nNodesE, &
    iE_B0, iE_E0
  USE MeshModule, ONLY: &
    MeshE
  USE InitializationModule, ONLY: &
    BB5, BB6, AA, nRK, nRKOrder
  USE InitializationModule, ONLY: &
    YIdentity, Hvf, &
    Enu, &
    fMatrixOsc, &
    SMatrixOsc, &
    nF, nM, nY, nS, nE_G
  USE OscillationsUtilsModule, ONLY: &
    Im, B, W, JInverse, CSI, &
    Eigenvalues, EigenvectorMatrix

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: EvolveOscillations
  PUBLIC :: Diagonalize
  PUBLIC :: InitializeOscArrays
  PUBLIC :: FinalizeOscArrays

  ! --- HAMILTONIANS & CO. --- !
  COMPLEX(DP), ALLOCATABLE :: Hmsw(:,:,:,:) !MSW Hamiltonian in flavor basis
  COMPLEX(DP), ALLOCATABLE :: Hf(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: U0(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: pm0(:,:,:,:)

  REAL(DP), ALLOCATABLE :: k0(:,:,:)

  ! --- CONSTANTS USED IN OSCILLATION CALCULATIOnS --- ! 
  REAL(DP), PARAMETER :: eV_to_erg = 1.60218d-12 !erg
  REAL(DP), PARAMETER :: hbar = 1.05457266d-27 ! erg s
  REAL(DP), PARAMETER :: clite = 2.99792458d10 !cm/s
  REAL(DP), PARAMETER :: hbarc = hbar*clite
  REAL(DP), PARAMETER :: GF = 1.1663787d-5/ (1e9*1e9*eV_to_erg**2) &
                            * hbarc**3 !erg cm^3
  REAL(DP), PARAMETER :: Mp = 1.6726219e-24 ! g
  
CONTAINS

  SUBROUTINE EvolveOscillations( R, t, dt )

    REAL(DP), INTENT(INOUT) :: t, dt, R

    LOGICAL     :: Reloop
    REAL(DP)    :: Yerror, MaxError, TimeTmp
    REAL(DP)    :: Ks(nRK,nM,nE_G,nS,nY)
    REAL(DP)    :: Ytmp(nM,nE_G,nS,nY), Y(nM,nE_G,nS,nY)
    COMPLEX(DP) :: SMSW(nF,nF), SSI(nF,nF), SNew(nF,nF)
    
    INTEGER :: k,l,m,iN_E,iY,x

    REAL(DP) :: Accuracy = 1d-12
    REAL(DP) :: Increase = 1.1

    Ks(:,:,:,:,:) = 0.0d0

    DO m = 1,nM
      DO iN_E = 1,nE_G

        Y(m,iN_E,:,:) = YIdentity

      END DO
    END DO

    CALL Get_P
    
    Reloop = .true.

    DO WHILE(Reloop .eqv. .true.)

      DO k = 1,nRK

        Ytmp = Y
        
        DO l = 1,k-1
          !$OMP PARALLEL DO COLLAPSE(4) &
          !$OMP DEFAULT(SHARED) PRIVATE( m,iN_E,x,iY )
          DO m = 1,nM
            DO iN_E = 1,nE_G
              DO x = 1,nS
                DO iY = 1,nY

                  Ytmp(m,iN_E,x,iY) = Ytmp(m,iN_E,x,iY) + AA(k,l) * Ks(l,m,iN_E,x,iY)

                END DO
              END DO
            END DO
          END DO
          !$OMP END PARALLEL DO
        END DO

        CALL RK_step( Ytmp, Ks(k,:,:,:,:), pm0, dt, R )
        
      END DO
      
      ! Do the rest of the RK and evaluate error
      MaxError = 0.0d0
      Ytmp = Y

      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(m,iN_E,x,iY,Yerror)
      !$OMP DO COLLAPSE(2) REDUCTION(Max:MaxError)
      DO m = 1,nM
        DO iN_E = 1,nE_G
          DO x=1,nS
            DO iY=1,nY

              Yerror = 0.0d0

              DO k=1,nRK

                Ytmp(m,iN_E,x,iY) = Ytmp(m,iN_E,x,iY) + BB5(k) * Ks(k,m,iN_E,x,iY)
                Yerror = Yerror + (BB5(k) - BB6(k)) * Ks(k,m,iN_E,x,iY)

              END DO

              MaxError = MAX(MaxError,ABS(Yerror))

            END DO
          END DO
        END DO
      END DO
      !$OMP END DO
      !$OMP END PARALLEL

      TimeTmp = t + dt
      
      IF(MaxError .gt. Accuracy) THEN

        dt = dt * 0.9 *( (Accuracy/MaxError)**(1.0d0/(nRKOrder-1.0d0)) )
        Reloop=.true.
      
      ELSE

        dt = dt * Increase
        Reloop = .false.

        IF(MaxError .gt. 0.0d0) &
          dt = dt * MIN( 1.0, ((Accuracy/MaxError)**(1.0d0/5.0d0))/Increase)
          Y = Ytmp

      END IF
    END DO
    
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(m,iN_E,SMSW,SSI,SNew)
    !$OMP DO COLLAPSE(2)
    DO m = 1,nM
      DO iN_E = 1,nE_G

        SMSW = MATMUL(W(Y(m,iN_E,1,:)),B(Y(m,iN_E,1,:)))
        SSI  = MATMUL(W(Y(m,iN_E,2,:)),B(Y(m,iN_E,2,:)))
    
        SNew = MATMUL(SMSW,SSI)
        SMatrixOsc(m,iN_E,:,:) = MATMUL(MATMUL(U0(m,iN_E,:,:),SNew), &
                      CONJG(TRANSPOSE(U0(m,iN_E,:,:))))
        fMatrixOsc(m,iN_E,:,:) =  &
            MATMUL(MATMUL(SMatrixOsc(m,iN_E,:,:),fMatrixOsc(m,iN_E,:,:)), &
                      CONJG(TRANSPOSE((SMatrixOsc(m,iN_E,:,:)))))
        Y(m,iN_E,:,:) = YIdentity !not sure it's necessary

      END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
    
    t = TimeTmp

  END SUBROUTINE EvolveOscillations

  SUBROUTINE RK_step( Y ,K, pm0, dt, R)

    REAL(DP),    INTENT(IN)    :: dt
    REAL(DP),    INTENT(IN)    :: Y(nM,nE_G,nS,nY)
    COMPLEX(DP), INTENT(IN)    :: pm0(nM,nE_G,NF,NF) 
    REAL(DP),    INTENT(INOUT) :: R
    REAL(DP),    INTENT(INOUT) :: K(nM,nE_G,nS,nY)
    ! pm0 is the integrand of the SI potential in the mass basis
    
    COMPLEX(DP) :: VfSI(nM,NF,NF)
    COMPLEX(DP) :: VfSIE(NF,NF) !it's the integrand of the SI potential in
    !the flavor basis 
    COMPLEX(DP) :: Sfm(NF,NF) !it transforms S to f to calculate the potential,
    !and THEN it transforms again with U0 to go to the flavor basis
    COMPLEX(DP) :: UWBW(nM,nE_G,NF,NF)
    COMPLEX(DP) :: BSI(nM,nE_G,NF,NF) !S_{si} = W(Y_{si}) * B(Y_{si}), so BSI is
    !one of the two contributions to S_{si}
    INTEGER :: m,iN_E,iY,l

    !Now variables that I'm not sure what they are
    COMPLEX(DP) :: Ha(NF,NF)
    COMPLEX(DP) :: HB0,HB1
    REAL(DP) :: dvdr(4)
    REAL(DP) :: JI(3,4)

    VfSI(:,:,:) = Zero
   
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(m,iN_E,Sfm,VfSIE)
    !$OMP DO COLLAPSE(2)
    DO m = 1,nM
      
      DO iN_E = 1,nE_G
        
        !setup transformation 
        UWBW(m,iN_E,:,:) = MATMUL(MATMUL(U0(m,iN_E,:,:)  , W(Y(m,iN_E,1,:))), &
                               MATMUL(B(Y(m,iN_E,1,:)), W(Y(m,iN_E,2,:))))

        !For 1=msw only 4 and 5 are relevant (i.3. B(Y(msw)) = I)
        DO iY = 1,4
        
          K(m,iN_E,1,iY) = 0.0d0
        
        END DO
       
        K(m,iN_E,1,5) = k0(m,iN_E,1) * dt / hbar / TwoPi
        K(m,iN_E,1,6) = k0(m,iN_E,2) * dt / hbar / TwoPi
        !DONE WITH MSW!

        BSI(m,iN_E,:,:) = B(Y(m,iN_E,2,:))
        
        !Evaluate Si potential
        Sfm = MATMUL(UWBW(m,iN_E,:,:),BSI(m,iN_E,:,:))
        VfSIE = MATMUL(MATMUL(Sfm,pm0(m,iN_E,:,:)),CONJG(TRANSPOSE(Sfm)))

        !Below is the "MINus" part of eq. 9 from Richers 2019
        IF(m .eq. 2) VfSIE = - CONJG(VfSIE)
        
        !$OMP CRITICAL        
        VfSI(1,:,:) = VfSI(1,:,:) + VfSIE(:,:)
        !$OMP END CRITICAL
      
      END DO
    
    END DO
    !$OMP END DO
    !$OMP END PARALLEL

    VfSI(2,:,:) = - CONJG(VfSI(1,:,:))
    VfSI = VfSI * CSI( R )
    !VfSI(:,:,:) = 0.0d0 
    
    !Now the part that I DOn't REALly understand
    
    !$OMP PARALLEL DEFAULT(SHARED) &
    !$OMP PRIVATE(m,iN_E,Ha,HB0,HB1,dvdr,JI)
    !$OMP DO COLLAPSE(2) PRIVATE(m,iN_E,Ha,HB0,HB1,dvdr,JI)
    DO m = 1,nM
      
      DO iN_E = 1,nE_G
        
        Ha = MATMUL(MATMUL(CONJG(TRANSPOSE(UWBW(m,iN_E,:,:))),VfSI(m,:,:)),UWBW(m,iN_E,:,:))

        K(m,iN_E,2,5) = dt * REAL(Ha(1,1)) / hbar / TwoPi
        K(m,iN_E,2,6) = dt * REAL(Ha(2,2)) / hbar / TwoPi

        HB0 = -Im/hbar*( Ha(1,2)*BSI(m,iN_E,2,1) )
        HB1 = -Im/hbar*( Ha(1,2)*BSI(m,iN_E,2,2) )

        dvdr(1) = REAL(HB1)
        dvdr(2) = AIMAG(HB1)
        dvdr(3) = REAL(HB0)
        dvdr(4) = AIMAG(HB0)

        JI = JInverse(Y(m,iN_E,2,:))

        DO iY = 1,3
          
          K(m,iN_E,2,iY) = Zero
          
          DO l = iY,4

            K(m,iN_E,2,iY) = K(m,iN_E,2,iY) + JI(iY,l)*dvdr(l)
          
          END DO
          
          K(m,iN_E,2,iY) = K(m,iN_E,2,iY) * dt
        
        END DO

        K(m,iN_E,2,4) = Zero
      
      END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
 
  END SUBROUTINE RK_step

  SUBROUTINE Get_P

    REAL(DP)    :: nuHz, dE, dnuHz
    INTEGER     :: m, iN_E, iE, iNodeE

    DO m = 1,nM 
      
      iN_E = 0

      DO iE = iE_B0,iE_E0
        DO iNodeE = 1,nNodesE
          
          iN_E = iN_E + 1

          dE = MeshE % Width(iE) / nNodesE / MeV
          dnuHz = dE * 1e6 * eV_to_erg / (TwoPi*hbar)
          
          nuHz = Enu(iN_E) * 1e6 * eV_to_erg / (TwoPi*hbar)
          pm0(m,iN_E,:,:) = &
              MATMUL(MATMUL(CONJG(TRANSPOSE(U0(m,iN_E,:,:))), &
              fMatrixosc(m,iN_E,:,:)), &
              U0(m,iN_E,:,:)) * &
              SQRT(2.0d0)*GF*4.0d0*pi*nuHz**2*dnuHz/(clite**3)
         
          END DO
        END DO

    END DO

  END SUBROUTINE Get_P

  SUBROUTINE Diagonalize( Rho, Ye )
    
    REAL(DP), INTENT(IN) :: Rho, Ye

    INTEGER :: m, iN_E

    Hmsw(:,:,:,:) = Zero                            
    Hmsw(1,:,1,1) = SQRT( Two ) * GF * Rho * Ye / Mp
    Hmsw(2,:,1,1) = - CONJG( Hmsw(1,:,1,1) )

    Hf(:,:,:,:) = Hvf(:,:,:,:) + Hmsw(:,:,:,:)

    DO m = 1,nM
      DO iN_E = 1,nE_G

         k0(m,iN_E,:)   = Eigenvalues       ( Hf(m,iN_E,:,:) )
         U0(m,iN_E,:,:) = EigenvectorMatrix( Hf(m,iN_E,:,:), k0(m,iN_E,:), iN_E )

         IF(m == 2) U0(m,iN_E,:,:) = CONJG( U0(m,iN_E,:,:) )

      END DO
    END DO


  END SUBROUTINE Diagonalize

  SUBROUTINE InitializeOscArrays

    ALLOCATE( Hmsw(nM,nE_G,nF,nF) )
    ALLOCATE( Hf  (nM,nE_G,nF,nF) )
    ALLOCATE( U0  (nM,nE_G,nF,nF) )
    ALLOCATE( pm0 (nM,nE_G,nF,nF) )
    ALLOCATE( k0  (nM,nE_G,nF)    )

  END SUBROUTINE InitializeOscArrays

  SUBROUTINE FinalizeOscArrays

    DEALLOCATE( Hmsw, Hf, U0, pm0, k0 )

  END SUBROUTINE FinalizeOscArrays

END MODULE OscillationsModule
