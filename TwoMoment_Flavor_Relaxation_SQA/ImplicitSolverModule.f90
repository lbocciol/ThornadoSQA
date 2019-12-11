MODULE ImplicitSolverModule 

  USE KindModule, ONLY: &
    DP, Zero, Half, One
  USE UnitsModule, ONLY: &
    MeV, Centimeter, Second, &
    BoltzmannConstant
  USE ProgramHeaderModule, ONLY: &
    nZ, nX, nDOF, nDOFX, &
    nDOFE, nNodesE, &
    iZ_E0, iZ_B0, iE_B0, &
    iE_E0, iX_E0, iX_B0
  USE ReferenceElementModule, ONLY: &
    NodeNumberTable
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE UtilitiesModule, ONLY: &
    NodeNumber
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    ComputeNeutrinoOpacities_EC_Points
  USE RadiationFieldsModule, ONLY: &
    nSpecies, nCR, uCR, &
    iCR_N, iCR_G1, iCR_G2, iCR_G3, &
    iNuE, iNuE_Bar, iNuX, iNuX_Bar
  USE MeshModule, ONLY: &
    MeshE, NodeCoordinate
  USE FluidFieldsModule, ONLY: &                          
    uPF, uAF, nPF, nAF, &
    iPF_D, iAF_T, iAF_Ye, &                             
    iAF_Me, iAF_Mn, iAF_Mp
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    FermiDirac
  USE InitializationModule, ONLY: &
    nX_G, nE_G, Energies,&
    SigmaOsc, ChiOsc, EtaOsc

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: EvolveImplicit 
  PUBLIC :: ComputeEmission

  REAL(DP), ALLOCATABLE :: Chi_Temp(:,:,:), Chi(:,:,:)
  REAL(DP), ALLOCATABLE :: Eta(:,:,:)
  REAL(DP), ALLOCATABLE :: PF_N(:,:) , AF_N(:,:) 
  REAL(DP), ALLOCATABLE :: CR_New(:,:,:,:), CR_old(:,:,:,:)
  REAL(DP), ALLOCATABLE :: CR_check(:,:,:,:)

CONTAINS

  SUBROUTINE EvolveImplicit &
          ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, dt, U_R, dU_R )

    REAL(DP), INTENT(IN) :: dt
    INTEGER,  INTENT(in)    :: &
      iZ_B0(4), iZ_E0(4), iZ_B1(4), iZ_E1(4)
    REAL(DP), INTENT(in)    :: &
      U_R (1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)
    REAL(DP), INTENT(out) :: &
      dU_R(1:nDOF ,iZ_B1(1):iZ_E1(1),iZ_B1(2):iZ_E1(2),iZ_B1(3):iZ_E1(3),iZ_B1(4):iZ_E1(4),1:nCR,1:nSpecies)

    INTEGER  :: iX1, iX2, iX3, iPF, iAF, iCR, iS, iE, iN_E, iN_X
    INTEGER  :: iNode, iNodeX, iNodeE, iNodeX1, iNodeX2, iNodeX3

    LOGICAL  :: Recycle
    REAL(DP) :: Error

    CALL InitializeArrays
    
    !$OMP PARALLEL DO PRIVATE( iPF, iAF, iN_X, iX1, iX2, iX3 )
    DO iN_X = 1, nX_G

      DO iPF = 1, nPF
        iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
        iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
        iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
        iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

        PF_N(iN_X,iPF) = uPF(iNodeX,iX1,iX2,iX3,iPF)
      END DO

      DO iAF = 1, nAF
        iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
        iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
        iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
        iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

        AF_N(iN_X,iAF) = uAF(iNodeX,iX1,iX2,iX3,iAF)
      END DO

    END DO
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DEFAULT(SHARED) &
    !$OMP PRIVATE( iE,iNodeE,iNodeX,iNode,iX1,iX2,iX3 )
    !$OMP DO COLLAPSE(4)
    DO iS = 1, nSpecies
      DO iCR = 1, nCR
        DO iN_X = 1, nX_G
          DO iN_E = 1, nE_G

            iE     = MOD( (iN_E-1) / nDOFE, nZ(1) ) + iE_B0
            iNodeE = MOD( (iN_E-1)        , nDOFE ) + 1

            iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
            iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
            iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
            iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

            iNode  = iNodeE &
                     + ( iNodeX - 1 ) * nDOFE

            CR_Old(iN_E,iN_X,iS,iCR) = U_R(iNode,iE,iX1,iX2,iX3,iCR,iS)

          END DO
        END DO
      END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
    
    DO iS = 1,2 

      CALL ComputeNeutrinoOpacities_EC_Points &   
             ( 1, nE_G, 1, nX_G, &                
               Energies (:), &                         
               PF_N(:,iPF_D ), &                  
               AF_N(:,iAF_T ), &                  
               AF_N(:,iAF_Ye), &                  
               iS, Chi_Temp(:,:,iS) )                  

    END DO

    CALL ComputeEmission(  &
        Energies(:), Chi_Temp(:,:,:), &
        Chi(:,:,:), Eta(:,:,:), nE_G, nX_G )
   
    !Calculate Contribution from oscillations
    DO iN_X = 1, nX_G
      iN_E = 0
      DO iE = iE_B0,iE_E0
        DO iNodeE = 1,nNodesE

          iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
          iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
          iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
          iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1
          iNodeX1 = NodeNumberTableX(1,iNodeX)
          iNode = NodeNumber(iNodeE, iNodeX1, 1, 1 )
          iN_E = iN_E + 1

          DO iCR = 1,nCR

              EtaOsc(iN_E,iN_X,iCR,iNuE) = SigmaOsc(iN_E,iN_X,1) * uCR(iNode,iE,iX1,iX2,iX3,iCR,iNuX)
              ChiOsc(iN_E,iN_X,iNuE)     = SigmaOsc(iN_E,iN_X,1)

              EtaOsc(iN_E,iN_X,iCR,iNuX) = SigmaOsc(iN_E,iN_X,1) * uCR(iNode,iE,iX1,iX2,iX3,iCR,iNuE)
              ChiOsc(iN_E,iN_X,iNuX)     = SigmaOsc(iN_E,iN_X,1)

              EtaOsc(iN_E,iN_X,iCR,iNuE_Bar) = SigmaOsc(iN_E,iN_X,2) * uCR(iNode,iE,iX1,iX2,iX3,iCR,iNuX_Bar)
              ChiOsc(iN_E,iN_X,iNuE_Bar)     = SigmaOsc(iN_E,iN_X,2)

              EtaOsc(iN_E,iN_X,iCR,iNuX_Bar) = SigmaOsc(iN_E,iN_X,2) * uCR(iNode,iE,iX1,iX2,iX3,iCR,iNuE_Bar)
              ChiOsc(iN_E,iN_X,iNuX_Bar)     = SigmaOsc(iN_E,iN_X,2)

          END DO

        END DO
      END DO
    END DO


    CALL SolveThisIteration( dt, CR_New)

    !$OMP PARALLEL DEFAULT(SHARED) &
    !$OMP PRIVATE( iN_E,iNodeE,iNodeX,iN_X )
    !$OMP DO COLLAPSE(7)
    DO iS = 1, nSpecies
      DO iCR = 1, nCR
        DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
        DO iX1 = iX_B0(1), iX_E0(1)
          DO iE = iE_B0, iE_E0
            DO iNode = 1, nDOF

              iNodeX = MOD( (iNode-1) / nNodesE, nDOFX   ) + 1
              iNodeE = MOD( (iNode-1)          , nNodesE ) + 1

              iN_X = iNodeX &
                     + ( iX1 - iX_B0(1) ) * nDOFX &
                     + ( iX2 - iX_B0(2) ) * nDOFX * nX(1) &
                     + ( iX3 - iX_B0(3) ) * nDOFX * nX(1) * nX(2)
              iN_E = iNodeE &
                     + ( iE  - iE_B0    ) * nDOFE

              dU_R(iNode,iE,iX1,iX2,iX3,iCR,iS) = &
                  ( CR_New(iN_E,iN_X,iS,iCR) - CR_Old(iN_E,iN_X,iS,iCR) ) / dt

            END DO
          END DO
        END DO
        END DO
        END DO
      END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
    
    CALL FinalizeArrays

  END SUBROUTINE EvolveImplicit

  SUBROUTINE SolveThisIteration( dt, CR_New_Local )
    
    REAL(DP), INTENT(IN)  :: dt
    REAL(DP), INTENT(OUT) :: CR_New_Local(nE_G,nX_G,nSpecies,nCR)

    INTEGER :: iS, iN_E, iN_X
   
    !$OMP PARALLEL DEFAULT(SHARED)
    !$OMP DO COLLAPSE(3)
    DO iN_E = 1,nE_G
      DO iN_X = 1,nX_G
        DO iS = 1,nSpecies
    
            ! Add contribution to Chi from Oscillations
            Chi(iN_E,iN_X,iS) = Chi(iN_E,iN_X,iS) + ChiOsc(iN_E,iN_X,iS)
            
            CR_New_Local(iN_E,iN_X,iS,iCR_N) = &
                ( CR_Old(iN_E,iN_X,iS,iCR_N) + dt * Eta(iN_E,iN_X,iS) & 
                + dt * EtaOsc(iN_E,iN_X,iCR_N,iS) ) &
                / ( One + Chi(iN_E,iN_X,iS) * dt )

            CR_New_Local(iN_E,iN_X,iS,iCR_G1) = &
                (   CR_Old(iN_E,iN_X,iS,iCR_G1) &
                + dt * EtaOsc(iN_E,iN_X,iCR_G1,iS) ) &
                / ( One + Chi(iN_E,iN_X,iS) * dt )

            CR_New_Local(iN_E,iN_X,iS,iCR_G2) = &
                (   CR_Old(iN_E,iN_X,iS,iCR_G2) &
                + dt * EtaOsc(iN_E,iN_X,iCR_G2,iS) ) &
                / ( One + Chi(iN_E,iN_X,iS) * dt )

            CR_New_Local(iN_E,iN_X,iS,iCR_G3) = &
                (   CR_Old(iN_E,iN_X,iS,iCR_G3) &
                + dt * EtaOsc(iN_E,iN_X,iCR_G3,iS) ) &
                / ( One + Chi(iN_E,iN_X,iS) * dt )

        END DO
      END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL

  END SUBROUTINE SolveThisIteration

  SUBROUTINE ComputeEmission(E_N, Chi_In, &
                     Chi, Eta, nE_G, nX_G )

    INTEGER,  INTENT(IN)  :: nE_G, nX_G
    REAL(DP), INTENT(IN)  :: E_N(nE_G)
    REAL(DP), INTENT(IN)  :: Chi_In(nE_G,nX_G,2)
    REAL(DP), INTENT(OUT) :: Chi(nE_G,nX_G,nSpecies)
    REAL(DP), INTENT(OUT) :: Eta(nE_G,nX_G,nSpecies)

    INTEGER  :: iS, iN_E, iN_X
    INTEGER  :: iX1, iX2, iX3, iNodeX
    REAL(DP) :: Mnu, kT, FD

    !$OMP PARALLEL DO DEFAULT( SHARED ) &
    !$OMP PRIVATE( iNodeX,iX1,iX2,iX3,iS,Mnu,kT,iN_E,FD )
    DO iN_X = 1, nX_G

      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

      DO iS = 1, nSpecies

        IF ( iS .EQ. iNuE ) THEN

          Mnu = + ( uAF(iNodeX,iX1,iX2,iX3,iAF_Me) &
                + uAF(iNodeX,iX1,iX2,iX3,iAF_Mp) &
                - uAF(iNodeX,iX1,iX2,iX3,iAF_Mn) )

        ELSE IF ( iS .EQ. iNuE_Bar ) THEN

          Mnu = - ( uAF(iNodeX,iX1,iX2,iX3,iAF_Me) &
                      + uAF(iNodeX,iX1,iX2,iX3,iAF_Mp) &
                      - uAF(iNodeX,iX1,iX2,iX3,iAF_Mn) )
        ELSE

          Mnu = Zero

        END IF

        kT = BoltzmannConstant &
             * uAF(iNodeX,iX1,iX2,iX3,iAF_T)

        DO iN_E = 1, nE_G
            
          FD = FermiDirac( E_N(iN_E), Mnu, kT )

          IF ( iS == iNuE ) THEN
          
            Eta(iN_E,iN_X,iS) = Chi_In(iN_E,iN_X,1) * FD
            Chi(iN_E,iN_X,iS) = Chi_In(iN_E,iN_X,1)

          ELSE IF ( iS == iNuE_Bar ) THEN

            Eta(iN_E,iN_X,iS) = Chi_In(iN_E,iN_X,2) * FD
            Chi(iN_E,iN_X,iS) = Chi_In(iN_E,iN_X,2)

          ELSE 

            Eta(iN_E,iN_X,iS) = Zero
            Chi(iN_E,iN_X,iS) = Zero
          
          END IF

        END DO

      END DO

    END DO
    !$OMP END PARALLEL DO

  END SUBROUTINE ComputeEmission

  SUBROUTINE InitializeArrays

    ALLOCATE( PF_N (nX_G,nPF) )
    ALLOCATE( AF_N (nX_G,nAF) )
    ALLOCATE( Chi_Temp(nE_G,nX_G,2) )
    ALLOCATE( Chi  (nE_G,nX_G,nSpecies) )
    ALLOCATE( Eta  (nE_G,nX_G,nSpecies) )
 
    ALLOCATE( CR_New  (nE_G,nX_G,nSpecies,nCR) )
    ALLOCATE( CR_Old  (nE_G,nX_G,nSpecies,nCR) )
    ALLOCATE( CR_check(nE_G,nX_G,nSpecies,nCR) )
    
  END SUBROUTINE InitializeArrays

  SUBROUTINE FinalizeArrays
   
    INTEGER :: iN_X, iN_E
    INTEGER :: iX1, iX2, iX3, iNodeX
    
    DEALLOCATE( PF_N, AF_N,  &
                Chi_Temp, Chi, Eta,       &
                CR_New, CR_old, CR_check )
    
  END SUBROUTINE

END MODULE ImplicitSolverModule 
