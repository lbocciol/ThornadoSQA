MODULE ThornadoSQAInterfaceModule

  USE KindModule, ONLY: &
    DP, One, Zero, Five, &
    Half, Two, Pi, TwoPi
  USE UnitsModule, ONLY: &
    MeV, Gram, Centimeter, &
    Second, Kilometer, Erg, &
    BoltzmannConstant, &
    PlanckConstant
  USE PhysicalConstantsModule, ONLY: &
    SpeedOfLightCGS
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nX, nDOFX, &
    nNodesX, nNodesE, &
    iX_B0, iX_E0, &
    iE_B0, iE_E0
  USE MeshModule, ONLY: &
    MeshE, MeshX, &
    NodeCoordinate
  USE ReferenceElementModuleX, ONLY: & 
    NodeNumberTableX
<<<<<<< HEAD
  USE ReferenceElementModuleE, ONLY: &
    WeightsE
=======
  USE UtilitiesModule, ONLY: &
    NodeNumber
>>>>>>> 3b00b0779dac05edbb208cd4135870aa318f3008
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    ComputeNeutrinoOpacities_EC_Points
  USE NeutrinoOpacitiesModule, ONLY: &
    unitsEC
  USE FluidFieldsModule, ONLY: &
    uAF, iAF_Ye, iAF_T, nAF, &
    uPF, iPF_D, nPF, &
    iAF_Me, iAF_Mp, iAF_Mn
  USE RadiationFieldsModule, ONLY: &
    nSpecies, uCR, nCR, iCR_N, &
    iCR_G1, iNuE, iNuE_Bar, &
    iNuX, iNuX_Bar
  USE InitializationModule, ONLY: &
    IdentityC, Energies,  &
    fMatrixOsc, Psi0_loc, &
    Psi1_loc, &
    SMatrixOsc, SigmaOsc, &
    nF, nM, nE_G, nX_G,   &
    R_Shock, Rnu
  USE OscillationsModule, ONLY: &
    EvolveOscillations, &
    Diagonalize, &
    InitializeOscArrays, &
    FinalizeOscarrays
  USE OscillationsUtilsModule, ONLY: &
    CSI
  USE IntegrationModule, ONLY: &
    Update_My_IMEX_PDARS
  USE ImplicitSolverModule, ONLY: &
    ComputeEmission
  USE InputOutputRelaxationModule, ONLY: &
    WritefOscillations

  IMPLICIT NONE 
  PRIVATE
    
  PUBLIC :: SQADriver
  PUBLIC :: IntegrationDriver
  PUBLIC :: Get_Profile
  PUBLIC :: GetShockRadius
 
  ! These are public in case you want to use python wrapping
  PUBLIC :: OscillationsDriver
  PUBLIC :: InitializefMatrixOsc
  PUBLIC :: ResetSMatrix
  PUBLIC :: CalculateOpacitiesOsc
  PUBLIC :: SetNonDiagonalEl
  PUBLIC :: OscillationsInterface
  PUBLIC :: InitializeOscInterface
  PUBLIC :: FinalizeOscInterface

  REAL(DP) :: dt_loc = 1.0d-13
  REAL(DP) :: MaxTDump = 1d-9

  COMPLEX(DP), ALLOCATABLE :: SMatrixTotal(:,:,:,:)

CONTAINS

  SUBROUTINE SQADriver

    LOGICAL :: DoOscillations
    LOGICAL :: InitializeFirstZone

    INTEGER  :: iN_X, iX1, iX2, iX3, iNodeX, iNodeX1

    REAL(DP) :: tmax
    REAL(DP) :: R,dr
    
    INTEGER  :: iPF, iAF, iS
    REAL(DP) :: Chi_Temp(nE_G, nX_G, nSpecies)
    REAL(DP) :: Eta(nE_G, nX_G, nSpecies)
    REAL(DP) :: Chi(nE_G, nX_G, nSpecies)
    REAL(DP) :: PF_N(nX_G,nPF), AF_N(nX_G,nAF)

    DoOscillations = .FALSE.
    InitializeFirstZone = .TRUE.

    dt_loc = 1.0d-13 !Seconds
    CALL ResetSmatrix
   
    iN_X = 0

    DO iX1 = iX_B0(1),iX_E0(1)
    DO iX2 = iX_B0(2),iX_E0(2)
    DO iX3 = iX_B0(3),iX_E0(3)
      
      DO iNodeX = 1,nDOFX
        
        iN_X = iN_X + 1
        iNodeX1 = NodeNumberTableX(1,iNodeX)

        IF ( iN_X == nX_G ) CONTINUE
 
        ! This only works with 2 nodes

        R = NodeCoordinate( MeshX(1), iX1, iNodeX1 ) - &
              MeshX(1) % Width(iX1) / 4.0_DP
        R = R / Kilometer

        ! FOR 1D ONLY. FIND FIRST ZONE WHERE YOU WANT TO 
        ! START YOUR CALCULATIONS (E.G. 10 KM OUTSIDE
        ! NEUTRINOSPHERE      
        IF ( R > Rnu + 10.0_DP .AND. .NOT. DoOscillations ) THEN
  
          DoOscillations = .TRUE.
  
        END IF
  
        IF ( R > R_Shock / Kilometer ) THEN
  
          DoOscillations = .FALSE.
  
        END IF
  
        IF ( DoOscillations ) THEN
  
          IF ( InitializeFirstZone ) THEN
            
            CALL InitializefMatrixOsc( iNodeX,iX1,iX2,iX3, SMatrixTotal(:,:,:,:), nM, nE_G, nF )
  
            InitializeFirstZone = .FALSE.
  
          ELSE
  
            ! iX1-1 because you initialize fMatrixOsc by getting the 
            ! spectrum coming from the previous zone
            CALL InitializefMatrixOsc( iNodeX,iX1-1,iX2,iX3, SMatrixTotal(:,:,:,:), nM, nE_G, nF )
  
          END IF
  
          dr = MeshX(1) % Width(iX1) / nNodesX(1) / Centimeter
  
          tmax = dr / SpeedOfLightCGS
          
          CALL OscillationsDriver( R, tmax, &
             uPF(iNodeX,iX1,iX2,iX3,iPF_D ), &
             uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) )
  
          CALL CalculateOpacitiesOsc( dr, iN_X, iNodeX1, iX1, iX2, iX3, &
                      ABS( SMatrixTotal(:,:,1,1) )**2, nM, nE_G )
          
          WRITE(*,'(A16,4i4)') 'Done with Zone:',iNodeX1,iX1,iX2,iX3
  
        END IF
       
      END DO
    END DO
    END DO
    END DO

    DEALLOCATE( SMatrixTotal )

  END SUBROUTINE SQADriver

  SUBROUTINE OscillationsInterface( R, tOsc, dt_osc, &
                    f, S, nMi, nEi, nFi )

    REAL(DP), INTENT(INOUT) :: tOsc, dt_osc, R
    
    INTEGER,     INTENT(IN)  :: nMi, nEi, nFi
    COMPLEX(KIND=8), INTENT(OUT) :: f(nMi,nEi,nFi,nFi)
    COMPLEX(KIND=8), INTENT(OUT) :: S(nMi,nEi,nFi,nFi)

    CALL EvolveOscillations( R, tOsc, dt_osc )

    !IF (tOsc <= MaxTDump ) CALL WritefOscillations( tOsc )
    
    f = fMatrixOsc
    S = SMatrixOsc

  END SUBROUTINE OscillationsInterface

  SUBROUTINE InitializeOscInterface( Rho_in, Ye_in )

    REAL(DP), INTENT(IN) :: Rho_in, Ye_in
    REAL(DP) :: Rho, Ye
    
    Rho = Rho_in / Gram * Centimeter**3
    Ye  = Ye_in

    CALL InitializeOscArrays

    CALL Diagonalize( Rho, Ye )

  END SUBROUTINE InitializeOscInterface

  SUBROUTINE FinalizeOscInterface

    CALL FinalizeOscArrays

  END SUBROUTINE FinalizeOscInterface


  SUBROUTINE OscillationsDriver( R, TimeBlockEnd, &
                        Rho_in, Ye_in )             

    REAL(DP), INTENT(INOUT) :: TimeBlockEnd, R
    REAL(DP), INTENT(IN)    :: Rho_in, Ye_in

    REAL(DP) :: Rho, Ye
    REAL(DP) :: TimeBlock, t_wrts
    REAL(DP) :: dt_wrts = 1.0d-8
    LOGICAL  :: wrts, BlockDone
    
    INTEGER :: m, iN_E
    CHARACTER(256) :: rmstring

    Rho      = Rho_in / Gram * Centimeter**3
    Ye       = Ye_in

    CALL InitializeOscArrays

    CALL Diagonalize( Rho, Ye )

    BlockDone = .FALSE.
    TimeBlock = Zero
    t_wrts = dt_wrts

    DO WHILE( BlockDone .EQV. .FALSE. )

      IF( TimeBlock + dt_loc >= TimeBlockEnd ) THEN

        dt_loc = TimeBlockEnd - TimeBlock
        BlockDone = .TRUE.

      END IF

      IF ( TimeBlock > t_wrts ) THEN

        t_wrts = t_wrts + dt_wrts
        WRITE(*,*) TimeBlockEnd, TimeBlock, dt_loc
        WRITE(*,*) 'loc',SMatrixOsc(1,1,1,1)
        WRITE(*,*) 'tot',SMatrixTotal(1,1,1,1)

      END IF

      CALL EvolveOscillations( R, TimeBlock, dt_loc )
  
      !IF (TimeBlock <= MaxTDump ) CALL WritefOscillations( TimeBlock )

      ! --- Get S Total
      DO m = 1,nM
        DO iN_E = 1, nE_G

          SmatrixTotal(m,iN_E,:,:) = &
              MatMul( SMatrixTotal(m,iN_E,:,:), &
                      SMatrixOsc(m,iN_E,:,:) )
        
        END DO
      END DO

    END DO

    CALL FinalizeOscArrays
    
  END SUBROUTINE OscillationsDriver                 

  SUBROUTINE IntegrationDriver( dt )

    REAL(DP), INTENT(INOUT) :: dt

    CALL Update_My_IMEX_PDARS &
       ( dt, uCR, &
         Explicit_Option = .TRUE., &
         Implicit_Option = .TRUE., &
         SingleStage_Option = .FALSE., &
         CallFromThornado_Option = .TRUE. )

  END SUBROUTINE IntegrationDriver

  SUBROUTINE InitializefMatrixOsc( iNodeX,iX1,iX2,iX3, STotal, nMi, nEi, nFi )

    INTEGER, INTENT(IN) :: iNodeX,iX1,iX2,iX3, nMi, nEi, nFi
    COMPLEX(KIND=8), INTENT(IN) :: STotal(nMi,nEi,nFi,nFi)

    INTEGER  :: iS, m, iE, iN_E
    INTEGER  :: iNodeE, iNode, iNodeX1

    REAL(DP) :: fPinched(nM,nE_G,nF,nF)
    REAL(DP) :: AvgE(nSpecies), AvgE2(nSpecies), Norm(nSpecies)
    REAL(DP) :: LumE(nSpecies), LumE2(nSpecies)
    REAL(DP) :: dE
    REAL(DP) :: alpha(nSpecies)
    REAL(DP) :: f_loc(nE_G,nSpecies)
    REAL(DP) :: R, Mnu, kT
    CHARACTER(3) :: FileNUmber
    REAL(DP) :: Dummy

    iNodeX1 = NodeNumberTableX(1,iNodeX)
    R = NodeCoordinate( MeshX(1), iX1, iNodeX1 ) / Kilometer 

    iN_E = 0

    DO iE = iE_B0, iE_E0

      DO iNodeE = 1,nNodesE
     
        iN_E = iN_E + 1
        
        iNode = NodeNumber(iNodeE, iNodeX, 1, 1 )
        
        Psi0_loc(iN_E,:) = uCR(iNode,iE,iX1,iX2,iX3,iCR_N,:)
        Psi1_loc(iN_E,:) = uCR(iNode,iE,iX1,iX2,iX3,iCR_G1,:)

      END DO

    END DO

    AvgE(:)  = Zero
    AvgE2(:) = Zero
    Norm(:)  = Zero
    LumE(:)  = Zero

    DO iS = 1,nSpecies
        
      iN_E = 0

      DO iE = iE_B0,iE_E0
    
        dE = MeshE % Width(iE)
        
        DO iNodeE = 1,nNodesE

          iN_E = iN_E + 1
               
          Norm(iS)  = Norm(iS)  + dE * WeightsE(iNodeE) &
              * Energies(iN_E)**2 * Psi0_loc(iN_E,iS)
          
          AvgE(iS)  = AvgE(iS)  + dE * WeightsE(iNodeE) &
              * Energies(iN_E)**3 * Psi0_loc(iN_E,iS)

          AvgE2(iS) = AvgE2(iS) + dE * WeightsE(iNodeE) &
              * Energies(iN_E)**4 * Psi0_loc(iN_E,iS)
          
          LumE(iS)  = LumE(iS)  + dE * WeightsE(iNodeE) &
              * Energies(iN_E)**3 * Psi1_loc(iN_E,iS) * &
              2.0_DP * Pi / Erg**4

        END DO
        
      END DO
    
    END DO

    AvgE  = AvgE  / Norm
    AvgE2 = AvgE2 / Norm
    LumE  = LumE * 4.0_DP * Pi * (R*1.0d5)**2 * SpeedOfLightCGS / &
            (PlanckConstant / (Erg*Second) * SpeedOfLightCGS)**3
    
    fPinched(:,:,:,:) = Zero

    DO iS = 1,nSpecies

      alpha(iS) = - ( Two - SQRT( AvgE2(iS) ) / AvgE(iS) ) / &
                    ( One - SQRT( AvgE2(iS) ) / AvgE(iS) )

      kT = BoltzmannConstant &
               * AvgE(iS) / 3.15d0

      kT = BoltzmannConstant &
            * uAF(iNodeX,iX1,iX2,iX3,iAF_T)
    
      IF( iS .EQ. iNuE )THEN

        Mnu = + ( uAF(iNodeX,iX1,iX2,iX3,iAF_Me) &
                    + uAF(iNodeX,iX1,iX2,iX3,iAF_Mp) &
                    - uAF(iNodeX,iX1,iX2,iX3,iAF_Mn) )

      ELSE IF( iS .EQ. iNuE_Bar )THEN

        Mnu = - ( uAF(iNodeX,iX1,iX2,iX3,iAF_Me) &
                    + uAF(iNodeX,iX1,iX2,iX3,iAF_Mp) &
                    - uAF(iNodeX,iX1,iX2,iX3,iAF_Mn) )

      ELSE

        Mnu = Zero

      END IF

      DO iN_E = 1,nE_G

        !This is the pinched spectrum from Keil 2003
        f_loc(iN_E,iS) = ( Energies(iN_E) / AvgE(iS) )**alpha(iS) * &
            EXP( - (alpha(iS) + One) * Energies(iN_E) / AvgE(iS) ) 
    
        !Now normalize this spectrum to 1
        f_loc(iN_E,iS) = f_loc(iN_E,iS) * ( alpha(iS)+1.0_DP )**( alpha(iS)+1.0_DP ) / &
            ( AvgE(iS) / Erg )/ GAMMA( alpha(iS)+1.0_DP )
        
        !Now multiply by Lum / AvgE / 4piR^2c as in Duan 2006 and get units of
        !1/erg/cm^3. Also, this has an extra factor of 2pi because we have integrated
        ! over the azimuthal angle
        f_loc(iN_E,iS) = f_loc(iN_E,iS) * LumE(iS) / ( AvgE(iS) / Erg ) / &
            ( 2.0_DP * Pi * SpeedOfLightCGS * ( Rnu*1.0d5 )**2 )

        !f_loc(iN_E,iS) =  ( One / ( One + EXP( ( Energies(iN_E) - Mnu ) / kT ) ) )

      END DO
     
      DO iN_E = 1,nE_G

        !f_loc(iN_E,iS) = f_loc(iN_E,iS) / Norm(iS)
        
        IF( iS == iNuE ) fPinched(1,iN_E,1,1) = f_loc(iN_E,iS)
        IF( iS == iNuX ) fPinched(1,iN_E,2,2) = f_loc(iN_E,iS)

        IF( iS == iNuE_Bar ) fPinched(2,iN_E,1,1) = f_loc(iN_E,iS)
        IF( iS == iNuX_Bar ) fPinched(2,iN_E,2,2) = f_loc(iN_E,iS)

      END DO

    END DO
    
!    DO iN_E = 1,nE_G
!        WRITE(FileNUmber,'(I3.3)') iN_E
!        OPEN(UNIT=666,FILE='pinch' // FileNumber, STATUS='unknown', &
!            FORM = 'formatted', POSITION='append')
!        OPEN(UNIT=777,FILE='psi' // FileNumber, STATUS='unknown', &
!            FORM = 'formatted', POSITION='append')
! 
!        WRITE(666,'(10ES22.11E3)') R, Energies(iN_E) / MeV, &
!                             fPinched(1,iN_E,1,1), &
!                             fPinched(2,iN_E,1,1), &
!                             fPinched(1,iN_E,2,2), &
!                             fPinched(2,iN_E,2,2), &
!                             alpha(1),alpha(2),alpha(3),alpha(4)
!        
!         WRITE(777,'(6ES22.11E3)') R, Energies(iN_E) / MeV, &
!                             Psi0_loc(iN_E,1), &
!                             Psi0_loc(iN_E,2), &
!                             Psi0_loc(iN_E,3), &
!                             Psi0_loc(iN_E,4)
!        CLOSE(666)
!        CLOSE(777)
!    END DO

    DO m = 1,nM
      DO iN_E = 1,nE_G

        !CALL SetNonDiagonalEl( &
        !    fPinched(m,iN_E,1,1), &
        !    fPinched(m,iN_E,2,2), &
        !    fPinched(m,iN_E,1,2), &
        !    Mixing_Optional = Zero )

        !CALL SetNonDiagonalEl( &
        !    fPinched(m,iN_E,1,1), &
        !    fPinched(m,iN_E,2,2), &
        !    fPinched(m,iN_E,2,1), &
        !    Mixing_Optional = Zero )
        fPinched(m,iN_E,1,2) = Zero
        fPinched(m,iN_E,2,1) = Zero

      END DO
    END DO

    DO m = 1,nM
      DO iN_E = 1,nE_G

        fMatrixOsc(m,iN_E,:,:) = MATMUL( MATMUL( STotal(m,iN_E,:,:), &
          fPinched(m,iN_E,:,:) ), CONJG(TRANSPOSE( STotal(m,iN_E,:,:) )) )
        
      END DO
    END DO

  END SUBROUTINE InitializefMatrixOsc


  SUBROUTINE ResetSMatrix

    INTEGER :: m, iN_E

    ALLOCATE( SMatrixTotal(nM,nE_G,nF,nF) )

    DO m = 1,nM
      DO iN_E = 1,nE_G

        SMatrixOsc(m,iN_E,:,:) = IdentityC
        SMatrixTotal(m,iN_E,:,:) = IdentityC

      END DO
    END DO

  END SUBROUTINE ResetSMatrix


  SUBROUTINE CalculateOpacitiesOsc( ZoneWidth, iN_X, iNodeX1, iX1, iX2, iX3, &
                    See2, nMi, nEi )

    INTEGER,  INTENT(IN) :: nMi, nEi
    REAL(DP), INTENT(IN) :: See2(nMi, nEi)
    
    REAL(DP), INTENT(IN) :: ZoneWidth
    INTEGER,  INTENT(IN) :: iN_X, iNodeX1, iX1, iX2, iX3
    
    REAL(DP) :: TransProb(nE_G)
    INTEGER  :: iE, iN_E, m, iNode, iNodeE, iS, iCR

    DO iN_E = 1,nE_G
      DO m = 1,nM
        
        TransProb(iN_E) = One - MIN( See2(m,iN_E) , One )
        TransProb(iN_E) = 1.0d-3
        SigmaOsc(iN_E,iN_X,m) = TransProb(iN_E) * SpeedOfLightCGS / ZoneWidth
    
        WRITE(*,*) iN_E, See2(m,iN_E), REAL(fMatrixOsc(m,iN_E,1,1))
      END DO
    END DO
    
    ! Now convert to code units
    SigmaOsc(:,:,:) = SigmaOsc(:,:,:) * One / Second
    
  END SUBROUTINE CalculateOpacitiesOsc

  SUBROUTINE GetShockRadius( R_Shock_Loc)

    REAL(DP), INTENT(OUT) :: R_Shock_Loc

    R_Shock_Loc = R_Shock

  END SUBROUTINE GetShockRadius

  SUBROUTINE Get_Profile( iNodeX, iX1, iX2, iX3, &
                          Rnu_Out, R, dR, Rho, Ye )

    INTEGER,  INTENT(IN)  :: iX1, iX2, iX3, iNodeX
    REAL(DP), INTENT(OUT) :: Rnu_Out, R, dR, Ye, Rho
    
    INTEGER :: iNodeX1

    iNodeX1 = NodeNumberTableX(1,iNodeX) 
    
    Rnu_Out = Rnu
    R   = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
    dR  = MeshX(1) % Width(iX1) 
    Ye  = uAF(iNodeX,iX1,iX2,iX3,iAF_Ye)
    Rho = uPF(iNodeX,iX1,iX2,iX3,iPF_D)
    
    WRITE(*,*) 'CSI', CSI( R )
  
  END SUBROUTINE Get_Profile

  SUBROUTINE SetNonDiagonalEl( El_11, El_22, El_12,  &
                        Mixing_Optional )

    REAL(DP), INTENT(IN)  :: El_11, El_22
    REAL(DP), INTENT(OUT) :: El_12
    REAL(DP), INTENT(IN), OPTIONAL :: Mixing_Optional

    REAL(DP) :: Mixing
    REAL(DP) :: Lmax
    REAL(DP) :: El_z, El_x
    REAL(DP) :: Trace

    Mixing = One
    IF( PRESENT(Mixing_Optional) ) &
      Mixing = Mixing_Optional

    Trace = El_11 + El_22
    Lmax = MIN( Half * Trace, One - Half * Trace )
    El_z = Half*( El_11 - El_22 )

    El_x = SQRT(Lmax**2 - El_z**2)

    El_12 = El_x * Mixing

  END SUBROUTINE SetNonDiagonalEl

END MODULE ThornadoSQAInterfaceModule
