MODULE InitializationModule

  USE KindModule, ONLY: &
    DP, Zero, Half, One, Two, &
    Five, Pi
  USE UnitsModule, ONLY: &
    Gram, &
    Centimeter, &
    Kilometer,  &
    Kelvin, MeV, &
    BoltzmannConstant
  USE ProgramHeaderModule, ONLY: &       
    iX_B0, iX_E0, &                      
    iE_B0, iE_E0, &
    iZ_B0, iZ_E0, &
    nDOF, nDOFX,  &
    nDOFE, nNodesE
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE ReferenceElementModule, ONLY: &  
    NodeNumbersX, &                    
    NodeNumberTable                    
  USE MeshModule, ONLY: &
    MeshE, MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &                
    CoordinateSystem, &                            
    uGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33  
  USE FluidFieldsModule, ONLY: &
    uPF, iPF_D, &
    uAF, iAF_P, iAF_T, iAF_Ye, iAF_S, iAF_E, &
    iAF_Me, iAF_Mp, iAF_Mn, iAF_Xp, iAF_Xn, &
    iAF_Xa, iAF_Xh, iAF_Gm
  USE RadiationFieldsModule, ONLY: &                 
    nSpecies, iNuE, iNuE_Bar, &
    uPR, nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3, &       
    uCR, nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3          
  USE EquationOfStateModule_TABLE, ONLY: &
    ApplyEquationOfState_TABLE
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    FermiDirac
  USE TwoMoment_UtilitiesModule, ONLY: &  
    ComputeConserved_TwoMoment            
  USE ReadProfileModule, ONLY: &
    ReadChimeraMeshDimensions, &
    ReadChimeraRadiationDimensions, &
    ReadChimeraEnergyGrid, &
    Read1DChimeraProfile, &
    ReadChimeraMoments, &
    ReadGR1DProfile
  USE InterpolationModule, ONLY: &
    LinearInterpolation1D, &
    LinearInterpolation2D, &
    QuadraticInterpolation1D

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_Relaxation_SQA
  PUBLIC :: InitializeOscillations
  PUBLIC :: CofactorMatrix

  ! --- OSCILLATIONS RELATED QUANTITIES --- !
  ! --- DIMENSIONS OF ARRAYS --- !
  INTEGER,            PUBLIC :: nX_G, nE_G
  INTEGER,            PUBLIC :: nF
  INTEGER, PARAMETER, PUBLIC :: nM = 2
  INTEGER, PARAMETER, PUBLIC :: nY = 6
  INTEGER, PARAMETER, PUBLIC :: nS = 2

  REAL(DP),    ALLOCATABLE, PUBLIC :: Enu(:)
  REAL(DP),    ALLOCATABLE, PUBLIC :: IdentityR(:,:)
  COMPLEX(DP), ALLOCATABLE, PUBLIC :: YIdentity(:,:)
  COMPLEX(DP), ALLOCATABLE, PUBLIC :: IdentityC(:,:)
  COMPLEX(DP), ALLOCATABLE, PUBLIC :: AV(:,:,:), CV(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE, PUBLIC :: Hvf(:,:,:,:) !vacuum hamiltonian in flavor basis
  COMPLEX(DP),              PUBLIC :: PMNS(2,2)
  COMPLEX(DP), ALLOCATABLE, PUBLIC :: fMatrixOsc(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE, PUBLIC :: SMatrixOsc(:,:,:,:)

  REAL(DP),    ALLOCATABLE, PUBLIC :: Energies(:)

  ! --- RUNGE KUTTA PARAMETERS --- !
  INTEGER,  PARAMETER, PUBLIC :: nRK = 6
  REAL(DP), PARAMETER, PUBLIC :: nRKOrder = Five
  REAL(DP),            PUBLIC :: AA(6,5)
  REAL(DP), PARAMETER, PUBLIC :: BB5(6)=(/ 37.0d0/378.0d0 , 0.0d0 , 250.0d0/621.0d0 , &
                               125.0d0/594.0d0 , 0.0d0 , 512.0d0/1771.0d0 /)
  REAL(DP), PARAMETER, PUBLIC :: BB6(6)=(/ 2825.0d0/27648.0d0 , 0.0d0 , 18575.0d0/48384.0d0 , &
                               13525.0d0/55296.0d0 , 277.0d0/14336.0d0 , 1.0d0/4.0d0 /)
  REAL(DP), PARAMETER :: b2(5)=(/ 1.0d0/5.0d0 , 0.0d0 , 0.0d0 , 0.0d0 , 0.0d0/)
  REAL(DP), PARAMETER :: b3(5)=(/ 3.0d0/40.0d0 , 9.0d0/40.0d0 , 0.0d0 , 0.0d0 , 0.0d0/)
  REAL(DP), PARAMETER :: b4(5)=(/ 3.0d0/10.0d0 , -9.0d0/10.0d0 , 6.0d0/5.0d0 , 0.0d0 , 0.0d0 /)
  REAL(DP), PARAMETER :: b5(5)=(/ -11.0d0/54.0d0 , 5.0d0/2.0d0 , &
                            -70.0d0/27.0d0 , 35.0d0/27.0d0 ,0.0d0/)
  REAL(DP), PARAMETER :: b6(5)=(/ 1631.0d0 /55296.0d0 , 175.0d0/512.0d0 , 575.0d0/13824.0d0 , &
                              44275.0d0/110592.0d0 , 253.0d0/4096.0d0 /)


  ! --- Other useful quantities --- !
  REAL(DP),              PUBLIC :: R_Shock
  REAL(DP), ALLOCATABLE, PUBLIC :: ChiOsc(:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: EtaOsc(:,:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: SigmaOsc(:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC :: Psi0_loc(:,:)

  REAL(DP), ALLOCATABLE :: R_Ch(:)

CONTAINS

  SUBROUTINE InitializeFields_Relaxation_SQA( TimeSlice )
    
    REAL(DP), INTENT(INOUT) :: TimeSlice

    CALL InitializeFluidFields_Relaxation( TimeSlice )
    CALL InitializeRadiationFields_Relaxation( TimeSlice )
    !CALL InitializeRadiationFieldsFD_Relaxation
    
  END SUBROUTINE InitializeFields_Relaxation_SQA


  SUBROUTINE InitializeFluidFields_Relaxation( TimeSlice )

    REAL(DP), INTENT(INOUT)  :: TimeSlice

    INTEGER  :: iX1, iX2, iX3, iNodeX, iNodeX1
    REAL(DP) :: R, Y_New
    
    LOGICAL               :: ProfileFromGR1D
    INTEGER               :: nX1, nX2, nX3
    INTEGER               :: FileNumberMax
    REAL(DP), ALLOCATABLE :: D_Ch(:), T_Ch(:), Ye_Ch(:)
    
    ProfileFromGR1D = .FALSE.

    IF ( ProfileFromGR1D ) THEN
      
      nX1 = 600
      ALLOCATE(R_Ch(nX1),D_Ch(nX1),T_Ch(nX1),Ye_Ch(nX1))
      
      CALL ReadGR1DProfile(R_Ch, D_Ch, T_Ch, Ye_Ch, & 
            nX1, TimeSlice)

      R_Ch = R_Ch * Centimeter
      D_Ch = D_Ch * Gram / Centimeter ** 3
      T_Ch = T_Ch * MeV

    ELSE
      
      CALL ReadChimeraMeshDimensions( nX1, nX2, nX3 )

      ALLOCATE(R_Ch(nX1),D_Ch(nX1),T_Ch(nX1),Ye_Ch(nX1)) 
      
      FileNumberMax = 3000
      
      CALL Read1DChimeraProfile(R_Ch, D_Ch, T_Ch, Ye_Ch, &
            [nX1,nX2,nX3], R_Shock, TimeSlice, FileNumberMax) 
    
      R_Shock = R_Shock * Centimeter
      R_Ch = R_Ch * Centimeter
      D_Ch = D_Ch * Gram / Centimeter ** 3
      T_Ch = T_Ch * Kelvin

    END IF

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        R = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
    
        WRITE(*,*) R / Kilometer
        ! Interpolate Chimera Data to R and set variables
         
        CALL QuadraticInterpolation1D( R_Ch, D_Ch, nX1, R,  &
                    uPF(iNodeX,iX1,iX2,iX3,iPF_D) )

        CALL QuadraticInterpolation1D( R_Ch, T_Ch, nX1, R,  &
                    uAF(iNodeX,iX1,iX2,iX3,iAF_T) )

        CALL QuadraticInterpolation1D( R_Ch, Ye_Ch, nX1, R, &
                    uAF(iNodeX,iX1,iX2,iX3,iAF_Ye) )
    
      END DO

      CALL ApplyEquationOfState_TABLE &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uAF(:,iX1,iX2,iX3,iAF_T ), &
               uAF(:,iX1,iX2,iX3,iAF_Ye), uAF(:,iX1,iX2,iX3,iAF_P ), &
               uAF(:,iX1,iX2,iX3,iAF_S ), uAF(:,iX1,iX2,iX3,iAF_E ), &
               uAF(:,iX1,iX2,iX3,iAF_Me), uAF(:,iX1,iX2,iX3,iAF_Mp), &
               uAF(:,iX1,iX2,iX3,iAF_Mn), uAF(:,iX1,iX2,iX3,iAF_Xp), &
               uAF(:,iX1,iX2,iX3,iAF_Xn), uAF(:,iX1,iX2,iX3,iAF_Xa), &
               uAF(:,iX1,iX2,iX3,iAF_Xh), uAF(:,iX1,iX2,iX3,iAF_Gm) )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFluidFields_Relaxation

  SUBROUTINE InitializeRadiationFields_Relaxation( TimeSlice )

    REAL(DP), INTENT(INOUT) :: TimeSlice

    INTEGER  :: iE, iX1, iX2, iX3, iS
    INTEGER  :: iNodeX, iNodeX1, iNodeX2, iNodeX3
    INTEGER  :: iNode, iNodeE
    INTEGER  :: iN_E
    INTEGER  :: nZ(4), nX(3)
    
    REAL(DP) :: Gm_dd_11(nDOF) 
    REAL(DP) :: Gm_dd_22(nDOF) 
    REAL(DP) :: Gm_dd_33(nDOF) 

    REAL(DP) :: E, R
   
    !Chimera stuff
    INTEGER  :: nX1, nX2, nX3
    INTEGER  :: nE_Ch, nS_Ch
    INTEGER  :: FileNumberMax
    REAL(DP), ALLOCATABLE :: Psi0(:,:,:)
    REAL(DP), ALLOCATABLE :: Psi1(:,:,:)
    REAL(DP), ALLOCATABLE :: E_Ch(:)


    ! -- Create Array With Flattened Energies --- !
    nZ = iZ_E0 - iZ_B0 + 1
    nX = nZ(2:4)
    nE_G = nNodesE * nZ(1)
    nX_G = nDOFX * PRODUCT( nX )

    ALLOCATE( Energies( nE_G ) )

    DO iN_E = 1, nE_G

      iE     = MOD( (iN_E-1) / nDOFE, nZ(1) ) + iE_B0
      iNodeE = MOD( (iN_E-1)        , nDOFE ) + 1

      Energies(iN_E) = NodeCoordinate( MeshE, iE, iNodeE )

    END DO

    ! Read First Moment From Chimera
    
    CALL ReadChimeraMeshDimensions( nX1, nX2, nX3 )

    CALL ReadChimeraRadiationDimensions( nE_Ch, nS_Ch )

    ALLOCATE( E_Ch(nE_Ch) )

    CALL ReadChimeraEnergyGrid( E_Ch, nE_Ch )
    
    E_Ch = E_Ch * MeV    
    
    ALLOCATE( Psi0(nE_Ch, nS_Ch, nX1) )
    ALLOCATE( Psi1(nE_Ch, nS_Ch, nX1) )
    
    FileNumberMax = 3000
    
    CALL ReadChimeraMoments( Psi0, Psi1, [nE_Ch,nS_Ch,nX1,nX2,nX3], &
                nS_Ch, TimeSlice, FileNumberMax )
    
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iE = iE_B0, iE_E0
        
        DO iNode = 1, nDOF

          iNodeE  = NodeNumberTable(1,iNode)

          iNodeX1 = NodeNumberTable(2,iNode) 
          
          R = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
          
          E = NodeCoordinate( MeshE, iE, iNodeE )
            
          DO iS = 1,nSpecies
            
            CALL LinearInterpolation2D( E_Ch, R_Ch, Psi0(:,iS,:), &
                              nE_Ch, nX1, E, R, &
                              uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS) )
            
            CALL LinearInterpolation2D( E_Ch, R_Ch, Psi1(:,iS,:) / 3.0_DP, &
                              nE_Ch, nX1, E, R, &
                              uPR(iNode,iE,iX1,iX2,iX3,iPR_I1,iS) )

            !uPR(iNode,iE,iX1,iX2,iX3,iPR_D ,iS) = 1d-149

            !uPR(iNode,iE,iX1,iX2,iX3,iPR_I1,iS) = Zero

            uPR(iNode,iE,iX1,iX2,iX3,iPR_I2,iS) = Zero

            uPR(iNode,iE,iX1,iX2,iX3,iPR_I3,iS) = Zero

        
          END DO
          
        END DO

        DO iS = 1,nSpecies

          Gm_dd_11 &
            = uGF(NodeNumbersX,iX1,iX2,iX3,iGF_Gm_dd_11)

          Gm_dd_22 &
            = uGF(NodeNumbersX,iX1,iX2,iX3,iGF_Gm_dd_22)

          Gm_dd_33 &
            = uGF(NodeNumbersX,iX1,iX2,iX3,iGF_Gm_dd_33)

          CALL ComputeConserved_TwoMoment &
               ( uPR(:,iE,iX1,iX2,iX3,iPR_D, iS), &
                 uPR(:,iE,iX1,iX2,iX3,iPR_I1,iS), &
                 uPR(:,iE,iX1,iX2,iX3,iPR_I2,iS), &
                 uPR(:,iE,iX1,iX2,iX3,iPR_I3,iS), &
                 uCR(:,iE,iX1,iX2,iX3,iCR_N, iS), &
                 uCR(:,iE,iX1,iX2,iX3,iCR_G1,iS), &
                 uCR(:,iE,iX1,iX2,iX3,iCR_G2,iS), &
                 uCR(:,iE,iX1,iX2,iX3,iCR_G3,iS), &
                 Gm_dd_11, Gm_dd_22, Gm_dd_33 )
        
        END DO

      END DO

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeRadiationFields_Relaxation


  SUBROUTINE InitializeRadiationFieldsFD_Relaxation
    
    INTEGER  :: nZ(4), nX(3)
    INTEGER  :: iN_E, iE, iX1, iX2, iX3, iS
    INTEGER  :: iNode, iNodeE
    REAL(DP) :: kT(nDOF)
    REAL(DP) :: Mnu(nDOF), E
    REAL(DP) :: Gm_dd_11(nDOF) 
    REAL(DP) :: Gm_dd_22(nDOF) 
    REAL(DP) :: Gm_dd_33(nDOF) 

    ! -- Create Array With Flattened Energies --- !
    nZ = iZ_E0 - iZ_B0 + 1
    nX = nZ(2:4)
    nE_G = nNodesE * nZ(1)
    nX_G = nDOFX * PRODUCT( nX )

    ALLOCATE( Energies( nE_G ) )

    DO iN_E = 1, nE_G

      iE     = MOD( (iN_E-1) / nDOFE, nZ(1) ) + iE_B0
      iNodeE = MOD( (iN_E-1)        , nDOFE ) + 1

      Energies(iN_E) = NodeCoordinate( MeshE, iE, iNodeE )

    END DO


    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iS = 1, nSpecies
  
      Gm_dd_11 &                                        
        = uGF(NodeNumbersX,iX1,iX2,iX3,iGF_Gm_dd_11)    
                                                        
      Gm_dd_22 &                                        
        = uGF(NodeNumbersX,iX1,iX2,iX3,iGF_Gm_dd_22)    
                                                        
      Gm_dd_33 &                                        
        = uGF(NodeNumbersX,iX1,iX2,iX3,iGF_Gm_dd_33)
  
      kT = BoltzmannConstant &
             * uAF(NodeNumbersX,iX1,iX2,iX3,iAF_T)
  
      IF( iS .EQ. iNuE )THEN
  
        Mnu = + ( uAF(NodeNumbersX,iX1,iX2,iX3,iAF_Me) &
                  + uAF(NodeNumbersX,iX1,iX2,iX3,iAF_Mp) &
                  - uAF(NodeNumbersX,iX1,iX2,iX3,iAF_Mn) )
  
      ELSEIF( iS .EQ. iNuE_Bar )THEN
  
        Mnu = - ( uAF(NodeNumbersX,iX1,iX2,iX3,iAF_Me) &
                  + uAF(NodeNumbersX,iX1,iX2,iX3,iAF_Mp) &
                  - uAF(NodeNumbersX,iX1,iX2,iX3,iAF_Mn) )
      
      ELSE
  
        Mnu(:) = Zero
      
      END IF
  
      DO iE = iE_B0, iE_E0
        DO iNode = 1, nDOF
  
          iNodeE = NodeNumberTable(1,iNode)
  
          E = NodeCoordinate( MeshE, iE, iNodeE )
          
          uPR(iNode,iE,iX1,iX2,iX3,iPR_D,iS) &
            = FermiDirac( E, Mnu(iNode), kT(iNode) )
        
          ! --- Set First moment to zero --- !

          uPR(iNode,iE,iX1,iX2,iX3,iPR_I1,iS) &
             = Zero
  
          uPR(iNode,iE,iX1,iX2,iX3,iPR_I2,iS) &
             = Zero
  
          uPR(iNode,iE,iX1,iX2,iX3,iPR_I3,iS) &
             = Zero
  
          ! --- Compute conserved moments --- !
  
          CALL ComputeConserved_TwoMoment &
                   ( uPR(:,iE,iX1,iX2,iX3,iPR_D, iS), &
                     uPR(:,iE,iX1,iX2,iX3,iPR_I1,iS), &
                     uPR(:,iE,iX1,iX2,iX3,iPR_I2,iS), &
                     uPR(:,iE,iX1,iX2,iX3,iPR_I3,iS), &
                     uCR(:,iE,iX1,iX2,iX3,iCR_N, iS), &
                     uCR(:,iE,iX1,iX2,iX3,iCR_G1,iS), &
                     uCR(:,iE,iX1,iX2,iX3,iCR_G2,iS), &
                     uCR(:,iE,iX1,iX2,iX3,iCR_G3,iS), &
                     Gm_dd_11, Gm_dd_22, Gm_dd_33 )

        END DO
      END DO

    END DO
    END DO
    END DO
    END DO


  END SUBROUTINE InitializeRadiationFieldsFD_Relaxation

  SUBROUTINE InitializeOscillations

    INTEGER :: k

    ALLOCATE( fMatrixOsc(nM,nE_G,nF,nF) )
    ALLOCATE( SMatrixOsc(nM,nE_G,nF,nF) )
    ALLOCATE( Enu(nE_G) )
    ALLOCATE( Psi0_loc(nE_G,nSpecies) )
    ALLOCATE( SigmaOsc(nE_G,nX_G,nM) )
    ALLOCATE( EtaOsc(nE_G,nX_G,nCR,nSpecies) )
    ALLOCATE( ChiOsc(nE_G,nX_G,nSpecies) )

    ChiOsc = Zero
    EtaOsc = Zero
    Enu = Energies / MeV

    ! Initialize array for RK
    DO k = 1,nRK-1

      AA(1,k) = Zero
      AA(2,k) = b2(k)
      AA(3,k) = b3(k)
      AA(4,k) = b4(k)
      AA(5,k) = b5(k)
      AA(6,k) = b6(k)

    END DO

    CALL InitializeMatricesOscillations

  END SUBROUTINE InitializeOscillations

  SUBROUTINE InitializeMatricesOscillations

    INTEGER :: x
    INTEGER :: jF, iE, m

    REAL(DP), PARAMETER :: eV_to_erg = 1.60218d-12
    REAL(DP), PARAMETER :: clite = 2.99792458d10 !cm/s
    REAL(DP), PARAMETER :: m1 = 0.0d0 * eV_to_erg**2/clite**4 ! g,
    REAL(DP), PARAMETER :: dm21 = 2.43e-3 * eV_to_erg**2/clite**4 !g
    REAL(DP), PARAMETER :: theta12V = 9.0d0 * Pi/180.0d0 ! rad
    REAL(DP), PARAMETER :: c12V = COS(theta12V)
    REAL(DP), PARAMETER :: s12V = SIN(theta12V)
    REAL(DP), PARAMETER :: SIN2thetaW = 0.23122

    REAL(DP) :: kV(nE_G,nF)
    REAL(DP) :: Hvm(nF,nF)

    ALLOCATE( IdentityR(nF,nF) )
    ALLOCATE( IdentityC(nF,nF) )
    ALLOCATE( YIdentity(nS,nY) )

    ALLOCATE( Hvf(nM,nE_G,nF,nF) )

    ALLOCATE( AV(nE_G,nF,nF)    )
    ALLOCATE( CV(nE_G,nF,nF,nF) )

    DO x = 1,nS

      YIdentity(x,:) = (/ Pi/2.0d0, Pi/2.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0 /)

    END DO

    IdentityR(:,:) = Zero
    IdentityC(:,:) = Zero

    DO jF = 1,nF

      IdentityR(jF,jF) = One
      IdentityC(jF,jF) = One

    END DO

    PMNS(1,1) = c12V
    PMNS(1,2) = s12V
    PMNS(2,1) = -s12V
    PMNS(2,2) = c12V

    !Initialize vector of mass Eigenvalues
    DO iE = 1,nE_G

      kV(iE,1) = m1**2* clite**4 /(Two * Enu(iE)*1e6*eV_to_erg)
      kV(iE,2) = (kV(iE,1) + dm21) * clite**4 /(Two * Enu(iE)*1e6*ev_to_erg)

    END DO

    !Initialize Hamiltonians
    Hvm(:,:) = Zero

    DO iE = 1,nE_G

      !Vacuum Hamiltonian, constant, no need to reevaluate it
      Hvm(1,1) = kV(iE,1)
      Hvm(2,2) = kV(iE,2)
      Hvf(1,iE,:,:) = Matmul(Matmul(PMNS,Hvm(:,:)),Conjg(Transpose(PMNS)))
      Hvf(2,iE,:,:) = Conjg(Hvf(1,iE,:,:))

    END DO

    DO m = 1,nM
      DO iE = 1,nE_G

        CV(iE,:,:,:) = CofactorMatrix(Hvf(m,iE,:,:), kV(iE,:))
        AV(iE,:,:)   = Evaluate_AV    (Hvf(m,iE,:,:), kV(iE,:))

      END DO
    END DO

  END SUBROUTINE InitializeMatricesOscillations

  FUNCTION CofactorMatrix(M,Lambda) RESULT(CM)

    COMPLEX(DP), INTENT(IN)  :: M(2,2)
    REAL(DP),    INTENT(IN)  :: Lambda(2)

    COMPLEX(DP) :: CM(2,2,2)

    !local
    INTEGER :: j

    DO j=1,2

        CM(j,1,1) = M(2,2) - Lambda(j)
        CM(j,1,2) = -M(2,1)
        CM(j,2,1) = CONJG(CM(j,1,2))
        CM(j,2,2) = M(1,1) - Lambda(j)

    END DO

  END FUNCTION CofactorMatrix

  FUNCTION Evaluate_AV(M,Lambda) RESULT(AV_loc)

    COMPLEX(DP), INTENT(IN)  :: M(2,2)
    REAL(DP),    INTENT(IN)  :: Lambda(2)

    COMPLEX(DP) :: Av_loc(2,2)

    !local
    INTEGER :: j
    REAL(DP) :: Delta
    REAL(DP) :: re2, rmu2
    COMPLEX(DP) :: CM(2,2,2)

    CM = CofactorMatrix(M,Lambda)

    DO j=1,2

        IF(j == 1) Delta = ( lambda(2) - lambda(1) )
        IF(j == 2) Delta = ( lambda(1) - lambda(2) )

        re2  = Delta * REAL( CM(j,1,1) )
        rmu2 = Delta * REAL( CM(j,2,2) )

        IF( ABS(PMNS(1,j)) > ABS(PMNS(2,j)) ) THEN

            AV_loc(j,1)=REAL(PMNS(1,j) * SQRT(re2) / CM(j,1,1))
            AV_loc(j,2)=REAL(PMNS(1,j) * SQRT(rmu2) / CM(j,2,1))

        END IF

        IF( ABS(PMNS(2,j)) > ABS(PMNS(1,j)) ) THEN

            AV_loc(j,1)=REAL(PMNS(2,j)*SQRT(re2) / CM(j,1,2))
            AV_loc(j,2)=REAL(PMNS(2,j)*SQRT(rmu2) / CM(j,2,2))

        END IF
    END DO

  END FUNCTION Evaluate_AV


END MODULE InitializationModule
