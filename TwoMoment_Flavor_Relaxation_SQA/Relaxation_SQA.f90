PROGRAM Relaxation_SQA

  USE KindModule, ONLY: &
    DP, Zero, One, Two, &
    Pi, TwoPi, SqrtTiny, &
    Half, Five
  USE UnitsModule, ONLY: &
    Kilometer, MeV, Millisecond, &
    Second, Centimeter
  USE PhysicalConstantsModule, ONLY: &
    SpeedOfLightCGS
  USE ProgramHeaderModule, ONLY: &
    nZ, nX, nNodesZ, &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    iE_B0, iE_E0, iE_B1, iE_E1, &
    iZ_B0, iZ_E0, iZ_B1, iZ_E1
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE TimersModule, ONLY: &
    InitializeTimers, &
    FinalizeTimers, &
    TimersStart, &
    TimersStop, &
    Timer_InputOutput, &
    Timer_Evolve
  USE ReferenceElementModuleX, ONLY: &
    InitializeReferenceElementX, &
    FinalizeReferenceElementX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    InitializeReferenceElementX_Lagrange, &
    FinalizeReferenceElementX_Lagrange
  USE ReferenceElementModuleE, ONLY: &
    InitializeReferenceElementE, &
    FinalizeReferenceElementE
  USE ReferenceElementModuleE_Lagrange, ONLY: &
    InitializeReferenceElementE_Lagrange, &
    FinalizeReferenceElementE_Lagrange
  USE ReferenceElementModule, ONLY: &
    InitializeReferenceElement, &
    FinalizeReferenceElement
  USE ReferenceElementModule_Lagrange, ONLY: &
    InitializeReferenceElement_Lagrange, &
    FinalizeReferenceElement_Lagrange
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE GeometryComputationModuleE, ONLY: &
    ComputeGeometryE
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE MeshModule, ONLY: &
    MeshX
  USE FluidFieldsModule, ONLY: &
    uCF, uPF, uAF, &
    iPF_D, iAF_T, iAF_Ye
  USE RadiationFieldsModule, ONLY: &
    uCR, uPR, iPR_D
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF, FileNumber
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState
  USE OpacityModule_TABLE, ONLY: &
    InitializeOpacities_TABLE, &
    FinalizeOpacities_TABLE
  USE NeutrinoOpacitiesModule, ONLY: &
    CreateNeutrinoOpacities, &
    DestroyNeutrinoOpacities
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    ComputeNeutrinoOpacities  
  USE TwoMoment_ClosureModule, ONLY: &
    InitializeClosure_TwoMoment
  USE TwoMoment_PositivityLimiterModule, ONLY: &
    InitializePositivityLimiter_TwoMoment, &    
    FinalizePositivityLimiter_TwoMoment, &
    ApplyPositivityLimiter_TwoMoment
  USE IntegrationModule, ONLY: &
    Update_My_IMEX_PDARS
  USE InitializationModule, ONLY: &
    InitializeFields_Relaxation_SQA, &
    InitializeOscillations, &
    fMatrixOsc, SMatrixOsc, &
    nF, nE_G, nX_G, &
    R_Shock
  USE ThornadoSQAInterfaceModule, ONLY: &
    SQADriver, &
    IntegrationDriver
  USE InputOutputRelaxationModule, ONLY: &
    InitializeFromRestart

  IMPLICIT NONE

  LOGICAL  :: wrt
  LOGICAL  :: Restart
  LOGICAL  :: StartSQA, StopSQA

  INTEGER  :: RestartNumber
  INTEGER  :: iCycle, iCycleD
  INTEGER  :: nNodes, nSpecies
  INTEGER  :: nE, bcE, bcX(3)
  INTEGER  :: swE, swX(3)

  REAL(DP) :: eL, eR, xL(3), xR(3)
  REAL(DP) :: ZoomE, ZoomX(3)
  REAL(DP) :: t, dt_block, t_end
  REAL(DP) :: tmax
  REAL(DP) :: dt, t_SQA
  REAL(DP) :: t_wrt, dt_wrt
  REAL(DP) :: TimeSlice = 0.3d0 !Seconds

  Restart = .FALSE.
  RestartNumber = 101

  nNodes = 2
  nSpecies = 4
  nF = 2

  nX  = [ 128, 1, 1 ]
  xL  = [ 0.0d0 * Kilometer, 0.0_DP, 0.0_DP ]
  xR  = [ 500.0d0 * Kilometer, Pi,     TwoPi  ]
  bcX = [ 32, 0, 0 ]
  swX = [ 1, 1, 1 ]

  ZoomX = [ 1.01413296311415_DP, 1.0_DP, 1.0_DP ]
  ZoomX = [ 1.0_DP, 1.0_DP, 1.0_DP ]  

  nE  = 16
  eL  = 0.0d0 * MeV
  eR  = 2.5d2 * MeV
  bcE = 0
  swE = 0
  ZoomE = 1.09791_DP

  t      = 0.0_DP
  dt_wrt = 2.0d-2 * MilliSecond
  t_SQA  = 1.0d-9  * MilliSecond
  t_end  = 2.0d0  * MilliSecond 

  iCycleD = 1

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'Relaxation_SQA', &
           nX_Option &
             = nX, &
           swX_Option &
             = swX, &
           bcX_Option &
             = bcX, &
           xL_Option &
             = xL, &
           xR_Option &
             = xR, &
           ZoomX_Option &
             = ZoomX, &
           nE_Option &
             = nE, &
           swE_Option &
             = swE, &
           bcE_Option &
             = bcE, &
           eL_Option &
             = eL, &
           eR_Option &
             = eR, &
           ZoomE_Option &
             = ZoomE, &
           nNodes_Option &
             = nNodes, &
           CoordinateSystem_Option &
             = 'SPHERICAL', &
           ActivateUnits_Option &
             = .TRUE., &
           nSpecies_Option &
             = nSpecies, &
           BasicInitialization_Option &
             = .TRUE. )

  ! --- Position Space Reference Element and Geometry ---

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL ComputeGeometryX &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

  ! --- Energy Space Reference Element and Geometry ---

  CALL InitializeReferenceElementE

  CALL InitializeReferenceElementE_Lagrange

  CALL ComputeGeometryE &
         ( iE_B0, iE_E0, iE_B1, iE_E1, uGE )

  ! --- Phase Space Reference Element ---

  CALL InitializeReferenceElement

  CALL InitializeReferenceElement_Lagrange

  ! --- Initialize Equation of State ---

  CALL InitializeEquationOfState &
         ( EquationOfState_Option &
             = 'TABLE', &
           EquationOfStateTableName_Option &
             = 'EquationOfStateTable.h5' )
  
  ! --- Initialize and create Opacities

  CALL InitializeOpacities_TABLE &                         
         ( OpacityTableName_EmAb_Option &                  
             = 'wl-Op-SFHo-15-25-50-E40-B85-AbEm.h5', &    
           OpacityTableName_Iso_Option  &                  
             = 'wl-Op-SFHo-15-25-50-E40-B85-Iso.h5', &     
!           OpacityTableName_NES_Option &                   
!             = 'wl-Op-SFHo-15-25-50-E40-B85-NES.h5', &     
!           OpacityTableName_Pair_Option &                  
!             = 'wl-Op-SFHo-15-25-50-E40-B85-Pair.h5', &    
           Verbose_Option = .TRUE. )                       

  CALL InitializePositivityLimiter_TwoMoment &     
       ( Min_1_Option = 0.0d0 + SqrtTiny, &      
         Max_1_Option = 1.0d0 - EPSILON(1.0d0), &
         Min_2_Option = 0.0d0 + SqrtTiny, &      
         UsePositivityLimiter_Option &           
           = .TRUE. )

  CALL CreateNeutrinoOpacities( nZ, nNodesZ, nF )

  ! --- Initialize Moment Closure ---

  CALL InitializeClosure_TwoMoment

  ! --- Create Profile and initialize Moments ---
  
  CALL InitializeFields_Relaxation_SQA( TimeSlice )

  CALL InitializeOscillations
  
  IF( Restart ) THEN

    CALL InitializeFromRestart( RestartNumber, &
        t, nX, swX, nE, swE )

    FileNumber = RestartNumber + 1
    
  ELSE

    t = 0.0_DP

  END IF
 
  CALL WriteFieldsHDF &
         ( Time = 0.0_DP, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE. )

  CALL ApplyPositivityLimiter_TwoMoment &
      ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, uGF, uCR )

  ! --- Evolve ---

  WRITE(*,*)
  WRITE(*,'(A6,A,ES8.2E2,A8,ES8.2E2)') &
    '', 'Evolving from t = ', t / Millisecond, &
    ' to t = ', t_end / Millisecond
  WRITE(*,*)

  iCycle   = 0
  t_wrt    = dt_wrt
  wrt      = .FALSE.
  StartSQA = .FALSE.
  StopSQA  = .FALSE.
  DO WHILE( t < t_end )
    
    iCycle = iCycle + 1

    dt = Half * MINVAL( MeshX(1) % Width(iX_B0(1):iX_E0(1)) ) &
          / ( 2.0_DP * DBLE( nNodes - 1 ) + 1.0_DP )

    IF( t + dt > t_SQA ) THEN
      
      StartSQA = .TRUE.

    END IF

    IF( t + dt > t_wrt )THEN

      dt    = t_wrt - t
      t_wrt = t_wrt + dt_wrt
      wrt   = .TRUE.

    END IF

    IF( MOD( iCycle, iCycleD ) == 0 )THEN
   
      WRITE(*,'(A8,A8,I8.8,A2,A4,ES12.6E2,A5,A5,ES12.6E2,A4)') &
          '', 'Cycle = ', iCycle, &
          '', 't = ',  t  / Second, &
          ' (s) ', 'dt = ', dt / Second, ' (s)'

    END IF

    CALL TimersStart( Timer_Evolve )

    IF( t + dt > t_end ) THEN

      dt = t_end - t

    END IF

    IF ( StartSQA .AND. .NOT. StopSQA ) THEN
    
      WRITE(*,*) 'Starting SQA oscillations'
      CALL SQADriver
      StopSQA = .TRUE.

    END IF
    
    CALL IntegrationDriver( dt )

    t = t + dt

    CALL TimersStop( Timer_Evolve )

    IF( wrt )THEN

      CALL TimersStart( Timer_InputOutput )
      
      CALL WriteFieldsHDF &
             ( Time = t, &
               WriteGF_Option = .FALSE., &
               WriteFF_Option = .FALSE., &
               WriteRF_Option = .TRUE., &
               WriteOP_Option = .FALSE. )

      CALL TimersStop( Timer_InputOutput )

      wrt = .FALSE.

    END IF

  END DO

  ! --- Write Final Solution ---

  CALL TimersStart( Timer_InputOutput )

  CALL WriteFieldsHDF &
         ( Time = t, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE., &
           WriteOP_Option = .TRUE. )

  CALL TimersStop( Timer_InputOutput )

  WRITE(*,*)
  WRITE(*,'(A6,A,I6.6,A,ES12.6E2,A)') &
    '', 'Finished ', iCycle, ' Cycles in ', Timer_Evolve, ' s'
  WRITE(*,*)

  CALL FinalizeTimers

  CALL FinalizeEquationOfState

  CALL FinalizeOpacities_TABLE
    
  CALL DestroyNeutrinoOpacities

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeReferenceElementE

  CALL FinalizeReferenceElementE_Lagrange

  CALL FinalizeReferenceElement

  CALL FinalizeReferenceElement_Lagrange
  
  CALL FinalizePositivityLimiter_TwoMoment

  CALL FinalizeProgram

END PROGRAM Relaxation_SQA
