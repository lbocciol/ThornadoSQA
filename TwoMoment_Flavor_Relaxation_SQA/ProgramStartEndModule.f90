MODULE ProgramStartEndModule 

  USE KindModule, ONLY: &
    DP, SqrtTiny
  USE UnitsModule, ONLY: &
    Kilometer, MeV, Second
  USE PhysicalConstantsModule, ONLY: &
    SpeedOfLightCGS
  USE ProgramHeaderModule, ONLY: &
    nZ, nNodesZ, &
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
    Timer_InputOutput
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
  USE RadiationFieldsModule, ONLY: &
    uCR
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF, FileNumber
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState
  USE OpacityModule_TABLE, ONLY: &
    InitializeOpacities_TABLE, &
    FinalizeOpacities_TABLE
  USE TwoMoment_PositivityLimiterModule, ONLY: &
    InitializePositivityLimiter_TwoMoment, &    
    FinalizePositivityLimiter_TwoMoment, &      
    ApplyPositivityLimiter_TwoMoment            
  USE NeutrinoOpacitiesModule, ONLY: &
    CreateNeutrinoOpacities, &
    DestroyNeutrinoOpacities
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    ComputeNeutrinoOpacities  
  USE TwoMoment_ClosureModule, ONLY: &
    InitializeClosure_TwoMoment
  USE InitializationModule, ONLY: &
    InitializeFields_Relaxation_SQA, &
    InitializeOscillations, &
    nF
  USE InputOutputRelaxationModule, ONLY: &
    InitializeFromRestart

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeUnits
  PUBLIC :: ProgramStart
  PUBLIC :: ProgramEnd
  PUBLIC :: DumpFields

CONTAINS

  SUBROUTINE InitializeUnits( Kilometer_Unit, &
            Second_Unit, MeV_Unit, SpeedOfLightCGS_Unit )

    REAL(DP), INTENT(OUT) :: Kilometer_Unit
    REAL(DP), INTENT(OUT) :: Second_Unit
    REAL(DP), INTENT(OUT) :: MeV_Unit
    REAL(DP), INTENT(OUT) :: SpeedOfLightCGS_Unit

    Kilometer_Unit       = Kilometer
    Second_Unit          = Second
    MeV_Unit             = MeV
    SpeedOfLightCGS_Unit = SpeedOfLightCGS
    
  END SUBROUTINE InitializeUnits

  SUBROUTINE ProgramStart( nNodes, nSpecies, nF_in, &
               nX, xL, xR, bcX, swX, ZoomX, &
               nE, eL, eR, bcE, swE, ZoomE, &
               TimeSlice, t, &
               Restart_Option, RestartNumber_Option )
  
    INTEGER, INTENT(IN)  :: nNodes, nSpecies, nF_in
    INTEGER, INTENT(IN)  :: nE, bcE, bcX(3), swX(3), swE
    INTEGER, INTENT(IN)  :: nX(3)

    REAL(DP), INTENT(IN) :: eL, eR, xL(3), xR(3)
    REAL(DP), INTENT(IN) :: ZoomX(3), ZoomE
    REAL(DP), INTENT(INOUT) :: TimeSlice, t
    
    LOGICAL, INTENT(IN), OPTIONAL :: Restart_Option
    INTEGER, INTENT(IN), OPTIONAL :: RestartNumber_Option
   
    LOGICAL :: Restart
    INTEGER :: RestartNumber

    IF ( PRESENT( Restart_Option ) ) THEN
      Restart = Restart_Option
    ELSE
      Restart = .FALSE.
    END IF

    IF ( PRESENT( RestartNumber_Option ) ) THEN
      RestartNumber = RestartNumber_Option
    ELSE
      RestartNumber = 0
    END IF

    nF = nF_in

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
  

  END SUBROUTINE ProgramStart

  SUBROUTINE DumpFields( t, WriteGF, WriteFF, &
                            WriteRF, WriteOp ) 

    REAL(DP), INTENT(IN) :: t
    LOGICAL,  INTENT(IN) :: WriteGF
    LOGICAL,  INTENT(IN) :: WriteFF
    LOGICAL,  INTENT(IN) :: WriteRF
    LOGICAL,  INTENT(IN) :: WriteOP

    CALL TimersStart( Timer_InputOutput )

    CALL WriteFieldsHDF &
           ( Time = t, &
             WriteGF_Option = .FALSE., &
             WriteFF_Option = .FALSE., &
             WriteRF_Option = .TRUE., &
             WriteOP_Option = .FALSE. )
  
    CALL TimersStop( Timer_InputOutput )
  
  END SUBROUTINE DumpFields

  SUBROUTINE ProgramEnd( t )

    REAL(DP), INTENT(INOUT) :: t
  
    ! --- Write Final Solution ---
  
    CALL TimersStart( Timer_InputOutput )
  
    CALL WriteFieldsHDF &
           ( Time = t, &
             WriteGF_Option = .TRUE., &
             WriteFF_Option = .TRUE., &
             WriteRF_Option = .TRUE., &
             WriteOP_Option = .TRUE. )
  
    CALL TimersStop( Timer_InputOutput )
  
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

  END SUBROUTINE ProgramEnd
 
END MODULE ProgramStartEndModule
