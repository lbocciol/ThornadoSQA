MODULE ReadProfileModule

  USE KindModule, ONLY: &
    DP
  
  USE HDF5

  IMPLICIT NONE
  PRIVATE
  
  PUBLIC :: ReadGR1DProfile
  PUBLIC :: ReadChimeraMeshDimensions
  PUBLIC :: ReadChimeraRadiationDimensions
  PUBLIC :: ReadChimeraEnergyGrid
  PUBLIC :: Read1DChimeraProfile
  PUBLIC :: ReadChimeraMoments

CONTAINS
    
  SUBROUTINE ReadGR1DProfile(R_GR1D, D_GR1D, T_GR1D, Ye_GR1D, n0, &
                            TimeSlice)
    
    REAL(DP), INTENT(INOUT) :: R_GR1D(n0), D_GR1D(n0)
    REAL(DP), INTENT(INOUT) :: T_GR1D(n0), Ye_GR1D(n0)
    REAL(DP), INTENT(INOUT) :: TimeSlice
    INTEGER, INTENT(IN)     :: n0
    
    INTEGER :: istate
    REAL(DP) :: Tbounce

    OPEN(UNIT = 101, FILE = "tbounce.dat", STATUS = "old", &
      FORM = 'formatted', IOSTAT = istate)
      READ(101,'(E27.18)') Tbounce
    CLOSE(101)

    TimeSlice = TimeSlice + Tbounce
    
    WRITE(*,*)

    CALL ReadGR1DProfileXg(R_GR1D, D_GR1D, n0, 4,      &
                           TimeSlice, "rho.xg")
    
    CALL ReadGR1DProfileXg(R_GR1D, T_GR1D, 600, 4,      &
                           TimeSlice, "temperature.xg")

    CALL ReadGR1DProfileXg(R_GR1D, Ye_GR1D, 600, 4,     &
                           TimeSlice, "ye.xg")
    
  END SUBROUTINE ReadGR1DProfile

  SUBROUTINE ReadChimeraMeshDimensions( nX1, nX2, nX3 )

    INTEGER, INTENT(OUT) :: nX1, nX2, nX3

    INTEGER :: nX(3)
    INTEGER :: HDFERR

    CHARACTER(256) :: FileName
    CHARACTER(256) :: SetName

    WRITE(*,*)
    WRITE(*,'(A4,A)') '','Reading Mesh Dimensions...'

    FileName = "CHIMERA/SFHo/E15-1DSFHo-Frames-00/" // &
               "chimera_00001_grid_1_01.h5"

    CALL H5OPEN_F( HDFERR )

    SetName = "mesh/array_dimensions"
    CALL ReadHDFVectorInt1D( FileName, SetName, nX, 3 )

    CALL H5CLOSE_F( HDFERR )

    nX1 = nX(1)
    nX2 = nX(2)
    nX3 = nX(3)

  END SUBROUTINE ReadChimeraMeshDimensions
    
  SUBROUTINE ReadChimeraRadiationDimensions( nE, nS )

    INTEGER, INTENT(OUT) :: nE, nS

    INTEGER :: nR(2)
    INTEGER :: HDFERR

    CHARACTER(256) :: FileName
    CHARACTER(256) :: SetName

    WRITE(*,'(A4,A)') '','Reading Radiation Dimensions...' 

    FileName = "CHIMERA/SFHo/E15-1DSFHo-Frames-00/" // &
               "chimera_00001_grid_1_01.h5"

    CALL H5OPEN_F( HDFERR )

    SetName = "radiation/raddim"
    CALL ReadHDFVectorInt1D( FileName, SetName, nR, 2 )

    CALL H5CLOSE_F( HDFERR )
    
    nE = nR(1)
    nS = nR(2)

  END SUBROUTINE ReadChimeraRadiationDimensions

  SUBROUTINE ReadChimeraEnergyGrid( E_Ch, nE_Ch)

    REAL(DP), INTENT(OUT)   :: E_Ch(nE_Ch)
    INTEGER,  INTENT(INOUT) :: nE_Ch

    INTEGER :: HDFERR

    CHARACTER(256) :: FileName
    CHARACTER(256) :: SetName

    WRITE(*,'(A4,A)') '','Reading Chimera Energy Grid...'

    FileName = "CHIMERA/SFHo/E15-1DSFHo-Frames-00/" // &
               "chimera_00001_grid_1_01.h5"

    CALL H5OPEN_F( HDFERR )

    SetName = "radiation/unui"
    CALL ReadHDFVectorDouble1D( FileName, SetName, E_Ch, nE_Ch )

    CALL H5CLOSE_F( HDFERR )

  END SUBROUTINE ReadChimeraEnergyGrid

  SUBROUTINE Read1DChimeraProfile( R_Ch, D_Ch, T_Ch, Ye_Ch, &
                    nX, R_Shock, TimeSlice, iMax )

    INTEGER,  INTENT(IN)  :: nX(3)
    REAL(DP), INTENT(OUT) :: R_Ch(:), D_Ch(:)
    REAL(DP), INTENT(OUT) :: T_Ch(:), Ye_Ch(:)
    REAL(DP), INTENT(OUT) :: R_Shock

    INTEGER, INTENT(IN)  :: iMax
    REAL(DP), INTENT(IN) :: TimeSlice              

    INTEGER        :: HDFERR
    INTEGER        :: iFile
    CHARACTER(5)   :: FileNumberString 
    CHARACTER(20)  :: DirectoryString
    CHARACTER(256) :: FileName
    CHARACTER(256) :: SetName

    REAL(DP)       :: Time, Tbounce

    REAL(DP), ALLOCATABLE :: Temp3D(:,:,:)
    REAL(DP)              :: Temp1D(1) 
    
    WRITE(*,*)
    WRITE(*,'(A4,A)') '','Reading Chimera Profiles...'

    CALL H5OPEN_F( HDFERR )

    DO iFile = 1,iMax

      WRITE(DirectoryString , FMT = '(A18,I2.2)') &
          'E15-1DSFHo-Frames-', (iFile / 1000) * 10
      WRITE(FileNumberString, FMT = '(I5.5)') iFile

      FileName = "CHIMERA/SFHo/" // DirectoryString // "/chimera_" &
            // FileNumberString // "_grid_1_01.h5"

      SetName = "mesh/time"
      CALL ReadHDFScalar(FileName, SetName, Time)
    
      SetName = "mesh/t_bounce"
      CALL ReadHDFScalar(FileName, SetName, Tbounce)

      IF ( Time - Tbounce >= TimeSlice .AND. Tbounce > 0 ) &
        EXIT

    END DO
    
    ALLOCATE( Temp3D(nX(1),nX(2),nX(3)) )

    SetName = "mesh/x_cf"
    CALL ReadHDFVectorDouble3D( FileName, SetName, Temp3D, nX )
    R_Ch = Temp3D(:,1,1)
    
    SetName = "fluid/rho_c"
    CALL ReadHDFVectorDouble3D( FileName, SetName, Temp3D, nX )
    D_Ch = Temp3D(:,1,1)

    SetName = "fluid/t_c"
    CALL ReadHDFVectorDouble3D( FileName, SetName, Temp3D, nX )
    T_Ch = Temp3D(:,1,1)

    SetName = "fluid/ye_c"
    CALL ReadHDFVectorDouble3D( FileName, SetName, Temp3D, nX )
    Ye_Ch = Temp3D(:,1,1)
    
    SetName = "analysis/r_shock"
    CALL ReadHDFVectorDouble1D( FileName, SetName, Temp1D, 1 )
    R_Shock = Temp1D(1)

    CALL H5CLOSE_F( HDFERR )

    DEALLOCATE( Temp3D )

  END SUBROUTINE Read1DChimeraProfile

  SUBROUTINE ReadChimeraMoments( Psi0, Psi1, nMoments, &
                nWeights, TimeSlice, iMax )

    INTEGER,  INTENT(IN)    :: nMoments(5), nWeights
    REAL(DP), INTENT(INOUT) :: Psi0(:,:,:), Psi1(:,:,:)

    INTEGER, INTENT(IN)  :: iMax
    REAL(DP), INTENT(IN) :: TimeSlice

    INTEGER        :: HDFERR
    INTEGER        :: iFile, iS
    CHARACTER(5)   :: FileNumberString
    CHARACTER(20)  :: DirectoryString
    CHARACTER(256) :: FileName
    CHARACTER(256) :: SetName

    REAL(DP)       :: Time, Tbounce

    REAL(DP), ALLOCATABLE :: Temp(:,:,:,:,:)
    INTEGER,  ALLOCATABLE :: Weights(:)

    WRITE(*,*)
    WRITE(*,'(A4,A)') '','Reading Chimera Moments...'

    CALL H5OPEN_F( HDFERR )

    DO iFile = 1,iMax

      WRITE(DirectoryString , FMT = '(A18,I2.2)') &
          'E15-1DSFHo-Frames-', (iFile / 1000) * 10
      WRITE(FileNumberString, FMT = '(I5.5)') iFile

      FileName = "CHIMERA/SFHo/" // DirectoryString // "/chimera_" &
            // FileNumberString // "_grid_1_01.h5"

      SetName = "mesh/time"
      CALL ReadHDFScalar(FileName, SetName, Time)

      SetName = "mesh/t_bounce"
      CALL ReadHDFScalar(FileName, SetName, Tbounce)

      IF ( Time - Tbounce >= TimeSlice .AND. Tbounce > 0 ) &
        EXIT

    END DO
    
    ALLOCATE( Temp(nMoments(1), &
                   nMoments(2), &
                   nMoments(3), &
                   nMoments(4), &
                   nMoments(5)) )

    ALLOCATE( Weights( nWeights ) )

    SetName = "radiation/stwt"
    CALL ReadHDFVectorInt1D( FileName, SetName, Weights, nWeights )

    SetName = "radiation/psi0_c"
    CALL ReadHDFVectorDouble5D( FileName, SetName, Temp, nMoments )
    
    DO iS = 1,nWeights
      
      Psi0(:,iS,:) = Temp(:,iS,:,1,1) / Weights(iS)

    END DO

    SetName = "radiation/psi1_e"
    CALL ReadHDFVectorDouble5D( FileName, SetName, Temp, nMoments )

    DO iS = 1,nWeights
      
      Psi1(:,iS,:) = Temp(:,iS,:,1,1) / Weights(iS)

    END DO
    
    DEALLOCATE( Temp, Weights )

    CALL H5CLOSE_F( HDFERR )

  END SUBROUTINE ReadChimeraMoments

  SUBROUTINE ReadHDFScalar(FileName, SetName, Variable)
    
    CHARACTER(256), INTENT(IN)  :: FileName
    CHARACTER(256), INTENT(IN)  :: SetName
    REAL(DP),       INTENT(OUT) :: Variable
    
    INTEGER(HID_T)              :: FILE_ID
    INTEGER(HID_T)              :: DSET_ID
    INTEGER                     :: HDFERR
    INTEGER(HSIZE_T)            :: DIMS(1)

    DIMS(1) = 1 

    ! --- Open File --- !
    CALL H5FOPEN_F( TRIM(ADJUSTL( FileName )), &
           H5F_ACC_RDONLY_F, FILE_ID, HDFERR)
    
    CALL H5DOPEN_F( FILE_ID, TRIM(ADJUSTL( SetName )), &
              DSET_ID, HDFERR)
    
    ! --- Read Data --- !
    CALL H5DREAD_F( DSET_ID, H5T_NATIVE_DOUBLE, & 
              Variable, DIMS, HDFERR )
    
    ! --- Close File --- !
    CALL H5DCLOSE_F( DSET_ID, HDFERR )
      
    CALL H5FCLOSE_F( FILE_ID, HDFERR )

  END SUBROUTINE ReadHDFScalar 

  SUBROUTINE ReadHDFVectorInt1D(FileName, SetName, Variable, nArray)

    CHARACTER(256), INTENT(IN)  :: FileName
    CHARACTER(256), INTENT(IN)  :: SetName
    INTEGER,        INTENT(IN)  :: nArray
    
    INTEGER,        INTENT(OUT) :: Variable(:)

    INTEGER(HID_T)              :: FILE_ID
    INTEGER(HID_T)              :: DSET_ID
    INTEGER                     :: HDFERR
    INTEGER(HSIZE_T)            :: DIMS(1)

    DIMS(1) = nArray
    
    ! --- Open File --- !
    CALL H5FOPEN_F( TRIM(ADJUSTL( FileName )), &
           H5F_ACC_RDONLY_F, FILE_ID, HDFERR)

    CALL H5DOPEN_F( FILE_ID, TRIM(ADJUSTL( SetName )), &
              DSET_ID, HDFERR)
    
    ! --- Read Data --- !
    CALL H5DREAD_F( DSET_ID, H5T_NATIVE_INTEGER, &
              Variable, DIMS, HDFERR )

    ! --- Close File --- !
    CALL H5DCLOSE_F( DSET_ID, HDFERR )

    CALL H5FCLOSE_F( FILE_ID, HDFERR )
 
  END SUBROUTINE ReadHDFVectorInt1D

  SUBROUTINE ReadHDFVectorDouble1D(FileName, SetName, Variable, nArray)

    CHARACTER(256), INTENT(IN)  :: FileName
    CHARACTER(256), INTENT(IN)  :: SetName
    INTEGER,        INTENT(IN)  :: nArray

    REAL(DP),       INTENT(OUT) :: Variable(:)

    INTEGER(HID_T)              :: FILE_ID
    INTEGER(HID_T)              :: DSET_ID
    INTEGER                     :: HDFERR
    INTEGER(HSIZE_T)            :: DIMS(1)

    DIMS(1) = nArray

    ! --- Open File --- !
    CALL H5FOPEN_F( TRIM(ADJUSTL( FileName )), &
           H5F_ACC_RDONLY_F, FILE_ID, HDFERR)

    CALL H5DOPEN_F( FILE_ID, TRIM(ADJUSTL( SetName )), &
              DSET_ID, HDFERR)

    ! --- Read Data --- !
    CALL H5DREAD_F( DSET_ID, H5T_NATIVE_DOUBLE, &
              Variable, DIMS, HDFERR )

    ! --- Close File --- !
    CALL H5DCLOSE_F( DSET_ID, HDFERR )

    CALL H5FCLOSE_F( FILE_ID, HDFERR )

  END SUBROUTINE ReadHDFVectorDouble1D


  SUBROUTINE ReadHDFVectorDouble3D( FileName, SetName, &
                        Variable, nArray )

    CHARACTER(256), INTENT(IN) :: FileName
    CHARACTER(256), INTENT(IN) :: SetName
    INTEGER,        INTENT(IN) :: nArray(3)
    
    REAL(DP),      INTENT(OUT) :: Variable(:,:,:)
    
    INTEGER(HID_T)             :: FILE_ID
    INTEGER(HID_T)             :: DSET_ID
    INTEGER(HSIZE_T)           :: DIMS(3)
    INTEGER                    :: HDFERR

    DIMS = nArray 
    
    ! --- Open File --- !
    CALL H5FOPEN_F( TRIM(ADJUSTL( FileName )), &
           H5F_ACC_RDONLY_F, FILE_ID, HDFERR)

    CALL H5DOPEN_F( FILE_ID, TRIM(ADJUSTL( SetName )), &
              DSET_ID, HDFERR)

    ! --- Read Data --- !
    CALL H5DREAD_F( DSET_ID, H5T_NATIVE_DOUBLE, &
              Variable, DIMS, HDFERR )
    
    ! --- Close File --- !
    CALL H5DCLOSE_F( DSET_ID, HDFERR )

    CALL H5FCLOSE_F( FILE_ID, HDFERR )

  END SUBROUTINE ReadHDFVectorDouble3D


  SUBROUTINE ReadHDFVectorDouble5D( FileName, SetName, & 
                        Variable, nArray )

    CHARACTER(256), INTENT(IN) :: FileName
    CHARACTER(256), INTENT(IN) :: SetName
    INTEGER,        INTENT(IN) :: nArray(5)

    REAL(DP),      INTENT(OUT) :: Variable(:,:,:,:,:)

    INTEGER(HID_T)             :: FILE_ID
    INTEGER(HID_T)             :: DSET_ID
    INTEGER(HSIZE_T)           :: DIMS(5)
    INTEGER                    :: HDFERR

    DIMS = nArray

    ! --- Open File --- !
    CALL H5FOPEN_F( TRIM(ADJUSTL( FileName )), &
           H5F_ACC_RDONLY_F, FILE_ID, HDFERR)

    CALL H5DOPEN_F( FILE_ID, TRIM(ADJUSTL( SetName )), &
              DSET_ID, HDFERR)

    ! --- Read Data --- !
    CALL H5DREAD_F( DSET_ID, H5T_NATIVE_DOUBLE, &
              Variable, DIMS, HDFERR )

    ! --- Close File --- !
    CALL H5DCLOSE_F( DSET_ID, HDFERR )

    CALL H5FCLOSE_F( FILE_ID, HDFERR )

  END SUBROUTINE ReadHDFVectorDouble5D


  SUBROUTINE ReadGR1DProfileXg(R_GR1D, U_GR1D, n0, nGhost, &
                               TimeSlice, FileName)

    REAL(DP), INTENT(INOUT)       :: R_GR1D(n0)
    REAL(DP), INTENT(INOUT)       :: U_GR1D(n0)
    INTEGER, INTENT(IN)           :: n0, nGhost
    REAL(DP), INTENT(IN)          :: TimeSlice
    CHARACTER(len=*), INTENT(IN)  :: FileName

    CHARACTER(len=200) :: TimeString
    REAL(DP) :: Time 
    INTEGER  :: i, istate

    111 FORMAT(1P20E18.9)
    
    WRITE(*,'(A4,2A)') '','Reading from GR1D: ',TRIM(ADJUSTL(FileName))
    
    OPEN(UNIT=999, FILE = TRIM(ADJUSTL(FileName)), STATUS = "old", &
        FORM = 'formatted', IOSTAT = istate) 

    DO WHILE (istate .ge. 0 )
      
      READ(999,"(A8,E30.9)") TimeString, Time
      
      !Skip ghost
      DO i = 1,nGhost
        READ(999,*)
      END DO

      DO i = 1,n0 
        READ(999,111) R_GR1D(i),U_GR1D(i)
      END DO
      
      !Skip ghost and 2 blank spaces
      DO i = 1,nGhost + 2
        READ(999,*) 
      END DO
      
      IF (Time .ge. TimeSlice) &
        EXIT
    END DO
    CLOSE(999)

  END SUBROUTINE ReadGR1DProfileXg

  SUBROUTINE ReadGR1DProfileDat(TimeSlice, Var, FileName)            
                                                                      
    REAL(DP), INTENT(INOUT)       :: Var
    REAL(DP), INTENT(IN)          :: TimeSlice         
    CHARACTER(len=*), INTENT(IN)  :: FileName                           
                                                                        
    REAL(DP) :: Time                                                    
    INTEGER  :: istate                                               
                                                                        
    111 FORMAT(1P20E18.9)                                               
    WRITE(*,'(A8,A)') "Reading",TRIM(ADJUSTL(FileName))                                                                                     
    OPEN(UNIT=999, FILE = TRIM(ADJUSTL(FileName)), STATUS = "old", &    
        FORM = 'formatted', IOSTAT = istate)                            
                                                                        
    DO WHILE (istate .ge. 0 )                                           
                                                                      
      READ(999,111) Time, Var                             
                                                           
      IF (Time .ge. TimeSlice) &
        EXIT
    END DO
    CLOSE(999)

  END SUBROUTINE ReadGR1DProfileDat

END MODULE ReadProfileModule 
