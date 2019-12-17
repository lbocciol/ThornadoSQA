MODULE InputOutputRelaxationModule 

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    UnitsDisplay
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    iX_B0, iX_E0, &
    nDOF, nNodesE, nNodesX
  USE ReferenceElementModule, ONLY: &
    NodeNumberTable
  USE RadiationFieldsModule, ONLY: &
    uPR, uCR, nCR, iPR_D, iCR_N, &  
    iPR_I1, iCR_G1, &
    iPR_I2, iCR_G2, &
    iPR_I3, iCR_G3, &
    nSpecies
  USE InitializationModule, ONLY: &
    nX_G, nE_G, nM, nF, &
    ChiOsc, EtaOsc, &
    fMatrixfOsc

  USE HDF5

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: WritefOscillations

  CHARACTER(9),  PARAMETER :: &
    OutputDirectory = '../Output'
  CHARACTER(13), PARAMETER :: &
    ImPartSuffix = 'ImaginaryPart'
  CHARACTER(15), PARAMETER :: &
    RadiationSuffix = 'RadiationFields'
  CHARACTER(12),  PARAMETER :: &
    OpacitySuffixOsc   = 'OpacitiesOsc'
  CHARACTER(12),  PARAMETER :: &
    OpacitySuffixStd   = 'OpacitiesStd'
  CHARACTER(10),  PARAMETER :: &
    fMatrixfOscSuffix   = 'fMatrixfOsc'

  INTEGER :: FileNumberfOsc = 0  
  INTEGER :: FileNumberOp   = 111111  
  INTEGER :: FileNumberRestart = 0
  !INTEGER :: FileNumberRestart = FileNumber !This way you don't erase previous files

  INTEGER :: HDFERR

CONTAINS

  SUBROUTINE InitializeFromRestart( RestartNumber, Time, &
          nX, swX, nE, swE) 

    INTEGER, INTENT(IN) :: nX(3), swX(3)
    INTEGER, INTENT(IN) :: nE,    swE
    INTEGER, INTENT(IN) :: RestartNumber

    REAL(DP), INTENT(INOUT) :: Time
    
    REAL(DP), ALLOCATABLE :: Dummy(:,:,:,:)

    INTEGER(HSIZE_T) :: DATASIZE(1)
    INTEGER(HID_T)   :: FILE_ID
    INTEGER(HID_T)   :: DSPACE_ID
    INTEGER(HID_T)   :: DSET_ID
    INTEGER(HID_T)   :: DIMS1(1), DIMS3(3), DIMS4(4)

    INTEGER :: iS

    CHARACTER(256) :: FileName
    CHARACTER(256) :: DatasetName 
    CHARACTER(6)   :: FileNumberString
    CHARACTER(10)  :: SpeciesString
   
    FileNumberRestart = RestartNumber
    WRITE( FileNumberString, FMT='(i6.6)') FileNumberRestart

    CALL H5OPEN_F( HDFERR )
    
    ! --- Open Radiation File --- !
    FileName &
      = OutputDirectory // '/' // &
        TRIM( ProgramName ) // '_' // &
        RadiationSuffix // '_' // &
        FileNumberString // '.h5'

    CALL H5FOPEN_F( TRIM( FileName ), &
        H5F_ACC_RDONLY_F, FILE_ID, HDFERR )

    IF( HDFERR < 0 ) THEN

      WRITE(*,'(A12,A1,I6.6,A1,A9)') &
          'Restart File','',RestartNumber,'','Not Found'
      STOP

    END IF

    DIMS4(1) = nE*nNodesE
    DIMS4(2) = nX(1)*nNodesX(1) 
    DIMS4(3) = nX(2)*nNodesX(2)
    DIMS4(4) = nX(3)*nNodesX(3)
   
    ALLOCATE( Dummy(DIMS4(1),DIMS4(2),DIMS4(3),DIMS4(4)) )

    DO iS = 1,nSpecies
      
      WRITE( SpeciesString, FMT='(A8,I2.2)') 'Species_',iS
      
      ! --- Read Number Density --- !
      DataSetName = '/Radiation Fields/' // SpeciesString // '/Primitive/Lagrangian Number Density'
    
      CALL H5DOPEN_F( FILE_ID, DatasetName, &
          DSET_ID, HDFERR )

      CALL H5DREAD_F( DSET_ID, H5T_NATIVE_DOUBLE, &
          Dummy, DIMS4, HDFERR )
    
      uPR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iPR_D,iS) = &
          Field4D_Inverse( Dummy, [ nE, nX(1), nX(2), nX(3) ], &
              [ nNodesE, nNodesX(1), nNodesX(2), nNodesX(3) ], &
              nDOF, NodeNumberTable )
    
      CALL H5DCLOSE_F( DSET_ID, HDFERR )
    
      DataSetName = '/Radiation Fields/' // SpeciesString // '/Conserved/Eulerian Number Density'

      CALL H5DOPEN_F( FILE_ID, DatasetName, &
          DSET_ID, HDFERR )

      CALL H5DREAD_F( DSET_ID, H5T_NATIVE_DOUBLE, &
          Dummy, DIMS4, HDFERR )

      uCR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iCR_N,iS) = &
          Field4D_Inverse( Dummy, [ nE, nX(1), nX(2), nX(3) ], &
              [ nNodesE, nNodesX(1), nNodesX(2), nNodesX(3) ], &
              nDOF, NodeNumberTable )
       
      CALL H5DCLOSE_F( DSET_ID, HDFERR )

      ! --- Read Flux (1) --- !
      DataSetName = '/Radiation Fields/' // SpeciesString // '/Primitive/Lagrangian Number Flux Density (1)'

      CALL H5DOPEN_F( FILE_ID, DatasetName, &
          DSET_ID, HDFERR )

      CALL H5DREAD_F( DSET_ID, H5T_NATIVE_DOUBLE, &
          Dummy, DIMS4, HDFERR )

      uPR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iPR_I1,iS) = &
          Field4D_Inverse( Dummy, [ nE, nX(1), nX(2), nX(3) ], &
              [ nNodesE, nNodesX(1), nNodesX(2), nNodesX(3) ], &
              nDOF, NodeNumberTable )


      CALL H5DCLOSE_F( DSET_ID, HDFERR )

      DataSetName = '/Radiation Fields/' // SpeciesString // '/Conserved/Eulerian Number Flux Density (1)'

      CALL H5DOPEN_F( FILE_ID, DatasetName, &
          DSET_ID, HDFERR ) 

      CALL H5DREAD_F( DSET_ID, H5T_NATIVE_DOUBLE, &
          Dummy, DIMS4, HDFERR )

      uCR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iCR_G1,iS) = &
          Field4D_Inverse( Dummy, [ nE, nX(1), nX(2), nX(3) ], &
              [ nNodesE, nNodesX(1), nNodesX(2), nNodesX(3) ], &
              nDOF, NodeNumberTable )

      CALL H5DCLOSE_F( DSET_ID, HDFERR )

      ! --- Read Flux (2) --- !
      DataSetName = '/Radiation Fields/' // SpeciesString // '/Primitive/Lagrangian Number Flux Density (2)'

      CALL H5DOPEN_F( FILE_ID, DatasetName, &
          DSET_ID, HDFERR )

      CALL H5DREAD_F( DSET_ID, H5T_NATIVE_DOUBLE, &
          Dummy, DIMS4, HDFERR )

      uPR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iPR_I2,iS) = &
          Field4D_Inverse( Dummy, [ nE, nX(1), nX(2), nX(3) ], &
              [ nNodesE, nNodesX(1), nNodesX(2), nNodesX(3) ], &
              nDOF, NodeNumberTable )

      CALL H5DCLOSE_F( DSET_ID, HDFERR )

      DataSetName = '/Radiation Fields/' // SpeciesString // '/Conserved/Eulerian Number Flux Density (2)'

      CALL H5DOPEN_F( FILE_ID, DatasetName, &
          DSET_ID, HDFERR )

      CALL H5DREAD_F( DSET_ID, H5T_NATIVE_DOUBLE, &
          Dummy, DIMS4, HDFERR )

      uCR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iCR_G2,iS) = &
          Field4D_Inverse( Dummy, [ nE, nX(1), nX(2), nX(3) ], &
              [ nNodesE, nNodesX(1), nNodesX(2), nNodesX(3) ], &
              nDOF, NodeNumberTable )

      CALL H5DCLOSE_F( DSET_ID, HDFERR )

      ! --- Read Flux (3) --- !
      DataSetName = '/Radiation Fields/' // SpeciesString // '/Primitive/Lagrangian Number Flux Density (3)'

      CALL H5DOPEN_F( FILE_ID, DatasetName, &
          DSET_ID, HDFERR )

      CALL H5DREAD_F( DSET_ID, H5T_NATIVE_DOUBLE, &
          Dummy, DIMS4, HDFERR )

      uPR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iPR_I3,iS) = &
          Field4D_Inverse( Dummy, [ nE, nX(1), nX(2), nX(3) ], &
              [ nNodesE, nNodesX(1), nNodesX(2), nNodesX(3) ], &
              nDOF, NodeNumberTable )

      CALL H5DCLOSE_F( DSET_ID, HDFERR )

      DataSetName = '/Radiation Fields/' // SpeciesString // '/Conserved/Eulerian Number Flux Density (3)'

      CALL H5DOPEN_F( FILE_ID, DatasetName, &
          DSET_ID, HDFERR )

      CALL H5DREAD_F( DSET_ID, H5T_NATIVE_DOUBLE, &
          Dummy, DIMS4, HDFERR )

      uCR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iCR_G3,iS) = &
          Field4D_Inverse( Dummy, [ nE, nX(1), nX(2), nX(3) ], &
              [ nNodesE, nNodesX(1), nNodesX(2), nNodesX(3) ], &
              nDOF, NodeNumberTable )

      CALL H5DCLOSE_F( DSET_ID, HDFERR )

    END DO
    
    ! --- Read Time --- !
    DataSetName = '/Time'
    DIMS1(1) = 1

    CALL H5DOPEN_F( FILE_ID, DatasetName, &
        DSET_ID, HDFERR )

    CALL H5DREAD_F( DSET_ID, H5T_NATIVE_DOUBLE, &
        Time, DIMS1, HDFERR )
    
    ASSOCIATE( U => UnitsDisplay )

      Time = Time * U % TimeUnit 

    END ASSOCIATE
    
    CALL H5DCLOSE_F( DSET_ID, HDFERR )
    
    CALL H5FCLOSE_F( FILE_ID, HDFERR )

    CALL H5CLOSE_F( HDFERR )
   
    DEALLOCATE( Dummy )

    FileNumberRestart = FileNumberRestart + 1

  END SUBROUTINE InitializeFromRestart

  FUNCTION Field4D_Inverse( F, nX, nN, nDOF, Tab )

    INTEGER,  INTENT(in) :: &
      nX(4), nN(4), nDOF, Tab(4,nDOF)
    REAL(DP), INTENT(in) :: &
      F(nX(1)*nN(1),nX(2)*nN(2),nX(3)*nN(3),nX(4)*nN(4))
    REAL(DP) :: &
      Field4D_Inverse(nDOF,nX(1),nX(2),nX(3),nX(4))

    INTEGER :: iX1, iX2, iX3, iX4
    INTEGER :: iN1, iN2, iN3, iN4, iNode
    
    DO iX4 = 1, nX(4)
      DO iX3 = 1, nX(3)
        DO iX2 = 1, nX(2)
          DO iX1 = 1, nX(1)

            DO iNode = 1, nDOF

              iN1 = Tab(1,iNode)
              iN2 = Tab(2,iNode)
              iN3 = Tab(3,iNode)
              iN4 = Tab(4,iNode)

              Field4D_Inverse(iNode,iX1,iX2,iX3,iX4) &
                = F( (iX1-1)*nN(1)+iN1, (iX2-1)*nN(2)+iN2, &
                  (iX3-1)*nN(3)+iN3, (iX4-1)*nN(4)+iN4 )

            END DO
          END DO
        END DO
      END DO
    END DO

    RETURN

  END FUNCTION Field4D_Inverse


  SUBROUTINE WritefOscillations( Time )

    REAL(DP), INTENT(IN) :: Time

    CHARACTER(6)   :: FileNumberString
    CHARACTER(256) :: FileName
    CHARACTER(256) :: DatasetName
    INTEGER(HID_T) :: FILE_ID
    
    REAL(DP)       :: fToBeWritten(nM,nE_G,nF,nF,2)
    
    WRITE( FileNumberString, FMT='(i6.6)') FileNumberfOsc

    FileName &
      = OutputDirectory // '/' // &
        TRIM( ProgramName ) // '_' // &
        fMatrixfOscSuffix // '_' // &
        FileNumberString // '.h5'

    CALL H5OPEN_F( HDFERR )

    CALL H5FCREATE_F( TRIM( FileName ), H5F_ACC_TRUNC_F, FILE_ID, HDFERR )

    DatasetName = '/Time'

    CALL WriteDataset1DHDF &
           ( [ Time ], DatasetName, FILE_ID )

    fToBeWritten(:,:,:,:,1) = REAL ( fMatrixfOsc(:,:,:,:) )
    fToBeWritten(:,:,:,:,2) = AIMAG( fMatrixfOsc(:,:,:,:) )

    DatasetName = '/fMatrixfOsc'

    CALL WriteDataset5DHDF &
          ( fToBeWritten, DatasetName, FILE_ID )

    CALL H5FCLOSE_F( FILE_ID, HDFERR )

    CALL H5CLOSE_F( HDFERR )
    
    FileNumberfOsc = FileNumberfOsc + 1
      
  END SUBROUTINE WritefOscillations

  SUBROUTINE WriteDataset5DHDF( Dataset, DatasetName, FILE_ID )

    REAL(DP),         INTENT(in) :: Dataset(:,:,:,:,:)
    CHARACTER(LEN=*), INTENT(in) :: DatasetName
    INTEGER(HID_T),   INTENT(in) :: FILE_ID

    INTEGER(HSIZE_T) :: DATASIZE(5)
    INTEGER(HID_T)   :: DATASPACE_ID
    INTEGER(HID_T)   :: DATASET_ID

    DATASIZE = SHAPE( Dataset )

    CALL H5SCREATE_F( H5S_SIMPLE_F, DATASPACE_ID, HDFERR )

    CALL H5SSET_EXTENT_SIMPLE_F &
           ( DATASPACE_ID, 5, DATASIZE, DATASIZE, HDFERR )

    CALL H5DCREATE_F &
           ( FILE_ID, TRIM( DatasetName ), H5T_NATIVE_DOUBLE, &
             DATASPACE_ID, DATASET_ID, HDFERR )

    CALL H5DWRITE_F &
           ( DATASET_ID, H5T_NATIVE_DOUBLE, Dataset, DATASIZE, HDFERR )

    CALL H5SCLOSE_F( DATASPACE_ID, HDFERR )

    CALL H5DCLOSE_F( DATASET_ID, HDFERR )

  END SUBROUTINE WriteDataset5DHDF


  SUBROUTINE WriteDataset3DHDF( Dataset, DatasetName, FILE_ID )

    REAL(DP),         INTENT(in) :: Dataset(:,:,:)
    CHARACTER(LEN=*), INTENT(in) :: DatasetName
    INTEGER(HID_T),   INTENT(in) :: FILE_ID

    INTEGER(HSIZE_T) :: DATASIZE(3)
    INTEGER(HID_T)   :: DATASPACE_ID
    INTEGER(HID_T)   :: DATASET_ID

    DATASIZE = SHAPE( Dataset )

    CALL H5SCREATE_F( H5S_SIMPLE_F, DATASPACE_ID, HDFERR )

    CALL H5SSET_EXTENT_SIMPLE_F &
           ( DATASPACE_ID, 5, DATASIZE, DATASIZE, HDFERR )

    CALL H5DCREATE_F &
           ( FILE_ID, TRIM( DatasetName ), H5T_NATIVE_DOUBLE, &
             DATASPACE_ID, DATASET_ID, HDFERR )

    CALL H5DWRITE_F &
           ( DATASET_ID, H5T_NATIVE_DOUBLE, Dataset, DATASIZE, HDFERR )

    CALL H5SCLOSE_F( DATASPACE_ID, HDFERR )

    CALL H5DCLOSE_F( DATASET_ID, HDFERR )

  END SUBROUTINE WriteDataset3DHDF


  SUBROUTINE WriteDataset1DHDF( Dataset, DatasetName, FILE_ID )

    REAL(DP),         INTENT(in) :: Dataset(:)
    CHARACTER(LEN=*), INTENT(in) :: DatasetName
    INTEGER(HID_T),   INTENT(in) :: FILE_ID

    INTEGER(HSIZE_T) :: DATASIZE(1)
    INTEGER(HID_T)   :: DATASPACE_ID
    INTEGER(HID_T)   :: DATASET_ID

    DATASIZE = SIZE( Dataset )

    CALL H5SCREATE_F( H5S_SIMPLE_F, DATASPACE_ID, HDFERR )

    CALL H5SSET_EXTENT_SIMPLE_F &
           ( DATASPACE_ID, 1, DATASIZE, DATASIZE, HDFERR )

    CALL H5DCREATE_F &
           ( FILE_ID, TRIM( DatasetName ), H5T_NATIVE_DOUBLE, &
             DATASPACE_ID, DATASET_ID, HDFERR )

    call H5DWRITE_F &
           ( DATASET_ID, H5T_NATIVE_DOUBLE, Dataset, DATASIZE, HDFERR )

    call H5SCLOSE_F( DATASPACE_ID, HDFERR )

    call H5DCLOSE_F( DATASET_ID, HDFERR )

  END SUBROUTINE WriteDataset1DHDF



  SUBROUTINE WriteDataset2DHDF( Dataset, DatasetName, FILE_ID )

    REAL(DP),         INTENT(in) :: Dataset(:,:)
    CHARACTER(LEN=*), INTENT(in) :: DatasetName
    INTEGER(HID_T),   INTENT(in) :: FILE_ID

    INTEGER(HSIZE_T) :: DATASIZE(2)
    INTEGER(HID_T)   :: DATASPACE_ID
    INTEGER(HID_T)   :: DATASET_ID

    DATASIZE = SHAPE( Dataset )

    CALL H5SCREATE_F( H5S_SIMPLE_F, DATASPACE_ID, HDFERR )

    CALL H5SSET_EXTENT_SIMPLE_F &
           ( DATASPACE_ID, 2, DATASIZE, DATASIZE, HDFERR )

    CALL H5DCREATE_F &
           ( FILE_ID, TRIM( DatasetName ), H5T_NATIVE_DOUBLE, &
             DATASPACE_ID, DATASET_ID, HDFERR )

    CALL H5DWRITE_F &
           ( DATASET_ID, H5T_NATIVE_DOUBLE, Dataset, DATASIZE, HDFERR )

    CALL H5SCLOSE_F( DATASPACE_ID, HDFERR )

    CALL H5DCLOSE_F( DATASET_ID, HDFERR )

  END SUBROUTINE WriteDataset2DHDF


  SUBROUTINE CreateGroupHDF( FileName, GroupName, FILE_ID )

    CHARACTER(len=*), INTENT(in) :: FileName
    CHARACTER(len=*), INTENT(in) :: GroupName
    INTEGER(HID_T),   INTENT(in) :: FILE_ID

    INTEGER(HID_T) :: GROUP_ID

    CALL H5GCREATE_F( FILE_ID, TRIM( GroupName ), GROUP_ID, HDFERR )

    CALL H5GCLOSE_F( GROUP_ID, HDFERR )

  END SUBROUTINE CreateGroupHDF


END MODULE InputOutputRelaxationModule 
