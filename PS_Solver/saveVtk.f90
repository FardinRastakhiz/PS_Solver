
    MODULE saveData

  
    IMPLICIT NONE
    SAVE
    !YOU MAY ADD DECLARATIONS HERE, TO BE USED BY SUBPROGRAMS BY ADDING "USE Input_data" TO THEM

    !*********************SUBROUTINMES INCLUDED IN THIS MODULE********************************
    CONTAINS
    
    subroutine saveVtk_int(openFileID, firstFlag, firstVtkType, CN,   Npoint, CONF,      NLink, DataType, dataVal, dataName)

   ! REAL*8,DIMENSION(:,:)::CN(3,NPoint),CELLS_vtk(NLink,3)
   ! Integer,DIMENSION(:,:)::CONF(2,NLink)
    REAL*8::CN(3,NPoint),CELLS_vtk(NLink,3)
    Integer::CONF(2,NLink)
    
    Integer:: dataVal(:)
    !DOUBLEPRECISION,DIMENSION(:),allocatable::dataVal
    !character(len=10) :: dataName
    CHARACTER(*):: dataName
    CHARACTER(*):: DataType
    INTEGER::firstFlag,I,J,numData,firstVtkType, openFileID, Npoint, Nlink
    LOGICAL:: vtkAscciFormat
    LOGICAL:: vtkPolydataFormat

    character :: buffer*100, lf*1, str1*8, str2*8, str3*8
    integer   :: ivtk = 9, int
                         
    lf = char(10) ! line feed character

    vtkAscciFormat = .TRUE.
    !vtkAscciFormat = .FALSE.
    vtkPolydataFormat =  .TRUE.
    !allocate(character(len=10) :: dataName)
    !ALLOCATE(dataVal(numData))
    CELLS_vtk(:,1)=2
    CELLS_vtk(:,2)=dble(CONF(1,:)-1)
    CELLS_vtk(:,3)=dble(CONF(2,:)-1)

    IF(firstFlag .EQ. 1)THEN
        IF(vtkAscciFormat)THEN
            IF(vtkPolydataFormat)THEN
                write (openFileID,'(A26)')"# vtk DataFile Version 3.0"
                write (openFileID,'(A12)')"# vtk output"
                write (openFileID,'(A5)')"ASCII"
                write (openFileID,'(A25)') " "
                write (openFileID,'(A16)') "DATASET POLYDATA"
                write (openFileID,'(A6,I15,A7)') "POINTS", NPoint, " double"
            ELSE
                write (openFileID,'(A26)')"# vtk DataFile Version 2.0"
                write (openFileID,'(A25)')"# SATURATION DISTRIBUTION"
                write (openFileID,'(A5)')"ASCII"
                write (openFileID,'(A25)') " "
                write (openFileID,'(A25)') "DATASET UNSTRUCTURED_GRID"
                write (openFileID,'(A6,I15,A7)') "POINTS", NPoint, " double"
            ENDIF
        ELSE
            CLOSE(openFileID)
            open(unit=openFileID,file='vtk_poreInlet_initial.vtk',form='binary',convert='BIG_ENDIAN')
            !open(unit=openFileID,file=savePath//'vtk_poreInlet_initial.vtk',access='stream')
            !open(unit=openFileID,file=savePath//'vtk_poreInlet_initial.vtk',form='binary')
            !         open(unit=openFileID,file=savePath//'vtk_poreInlet_initial.vtk', form='unformatted',&
            !recordtype='STREAM_LF',&
            !action='write'                          ,&
            !convert='BIG_ENDIAN'                    ,&
            !access='sequential'                      )

            buffer = '# vtk DataFile Version 3.0'//lf  ; write(openFileID) trim(buffer)
            buffer = 'vtk output'//lf                  ; write(openFileID) trim(buffer)
            buffer = 'BINARY'//lf                      ; write(openFileID) trim(buffer)
            buffer = 'DATASET UNSTRUCTURED_GRID'//lf    ; write(openFileID) trim(buffer)
            buffer = 'POINTS '//intToCharConvert(NPoint)//' double'//lf  ; write(openFileID) trim(buffer)
            !CLOSE(openFileID)
        ENDIF

        IF(vtkAscciFormat)THEN
            DO I = 1,NPoint
                DO J = 1,3
                    !set very small number to zero
                    IF( ABS( CN(J,I) ) .LT. 1d-20)THEN
                        CN(J,I) = 0D0
                    ENDIF
                ENDDO
                write (openFileID,'(3G20.7)') CN(:,I)
            ENDDO
        ELSE
            DO I = 1,NPoint
                DO J = 1,3
                    IF( ABS( CN(J,I) ).LT.1D-20)THEN
                        CN(J,I) = 0D0
                    ENDIF
                ENDDO
                write (openFileID) CN(:,I)
            ENDDO

            ! write (openFileID) (CN(:,I), I=1,NPoint)
            ! write (openFileID) (CELLS_vtk(I,1), I=1,10), lf
            ! write (openFileID) (CN(:,I),real(0.,kind=8) , I=1,NPoint)
            ! CELLS_vtk(1:NPoint,1) = CN(1,:)
            ! write(openFileID)(CELLS_vtk(I,1) , i=1,NPoint)
        ENDIF
        !CLOSE(openFileID)
        !write (openFileID,*) " "
        !write (openFileID,'(A5,2I10)') "CELLS", NLink, 3*NLink


        if(vtkAscciFormat)then
            IF(vtkPolydataFormat)THEN
                write (openFileID,*) " "
                write (openFileID,'(A5,2I15)') "LINES", NLink, 3*NLink    
            ELSE
                write (openFileID,*) " "
                write (openFileID,'(A5,2I15)') "CELLS", NLink, 3*NLink
            ENDIF
        else
            buffer = ' '//lf    ; write(openFileID) trim(buffer)
            buffer = 'CELLS '//intToCharConvert(NLink)//' '//intToCharConvert(3*NLink)//lf  ; write(openFileID) trim(buffer)
        endif

        
        IF(vtkAscciFormat)THEN
            DO I = 1,NLink
                write (openFileID,'(A2,2I15)') "2 ", CONF(1,I)-1, CONF(2,I)-1
            ENDDO
        ELSE
            !write (openFileID) ( 2+CONF(1,I)*0, I = 1,NLink), ( CONF(1,I)-1, I = 1,NLink), (CONF(2,I)-1, I = 1,NLink)
            ! write (openFileID) (CELLS_vtk(I,:), I=1,NLink)
            !write (openFileID) (CN(:,I), I=1,NPoint)
            DO I = 1,NLink
                write (openFileID) CELLS_vtk(I,:)
            ENDDO
        ENDIF

        !CLOSE(openFileID)
        IF(vtkPolydataFormat .EQ. .FALSE.)THEN
            ! write (openFileID,*) " "
            !write (openFileID,'(A10,I10)') "CELL_TYPES", NLink
            IF(vtkAscciFormat)THEN
                write (openFileID,*) " "
                write (openFileID,'(A10,I15)') "CELL_TYPES", NLink
            ELSE
                buffer = ' '//lf    ; write(openFileID) trim(buffer)
                buffer = 'CELL_TYPES '//intToCharConvert(NLink)//lf  ; write(openFileID) trim(buffer)
            ENDIF

            IF(vtkAscciFormat)THEN
                DO I = 1,NLink
                    write (openFileID,'(A1)') "3"
                ENDDO
            ELSE
                CELLS_vtk(:,1)=3
                !write (openFileID) (3 + CONF(1,I)*0,I = 1,NLink)
                ! write (openFileID) CELLS_vtk(:,1)
                DO I = 1,NLink
                    write (openFileID) CELLS_vtk(I,1)
                ENDDO
            ENDIF
        ENDIF
    ENDIF

    !CLOSE(openFileID)
                        
    IF(DataType .EQ. 'nnodeData')THEN
        IF(vtkAscciFormat)THEN
        write (openFileID,*) " "
        !add porebody properties
        IF(firstVtkType.EQ.1)THEN
            write (openFileID,'(A10,I15)') "POINT_DATA", NPoint
        ENDIF
        write (openFileID,'(A8,A,A7)') "SCALARS ",dataName, " int"
        write (openFileID,'(A20)') "LOOKUP_TABLE default"
        ELSE
        !IF(0)THEN
            write (openFileID) " "
            !add porebody properties
            IF(firstVtkType.EQ.1)THEN
                write (openFileID) "POINT_DATA", NPoint
            ENDIF
            write (openFileID) "SCALARS ",dataName, " int"
            write (openFileID) "LOOKUP_TABLE default"
        ENDIF

        IF(vtkAscciFormat)THEN
            DO I = 1,NPoint
                write (openFileID,*) dataVal(I)
            ENDDO
        ELSE
            write (openFileID) (dataVal(I),I = 1,NPoint)
        ENDIF
    
      ELSEIF(DataType .EQ. 'ntubeData')THEN
        IF(vtkAscciFormat)THEN
            write (openFileID,*) " "
            IF(firstVtkType.EQ.1)THEN
                write (openFileID,'(A9,I15)') "CELL_DATA", NLink
            ENDIF
            write (openFileID,'(A8,A,A7)') "SCALARS ",dataName, " int"
            write (openFileID,'(A20)') "LOOKUP_TABLE default"
        ELSE
            buffer = ' '//lf    ; write(openFileID) trim(buffer)
            IF(firstVtkType.EQ.1)THEN
                buffer = 'CELL_DATA '//intToCharConvert(NLink)//lf  ; write(openFileID) trim(buffer)
            ENDIF
            buffer = 'SCALARS '//dataName//' '//'int '//lf  ; write(openFileID) trim(buffer)
            buffer = 'LOOKUP_TABLE default '//lf  ; write(openFileID) trim(buffer)
        ENDIF

        IF(vtkAscciFormat)THEN
            DO I = 1,NLink
                write (openFileID,*) dataVal(I)
            ENDDO
        ELSE
            !write (openFileID) (dataVal(I),I = 1,NLink)
            DO I = 1,NLink
                write (openFileID) dble(dataVal(I))
            ENDDO
        ENDIF
      ELSE
        WRITE(*         , '(2A)' ) ' !!!ERROR: WRONG DATA TYPE FOR VTK: ', DataType
    ENDIF

    !CLOSE(openFileID)

    END subroutine saveVtk_int

    
    
    
    subroutine saveVtk_real(openFileID, firstFlag, firstVtkType, CN,   Npoint, CONF,      NLink, DataType, dataVal, dataName)

   ! REAL*8,DIMENSION(:,:)::CN(3,NPoint),CELLS_vtk(NLink,3)
   ! Integer,DIMENSION(:,:)::CONF(2,NLink)
    REAL*8::CN(3,NPoint),CELLS_vtk(NLink,3)
    Integer::CONF(2,NLink)
    
    real*8:: dataVal(:)
    !DOUBLEPRECISION,DIMENSION(:),allocatable::dataVal
    !character(len=10) :: dataName
    CHARACTER(*):: dataName
    CHARACTER(*):: DataType
    INTEGER::firstFlag,I,J,numData,firstVtkType, openFileID, Npoint, Nlink
    LOGICAL:: vtkAscciFormat
    LOGICAL:: vtkPolydataFormat

    character :: buffer*100, lf*1, str1*8, str2*8, str3*8
    integer   :: ivtk = 9, int
                         
    lf = char(10) ! line feed character

    vtkAscciFormat = .TRUE.
    !vtkAscciFormat = .FALSE.
    vtkPolydataFormat =  .TRUE.
    !allocate(character(len=10) :: dataName)
    !ALLOCATE(dataVal(numData))
    CELLS_vtk(:,1)=2
    CELLS_vtk(:,2)=dble(CONF(1,:)-1)
    CELLS_vtk(:,3)=dble(CONF(2,:)-1)

    IF(firstFlag .EQ. 1)THEN
        IF(vtkAscciFormat)THEN
            IF(vtkPolydataFormat)THEN
                write (openFileID,'(A26)')"# vtk DataFile Version 3.0"
                write (openFileID,'(A12)')"# vtk output"
                write (openFileID,'(A5)')"ASCII"
                write (openFileID,'(A25)') " "
                write (openFileID,'(A16)') "DATASET POLYDATA"
                write (openFileID,'(A6,I15,A7)') "POINTS", NPoint, " double"
            ELSE
                write (openFileID,'(A26)')"# vtk DataFile Version 2.0"
                write (openFileID,'(A25)')"# SATURATION DISTRIBUTION"
                write (openFileID,'(A5)')"ASCII"
                write (openFileID,'(A25)') " "
                write (openFileID,'(A25)') "DATASET UNSTRUCTURED_GRID"
                write (openFileID,'(A6,I15,A7)') "POINTS", NPoint, " double"
            ENDIF
        ELSE
            CLOSE(openFileID)
            open(unit=openFileID,file='vtk_poreInlet_initial.vtk',form='binary',convert='BIG_ENDIAN')
            !open(unit=openFileID,file=savePath//'vtk_poreInlet_initial.vtk',access='stream')
            !open(unit=openFileID,file=savePath//'vtk_poreInlet_initial.vtk',form='binary')
            !         open(unit=openFileID,file=savePath//'vtk_poreInlet_initial.vtk', form='unformatted',&
            !recordtype='STREAM_LF',&
            !action='write'                          ,&
            !convert='BIG_ENDIAN'                    ,&
            !access='sequential'                      )

            buffer = '# vtk DataFile Version 3.0'//lf  ; write(openFileID) trim(buffer)
            buffer = 'vtk output'//lf                  ; write(openFileID) trim(buffer)
            buffer = 'BINARY'//lf                      ; write(openFileID) trim(buffer)
            buffer = 'DATASET UNSTRUCTURED_GRID'//lf    ; write(openFileID) trim(buffer)
            buffer = 'POINTS '//intToCharConvert(NPoint)//' double'//lf  ; write(openFileID) trim(buffer)
            !CLOSE(openFileID)
        ENDIF

        IF(vtkAscciFormat)THEN
            DO I = 1,NPoint
                DO J = 1,3
                    !set very small number to zero
                    IF( ABS( CN(J,I) ) .LT. 1d-20)THEN
                        CN(J,I) = 0D0
                    ENDIF
                ENDDO
                write (openFileID,'(3G20.7)') CN(:,I)
            ENDDO
        ELSE
            DO I = 1,NPoint
                DO J = 1,3
                    IF( ABS( CN(J,I) ).LT.1D-20)THEN
                        CN(J,I) = 0D0
                    ENDIF
                ENDDO
                write (openFileID) CN(:,I)
            ENDDO

            ! write (openFileID) (CN(:,I), I=1,NPoint)
            ! write (openFileID) (CELLS_vtk(I,1), I=1,10), lf
            ! write (openFileID) (CN(:,I),real(0.,kind=8) , I=1,NPoint)
            ! CELLS_vtk(1:NPoint,1) = CN(1,:)
            ! write(openFileID)(CELLS_vtk(I,1) , i=1,NPoint)
        ENDIF
        !CLOSE(openFileID)
        !write (openFileID,*) " "
        !write (openFileID,'(A5,2I10)') "CELLS", NLink, 3*NLink


        if(vtkAscciFormat)then
            IF(vtkPolydataFormat)THEN
                write (openFileID,*) " "
                write (openFileID,'(A5,2I15)') "LINES", NLink, 3*NLink    
            ELSE
                write (openFileID,*) " "
                write (openFileID,'(A5,2I15)') "CELLS", NLink, 3*NLink
            ENDIF
        else
            buffer = ' '//lf    ; write(openFileID) trim(buffer)
            buffer = 'CELLS '//intToCharConvert(NLink)//' '//intToCharConvert(3*NLink)//lf  ; write(openFileID) trim(buffer)
        endif

        
        IF(vtkAscciFormat)THEN
            DO I = 1,NLink
                write (openFileID,'(A2,2I15)') "2 ", CONF(1,I)-1, CONF(2,I)-1
            ENDDO
        ELSE
            !write (openFileID) ( 2+CONF(1,I)*0, I = 1,NLink), ( CONF(1,I)-1, I = 1,NLink), (CONF(2,I)-1, I = 1,NLink)
            ! write (openFileID) (CELLS_vtk(I,:), I=1,NLink)
            !write (openFileID) (CN(:,I), I=1,NPoint)
            DO I = 1,NLink
                write (openFileID) CELLS_vtk(I,:)
            ENDDO
        ENDIF

        !CLOSE(openFileID)
        IF(vtkPolydataFormat .EQ. .FALSE.)THEN
            ! write (openFileID,*) " "
            !write (openFileID,'(A10,I10)') "CELL_TYPES", NLink
            IF(vtkAscciFormat)THEN
                write (openFileID,*) " "
                write (openFileID,'(A10,I15)') "CELL_TYPES", NLink
            ELSE
                buffer = ' '//lf    ; write(openFileID) trim(buffer)
                buffer = 'CELL_TYPES '//intToCharConvert(NLink)//lf  ; write(openFileID) trim(buffer)
            ENDIF

            IF(vtkAscciFormat)THEN
                DO I = 1,NLink
                    write (openFileID,'(A1)') "3"
                ENDDO
            ELSE
                CELLS_vtk(:,1)=3
                !write (openFileID) (3 + CONF(1,I)*0,I = 1,NLink)
                ! write (openFileID) CELLS_vtk(:,1)
                DO I = 1,NLink
                    write (openFileID) CELLS_vtk(I,1)
                ENDDO
            ENDIF
        ENDIF
    ENDIF

    !CLOSE(openFileID)
                        
    IF(DataType .EQ. 'nnodeData')THEN
        IF(vtkAscciFormat)THEN
        write (openFileID,*) " "
        !add porebody properties
        IF(firstVtkType.EQ.1)THEN
            write (openFileID,'(A10,I15)') "POINT_DATA", NPoint
        ENDIF
        write (openFileID,'(A8,A,A7)') "SCALARS ",dataName,  " double"
        write (openFileID,'(A20)') "LOOKUP_TABLE default"
        ELSE
        !IF(0)THEN
            write (openFileID) " "
            !add porebody properties
            IF(firstVtkType.EQ.1)THEN
                write (openFileID) "POINT_DATA", NPoint
            ENDIF
            write (openFileID) "SCALARS ",dataName,  " double"
            write (openFileID) "LOOKUP_TABLE default"
        ENDIF

        IF(vtkAscciFormat)THEN
            DO I = 1,NPoint
                write (openFileID,*) dataVal(I)
            ENDDO
        ELSE
            write (openFileID) (dataVal(I),I = 1,NPoint)
        ENDIF
    
      ELSEIF(DataType .EQ. 'ntubeData')THEN
        IF(vtkAscciFormat)THEN
            write (openFileID,*) " "
            IF(firstVtkType.EQ.1)THEN
                write (openFileID,'(A9,I15)') "CELL_DATA", NLink
            ENDIF
            write (openFileID,'(A8,A,A7)') "SCALARS ",dataName, " int"
            write (openFileID,'(A20)') "LOOKUP_TABLE default"
        ELSE
            buffer = ' '//lf    ; write(openFileID) trim(buffer)
            IF(firstVtkType.EQ.1)THEN
                buffer = 'CELL_DATA '//intToCharConvert(NLink)//lf  ; write(openFileID) trim(buffer)
            ENDIF
            buffer = 'SCALARS '//dataName//' '//'int '//lf  ; write(openFileID) trim(buffer)
            buffer = 'LOOKUP_TABLE default '//lf  ; write(openFileID) trim(buffer)
        ENDIF

        IF(vtkAscciFormat)THEN
            DO I = 1,NLink
                write (openFileID,*) dataVal(I)
            ENDDO
        ELSE
            !write (openFileID) (dataVal(I),I = 1,NLink)
            DO I = 1,NLink
                write (openFileID) dble(dataVal(I))
            ENDDO
        ENDIF
      ELSE
        WRITE(*         , '(2A)' ) ' !!!ERROR: WRONG DATA TYPE FOR VTK: ', DataType
    ENDIF

    !CLOSE(openFileID)

    END subroutine saveVtk_real
  
    
    
    FUNCTION intToCharConvert(intNum) RESULT (intToCharConvertWasDone)
    IMPLICIT NONE
    INTEGER:: intNum
    character(len=15) :: intToCharConvertWasDone
    !write (intToCharConvertWasDone, '(G15.5)') intNum
    write (intToCharConvertWasDone, '(I10)') intNum
    intToCharConvertWasDone = adjustl(intToCharConvertWasDone)
    intToCharConvertWasDone = trim(intToCharConvertWasDone)
    END FUNCTION intToCharConvert
    
    end MODULE saveData