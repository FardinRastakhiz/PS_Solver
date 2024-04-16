
    Module solveMatrix
    use, intrinsic :: iso_fortran_env, only : wp => real64
    use iso_c_binding
    use matrixBuilder
    implicit none


    !!!note: in the project folder i have put two matrices for examples. a larger size (with name append of _L) and a medium size (with name append of _M).

    !medium size
    !character(len=30) :: matrixFile = "Matrix_M", bFile = "Pknown_M"
    !integer, parameter :: n = 11720, nnz = 27839     ! medium size (with name append of _M)

    !larger size
    !  integer, parameter :: n = 997350, nnz = 26218238     !larger size (with name append of _L)
    !   character(len=30) :: matrixFile = "Matrix_L", bFile = "Pknown_L"


    integer, allocatable :: rowIndex(:), colIndex(:)
    real(c_double), allocatable :: Avalues(:), b(:), x(:)
    integer :: i, j, error
    integer, allocatable :: ia(:), ja(:)

    !parameters
    !
    !general
    doubleprecision, parameter :: pi = 3.14159265359
    !fluid
    doubleprecision, parameter :: mu_w = 0.0010e4/60.0_8 !Pa.min

    contains

    !subroutine readSparseMatrix(filename, nnz, rowIndex, colIndex, Avalues)
    !!this subroutine reads 'Matrix' from the file and saves it in rowIndex, colIndex, and Avalues vectors
    !character(len=*), intent(in) :: filename
    !integer, intent(in) :: nnz
    !integer :: rowIndex(nnz), colIndex(nnz)
    !real(wp) :: Avalues(nnz)
    !integer :: i, unit
    !unit  = 1001
    !open(unit, file=filename, action="read", status="old")
    !do i = 1, nnz
    !    read(unit, *) rowIndex(i), colIndex(i), Avalues(i)
    !end do
    !close(unit)
    !end subroutine

    subroutine readVector(filename, n, b)
    !this subroutine reads 'Pknown' from the file and saves it in 'b' vector
    character(len=*), intent(in) :: filename
    integer, intent(in) :: n
    doubleprecision, intent(out) :: b(n)
    integer :: i, unit
    open(newunit=unit, file=filename, action="read", status="old")
    do i = 1, n
        read(unit, *) b(i)
    end do
    close(unit)
    end subroutine

    subroutine readVectorInt (filename, n, b)
    !this subroutine reads 'Pknown' from the file and saves it in 'b' vector
    character(len=*), intent(in) :: filename
    integer, intent(in) :: n
    integer, intent(out) :: b(n)
    integer :: i, unit
    open(newunit=unit, file=filename, action="read", status="old")
    do i = 1, n
        read(unit, *) b(i)
    end do
    close(unit)
    end subroutine

    end Module solveMatrix


    program Main
    use iso_c_binding
    use solveMatrix
    use saveData
    implicit none

    integer :: idx, errCode
    integer :: num_threads
    integer(c_int) :: status
    type(c_ptr) :: dll_handle
    !read network info
    integer :: nnode,ntube, ios
    character(200) :: line
    character(len=20) :: filename
    integer :: iunit

    !pore locations
    doubleprecision,dimension(:,:),allocatable:: pore_loc
    !pore connections
    integer,dimension(:,:),allocatable:: conf
    !pore sizes
    doubleprecision,dimension(:),allocatable:: porer,piper,pipel
    !pressure
    doubleprecision,dimension(:),allocatable:: edgeCond
    !boundary pores
    integer,dimension(:),allocatable:: porebnd

    !re-indexing
    integer,dimension(:),allocatable:: index_node_in_act, index_edg_in_origin, porebnd_act
    integer,dimension(:,:),allocatable:: conf_act
    doubleprecision,dimension(:,:),allocatable:: pore_loc_act
    integer :: count, new_count, nnode_act, ntube_act
    logical :: include_connection

    !solving pressure
    integer :: nnz
    integer,dimension(:),allocatable:: ai, aj
    doubleprecision,dimension(:),allocatable:: pr_act, bc, rhs
    doubleprecision:: pr_inlet, pr_outlet

    !calucing Q
    doubleprecision,dimension(:),allocatable:: p_origin, delp, q
    doubleprecision:: Q_tot_inlet, Q_tot_outlet

    !calucing concentration
    logical::  PASSTO
    integer :: inlet_edge, outlet_edge,II
    doubleprecision,dimension(:),allocatable:: AREAV, VOLUMEV, POREV, Q_tot_node
    doubleprecision,dimension(:),allocatable:: COEFB,COEF1,COEF2,COEF3,COEF4
    doubleprecision,dimension(:),allocatable:: conc_node,  conc_edje
    integer,dimension(:),allocatable::  UPnode,DownNODE
    doubleprecision:: DT, Dmol, lDiff, conc_inlet, DtPore, DtPipe, DtDiff




    ! DSDBCG solver
    INTEGER  ISYM,ITOL,ITMAX,NU,ITER,IERR,LENW,LENIW !DECLARATION
    doubleprecision TOL, ERR
    integer,allocatable::IWORK(:)
    doubleprecision,dimension(:),allocatable::RWORK

    !save as vtk
    INTEGER::firstFlag,firstVtkType

        ! c++ solvers parameters
    ! below codes can be defined somewhere else
    integer :: solverLibrary = 2 ! acceptable values : 1 = PETSC_CPU , 2 = PETSC_GPU , 3 = ViennaCL_GPU 
    integer :: useOpenMp = 1 ! acceptable values : 0 = dont use openMP , 1 = use openMP
    integer :: numOfThreads = 8 ! acceptable values : any positive integer value
    integer :: platform = 2 ! acceptable values should be extracted from ses_get_devices function
    integer :: device = 0 ! acceptable values should be extracted from ses_get_devices function
    integer :: iterationCount = -1 ! acceptable values : -1 = default , any positive integer value
    real :: precision = -1 ! acceptable values : -1 = default , any positive double value
    integer :: preconditioner = 7 ! acceptable values : 1 = jacobi, 2 = ILU, 3 = gasm, 4 = icc, 5 = ksp, 6 = bjacobi, 7 = sor, 8 = asm, 9 = cholesky
    integer :: returnValue
    type(c_ptr) :: solver_pointer
    
    ! c++ solvers interfaces
    interface
        function ses_build_initial_guess(numRows, numRowsAct, locX , locY , locZ, locActX , locActY , locActZ, bnd, x) result(ok) bind(C, name="ses_build_initial_guess")
            use iso_c_binding
            integer(c_int) :: ok
            integer(c_int), value :: numRows , numRowsAct
            integer(c_int), dimension(*) :: bnd
            doubleprecision, dimension(*) :: locX , locY , locZ , locActX , locActY , locActZ
            doubleprecision, dimension(*), intent(out) :: x
        end function ses_build_initial_guess
    end interface
    interface
        function ses_solve_pressure_cpu(numRows, numNonzero, rowIndices, colIndices, values, b, x, iterationCount, precision, useOpenMp, numOfThreads, preconditioner) bind(C, name="ses_solve_pressure_cpu")
            use iso_c_binding
            integer(c_int), value :: numRows, numNonzero, iterationCount, useOpenMp, numOfThreads, preconditioner
            integer(c_int), dimension(*) :: rowIndices, colIndices
            real(c_double), dimension(*) :: values, b
            real(c_double), dimension(*), intent(out) :: x
            real(c_double), value :: precision
            type(c_ptr) :: ses_solve_pressure_cpu
        end function ses_solve_pressure_cpu
    end interface
    
    interface
        function ses_solve_pressure_gpu(numRows, numNonzero, rowIndices, colIndices, values, b, x, solverLibrary, iterationCount, precision, platform , device, preconditioner) bind(C, name="ses_solve_pressure_gpu")
            use iso_c_binding
            integer(c_int), value :: numRows, numNonzero, solverLibrary, iterationCount, platform, device, preconditioner
            integer(c_int), dimension(*) :: rowIndices, colIndices
            real(c_double), dimension(*) :: values, b
            real(c_double), dimension(*) , intent(out) :: x
            real(c_double), value :: precision
            type(c_ptr) :: ses_solve_pressure_gpu
        end function ses_solve_pressure_gpu
    end interface
    
    interface
        function ses_solve_begin_density_cpu(numRows, numNonzero, rowIndices, colIndices, values, b, x, iterationCount, precision, useOpenMp, numOfThreads, preconditioner) bind(C, name="ses_solve_begin_density_cpu")
            use iso_c_binding
            integer(c_int), value :: numRows, numNonzero, iterationCount, useOpenMp, numOfThreads, preconditioner
            integer(c_int), dimension(*) :: rowIndices, colIndices
            real(c_double), dimension(*) :: values, b
            real(c_double), dimension(*), intent(out) :: x
            real(c_double), value :: precision
            type(c_ptr) :: ses_solve_begin_density_cpu
        end function ses_solve_begin_density_cpu
    end interface
    
    interface
        function ses_solve_begin_density_gpu(numRows, numNonzero, rowIndices, colIndices, values, b, x, solverLibrary, iterationCount, precision, platform , device, preconditioner) bind(C, name="ses_solve_begin_density_gpu")
            use iso_c_binding
            integer(c_int), value :: numRows, numNonzero, solverLibrary, iterationCount, platform, device, preconditioner
            integer(c_int), dimension(*) :: rowIndices, colIndices
            real(c_double), dimension(*) :: values, b
            real(c_double), dimension(*) , intent(out) :: x
            real(c_double), value :: precision
            type(c_ptr) :: ses_solve_begin_density_gpu
        end function ses_solve_begin_density_gpu
    end interface
    
    interface
        function ses_solve_next(solver_container, b, x, iterationCount, precision) bind(C, name="ses_solve_next")
            use iso_c_binding
            integer(c_int), value :: iterationCount
            real(c_double), dimension(*) :: b
            real(c_double), dimension(*) , intent(out) :: x
            real(c_double), value :: precision
            type(c_ptr), value :: solver_container
            integer(c_int) :: ok
        end function ses_solve_next
    end interface
    
    
    interface
        function ses_write_devices_to_file() result(ok) bind(C, name="ses_write_devices_to_file")
            use iso_c_binding
            integer(c_int) :: ok
        end function ses_write_devices_to_file
    end interface
    
    ! Allocate rowIndex, colIndex, Avalues, b, x, ia, and ja arrays
    !    allocate(rowIndex(nnz), colIndex(nnz), Avalues(nnz), b(n), x(n), ia(n+1), ja(nnz))

    ! Read the sparse matrix A and vector b from the text files
    !   call readSparseMatrix(matrixFile, nnz, rowIndex, colIndex, Avalues)
    !   call readVector(bFile, n, b)

    !   x = solve_matrix(n, nnz, rowIndex, colIndex, Avalues, b, x)
    ! Wait for user input before exiting
    !   write(*,*) "Press Enter to exit"
    !    read(*,*)
    !
    !
    !get network size:

    ! open the file


    !nnode
    !
    !
    !number of rows in the 'PORE_LOC.txt' is qual to number of nodes (i.e., nnode)
    !this is location of each node ( each row is x,y,z coordinate of one node). so, numebr of lines in this file
    !is eqaul to the number of nodes (i.e., nnode)

    !first, open the 'PORE_LOC.txt' and read the number of files to obtain 'nnode', needed to allocated 'pore_loc' variable
    iunit=100
    open(unit=iunit, file='PORE_LOC.txt', status='old')
    nnode = 0
    do  !read through the 'PORE_LOC.txt' file line by line
        read(iunit,'(A)',iostat=ios) line
        if (ios /= 0) exit  ! exit the loop at end of file
        nnode = nnode + 1
    end do
    close(iunit)   ! close the file

    !now read the content of pore locations
    allocate(pore_loc(3,nnode))
    open(unit=iunit, file='PORE_LOC.txt', status='old')
    do i = 1, nnode  !read through the 'PORE_LOC.txt' file line by line
        read(iunit,*) pore_loc(:,I) !pore_loc(:,i) is x,y,z of ith node
    end do
    close(iunit)   ! close the file
    print *, 'number of pore bodies (nodes) = ', nnode

    !ntune
    !
    !
    !number of rows in the 'CONF.txt' is qual to number of edges/tubes/pore throats (i.e. ntube)
    !each line of 'CONF.txt' shows the two nodes that are connedcted together
    !first, open the 'CONF.txt' and read the number of files to obtain 'ntube', needed to allocated 'conf' variable
    iunit=101
    open(unit=iunit, file='CONF.txt', status='old')
    ntube = 0
    do !read through the 'PORE_LOC.txt' file line by line
        read(iunit,'(A)',iostat=ios) line
        if (ios /= 0) exit  ! exit the loop at end of file
        ntube = ntube + 1
    end do
    close(iunit)   ! close the file
    print *, 'number of pore throats (ntube) = ', ntube

    !now read content of pore connections (i.e., variable 'conf')
    allocate(conf(2,ntube))
    open(unit=iunit, file='CONF.txt', status='old')
    do i = 1, ntube  !read through the 'PORE_LOC.txt' file line by line
        read(iunit,*) conf(:,I)
    end do
    close(iunit)   ! close the file


    !allocate pore sizes
    Allocate(porer(nnode),piper(ntube),pipel(ntube))

    !read pore sizes
    !
    !
    call readVector('PORER.TXT', nnode, porer) ! for eahc node (i.e,. pore) we have one radius.
    call readVector('PIPER.TXT', ntube, piper)
    call readVector('PIPEL.TXT', ntube, pipel)
    porer = porer / 1000.0d0  !um to mm
    piper = piper / 1000.0d0  !um to mm
    pipel = pipel / 1000.0d0  !um to mm

    !check negative values
    do i = 1, ntube
        if (piper(i) .lt. 1e-10) then
            print*, 'small or nagative piper', i, piper(i)
            piper(i) = 1.0
        end if
        if (pipel(i) .lt. 1e-10) then
            print*, 'small or nagative pipel', i, pipel(i)
            pipel(i) = 1.0
        end if
    end do
    porer = porer / 1000.0d0  !um to mm

    do i = 1, nnode
        if (porer(i) .lt. 1e-10) then
            print*, 'small or nagative porer', i, porer(i)
            porer(i) = 1.0
        end if
    end do


    !at this point the network geometry is completed.

    !now we need to assign the boudaries (inlet and outlet). this data will be saved in 'porebnd' varible
    !Flags in portebnd: 1 is inlet, 2 is outlet, and 0 is a normal (active) node that we want to calculate its pressure
    !allocate boudary pores
    Allocate(porebnd(nnode))
    call readVectorInt('PORE_INLET.txt', nnode, porebnd)

    !now proceed with preprartion for solving for pressure

    !First, determine the size of new_nodes by counting
    count = 0
    do i = 1, ntube
        include_connection = porebnd(conf(1, i)) == 0 .and. porebnd(conf(2, i)) == 0
        if (include_connection) then
            count = count + 1
        end if
    end do

    ! Now, allocate the new_nodes array based on the determined size
    allocate(conf_act(2, count), index_node_in_act(nnode), index_edg_in_origin(count))

    ! Populate the new_nodes array
    idx = 0
    new_count = 0
    index_node_in_act = 0
    index_edg_in_origin = 0
    !do i = 1, ntube
    !    include_connection = porebnd(conf(1, i)) == 0 .and. porebnd(conf(2, i)) == 0
    !    if (include_connection) then
    !        new_count = new_count + 1
    !        index_edg_in_origin(new_count) = i
    !        do j = 1, 2
    !            if (index_node_in_act(conf(j, i)) == 0) then
    !                idx = idx + 1
    !                index_node_in_act(conf(j, i)) = idx
    !            end if
    !            conf_act(j, new_count) = index_node_in_act(conf(j, i))
    !        end do
    !    end if
    !end do

    do i = 1, ntube
        include_connection = porebnd(conf(1, i)) == 0
        if (include_connection) then
            if (index_node_in_act(conf(1, i)) == 0) then
                idx = idx + 1
                index_node_in_act(conf(1, i)) = idx
            end if
        end if

        include_connection = porebnd(conf(2, i)) == 0
        if (include_connection) then
            if (index_node_in_act(conf(2, i)) == 0) then
                idx = idx + 1
                index_node_in_act(conf(2, i)) = idx
            end if
        end if

        !add to  conf_act if both end nodes are active
        include_connection = porebnd(conf(1, i)) == 0 .and. porebnd(conf(2, i)) == 0
        if (include_connection) then
            new_count = new_count + 1
            index_edg_in_origin(new_count) = i
            conf_act(:, new_count) = index_node_in_act(conf(:, i))
        end if

    end do



    nnode_act = idx
    ntube_act = new_count

    ! Print new_nodes array to check
    do i = 1, new_count
        !      print *, "new_nodes(", i, ", :) = ", conf_act(1, i), conf_act(2, i)
    end do

    allocate(pore_loc_act(3,nnode_act), porebnd_act(nnode_act))
    j = 0
    do i = 1, nnode
        if(index_node_in_act(i) .gt. 0)then
            j = j +1
            pore_loc_act(:,index_node_in_act(i)) = pore_loc (:,i)
            porebnd_act(index_node_in_act(i))    = porebnd(i)
        endif
    end do

    !lets test this by visualizing the two graphs
    open(unit=10, file='myVtkOriginal.vtk')
    firstFlag  =  1
    firstVtkType =   1
    call saveVtk_int(10, firstFlag, firstVtkType, pore_loc,   nnode, CONF, ntube, 'nnodeData', porebnd, 'porebnd')
    close(10)

    open(unit=10, file='myVtkFinal.vtk')  !mena with _act
    firstFlag  =  1
    firstVtkType =   1
    call saveVtk_int(10, firstFlag, firstVtkType, pore_loc_act,   nnode_act, conf_act, ntube_act, 'nnodeData', porebnd_act, 'porebnd_act')
    close(10)

    !prepare to calculate pressure
    !
    !calculate conductance
    Allocate(edgeCond(ntube))
    edgeCond = 0d0
    edgeCond = Pi * piper**4 / (8.0D0 * mu_w * pipel)
    do i = 1, ntube
        if (edgeCond(i) .lt. 1e-10) then
            print*, 'small or nagative edgeCond', i, edgeCond(i)
            edgeCond(i) = 0.0d0
        end if
    end do


    !set BCs:
    pr_inlet = 1.0d0
    pr_outlet = 0.d0
    !nnz
    nnz = nnode_act + ntube_act

    allocate(rhs(nnode_act), ai(nnz),aj(nnz),bc(nnz), pr_act(nnode_act) )
    bc=0.D0
    ai=0
    aj=0
    pr_act = 0 !to be solved for

    !add the -sum of throat conductances to the main diagonal
    !
    !add index of main digonals (to be used to add -sums to the corresponding BC)
    DO I=1,nnode_act
        ai(I) = I
        aj(I) = I
    ENDDO

    !PUT (-SUMMATION) ON DIAGONAL. FILL MNODE POSITIONS
    DO I=1,Ntube !EFECT OF RP ON RHS
        IF( poreBnd( CONF(1,i)) .EQ. 0 )THEN
            bc(  index_node_in_act( CONF(1,I))  ) = bc(  index_node_in_act( CONF(1,I))  )  - edgeCond(I)
        ENDIF

        IF( poreBnd( CONF(2,i)) .EQ. 0 )THEN
            bc(  index_node_in_act( CONF(2,I))  ) = bc(  index_node_in_act( CONF(2,I))  )  - edgeCond(I)
        ENDIF
    ENDDO


    count = nnode_act  !the numbber of filled positions so far

    !A * x = b
    !in our code the corresponding varibles are:
    !Bc * pr_act = rhs

    !NOW PUT INDIVIDUAL CONECTIONS: 2*NTUBE-Ltb-RT
    !
    !SYMMETRIC MATRIX: WE NEED ONLY PUT cond FOR one direction (e.g., only in FORWARD numbering CONNECTIONS)
    do i=1,ntube_act  !active tubes in Ai and Aj (not sum terms)
        count = count + 1
        Ai(count) = conf_act(2,i)
        Aj(count) = conf_act(1,i)
        BC(count) = edgeCond( index_edg_in_origin(i) )
    enddo

    !fill the rhs
    rhs = 0.0D0
    DO I=1,Ntube !EFECT OF RP ON RHS

        if(poreBnd(CONF(1,I)) .gt. 0 .and. poreBnd(CONF(2,I)) .gt. 0)then
            !this is a tube that both sides of it are bnd pores (e.g., both are inlet or outlet). so we dont need to do anything
            !do nothing
        else
            IF(poreBnd(CONF(1,I)).eq.1)THEN
                !node "CONF(1,I)"is boundary node, so in RHS add "-pr_inlet*COND(I)" to the node connected to each, which is "CONF(2,I)"
                rhs(  index_node_in_act( CONF(2,I))  ) = rhs(  index_node_in_act( CONF(2,I))  )  - pr_inlet*edgeCond(I)
            ELSEIF(poreBnd(CONF(1,I)).eq.2)THEN
                rhs(  index_node_in_act( CONF(2,I))  ) = rhs(  index_node_in_act( CONF(2,I))  )  - pr_outlet*edgeCond(I)
            ENDIF

            IF(poreBnd(CONF(2,I)).eq.1)THEN
                rhs(  index_node_in_act( CONF(1,I))  ) = rhs(  index_node_in_act( CONF(1,I))  )  - pr_inlet*edgeCond(I)
            ELSEIF(poreBnd(CONF(2,I)).eq.2)THEN
                rhs(  index_node_in_act( CONF(1,I))  ) = rhs(  index_node_in_act( CONF(1,I))  )  - pr_outlet*edgeCond(I)
            ENDIF
        endif

    ENDDO


    ISYM = 1  !is symmetric matrix (1 = yes)
    ITOL=1
    ITMAX=5000
    ITER=30000+1;IERR=0;IUNIT=0
    NU=nnz-nnode_act !Non-Zeros in the upper + lower triangular matrix
    LENW=1+nnode_act*(60+7)+60*(60+3)+NU !LENW >= 8*N.
    LENW=10*nnode_act
    LENIW=NU+4*nnode_act+32 !LENIW >= 10
    LENIW=10*nnode_act
    allocate(RWORK(LENW),IWORK(LENIW))
    RWORK=0.1;IWORK=1
    ERR=0.00001
    TOL = 10e-18 !tollerance
    !DSDBCG is the Diagonally Scaled BiConjugate Gradient Sparse Ax=b Solver.
    ! nnode_act is number of active nodes (i.e., n in standard Ax=b)
    !pr_act is the unknown pressure (i.e., x in Ax=b)
    !ai,aj,bc is sparse represention of coefficient matrxi (sparese of A in Ax=b)
    !ai = row index, aj=col index, bc = values
    !call DSDBCG (nnode_act,rhs,pr_act,nnz,ai,aj,bc,ISYM,ITOL, &  !MAIN SOLVER
    !    TOL,ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW)
    !IF(IERR.EQ.0)THEN
    !    WRITE(*  ,'(A45)')," solverReturnFlag: IERR = 0 => All went well."
    !else
    !    WRITE(*  ,'(A45)')," solverReturnFlag: IERR /= 0 => Not Succeed!!"
    !ENDIF
    
    
    if(solverLibrary.eq.1)THEN
        solver_pointer = ses_solve_pressure_cpu(nnode_act, nnz, ai, aj, bc, rhs, pr_act, iterationCount, precision, useOpenMp, numOfThreads, preconditioner)
    endif
    
    ! solve using gpu
    if(solverLibrary.eq.2 .or. solverLibrary.eq.3)THEN
        ! here first we need to take platforms and devices
        returnValue = ses_write_devices_to_file()
        ! choose your platform and device by index - indices start from 0
        solver_pointer = ses_solve_pressure_gpu(nnode_act, nnz, ai, aj, bc, rhs, pr_act, solverLibrary, iterationCount, precision, platform , device, preconditioner)
        
    endif
    WRITE(*, *)," Min pr= ", minval(pr_act), " Min pr= ", maxval(pr_act)


    !calculate fluxes
    !
    !first make the full pressure vector (i.e., add boundary pressures)
    allocate(p_origin(nnode))
    p_origin=-999
    J=0
    DO I = 1,NNODE
        IF(poreBnd( I ) .EQ. 0)THEN
            J=J+1
            p_origin(I) = pr_act(J)
        ELSEIF(poreBnd( I ) .EQ. 1)THEN !This is a inlet boundary node
            p_origin(I) = pr_inlet
        ELSE
            p_origin(I) = pr_outlet
        ENDIF
    ENDDO

    !caclulate pressure difference for each edge
    allocate(delp(Ntube),q(Ntube))
    delp=0.0
    q=0.0D0
    DO I = 1,NTUBE
        delp(I)=DABS(   p_origin(CONF(1,I)) - p_origin(CONF(2,I))  )
    ENDDO
    !  delp(:)=DABS(   p_origin(CONF(1,:)) - p_origin(CONF(2,:))  )


    !calculate volumetric flux:
    q(:) = delp(:)  *  edgeCond(:)

    !calculte total inflow, and total outflow. in theory these two should be the same
    Q_tot_inlet = 0d0
    Q_tot_outlet = 0d0
    DO I = 1,NTUBE
        if(poreBnd( CONF(1,I) ) .EQ. 1 .and. poreBnd( CONF(2,I) ) .EQ. 2)then
            WRITE(* ,'(A)'), " error!!!! "; pause
        endif
        if(poreBnd( CONF(2,I) ) .EQ. 1 .and. poreBnd( CONF(1,I) ) .EQ. 2)then
            WRITE(* ,'(A)'), " error!!!! "; pause
        endif

        if(poreBnd( CONF(2,I) ) .EQ. 1 .and. poreBnd( CONF(1,I) ) .EQ. 1)then
            WRITE(* ,'(A)'), " INLET-INLET!!!! "; !pause
        endif
        if(poreBnd( CONF(2,I) ) .EQ. 2 .and. poreBnd( CONF(1,I) ) .EQ. 2)then
            WRITE(* ,'(A)'), " OUTLET-OUTLET!!!! "; !pause
        endif

        IF(  poreBnd( CONF(1,I) ) .EQ. 1 .OR. poreBnd( CONF(2,I) ) .EQ. 1)THEN
            Q_tot_inlet = Q_tot_inlet + q(I)
        ELSEIF( poreBnd( CONF(1,I) ) .EQ. 2 .OR. poreBnd(CONF(2,I) ) .EQ. 2)THEN
            Q_tot_outlet = Q_tot_outlet + q(I)
        ENDIF
    ENDDO
    WRITE(* ,'(A, E)'), " Total inflow  = ", ( Q_tot_inlet )
    WRITE(* ,'(A, E)'), " Total outflow = ", ( Q_tot_outlet )

    WRITE(* ,'(A, E, A)'), " Percentage Error (using: inflow - outflow) = ", ( abs(Q_tot_inlet-Q_tot_outlet)/abs(Q_tot_inlet) ), '%'

    !SOLUTE TRANSPORT LOOP

    !rhs
    !
    !
    deallocate(rhs, ai,aj,bc )
    allocate(rhs(nnode))
    rhs = 0d0
    allocate(conc_node(nnode), conc_edje(Ntube))
    conc_node = 0d0
    conc_edje = 0d0
    allocate(AREAV(Ntube), VOLUMEV(Ntube),Q_tot_node(nnode))
    allocate(COEFB(Ntube),COEF1(Ntube),COEF2(Ntube))

    DT = 0.000000001
    Dmol = 1e-6
    lDiff = 0.5d0
    conc_inlet = 1.0d0

    AREAV = pi * piper**2
    VOLUMEV = pi * piper**2 * pipel
    POREV = 4.0d0/3.0 * pi * porer**3

    !needed by the loop
    !dt
    Q_tot_node = 0d0
    DO I=1,NTUBE
        IF(poreBnd(CONF(1,I)).EQ.1)THEN !so CONF(1,I)) is the upstream pore
            Q_tot_node(CONF(1,I))=Q_tot_node(CONF(1,I))+ Q(I)
        ENDIF
        IF(poreBnd(CONF(2,I)).EQ.1)THEN
            Q_tot_node(CONF(2,I))=Q_tot_node(CONF(2,I))+ Q(I)
        ENDIF
    ENDDO
    !now we fo with all other pores
    DO I=1,NTUBE
        !add Qtot to the downstream pore (ie., pore with lower pressure)
        IF(p_origin(CONF(1,I)) .GT. p_origin(CONF(2,I)))THEN
            Q_tot_node(CONF(2,I))= &
                Q_tot_node(CONF(2,I))+(q(I))
        ELSE
            Q_tot_node(CONF(1,I))= &
                Q_tot_node(CONF(1,I))+(q(I))
        ENDIF
    ENDDO
    
    DtPore=MINVAL(POREV/Q_tot_node)*0.5 ![min]
    DtPipe=MINVAL(VOLUMEV/Q)*0.5
    DtDiff=0.5D0*lDiff**2/Dmol
    Dt=DtPore
    IF(DtPipe.LT.Dt)THEN
        Dt=DtPipe
    ENDIF
    IF(DtDiff.LT.Dt)THEN
        Dt=DtDiff
    ENDIF
    Dt = 100 * Dt
    IF(Dt .LT. 0)THEN
        PRINT*,"ERROR! NEGATIVE Dt. Dt = ",Dt
    ENDIF
    WRITE(* ,'(A, E)'), " Dt = ", ( Dt )
    
    COEFB=0.;COEF1=0.;COEF2=0.
    DO I=1,NTUBE !FIIL TWO COEF. THAT HAVE DIMENTION OF NTUBE
        COEFB(I)=1.0D0+ (Q(I)*DT/VOLUMEV(I))+ &
            + (Dmol*AREAV(I)*2.D0*DT / (VOLUMEV(I)*lDiff))

        COEF1(I)=   ( Q(I)*DT ) &
            /(VOLUMEV(I)*COEFB(I))

        COEF2(I)=(Dmol*AREAV(I)*DT) / (VOLUMEV(I)*lDiff*COEFB(I))
    ENDDO

    allocate(UPNODE(Ntube), DownNODE(Ntube))
    !UPNODE  (NTUBE) DETERMINES THE UPSTREAM NODE IN EACH TUBE
    DO I=1,NTUBE
        IF (p_origin(CONF(1,I)).GT.p_origin(CONF(2,I)) )THEN
            UPNODE(I)  =CONF(1,I)
            DownNODE(I)=CONF(2,I)
        ELSE
            UPNODE(I)  =CONF(2,I)
            DownNODE(I)=CONF(1,I)
        ENDIF
    ENDDO

    nnz = 0
    call soluteCoefMatrix(AI,AJ,BC,nnode,ntube,poreBnd,CONF,Q,POREV,Dmol, VOLUMEV,AREAV, DT, p_origin,Q_tot_node, nnz)

    !solver setting
    open(11,FILE='error-record.dat')
    ISYM=0; ITOL=1;ITMAX=1000
    ITER=5000+1;IERR=0;IUNIT=0
    LENW=1+(NNODE)*(60+7)+60*(60+3)+nnz
    LENIW=nnz+(4*NNODE)+32
    deallocate(RWORK,IWORK)
    allocate(RWORK(LENW),IWORK(LENIW))
    RWORK=0.1;IWORK=1
    TOL=1D-12; ERR=0.001

    
    !note: this loop is now 10000. Normally it is much higher 
    DO I = 1,20
        
        DO J=1,NNODE
            IF(poreBnd(J).EQ.1)THEN
                conc_node(J) = conc_inlet
            ENDIF
        ENDDO
                    
        call soluteCoefRHS(rhs, nnode,ntube,poreBnd,CONF,Q,POREV,Dmol, VOLUMEV,AREAV, DT, p_origin,conc_node,conc_edje,PASSTO,conc_inlet)

        !call DSDBCG (nnode,rhs,conc_node,nnz,ai,aj,bc,ISYM,ITOL, &
        !    TOL,ITMAX, ITER, ERR, IERR, IUNIT, RWORK,  LENW, IWORK, LENIW)
        !IF(IERR.EQ.0)THEN
        !    WRITE(*  ,'(A45)')," solverReturnFlag: IERR = 0 => All went well."
        !else
        !    WRITE(*  ,'(A45)')," solverReturnFlag: IERR /= 0 => Not Succeed!!"
        !ENDIF
        IF(I.EQ.1) THEN
            if(solverLibrary.eq.1)THEN
                solver_pointer = ses_solve_begin_density_cpu(nnode, nnz, ai, aj, bc, rhs, conc_node, iterationCount, precision, useOpenMp, numOfThreads, preconditioner)
            endif
    
            ! solve using gpu
            if(solverLibrary.eq.2 .or. solverLibrary.eq.3)THEN
                ! here first we need to take platforms and devices
                returnValue = ses_write_devices_to_file()
                ! choose your platform and device by index - indices start from 0
                solver_pointer = ses_solve_begin_density_gpu(nnode, nnz, ai, aj, bc, rhs, conc_node, solverLibrary, iterationCount, precision, platform , device, preconditioner)
        
            endif
        ELSE
            returnValue = ses_solve_next(solver_pointer, rhs, conc_node,iterationCount, precision)
        ENDIF
        WRITE(*, *)," Min conc_node= ", minval(conc_node), " Max conc_node= ", maxval(conc_node)
        
        DO J=1,NTUBE
            conc_edje(J)=(conc_edje(J)/COEFB(J))+ &
                ( COEF1(J)+COEF2(J) )*conc_node(UPNODE(J)) + &
                ( COEF2(J)*conc_node(DownNODE(J)) )
        ENDDO
    ENDDO


  ! Wait for user input before exiting
    write(*,*) "Press Enter to exit"
    read(*,*)
    end program Main
