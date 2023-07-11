
    Module solveMatrix
    use, intrinsic :: iso_fortran_env, only : wp => real64
    use iso_c_binding
    implicit none

    
    !!!note: in the project folder i have put two matrices for examples. a larger size (with name append of _L) and a medium size (with name append of _M). 
    
    !medium size
    !character(len=30) :: matrixFile = "Matrix.txt", bFile = "Pknown.txt"
    !integer, parameter :: n = 22591, nnz = 67068     ! medium size (with name append of _M)
    
    !larger size
    integer, parameter :: n = 997350, nnz = 26218238     !larger size (with name append of _L)
    character(len=30) :: matrixFile = "Matrix_L", bFile = "Pknown_L"
    
    
    integer, allocatable :: rowIndex(:), colIndex(:)
    real(c_double), allocatable :: Avalues(:), b(:), x(:)
    integer :: i, j, error
    integer, allocatable :: ia(:), ja(:)

    contains

    subroutine readSparseMatrix(filename, nnz, rowIndex, colIndex, Avalues)
    !this subroutine reads 'Matrix' from the file and saves it in rowIndex, colIndex, and Avalues vectors
    character(len=*), intent(in) :: filename
    integer, intent(in) :: nnz
    integer :: rowIndex(nnz), colIndex(nnz)
    real(wp) :: Avalues(nnz)
    integer :: i, unit
    unit  = 1001
    open(unit, file=filename, action="read", status="old")
    do i = 1, nnz
        read(unit, *) rowIndex(i), colIndex(i), Avalues(i)
    end do
    close(unit)
    end subroutine

    subroutine readVector(filename, n, b)
    !this subroutine reads 'Pknown' from the file and saves it in 'b' vector
    character(len=*), intent(in) :: filename
    integer, intent(in) :: n
    real*8, intent(out) :: b(n)
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
    implicit none

    integer :: idx, errCode
    integer :: num_threads
    integer(c_int) :: status
    type(c_ptr) :: dll_handle
    interface
        function solve_matrix(numRows, numNonzero, rowIndices, colIndices, values, b, x) bind(C, name="solve_matrix")
            use iso_c_binding
            integer(c_int), value :: numRows, numNonzero
            integer(c_int), dimension(*) :: rowIndices, colIndices
            real(c_double), dimension(*) :: values, b, x
            integer(c_int) :: solve_matrix
        end function solve_matrix
    end interface
    ! Allocate rowIndex, colIndex, Avalues, b, x, ia, and ja arrays
    allocate(rowIndex(nnz), colIndex(nnz), Avalues(nnz), b(n), x(n), ia(n+1), ja(nnz))

    ! Read the sparse matrix A and vector b from the text files
    call readSparseMatrix(matrixFile, nnz, rowIndex, colIndex, Avalues)
    call readVector(bFile, n, b)

    x = solve_matrix(n, nnz, rowIndex, colIndex, Avalues, b, x)
    ! Wait for user input before exiting
    write(*,*) "Press Enter to exit"
    read(*,*)

    end program Main
