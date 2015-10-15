#include "fintrf.h"


subroutine mexFunction(nlhs, plhs, nrhs, prhs)

    use gdma_driver
    implicit none
    
    interface
        function ReadField(structPtr, fieldName, refM, refN)
            real*8, allocatable :: ReadField(:, :)
            mwPointer :: structPtr
            character*(*) :: fieldName
            integer refM, refN
        end function ReadField
        function Vector(tempMatrix)
            real*8, allocatable :: Vector(:)
            real*8, allocatable :: tempMatrix(:, :)
        end function Vector
        function Scalar(tempMatrix)
            real*8 Scalar
            real*8, allocatable :: tempMatrix(:, :)
        end function Scalar
    end interface
    
    mwPointer :: plhs(*), prhs(*)
    integer :: nlhs, nrhs
    
    mwPointer mxGetPr
    mwPointer mxCreateDoubleMatrix
    type(gdma_input) input_args
    real(kind(1d0)), dimension(:,:), allocatable :: multipoles
    
    integer numSites, numShells, numBasisFunc, numPrims
    mwPointer size_m, size_n
    
    if(nrhs .ne. 1) & 
        call mexErrMsgIdAndTxt("mexFunction:nrhs", "Expect only 1 input.")
        
    ! sites info
    input_args%limit = Vector(ReadField(prhs(1), 'limit', -1, -1))
    numSites = size(input_args%limit)
    input_args%nucleiCharges = Vector(ReadField(prhs(1), 'nucleiCharges', numSites, 1))
    input_args%xyzSites = ReadField(prhs(1), 'xyzSites', 3, numSites)
    
    ! shells info
    input_args%shellNfuncs = Vector(ReadField(prhs(1), 'shellNfuncs', -1, -1))
    numShells = size(input_args%shellNfuncs)
    input_args%shellNprims = Vector(ReadField(prhs(1), 'shellNprims', numShells, 1))
    input_args%shell2atom = Vector(ReadField(prhs(1), 'shell2atom', numShells, 1))
    
    ! primitives info
    numPrims = sum(input_args%shellNprims)
    input_args%primExps = Vector(ReadField(prhs(1), 'primExps', numPrims, 1))
    input_args%primCoefs = Vector(ReadField(prhs(1), 'primCoefs', numPrims, 1))
    
    numBasisFunc = sum(input_args%shellNfuncs)
    input_args%density = ReadField(prhs(1), 'density', numBasisFunc, numBasisFunc)
    
    input_args%bigexp = Scalar(ReadField(prhs(1), 'bigexp', 1, 1))
    
    call gdma_driver_routine(multipoles, input_args)
    
    size_m = size(multipoles, 1)
    size_n = size(multipoles, 2)
    plhs(1) = mxCreateDoubleMatrix(size_m, size_n, 0)
    call mxCopyReal8ToPtr(multipoles, mxGetPr(plhs(1)), size_m*size_n)
    deallocate(input_args%limit)
    deallocate(input_args%nucleiCharges)
    deallocate(input_args%xyzSites)
    deallocate(input_args%shellNfuncs)
    deallocate(input_args%shellNprims)
    deallocate(input_args%shell2atom)
    deallocate(input_args%primExps)
    deallocate(input_args%primCoefs)
    deallocate(input_args%density)
    deallocate(multipoles)
end subroutine mexFunction


function Vector(tempMatrix)
    real*8, allocatable :: Vector(:)
    real*8, allocatable :: tempMatrix(:, :)
    
    Vector = tempMatrix(:, 1)
end function Vector


function Scalar(tempMatrix)
    real*8 :: Scalar
    real*8, allocatable :: tempMatrix(:, :)
    
    Scalar = tempMatrix(1, 1)
end function Scalar


function ReadField(structPtr, fieldName, refM, refN) ! returns a matrix
    real*8, allocatable :: ReadField(:, :)
    mwPointer :: structPtr
    character*(*) :: fieldName
    integer :: refM, refN
    
    integer dimM, dimN
    integer*4 mxIsStruct
    integer*4 mxIsClass
    mwPointer mxGetPr
    mwPointer mxGetField
    mwPointer mxGetProperty
    mwPointer mxGetM
    mwPointer mxGetN
    mwPointer size_m, size_n
    mwPointer tempMwPointer
    
    if(mxIsStruct(structPtr) .eq. 1 ) then
        tempMwPointer = mxGetField(structPtr, 1, fieldName)
    else
        tempMwPointer = mxGetProperty(structPtr, 1, fieldName)
    end if
    
    ! field existence check
    if(tempMwPointer .eq. 0) &
        call mexErrMsgIdAndTxt("ReadField:structPtr", "Field "//fieldName//" does not exist.")
    
    size_m = mxGetM(tempMwPointer)
    size_n = mxGetN(tempMwPointer)
    
    if(min0(size_m, size_n) .eq. 1) then ! vector
        dimM = max0(size_m, size_n) ! make sure we output a column
        dimN = 1;
    else
        dimM = size_m
        dimN = size_n
    end if
    
    if(refM .ne. -1) then ! if -1 then no dimension check
        call DimensionCheck(fieldName, dimM, refM)
        call DimensionCheck(fieldName, dimN, refN)
    end if
    
    allocate(ReadField(dimM, dimN))
    
    call mxCopyPtrToReal8(mxGetPr(tempMwPointer), ReadField, size_m * size_n)
    
    if(isnan(sum(ReadField))) then
        call mexErrMsgIdAndTxt("ReadField:structPtr", "Field "//fieldName//" has NaN in it.")
    end if
end function ReadField


subroutine DimensionCheck(fieldName, dim_, dim_ref)
    character*(*) :: fieldName
    integer dim_, dim_ref
    if(dim_ .ne. dim_ref) then
        call mexErrMsgIdAndTxt("DimensionCheck:dim_", "Field "//fieldName//" has wrong dimension.")
    end if
    return
end subroutine DimensionCheck


