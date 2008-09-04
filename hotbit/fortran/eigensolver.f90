subroutine geig(A,B,n,ev,vec)
implicit none
real(8), intent(in) :: A(n,n)
real(8), intent(in) :: B(n,n)
integer, intent(in) :: n
real(8), intent(out) :: ev(n)
real(8), intent(out) :: vec(n,n)
integer :: liwork, lwork, info, i
integer, allocatable, save :: iwork(:)
real(8), allocatable, save :: work(:), B2(:,:)

liwork = 3+5*n 
lwork  = 1+6*n+2*n**2 

if( allocated(iwork) .and. size(B2(:,1))/=n ) then
    ! deallocate if problem size changes
    deallocate(iwork,work,B2)
end if

if( .not. allocated(iwork) ) then
    allocate( iwork(liwork),work(lwork),B2(n,n) )
end if

vec=A
B2=B
call dsygvd(1,'V','L',n,vec,n,B2,n,ev,work,lwork,iwork,liwork,info)

if(info/=0) then
    write(*,*) "***********************************************************"
    write(*,*) "Parameter INFO from dsygvd:",info
    write(*,*) "matrix size N=",n
    if ( info<0 ) then
        write(*,*) "The info-th argument for dsygvd had an illegal value."
    else if( info>0 .and. info<=n ) then
        write(*,*) "The algorithm failed to compute an eigenvalue while working"
        write(*,*) "on the submatrix lying in rows and columns INFO/(N+1) through mod(INFO,N+1)"
    else if( info>n ) then
        write(*,*) "If INFO = N + i, for 1 <= i <= N, then the leading minor of order i of "
        write(*,*) "S is not positive definite. The factorization of S could not be completed and"
        writE(*,*) "no eigenvalues or eigenvectors were computed. (S=overlap matrix)"
    end if
    write(*,*) "***********************************************************"
    write(*,*) "H="
    do i=1,min(4,n)
        write(*,*) A(i,1:min(4,n))
    end do
    write(*,*) '... and so on.'
    write(*,*) "S="
    do i=1,min(4,n)
        write(*,*) B(i,1:min(4,n))
    end do
    write(*,*) '... and so on.'
    stop 'Error in LAPACK diagonalization routine.'
end if
end subroutine geig