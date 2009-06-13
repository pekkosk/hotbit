! return density matrix, where elements only with non-zero S are calculated
subroutine fortran_rho0(wf,occ,nr_ia_orbitals,ia_orbitals,norb,rho)
implicit none
integer, intent(in) :: norb
real(8), intent(in) :: wf(0:norb-1,0:norb-1)
real(8), intent(in) :: occ(0:norb-1)
integer, intent(in) :: nr_ia_orbitals(0:norb-1)
integer, intent(in) :: ia_orbitals(0:norb-1,0:norb-1)
real(8), intent(out) :: rho(0:norb-1,0:norb-1)
integer :: i,j,k,high,fully
real(8) :: wft(0:norb-1,0:norb-1)
real(8), external ::ddot
rho=0d0
wft=transpose(wf)

fully=-1
do i=0, norb-1
    if( occ(i)>1E-20 ) high=i
    if( abs(occ(i)-2d0)<1E-20 ) fully=i
end do
 

do i=0,norb-1
    do k=0,nr_ia_orbitals(i)-1
        j=ia_orbitals(i,k)  
        ! first fully occupied part...
        if( fully>=0 ) rho(i,j)=2*ddot(fully+1,wft(0:fully,i),1,wft(0:fully,j),1)
        ! ...then partly occupied.
        rho(i,j)=rho(i,j)+sum( occ(fully+1:high)*wft(fully+1:high,i)*wft(fully+1:high,j) )                
        rho(j,i)=rho(i,j)
    end do
end do
end subroutine fortran_rho0


! return density matrix weighted by energies only with non-zero S are calculated
subroutine fortran_rhoe0(wf,occ,e,nr_ia_orbitals,ia_orbitals,norb,rhoe)
implicit none
integer, intent(in) :: norb
real(8), intent(in) :: wf(0:norb-1,0:norb-1)
real(8), intent(in) :: occ(0:norb-1)
real(8), intent(in) :: e(0:norb-1)
integer, intent(in) :: nr_ia_orbitals(0:norb-1)
integer, intent(in) :: ia_orbitals(0:norb-1,0:norb-1)
real(8), intent(out) :: rhoe(0:norb-1,0:norb-1)
real(8) :: wft(0:norb-1,0:norb-1)
integer :: i,j,k,mx
rhoe=0d0
wft=transpose(wf)

do i=norb-1,0,-1
    if( occ(i)>1E-15 ) then
        mx=i
        exit
    end if
end do    

do i=0,norb-1
    do k=0,nr_ia_orbitals(i)-1
        j=ia_orbitals(i,k)  
        rhoe(i,j)=sum( e(0:mx)*occ(0:mx)*wft(0:mx,i)*wft(0:mx,j) )
        rhoe(j,i)=rhoe(i,j)
    end do
end do
end subroutine fortran_rhoe0




! return density matrix
subroutine fortran_rho(wf,occ,norb,rho)
implicit none
integer, intent(in) :: norb
real(8), intent(in) :: wf(0:norb-1,0:norb-1)
real(8), intent(in) :: occ(0:norb-1)
real(8), intent(out) :: rho(0:norb-1,0:norb-1)
integer :: i,j,mx
real(8) :: wft(0:norb-1,0:norb-1)
rho=0d0
wft=transpose(wf)

do i=norb-1,0,-1
    if( occ(i)>1E-15 ) then
        mx=i
        exit
    end if
end do    

do i=0,norb-1
do j=i,norb-1
    rho(i,j)=sum( occ(0:mx)*wft(0:mx,i)*wft(0:mx,j) )
    rho(j,i)=rho(i,j)
end do
end do
end subroutine fortran_rho



! return complex density matrix
subroutine fortran_rhoc(wf,occ,norb,nk,rho)
implicit none
integer,    intent(in) :: norb
integer,    intent(in) :: nk
complex(8), intent(in)  :: wf(0:nk-1,0:norb-1,0:norb-1)
real(8),    intent(in)  :: occ(0:nk-1,0:norb-1)
complex(8), intent(out) :: rho(0:nk-1,0:norb-1,0:norb-1)
integer :: i,j,k,mx
complex(8) :: wft(0:norb-1,0:norb-1)


first: do i=norb-1,0,-1
    do k=0,nk-1
    if( occ(k,i)>1E-15 ) then
        mx=i
        exit first
    end if
    end do
end do first

rho=0d0
do k=0,nk-1
   ! because the first index must be faster:
   wft=transpose(wf(k,:,:))
   do i=0,norb-1
   do j=i,norb-1
       rho(k,i,j)=sum( occ(k,0:mx)*wft(0:mx,i)*conjg(wft(0:mx,j)) )
       if(i/=j) then
          rho(k,j,i)=conjg(rho(k,i,j))
       end if
   end do
   end do
end do
end subroutine fortran_rhoc



! return density matrix weighted by energies
subroutine fortran_rhoe(wf,occ,e,norb,rhoe)
implicit none
integer, intent(in) :: norb
real(8), intent(in) :: wf(0:norb-1,0:norb-1)
real(8), intent(in) :: occ(0:norb-1),e(0:norb-1)
real(8), intent(out) :: rhoe(0:norb-1,0:norb-1)
real(8) :: wft(0:norb-1,0:norb-1)
integer :: i,j,mx
rhoe=0d0
wft=transpose(wf)

do i=norb-1,0,-1
    if( occ(i)>1E-15 ) then
        mx=i
        exit
    end if
end do    

do i=0,norb-1
do j=i,norb-1
    rhoe(i,j)=sum( e(0:mx)*occ(0:mx)*wft(0:mx,i)*wft(0:mx,j) )
    rhoe(j,i)=rhoe(i,j)
end do
end do
end subroutine fortran_rhoe


! return complex energy-weighted density matrix
subroutine fortran_rhoec(wf,occ,e,norb,nk,rho)
implicit none
integer,    intent(in) :: norb,nk
complex(8), intent(in)  :: wf(0:nk-1,0:norb-1,0:norb-1)
real(8),    intent(in)  :: occ(0:nk-1,0:norb-1)
real(8),    intent(in)  :: e(0:nk-1,0:norb-1)
complex(8), intent(out) :: rho(0:nk-1,0:norb-1,0:norb-1)
integer :: i,j,k,mx
complex(8) :: wft(0:norb-1,0:norb-1)

rho=0d0
first: do i=norb-1,0,-1
    do k=0,nk-1
    if( occ(k,i)>1E-15 ) then
        mx=i
        exit first
    end if
    end do
end do first

do k=0,nk-1
   ! because the state index must be faster, it must come first
   wft=transpose(wf(k,:,:))
   do i=0,norb-1
   do j=i,norb-1
       rho(k,i,j)=sum( e(k,0:mx)*occ(k,0:mx)*wft(0:mx,i)*conjg(wft(0:mx,j)) )
       rho(k,j,i)=conjg(rho(k,i,j))
   end do
   end do
end do
end subroutine fortran_rhoec


! matrix multiplication a*b, where a is a symmetric matrix
subroutine symmetric_matmul(a,b,c,n)
implicit none
integer, intent(in) :: n
real(8), intent(in) :: a(0:n-1,0:n-1),b(0:n-1,0:n-1)
real(8), intent(out) :: c(0:n-1,0:n-1)
call dsymm('L','U',n,n,1d0,a,n,b,n,0,c,n)
end subroutine symmetric_matmul


! return diagonal of matrix multiplication a*b
subroutine matmul_diagonal(a,b,c,n)
implicit none
integer, intent(in) :: n
real(8), intent(in) :: a(0:n-1,0:n-1),b(0:n-1,0:n-1)
real(8), intent(out) :: c(0:n-1)
integer :: i
do i=0,n-1
    c(i)=sum( a(i,:)*b(:,i) )
end do
end subroutine matmul_diagonal



! Return the band-structure energy
subroutine fortran_fbs(rho,rhoe,dH,dS,norbs,indices,f,norb,nat)
implicit none
integer, intent(in) :: norb
integer, intent(in) :: nat
real(8), intent(in) :: rho(0:norb-1,0:norb-1)
real(8), intent(in) :: rhoe(0:norb-1,0:norb-1)
real(8), intent(in) :: dH(0:norb-1,0:norb-1,0:2)
real(8), intent(in) :: dS(0:norb-1,0:norb-1,0:2)
integer, intent(in) :: norbs(0:nat-1)
integer, intent(in) :: indices(0:nat-1,0:8)
real(8), intent(out) :: f(0:nat-1,0:2)
real(8) :: diff(0:norb-1)
integer :: a,i,noi

f=0d0
do a=0,2
    do i=0,norb-1
        ! diff is the diagonal of difference rho*dH - rhoe*dS
        diff(i)=sum(rho(i,:)*dH(:,i,a))-sum(rhoe(i,:)*dS(:,i,a))
    end do

    do i=0,nat-1
        noi=norbs(i)
        f(i,a)=sum( diff(indices(i,0:noi-1)) )
    end do
end do
end subroutine fortran_fbs    



! Return the band-structure energy (complex version)
! subroutine fortran_fbsc(rho,rhoe,dH,dS,norbs,indices,wk,f,norb,nat,nk)
! implicit none
! integer, intent(in) :: norb
! integer, intent(in) :: nat
! integer, intent(in) :: nk
! complex(8), intent(in) :: rho(0:nk-1,0:norb-1,0:norb-1)
! complex(8), intent(in) :: rhoe(0:nk-1,0:norb-1,0:norb-1)
! complex(8), intent(in) :: dH(0:nk-1,0:norb-1,0:norb-1,0:2)
! complex(8), intent(in) :: dS(0:nk-1,0:norb-1,0:norb-1,0:2)
! integer, intent(in) :: norbs(0:nat-1)
! integer, intent(in) :: indices(0:nat-1,0:8)
! real(8), intent(in) :: wk(0:nk-1)
! real(8), intent(out) :: f(0:nat-1,0:2)
! complex(8) :: diff(0:norb-1)
! integer :: a,i,noi,ik
! complex(8) :: f2(0:nat-1,0:2)
! 
! f2=0d0
! do ik=0,nk-1
!   do a=0,2
!       do i=0,norb-1
!         ! diff is the diagonal of difference rho*dH - rhoe*dS
!         diff(i) = sum(rho(ik,i,:)*dH(ik,:,i,a))-sum(rhoe(ik,i,:)*dS(ik,:,i,a))
!       end do
! 
!       do i=0,nat-1
!         noi = norbs(i)
!         f2(i,a) = f2(i,a) + wk(ik)*sum( diff(indices(i,0:noi-1)) )
!       end do
!   end do
! end do
! f = 2*real(f2)
! end subroutine fortran_fbsc      




subroutine fortran_fbsc(rho,rhoe,dH,dS,norbs,indices,wk,f,norb,nat,nk)
implicit none
integer, intent(in) :: norb
integer, intent(in) :: nat
integer, intent(in) :: nk
complex(8), intent(in) :: rho(0:nk-1,0:norb-1,0:norb-1)
complex(8), intent(in) :: rhoe(0:nk-1,0:norb-1,0:norb-1)
complex(8), intent(in) :: dH(0:nk-1,0:norb-1,0:norb-1,0:1,0:2)
complex(8), intent(in) :: dS(0:nk-1,0:norb-1,0:norb-1,0:1,0:2)
integer, intent(in) :: norbs(0:nat-1)
integer, intent(in) :: indices(0:nat-1,0:8)
real(8), intent(in) :: wk(0:nk-1)
real(8), intent(out) :: f(0:nat-1,0:2)
complex(8) :: diag(0:norb-1)
integer :: a,i,noi,ik
complex(8) :: f2(0:nat-1,0:2)

f2=0d0
do ik=0,nk-1
do a=0,2
      do i=0,norb-1
        ! diff is the diagonal of rho*dH2 + dH1*rho-rhoe*dS2-dS1*rhoe
        diag(i) = sum(  rho(ik,i,:)*dH(ik,:,i,1,a) + dH(ik,i,:,0,a)*rho(ik,:,i) &
                      -rhoe(ik,i,:)*dS(ik,:,i,1,a) - dS(ik,i,:,0,a)*rhoe(ik,:,i) )
      end do

      do i=0,nat-1
        noi = norbs(i)
        f2(i,a) = f2(i,a) - wk(ik)*sum( diag(indices(i,0:noi-1)) )
      end do
end do
end do
if( any(abs(imag(f2))>1E-10) ) then
  stop "complex forces"
end if
f = real(f2)
end subroutine fortran_fbsc           
