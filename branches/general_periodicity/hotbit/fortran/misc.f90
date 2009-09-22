! return density matrix, where elements only with non-zero S are calculated
subroutine fortran_rho0(wf,occ,nr_ia_orbitals,ia_orbitals,norb,rho)
implicit none
integer, parameter :: dp = SELECTED_REAL_KIND(12)
integer, intent(in) :: norb
real(dp), intent(in) :: wf(0:norb-1,0:norb-1)
real(dp), intent(in) :: occ(0:norb-1)
integer, intent(in) :: nr_ia_orbitals(0:norb-1)
integer, intent(in) :: ia_orbitals(0:norb-1,0:norb-1)
real(dp), intent(out) :: rho(0:norb-1,0:norb-1)
integer :: i,j,k,high,fully
real(dp) :: wft(0:norb-1,0:norb-1)
real(dp), external ::ddot
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
integer, parameter :: dp = SELECTED_REAL_KIND(12)
integer, intent(in) :: norb
real(dp), intent(in) :: wf(0:norb-1,0:norb-1)
real(dp), intent(in) :: occ(0:norb-1)
real(dp), intent(in) :: e(0:norb-1)
integer, intent(in) :: nr_ia_orbitals(0:norb-1)
integer, intent(in) :: ia_orbitals(0:norb-1,0:norb-1)
real(dp), intent(out) :: rhoe(0:norb-1,0:norb-1)
real(dp) :: wft(0:norb-1,0:norb-1)
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
integer, parameter :: dp = SELECTED_REAL_KIND(12)
integer, intent(in) :: norb
real(dp), intent(in) :: wf(0:norb-1,0:norb-1)
real(dp), intent(in) :: occ(0:norb-1)
real(dp), intent(out) :: rho(0:norb-1,0:norb-1)
integer :: i,j,mx
real(dp) :: wft(0:norb-1,0:norb-1)
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
integer, parameter :: dp = SELECTED_REAL_KIND(12)
integer,    intent(in) :: norb
integer,    intent(in) :: nk
complex(2*dp), intent(in)  :: wf(0:nk-1,0:norb-1,0:norb-1)
real(dp),    intent(in)  :: occ(0:nk-1,0:norb-1)
complex(2*dp), intent(out) :: rho(0:nk-1,0:norb-1,0:norb-1)
integer :: i,j,k,mx
complex(2*dp) :: wft(0:norb-1,0:norb-1)


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
integer, parameter :: dp = SELECTED_REAL_KIND(12)
integer, intent(in) :: norb
real(dp), intent(in) :: wf(0:norb-1,0:norb-1)
real(dp), intent(in) :: occ(0:norb-1),e(0:norb-1)
real(dp), intent(out) :: rhoe(0:norb-1,0:norb-1)
real(dp) :: wft(0:norb-1,0:norb-1)
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
integer, parameter :: dp = SELECTED_REAL_KIND(12)
integer,    intent(in) :: norb,nk
complex(2*dp), intent(in)  :: wf(0:nk-1,0:norb-1,0:norb-1)
real(dp),    intent(in)  :: occ(0:nk-1,0:norb-1)
real(dp),    intent(in)  :: e(0:nk-1,0:norb-1)
complex(2*dp), intent(out) :: rho(0:nk-1,0:norb-1,0:norb-1)
integer :: i,j,k,mx
complex(2*dp) :: wft(0:norb-1,0:norb-1)

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





