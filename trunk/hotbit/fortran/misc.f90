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
!real(8) :: wft(0:norb-1,0:norb-1) this is not needed anymore
rho=0d0
!wft=transpose(wf) 

do i=norb-1,0,-1
    if( occ(i)>1E-15 ) then
        mx=i
        exit
    end if
end do    

do i=0,norb-1
do j=i,norb-1
    rho(i,j)=sum( occ(0:mx)*wf(0:mx,i)*wf(0:mx,j) )
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
!complex(8) :: wft(0:norb-1,0:norb-1) not needed anymore


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
   !wft=transpose(wf(k,:,:))
   do i=0,norb-1
   do j=i,norb-1
       rho(k,i,j)=sum( occ(k,0:mx)*wf(k,0:mx,i)*conjg(wf(k,0:mx,j)) )
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
!real(8) :: wft(0:norb-1,0:norb-1)
integer :: i,j,mx
rhoe=0d0
!wft=transpose(wf)

do i=norb-1,0,-1
    if( occ(i)>1E-15 ) then
        mx=i
        exit
    end if
end do    

do i=0,norb-1
do j=i,norb-1
    rhoe(i,j)=sum( e(0:mx)*occ(0:mx)*wf(0:mx,i)*wf(0:mx,j) )
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
!complex(8) :: wft(0:norb-1,0:norb-1)

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
   !wft=transpose(wf(k,:,:))
   do i=0,norb-1
   do j=i,norb-1
       rho(k,i,j)=sum( e(k,0:mx)*occ(k,0:mx)*wf(k,0:mx,i)*conjg(wf(k,0:mx,j)) )
       rho(k,j,i)=conjg(rho(k,i,j))
   end do
   end do
end do
end subroutine fortran_rhoec



! check that any atom in cell (n1,n2,n3) interacts with atoms
! in cell (0,0,0)
! TODO invent a better name for this subroutine
subroutine fortran_doublefor(rn,r,cut,natoms,bool)
implicit none
integer, intent(in) :: natoms
real(8), intent(in) :: rn(0:2,0:natoms-1)
real(8), intent(in) :: r(0:2,0:natoms-1)
real(8), intent(in) :: cut(0:natoms-1,0:natoms-1)
integer, intent(out) :: bool
integer :: i,j
real(8) :: drx,dry,drz, dist, cu

bool = 0
do i=0, natoms-1
  if (bool == 1) exit
  do j=0, natoms-1
    drx = rn(0,i)-r(0,j)
    dry = rn(1,i)-r(1,j)
    drz = rn(2,i)-r(2,j)
    dist = drx**2+dry**2+drz**2
    cu = cut(j,i)
    if (dist .LE. cu) then
      bool=1
      exit
    endif
  end do
end do
end subroutine fortran_doublefor




subroutine geometry(rn,r,cut,natoms,bool)
implicit none
integer, intent(in) :: natoms
real(8), intent(in) :: rn(0:2,0:natoms-1)
real(8), intent(in) :: r(0:2,0:natoms-1)
real(8), intent(in) :: cut(0:natoms-1,0:natoms-1)
integer, intent(out) :: bool
integer :: i,j
real(8) :: drx,dry,drz, dist, cu

bool = 0
do i=0, natoms-1
  if (bool == 1) exit
  do j=0, natoms-1
    drx = rn(0,i)-r(0,j)
    dry = rn(1,i)-r(1,j)
    drz = rn(2,i)-r(2,j)
    dist = drx**2+dry**2+drz**2
    cu = cut(j,i)
    if (dist .LE. cu) then
      bool=1
      exit
    endif
  end do
end do
end subroutine geometry



! 
! ! subroutine matrix_blocks(ht,st,dht,dst,DT,R,phase,o1i,o1j,noi,noj,nk,norb,H0) !,S,dH0,dS)
! ! subroutine matrix_blocks(o1i,o1j,noi,noj,nk,norb,ht,st,dht,dst,DT,R,phase,H0) !,S,dH0,dS)
! subroutine matrix_blocks(norb,o1i,o1j,noi,noj,nk,ht,DT,R,phase,H0) !,S,dH0,dS)
! !subroutine matrix_blocks(norb,noi,noj,nk,ht,DT,R,phase,H0) !,S,dH0,dS)
! ! subroutine matrix_blocks(norb,nk,noi,noj,ht) !,S,dH0,dS)
! implicit none
! integer, intent(in) :: norb
! integer, intent(in) :: o1i
! integer, intent(in) :: o1j
! integer, intent(in) :: noi
! integer, intent(in) :: noj
! integer, intent(in) :: nk
! real(8), intent(in) :: ht(:,:)
! ! real, intent(in) :: st(0:noi-1,0:noj-1)
! ! real, intent(in) :: dht(0:noi-1,0:noj-1,0:2)
! ! real, intent(in) :: dst(0:noi-1,0:noj-1,0:2)
! real(8), intent(in) :: DT(0:8,0:8)
! real(8), intent(in) :: R(0:2,0:2)
! 
! complex(8), intent(in) :: phase(0:nk-1)
! complex(8), intent(out) :: H0(0:nk-1,0:norb-1,0:norb-1)
! ! complex(8), intent(out) :: S(0:nk-1,0:norb-1,0:norb-1)
! ! complex(8), intent(out) :: dH0(0:nk-1,0:norb-1,0:norb-1,0:2)
! ! complex(8), intent(out) :: dS(0:nk-1,0:norb-1,0:norb-1,0:2)
! complex(8) :: hblock(0:nk-1,0:noi-1,0:noj-1)
! ! complex(8) :: hblock(0:nk-1,0:noi-1,0:noj-1)
! integer :: i,j,k,p
! integer :: a,b,c,d
! 
! 
! do j=0,noj-1
! do i=0,noi-1
! do k=0,nk-1
!     hblock(k,i,j) = phase(k)*sum( ht(i,0:noj-1)*DT(0:noj-1,j) )
! end do
! end do
! end do
! a = o1i
! b = o1i+noi
! c = o1j
! d = o1j+noj
! H0(:,a:b,c:d) = H0(:,a:b,c:d) + hblock
! H0(:,c:d,a:b) = conjugate( transpose(hblock) )
! 
! end subroutine matrix_blocks
! 
