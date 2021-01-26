  program schrodinger_1d
  implicit none
  !---------------------!
  ! error for bisection !
  !---------------------!
  real(8), parameter :: eps=1.d-6

  integer :: i, n
  real(8) :: L, psi0, psi1, E0, dE, k, d, V, E3
  real(8) :: f, mpsi0, mpsi1, npsi, dx, E ! f=2*[V(x)-E]; mpsi - modified psi, mpsi=psi*[1-dx**2*f/12]
  common/block/f, mpsi0, mpsi1, npsi, dx, E
  real(8), allocatable :: psi(:)  

  L=1.d0
  dE=1.d-1
  E0=5.d-1
  psi0=0.d0
  psi1=1.d-4
  n=10000
  dx=L/dble(2*n)
  allocate(psi(-n:n))
  psi(-n)=psi0
  psi(-n+1)=psi1
  f=2.d0*(V(-L)-E0)  
  mpsi0=psi0*(1.d0-dx**2*f/12.d0)
  f=2.d0*(V(-L+dx)-E0)  
  mpsi1=psi1*(1.d0-dx**2*f/12.d0)
  
  call integrate(psi, E0, n)
  do
    k=psi(n)
    E=E0+dE
    call integrate(psi, E, n)
    if (psi(n)*k < 0) then 
      print*, E
      exit
    endif
    E0=E
    print*,E
  enddo

  do
    E3=(E+E0)/2
    d=psi(n)
    call integrate(psi, E3, n)
    k=psi(n)
    if (abs(d)<=eps) exit
    if (k*d<=0) then 
      E=E3
    else
      E0=E3
    endif
  enddo


  end program



  subroutine integrate(psi, E0, n)
  implicit none

  integer :: i, n
  real(8) :: psi(-n:n), E0, x
  real(8) :: f, mpsi0, mpsi1, npsi, dx, E
  common/block/f, mpsi0, mpsi1, npsi, dx, E
  do i=-n+2,n
    x=dble(i)*dx                          !--x for calculating the potential--!
    call numerov(x)
    psi(i)=npsi  
  enddo
  call normalize(psi, n, dx)
  end subroutine integrate


  subroutine numerov(x)
  implicit none

  real(8) :: f, mpsi0, mpsi1, npsi, dx, E
  common/block/f, mpsi0, mpsi1, npsi, dx, E
  real(8) :: x, V, mpsi2
  
  f=2.d0*(V(x)-E)                       !-----------------------------------!
  mpsi2=2.d0*mpsi1-mpsi0+dx*dx*f*npsi   !------calculate modified psi-------!
  mpsi0=mpsi1; mpsi1=mpsi2              !--------new mpsi1 and mpsi0--------!
  npsi=mpsi1/(1.d0-dx*dx*f/12.d0)       !calculate normal psi for this step-!

  end subroutine numerov



  subroutine normalize(psi, n, dx)
  implicit none

  integer :: i, n
  real(8) :: psi(-n:n), norm, dx
  norm=psi(-n)**2+psi(n)**2
  do i=-n+1,n-3,2
    norm=norm+4.d0*psi(i)**2+2.d0*psi(i+1)**2
  enddo
  norm=norm+4.d0*psi(n-1)**2
  norm=1.d0/sqrt(norm*dx/3.d0)
  do i=-n,n
    psi(i)=psi(i)*norm
  enddo
  end subroutine normalize

  !----------------!
  ! Potential V(x) !
  !----------------!
  real(8) function V(x)
    implicit none 

    real(8) :: x
    V=0.d0
  end function V
