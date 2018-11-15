Program ExoTP
  Implicit NONE

  Real*8,dimension(:),allocatable::Ep,Tint
  Integer::Np,i,j,diviseur,Nt
  real*8::Tau,dt,T,R,S,Tf,sigma,dtau,RTdtau,coeff
  REAL*8, PARAMETER :: Pi = 3.14159267

  Tau=5.d0
  dt=1.d0/100
  Np=500000
  T=300
  R=280
  Tf = 20.d0
  S=0.d0
  Nt = int(Tf/dt)
  dtau=dt/Tau
  RTdTau=R*T*dtau
  diviseur=int(Nt/100)
  coeff=1.d0/(1.d0+2.d0*dtau)

  Allocate(Ep(1:Np),Tint(1:Nt))

  do i=1,Np
    Ep(i) = rand(0) * (100000.d0) + 100000.d0
    S=S+Ep(i)
  end do

  ! open(unit=10,file='data')
  open(unit=11,file='Hist')

  Tint(1)=S/(R*Np)

  ! write(10,*) Tint(0), '0'

  do j=1,int(Tf/dt)
    S=0.d0

    do i=2,Np
      sigma= rn_std_normal_dist()
      Ep(i)=coeff*(Ep(i)+RTdTau*(1.d0+sigma**2)+2.d0*sigma*sqrt(RTdTau*Ep(i)))
      S=S+Ep(i)
    end do

    Tint(j)=S/(R*Np)
    ! write(10,*) Tint(j), j
    if (mod(j,diviseur)==0) then
      write(11,*) j,Tint(j)
      ! print*,j
    endif
  end do

  Deallocate(Ep,Tint)

contains

real(kind=8) function rn_std_normal_dist()

real*8 :: s = 0.449871d0, t = -0.386595d0, a = 0.19600d0, b = 0.25472d0, &
          r1 = 0.27597d0, r2 = 0.27846d0, u, v, x, y, q

do

  call random_number(u)
  call random_number(v)
  v = 1.7156d0 * (v - 0.5d0)

  x = u - s
  y = ABS(v) - t
  q = x**2 + y*(a*y - b*x)

  if (q < r1) exit
  if (q > r2) cycle
  if (v**2 < -4.0*log(u)*u**2) exit

end do

rn_std_normal_dist = v/u

end function rn_std_normal_dist

end program

! sigma : https://www.science-emergence.com/Articles/G%C3%A9n%C3%A9rer-des-nombres-al%C3%A9atoires-depuis-une-loi-normale-en-fortran-90/
