Program ExoTP
  Implicit NONE

  Real*8,dimension(:),allocatable::Ep,Tint
  Integer::Np,i,j,diviseur,k
  real*8::Tau,dt,T,R,S,Tf,sigma,loinorm
  REAL*8, PARAMETER :: Pi = 3.14159267
  character(len=1)::str

  Tau=5
  dt=1.d0/100
  Np=500000
  T=300
  R=280

  Tf = 4.d0*Tau

  Allocate(Ep(1:Np),Tint(0:int(Tf/dt)))

  do i=1,Np
    Ep(i) = rand(0) * (100000.d0) + 100000.d0
    S=S+Ep(i)
  end do

  open(unit=10,file='data')
  open(unit=11,file='Hist')

  Tint(0)=S/(R*Np)

  write(10,*) Tint(0), '0'
  do j=1,int(Tf/dt)
    S=0.d0

    do i=1,Np
      sigma= rn_std_normal_dist()
      Ep(i)=1.d0/(1.d0+2.d0*dt/Tau)*(Ep(i)+R*T*(dt/tau)*(1.D0+sigma**2)+2.d0*sigma*sqrt((dt/Tau)*R*T*Ep(i)))
      S=S+Ep(i)
    end do

    Tint(j)=S/(R*Np)
    write(10,*) Tint(j), j
    diviseur=int(Tf/(dt*100))
    if (mod(j,diviseur)==0) then
      write(11,*) j,Tint(j)
      print*,j
    endif
  end do

  Deallocate(Ep,Tint)

contains

real(kind=8) function rn_std_normal_dist()

real :: half = 0.5

real :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472, &
        r1 = 0.27597, r2 = 0.27846, u, v, x, y, q
      x=rand(0) * (200000.d0+200000.d0) - 200000.d0

do

  call random_number(u)
  call random_number(v)
  v = 1.7156 * (v - half)

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
