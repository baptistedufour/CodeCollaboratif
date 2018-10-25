Program ExoTP
  Implicit NONE

  Real*8,dimension(:),allocatable::Ep,Tint
  Integer::Np,i,j,diviseur,k
  real*8::Tau,dt,T,R,S,Tf,temps,sigma,x
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
    Ep(i) = rand(0) * (200000-100000) + 100000
    S=S+Ep(i)
  end do

  open(unit=10,file='data')
  open(unit=11,file='Hist')

  Tint(0)=S/(R*Np)

  write(10,*) Tint(0), '0'

  j=1
  temps=0.d0

  do while (temps<Tf)
    S=0.d0
    do i=1,Np
      x=rand(0) * (200000+200000) - 200000
      sigma=exp((-x**2)/2)*1/sqrt(2*Pi)
      Ep(i)=1.d0/(1+2.d0*dt/Tau)*(Ep(i)+R*T*(dt/tau)*(1+sigma**2)+2*sigma*sqrt((dt/Tau)*R*T*Ep(i)))
      S=S+Ep(i)
    end do
    Tint(j)=S/(R*Np)
    write(10,*) Tint(j), j
    diviseur=int(Tf/(dt*100))
    if (mod(j,diviseur)==0) then
      write(11,*) j,Tint(j)
    endif
    j=j+1
    temps=temps+dt
  end do

  Deallocate(Ep,Tint)



end program

! sigma : https://www.science-emergence.com/Articles/G%C3%A9n%C3%A9rer-des-nombres-al%C3%A9atoires-depuis-une-loi-normale-en-fortran-90/
