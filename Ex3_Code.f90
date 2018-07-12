!--------------------------------------------------------------------------------------------------
module Constants
!--------------------------------------------------------------------------------------------------
! Module to keep constants in
!--------------------------------------------------------------------------------------------------
	implicit none
	integer,parameter :: N=5										! Number of bodies
	integer,parameter :: Ndim=3*2*N							! Dimension of solution vector y
	integer,parameter :: N1=Ndim/4,N2=Ndim/2		! Integers used in iterative loops
	real*8,parameter  :: pi=acos(-1.d0)					! The value of pi
	real*8,parameter  :: G=4.d0*pi*pi						! Gravitational constant [AU^3/yr^2/Msolar]
	real*8,dimension(N) :: M										! Store masses in an array
	real*8  :: Mtot															! Total mass of system [Msolar]
	integer,parameter :: Nt=100000							! Number of time steps to do
	real*8,parameter :: alpha=1d-3							! Dimensionless parameter for adjusted timestep
	
end module Constants

!--------------------------------------------------------------------------------------------------
integer function vari(i)
!--------------------------------------------------------------------------------------------------
! Function to return a varying integer to access desired part of array
! Doing this to tidy up code as this varying integer is regularly used
!--------------------------------------------------------------------------------------------------
	implicit none
	integer,intent(in) :: i		! Takes in iteration integer i
	
	vari = 3*(i-1)						! Outputted integer varies with inputted integer i

end function vari

!--------------------------------------------------------------------------------------------------
integer function N2vari(i)
!--------------------------------------------------------------------------------------------------
! Same as above function vari(i), but adds on N2=Ndim/2=3*N
! This is done to access the second half of an array (e.g. velocities in the y-array)
! Doing this to tidy up code as this varying integer is regularly used
!--------------------------------------------------------------------------------------------------
	use Constants, only: N2
	implicit none
	integer,intent(in) :: i		! Takes in iteration integer i
	
	N2vari = N2+3*(i-1)				! Outputted integer varies with inputted integer i

end function N2vari

!--------------------------------------------------------------------------------------------------
program NBody
!--------------------------------------------------------------------------------------------------
! The main program
!--------------------------------------------------------------------------------------------------
	use Constants, only: Ndim
	implicit none
	real*8,dimension(Ndim) :: y0,y	! Solution vector y (initial and final)
	real*8 :: f(Ndim)
	real*8,dimension(3,Ndim/3/2) :: r
	
	call Set_Masses
	call init_cond(y0)							! Setup initial values for y
	call Validity(y0,y,f,r)

end program NBody

!--------------------------------------------------------------------------------------------------
subroutine Set_Masses
!--------------------------------------------------------------------------------------------------
! Subroutine to define the masses of each body and store them in an array for later use
!--------------------------------------------------------------------------------------------------
	use Constants, only: M,Mtot
	implicit none
	! Define the mass of each body in the simulation
	real*8,parameter  :: Ms=1.d0			! Mass of sun [Msolar]
	real*8,parameter  :: Mj=9.543d-4	! Mass of Jupiter [Msolar]
	
	M = (/Ms,Mj,2.85716656d-4,2.41315168d-8,5d0/)	! Store each mass in an array for later use
	Mtot = sum(M)											! Total mass of the system [Msolar]
	
end subroutine Set_Masses

!--------------------------------------------------------------------------------------------------
subroutine init_cond(y0)
!--------------------------------------------------------------------------------------------------
! Subroutine to set the initial conditions
!--------------------------------------------------------------------------------------------------
	use Constants
	implicit none
	real*8,dimension(Ndim),intent(out) :: y0	! Solution vector y
	real*8,parameter  :: dj=5.20336301d0			! Jupiter's initial orbital distance [AU]
	real*8,parameter  :: vj=2.7531151d0				! Jupiter's initial velocity [AU/yr]
	integer,external :: vari,N2vari
	real*8  :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5				! Initial positions
	real*8  :: vx1,vy1,vz1,vx2,vy2,vz2,vx3,vy3,vz3,vx4,vy4,vz4,vx5,vy5,vz5				! Initial velocities
	real*8, dimension(3) :: rc,vc						! Co-moving frame vectors
	integer :: i
	
	! Set initial positions of each body
	x1  = 0.d0;	y1  = 0.d0;	z1  = 0.d0				! Sun	
	x2  = dj;		y2  = 0.d0;	z2  = 0.d0				! Jupiter
	x3  = 9.582d0; y3  = 0d0;  z3 = 0d0
	x4  = dj; y4  = 4.48602642d-2;  z4 = 0d0
	x5  = -8d0; y5  = 3d0;  z5 = 0d0
	
	! Set initial velocities of each body
	vx1 = 0.d0;	vy1 = 0.d0;				 vz1 = 0.d0				! Sun
	vx2 = 0.d0;	vy2 = vj;					 vz2 = 0.d0				! Jupiter
	vx3 = 0.d0; vy3 = 2.0419478d0; vz3 = 0.d0
	vx4 = 0.210945021d-1; vy4 = vj; vz4 = 0.d0
	vx5 = 1d0; vy5 = -5d-1; vz5 = 0.d0
	! Set the initial values for the solution vector y(t=0)
	y0 = (/&
					x1,y1,z1,&
					x2,y2,z2,&
					x3,y3,z3,&
					x4,y4,z4,&
					x5,y5,z5,&
					vx1,vy1,vz1,&
					vx2,vy2,vz2,&
					vx3,vy3,vz3, &
					vx4,vy4,vz4, &
					vx5,vy5,vz5 &
											/)
	
	rc = 0.d0;vc = 0.d0
	
	! Transform to co-moving frame
	! Generate the comoving frame coordinates and velocities
	do i=1,N
		rc = rc+(M(i)/Mtot)*y0(1+vari(i):3+vari(i))			! Positions
		vc = vc+(M(i)/Mtot)*y0(1+N2vari(i):3+N2vari(i))	! Velocities
	enddo
	
	! Correct the solution vector
	do i=1,N
		y0(1+vari(i):3+vari(i)) = y0(1+vari(i):3+vari(i)) - rc	! Position correction
		y0(1+N2vari(i):3+N2vari(i)) = y0(1+N2vari(i):3+N2vari(i)) - vc	! Velocity correction
	enddo
	
end subroutine init_cond

!----------------------------------------------------------------------------
subroutine RK4(y, t, tstep,calls,r,v)
!----------------------------------------------------------------------------
! ***  fourth-order Runge-Kutta integrator                                ***
!----------------------------------------------------------------------------
  use Constants, only : N,Ndim
  implicit none
  real*8,intent(inout) :: y(Ndim)      ! the solution vector at t, y(t)
                                       ! on return, it contains y(t+dt)
  real*8,intent(in) :: t,tstep         ! time t, and time step to be integrated
  integer,intent(inout) :: calls
  real*8,dimension(3,N),intent(inout) :: r,v
  real*8 :: h                          ! half time step
  real*8,dimension(Ndim) :: f1, f2, f3, f4, ytmp1, ytmp2, ytmp3

  h=0.5d0*tstep

  call RHS(y,f1,calls,r,v)             ! calulate f(y,t)
  ytmp1 = y + h*f1                     ! half step

  call RHS(ytmp1,f2,calls,r,v)       ! calculate f(ytmp1,t+h)
  ytmp2 = y + h*f2                     ! half step in another way

  call RHS(ytmp2,f3,calls,r,v)       ! calculate f(ytmp2,t+h)
  ytmp3 = y + tstep*f3                 ! full step

  call RHS(ytmp3,f4,calls,r,v)   ! calculate f(ytmp3,t+tstep)

  y = y + (f1 + 2.d0*f2 + 2.d0*f3 + f4)*tstep/6.d0    ! Runge-Kutta recipe

end subroutine RK4

!--------------------------------------------------------------------------------------------------
subroutine RHS(y,f,calls,r,v)
!--------------------------------------------------------------------------------------------------
! Subroutine to calculate the r.h.years^2 to s^2s. vector f
!--------------------------------------------------------------------------------------------------
	use Constants, only : N,Ndim,N1,N2,M
	implicit none
	real*8,intent(in)  :: y(Ndim)			! Solution vector y is input
	real*8,intent(out) :: f(Ndim)			! r.h.s. vector f is output
	integer,intent(out) :: calls			! no. of calls to this subroutine
	real*8,dimension(3,N),intent(inout) :: r,v
	integer,external   :: vari,N2vari
	real*8,dimension(N2) :: Fo				! Array of forces on each object
	integer :: i
	
	calls=calls+1
	
	call Forces(y,Fo,r,v)
	
	do i=1,N
		f(1+vari(i):3+vari(i)) = y(1+N2vari(i):3+N2vari(i))				! Velocities
		f(1+N2vari(i):3+N2vari(i)) = Fo(1+vari(i):3+vari(i))/M(i)	! Accelerations
	enddo

end subroutine RHS

!--------------------------------------------------------------------------------------------------
subroutine Forces(y,Fo,r,v)
!--------------------------------------------------------------------------------------------------
! Subroutine to calculate the gravitational forces present
! Force calculated in [Msolar*AU/yr^2]
!--------------------------------------------------------------------------------------------------
	use Constants
	implicit none
	integer,external :: vari,N2vari
	real*8,external :: mag_rirj
	real*8,intent(in),dimension(Ndim) :: y	! Solution vector y contains positions
	real*8,intent(out) :: Fo(N2)						! Array of forces on each object
	real*8,dimension(3,N),intent(out) :: r,v
	integer :: i,j

	do i=1,N
		r(:,i) = y(1+vari(i):3+vari(i))
		v(:,i)  = y(1+N2vari(i):3+N2vari(i))
	enddo

	! Calculate the resulting forces on the sun and on jupiter [Msolar*AU/yr^2]
	Fo = 0.d0
	do j=1,N
		do i=1,N
			if (i==j) cycle
			Fo(1+vari(j):3+vari(j)) = Fo(1+vari(j):3+vari(j))+G*M(j)*M(i)*(r(:,i)-r(:,j))/mag_rirj(r,i,j)**3
		enddo
	enddo

end subroutine Forces

!--------------------------------------------------------------------------------------------------
subroutine Orbit(y0,y,f,ttot,r,v,calls)
!--------------------------------------------------------------------------------------------------
! Subroutine to calculate the coordinates of the orbits and write them to a data file to plot
!--------------------------------------------------------------------------------------------------
	use Constants
	implicit none
	real*8,intent(inout),dimension(Ndim) :: y0,y,f
	real*8,dimension(3,N),intent(inout) :: r,v
	real*8,intent(in) :: ttot				! Total time to integrate over
	integer,intent(out) :: calls		! No. of calls to subroutine RHS
	real*8  :: dt										! Timestep [yr]
	integer :: i=0.d0
	real*8 :: t=0.d0
	
	y=y0	 						! Initialise position vector y
	call RHS(y0,f,calls,r,v)	! Initialise f
	
	open(1,file="Orbit.dat")
	write(1,*) y0(1:N2)							! Write initial coordinates to file
	do while (t<=ttot)
		i=i+1
		t=t+dt
		call Timestep(r,v,dt)					! Adjust the timestep
		call RK4(y, t, dt,calls,r,v)
		write(1,*) y(1:N2)						! Write coordinates to file
		if (mod(i,10000)==0) then
			print '(F6.2,"%")', t/ttot *100.d0	
			!write(1,*) y(1:N2)						! Write coordinates to file
		endif
	enddo
	close(1)

	print '(A)', 'File has been written' ! To know subroutine ran successfully

end subroutine Orbit

!--------------------------------------------------------------------------------------------------
subroutine Timestep(r,v,dt)
!--------------------------------------------------------------------------------------------------
! Subroutine to adjust the timestep during the simulation
!--------------------------------------------------------------------------------------------------
	use Constants, only: alpha,N
	implicit none
	real*8,intent(in),dimension(3,N) :: r,v	! Position and velocity vectors of each body
	real*8,intent(out) :: dt								! Adjusted timestep 
	real*8,external :: mag_rirj,mag_vivj
	real*8,dimension(N,N) :: magr,magv,ratio	! Arrays to be used in this subroutine
	integer :: i,j 
	 
	magr=0.d0;magv=0.d0;ratio=0.d0
	do i=1,N
		do j=1,N
			if (j==i) cycle
			magr(j,i) = mag_rirj(r,i,j)
			magv(j,i) = mag_vivj(v,i,j)
			ratio(j,i) = magr(j,i)/magv(j,i)
		enddo
	enddo
	dt = alpha*minval(ratio,mask=ratio>0.d0)

end subroutine Timestep

!--------------------------------------------------------------------------------------------------
real*8 function mag_rirj(r,i,j)
!--------------------------------------------------------------------------------------------------
! Function to calculate the magnitude of ri-rj, where ri and rj are position vectors for bodies i and j
!--------------------------------------------------------------------------------------------------
	use Constants, only: N
	implicit none
	real*8,dimension(3,N),intent(in) :: r	
	integer,intent(in) :: i,j		! To select which body (e.g. 1=Sun,2=Jupiter)
	real*8,dimension(3) :: rirj
	
	rirj = r(:,i)-r(:,j)
	mag_rirj = sqrt(rirj(1)**2+rirj(2)**2+rirj(3)**2)
	
end function mag_rirj

!--------------------------------------------------------------------------------------------------
real*8 function mag_vivj(v,i,j)
!--------------------------------------------------------------------------------------------------
! Function to calculate the magnitude of vi-vj, where vi and vj are velocity vectors for bodies i and j
!--------------------------------------------------------------------------------------------------
	use Constants, only: N
	implicit none
	real*8,dimension(3,N),intent(in) :: v	
	integer,intent(in) :: i,j		! To select which body (e.g. 1=Sun,2=Jupiter)
	real*8,dimension(3) :: vivj
	
	vivj = v(:,i)-v(:,j)
	mag_vivj = sqrt(vivj(1)**2+vivj(2)**2+vivj(3)**2)
	
end function mag_vivj

!--------------------------------------------------------------------------------------------------
real*8 function magv2(y,i)
!--------------------------------------------------------------------------------------------------
! Function to calculate the magnitude squared of velocity vector
!--------------------------------------------------------------------------------------------------
	use Constants, only: Ndim
	implicit none
	real*8,dimension(Ndim),intent(in) :: y		! Make use of velocities stored in y
	integer,intent(in) :: i		! To select which body (e.g. 1=Sun,2=Jupiter)
	integer,external :: vari,N2vari
	
	magv2 = y(1+N2vari(i))**2 + y(2+N2vari(i))**2 + y(3+N2vari(i))**2
	
end function magv2

!--------------------------------------------------------------------------------------------------
subroutine Validity(y0,y,f,r)
!--------------------------------------------------------------------------------------------------
! Subroutine to check the validity of the simulation
! Will check the deviations of Etot,Lz, and position after one period compared to initial position
!--------------------------------------------------------------------------------------------------
	use Constants
	implicit none
	real*8,intent(inout),dimension(Ndim) :: y0,y,f
	real*8,dimension(3,N),intent(inout) :: r
	integer :: calls=0
	real*8,dimension(3,N) :: r0
	real*8,external :: mag_rirj,magv2
	integer,external ::vari,N2vari
	real*8  :: ttot
	real*8  :: Ekin0,Ekin,Epot0,Epot,Etot0,Etot
	real*8  :: a,P,pos0,pos
	integer :: i,j
	real*8,dimension(N2) :: Fo					! Array of forces on each object
	real*8,dimension(3,N) :: v0,v
	real*8 :: Lx0,Lx,Ly0,Ly,Lz0,Lz
	real*8,dimension(3) :: L0,L
	real*8 :: magL0,magL

	do i=1,N
		r0(:,i) = y0(1+vari(i):3+vari(i))
		v0(:,i) = y0(1+N2vari(i):3+N2vari(i))
	enddo
	
	Ekin0=0.d0;Epot0=0.d0
	do i=1,N
		Ekin0 = Ekin0+0.5d0*M(i)*magv2(y0,i)
		if (i==N) exit
		do j=(i+1),N
			Epot0 = Epot0 - G*M(i)*M(j)/mag_rirj(r0,i,j)
		enddo
	enddo
	
	Etot0 = Ekin0+Epot0
	
	! Calculate the semi-major axis of the orbit of Jupiter
	a = G*(Mtot-M(2))*M(2)/(2.d0*abs(Etot0))
	! Calculate the orbital period of Jupiter
	P = 2.d0*pi*sqrt(a**3/(G*(Mtot)))
	
	! Set the total length of time to integrate over
	ttot = 10d0*12d0!*P
	
	call Orbit(y0,y,f,ttot,r,v,calls)
	
	! Kinetic and potential energies
	Ekin=0.d0;Epot=0.d0
	do i=1,N
		Ekin = Ekin+0.5d0*M(i)*magv2(y,i)
		if (i==N) exit
		do j=(i+1),N
			Epot = Epot - G*M(i)*M(j)/mag_rirj(r,i,j)
		enddo
	enddo
	
	! Final total energy
	Etot = Ekin+Epot
		
	! Velocity vectors of each body
	do i=1,N
		v(:,i)  = y(1+N2vari(i):3+N2vari(i))
	enddo
	
	! Angular momentum of system
	L0=0.d0;L=0.d0
	do i=1,N
		! Initial
		Lx0 = r0(2,i)*v0(3,i)-r0(3,i)*v0(2,i)
		Ly0 = r0(3,i)*v0(1,i)-r0(1,i)*v0(3,i)
		Lz0 = r0(1,i)*v0(2,i)-r0(2,i)*v0(1,i)
		! Final
		Lx = r(2,i)*v(3,i)-r(3,i)*v(2,i)
		Ly = r(3,i)*v(1,i)-r(1,i)*v(3,i)
		Lz = r(1,i)*v(2,i)-r(2,i)*v(1,i)
		
		L0 = L0 + M(i)*(/Lx0,Ly0,Lz0/)
		L = L + M(i)*(/Lx,Ly,Lz/)
	enddo
	
	! Initial and final angular momentum of system
	magL0 = sqrt(L0(1)**2+L0(2)**2+L0(3)**2)
	magL = sqrt(L(1)**2+L(2)**2+L(3)**2)
	
	! Initial and final position of Jupiter
	pos0 = sqrt(r0(1,2)**2+r0(2,2)**2+r0(3,2)**2)
	pos = sqrt(r(1,2)**2+r(2,2)**2+r(3,2)**2)
	
	! Output the results
	print '(2(A,ES18.10,10x)/x)', &
											'a     = ',a, &
											'P     = ',P, &
											'Eki   = ',Ekin0, &
											'Ekf   = ',Ekin, &
											'Epi   = ',Epot0, &
											'Epf   = ',Epot, &
											'Etoti = ',Etot0, &
											'Etotf = ',Etot, &
											'magLi = ',magL0, &
											'magLf = ',magL															
	print '(A,3ES18.10/x)', &			
											'posi  = ',r0(:,2), &
											'posf  = ',r(:,2)
	print '(A,ES18.10/x)', &
											'Q1    = ',Etot/Etot0 -1.d0, &	! Energy Conservation
											'Q2    = ',magL/magL0 -1.d0, &	! L Conservation
											'Q3    = ',pos/pos0 - 1.d0			! Deviation from initial position
	print '(A,I10,A)',	'The RHS subroutine was called ',calls, ' times.'
	
end subroutine Validity
