!--------------------------------------------------------------------------------------------------
module Constants
!--------------------------------------------------------------------------------------------------
! Module to keep constants in
!--------------------------------------------------------------------------------------------------
	implicit none
	integer,parameter :: N=2					! Number of bodies
	integer,parameter :: sel=1d1					! For selecting when to write to file
	character(len=99),parameter :: xyz="Binary_Stars.dat"	! Name of file to save coordinates to
	character(len=99),parameter :: time="T_Binary_Stars.dat"		! Name of file to save times to
	real*8,parameter :: alpha=1d-3			! Dimensionless parameter for adjusted timestep
	real*8,parameter :: L=1d2					! Arbitrary length [Mpc] of box for periodic B.C.
	integer,parameter :: Ndim=3*2*N			! Dimension of solution vector y
	integer,parameter :: N3=3*N				! Integer used in iterative loops
	real*8 :: M(N)									! Store masses in an array
	real*8 :: Mtot									! Total mass of system [Msolar]
	real*8,parameter :: pi=acos(-1.d0)		! The value of pi
	real*8,parameter :: pc=4.84814d-6		! pc per AU
	real*8,parameter :: Mega=1d6				! Mega prefix (10^6)
	real*8,parameter :: G=4.d0*pi*pi	! Gravitational constant [Mpc^3/yr^2/Msolar]
	real*8,dimension(3,27) :: offset			! For mirror particles
	
end module Constants

!--------------------------------------------------------------------------------------------------
integer function vari(i)
!--------------------------------------------------------------------------------------------------
! Function to return a varying integer to access desired part of array
! Doing this to tidy up code as this varying integer is regularly used
!--------------------------------------------------------------------------------------------------
	implicit none
	integer,intent(in) :: i			! Takes in iteration integer i
	
	vari = 3*(i-1)						! Outputted integer varies with inputted integer i

end function vari

!--------------------------------------------------------------------------------------------------
integer function N3vari(i)
!--------------------------------------------------------------------------------------------------
! Same as above function vari(i), but adds on N3=3*N
! This is done to access the second half of an array (e.g. velocities in the y-array)
! Doing this to tidy up code as this varying integer is regularly used
!--------------------------------------------------------------------------------------------------
	use Constants, only: N3
	implicit none
	integer,intent(in) :: i			! Takes in iteration integer i
	
	N3vari = N3+3*(i-1)				! Outputted integer varies with inputted integer i

end function N3vari

!--------------------------------------------------------------------------------------------------
program NBody
!--------------------------------------------------------------------------------------------------
! The main program
!--------------------------------------------------------------------------------------------------
	use Constants, only: xyz,time,Ndim,L,offset
	implicit none
	real*8,dimension(Ndim) :: y0,y	! Solution vector y (initial, final)
	real*8 :: f(Ndim)						! RHS vector f
	real*8 :: O=0d0	! Dummy variable to store the value 0
	
	call random_seed()					! Set the seed for the RNG used later
	
	offset = reshape((/O,O,O,  L,O,O,  -L,O,O,   O,L,O, O,-L,O, O,O,L, O,O,-L,&
							 L,L,O,  L,-L,O,  L,O,L,   L,O,-L,&
							-L,L,O, -L,-L,O, -L,O,L,  -L,O,-L,&
							 O,L,L,  O,L,-L,  O,-L,L,  O,-L,-L,&
							 L,L,L,  L,-L,L,  L,L,-L,  L,-L,-L,&
							-L,L,L, -L,-L,L, -L,L,-L, -L,-L,-L/), (/3,27/))
	
	call Set_Masses						! Setup the array of masses
	call init_cond(y0)					! Setup initial values for y
	call Orbit(y0,y,f)					! Determine the orbital paths of each body

end program NBody

!--------------------------------------------------------------------------------------------------
subroutine Set_Masses
!--------------------------------------------------------------------------------------------------
! Subroutine to define the masses of each body and store them in an array for later use
!--------------------------------------------------------------------------------------------------
	use Constants, only: N,M,Mtot
	implicit none
	real*8 :: m1,m2		! Masses [Msolar]
	integer :: i
	real*8,external :: RNG
	
	! Store each mass in an array for later use
	M = (/1d1,-5d0/)
	! Total mass of the system [Msolar]
	Mtot = sum(M)
	
end subroutine Set_Masses

!--------------------------------------------------------------------------------------------------
subroutine init_cond(y0)
!--------------------------------------------------------------------------------------------------
! Subroutine to set the initial conditions
!--------------------------------------------------------------------------------------------------
	use Constants, only: N,Ndim,N3,M,Mtot,L
	implicit none
	real*8,dimension(Ndim),intent(out) :: y0	! Solution vector y
	integer,external :: vari,N3vari
	real*8,external  :: RNG										! Use RNG to set initial conditions
	real*8  :: rmin,rmax,vmin,vmax				! For use with the RNG
	real*8, dimension(3) :: rc,vc					! Co-moving frame vectors
	integer :: i

	y0 = (/	-L/4d0,L/2d1,L/2d1, -L/4d0,-L/2d1,-L/2d1,&
				L/8d0, 0d0,  0d0, 	L/10d0,0d0,   0d0		/)

end subroutine init_cond

!--------------------------------------------------------------------------------------------------
real*8 function RNG(minimum,maximum)
!--------------------------------------------------------------------------------------------------
! RNG function that returns a random number in the interval [minimum,maximum)
! Since the intrinsic RNG function only works with the interval [0,1), must scale to produce a pseudorandom number within any specified interval
! RNG scaling taken from http://infohost.nmt.edu/tcc/help/lang/fortran/scaling.html
!--------------------------------------------------------------------------------------------------
	implicit none
	real*8,intent(in) :: minimum,maximum	! Desired range
	real*8 :: num	! Local dummy variables
	
	! The RNG
	call random_number(num)
	! Scale the randomly generated number
	RNG = (num*(maximum-minimum))+minimum

end function RNG

!--------------------------------------------------------------------------------------------------
subroutine RK4(y,t,tstep,r,v)
!--------------------------------------------------------------------------------------------------
! ***  fourth-order Runge-Kutta integrator                                ***
!--------------------------------------------------------------------------------------------------
  use Constants, only: N,Ndim
  implicit none
  real*8,intent(inout) :: y(Ndim)      ! the solution vector at t, y(t)
                                       ! on return, it contains y(t+dt)
  real*8,intent(in) :: t,tstep         ! time t, and time step to be integrated
  real*8,dimension(3,N),intent(inout) :: r,v
  real*8 :: h                          ! half time step
  real*8,dimension(Ndim) :: f1, f2, f3, f4, ytmp1, ytmp2, ytmp3

  h=0.5d0*tstep

  call RHS(y,f1,r,v)             ! calulate f(y,t)
  ytmp1 = y + h*f1                     ! half step

  call RHS(ytmp1,f2,r,v)       	! calculate f(ytmp1,t+h)
  ytmp2 = y + h*f2                     ! half step in another way

  call RHS(ytmp2,f3,r,v)       	! calculate f(ytmp2,t+h)
  ytmp3 = y + tstep*f3                 ! full step

  call RHS(ytmp3,f4,r,v)   		! calculate f(ytmp3,t+tstep)

  y = y + (f1 + 2.d0*f2 + 2.d0*f3 + f4)*tstep/6.d0    ! Runge-Kutta recipe

end subroutine RK4

!--------------------------------------------------------------------------------------------------
subroutine Boundary(y)
!--------------------------------------------------------------------------------------------------
! Subroutine to apply a periodic boundary condition to contain particles within a box of length L
!--------------------------------------------------------------------------------------------------
	use Constants, only: N,Ndim,L
	implicit none
	real*8,intent(inout) :: y(Ndim)	! Solution vector is checked and possibly adjusted
	integer :: i=0
	integer,external :: vari
	real*8 :: BC=L/2d0					! Boundary Condition (edge of box)
	
	do i=1,N
		! Check solution vector y against the boundary conditions for each particle
		if (y(1+vari(i))>BC) then			! Check x-coordinate
			y(1+vari(i))=y(1+vari(i))-L	! Make adjustment to put inside box
		elseif (y(1+vari(i))<-BC) then
			y(1+vari(i))=y(1+vari(i))+L
		endif
		
		if (y(2+vari(i))>BC) then			! Check y-coordinate
			y(2+vari(i))=y(2+vari(i))-L	! Make adjustment to put inside box
		elseif (y(2+vari(i))<-BC) then
			y(2+vari(i))=y(2+vari(i))+L
		endif
		
		if (y(3+vari(i))>BC) then			! Check z-coordinate
			y(3+vari(i))=y(3+vari(i))-L	! Make adjustment to put inside box
		elseif (y(3+vari(i))<-BC) then
			y(3+vari(i))=y(3+vari(i))+L
		endif
	enddo

end subroutine Boundary

!--------------------------------------------------------------------------------------------------
subroutine RHS(y,f,r,v)
!--------------------------------------------------------------------------------------------------
! Subroutine to calculate the RHS vector f
!--------------------------------------------------------------------------------------------------
	use Constants, only: N,Ndim,N3,M
	implicit none
	real*8,intent(in)  :: y(Ndim)			! Solution vector y is input
	real*8,intent(out) :: f(Ndim)			! RHS vector f is output
	real*8,dimension(3,N),intent(inout) :: r,v
	integer,external   :: vari,N3vari
	real*8,dimension(N3) :: a				! Array of accelerations on each object
	integer :: i
	
	! Calculate the forces on each body
	call Accelerations(y,a,r,v)
	
	! Determine the RHS vector f
	do i=1,N
		f(1+vari(i):3+vari(i)) = y(1+N3vari(i):3+N3vari(i))			! Velocities [Mpc/yr]
		f(1+N3vari(i):3+N3vari(i)) = a(1+vari(i):3+vari(i))	! Accelerations [Mpc/yr^2]
	enddo

end subroutine RHS

!--------------------------------------------------------------------------------------------------
subroutine Accelerations(y,a,r,v)
!--------------------------------------------------------------------------------------------------
! Subroutine to calculate the gravitational accelerations of each body
! Acceleration calculated in [Mpc/yr^2]
!--------------------------------------------------------------------------------------------------
	use Constants, only: N,Ndim,N3,G,M,L,offset
	implicit none
	integer,external :: vari,N3vari
	real*8,external :: mag_rirj
	real*8,intent(in),dimension(Ndim) :: y		! Solution vector y contains positions [Mpc]
	real*8,intent(out) :: a(N3)					! Array of accelerations on each object [Mpc/yr^2]
	real*8,dimension(3,N),intent(out) :: r,v	! Position and velocity vectors
	integer :: i,j,k
	real*8 :: inert	! To account for sign of inertial mass

	
	! Store the positions and velocities of each body for easy reference
	do i=1,N
		r(:,i) = y(1+vari(i):3+vari(i))			! (x,y,z) of Body i [Mpc]
		v(:,i)  = y(1+N3vari(i):3+N3vari(i))	! (vx,vy,vz) of Body i [Mpc/yr]
	enddo

	! Calculate the resulting gravitational accelerations of each body [Msolar*Mpc/yr^2]
	! Calculated via force on body j due to each other body i/=j 
	a = 0.d0
	do j=1,N
		if (M(j)<0d0) inert=-1d0
		if (M(j)>0d0) inert=+1d0
		do i=1,N
			if (i==j) cycle
			do k=1,27
				a(1+vari(j):3+vari(j)) = a(1+vari(j):3+vari(j))+inert*G*M(i)*((r(:,i)+offset(:,k))-r(:,j))/mag_rirj(r,i,j,k)**3
			enddo
		enddo
	enddo
	
end subroutine Accelerations

!--------------------------------------------------------------------------------------------------
subroutine Orbit(y0,y,f)
!--------------------------------------------------------------------------------------------------
! Subroutine to calculate the orbital paths of each body and write them to a data file to plot
!--------------------------------------------------------------------------------------------------
	use Constants, only: N,sel,xyz,time,Ndim,N3,L,Mtot,G
	implicit none
	real*8,intent(inout),dimension(Ndim) :: y0,y,f
	real*8,dimension(3,N) :: r,v				! Position and velocity vectors
	real*8 :: rho									! Density [Msolar/Mpc^3]
	real*8 :: tdyn									! The dynamical timescale [yr]
	real*8 :: ttot									! Total time to integrate over [yr]
	real*8 :: dt									! Timestep [yr]
	integer :: i=0
	real*8 :: t=0d0								! Time at current iteration [yr]
	
	y=y0	 											! Initialise position vector y
	call RHS(y0,f,r,v)							! Initialise f
	
	! Set the total length of time to integrate over
	rho  = Mtot/L**3			![Msolar/Mpc^3]
	tdyn = 1d0/sqrt(rho*G)	![yrs]
	ttot = 5d0*tdyn			![yrs]
	open(1,file=trim(xyz))
	open(2,file=trim(time))
	write(2,*) tdyn
	do while (t<=ttot)
		if (mod(i,sel)==0) then
			print '(F7.3,"%")', t/ttot *100.d0	! Display % complete
			write(1,*) y(1:N3)	! Write coordinates to file at certain timesteps
			write(2,*) t
		endif
		call Timestep(r,v,dt)					! Adjust the timestep
		call RK4(y,t,dt,r,v)						! Integrate y using the RK4 technique
		call Boundary(y)							! Check against boundary conditions
		t=t+dt
		i=i+1
	enddo
	close(2)
	close(1)

	print '(A)', 'File has been written'	! To indicate completion
	print '(A,ES16.10,A)', 'rho  = ',rho,' Msolar/AU^3',&
								  'G    = ',G,' AU^3/yr^2/Msolar',&
								  'tdyn = ',tdyn,' yr',&
								  'ttot = ',ttot,' yr'

end subroutine Orbit

!--------------------------------------------------------------------------------------------------
subroutine Timestep(r,v,dt)
!--------------------------------------------------------------------------------------------------
! Subroutine to adjust the timestep during the simulation
!--------------------------------------------------------------------------------------------------
	use Constants, only: alpha,N
	implicit none
	real*8,intent(in),dimension(3,N) :: r,v	! Position and velocity vectors of each body
	real*8,intent(out) :: dt						! Adjusted timestep 
	real*8,external :: mag_rirj,mag_vivj
	real*8,dimension(N,N) :: magr,magv,ratio	! Arrays to be used in this subroutine
	integer :: i,j
	 
	magr=0.d0;magv=0.d0;ratio=0.d0;
	do i=1,N
		do j=1,N
			if (j==i) cycle
			magr(j,i) = mag_rirj(r,i,j,1)
			magv(j,i) = mag_vivj(v,i,j)
			ratio(j,i) = magr(j,i)/magv(j,i)
		enddo
	enddo
	dt = alpha*minval(ratio,mask=ratio>0.d0) ! Mask ensures diagonal of matrix is ignored

end subroutine Timestep

!--------------------------------------------------------------------------------------------------
real*8 function mag_rirj(r,i,j,k)
!--------------------------------------------------------------------------------------------------
! Function to calculate the magnitude of ri-rj, where ri and rj are position vectors for bodies i and j
!--------------------------------------------------------------------------------------------------
	use Constants, only: N,L,offset
	implicit none
	real*8,dimension(3,N),intent(in) :: r	
	integer,intent(in) :: i,j,k	! Indices
	real*8,dimension(3) :: rirj
		
	rirj = (r(:,i)+offset(:,k))-r(:,j)
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
	integer,intent(in) :: i,j		! Index, to select the relevant bodies
	real*8,dimension(3) :: vivj
	
	vivj = v(:,i)-v(:,j)
	mag_vivj = sqrt(vivj(1)**2+vivj(2)**2+vivj(3)**2)
	
end function mag_vivj

!--------------------------------------------------------------------------------------------------
real*8 function magv2(y,i)
!--------------------------------------------------------------------------------------------------
! Function to calculate the magnitude squared of the velocity vector for body i
!--------------------------------------------------------------------------------------------------
	use Constants, only: Ndim
	implicit none
	real*8,dimension(Ndim),intent(in) :: y	! Make use of velocities stored in y
	integer,intent(in) :: i						! Index, to select the relevant bodies
	integer,external :: vari,N3vari
	
	magv2 = y(1+N3vari(i))**2 + y(2+N3vari(i))**2 + y(3+N3vari(i))**2
	
end function magv2
