!--------------------------------------------------------------------------------------------------
module Constants
!--------------------------------------------------------------------------------------------------
! Module to keep constants in
!--------------------------------------------------------------------------------------------------
	implicit none
	integer,parameter :: N=2					! Number of bodies
	integer,parameter :: k=1					! Locate the kth nearest particle
	integer,parameter :: Np=2d3				! Number of data points to write
	real*8,parameter :: Ntdyn=2d0				! Number of dynamical timescales to iterate over
	character(len=99),parameter :: nam="Plot"	! Start of File Name
	character(len=99),parameter :: dat=trim(nam)//"_Data.dat"	! Name of file to save data to
	character(len=99),parameter :: mass=trim(nam)//"_Signs.dat"	! File to save mass signs to
	real*8,parameter :: alpha=1d-3			! Dimensionless parameter for adjusted time step
	real*8,parameter :: L=3d2					! Arbitrary length [Mpc] of box for periodic B.C.
	integer,parameter :: Ndim=3*2*N			! Dimension of solution vector y
	integer,parameter :: N3=3*N				! Integer used in iterative loops
	real*8,dimension(N) :: M,MP,MA,MI		! Store masses in an array
	real*8 :: Mtot									! Total mass of system [Msolar]
	real*8,parameter :: pi=acos(-1.d0)		! The value of pi
	real*8,parameter :: Mpc=4.84814d-12		! Mpc per AU
	real*8,parameter :: AU=1.495978700d8	! km per AU
	real*8,parameter :: yr=3.155760000d7	! s per yr
	real*8,parameter :: G=(Mpc**3)*4.d0*pi*pi	! Gravitational constant [Mpc^3/yr^2/Msolar]
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
	use Constants, only: Ndim,L,offset
	implicit none
	real*8,dimension(Ndim) :: y0,y	! Solution vector y (initial, final)
	real*8 :: f(Ndim)						! RHS vector f
	! Dummy Variables
	real*8 :: O=0d0						! Just 0
	integer,allocatable :: S(:)		! For the RNG seed
	integer :: n							! For the RNG seed
	
	! Set the seed for the RNG used later
	call random_seed(size=n)			! Obtain seed size
	allocate(S(n))							! Ready the array to store the seed
	open(1,file="RNG_Seed.dat")	
	read(1,*) S								! Load seed from file
	close(1)
	call random_seed(PUT=S)				! Assign seed
	
	! Define an array to use when generating mirror copies of the box
	offset = reshape((/O,O,O,  L,O,O,  -L,O,O,   O,L,O, O,-L,O, O,O,L, O,O,-L,&
							 L,L,O,  L,-L,O,  L,O,L,   L,O,-L,&
							-L,L,O, -L,-L,O, -L,O,L,  -L,O,-L,&
							 O,L,L,  O,L,-L,  O,-L,L,  O,-L,-L,&
							 L,L,L,  L,-L,L,  L,L,-L,  L,-L,-L,&
							-L,L,L, -L,-L,L, -L,L,-L, -L,-L,-L/), (/3,27/))
	
	call Set_Masses						! Setup the array of masses		
	call init_cond(y0)					! Setup initial values for y
	call Orbit(y0,y,f)					! Determine the orbital paths of each body

	! Call python to plot the coordinate data
!	call SYSTEM("ipython Spacetime.py")

end program NBody

!--------------------------------------------------------------------------------------------------
subroutine Set_Masses
!--------------------------------------------------------------------------------------------------
! Subroutine to define the masses of each body and store them in an array for later use
!--------------------------------------------------------------------------------------------------
	use Constants, only: N,mass,L,M,MP,MA,MI,Mtot
	implicit none
	real*8 :: vol					! Volume occupied by each particle
	integer :: i					! Iteration integer
	real*8,external :: RNG		! The RNG function being used
	real*8,dimension(N) :: s	! Dummy variable
	
	vol = (L*L*L)/dble(N)
	
	! Assign a value to each mass [Msolar]
	M = vol*1d12
	! Set half the masses negative
!	M(1+N/2:) = -M(1+N/2:)
	! Assign Masses according to desired signs
	! To not violate equivalency principle, MP and MI must have same sign
	MP = M				! Passive Gravitational Mass
	MA = M				! Active Gravitational Mass
	! Use dsign to automatically ensure equivalency principle is obeyed
	MI = dsign(M,MP)	! Inertial Mass
	! Store the sign of each mass for python to reference
	s=dsign(1d0,M)


	! Total mass of the system [Msolar]
	Mtot = sum(abs(M))
	
	! Write a file for python to read for signs of masses
	open(1,file=trim(mass))
	write(1,*) s
	close(1)
	
end subroutine Set_Masses

!--------------------------------------------------------------------------------------------------
subroutine init_cond(y0)
!--------------------------------------------------------------------------------------------------
! Subroutine to set the initial conditions
!--------------------------------------------------------------------------------------------------
	use Constants, only: N,Ndim,N3,M,Mtot,L
	implicit none
	real*8,dimension(Ndim),intent(out) :: y0	! Solution vector y
	integer,external :: vari,N3vari				! Iteration functions being used
	real*8,external  :: RNG							! Use RNG to set initial conditions
	real*8  :: rmin,rmax,vmin,vmax				! Lower/Upper bounds for RNG
	real*8, dimension(3) :: rc,vc					! Co-moving frame vectors
	integer :: i										! Iteration integer
	
	! Set the min and max desired values
	rmin=-L/2d0; rmax=L/2d0
	vmin=0d0; vmax=vmin
	! Initialise the co-moving frame vectors
	rc = 0.d0;vc = 0.d0
	! Pseudo-randomly set initial values for the solution vector y(t=0)
	do i=1,N
		! Select desired distribution
		goto 20
		! Seperate +ve and -ve Mass particles
10		if (M(i) .gt. 0d0) then 
			y0(1+vari(i)) = RNG(0d0,rmax)		! x-positions
			y0(2+vari(i)) = RNG(0d0,rmax)		! y-positions
			y0(3+vari(i)) = RNG(0d0,rmax)		! z-positions
			! Options for desired seperation			
!			y0(1+vari(i)) = RNG(rmin,rmax)	! x-positions
!			y0(2+vari(i)) = RNG(rmin,rmax)	! y-positions
!			y0(3+vari(i)) = RNG(rmin,rmax)	! z-positions
		else
			y0(1+vari(i)) = RNG(rmin,0d0)		! x-positions
			y0(2+vari(i)) = RNG(rmin,0d0)		! y-positions
			y0(3+vari(i)) = RNG(rmin,0d0)		! z-positions
			! Options for desired seperation			
!			y0(1+vari(i)) = RNG(rmin,rmax)	! x-positions
!			y0(2+vari(i)) = RNG(rmin,rmax)	! y-positions
!			y0(3+vari(i)) = RNG(rmin,rmax)	! z-positions
		endif
		goto 40
		
		! Mix +ve and -ve mass particles together
20		y0(1+vari(i)) = RNG(rmin,rmax)		! x-positions
		y0(2+vari(i)) = RNG(rmin,rmax)		! y-positions
		y0(3+vari(i)) = RNG(rmin,rmax)		! z-positions
		goto 40
		
		! Fix positions of particles
30		goto 40
				
		! Assign initial velocities
40		y0(1+N3vari(i)) = RNG(vmin,vmax)		! x-velocities
		y0(2+N3vari(i)) = RNG(vmin,vmax)		! y-velocities
		y0(3+N3vari(i)) = RNG(vmin,vmax)		! z-velocities

		! Generate the co-moving frame coordinates and velocities
!		rc = rc+(abs(M(i))/Mtot)*y0(1+vari(i):3+vari(i))				! Positions  [Mpc]
!		vc = vc+(abs(M(i))/Mtot)*y0(1+N3vari(i):3+N3vari(i))			! Velocities [Mpc/yr]
	enddo
	! Transform to co-moving frame
!	do i=1,N
!		y0(1+vari(i):3+vari(i))=y0(1+vari(i):3+vari(i))-rc				! Position correction [Mpc]
!		y0(1+N3vari(i):3+N3vari(i))=y0(1+N3vari(i):3+N3vari(i))-vc	! Velocity correction [Mpc/yr]
!	enddo
	
end subroutine init_cond

!--------------------------------------------------------------------------------------------------
real*8 function RNG(minimum,maximum)
!--------------------------------------------------------------------------------------------------
! RNG function that returns a random number in the interval [minimum,maximum)
! Since the intrinsic RNG function only works with the interval [0,1), must scale to produce a pseudo-random number within any specified interval
!--------------------------------------------------------------------------------------------------
	implicit none
	real*8,intent(in) :: minimum,maximum	! Desired range
	real*8 :: num									! Local dummy variable
	
	! The RNG
	call random_number(num)
	! Scale the output to lie within the desired interval
	RNG = (num*(maximum-minimum))+minimum

end function RNG

!--------------------------------------------------------------------------------------------------
subroutine RK4(y,t,tstep,r,v,z)
!--------------------------------------------------------------------------------------------------
! ***  fourth-order Runge-Kutta integrator                                ***
!--------------------------------------------------------------------------------------------------
  use Constants, only: N,Ndim
  implicit none
  real*8,intent(inout) :: y(Ndim)      ! the solution vector at t, y(t)
                                       ! on return, it contains y(t+dt)
  real*8,intent(in) :: t,tstep         ! time t, and time step to be integrated
  real*8,dimension(3,N),intent(inout) :: r,v
  real*8,intent(in) :: z					! Redshift
  real*8 :: h                          ! half time step
  real*8,dimension(Ndim) :: f1, f2, f3, f4, ytmp1, ytmp2, ytmp3

  h=0.5d0*tstep

  call RHS(y,f1,r,v,z)             		! calculate f(y,t)
  ytmp1 = y + h*f1                     ! half step

  call RHS(ytmp1,f2,r,v,z)       		! calculate f(ytmp1,t+h)
  ytmp2 = y + h*f2                     ! half step in another way

  call RHS(ytmp2,f3,r,v,z)       		! calculate f(ytmp2,t+h)
  ytmp3 = y + tstep*f3                 ! full step

  call RHS(ytmp3,f4,r,v,z)   				! calculate f(ytmp3,t+tstep)

  y = y + (f1 + 2.d0*f2 + 2.d0*f3 + f4)*tstep/6.d0    ! Runge-Kutta recipe

end subroutine RK4

!--------------------------------------------------------------------------------------------------
subroutine Boundary(y)
!--------------------------------------------------------------------------------------------------
! Subroutine to apply a periodic boundary condition to contain particles within a box of length L
!--------------------------------------------------------------------------------------------------
	use Constants, only: N,Ndim,L
	implicit none
	real*8,intent(inout) :: y(Ndim)		! Solution vector is checked and possibly adjusted
	integer :: i
	integer,external :: vari
	real*8 :: BC=L/2d0						! Boundary Condition (edge of box)
	integer :: xind,yind,zind				! Dummy variables
	
	do i=1,N
		! Position indeces to reference
		xind = 1+vari(i)
		yind = 2+vari(i)
		zind = 3+vari(i)
		! Check solution vector y against the boundary conditions for each particle
		! Check x-coordinate
		do while (y(xind) .gt. BC)
			y(xind)=y(xind)-L
		enddo
		do while (y(xind) .lt.-BC)
			y(xind)=y(xind)+L
		enddo
		! Check y-coordinate
		do while (y(yind) .gt. BC)
			y(yind)=y(yind)-L
		enddo
		do while (y(yind) .lt. -BC)
			y(yind)=y(yind)+L
		enddo
		! Check z-coordinate
		do while (y(zind) .gt. BC)
			y(zind)=y(zind)-L
		enddo
		do while (y(zind) .lt. -BC)
			y(zind)=y(zind)+L
		enddo
	enddo

end subroutine Boundary

!--------------------------------------------------------------------------------------------------
subroutine RHS(y,f,r,v,z)
!--------------------------------------------------------------------------------------------------
! Subroutine to calculate the RHS vector f
!--------------------------------------------------------------------------------------------------
	use Constants, only: N,Ndim,N3
	implicit none
	real*8,intent(in)  :: y(Ndim)			! Solution vector y is input
	real*8,intent(out) :: f(Ndim)			! RHS vector f is output
	real*8,dimension(3,N),intent(inout) :: r,v
	real*8,intent(in)  :: z					! Redshift
	real*8,dimension(N3) :: acc			! Array of accelerations on each object
	
	! Calculate the acceleration of each body
	call Accelerations(y,acc,r,v,z)
	
	! Set the RHS vector f
	f(1:N3)  = y(1+N3:)	! Velocities
	f(1+N3:) = acc			! Accelerations
	
end subroutine RHS

!--------------------------------------------------------------------------------------------------
subroutine Accelerations(y,acc,r,v,z)
!--------------------------------------------------------------------------------------------------
! Subroutine to calculate the gravitational accelerations of each body
! Acceleration calculated in [Mpc/yr^2]
!--------------------------------------------------------------------------------------------------
	use Constants, only: N,Ndim,N3,G,MP,MA,MI,L,offset
	implicit none
	real*8,intent(in),dimension(Ndim) :: y		! Solution vector y contains positions [Mpc]
	real*8,intent(out) :: acc(N3)					! Array of accelerations on each object [Mpc/yr^2]
	real*8,dimension(3,N),intent(out) :: r,v	! Position and velocity vectors
	real*8,intent(in) :: z							! Redshift
	integer,external :: vari,N3vari
	real*8,external :: mag_rirj,H,a
	integer :: i,j,k									! Iteration integers
	real*8 :: GMj										! Ratio of MP/MI for body j, multiplied by G
	real*8 :: soft										! Softening Parameter
	real*8 :: magrirj2								! The magnitude squared of ri-rj
	real*8 :: dist3									! The softened distance between particles cubed
	real*8,dimension(3) :: rirj					! Vector pointing from rj to ri (ri-rj)
	real*8,dimension(3) :: grav					! Gravitational acceleration term
	real*8,dimension(3) :: drag					! Hubble Drag term

	
	! Store the positions and velocities of each body for easy reference
	do i=1,N
		r(:,i) = y(1+vari(i):3+vari(i))			! (x,y,z) of Body i [Mpc]
		v(:,i) = y(1+N3vari(i):3+N3vari(i))		! (vx,vy,vz) of Body i [Mpc/yr]
	enddo

	soft = 0.98*dble(N)**(-0.28)

	! Calculate the resulting gravitational accelerations of each body [Msolar*Mpc/yr^2]
	! Calculated via force on body j due to each other body i!=j 
	acc = 0.d0
	do j=1,N
		GMj = MP(j)/MI(j)*G
		grav=0d0
		do i=1,N
			if (i .eq. j) cycle
			do k=1,27
				rirj = (r(:,i)+offset(:,k))-r(:,j)
				magrirj2 = mag_rirj(r,i,j,k)*mag_rirj(r,i,j,k)
				if(magrirj2 .le. soft) then
					dist3 = (magrirj2+soft*soft)**(1.5)	! (Mag(ri-rj)^2 + Epsilon^2)^(3/2)
				else
					dist3 = magrirj2**1.5
				endif
				grav = grav+MA(i)*rirj/dist3
			enddo
		enddo
		grav = grav*GMj/a(z)/a(z)
		drag = -2d0*H(z)*v(:,j)
		acc(1+vari(j):3+vari(j)) = acc(1+vari(j):3+vari(j))+grav+drag
	enddo
	
end subroutine Accelerations

!--------------------------------------------------------------------------------------------------
subroutine Orbit(y0,y,f)
!--------------------------------------------------------------------------------------------------
! Subroutine to calculate the orbital paths of each body and write them to a data file to plot
!--------------------------------------------------------------------------------------------------
	use Constants, only: N,Np,Ntdyn,dat,mass,Ndim,N3,L,Mtot,G
	implicit none
	real*8,intent(inout),dimension(Ndim) :: y0,y,f
	real*8,external :: dz,a
	real*8,dimension(3,N) :: r,v				! Position and velocity vectors
	real*8 :: rhobar								! Average density [Msolar/Mpc^3]
	real*8,dimension(N) :: rho					! Local density [Msolar/Mpc^3]
	real*8 :: tdyn									! The dynamical timescale [yr]
	real*8 :: ttot									! Total time to integrate over [yr]
	real*8 :: dt									! Time step [yr]
	integer:: i=0
	real*8 :: t=0d0								! Time at current iteration [yr]
	real*8 :: dtmax,dtmin						! Timestep upper/lower bounds
	real*8 :: z=1.1d3								! Redshift
	real*8,dimension(N) :: delta				! Ratio between rho and rhobar
	real*8,dimension(N,N) :: magr				! Store distances to each other body
	
	y=y0/a(z)										! Initialise position vector y
	! Make sure every particle is IN the box
	call boundary(y)
	call RHS(y,f,r,v,z)							! Initialise f
	
	! Set the total length of time to integrate over
	rhobar  = Mtot/L/L/L			![Msolar/Mpc^3]
	tdyn = 1d0/sqrt(rhobar*G)	![yrs]
	ttot = Ntdyn*tdyn				![yrs]
	dtmax= tdyn/1d5				![yrs]
	dtmin= 1d2						![yrs]
								  
	open(1,file=trim(dat))
	write(1,*) N,L,tdyn
	do while (t .le. ttot)
		call Timestep(r,v,dt,dtmin,magr)		! Adjust the time step
		if (dt .gt. dtmax) dt=dtmax
		if (dt .lt. dtmin) dt=dtmin
		if (t+dt .gt. ttot .AND. t .ne. ttot) dt=ttot-t
		if (t .ge. i*ttot/Np) then				! Write per specified time interval
			call LocalDensity(r,rho,magr)		! Calculate local densities
			delta=rho/rhobar
			print '(F7.3,"%",A,F5.3,A)', t/ttot*1d2,' - t = ',t/tdyn,' tdyn'	! Display % complete
			write(1,*) y(1:N3),t,dt,delta		! Write data to file
			if (mod(i,Np/10) .eq. 0 .AND. i .ne. 0) call SYSTEM("python Spacetime.py")
			i=i+1
		endif
		call RK4(y,t,dt,r,v,z)					! Integrate y using the RK4 technique
		call Boundary(y)							! Check against boundary conditions
		t=t+dt
		z=z+dz(z,dt)
	enddo
	close(1)
	
	! Indicate completion
	print '(A)', 'Files Created:',&
					 trim(dat),&
					 trim(mass)
	print*
	print '(A,ES16.10,A)', 'rhobar = ',rhobar,' Msolar/Mpc^3',&
								  'G      = ',G,		' Mpc^3/yr^2/Msolar',&
								  'tdyn   = ',tdyn,	' yr',&
								  'ttot   = ',ttot,	' yr'
	
end subroutine Orbit

!--------------------------------------------------------------------------------------------------
subroutine Timestep(r,v,dt,dtmin,magr)
!--------------------------------------------------------------------------------------------------
! Subroutine to adjust the time step during the simulation
!--------------------------------------------------------------------------------------------------
	use Constants, only: alpha,N
	implicit none
	real*8,intent(in),dimension(3,N) :: r,v	! Position and velocity vectors of each body
	real*8,intent(out) :: dt						! Adjusted time step 
	real*8,intent(in) :: dtmin						! Lower timestep bound
	real*8,intent(out),dimension(N,N) :: magr	! Dummy variable
	real*8,external :: mag_rirj,mag_vivj
	real*8 :: magv										! Dummy variable
	real*8,dimension(N,N) :: ratio				! Array to be used in this subroutine
	integer :: i,j
	 
	magr=0.d0;magv=0.d0;ratio=0.d0;
	do i=1,N
		do j=1,N
			if (j .eq. i) cycle
			magr(j,i) = mag_rirj(r,i,j,1)
			magv = mag_vivj(v,i,j)
			if (magv .eq. 0d0) then
				ratio(j,i) = dtmin	! Prevent a division by 0
			else
				ratio(j,i) = magr(j,i)/magv
			endif
		enddo
	enddo
	dt = alpha*minval(ratio,mask=ratio .gt. 0.d0) ! Mask ensures diagonal of matrix is ignored

end subroutine Timestep

!--------------------------------------------------------------------------------------------------
subroutine LocalDensity(r,rho,magr)
!--------------------------------------------------------------------------------------------------
! Subroutine that returns the "local density" of each body
! This is represented by a sphere containing the k nearest bodies to each body
!--------------------------------------------------------------------------------------------------
	use Constants, only: N,k,M,pi
	implicit none
	real*8,dimension(3,N),intent(in) :: r	! Position vectors of each body
	real*8,dimension(N),intent(out) :: rho	! Local density for each particle
	real*8,dimension(N,N),intent(in) :: magr	! Dummy variable
	real*8 :: Mloc							! Local Total Mass
	real*8 :: rsphere						! Radius of sphere containing the nth nearest body
	real*8 :: V								! Local Volume
	real*8,allocatable :: search(:),mass(:)	! Dummy arrays for minimum search
	real*8 :: ind							! Dummy variable to store an index
	integer :: i,j							! Iteration integers

	! Search for kth nearest body
	do i=1,N
		! Initialise the search for body i's kth nearest body
		search=(/magr(1:i-1,i),magr(i+1:,i)/)
		mass=(/M(1:i-1),M(i+1:)/)
		Mloc=M(i)
		do j=1,k
			! Index of current minimum
			ind=MINLOC(search,1)
			! Add this body's mass to the local total
			Mloc=Mloc+mass(nint(ind))
			if (j .eq. k) exit	! Don't remove the kth element from search list
			! Remove this minimum from search to find next minimum
			search=(/search(1:nint(ind)-1),search(nint(ind)+1:)/)
			mass=(/mass(1:nint(ind)-1),mass(nint(ind)+1:)/)
		enddo
		! Calculate the local density for body i
		rsphere=MINVAL(search)
		V=4d0/3d0*pi*rsphere*rsphere*rsphere	! Volume of a sphere
		rho(i)=Mloc/V
	enddo
		
end subroutine LocalDensity

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
	mag_rirj = sqrt(rirj(1)*rirj(1)+rirj(2)*rirj(2)+rirj(3)*rirj(3))
	
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
	mag_vivj = sqrt(vivj(1)*vivj(1)+vivj(2)*vivj(2)+vivj(3)*vivj(3))
	
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
	
	magv2=y(1+N3vari(i))*y(1+N3vari(i))+y(2+N3vari(i))*y(2+N3vari(i))+y(3+N3vari(i))*y(3+N3vari(i))
	
end function magv2

!--------------------------------------------------------------------------------------------------
real*8 function H(z)
!--------------------------------------------------------------------------------------------------
! Returns the Hubble parameter in [yr^-1] at the inputted redshift z
! Require H in these units for calculating dz from dt
!--------------------------------------------------------------------------------------------------
	use Constants,only:Mpc,AU,yr
	implicit none
	real*8,intent(in) :: z						! Cosmological redshift
	real*8,external  :: a						! Scale factor at redshift z
	real*8,parameter :: H0=7d1*Mpc/AU*yr	! Hubble Constant at present day [yr^-1]
	real*8,parameter :: OmL=0.7				! Dark energy density parameter
	real*8,parameter :: OmM=0.3				! Matter density parameter
	real*8,parameter :: Omk=1d0-OmL-OmM		! Curvature term
	
	H=H0*SQRT(OmL+(OmM*a(z)+Omk)*a(z)*a(z))	

end function H

!--------------------------------------------------------------------------------------------------
real*8 function dz(z,dt)
!--------------------------------------------------------------------------------------------------
! Function to calculate the change in redshift dz over the timestep dt
!--------------------------------------------------------------------------------------------------
	implicit none
	real*8,intent(in) :: z,dt	! Current redshift and timestep
	real*8,external :: H,a		! Hubble parameter, scale factor
	
	dz = -H(z)*dt/a(z)			! Change in redshift
	
end function dz

!--------------------------------------------------------------------------------------------------
real*8 function a(z)
!--------------------------------------------------------------------------------------------------
! Function to calculate the scale factor a at redshift z
!--------------------------------------------------------------------------------------------------
	implicit none
	real*8,intent(in) :: z		! Redshift is input
	
	a = 1d0/(1d0+z)

end function

