program test
	implicit none
	real*8,dimension(10) :: s,M
	integer :: i
	
	do i=1,10
		call random_number(M(i))
		if (mod(i,2)==0) M(i)=-M(i)
	enddo
	
	s=dsign(1d0,M)
	print'(5F5.1)',s

end program test	
