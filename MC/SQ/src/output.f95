module output
 use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
 
 implicit none
 
 
 contains
 
 subroutine OUTPUTPOS(iCycle,NPart,TYP,X,Y)
	implicit none
	
	integer, intent(in) :: iCycle   ! Input, Cycle we are on
	integer, intent(in) :: NPart	! Input, Total number of particles
	
	integer, intent(in) :: TYP(:)	! Input, Type of Each Particle
	
	real(dp), intent(in) :: X(:) 	! Input, X Positions
 	real(dp), intent(in) :: Y(:)	! Input, Y Positions
	
	character(len=1024) :: F1		! Internal, Format description
	character(len=1024) :: OUTFILE  ! Internal, Name of File, changes based on iCycle
	character(len=1024) :: iOUT     ! Internal, iCycle converted to string
	
	character(len=10) :: NOUT        ! Internal, convert to string before writing to file
	character(len=10) :: TOUT	    ! Internal, convert to string before writing to file
	character(len=15) :: XOUT	    ! Internal, convert to string before writing to file
	character(len=15) :: YOUT        ! Internal, convert to string before writing to file
	
	integer :: iPart				! Internal, loop variable
	
	F1 = '(I20)'
	
	write(iOUT,F1) iCycle
	
	OUTFILE = 'out/YK'//trim(adjustl(iOUT))//'.xyz'
	
	open(1, file = OUTFILE, status = 'new')
	
	write(NOUT,'(I10)') NPart
	write(1,'(a8)') adjustl(NOUT)
	write(1,*) 
	do iPart= 1,NPart
		!write(TOUT, '(i8)') TYP(iPart)
		write(XOUT, '(F10.5)') X(iPart)
		write(YOUT, '(F10.5)') Y(iPart)
		if (TYP(iPart) .eq. 1) then
			write(TOUT, '(a8)') 'O'
		else if ( TYP(iPart) .eq. 2) then
			write(TOUT, '(a8)') 'I'
		else
			write(TOUT, '(a8)') 'S'
		end if
		write(1,'(4a10)') adjustl(TOUT), adjustl(XOUT), adjustl(YOUT), '0.00000'
		

	end do
	close(1)	
	
 end subroutine OUTPUTPOS

 subroutine OUTPUTDAT(iCycle,Usys,Area,DX,DA,ProbD,ProbS,ProbA)
	implicit none
	
	integer, intent(in) :: iCycle		! Input, Cycle of output
	
	real(dp), intent(in) :: Usys        ! Input, Energy of system
	real(dp), intent(in) :: Area        ! Input, Area of system
	
	real(dp), intent(in) :: DX          ! Input, displacement move step size
	real(dp), intent(in) :: DA          ! Input, area move step size
	
	real(dp), intent(in) :: ProbD       ! Input, Acceptance Probability of displacement moves
	real(dp), intent(in) :: ProbS       ! Input, Acceptance Probability of swap moves
	real(dp), intent(in) :: ProbA       ! Input, Acceptance Probability of area moves
	
	logical :: ex						! Internal, whether file exists
	
	character(len=20) :: iOUT           ! Internal, convert number to string
	character(len=20) :: UOUT           ! Internal, convert number to string
	character(len=20) :: AOUT           ! Internal, convert number to string
	character(len=20) :: DXOUT          ! Internal, convert number to string
	character(len=20) :: DAOUT          ! Internal, convert number to string	
	character(len=20) :: PDOUT          ! Internal, convert number to string
	character(len=20) :: PSOUT          ! Internal, convert number to string
	character(len=20) :: PAOUT          ! Internal, convert number to string

	
	
	
	write(iOUT , '(I10)') iCycle
	write(UOUT , '(F20.5)') Usys
	write(AOUT , '(F10.5)') Area
	write(DXOUT , '(F10.5)') DX
	write(DAOUT , '(F10.5)') DA	
	write(PDOUT, '(F10.5)') ProbD
	write(PSOUT, '(F10.5)') ProbS
	write(PAOUT, '(F10.5)') ProbA

	
	
	
	
	INQUIRE(file = 'out/stats.dat', exist = ex)
	
	if (ex) then
		open(2, file = 'out/stats.dat', access = 'append', status = 'old')
		write(2,'(9a20)') adjustl(iOUT), adjustl(UOUT), adjustl(AOUT), adjustl(DXOUT), adjustl(DAOUT),  &
						  adjustl(PDOUT), adjustl(PSOUT), adjustl(PAOUT)
		close(2)	
	else
		open(2, file = 'out/stats.dat', status = 'new')
		write(2,*) "MCCycle Energy Area DX DA ProbD ProbS ProbA" 
		write(2,'(9a20)') adjustl(iOUT), adjustl(UOUT), adjustl(AOUT), adjustl(DXOUT), adjustl(DAOUT),  &
						  adjustl(PDOUT), adjustl(PSOUT), adjustl(PAOUT)
		close(2)	
	end if	
 end subroutine OUTPUTDAT
 
end module output