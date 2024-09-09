module celllists 
 use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
 implicit none

 real(dp), parameter :: Rcut = 2.25_dp ! Cutoff for YK potential
 
 
 
contains

 subroutine MAKECELL(NPart,SLength,HOC,LIS,CELLNEIGH,NCell1,X,Y)
	implicit none
	
	integer, intent(in) :: NPart ! Input, Number of particles
	
	real(dp), intent(in) :: X(:)           ! Input, X Coordinates of Particles
	real(dp), intent(in) :: Y(:)           ! Input, Y Coordinates of Particles
	real(dp), intent(in) :: SLength        ! Input, Length of Side of Box
	
	integer, intent(out) :: HOC(:)          ! Output, Head of Cell Lists
	integer, intent(out) :: LIS(:)          ! Output, Cell Lists
	integer, intent(out) :: CELLNEIGH(:,:)  ! Output, List of Neighbors for each Cell
	
	integer, intent(out) :: NCell1         ! Output, Number of Cells in Each Direction
	
	integer :: NCell 					   ! Internal, Total number of cells
	integer :: Cell                        ! Internal, Cell to which atom belongs to
	
	integer :: iHOC						   ! Internal, Loop variable for HOC initialization
	integer :: iLIS   					   ! Internal, Loop variable for LIS initialization
	integer :: iNei						   ! Internal, Loop variable for CELLNEIGH initialization

	integer :: iTop                        ! Internal, Loop variable for CELLNEIGH initialization
	integer :: iBot                        ! Internal, Loop variable for CELLNEIGH initialization
	integer :: iLeft                       ! Internal, Loop variable for CELLNEIGH initialization
	integer :: iRight                      ! Internal, Loop variable for CELLNEIGH initialization
	
	integer :: cTop                        ! Internal, Loop variable for CELLNEIGH initialization
	integer :: cBot                        ! Internal, Loop variable for CELLNEIGH initialization
	integer :: cLeft                       ! Internal, Loop variable for CELLNEIGH initialization
	integer :: cRight                      ! Internal, Loop variable for CELLNEIGH initialization	
	
	integer, allocatable :: Bot(:) 	       ! Internal, Cell neighbors calculation variable
	integer, allocatable :: Top(:)         ! Internal, Cell neighbors calculation variable
	integer, allocatable :: Left(:)        ! Internal, Cell neighbors calculation variable
	integer, allocatable :: Right(:)       ! Internal, Cell neighbors calculation variable
	
	! Compute number of cells per side
	NCell1 = int(SLength/Rcut)
	NCell = NCell1**2
	
	! Allocate arrays
	
	allocate(Top(NCell1-2))
	allocate(Bot(NCell1-2))
	allocate(Right(NCell1-2))
	allocate(Left(NCell1-2))

	
	! Initialize HOC
	do iHOC = 1,NCell
		HOC(iHOC) = 0
	end do
	
	! Compute and Save Lists:
	
	do iLIS = 1,NPart
		Cell = CEILING((X(iLIS) + SLength/2)*NCell1/SLength) + NCell1*(CEILING((Y(iLIS) + SLength/2)*NCell1/SLength) - 1)	
		LIS(iLIS) = HOC(Cell)	
		HOC(Cell) = iLIS
	end do
	
	! Initialize which cells are on top, bottom, right and left
	
	cBot = 1
	do iBot = 2,NCell1 - 1
		Bot(cBot) = iBot
		cBot = cBot + 1
	end do

	cTop = 1
	do iTop = NCell - NCell1 + 2,NCell-1
		Top(cTop) = iTop
		cTop = cTop + 1
	end do
	
	cLeft = 1
	do iLeft = NCell1 + 1,NCell - 2*NCell1 + 1,NCell1
		Left(cLeft) = iLeft
		cLeft = cLeft + 1
	end do
	
	cRight = 1
	do iRight = 2*NCell1,NCell - NCell1,NCell1
		Right(cRight) = iRight
		cRight = cRight + 1
	end do
	
	
	! Save neighbors of cells
	
	do iNei = 1,NCell
		! If bottom left corner
		if (iNei .eq. 1) then
			CELLNEIGH(iNei,1) = iNei
			CELLNEIGH(iNei,2) = NCell1
			CELLNEIGH(iNei,3) = iNei + 1
			CELLNEIGH(iNei,4) = iNei + NCell1
			CELLNEIGH(iNei,5) = NCell - NCell1 + 1
			CELLNEIGH(iNei,6) = iNei + NCell1 + 1
			CELLNEIGH(iNei,7) = 2*NCell1
			CELLNEIGH(iNei,8) = NCell - NCell1 + 2
			CELLNEIGH(iNei,9) = NCell
		! If bottom right corner
		else if (iNei .eq. NCell1) then
			CELLNEIGH(iNei,1) = iNei
			CELLNEIGH(iNei,2) = iNei - 1
			CELLNEIGH(iNei,3) = 1
			CELLNEIGH(iNei,4) = iNei + NCell1
			CELLNEIGH(iNei,5) = NCell
			CELLNEIGH(iNei,6) = NCell1 + 1
			CELLNEIGH(iNei,7) = iNei + NCell1 - 1
			CELLNEIGH(iNei,8) = NCell - NCell1 + 1
			CELLNEIGH(iNei,9) = NCell - 1
		! If top left corner
		else if (iNei .eq. NCell - NCell1 + 1) then
			CELLNEIGH(iNei,1) = iNei
			CELLNEIGH(iNei,2) = NCell
			CELLNEIGH(iNei,3) = iNei + 1
			CELLNEIGH(iNei,4) = 1
			CELLNEIGH(iNei,5) = iNei - NCell1
			CELLNEIGH(iNei,6) = 2
			CELLNEIGH(iNei,7) = NCell1
			CELLNEIGH(iNei,8) = iNei - NCell1 + 1
			CELLNEIGH(iNei,9) = iNei - 1
		! If top right corner
		else if (iNei .eq. NCell) then
			CELLNEIGH(iNei,1) = iNei
			CELLNEIGH(iNei,2) = iNei - 1
			CELLNEIGH(iNei,3) = NCell - NCell1 + 1
			CELLNEIGH(iNei,4) = NCell1
			CELLNEIGH(iNei,5) = iNei - NCell1
			CELLNEIGH(iNei,6) = 1
			CELLNEIGH(iNei,7) = NCell1 - 1
			CELLNEIGH(iNei,8) = NCell-2*NCell1+1
			CELLNEIGH(iNei,9) = iNei - NCell1 - 1
		! If bottom row middle   
		else if (ANY(Bot .eq. iNei) ) then
			CELLNEIGH(iNei,1) = iNei
			CELLNEIGH(iNei,2) = iNei - 1
			CELLNEIGH(iNei,3) = iNei + 1
			CELLNEIGH(iNei,4) = iNei + NCell1
			CELLNEIGH(iNei,5) = NCell - NCell1 + iNei
			CELLNEIGH(iNei,6) = iNei + NCell1 + 1
			CELLNEIGH(iNei,7) = iNei + NCell1 - 1
			CELLNEIGH(iNei,8) = NCell - NCell1 + iNei + 1
			CELLNEIGH(iNei,9) = NCell - NCell1 + iNei - 1
		! If top row
		else if (ANY(Top .eq. iNei) ) then
			CELLNEIGH(iNei,1) = iNei
			CELLNEIGH(iNei,2) = iNei - 1
			CELLNEIGH(iNei,3) = iNei + 1
			CELLNEIGH(iNei,4) = iNei - NCell1*(NCell1 - 1)
			CELLNEIGH(iNei,5) = iNei - NCell1
			CELLNEIGH(iNei,6) = iNei - NCell1*(NCell1 - 1) + 1
			CELLNEIGH(iNei,7) = iNei - NCell1*(NCell1 - 1) - 1
			CELLNEIGH(iNei,8) = iNei - NCell1 + 1
			CELLNEIGH(iNei,9) = iNei - NCell1 - 1
		! If left column
		else if (ANY(Left .eq. iNei) ) then
			CELLNEIGH(iNei,1) = iNei
			CELLNEIGH(iNei,2) = iNei + NCell1 -1
			CELLNEIGH(iNei,3) = iNei + 1
			CELLNEIGH(iNei,4) = iNei + NCell1
			CELLNEIGH(iNei,5) = iNei - NCell1
			CELLNEIGH(iNei,6) = iNei + NCell1 + 1
			CELLNEIGH(iNei,7) = iNei + 2*NCell1 - 1
			CELLNEIGH(iNei,8) = iNei - NCell1 + 1
			CELLNEIGH(iNei,9) = iNei - 1
		! If right column
		else if (ANY(Right .eq. iNei) ) then
			CELLNEIGH(iNei,1) = iNei
			CELLNEIGH(iNei,2) = iNei - 1
			CELLNEIGH(iNei,3) = iNei - NCell1 + 1
			CELLNEIGH(iNei,4) = iNei + NCell1
			CELLNEIGH(iNei,5) = iNei - NCell1
			CELLNEIGH(iNei,6) = iNei + 1
			CELLNEIGH(iNei,7) = iNei + NCell1 - 1
			CELLNEIGH(iNei,8) = iNei - 2*NCell1 + 1
			CELLNEIGH(iNei,9) = iNei - NCell1 - 1
		! If in center
		else
			CELLNEIGH(iNei,1) = iNei
			CELLNEIGH(iNei,2) = iNei - 1
			CELLNEIGH(iNei,3) = iNei + 1
			CELLNEIGH(iNei,4) = iNei + NCell1
			CELLNEIGH(iNei,5) = iNei - NCell1
			CELLNEIGH(iNei,6) = iNei + NCell1 + 1
			CELLNEIGH(iNei,7) = iNei + NCell1 - 1
			CELLNEIGH(iNei,8) = iNei - NCell1 + 1
			CELLNEIGH(iNei,9) = iNei - NCell1 - 1
		end if
	end do
 end subroutine MAKECELL
 
 ! Update cell list after displacement move if needed
 subroutine UPDATECELLDISP(Part,CellOld,CellNew,HOC,LIS)
	implicit none
	
	integer, intent(in) :: Part   		! Input, Particle ID
								  
	integer, intent(in) :: CellOld	  	! Input, Cell number old
	integer, intent(in) :: CellNew  	! Input, Cell number new
	
	integer, intent(inout) :: HOC(:) 	! Input-Output, Head of Cells to be updated
	integer, intent(inout) :: LIS(:) 	! Input-Output, Lists to be updated
	
	integer :: PartPrev           		! Internal, To track switch between lists
	integer :: PartNex            		! Internal, To track switch between lists
	
	integer :: iLIS               		! Internal, Loop variable over LIS
	integer :: iHOC				  		! Internal, Loop variable over HOC
		
	
	! If Part is HOC, change so that LIS(Part) is new HOC then make Part new HOC of CellNew
	if (ANY(HOC .eq. Part) ) then
		PartPrev = 0  ! Not used
		do iHOC = 1,size(HOC)
			if (HOC(iHOC) .eq. Part) then
				HOC(iHOC) = LIS(Part)
				exit
			end if		
		end do		

	else
	! If Part is not HOC, check which particle points to part
		do iLIS = 1,size(LIS)
			if (LIS(iLIS) .eq. Part) then
				PartPrev = iLIS
				exit
			end if
		end do	
		PartNex = LIS(Part)     ! Check which particle Part points to
		LIS(PartPrev) = PartNex ! Make PartPrev point to PartNex
	end if
	
	LIS(Part) = HOC(CellNew) ! Make Part point to current HOC of CellNew 
	HOC(CellNew) = Part      ! Set HOC of CellNew to Part	
 end subroutine UPDATECELLDISP
 
 ! Update cell list after swap move - commented out, no longer in use
! subroutine UPDATECELLSWAP(Part1,Part2,CellN1,CellN2,HOC,LIS)
!	implicit none
!	
!	integer, intent(in) :: Part1            ! Input, ID of Particle 1 swapped
!	integer, intent(in) :: Part2			! Input, ID of Particle 2 swapped
!
!	integer, intent(in) :: CellN1			! Input, Particle 1 New Cell. Old Cell will be CellN2
!	integer, intent(in) :: CellN2			! Input, Particle 2 New Cell. Old Cell will be CellN1
!	
!	integer, intent(inout) :: HOC(:)		! Input-Output, HOC of Cell lists
!	integer, intent(inout) :: LIS(:)		! Input-Output, Cell lists
!	
!	integer :: TempPrev1                    ! Internal, For storing previous/next particle in list before swap
!	integer :: TempNext1                    ! Internal, For storing previous/next particle in list before swap
!
!	integer :: TempPrev2                    ! Internal, For storing previous/next particle in list before swap
!	integer :: TempNext2                    ! Internal, For storing previous/next particle in list before swap
!	
!	integer :: iLIS                         ! Internal, Loop variable over list
!	
!	! If Part 1 was HOC of CellN2 and Part2 was HOC of CellN1
!	if (HOC(CellN2) .eq. Part1 .and. HOC(CellN1) .eq. Part2) then
!	! We need to swap HOC and which particles they point to
!		TempPrev1 = 0 ! Not used
!		TempPrev2 = 0 ! Not used 
!		
!		TempNext1 = LIS(Part1) 
!		TempNext2 = LIS(Part2)
!		
!		! Swap in cell lists
!		LIS(Part1) = TempNext2
!		LIS(Part2) = TempNext1
!		 
!		! Swap in HOC
!		HOC(CellN1) = Part1
!		HOC(CellN2) = Part2 
!	
!
!	! If Part1 was HOC of CellN2 and Part2 was not HOC of CellN1
!	else if (HOC(CellN2) .eq. Part1 .and. HOC(CellN1) .ne. Part2) then
!		TempPrev1 = 0 ! Not used
!		
!		! Find TempPrev of Part2
!		
!		do iLIS = 1,size(LIS)
!			if (LIS(iLIS) .eq. Part2) then
!				TempPrev2 = iLis
!				exit
!			end if
!		end do
!		
!		! Find next in lists
!		TempNext1 = LIS(Part1)
!		TempNext2 = LIS(Part2)
!		
!		! Update Lists
!		LIS(TempPrev2) = Part1
!		
!		LIS(Part1) = TempNext2
!		LIS(Part2) = TempNext1
!		
!		! Update HOC
!		HOC(CellN2) = Part2		
!		
!	! If Part 1 was not HOC of CellN2 and Part2 was HOC of CellN1
!	
!	else if (HOC(CellN2) .ne. Part1 .and. HOC(CellN1) .eq. Part2) then
!		TempPrev2 = 0 ! Not Used
!		
!		! Find TempPrev of Part1
!	
!		do iLIS = 1,size(LIS)
!			if (LIS(iLIS) .eq. Part1) then
!				TempPrev1 = iLis
!				exit
!			end if
!		end do
!		
!		! Find next in lists
!		TempNext1 = LIS(Part1)
!		TempNext2 = LIS(Part2)
!		
!		! Update Lists
!		LIS(TempPrev1) = Part2
!		
!		LIS(Part1) = TempNext2
!		LIS(Part2) = TempNext1
!		
!		! Update HOC
!		HOC(CellN1) = Part1
!	
!	! If Part1 was not HOC of CellN2 and Part2 was not HOC of CellN1
!	else if (HOC(CellN2) .ne. Part1 .and. HOC(CellN1) .ne. Part2) then
!		! Find TempPrevs
!	
!		do iLIS = 1,size(LIS)
!			if (LIS(iLIS) .eq. Part1) then
!				TempPrev1 = iLis
!				exit
!			end if
!		end do
!	
!		do iLIS = 1,size(LIS)
!			if (LIS(iLIS) .eq. Part2) then
!				TempPrev2 = iLis
!				exit
!			end if
!		end do
!		
!		! Find next in lists
!		TempNext1 = LIS(Part1)
!		TempNext2 = LIS(Part2)
!		
!		! Update Lists
!		LIS(TempPrev1) = Part2
!		LIS(TempPrev2) = Part1 
!		
!		LIS(Part1) = TempNext2
!		LIS(Part2) = TempNext1
!	else
!		print *, 'Error found when updating cell lists after swap.'
! 	end if
! end subroutine UPDATECELLSWAP
end module celllists