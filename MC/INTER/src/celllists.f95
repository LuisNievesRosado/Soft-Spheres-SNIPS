module celllists 
 use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
 implicit none

 real(dp), parameter :: Rcut = 2.25_dp ! Cutoff for YK potential
 
 
 
contains

 subroutine MAKECELL(NPart,SLengthx,Slengthy,HOC,LIS,CELLNEIGH,NCellx,NCelly,X,Y)
	implicit none
	
	integer, intent(in) :: NPart ! Input, Number of particles
	
	real(dp), intent(in) :: X(:)           ! Input, X Coordinates of Particles
	real(dp), intent(in) :: Y(:)           ! Input, Y Coordinates of Particles
	real(dp), intent(in) :: SLengthx       ! Input, Length of Side of Box x
	real(dp), intent(in) :: SLengthy       ! Input, Length of Side of Box y
	
	integer, intent(out) :: HOC(:)          ! Output, Head of Cell Lists
	integer, intent(out) :: LIS(:)          ! Output, Cell Lists
	integer, intent(out) :: CELLNEIGH(:,:)  ! Output, List of Neighbors for each Cell
	
	integer, intent(out) :: NCellx         ! Output, Number of Cells in x Direction
	integer, intent(out) :: NCelly         ! Output, Number of Cells in y Direction
	
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
	NCellx = int(SLengthx/Rcut)
	NCelly = int(SLengthy/Rcut)
	NCell = NCellx*NCelly
	
	
	! Allocate arrays

	allocate(Top(NCellx-2))
	allocate(Bot(NCellx-2))
	allocate(Right(NCelly-2))
	allocate(Left(NCelly-2))

	
	! Initialize HOC
	do iHOC = 1,NCell
		HOC(iHOC) = 0
	end do
	
	! Compute and Save Lists:
	
	do iLIS = 1,NPart
		Cell = CEILING((X(iLIS) + SLengthx/2)*NCellx/SLengthx) + NCellx*(CEILING((Y(iLIS) + SLengthy/2)*NCelly/SLengthy) - 1)	
		LIS(iLIS) = HOC(Cell)	
		HOC(Cell) = iLIS
	end do
	
	! Initialize which cells are on top, bottom, right and left
	
	cBot = 1
	do iBot = 2,NCellx - 1
		Bot(cBot) = iBot
		cBot = cBot + 1
	end do

	cTop = 1
	do iTop = NCell - NCellx + 2,NCell-1
		Top(cTop) = iTop
		cTop = cTop + 1
	end do
	
	cLeft = 1
	do iLeft = NCellx + 1,NCell - 2*NCellx + 1,NCellx
		Left(cLeft) = iLeft
		cLeft = cLeft + 1
	end do
	
	cRight = 1
	do iRight = 2*NCellx,NCell - NCellx,NCellx
		Right(cRight) = iRight
		cRight = cRight + 1
	end do
	

	! Save neighbors of cells
	
	do iNei = 1,NCell
		! If bottom left corner
		if (iNei .eq. 1) then
			CELLNEIGH(iNei,1) = iNei
			CELLNEIGH(iNei,2) = NCellx
			CELLNEIGH(iNei,3) = iNei + 1
			CELLNEIGH(iNei,4) = iNei + NCellx
			CELLNEIGH(iNei,5) = NCell - NCellx + 1
			CELLNEIGH(iNei,6) = iNei + NCellx + 1
			CELLNEIGH(iNei,7) = 2*NCellx
			CELLNEIGH(iNei,8) = NCell - NCellx + 2
			CELLNEIGH(iNei,9) = NCell
		! If bottom right corner
		else if (iNei .eq. NCellx) then
			CELLNEIGH(iNei,1) = iNei
			CELLNEIGH(iNei,2) = iNei - 1
			CELLNEIGH(iNei,3) = 1
			CELLNEIGH(iNei,4) = iNei + NCellx
			CELLNEIGH(iNei,5) = NCellx
			CELLNEIGH(iNei,6) = NCellx + 1
			CELLNEIGH(iNei,7) = iNei + NCellx - 1
			CELLNEIGH(iNei,8) = NCell - NCellx + 1
			CELLNEIGH(iNei,9) = NCell - 1
		! If top left corner
		else if (iNei .eq. NCell - NCellx + 1) then
			CELLNEIGH(iNei,1) = iNei
			CELLNEIGH(iNei,2) = NCell
			CELLNEIGH(iNei,3) = iNei + 1
			CELLNEIGH(iNei,4) = 1
			CELLNEIGH(iNei,5) = iNei - NCellx
			CELLNEIGH(iNei,6) = 2
			CELLNEIGH(iNei,7) = NCellx
			CELLNEIGH(iNei,8) = iNei - NCellx + 1
			CELLNEIGH(iNei,9) = iNei - 1
		! If top right corner
		else if (iNei .eq. NCell) then
			CELLNEIGH(iNei,1) = iNei
			CELLNEIGH(iNei,2) = iNei - 1
			CELLNEIGH(iNei,3) = NCell - NCellx + 1
			CELLNEIGH(iNei,4) = NCellx
			CELLNEIGH(iNei,5) = iNei - NCellx
			CELLNEIGH(iNei,6) = 1
			CELLNEIGH(iNei,7) = NCellx - 1
			CELLNEIGH(iNei,8) = NCell-2*NCellx+1
			CELLNEIGH(iNei,9) = iNei - NCellx - 1
		! If bottom row middle   
		else if (ANY(Bot .eq. iNei) ) then
			CELLNEIGH(iNei,1) = iNei
			CELLNEIGH(iNei,2) = iNei - 1
			CELLNEIGH(iNei,3) = iNei + 1
			CELLNEIGH(iNei,4) = iNei + NCellx
			CELLNEIGH(iNei,5) = NCell - NCellx + iNei
			CELLNEIGH(iNei,6) = iNei + NCellx + 1
			CELLNEIGH(iNei,7) = iNei + NCellx - 1
			CELLNEIGH(iNei,8) = NCell - NCellx + iNei + 1
			CELLNEIGH(iNei,9) = NCell - NCellx + iNei - 1
		! If top row middle
		else if (ANY(Top .eq. iNei) ) then
			CELLNEIGH(iNei,1) = iNei
			CELLNEIGH(iNei,2) = iNei - 1
			CELLNEIGH(iNei,3) = iNei + 1
			CELLNEIGH(iNei,4) = iNei - NCellx*(NCelly - 1)
			CELLNEIGH(iNei,5) = iNei - NCellx
			CELLNEIGH(iNei,6) = iNei - NCellx*(NCelly - 1) + 1
			CELLNEIGH(iNei,7) = iNei - NCellx*(NCelly - 1) - 1
			CELLNEIGH(iNei,8) = iNei - NCellx + 1
			CELLNEIGH(iNei,9) = iNei - NCellx - 1
		! If left column middle
		else if (ANY(Left .eq. iNei) ) then
			CELLNEIGH(iNei,1) = iNei
			CELLNEIGH(iNei,2) = iNei + NCellx -1
			CELLNEIGH(iNei,3) = iNei + 1
			CELLNEIGH(iNei,4) = iNei + NCellx
			CELLNEIGH(iNei,5) = iNei - NCellx
			CELLNEIGH(iNei,6) = iNei + NCellx + 1
			CELLNEIGH(iNei,7) = iNei + 2*NCellx - 1
			CELLNEIGH(iNei,8) = iNei - NCellx + 1
			CELLNEIGH(iNei,9) = iNei - 1
		! If right column middle
		else if (ANY(Right .eq. iNei) ) then
			CELLNEIGH(iNei,1) = iNei
			CELLNEIGH(iNei,2) = iNei - 1
			CELLNEIGH(iNei,3) = iNei - NCellx + 1
			CELLNEIGH(iNei,4) = iNei + NCellx
			CELLNEIGH(iNei,5) = iNei - NCellx
			CELLNEIGH(iNei,6) = iNei + 1
			CELLNEIGH(iNei,7) = iNei + NCellx - 1
			CELLNEIGH(iNei,8) = iNei - 2*NCellx + 1
			CELLNEIGH(iNei,9) = iNei - NCellx - 1
		! If in center
		else
			CELLNEIGH(iNei,1) = iNei
			CELLNEIGH(iNei,2) = iNei - 1
			CELLNEIGH(iNei,3) = iNei + 1
			CELLNEIGH(iNei,4) = iNei + NCellx
			CELLNEIGH(iNei,5) = iNei - NCellx
			CELLNEIGH(iNei,6) = iNei + NCellx + 1
			CELLNEIGH(iNei,7) = iNei + NCellx - 1
			CELLNEIGH(iNei,8) = iNei - NCellx + 1
			CELLNEIGH(iNei,9) = iNei - NCellx - 1
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
 
 
end module celllists