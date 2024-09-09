module verletlists
 use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
 
 implicit none
 
 contains
 
 ! Make verlet list based on cell lists
 subroutine MAKELIST(VList,NVList,NVMax,VCut,HOC,LIS,CELLNEIGH,NCellx,NCelly,NPart,SLengthx,Slengthy,X,Y)
	implicit none
	integer, intent(inout) :: VList(:,:)   ! Input-Output Verlet Lists
	integer, intent(inout) :: NVList(:)    ! Input-Output Number of Neighbors per Particle
	
	integer, intent(in) :: NVMax           ! Input, Maximum number of neighbors per particle
	real(dp), intent(in) :: VCut           ! Input, Verlet List Cutoff
	
	integer, intent(in) :: NPart           ! Input, Number of particles
	
	integer, intent(in) :: NCellx 		   ! Input, Number of Cells per direction x
	integer, intent(in) :: NCelly 		   ! Input, Number of Cells per direction y
	integer, intent(in) :: HOC(:)          ! Input, HOC of Cell Lists
	integer, intent(in) :: LIS(:)          ! Input, Cell Lists
	integer, intent(in) :: CELLNEIGH(:,:)  ! Input, Neighbor list of Cells
	
	real(dp), intent(in) :: SLengthx       ! Input, Side Length of Box x
	real(dp), intent(in) :: SLengthy       ! Input, Side Length of Box y	
	real(dp), intent(in) :: X(:)           ! Input, X Positions of Particles
	real(dp), intent(in) :: Y(:)           ! Input, Y Positions of Particles
	
	real(dp) :: VCut2					   ! Internal, VCut squared for convenience
	
	integer :: Cell						   ! Internal, Cell of particle
	
	integer :: iPart					   ! Internal, Loop variable, 1 to NPart
	integer :: iCell					   ! Internal, Loop variable, 1 to 9
	integer :: jCell					   ! Internal, Loop variable, Cell of Neighbor under consideration	
	integer :: jPart					   ! Internal, Loop variable, Neighbor under consideration	
	
	real(dp) :: XD  					   ! Internal, Distance in X Direction
	real(dp) :: YD  					   ! Internal, Distance in Y Direction	
	real(dp) :: R2						   ! Internal, Distance squared
	
	! Initialize NVList:
	NVList(:) = 0
	VList(:,:) = 0
	
	VCut2 = VCut**2
	part_loop: do iPart = 1,NPart
		! Find Cell of Particle
		Cell = CEILING((X(iPart) + SLengthx/2)*NCellx/SLengthx) + NCellx*(CEILING((Y(iPart) + SLengthy/2)*NCelly/SLengthy) - 1)
		cell_loop: do iCell = 1,9
			! Choose neighboring cell
			jCell = CELLNEIGH(Cell,iCell)
			jPart = HOC(jCell)
			neigh_loop: do while (jPart .ne. 0)
			    ! Check that particle is not the same and is not already in verlet list
				if (iPart .ne. jPart .and. .not. ANY(VList(iPart,:) .eq. jPart)) then 
					XD = X(iPart)-X(jPart) - SLengthx*NInt((X(iPart)-X(jPart))/SLengthx)  ! Distance adjusting for PBC
					YD = Y(iPart)-Y(jPart) - SLengthy*NInt((Y(iPart)-Y(jPart))/SLengthy)  ! Distance adjusting for PBC
				
					R2 = XD**2 + YD**2 ! Compute distance squared
					if (R2 .lt. VCut2) then
						NVList(iPart) = NVList(iPart) + 1
						NVList(jPart) = NVList(jPart) + 1
						VList(iPart,NVList(iPart)) = jPart
						VList(jPart,NVList(jPart)) = iPart
					end if					
				end if
				jPart = LIS(jPart)
			end do neigh_loop
		end do cell_loop		
	end do part_loop	
 end subroutine MAKELIST

 ! After displacement, if needed update verlet lists of particle that moved too much and neighbors
 subroutine UPTLISTDISP(Part,VList,NVList,VCut,SLengthx,Slengthy,NCellx,Ncelly,HOC,LIS,CELLNEIGH,X,Y)
	implicit none
	integer, intent(in) :: Part					! Input, Particle to change

	integer, intent(inout) :: VList(:,:)		! Input-Output, Verlet Lists
	integer, intent(inout) :: NVList(:)			! Input-Output, Number of Neighbors for Verlet Lists
	
	real(dp), intent(in) :: X(:)				! Input, X Position of Particles after update
	real(dp), intent(in) :: Y(:)				! Input, Y Position of Particles after update
	
	real(dp), intent(in) :: VCut                ! Input, Cutoff for verlet list
	
	real(dp), intent(in) :: SLengthx			! Input, Length of Box Size x
	real(dp), intent(in) :: SLengthy			! Input, Length of Box Size y	
	integer, intent(in) :: NCellx				! Input, Number of Cells per Side x
	integer, intent(in) :: NCelly				! Input, Number of Cells per Side y	
		
	integer, intent(in) :: HOC(:)				! Input, HOC of Cell Lists
	integer, intent(in) :: LIS(:)               ! Input, Cell Lists
	integer, intent(in) :: CELLNEIGH(:,:)		! Input, Neighbors of Cells
		
	integer :: Cell								! Internal, Cell of Particle
	integer :: jPart							! Internal, loop over neighbors
	integer :: jCell							! Internal, loop over cell neighbors	
	integer :: iOld								! Internal, loop over old neighbors
	integer :: jOld								! Internal, loop over neighbors of no longer neighbor
	
	real(dp) :: XD								! Internal, X Distance
	real(dp) :: YD								! Internal, Y Distance
	real(dp) :: R2								! Internal, Squared distance
	real(dp) :: VCut2						    ! Internal, VCut**2 for convenience
	
	integer, allocatable :: ListOld(:)			! Internal, Initial neighbor list of Part
	
	VCut2 = VCut**2
	! Save old neighborlist for checking
	allocate(ListOld(size(VList(Part,:))))
	ListOld = VList(Part,:)
	
	! Update verlet list of part and new neighbors
	Cell = CEILING((X(Part) + SLengthx/2)*NCellx/SLengthx) + NCellx*(CEILING((Y(Part) + SLengthy/2)*NCelly/SLengthy) - 1)
	NVList(Part) = 0
	VList(Part,:) = 0
	cell_loop: do jCell = 1,9
		jPart = HOC(CELLNEIGH(Cell,jCell))
		neigh_loop : do while (jPart .ne. 0)
			if (Part .ne. jPart) then
				XD = X(Part)-X(jPart) - SLengthx*NInt((X(Part)-X(jPart))/SLengthx)  ! Distance adjusting for PBC
				YD = Y(Part)-Y(jPart) - SLengthy*NInt((Y(Part)-Y(jPart))/SLengthy)  ! Distance adjusting for PBC
				
				R2 = XD**2 + YD**2
				
				if (R2 .le. VCut2) then 
					NVList(Part) = NVList(Part) + 1
					VList(Part,NVList(Part)) = jPart
					
					! Check if Part is already neighbor of jPart, if not add
					if (.not. ANY(VList(jPart,:) .eq. Part) ) then
						NVList(jPart) = NVList(jPart) + 1
						VList(jPart,NVList(jPart)) = Part
					end if
				end if
			end if
			jPart = LIS(jPart)
		end do neigh_loop
	end do cell_loop
	
	! Remove Part from neighborlists it doesn't belong to anymore
	old_loop: do iOld = 1,size(ListOld)
		if (.not. ANY(VList(Part,:) .eq. ListOld(iOld))) then
			NVList(ListOld(iOld)) = NVList(ListOld(iOld)) - 1  ! Keeps track, but they will no longer be the first NVList(i)
			oldNeigh_loop: do jOld = 1, size(VList(ListOld(iOld),:))
				if (VList(ListOld(iOld),jOld) .eq. Part) then
					VList(ListOld(iOld),jOld) = 0
					exit oldNeigh_loop
				end if				
			end do oldNeigh_loop
		end if
	end do old_loop
 end subroutine UPTLISTDISP 
 
 ! Update verlet list after trying a swap move, for energy calculation , commented out no longer in use
 !subroutine UPTLISTSWAPTEMP(Part1,Part2,VList1,VList2)
!	implicit none
!	integer, intent(in) :: Part1          			 ! Input, Particle 1 ID
!	integer, intent(in) :: Part2          			 ! Input, Particle 2 ID
!	
!	integer, intent(inout) :: VList1(:)   			 ! Input, Initial Verlet List Particle 1
!	integer, intent(inout) :: VList2(:)   			 ! Input, Initial Verlet List Particle 2	
!	
!	integer, allocatable :: TList1(size(VList1))     ! Internal, Temporary storage for verlet list 1
!	integer, allocatable :: TList2(size(VList2))     ! Internal, Temporary storage for verlet list 2
!	
!	integer :: iList 					  			 ! Internal, Loop variable for list 
!	
!	TList1 = VList1
!	TList2 = VList2
 !
!	! If they were neighbors, exchange from list 
!	do iList = 1,size(VList1)
!		if (TList1(iList) .eq. Part2) then
!			TList1(iList) = Part1
!		end if
!		if (TList2(iList) .eq. Part1) then
!			TList2(iList) = Part2
!		end if
!	end do
!	
!	VList1 = TList2
!	VList2 = TList1	
!	
 !end subroutine UPTLISTSWAPTEMP
 
 ! Update verlet list after a succesful swap move, commented out no longer in use
! subroutine UPTLISTSWAPREAL(VList,NVList,Part1,Part2,XO,YO)
!	implicit none
!	integer, intent(inout) :: VList(:,:)	 ! Input-Output, Verlet Lists
!	integer, intent(inout) :: NVList(:)	     ! Input-Output, Number of Neighbors per Particle
!	
!	integer, intent(in) :: Part1             ! Input, Particle 1 of Swap Move
!	integer, intent(in) :: Part2             ! Input, Particle 2 of Swap Move
!	
!	real(dp), intent(inout) :: XO 			 ! Input-Output, Initial Positions X
!	real(dp), intent(inout) :: YO			 ! Input-Output, Initial Positions Y
!	
!	integer :: TList1(size(VList(Part1,:)))  ! Internal, Temporary storage for verlet list 1
!	integer :: TList2(size(VList(Part2,:)))  ! Internal, Temporary storage for verlet list 2
!
!	integer :: TNList1						 ! Internal, Temporary storage for size of verlet list 1
!	integer :: TNList2					     ! Internal, Temporary storage for size of verlet list 2
!	
!	integer :: iList						 ! Internal, Loop variable for list 
!	integer :: iOther1						 ! Internal, loop over other neighbors of 1 after switch
!	integer :: iOther2						 ! Internal, loop over other neighbors of 2 after switch
!	integer :: jOther1						 ! Internal, loop over other neighbors of 1 after switch
!	integer :: jOther2						 ! Internal, loop over other neighbors of 2 after switch
!
!	
!	real(dp) :: TXO1                         ! Internal, Temporary storage for initial positions
!	real(dp) :: TYO1                         ! Internal, Temporary storage for initial positions
!	real(dp) :: TXO2                         ! Internal, Temporary storage for initial positions
!	real(dp) :: TYO2                         ! Internal, Temporary storage for initial positions
!	
!	
!	TNList1 = NVList(Part1)
!	TNList2 = NVList(Part2)
!	
!	NVList(Part1) = TNList2
!	NVList(Part2) = TNList1
!	
!	TList1 = VList(Part1,:)
!	TList2 = VList(Part2,:)
!	do iList = 1,size(TList1)
!		if (TList1(iList) .eq. Part2) then
!			TList1(iList) = Part1
!		end if
!		if (TList2(iList) .eq. Part1) then
!			TList2(iList) = Part2
!		end if
!	end do
!	
!	VList(Part1,:) = TList2
!	VList(Part2,:) = TList1
!	
!	TXO1 = XO(Part1)
!	TYO1 = YO(Part1)
!	
!	TXO2 = XO(Part2)
!	TYO2 = YO(Part2)
!		
!	XO(Part1) = TXO2	
!	YO(Part1) = TYO2
!	
!	XO(Part2) = TXO1	
!	YO(Part2) = TYO1	
!	
!	! Update other neighbor lists
!	
!	other_loop1: do iOther1 = 1,size(VList(Part1,:))
!		if (ANY(VList(VList(Part1,iOther1),:) .eq. Part1) .and. ANY(VList(VList(Part1,iOther1),:) .eq. Part2)) then
!			exit other_loop1! Neighbor of both Part 1 and Part 2, leave alone
!		else
!			otherneigh_loop1: do jOther1 = 1,size(VList(VList(Part1,iOther1),:))
!				if (VList(VList(Part1,iOther1),jOther1) .eq. Part2) then
!					VList(VList(Part1,iOther1),:jOther1 = Part1
!					exit otherneigh_loop1
!				end if
!			end do otherneigh_loop1
!		end if	
!	end do other_loop1
!	
!	other_loop2: do iOther2 = 1,size(VList(Part2,:))
!		if (ANY(VList(VList(Part2,iOther1),:) .eq. Part1) .and. ANY(VList(VList(Part2,iOther1),:) .eq. Part2)) then
!			exit other_loop2! Neighbor of both Part 1 and Part 2, leave alone
!		else
!			otherneigh_loop2: do jOther2 = 1,size(VList(VList(Part2,iOther2),:))
!				if (VList(VList(Part2,iOther2),jOther2) .eq. Part1) then
!					VList(VList(Part2,iOther2),:jOther2 = Part2
!					exit otherneigh_loop2
!				end if
!			end do otherneigh_loop2
!		end if	
!	end do other_loop2
! end subroutine UPTLISTSWAPREAL
 
end module verletlists