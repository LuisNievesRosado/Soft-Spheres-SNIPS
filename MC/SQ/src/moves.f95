module moves 
 use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
 use energy
 use celllists
 use verletlists
 
 implicit none

 contains

 subroutine DISPLACEMOVE(AceptD,Beta,DX,NPart,NVList,VList,VCut,NCell1,HOC,LIS,CELLNEIGH,SLength,TYP,X,Y,XO,YO)
	implicit none
	
	integer, intent(in) :: NCell1           ! Input, Number of Cells in each direction
	real(dp), intent(in) :: SLength			! Input, Side length of box
	
	real(dp), intent(in) :: Beta			! Input, 1/Temp
	
	real(dp), intent(in) :: DX				! Input, Displacement extent
	
	integer, intent(in) :: NPart			! Input, Total Number of Particles
	integer, intent(in) :: TYP(:)	        ! Input, Types of Particles

	integer, intent(in) :: CELLNEIGH(:,:)   ! Input, Neighbors of Cells
	
	real(dp), intent(in) :: VCut			! Input, Verlet list cut
	
	integer, intent(inout) :: HOC(:)		! Input-Output, HOC of Cell Lists
	integer, intent(inout) :: LIS(:)		! Input-Output, Cell Lists

	
	integer, intent(inout) :: VList(:,:)    ! Input-Output, Verlet List
	integer, intent(inout) :: NVList(:)     ! Input-Output, Verlet List
	
	real(dp), intent(inout)	:: X(:)    	    ! Input-Output, X Positions
	real(dp), intent(inout)	:: Y(:)    		! Input-Output, Y Positions

	real(dp), intent(inout)	:: XO(:)   		! Input-Output, Initial X Positions
	real(dp), intent(inout)	:: YO(:)   		! Input-Output, Initial Y Positions		

	integer, intent(out) :: AceptD			! Output, whether move was accepted or rejected
	
	integer :: CellI						! Internal, Cell of particle before move
	integer :: CellF						! Internal, Cell of particle after move
	
	real(dp) :: randC						! Internal, Random number for acceptance check
	real(dp) :: randN    					! Internal, Random number to choose particle
	real(dp) :: randX 						! Internal, Random number to move in X
	real(dp) :: randY 						! Internal, Random number to move in Y
	real(dp) ::	UPre						! Internal, Initial energy of particle
	real(dp) :: UPos						! Internal, post-move energy of particle

	real(dp) :: Check						! Internal, result of check
	
	real(dp) :: DO2							! Internal, distance from X0, Y0, squared
	real(dp) :: MaxD						! Internal, max distance before update
	
	real(dp) :: XProp                       ! Internal, proposed X location post displacement
	real(dp) :: YProp                       ! Internal, proposed Y location post displacement
	
	integer, allocatable :: VListP(:)		! Internal, verlet list of particle
	integer :: Part      					! Internal, Particle to Displace
	
	
	
	! Choose a particle at random
	call RANDOM_NUMBER(randN)
	
	Part = CEILING(NPart*randN) ! Ceiling guarantees 1 to NPart, equally
	
	allocate(VListP(size(VList(Part,:))))
	VListP = VList(Part,:)
	
	XProp = X(Part)
	YProp = Y(Part)
	
	CellI = CEILING((XProp + SLength/2)*NCell1/SLength) + NCell1*(CEILING((YProp + SLength/2)*NCell1/SLength) - 1)
	
	! Compute initial particle energy
	call ENERGYPART(UPre,Part,SLength,VListP,TYP,X,Y,XProp,YProp)
	
	! Propose move
	call RANDOM_NUMBER(randX)
	call RANDOM_NUMBER(randY)
	
	XProp = X(Part) + DX*(randX - 0.5_dp)
	YProp = Y(Part) + DX*(randY - 0.5_dp)
	
	if (XProp .gt. SLength/2) then
		XProp  = XProp - SLength
	else if (XProp .lt. -SLength/2) then
		XProp = XProp + SLength
	end if
	
	if (YProp .gt. SLength/2) then
		YProp  = YProp - SLength
	else if (YProp .lt. -SLength/2) then
		YProp = YProp + SLength
	end if	
	
	CellF = CEILING((XProp + SLength/2)*NCell1/SLength) + NCell1*(CEILING((YProp + SLength/2)*NCell1/SLength) - 1)
	
	! Compute new particle energy

	call ENERGYPART(UPos,Part,SLength,VListP,TYP,X,Y,XProp,YProp)
	 

	! Check if accepted
	Check = exp(-Beta*(UPos-UPre))	
	
	call RANDOM_NUMBER(randC)
	

	if (Check .ge. 1.0_dp) then
		AceptD = 1
	else if (randC .lt. Check) then
		AceptD = 1
	else
		AceptD = 0
	end if
	
	
	

	! If accepted, update everything
	if (AceptD .eq. 1) then
		! Update positions
		X(Part) = XProp
		Y(Part) = YProp
		
		! Check if cell changed, if it did, update cell lists
		
		if (CellI .ne. CellF) then
			call UPDATECELLDISP(Part,CellI,CellF,HOC,LIS)
		end if
	
		! Check if particle has moved too much, if it did, update verlet lists
		
		MaxD = ((VCut - Rcut)/2.0_dp)**2 
		DO2 = (XProp-XO(Part))**2 + (YProp-YO(Part))**2
		
		if (DO2 .ge. MaxD) then 
			call UPTLISTDISP(Part,VList,NVList,VCut,SLength,NCell1,HOC,LIS,CELLNEIGH,X,Y)
			! Reset intial position vector
			XO(Part) = XProp
			YO(Part) = YProp
		end if
		
	end if	
 end subroutine DISPLACEMOVE

 subroutine SWAPMOVE(AceptS,Beta,NPart,NPartO,NPartI,NPartS, &
					NVList,VList,NCell1,HOC,LIS,CELLNEIGH,SLength,TYP,X,Y,XO,YO)
	implicit none

	integer, intent(in) :: NPart			! Input, Total Number of Particles
	integer, intent(in) :: NPartO			! Input, Total Number of ISVO
	integer, intent(in) :: NPartI           ! Input, Total Number of ISV
	integer, intent(in) :: NPartS 			! Input, Total Number of SVPS

	integer, intent(in) :: NCell1			! Input, number of cells per side
	
	integer, intent(in) :: CELLNEIGH(:,:)   ! Input, cell neighbors
	
	real(dp), intent(in) :: Beta			! Input, 1/Temp
	real(dp), intent(in) :: SLength			! Input, length of Side
	
	integer, intent(inout) :: TYP(:)		! Input-Output, Type of particles
	
	integer, intent(inout) :: NVList(:)		! Input-Output, Number of Neighbors
	integer, intent(inout) :: VList(:,:)	! Input-Output, Verlet Lists
	
	integer, intent(inout) :: HOC(:)		! Input-Output, HOC of Cell Lists
	integer, intent(inout) :: LIS(:)		! Input-Output, Cell Lists
	
	real(dp), intent(inout)	:: X(:)    	    ! Input-Output, X Positions
	real(dp), intent(inout)	:: Y(:)    		! Input-Output, Y Positions

	real(dp), intent(inout)	:: XO(:)   		! Input-Output, Initial X Positions
	real(dp), intent(inout)	:: YO(:)   		! Input-Output, Initial Y Positions
	
	integer, intent(out) :: AceptS			! Output, whether move was accepted
	
	integer, allocatable :: TYPT(:)			! Internal, temporary storage for types
	
	integer :: TYP1							! Internal, type of part1 before swap
	integer :: TYP2							! Internal, type of part2 before swap
	
	integer :: Part1						! Internal, Particle 1 to Swap
	integer :: Part2						! Internal, Particle 2 to Swap
	
	!integer :: CellN1						! Internal, Cell of Part 1 Post-Move
	!integer :: CellN2						! Internal, Cell of Part 2 Post-Move
	integer, allocatable :: VList1(:)		! Internal, Verlet list of particle 1 
	integer, allocatable :: VList2(:)		! Internal, Verlet list of particle 2 
	
	integer, allocatable :: Part2L(:)		! Internal, Possible swaps after Part1 has been chosen
	
	real(dp) :: UPart1I						! Internal, Initial Energy Part 1
	real(dp) :: UPart2I						! Internal, Initial Energy Part 2
	
	real(dp) :: UPart1F						! Internal, Final Energy Part 1
	real(dp) :: UPart2F						! Internal, Final Energy Part 2
	
	real(dp) :: randN1						! Internal, random number for particle selection
	real(dp) :: randN2						! Internal, random number for particle selection
    real(dp) :: randT						! Internal, random number for particle selection
	
	real(dp) :: randC						! Internal, random number for acceptance selection
	real(dp) :: Check					    ! Internal, acceptance check variable 
	
	integer :: iPart						! Internal, Loop variable
	integer :: cPart						! Internal, Loop variable
	
	! Choose particle at Random
	call RANDOM_NUMBER(randN1)
	call RANDOM_NUMBER(randN2)
	
	Part1 = CEILING(NPart*randN1)

	cPart = 1
	if (TYP(Part1) .eq. 1) then
		allocate(Part2L(NPartI+NPartS))
		do iPart = 1,NPart
			if (TYP(iPart) .ne. 1) then
				Part2L(cPart) = iPart
				cPart = cPart + 1
			end if	
		end do
		Part2 = Part2L(CEILING((NPartI+NPartS)*randN2))
	else if (TYP(Part1) .eq. 2) then
		allocate(Part2L(NPartO+NPartS))	
		do iPart = 1,NPart
			if (TYP(iPart) .ne. 2) then
				Part2L(cPart) = iPart
				cPart = cPart + 1
			end if	
		end do
		Part2 = Part2L(CEILING((NPartO+NPartS)*randN2))
	else	
		allocate(Part2L(NPartO+NPartI))	
		do iPart = 1,NPart
			if (TYP(iPart) .ne. 3) then
				Part2L(cPart) = iPart
				cPart = cPart + 1
			end if	
		end do
		Part2 = Part2L(CEILING((NPartO+NPartI)*randN2))
	end if
		
	
	
	! ================================================================================!
	! Old code, assumed order TYP would always remain the same

	! Choose another particle, with a different type
	!if (TYP(Part1) .eq. 1) then
	!	! Choose between NPartO to NPart
	!	call RANDOM_NUMBER(randN2)
	!	Part2 = NPartO + CEILING((NPart - NPartO)*randN2)
	!else if (TYP(Part1) .eq. 2) then
	!	! Choose between 1 to NPartO or NPartO+NPartI + 1 to NPart
	!	call RANDOM_NUMBER(randN2)		
	!	call RANDOM_NUMBER(randT)
	!	! Choose whether to use TYP 1 or TYP 3 Particle, scaled to populations
	!	if (randT .lt. float(NPartO)/(float(NPartO)+float(NPartI))) then 
	!		! Use TYP 1; 1 to NPartO
	!		Part2 = CEILING(NPartO*randN2)
	!	else 
	!		! Use TYP 3 ; NPartO + NPartI + 1 to NPart
	!		Part2 = NPartO + NPartI + CEILING(NPartS*randN2)
	!	end 
	!else
	!	! Choose between 1 to NPart-NPartS
	!	call RANDOM_NUMBER(randN2)
	!	Part2 = CEILING((NPart - NPartS)*randN2)
	!end if
	
	! ================================================================================!
	
	
	! ================================================================================!
	! Find Positions and Neighbors
	
	allocate(VList1(size(VList(Part1,:))))
	allocate(VList2(size(VList(Part2,:))))
	
	VList1 = VList(Part1,:)
	VList2 = VList(Part2,:)
	
	! Find energies Pre Swap, need to feed X(Partn),Y(Partn) due to using same routine as disp)
	call ENERGYPART(UPart1I,Part1,SLength,VList1,TYP,X,Y,X(Part1),Y(Part1))
	call ENERGYPART(UPart2I,Part2,SLength,VList2,TYP,X,Y,X(Part2),Y(Part2))	
	
	! Propose swap, change types
	allocate(TYPT(size(TYP)))
	
	TYP1 = TYP(Part1)
	TYP2 = TYP(Part2)
	
	TYPT = TYP
	
	TYPT(Part1) = TYP2
	TYPT(Part2) = TYP1
	
	!============================================================================!
	! Commented out old method switching positons, required more work
	!X1F = X2I
	!Y1F = Y2I
	!
	!X2F = X1I
	!Y2F = Y2I
	
	
	! Propose swap, switch verlet lists for energy calculation
	
	!call UPTLISTSWAPTEMP(Part1,Part2,VList1,VList2)
	!============================================================================!
	
	! Find Energies Post Swap

	call ENERGYPART(UPart1F,Part1,SLength,VList1,TYPT,X,Y,X(Part1),Y(Part1))
	call ENERGYPART(UPart2F,Part2,SLength,VList2,TYPT,X,Y,X(Part2),Y(Part2))	
	
	Check = exp(-Beta*(UPart1F + UPart2F - UPart1I - UPart2I))
	
	! Check if move is accepted
	
	call RANDOM_NUMBER(randC)
	if (Check .ge. 1.0_dp) then
		AceptS = 1
	else if (randC .lt. Check) then
		AceptS = 1
	else
		AceptS = 0
	end if
	
	! If Accepted, Update types
	
	if (AceptS .eq. 1) then
		TYP = TYPT
	end if
	
	
	
	! ================================================================================!
	! Old Code, switching positions instead of types, more work that way
	!
	!if (AceptD .eq. 1) then
	!! Update CellList
	!
	!	CellN1 = int((X1F + SLength/2)*NCell1/SLength) + NCell1*int((Y1F + SLength/2)*NCell1/SLength)
	!	CellN2 = int((X2F + SLength/2)*NCell1/SLength) + NCell1*int((Y2F + SLength/2)*NCell1/SLength)
	!
	!	call UPDATECELLSWAP(Part1,Part2,CellN1,CellN2,HOC,LIS)
	!
	!! Update Verlet-list, along with initial positions
	!
	!	call UPTLISTSWAPREAL(VList,NVList,Part1,Part2,XO,YO)
	!
	!! Update Positions
	!
	!	X(Part1) = X1F
	!	Y(Part1) = Y1F
	!
	!	X(Part2) = X2F
	!	Y(Part2) = Y2F
	!	
	!! Update system energy
	!	Usys = Usys + UPart1F + UPart2F - UPart1I - UPart2I
	!end if
	! ================================================================================!
 end subroutine SWAPMOVE
	
 subroutine AREAMOVE(AceptA,Pres,Beta,DA,Area,SLength,NPart, & 
			NCell1,NVMax,VCut,NVList,VList,HOC,LIS,CELLNEIGH, &
			TYP,X,Y,XO,YO,HOCOUT,CELLNEIGHOUT)
	implicit none
	
	integer, intent(in) :: NPart				! Input, Number of Particles
	integer, intent(in) :: TYP(:)				! Input, Types of Particles
	integer, intent(in) :: NVMax				! Input, Max Number of Neighbors per particle
	
	real(dp), intent(in) :: VCut				! Input, Cutoff for Verlet 
	real(dp), intent(in) :: DA					! Input, Extent of area move
	real(dp), intent(in) :: Pres				! Input, Pressure of system
	real(dp), intent(in) :: Beta				! Input, 1/Temp
	
	integer, intent(inout) :: NCell1			! Input-Output, Number of Cells in Each Direction	
	
	integer, intent(inout) :: NVList(:)			! Input-Output, Number of Neighbors per particle
	integer, intent(inout) :: VList(:,:)		! Input-Output, Verlet Lists
	
	integer, intent(inout) :: HOC(:)			! Input-Output, HOC of Cell Lists
	integer, intent(inout) :: LIS(:)			! Input-Output, Cell lists
	integer, intent(inout) :: CELLNEIGH(:,:)	! Input-Output, Cell Neighbors
	
	real(dp), intent(inout) :: Area				! Input-Output, Area of System
	real(dp), intent(inout) :: SLength			! Input-Output, Side Length
	
	real(dp), intent(inout) :: X(:)				! Input-Output, X Positions
	real(dp), intent(inout) :: Y(:)				! Input-Output, Y Positions
	
	real(dp), intent(inout) :: XO(:)				! Input-Output, Initial X Positions
	real(dp), intent(inout) :: YO(:)				! Input-Output, Initial Y Positions	
	
	integer, intent(out) :: AceptA				! Output, whether move is accepted 
	
	integer, allocatable, intent(out) :: HOCOUT(:)			 ! Output, exit HOC, might have changed length
	integer, allocatable, intent(out) :: CELLNEIGHOUT(:,:)   ! Output, exit Cell lists neighbors, might have changed length
	
	
	real(dp) :: lnNA							! Internal, Log of new area
	real(dp) :: NewA							! Internal, new area
	real(dp) :: NewSLength						! Internal, new side length
	real(dp) :: OldU							! Internal, initial system energy	
	real(dp) :: NewU							! Internal, new system energy	
	
	real(dp) :: randA							! Internal, random number for area move
	real(dp) :: randC							! Internal, random number for acceptance
	
	real(dp) :: Check							! Internal, check for acceptance
	
	real(dp), allocatable :: XT(:)              ! Internal, temporary storage for positions
	real(dp), allocatable :: YT(:)              ! Internal, temporary storage for positions
										    
	real(dp), allocatable :: XOT(:)             ! Internal, temporary storage for initial positions
	real(dp), allocatable :: YOT(:)             ! Internal, temporary storage for initial positions
	
	integer :: NCell1T							! Internal, temporary storage for Number of Cells in Each Direction	
	integer, allocatable :: HOCT(:)				! Internal, temporary storage for HOC of Cell Lists
	integer, allocatable :: LIST(:)				! Internal, temporary storage for Cell lists
	integer, allocatable :: CELLNEIGHT(:,:)		! Internal, temporary storage for Cell Neighbors
	
	integer, allocatable :: NVListT(:)			! Internal, temporary storage for Number of Neighbors per particle
	integer, allocatable :: VListT(:,:)		    ! Internal, temporary storage for Verlet Lists
	
	integer :: iPart							! Internal, loop variable over particles

	
	logical :: ex								! Internal, whether file exists
	character(len=20) :: AOUT                   ! Internal, for area move troubleshooting
	character(len=20) :: NAOUT                  ! Internal, for area move troubleshooting
	character(len=20) :: COUT                   ! Internal, for area move troubleshooting
	character(len=20) :: UOUT                   ! Internal, for area move troubleshooting
	character(len=20) :: NUOUT                  ! Internal, for area move troubleshooting	
	
	! Compute initial system energy
	call ENERGYSYS(OldU,NPart,SLength,NVList,VList,TYP,X,Y)

	allocate(XT(NPart))
	allocate(YT(NPart))	
	allocate(XOT(NPart))
	allocate(YOT(NPart))	
	
	! Make area move - logarithmic move
	!call RANDOM_NUMBER(randA)
	
	!lnNA = log(Area) + (randA - 0.5_dp)*DA
	
	!NewA = exp(lnNA)
	!NewSLength = sqrt(NewA)
	
	
	! Make area move - direct move
	call RANDOM_NUMBER(randA)
	
	NewA = Area + (randA-0.5_dp)*DA
	NewSLength = sqrt(NewA)
	
	! Update all Positions
	
	do iPart = 1,NPart
		XT(iPart) = X(iPart)*NewSLength/SLength
		YT(iPart) = Y(iPart)*NewSLength/SLength
		XOT(iPart) = XO(iPart)*NewSLength/SLength
		YOT(iPart) = YO(iPart)*NewSLength/SLength	
	end do
	
	! Update cell Lists
	
	allocate(HOCT(int(NewSLength/Rcut)**2))
	allocate(LIST(NPart))
	allocate(CELLNEIGHT(int(NewSLength/Rcut)**2,9))	

	call MAKECELL(NPart,NewSLength,HOCT,LIST,CELLNEIGHT,NCell1T,XT,YT)

	! Update verlet Lists
	allocate(VListT(NPart,NVMax))
	allocate(NVListT(NPart))
	
	call MAKELIST(VListT,NVListT,NVMax,VCut,HOCT,LIST,CELLNEIGHT,NCell1T,NPart,NewSLength,XT,YT)
	! Compute new system energy

	call ENERGYSYS(NewU,NPart,NewSLength,NVListT,VListT,TYP,XT,YT)

	! Check if move is accepted - logarithmic move
	!Check = exp(-Beta*(NewU - OldU + Pres*(NewA - Area) - (NPart + 1.0_dp)*log(NewA/Area)/Beta))
	
	! Check if move is accepted - direct
	
	Check = exp(-Beta*(NewU - OldU + Pres*(NewA - Area)) + NPart*log(NewA/Area))
	!print *, -Beta*(NewU - OldU + Pres*(NewA - Area)) + NPart*log(NewA/Area), Check 
	call RANDOM_NUMBER(randC)
	if (Check .ge. 1.0_dp) then
		AceptA = 1
	else if (randC .lt. Check) then
		AceptA = 1
	else
		AceptA = 0
	end if
	
	! Output data for troubleshooting
	!INQUIRE(file = 'out/area.dat', exist = ex)
	!
	!write(AOUT , '(F20.5)') Area
	!write(NAOUT, '(F20.5)') NewA
	!write(COUT , '(F20.5)') Check
	!write(UOUT , '(F20.5)') OldU
	!write(NUOUT , '(F20.5)') NewU	
	!
	!if (ex) then
	!	open(2, file = 'out/area.dat', access = 'append', status = 'old')
	!	write(2,'(9a20)') adjustl(AOUT), adjustl(NAOUT), adjustl(COUT), adjustl(UOUT), adjustl(NUOUT)
	!	close(2)	
	!else
	!	open(2, file = 'out/area.dat', status = 'new')
	!	write(2,*) "Area NewArea Check OldU NewU" 
	!	write(2,'(9a20)') adjustl(AOUT), adjustl(NAOUT), adjustl(COUT), adjustl(UOUT), adjustl(NUOUT)
	!	close(2)	
	!end if	
	
	
	
	! If move is accepted, update everything
	if (AceptA .eq. 1) then
	
	! Geometry
		Area = NewA
		SLength = NewSLength

	! Positions, XO and YO too since we rebuilt verlet lists
		X = XT
		Y = YT
		XO = XOT
		YO = YOT
		
	! Cell Lists
	
		allocate(HOCOUT(size(HOCT)))
		allocate(CELLNEIGHOUT(size(HOCT),9))
		
		NCell1 = NCell1T
		HOCOUT = HOCT                ! Changes length, need special treatment outside move
		LIS = LIST
		CELLNEIGHOUT = CELLNEIGHT	 ! Changes length, needs special treatment outside move
		
	! Verlet Lists
	
		VList = VListT
		NVList = NVListT

	else 
	! Output HOCOUT and CELLNEIGHOUT as old ones, not used outside
		allocate(HOCOUT(size(HOC)))
		allocate(CELLNEIGHOUT(size(HOC),9))
		
		HOCOUT = HOC
		CELLNEIGHOUT = CELLNEIGH
	end if
	
 end subroutine AREAMOVE

end module moves