program main
 use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
 use celllists
 use energy 
 use verletlists
 use moves
 use output
 
 implicit none
 
 !===========================================================================!
 ! Variables to be read from input file
 !===========================================================================!
 
 integer :: NCycle ! Number of MC Cycles
 integer :: NTherm ! Number of thermalization cycles
 integer :: NOutD  ! Number of Cycles Before Data Output
 integer :: NOutP  ! Number of Cycles Before Position Output
 integer :: NPartO ! Number of ISVO Particles
 integer :: NPartI ! Number of ISV Particles
 integer :: NPartS ! Number of SPVS Particles
 integer :: NVMax  ! Maximum Number of Verlet Neighbors per Particle
 
 real(dp) :: Rho   ! Initial Surface Density of Particles
 real(dp) :: Temp  ! Desired Temperature of System
 real(dp) :: Pres  ! Desired Pressure of System

 real(dp) :: DX    ! Displacement Move Magnitude
 real(dp) :: DA    ! Area Move Magnitude
 
 real(dp) :: Vcut  ! Cutoff for verlet list
 
 !===========================================================================!
 ! Internal Variables
 !===========================================================================!
 
 integer :: NPart    ! Total Number of Particles
 integer :: NMove    ! Total Number of Moves per MC Cycle

 integer :: AceptD   ! 1 or 0. Tracks Whether Last Displacement Move Was Accepted
 integer :: AceptS   ! 1 or 0. Tracks Whether Last Swap Move Was Accepted
 integer :: AceptA   ! 1 or 0. Tracks Whether Last Area Move Was Accepted

 integer :: NTotD    ! Total number of Displacement Moves Attempted
 integer :: NTotS    ! Total number of Swap Moves Attempted
 integer :: NTotA    ! Total number of Area Moves Attempted

 integer :: ATotD    ! Total number of Displacement Moves Accepted
 integer :: ATotS    ! Total number of Swap Moves Accepted
 integer :: ATotA    ! Total number of Area Moves Accepted
 
 real(dp) :: ProbD   ! Acceptance Probability of Displacement Move 
 real(dp) :: ProbS   ! Acceptance Swap of Displacement Move
 real(dp) :: ProbA   ! Acceptance Area of Displacement Move
 
 
 real(dp) :: Area    ! Area of System
 
 real(dp) :: SLength ! Length of Box
 
 real(dp) :: Beta    ! 1/T
 real(dp) :: Usys    ! Total Energy of System 

 
 real(dp) :: MoveType ! 0-1 Random number, sets type of move
 
 integer, allocatable :: TYP(:) ! Particle Type 
 
 real(dp), allocatable :: X(:) ! Particle X Positions
 real(dp), allocatable :: Y(:) ! Particle Y Positions
 
 real(dp), allocatable :: XO(:) ! Particle X Initial Positions for Verlet List
 real(dp), allocatable :: YO(:) ! Particle Y Initial Positions for Verlet List
 
 integer, allocatable :: HOC(:) ! Head of Cell
 integer, allocatable :: LIS(:) ! Cell List
 
 integer, allocatable :: CELLNEIGH(:,:) ! List of Neighbors for each Cell
 
 integer, allocatable :: HOCOUT(:)         ! For area move update
 integer, allocatable :: CELLNEIGHOUT(:,:) ! For area move update
 
 
 integer :: NCell1   ! Number of cells in each direction
 integer :: IniCycle ! Set to 0 for output
 integer, allocatable :: VList(:,:)   ! Verlet List of neigbors
 integer, allocatable :: NVList(:)    ! Number of neighbors of each particle
 
 real(dp), dimension(14) :: INP ! Input file data storage

  
 ! Loop integers
 integer :: iPart    ! Loop varaible for particle setup, 1: NPart
 integer :: iPartO   ! Loop variable for particle set up 1:NPartO
 integer :: iPartI   ! Loop variable for particle set up NPartO + 1:NPartO + NPartI
 integer :: iPartS   ! Loop variable for particle set up NPartO + NPartI + 1: NPartO + NPartI + NPartS
 integer :: iRead    ! Loop variable for reading input file
 integer :: iCycle   ! Loop variable for outer loop 1:NCycle
 integer :: iMove    ! Loop variable for inner loop 1:NMove
 integer :: OUTCOUNTD! Loop variable for output
 integer :: OUTCOUNTP! Loop variable for output
 integer :: ACOUNT   ! Loop variable for area move adjustment
 integer :: DCOUNT   ! Loop variable for displacement move adjustment
 !===========================================================================!
 ! Read data from file 
 !===========================================================================!
 
 open(1, file='input.dat',status='old')
 print *, 'Reading simulation parameters'
 do iRead = 1,14
	read(1,*) INP(iRead)
 end do	
 close(1)
 
 NCycle = int(INP(1))
 NTherm = int(INP(2))
 NOutD = int(INP(3))
 NOutP = int(INP(4))
 NPartO = int(INP(5))
 NPartI = int(INP(6))
 NPartS = int(INP(7))
 NVMax = int(INP(8))
 
 Rho  = INP(9)
 Temp = INP(10)
 Pres = INP(11)

 DX = INP(12)
 DA = INP(13)
 Vcut = INP(14)
 
 !===========================================================================!
 ! Initial System Setup - Simulation
 !===========================================================================!
 NPart = NPartO + NPartI + NPartS  
 
 ! Based on Prajwal Recomendation
 NMove = NPart + NPart/10 + 2  ! NPart Disp Moves, NPart/10 Swap Moves, 2 Area Moves
 
 Beta = 1.0_dp/Temp  
 !===========================================================================!
 ! Initial System Setup - Box
 !===========================================================================!
 Area = NPart/Rho                  
 SLength = sqrt(Area)
 
 !===========================================================================!
 ! Initial System Setup - Particles
 !===========================================================================!
 
 ! Initialize Types of Particles
 allocate(TYP(NPart))
 
 do iPartO = 1, NPartO
	TYP(iPartO) = 1
 end do

 do iPartI = NPartO + 1, NPartO + NPartI
	TYP(iPartI) = 2
 end do
 
 do iPartS = NPartO + NPartI + 1, NPartO + NPartI + NPartS
	TYP(iPartS) = 3
 end do 
 

 
 ! Initialize Particle Locations
 allocate(X(NPart))
 allocate(Y(NPart)) 
 
 call RANDOM_NUMBER(X)
 call RANDOM_NUMBER(Y)
 
 X = X*SLength - SLength/2
 Y = Y*SLength - SLength/2
 
 ! Save Initial Locations
 allocate(XO(NPart))
 allocate(YO(NPart))

 do iPart = 1,NPart
	XO(iPart) = X(iPart)
	YO(iPart) = Y(iPart)
 end do

 !===========================================================================!
 ! Initial System Setup - Cell Lists
 !===========================================================================!
 
 allocate(HOC(int(SLength/Rcut)**2))
 allocate(LIS(NPart))
 allocate(CELLNEIGH(int(SLength/Rcut)**2,9))
 
 call MAKECELL(NPart,SLength,HOC,LIS,CELLNEIGH,NCell1,X,Y)
 !===========================================================================
 ! Initial System Setup - Neighbor Lists
 !===========================================================================!
 allocate(VList(NPart,NVMax))
 allocate(NVList(NPart))
 
 call MAKELIST(VList,NVList,NVMax,Vcut,HOC,LIS,CELLNEIGH,NCell1,NPart,SLength,X,Y)
 !===========================================================================!
 ! Initial System Setup - Energy
 !===========================================================================!
 call ENERGYSYS(Usys,NPart,SLength,NVLISt,VList,TYP,X,Y)
 
 
 !===========================================================================!
 ! Initial System Setup - Output Initial Files
 !===========================================================================!
 IniCycle = 0
 call OUTPUTPOS(IniCycle,NPart,TYP,X,Y)
 
 !===========================================================================!
 ! Initial System Setup - Initialize Counters
 !===========================================================================!
 
 NTotD = 0
 NTotS = 0
 NTotA = 0
 
 ATotD = 0
 ATotS = 0
 ATotA = 0
 
 AceptD = 0
 AceptS = 0
 AceptA = 0

 OUTCOUNTD = -NTherm
 OUTCOUNTP = -NTherm
 ACOUNT = 0
 DCOUNT = 0
 
 ! Start MC run
 outer_loop: do iCycle = 1,NCycle
 	OUTCOUNTP = OUTCOUNTP + 1
	OUTCOUNTD = OUTCOUNTD + 1
	ACOUNT = ACOUNT + 1
	DCOUNT = DCOUNT + 1
	inner_loop: do iMove = 1,NMove
		! Choose type of move
		call RANDOM_NUMBER(MoveType)
		! If Displacement Move - aim for 40% acceptance
		if (MoveType .le. NPart/float(NMove)) then
			NTotD = NTotD + 1
			call DISPLACEMOVE(AceptD,Beta,DX,NPart,NVList,VList,VCut,NCell1,HOC,LIS,CELLNEIGH,SLength,TYP,X,Y,XO,YO)
			ATotD = ATotD + AceptD
			ProbD = ATotD*100.0_dp/NTotD  ! Here we want float result	
		! If Swap Move
		else if (MoveType .gt. NPart/float(NMove) .and. MoveType .le. (NPart/10 + NPart)/float(NMove)) then
			! Check that it isn't a pure system. If so, skip swap move
			if (NPart .ne. NPartO .and. NPart .ne. NPartI .and. NPart .ne. NPartS) then
				NTotS = NTotS + 1
				call SWAPMOVE(AceptS,Beta,NPart,NPartO,NPartI,NPartS, & 
							NVList,VList,NCell1,HOC,LIS,CELLNEIGH,SLength,TYP,X,Y,XO,YO)
				ATotS = ATotS + AceptS
				ProbS = ATotS*100.0_dp/NTotS ! Here we want float result
			end if 	
		! If Area Move - aim for 30% acceptance - only do after some level of thermalization
		else if (MoveType .gt. (NPart/10 + NPart)/float(NMove) .and. iCycle .gt. 50) then
			NTotA = NTotA + 1
			call AREAMOVE(AceptA,Pres,Beta,DA,Area,SLength,NPart, & 
						NCell1,NVMax,VCut,NVList,VList,HOC,LIS,CELLNEIGH, &
						TYP,X,Y,XO,YO,HOCOUT,CELLNEIGHOUT)
			! If succesful, update HOC and CELLNEIGH here
			if (AceptA .eq. 1) then
				deallocate(HOC)
				deallocate(CELLNEIGH)
				
				allocate(HOC(size(HOCOUT)))
				allocate(CELLNEIGH(size(HOCOUT),9))
				HOC = HOCOUT
				CELLNEIGH = CELLNEIGHOUT
			end if
			ATotA = ATotA + AceptA
			ProbA = ATotA*100.0_dp/NTotA ! Here we want float result			
		end if
	end do inner_loop
	
	if (OUTCOUNTP .eq. NOutP) then
		!print *, 'Outputting data Cycle', iCycle
		OUTCOUNTP = 0
		call OUTPUTPOS(iCycle,NPart,TYP,X,Y)
	end if
	
	if (OUTCOUNTD .eq. NOutD) then
		!print *, 'Outputting data Cycle', iCycle
		OUTCOUNTD = 0
		call ENERGYSYS(Usys,NPart,SLength,NVLISt,VList,TYP,X,Y)
		call OUTPUTDAT(iCycle,Usys,Area,DX,DA,ProbD,ProbS,ProbA)
	end if	
	
	if (DCOUNT .eq. 100) then
		DCOUNT = 0
		if (ProbD .lt. 40.0_dp) then
			DX = DX*0.95_dp
		else if (ProbD .gt. 40.0_dp) then
			DX = DX*1.05_dp
		end if
		! Reset counters for displacement moves
		ProbD = 0.0_dp
		NTotD = 0
		ATotD = 0
		! Reset counters for swap moves
		ProbS = 0.0_dp
		ATotS = 0
		NTotS = 0
	end if
	
		if (ACOUNT .eq. 500) then
		ACOUNT = 0
		if (ProbA .lt. 30.0_dp) then
			DA = DA*0.95_dp
		else if (ProbA .gt. 30.0_dp) then
			DA = DA*1.05_dp
		end if
		! Reset counters for area moves
		ProbA = 0.0_dp
		NTotA = 0
		ATotA = 0
	end if
	
 end do outer_loop	
end program main
