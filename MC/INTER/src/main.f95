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
 integer :: NVMax  ! Maximum Number of Verlet Neighbors per Particle
 
 real(dp) :: IPerc ! Percentage of ISV in ISV/O Mix
 real(dp) :: SLengthx ! Length of side
 
 real(dp) :: Temp  ! Desired Temperature of System
 real(dp) :: Pres  ! Desired Pressure of System

 real(dp) :: DX    ! Displacement Move Magnitude
 real(dp) :: DA    ! Area Move Magnitude
 
 real(dp) :: Vcut  ! Cutoff for verlet list
 
 
 
 !===========================================================================!
 ! Internal Variables
 !===========================================================================!
 
 real(dp), parameter :: LI = 1.15_dp   ! ISV/O Distance
 real(dp), parameter :: LSx = 1.0_dp   ! SVPS Distance X
 real(dp), parameter :: LSy = 0.85_dp  ! SVPS Distance Y 
 
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
 
 real(dp) :: SLengthy ! Side length in y
 real(dp) :: Area     ! Area of System
 
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
 
 
 integer :: NCellx   ! Number of cells in each direction x
 integer :: NCelly   ! Number of cells in each direction y
 
 
 integer :: IniCycle ! Set to 0 for output
 integer, allocatable :: VList(:,:)   ! Verlet List of neigbors
 integer, allocatable :: NVList(:)    ! Number of neighbors of each particle
 
 real(dp), dimension(12) :: INP ! Input file data storage

  
  
 real(dp) :: dY      ! Initial distance from y = 0 
 
 integer :: NRI      ! Number of rows of ISV/O
 integer :: NRS      ! Number of rows of SVPS 
  
 integer :: NCI      ! Number of columns of ISV/O
 integer :: NCS      ! Number of columns of SVPS
 
 integer :: NPartI   ! Number of ISV/O
 integer :: NPartS   ! Number of SVPS
 
 
 real(dp) :: IOTYP   ! Choose between ISV and ISVO
 
 ! Loop integers
 integer :: PartC    ! Loop variable for particle setup
 integer :: iXI      ! Loop variable for particle setup
 integer :: iYI      ! Loop variable for particle setup
 integer :: iXS      ! Loop variable for particle setup
 integer :: iYS      ! Loop variable for particle setup
 integer :: iPart    ! Loop variable for particle setup
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
 do iRead = 1,12
	read(1,*) INP(iRead)
 end do	
 close(1)

 
 
 NCycle = int(INP(1))
 NTherm = int(INP(2))
 NOutD = int(INP(3))
 NOutP = int(INP(4))
 NVMax = int(INP(5))
 
 IPerc = INP(6)
 SLengthx = INP(7)
 Temp = INP(8)
 Pres = INP(9)

 DX = INP(10)
 DA = INP(11)
 Vcut = INP(12)
 

 !===========================================================================!
 ! Initial System Setup - Box
 !===========================================================================!
 SLengthy = 4*SLengthx
 Area = SLengthx*SLengthy 
 
 !===========================================================================!
 ! Initial System Setup - Particles
 !===========================================================================!
 
 dY = sqrt(LI*LSy)
 
 NCI = int(SLengthx/LI)-1
 NCS = int(SLengthx/LSx)-1
 
 NRI = int((SLengthy/2 - dY)/LI)
 NRS = int((SLengthy/2 - dY)/LSy)
 
 NpartI = NCI*NRI
 NPartS = NCS*NRS
 

 
 NPart = NPartI + NPartS
 
 ! Allocate Types of Particles
 allocate(TYP(NPart))
 
 ! Allocate Particle Locations
 allocate(X(NPart))
 allocate(Y(NPart)) 
 

 
 
 ! Initialize square region
 PartC = 1
 do iXI = 1, NCI	
	do iYI = 1, NRI
		X(PartC) = iXI*LI - Slengthx/2
		Y(PartC) = (iYI - 1)*LI + dY/2
		call RANDOM_NUMBER(IOTYP)
		if (IOTYP .gt. IPerc) then
			TYP(PartC) = 2
		else
			TYP(PartC) = 1
		end if	
		PartC = PartC + 1
	end do
 end do

 ! Initialize hex region
 
 do iXS = 1, NCS	
	do iYS = 1, NRS 	
		if (MOD(iYS,2) .eq. 0) then
			X(PartC) = iXS*LSx  - LSx/2 - Slengthx/2
		else
			X(PartC) = iXS*LSx - Slengthx/2
		end if	
		Y(PartC) = -(iYS - 1)*LSy - dY/2	
		TYP(PartC) = 3
		PartC = PartC + 1
	end do
 end do
  
  ! Save Initial Locations
 allocate(XO(NPart))
 allocate(YO(NPart))

 do iPart = 1,NPart
	XO(iPart) = X(iPart)
	YO(iPart) = Y(iPart)
 end do


  !===========================================================================!
 ! Initial System Setup - Simulation
 !===========================================================================!
 ! Based on Prajwal Recomendation
 NMove = NPart + NPart/10 + 2  ! NPart Disp Moves, NPart/10 Swap Moves, 2 Area Moves
 
 Beta = 1.0_dp/Temp 
 !===========================================================================!
 ! Initial System Setup - Cell Lists
 !===========================================================================!
 
 allocate(HOC(int(SLengthx/Rcut)*int(SLengthy/Rcut)))
 allocate(LIS(NPart))
 allocate(CELLNEIGH(int(SLengthx/Rcut)*int(SLengthy/Rcut),9))

 call MAKECELL(NPart,SLengthx,Slengthy,HOC,LIS,CELLNEIGH,NCellx,NCelly,X,Y)
 
 !===========================================================================
 ! Initial System Setup - Neighbor Lists
 !===========================================================================!
 allocate(VList(NPart,NVMax))
 allocate(NVList(NPart))
 
 call MAKELIST(VList,NVList,NVMax,Vcut,HOC,LIS,CELLNEIGH,NCellx,NCelly,NPart,SLengthx,Slengthy,X,Y)
 !===========================================================================!
 ! Initial System Setup - Energy
 !===========================================================================!
 call ENERGYSYS(Usys,NPart,SLengthx,Slengthy,NVLISt,VList,TYP,X,Y)
 
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
			call DISPLACEMOVE(AceptD,Beta,DX,NPart,NVList,VList,VCut,NCellx,Ncelly,HOC,LIS,CELLNEIGH,SLengthx,Slengthy,TYP,X,Y,XO,YO)
			ATotD = ATotD + AceptD
			ProbD = ATotD*100.0_dp/NTotD  ! Here we want float result	
		! If Swap Move
		else if (MoveType .gt. NPart/float(NMove) .and. MoveType .le. (NPart/10 + NPart)/float(NMove)) then
			NTotS = NTotS + 1
			call SWAPMOVE(AceptS,Beta,NPart,NVList,VList,NCellx,Ncelly,HOC,LIS,CELLNEIGH,SLengthx,Slengthy,TYP,X,Y,XO,YO)
			ATotS = ATotS + AceptS
			ProbS = ATotS*100.0_dp/NTotS ! Here we want float result
		! If Area Move - aim for 30% acceptance - only do after some level of thermalization
		else if (MoveType .gt. (NPart/10 + NPart)/float(NMove) .and. iCycle .gt. 50) then
			NTotA = NTotA + 1
			call AREAMOVE(AceptA,Pres,Beta,DA,Area,SLengthx,Slengthy,NPart, & 
						NCellx,Ncelly,NVMax,VCut,NVList,VList,HOC,LIS,CELLNEIGH, &
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
		call ENERGYSYS(Usys,NPart,SLengthx,Slengthy,NVLISt,VList,TYP,X,Y)
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
