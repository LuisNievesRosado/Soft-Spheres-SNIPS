module energy
 use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
 implicit none
 
 real(dp), parameter :: b = 6.0_dp      ! YK potential parameter 
 real(dp), parameter :: l = 2.0_dp      ! YK potential parameter 
 real(dp), parameter :: Rcut2 = 2.25_dp**2 ! Cutoff for YK potential, squared for convinience
 
 real(dp), parameter :: aO = 3.0_dp ! Interaction Parameter ISVO
 real(dp), parameter :: aI = 3.0_dp ! Interaction Parameter ISV
 real(dp), parameter :: aS = 3.7_dp ! Interaction Parameter SVPS 
 
 real(dp), parameter :: aOI = sqrt(aO*aI)! Interaction Parameter ISVO-ISV
 real(dp), parameter :: aOS = sqrt(aO*aS)! Interaction Parameter ISVO-SVPS
 real(dp), parameter :: aIS = sqrt(aI*aS)! Interaction Parameter ISV-SVPS
 
 real(dp), parameter :: fO = 1.12_dp  ! Size Parameter ISVO
 real(dp), parameter :: fI = 0.95_dp ! Size Parameter ISV
 real(dp), parameter :: fS = 0.699_dp ! Size Parameter SVPS 
 
 real(dp), parameter :: fOI = sqrt(fO*fI)! Size Parameter ISVO-ISV
 real(dp), parameter :: fOS = sqrt(fO*fS)! Size Parameter ISVO-SVPS
 real(dp), parameter :: fIS = sqrt(fI*fS)! Size Parameter ISV-SVPS
 
 real(dp), parameter :: eOI = 1.0_dp  ! Mix Parameter ISVO-ISV
 real(dp), parameter :: eOS = 1.075_dp ! Mix Parameter ISVO-SVPS
 real(dp), parameter :: eIS = 1.075_dp ! Mix Parameter ISV-SVPS
 
contains
 
 ! Subroutine to YK potential energy of pair of particles
 subroutine YKpairE(Upair,TY1,TY2,R)
	implicit none
	integer, intent(in) :: TY1 ! Input, Type of particle 1
	integer, intent(in) :: TY2 ! Input, Type of particle 2
	
	real(dp), intent(in) :: R  ! Input, Distance Between Particles
	
	real(dp), intent(out) :: Upair ! Output, Pair energy by YK Potential
	
	if (TY1 .eq. 1 .and. TY2 .eq. 1) then ! ISVO-ISVO
		Upair = exp(aO*(1.0_dp-R/fO) -b*(1.0_dp-R/fO)**l*log(R/fO))
	else if (TY1 .eq. 2 .and. TY2 .eq. 2) then ! ISV-ISV
		Upair = exp(aI*(1.0_dp-R/fI) -b*(1.0_dp-R/fI)**l*log(R/fI))
	else if (TY1 .eq. 3 .and. TY2 .eq. 3) then ! SVPS-SVPS
		Upair = exp(aS*(1.0_dp-R/fS) -b*(1.0_dp-R/fS)**l*log(R/fS))

	else if (TY1 .eq. 1 .and. TY2 .eq. 2 .or. TY1 .eq. 2 .and. TY2 .eq. 1) then ! ISVO-ISV
		Upair = eOI*exp(aOI*(1.0_dp-R/fOI) -b*(1.0_dp-R/fOI)**l*log(R/fOI))
	else if (TY1 .eq. 1 .and. TY2 .eq. 3 .or. TY1 .eq. 3 .and. TY2 .eq. 1) then ! ISVO-SVPS	
		Upair = eOS*exp(aOS*(1.0_dp-R/fOS) -b*(1.0_dp-R/fOS)**l*log(R/fOS))
	else if (TY1 .eq. 2 .and. TY2 .eq. 3 .or. TY1 .eq. 3 .and. TY2 .eq. 2) then ! ISV-SVPS
		Upair = eIS*exp(aIS*(1.0_dp-R/fIS) -b*(1.0_dp-R/fIS)**l*log(R/fIS))
	else
		print *, 'ERROR:Invalid type pairing found in energy calculation.'
		print *, TY1
		print *, TY2
	end if	 
 end subroutine YKpairE
 
 ! Subroutine to compute total energy of system. Used at start and after area move
 subroutine ENERGYSYS(Usys,NPart,SLength,NVList,VList,TYP,X,Y)
	implicit none
	integer, intent(in) :: NPart    		! Input, Total number of particles
	
	integer, intent(in) :: VList(:,:)       ! Input, Verlet List
	integer, intent(in) :: NVList(:)        ! Input, Number of Neighbors
	
	real(dp), intent(in) :: SLength			! Input, Box Length
	
	
	integer, intent(in) :: TYP(:) 			! Input, Type of each particle
	real(dp), intent(in) :: X(:)  			! Input, X locations
	real(dp), intent(in) :: Y(:)  			! Input, Y locations
	
	real(dp), intent(out) :: Usys 			! Output, Total energy of system
	
	
	real(dp) :: Upair            			! Internal, Where we will store pair energy
	real(dp) :: USUM             			! Internal, Where we will sum all the energies
	real(dp) :: XD               			! Internal, Distance in X Direction, adjusted for PBC
	real(dp) :: YD               			! Internal, Distance in Y Direction, adjusted for PBC
	real(dp) :: R				 			! Internal, Where we will keep track of pair distances
	
	integer :: iPart             			! Internal, Loop variable 1:NPart
	integer :: jPart			 			! Internal, Particle we're comparing i to
	integer :: jNeigh						! Internal, Loop variable over neighbors of i
	
	! Initialize energy
	USUM = 0
	part_loop: do iPart = 1,NPart
		neigh_loop: do jNeigh = 1,size(VList(iPart,:))
						jPart = VList(iPart,jNeigh)
						if (jPart .ne. 0) then
							XD = X(iPart)-X(jPart) - SLength*NInt((X(iPart)-X(jPart))/SLength)  ! Distance adjusting for PBC
							YD = Y(iPart)-Y(jPart) - SLength*NInt((Y(iPart)-Y(jPart))/SLength)  ! Distance adjusting for PBC					
							R = sqrt(XD**2 + YD**2) ! Compute distance
							call YKpairE(Upair,TYP(iPart),TYP(jPart),R) ! Compute pair energy
							USUM = USUM + Upair ! Track total energy
						end if
		end do neigh_loop
	end do part_loop
	Usys = USUM/2
	
	
	
	!==================================================================================================
	! Old code, using cell lists instead of verlet list. Now we use both, but here need only verlet.
	!==================================================================================================
	!part_loop: do iPart = 1,NPart
	!! Compute which cell we're in
	!	Cell = int((X[iPart] + SLength/2)*NCell1/SLength) + NCell1*int((Y[iPart] + SLength/2)*NCell1/SLength)
	!	
	!! Loop Over Neighboring Cells:
	!	cell_loop: do iCell = 1,9
	!! Loop Over Particles in Each Cell:	
	!		jCell = CELLNEIGH(Cell,iCell) ! Choose neighbor cell for j
	!		jPart = HOC(jCell) ! Initialize jPart for cell
	!		neigh_loop: do while (jPart .ne. 0)
	!			if (jPart .ne. iPart) then ! Check that its not the same particle
	!				XD = X(iPart)-X(jPart) - SLength*NInt((X(iPart)-X(jPart))/SLength)  ! Distance adjusting for PBC
	!				YD = Y(iPart)-Y(jPart) - SLength*NInt((Y(iPart)-Y(jPart))/SLength)  ! Distance adjusting for PBC
	!				
	!				R2 = XD**2 + YD**2 ! Compute distance squared
	!				if (R2 .le. Rcut2) then
	!					call YKpairE(Upair,TYP(iPart),TYP(jPart),sqrt(R2)) ! Compute pair energy
	!					USUM = USUM + Upair ! Track total energy
	!				end if
	!			end if	
	!			jPart = LIS(jPart)
	!		end do neigh_loop
	!	end do cell_loop
	!end do part_loop
	!
	!Usys = USUM/2
 end subroutine ENERGYSYS
 
 ! Subroutine to compute energy of particle. Used at every displacement and swap move
 subroutine ENERGYPART(Upart,Part,SLength,VListP,TYP,X,Y,XP,YP)
	implicit none
	
	integer, intent(in) :: Part           ! Input, Particle of interest	
	integer, intent(in) :: TYP(:)         ! Input, Type of each Particle
	integer, intent(in) :: VListP(:)	  ! Input, Verlet list, here only of particle of interest
	
	real(dp), intent(in) :: SLength		  ! Input, Side length of box
	
	real(dp), intent(in) :: XP		      ! Input, Position of particle in X
	real(dp), intent(in) :: YP		 	  ! Input, Position of particle in Y
	
	real(dp), intent(in) :: X(:)          ! Input, X Coordinate of Each Particle
	real(dp), intent(in) :: Y(:)          ! Input, Y Coordinate of Each Particle
	
	real(dp), intent(out) :: Upart		  ! Output, energy of particle
	
	real(dp) :: Upair            		  ! Internal, Where we will store pair energy
	real(dp) :: XD               		  ! Internal, Distance in X Direction, adjusted for PBC
	real(dp) :: YD               		  ! Internal, Distance in Y Direction, adjusted for PBC
	real(dp) :: R				 		  ! Internal, Real distance, for when we need to compute energy

	integer :: jPart			 		  ! Internal, Particle we're comparing Part to
	integer :: jNeigh            		  ! Internal, Loop variable over verlet list
	

	! Initialize energy
	Upart = 0
	do jNeigh=1,size(VListP)
		jPart = VListP(jNeigh)
		if (jPart .ne. 0) then
			XD = XP-X(jPart) - SLength*NInt((XP-X(jPart))/SLength)  ! Distance adjusting for PBC
			YD = YP-Y(jPart) - SLength*NInt((YP-Y(jPart))/SLength)  ! Distance adjusting for PBC					
			R = sqrt(XD**2 + YD**2) ! Compute distance
			call YKpairE(Upair,TYP(Part),TYP(jPart),R) ! Compute pair energy
			Upart = Upart + Upair ! Track total energy		
			
		end if
	end do
	! ================================================================================!
	! Commented out Cell Based Method, Now using verlet lists with cell lists
	
	!Cell = int((X[Part] + SLength/2)*NCell1/SLength) + NCell1*int((Y[Part] + SLength/2)*NCell1/SLength)
	!
	!! Initialize energy
	!USUM = 0
	!cell_loop: do iCell = 1,9
	!! Loop Over Particles in Each Cell:	
	!	jCell = CELLNEIGH(Cell,iCell) ! Choose neighbor cell for j
	!	jPart = HOC(jCell) ! Initialize jPart for cell
	!	neigh_loop: do while (jPart .ne. 0)
	!		if (jPart .ne. Part) then ! Check that its not the same particle
	!			XD = X(Part)-X(jPart) - SLength*NInt((X(Part)-X(jPart))/SLength)  ! Distance adjusting for PBC
	!			YD = Y(Part)-Y(jPart) - SLength*NInt((Y(Part)-Y(jPart))/SLength)  ! Distance adjusting for PBC
	!			
	!			R2 = XD**2 + YD**2 ! Compute distance squared
	!			if (R2 .le. Rcut2) then
	!				call YKpairE(Upair,TYP(Part),TYP(jPart),sqrt(R2)) ! Compute pair energy
	!				USUM = USUM + Upair ! Track total energy
	!			end if
	!		end if	
	!		jPart = LIS(jPart)
	!	end do neigh_loop
	!end do cell_loop	
	!
	!Upart = USUM
	! ================================================================================!
 end subroutine ENERGYPART
 
end module energy
 
 