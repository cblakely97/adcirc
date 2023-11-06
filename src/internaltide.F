!-------------------------------------------------------------------------------!
!
! ADCIRC - The ADvanced CIRCulation model
! Copyright (C) 1994-2023 R.A. Luettich, Jr., J.J. Westerink
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!-------------------------------------------------------------------------------!
C---------------------------------------------------------------------C
      MODULE INTERNALTIDE
C---------------------------------------------------------------------C
C     CPB 03/2023: This module was created to make the
C     apply2dinternalwavedrag subroutine in nodalattr.F more readable.
C     For now I have simply moved the calculation of the de-tided
C     velocities here. I plan to move more of the subroutine into this
C     module to hopefully make it a little clearer what is happening in
C     that subroutine.
C---------------------------------------------------------------------C
      USE GLOBAL, ONLY : DEBUG, ECHO, INFO, WARNING, ERROR,
     &          setMessageSource, unsetMessageSource, allMessage,
     &          scratchMessage, logMessage      
      IMPLICIT NONE
      ! stored velocity samples for moving average calculation
      REAL(8),ALLOCATABLE,PRIVATE :: UAV(:,:), VAV(:,:)        
      ! averaged velocities so we don't calculate them every single time
      ! step
      REAL(8),ALLOCATABLE :: UBar(:), VBar(:)
      ! tidal velocities
      REAL(8),ALLOCATABLE :: UTide(:), VTide(:)
      ! weights for filters
      REAL(8),ALLOCATABLE,PRIVATE :: wts(:)
      ! CPB: testing short-term harmonic resynthesis
      ! mass matrix and rhs load vector for short term analysis. since
      ! we are solving for both u and v directions hb will have size 
      ! 2 by 2xnconst by number of nodes
      REAL(8),ALLOCATABLE,PRIVATE :: HA(:,:), hb(:,:,:)
      ! vectors to store the solution to the least squares analysis
      REAL(8),ALLOCATABLE,PRIVATE :: c_i(:,:,:), s_i(:,:,:)
      REAL(8),PRIVATE :: LSQ_L  ! length of analysis
      REAL(8),PRIVATE :: filtL  ! length of filter used before least
                                ! squares analysis can be done
      REAL(8),PRIVATE :: Fs     ! sampling frequency
      INTEGER,PRIVATE :: NS     ! number of samples
      INTEGER,PRIVATE :: nconst ! number of constituents considered
      REAL(8),ALLOCATABLE :: freq(:) ! frequencies we want to use


C---------------------------------------------------------------------C
      CONTAINS
C---------------------------------------------------------------------C
!---------------------------------------------------------------------!
      SUBROUTINE SHORT_LSQ_HA(Uin,Vin, TimeStep)
!---------------------------------------------------------------------!
!     subroutine that takes the current y value and timetep and, if   !
! it is time, performs the least squares analysis.                    !
!---------------------------------------------------------------------!
      USE MESH, ONLY : NP
      USE GLOBAL, ONLY : DTDP
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: Uin(:), Vin(:)
      INTEGER,INTENT(IN) :: TimeStep
      INTEGER :: L, Lm ! for checking if it is time to update/solve
      ! first_call is my hacky-ass way of initializing. I hopefully will
      ! fix it before actually implementing in ADCIRC
      LOGICAL,SAVE :: first_call = .true.
      ! This next flag is also a bit hacky but I use it to tell the
      ! program when the first HA has been performed. 
      LOGICAL :: first_solve = .false.
      ! these next two flags serve to check if it is time to update the
      ! mass matrix and load vector (Ax=b) and also if it is time to
      ! solve the harmonic analysis
      LOGICAL :: update_flag = .false., solve_flag = .false.
      INTEGER,SAVE :: solve_counter = 0
      ! To avoid roundoff issues in calculating sines/cosines we will
      ! calculate a "harmonic analysis time" that is in units of seconds
      ! since the last harmonic analysis. We will calculate this using
      ! solve_counter
      REAL(8) :: HA_Time
      INTEGER :: ii
      !
      ! This should really be moved to an initialization routine. If I
      ! have time I will make one once I have it all working
      !
      IF ( first_call ) THEN
         CALL INIT_SHORT_LSQ_HA()
         first_call = .false.
      ENDIF
      !
      ! Calculate harmonic time before solve_counter is messed with
      !
      HA_time = TimeStep*DTDP
      !
      ! Check if it is time to update/solve
      !
      ! to check if it is time to update check if this timestep and the
      ! previous timestep are in different sampling periods
      L = FLOOR(TimeStep*DTDP/Fs)
      Lm = FLOOR((TimeStep-1)*DTDP/Fs)
      IF ( L.GT.Lm ) THEN
         update_flag = .true.
      ELSE
         update_flag = .false.
      ENDIF
      ! to check if it is time to solve we just use a counter to track
      ! how long it has been since we last solved the analysis
      L = solve_counter*DTDP ! seconds since last time we solved
      IF ( L.EQ.LSQ_L ) THEN
         solve_flag = .true.
         solve_counter = 0
      ELSE
         solve_flag = .false.
      ENDIF
      !
      ! Update LA25 filter if we haven't gotten to the LSQHA
      !
      IF (TimeStep*DTDP.LT.LSQ_L) THEN
         first_solve = .false.
         CALL UNTIDE(Uin,Vin,TimeStep)
      ELSE
         first_solve = .true.
      ENDIF
      ! 
      ! Update A matrix and b vector
      !
      CALL UPDATE_LSQ_HA(Uin,Vin,HA_Time,update_flag)
      !
      ! Solve the least squares analysis
      !
      CALL SOLVE_LSQ_HA(solve_flag)
      !
      ! Estimate output
      !
      CALL ESTIMATE_TIDAL_VELO(Uin,Vin,HA_time,first_solve)
      solve_counter = solve_counter + 1
      return
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      ENDSUBROUTINE SHORT_LSQ_HA
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      SUBROUTINE INIT_SHORT_LSQ_HA()
!---------------------------------------------------------------------!
!     subroutine that takes the current y value and timetep and, if   !
! it is time, performs the least squares analysis.                    !
!---------------------------------------------------------------------!
      USE MESH, ONLY : NP
      USE CONSTANTS, ONLY : PI
      IMPLICIT NONE
      REAL(8),PARAMETER :: cpd2radpsec = 2d0*pi/86400
      !
      ! Initialize everything we need for the least squares harmonic
      ! analysis
      !
      ! number of constituents considered
      nconst = 2
      ! frequencies used
      ALLOCATE( freq(nconst) )
      freq = (/1.9323, 1.0027/) ! diurnal (M2) and semidiurnal (K1)
      freq = freq*cpd2radpsec            ! convert to radians per second
      ! allocate mass matrix and RHS vector for each node
      ALLOCATE( HA(2*nconst,2*nconst), hb(2,2*nconst,NP),
     &          c_i(2,nconst,NP), s_i(2,nconst,NP) )
      HA = 0d0
      hb = 0d0
      c_i = 0d0
      s_i = 0d0
      ! 
      Fs = 720d0 ! sampling rate in seconds (12 min)
      filtL = 90000d0 ! 25 hour filter (seconds)
      LSQ_L = 172800d0 ! 2 days for each analysis (seconds)
      NS = FLOOR(filtL/Fs)
      ALLOCATE( UAV(NS,NP), VAV(NS,NP), 
     &          UBar(NP), VBar(NP),
     &          UTide(NP), VTide(NP) )
      UAV = 0d0
      VAV = 0d0
      UBar = 0d0
      VBar = 0d0
      UTide = 0d0
      VTide = 0d0
      RETURN
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      ENDSUBROUTINE INIT_SHORT_LSQ_HA
!---------------------------------------------------------------------!
C-----------------------------------------------------------------------
      SUBROUTINE UNTIDE(U_i,V_i,TimeStep)
C---------------------------------------------------------------------C
C     This subroutine takes as inputs the current time step velocity,
C     updates the average velocity vectors (if necessary) and outputs
C     the 25 hour lagged average velocity based on a lagged, 25 hour
C     filter with sampling frequency of 12 minutes. It replaces a good
C     chunk of code in the apply2dinternalwavedrag to make that
C     subroutine more readable as well as to eliminate some global
C     variables
C---------------------------------------------------------------------C
      USE MESH, ONLY : NP
      USE GLOBAL, ONLY: DTDP
      IMPLICIT NONE
      REAL(8),INTENT(IN),DIMENSION(:) :: U_i, V_i ! current velocity
      INTEGER,INTENT(IN) :: TimeStep    ! current time step
      LOGICAL,SAVE :: first_call = .TRUE.
      ! need indices to check if the current and previous time step are
      ! in the same 12 min windo2
      INTEGER :: L, Lm           
      INTEGER :: ii, kk ! for loops
      ! for populating the UAV, VAV matrices at the start of the run
      INTEGER,SAVE :: ISTA

      IF ( first_call ) THEN
         ista = 1
         first_call = .false.
      ENDIF

      L = FLOOR(TimeStep*DTDP/Fs)
      Lm = FLOOR((TimeStep-1)*DTDP/Fs)
      IF ( L.GT.Lm ) THEN
         IF ( ISTA.GT.NS ) THEN
            DO ii = 1,NP
               DO kk = 1,NS-1
                  UAV(ii,kk) = UAV(ii,kk+1)
                  VAV(ii,kk) = VAV(ii,kk+1)
               ENDDO
               UAV(ii,NS) = U_i(ii)
               VAV(ii,NS) = V_i(ii)
            ENDDO
            DO ii = 1,NP
               UBar(ii) = SUM(UAV(ii,1:NS))/DBLE(NS)
               VBar(ii) = SUM(VAV(ii,1:NS))/DBLE(NS)
            ENDDO
         ELSE
            DO ii = 1,NP
               UAV(ii,ISTA) = U_i(ii)
               VAV(ii,ISTA) = V_i(ii)
               UBar(ii) = SUM(UAV(ii,1:ISTA))/DBLE(ISTA)
               VBar(ii) = SUM(VAV(ii,1:ISTA))/DBLE(ISTA)
            ENDDO
            ISTA = ISTA + 1
         ENDIF
      ENDIF
      RETURN
C-----------------------------------------------------------------------
      END SUBROUTINE UNTIDE
C-----------------------------------------------------------------------
!---------------------------------------------------------------------!
      SUBROUTINE UPDATE_LSQ_HA(Uin,Vin,HA_Time,flag)
!---------------------------------------------------------------------!
!     Checks if it is time to update the mass matrix and load vector  !
!     for the harmonic analysis. If it is it then updates
!---------------------------------------------------------------------!
      USE MESH, ONLY : NP
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: Uin(:), Vin(:)
      REAL(8),INTENT(IN) :: HA_Time
      LOGICAL,INTENT(IN) :: flag
      INTEGER :: ii, jj ! for loops
      INTEGER :: i1, i2, j1, j2 ! for indices
      REAL(8) :: cosii, cosjj, sinii, sinjj
      ! return if it is not time to update
      IF ( .NOT.flag ) THEN
         RETURN
      ENDIF
      ! update LHS matrix
      DO ii = 1,nconst
         cosii = COS(freq(ii)*HA_Time)
         sinii = SIN(freq(ii)*HA_Time)
         i1 = 2*ii-1
         i2 = i1+1
         DO jj = ii,nconst
            cosjj = COS(freq(jj)*HA_Time)
            sinjj = SIN(freq(jj)*HA_Time)
            j1 = 2*jj-1
            j2 = j1+1
            HA(i1,j1) = HA(i1,j1) + cosii*cosjj
            HA(i1,j2) = HA(i1,j2) + cosii*sinjj
            HA(i2,j2) = HA(i2,j2) + sinii*sinjj
            IF (i2.LE.j1) THEN
               HA(i2,j1) = HA(i2,j1) + sinii*cosjj
            ENDIF
         ENDDO
         ! load vectors
         DO jj = 1,NP
            hb(1,i1,jj) = hb(1,i1,jj) + Uin(jj)*cosii
            hb(1,i2,jj) = hb(1,i2,jj) + Uin(jj)*sinii
            hb(2,i1,jj) = hb(2,i1,jj) + Vin(jj)*cosii
            hb(2,i2,jj) = hb(2,i2,jj) + Vin(jj)*sinii
         ENDDO
      ENDDO

      RETURN
!---------------------------------------------------------------------!
      ENDSUBROUTINE UPDATE_LSQ_HA
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      SUBROUTINE SOLVE_LSQ_HA(flag)
!---------------------------------------------------------------------!
!     Checks if it is time to solve for the harmonic constituents.    !
!     If it is then decomposes mass matrix using LU decomp and solves !
!     for all points of interest. 
!---------------------------------------------------------------------!
      USE MESH, ONLY : NP
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: flag
      INTEGER :: ii, jj, i1,i2 ! for loops
      INTEGER,ALLOCATABLE :: pivot(:)
      ! return if it is not time to solve
      IF ( .NOT.flag ) THEN
         RETURN
      ENDIF
      ! 
      ! Populate the lower left part of the mass matrix.
      !
      DO jj = 1,2*nconst
         DO ii = jj,2*nconst
            HA(ii,jj) = HA(jj,ii)
         ENDDO
      ENDDO
      ! 
      ! We will solve using LU Decomposition
      !
      ! Pivot matrix allows for non-diagonal pivot points. This ensures
      ! that the LU decomposition converges.
      ALLOCATE( pivot(2*nconst) )
      ! decompose matrix
      CALL LU_DECOMP(HA,pivot)
      !
      ! Solve for each point of interest
      !
      DO jj = 1,NP
         ! u-direction
         CALL LU_SOLVE(HA,hb(1,:,jj),pivot)
         DO ii = 1,nconst
            i1 = 2*ii-1
            i2 = i1+1
            c_i(1,ii,jj) = hb(1,i1,jj)
            s_i(1,ii,jj) = hb(1,i2,jj)
         ENDDO
         ! v-direction
         CALL LU_SOLVE(HA,hb(2,:,jj),pivot)
         DO ii = 1,nconst
            i1 = 2*ii-1
            i2 = i1+1
            c_i(2,ii,jj) = hb(2,i1,jj)
            s_i(2,ii,jj) = hb(2,i2,jj)
         ENDDO
      ENDDO
      HA = 0d0
      hb = 0d0
      DEALLOCATE( pivot )
      RETURN
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      ENDSUBROUTINE SOLVE_LSQ_HA
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      SUBROUTINE ESTIMATE_TIDAL_VELO(Uin,Vin,HA_time,first_solve)
!---------------------------------------------------------------------!
!     Returns an estimated tidal output. If we have performed a       !
!     harmonic analysis already then it estimates the output at       !
!     HA_time. If not, it estimates the tidal output using a lagged   !
!     25-hour average filter.
!---------------------------------------------------------------------!
      USE MESH, ONLY : NP
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: Uin(:),Vin(:)
      REAL(8),INTENT(IN) :: HA_Time
      LOGICAL,INTENT(IN) :: first_solve
      ! local
      INTEGER :: kk, iconst, nconst
      REAL(8) :: cost,sint
      UTide = 0d0
      VTide = 0d0
      IF (.NOT.first_solve) THEN
         ! use LA25 filter
         DO kk = 1,NP
            UTide(kk) = Uin(kk) - UBar(kk)
            VTide(kk) = Vin(kk) - VBar(kk)
         ENDDO
         RETURN
      ENDIF
      ! resynthesize
      DO kk = 1,NP
         DO iconst = 1,nconst
            cost = c_i(1,iconst,kk)*cos(freq(iconst)*HA_Time)
            sint = s_i(1,iconst,kk)*sin(freq(iconst)*HA_Time)
            UTide(kk) = UTide(kk) + cost + sint
            cost = c_i(2,iconst,kk)*cos(freq(iconst)*HA_Time)
            sint = s_i(2,iconst,kk)*sin(freq(iconst)*HA_Time)
            VTide(kk) = VTide(kk) + cost + sint
         ENDDO
      ENDDO
      RETURN
!---------------------------------------------------------------------!
      ENDSUBROUTINE ESTIMATE_TIDAL_VELO
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      SUBROUTINE LU_DECOMP(A,Pivot)
!---------------------------------------------------------------------!
!     Performs the LU decomposition of mass matrix A such that        !
!     L*U = A. Note that the matrices LU are stored in A such that:   !
!        A = | u_{1,1}, u_{1,2}, u_{1,3}, u_{1,4} |                   !
!            | l_{2,1}, u_{2,2}, u_{2,3}, u_{2,4} |                   !
!            | l_{3,1}, l_{3,2}, u_{3,3}, u_{3,4} |                   !
!            | l_{4,1}, l_{4,2}. l_{4,3}, u_{4,4} |                   !
!     and l_{i,i} = 1                                                 !
!     This decomposition can then be used to solve for every node.    !
!---------------------------------------------------------------------!
      IMPLICIT NONE
      ! 
      ! arguments
      ! 
      ! input/output matrix
      REAL(8),DIMENSION(:,:),INTENT(INOUT) :: A
      ! permutations
      INTEGER,DIMENSION(:),INTENT(OUT) :: Pivot
      !
      ! local variables
      !
      INTEGER :: i,j,k,n
      INTEGER,DIMENSION(:),ALLOCATABLE :: p
      REAL,DIMENSION(:),ALLOCATABLE :: Tmp ! temporary row
      n = size(A,1)
      ALLOCATE( Tmp(n), p(n) )
      !
      DO j = 1,n-1
         p = MAXLOC(ABS(A(J:n,j)))+j-1
         IF (p(1).NE.j) THEN
            Tmp(:) = A(j,:)
            A(j,:) = A(p(1),:)
            A(p(1),:) = Tmp(:)
            Pivot(j) = p(1)
         ELSE
            Pivot(j) = j
         ENDIF
         DO i = j+1,n
            A(i,j) = A(i,j)/A(j,j)
            A(i,j+1:n) = A(i,j+1:n) - A(i,j)*A(j,j+1:n)
         ENDDO
      ENDDO
      DEALLOCATE( Tmp, p )
      RETURN
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      END SUBROUTINE LU_DECOMP
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      SUBROUTINE LU_SOLVE(A,b,pivot)
!---------------------------------------------------------------------!
!     Solves the system of linear equations Ax=b given:               !
!        - LU decomposition of A (stored compactly)                   !
!        - RHS load vector b                                          !
!        - Vector of pivots if LU decomp did not pivot on diagonal    !
!     Returns:                                                        !
!        - Solution to system stored in b                             !
!     First solves Ly=b by forward substitution and then solves       !
!     Ux=y by backwards substitution.                                 !
!---------------------------------------------------------------------!
      IMPLICIT NONE
      ! 
      ! arguments
      ! 
      ! A-matrix decomposed into LU and stored compactly
      REAL(8),DIMENSION(:,:),INTENT(IN) :: A
      ! vector with stored pivots if they are not the diagonal
      INTEGER,DIMENSION(:),INTENT(IN) :: Pivot
      ! RHS Load vector
      REAL(8),DIMENSION(:),INTENT(INOut) :: b
      !
      ! local variables
      !
      INTEGER :: i,j,n
      REAL(8) :: Tmp ! temporary row
      n = size(A,1)
      !
      DO j = 1,n-1
         IF (Pivot(j).NE.j) THEN
            Tmp = b(j)
            b(j) = b(Pivot(j))
            b(Pivot(j)) = Tmp
         ENDIF
      ENDDO
      !
      ! Forward substitution
      !
      DO i = 2,n
         b(i) = b(i) - DOT_PRODUCT(A(i,1:i-1),b(1:i-1))
      ENDDO
      !
      ! Backward substitution
      !
      b(n) = b(n)/A(n,n)
      DO i = n-1,1,-1
         b(i) = (b(i) - DOT_PRODUCT(A(i,i+1:n),b(i+1:n)))/A(i,i)
      ENDDO
      RETURN
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      ENDSUBROUTINE LU_SOLVE
!---------------------------------------------------------------------!
      SUBROUTINE MunkHPFilter(U_i,V_i,TimeStep)
C---------------------------------------------------------------------C
C     This subroutine applies a high-pass filter derived from the
C     so-called Munk "Tide Killer" filter to the velocity field. The 
C     coefficients for the original, low pass filter can be found at: 
C
C        https://www.sonel.org/Filters-for-the-daily-mean-sea.html
C
C     The high-pass filter is derived from the normalized, nonrecursive
C     low-pass filter as:
C
C        W^{HP}_0 = 1-W^{LP}_0
C        W^{HP}_k = -W^{LP}_k   (k not equal to 0)
C
C     The filter is applied as:
C
C        y_n = \sum_{k=-m}^m W_k * x_{n+k}
C
C     Note that the output of the filter is lagged by 25 hours from the
C     current time step (which is fine since we are going to apply this
C     to the internal wave drag and 25 hours is approximately twice the
C     semi-diurnal period).
C
C     Written by: Coleman Blakely 3/2023
C---------------------------------------------------------------------C
      USE MESH, ONLY : NP
      USE GLOBAL, ONLY: DTDP
      IMPLICIT NONE
      REAL(8),INTENT(IN),DIMENSION(:) :: U_i, V_i ! current velocity
      INTEGER,INTENT(IN) :: TimeStep    ! current time step
      REAL(8) :: filtL = 49d0*3600d0        ! filter length (s)
      REAL(8) :: T = 60d0*60d0             ! sampling interval (s)
      INTEGER,SAVE :: NS = 1   ! number of samples
      LOGICAL,SAVE :: first_call = .TRUE.
      ! need indices to check if the current and previous time step are
      ! in the same 1 hour window
      INTEGER :: L, Lm           
      INTEGER :: ii, kk ! for loops
      ! for populating the UAV, VAV matrices at the start of the run
      INTEGER,SAVE :: ISTA

      IF ( first_call ) THEN
         first_call = .FALSE.
         NS = FLOOR(filtL/T)
         ista = 1
         ALLOCATE( UAV(NP,NS), VAV(NP,NS) )
         UAV = 0d0
         VAV = 0d0
         ALLOCATE( UBar(NP), VBar(NP) )
         UBar = 0d0
         VBar = 0d0
         CALL CalcMunkWeights()
      ENDIF
      ! check if this and the previous timestep are in the same hour
      L = FLOOR(TimeStep*DTDP/T)
      Lm = FLOOR((TimeStep-1)*DTDP/T)
      IF ( L.GT.Lm ) THEN
         ! check if we are >= 49 hours into the run
         IF ( ISTA.GT.NS ) THEN
            DO ii = 1,NP
               UAV(ii,1:NS-1) = UAV(ii,2:NS)
               VAV(ii,1:NS-1) = VAV(ii,2:NS)
               UAV(ii,NS) = U_i(ii)
               VAV(ii,NS) = V_i(ii)
            ENDDO
            UBar = 0d0
            VBar = 0d0
            DO ii = 1,NP
               DO kk = 1,NS
                  UBar(ii) = UBar(ii) + wts(kk)*UAV(ii,kk)
                  VBar(ii) = VBar(ii) + wts(kk)*VAV(ii,kk)
               ENDDO
            ENDDO
         ELSE
            ! if we do not have a 49 hr record yet just use a lagged
            ! average
            DO ii = 1,NP
               UAV(ii,ISTA) = U_i(ii)
               VAV(ii,ISTA) = V_i(ii)
               UBar(ii) = U_i(ii) - SUM(UAV(ii,1:ISTA))/DBLE(ISTA)
               VBar(ii) = V_i(ii) - SUM(VAV(ii,1:ISTA))/DBLE(ISTA)
            ENDDO
            ISTA = ISTA + 1
         ENDIF
      ENDIF
C---------------------------------------------------------------------C
      END SUBROUTINE MunkHPFilter
C---------------------------------------------------------------------C
C---------------------------------------------------------------------C
      SUBROUTINE CalcMunkWeights()
C---------------------------------------------------------------------C
C     This subroutine sets up the filter weights for use in the
C     high-pass filter derived from the so-called Munk "Tide Killer"
C     low-pass filter. 
C---------------------------------------------------------------------C
      IMPLICIT NONE
      REAL(8),DIMENSION(25) :: LPwts   ! original weights (one sided)
      INTEGER :: NS = 49                ! length of HP filter needed
      REAL(8) :: K                     ! for normalizing the LP filter
      INTEGER :: ii                     ! for loops
      !
      ! allocate weights
      ALLOCATE( wts(NS) )
      ! define low-pass filter weights (not normalized and one-sided)
      LPwts = (/395287,386839,370094,354118,338603,325633,314959,
     &          300054,278167,251492,234033,219260,208050,195518,
     &          180727,165525,146225,122665,101603,85349,72261,
     &          60772,47028,30073,13307 /)
      ! find sum to normalize low pass filter
      K = LPwts(1)
      DO ii = 2,25
         K = K + 2d0*LPwts(ii)
      END DO
      ! normalize low pass filter weights
      DO ii = 1,25
         LPwts(ii) = LPwts(ii)/K
      END DO
      ! turn one-sided low-pass filter into two-sided high-pass filter
      DO ii = 1,24
         ! -m to -1
         wts(ii) = -LPwts(26-ii)
         ! 1 to m
         wts(ii+25) = -LPwts(ii+1)
      END DO
      ! 0 (center point)
      wts(25) = 1 - LPwts(1)
C---------------------------------------------------------------------C
      END SUBROUTINE CalcMunkWeights
C---------------------------------------------------------------------C
C-----------------------------------------------------------------------
C     S U B R O U T I N E   C H E C K  _  E R R
C-----------------------------------------------------------------------
C     jgf49.17.02 Checks the return value from netCDF calls; if there
C     was an error, it writes the error message to the screen and to the
C     fort.16 file.
C-----------------------------------------------------------------------
      subroutine check_err(iret)
      USE SIZES, ONLY : myproc
      USE GLOBAL, ONLY : screenUnit
#ifdef CMPI
      USE MESSENGER, ONLY : MSG_FINI
#endif
#ifdef ADCNETCDF
      USE NETCDF
#endif
      IMPLICIT NONE

      INTEGER, intent(in) :: iret

      !call setMessageSource("check_err")
#ifdef ALL_TRACE
      call allMessage(DEBUG,"Enter.")
#endif
#ifdef ADCNETCDF
      if (iret .ne. NF90_NOERR) then
         call allMessage(ERROR,nf90_strerror(iret))
         call ITFRICTERMINATE()
      endif
#endif
#ifdef ALL_TRACE
      call allMessage(DEBUG,"Return.")
#endif
      !call unsetMessageSource()
C-----------------------------------------------------------------------
      end subroutine check_err
C-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     S U B R O U T I N E   I T F R I C T E R M I N A T E
!-----------------------------------------------------------------------
      SUBROUTINE ITFRICTERMINATE(NO_MPI_FINALIZE)
      USE SIZES, ONLY : myproc
#ifdef CMPI 
      USE MESSENGER
#endif
      IMPLICIT NONE 
      LOGICAL, OPTIONAL :: NO_MPI_FINALIZE
#ifdef ALL_TRACE
      REAL, ALLOCATABLE :: dummy(:)
#endif
C
      !call setMessageSource("BC3DTerminate")
#ifdef ALL_TRACE
      call allMessage(DEBUG,"Enter.")
#endif

      call allMessage(INFO,"ADCIRC Terminating.")

#ifdef ALL_TRACE
      ! intentionally create a segmentation fault so that we can get
      ! a stack trace to determine the line number of the netcdf call
      ! that went bad ... this assumes that the code was compiled with
      ! debugging symbols, bounds checking, and stack trace turned on.
      write(6,*) 'myproc = ',myProc      
      dummy(1) = 99.9d0
#endif

#ifdef CMPI
      subdomainFatalError = .true.
      IF (PRESENT(NO_MPI_FINALIZE)) THEN
        CALL MSG_FINI(NO_MPI_FINALIZE)
      ELSE
        CALL MSG_FINI()
      ENDIF
#endif
      CALL EXIT(1)
C
#ifdef ALL_TRACE
      call allMessage(DEBUG,"Return.") ! should be unreachable
#endif
      !call unsetMessageSource()
C---------------------------------------------------------------------C
      END SUBROUTINE ITFRICTERMINATE
C---------------------------------------------------------------------C
C---------------------------------------------------------------------C
      END MODULE INTERNALTIDE
C---------------------------------------------------------------------C
