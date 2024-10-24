!******************************************************
!*                                                    *
!* This module handles I/O of the program including   *
!* definition of inputs/outputs, reading inputs from  *
!* a data file, echo inputs, and write outputs        *
!* to data files. To interface with outside           *
!* environment, one needs to provide such             *
!* capability from the outside environment.           *
!*                                                    *
!*             Copyright (c) Wenbin Yu                *
!*         All rights reserved                        *
!*     Mechanical and Aerospace Engineering           *
!*                                                    *
!*           Utah State University                    *
!******************************************************

MODULE IO

USE GlobalDataFun
USE PrescribedCondition
USE TimeFunctionModule
IMPLICIT NONE

PRIVATE                              ! So everything is private, except declared by public

! functions called by the main program
!---------------------------------------------
PUBLIC Input,WriteError,Output

! variables needed by the main program
!----------------------------------------------
PUBLIC error,nkp,nelem,ndof_el,coord,member,pt_condition,material,sol_pt,sol_mb, &
     & frame,mb_condition,distr_fun,curvature,omega_a0,omega_a_tf,v_root_a0,v_root_a_tf,time_function,simu_time, &
	 & EIN,niter,nstep,ncond_mb,ntimefun,analysis_flag,init_cond,nev,eigen_val,eigen_vec_pt,eigen_vec_mb,&
	 & nmemb,ncond_pt,nmate,nframe,ndistrfun,ncurv

! files needed/generated
!--------------------------------------------------------------------------------
INTEGER,PARAMETER,PRIVATE:: CHAR_LEN=256
INTEGER,PARAMETER:: IN  =10     ! input file: inp_name
CHARACTER(CHAR_LEN)    :: inp_name

INTEGER,PARAMETER:: EIN =20     ! file for echoing the inputs: inp_name.ech
CHARACTER(CHAR_LEN+3)    :: ech_name

INTEGER,PARAMETER:: OUT =40     ! file for output: inp_name.out
CHARACTER(CHAR_LEN+3)    :: out_name
!--------------------------------------------------------------------------------
INTEGER,PARAMETER:: INIT =50     ! file for initial conditions: inp_name.ini
CHARACTER(CHAR_LEN+3)    :: init_name
!--------------------------------------------------------------------------------


!Private constants
!-----------------------------------------------


! Private variables
!---------------------------------------------------------------------------
INTEGER::nkp           ! number of key points
INTEGER::nelem         ! total number of elements
INTEGER::nmemb         ! number of members
INTEGER::nmate         ! number of cross-sectional properties sets
INTEGER::nframe        ! number of frames
INTEGER::ncond_pt      ! number of point conditions for concentrated loads and boundary conditions
INTEGER::ndistrfun     ! number of distributed functions
INTEGER::ncurv         ! number of initial curvatures/twists
INTEGER::analysis_flag ! 0: static analysis; 1: steady state response; 2: transient analysis; 3: eigenvalue analysis
INTEGER::nev           ! number of frequencies and modeshapes.


!Public integer variables
!---------------------------------------
INTEGER::ncond_mb      ! number of member conditions for distributed loads
INTEGER::ntimefun      ! number of time functions
INTEGER::niter         ! number of maximum iterations
INTEGER::nstep         ! number of time steps/load steps
INTEGER,ALLOCATABLE::member(:,:)   ! member property array: member(nmemb,MEMB_CONST)
INTEGER::ndof_el       ! dofs per element: 12 for static analysis, 18 for dynamic analysis
INTEGER::omega_a_tf(NDIM)    ! time function numbers for the angular velocity of frame a
INTEGER::v_root_a_tf(NDIM)   ! time function numbers for the velocity of the starting point of the first member


!Public real variables
!--------------------------------------------------------------------------------
REAL(DBL),ALLOCATABLE:: coord(:,:)       ! nodal coordinates: coord(nkp,NDIM)
REAL(DBL),ALLOCATABLE:: material(:,:,:)  ! flexibility matrix: (nmate,12,6)
REAL(DBL),ALLOCATABLE:: frame(:,:,:)     ! member frames: (nframe,3,3)
REAL(DBL),ALLOCATABLE:: distr_fun(:,:)   ! prescribed functions: (ndistrfun,6)
REAL(DBL),ALLOCATABLE:: curvature(:,:)   ! curvatures: (ncurv,NDIM)
REAL(DBL),ALLOCATABLE:: sol_pt(:,:,:)    ! solutions for points sol_pt(nstep,nkp,NDIM+NDOF_ND)
REAL(DBL),ALLOCATABLE:: sol_mb(:,:,:)    ! solutions for member sol_mb(nstep,nelem,NDIM+ndof_el): nelem: total number of elements
REAL(DBL)            :: simu_time(2)     ! start and end time of the simulation.
REAL(DBL)            :: omega_a0(NDIM)   ! the magnitude of angular velocity of frame a
REAL(DBL)            :: v_root_a0(NDIM)  ! the magnitude of linear velocity of the starting point of the first member
REAL(DBL),ALLOCATABLE:: init_cond(:,:)   ! initial conditions: init_cond(nelem,12);
                                         ! init_cond(nelem,1:6) for initial displacements/rotations
									     ! init_cond(nelem,7:12) for initial velocities
REAL(DBL),ALLOCATABLE::eigen_val(:,:),eigen_vec_pt(:,:,:),eigen_vec_mb(:,:,:) ! arrays for holding eigenvalues and eigenvectors

!Public derived types
!--------------------------------------------------------------------------------
TYPE(PrescriInf),ALLOCATABLE::pt_condition(:) ! prescribed information concentrated at nodes
TYPE(PrescriInf),ALLOCATABLE::mb_condition(:) ! prescribed information distributed along beam members
TYPE(TimeFunction),ALLOCATABLE::time_function(:) ! time functions


!Public character variables
!============================================================================
CHARACTER(300)::error         ! a character variable holding  error message
!=========================================================================================

CONTAINS
!=============================


!*************************************************************
!*                                                           *
!* To read and echo problem definition data                  *
!*															 *
!*************************************************************
SUBROUTINE Input

INTEGER:: i,j
INTEGER::tmp_no ! a temporary integer

!Get the input file name from the command line
!------------------------------------------------------------------------------
CALL GETARG(1,inp_name)
IF(TRIM(inp_name)=='') THEN
 error='Please provide an input file name, executing as GEBT input_file_name'
 RETURN
ENDIF
IF(FileOpen(IN, inp_name, 'OLD', 'READ',error))	 RETURN


! Create a file name for echoing the input data
!--------------------------------------------------
ech_name=TRIM(inp_name) // ".ech"
IF(FileOpen(EIN,  ech_name,'REPLACE','WRITE',error)) RETURN


! Input and echo analysis control parameters.
!---------------------------------------------------------
READ(IN,*,IOSTAT=in_stat) analysis_flag,niter,nstep
IF(IOError('read analysis conontrol parameters',error)) RETURN

CALL TitlePrint(EIN, 'Analysis Control Parameters')
WRITE(EIN,*) "Analysis Type                = ", analysis_flag
WRITE(EIN,*) "Number of Maximum Iterations = ", niter
WRITE(EIN,*) "Number of Time/Load steps    = ", nstep

IF(analysis_flag==0) THEN
	ndof_el=NDOF_ND
ELSE
    ndof_el=18
	READ(IN,*,IOSTAT=in_stat) omega_a0
	IF(IOError('read angular velocity of the global body-attached frame',error)) RETURN
    WRITE(EIN,*)"Angular Velocity of a frame"
    WRITE(EIN,*) '--------------------------------'
	CALL WriteVec(EIN,omega_a0)

	READ(IN,*,IOSTAT=in_stat) omega_a_tf
	IF(IOError('read time functions for angular velocity of the global body-attached frame',error)) RETURN
    WRITE(EIN,*)"Time Function for Angular Velocity of a frame"
    WRITE(EIN,*) '--------------------------------'
	CALL WriteVec(EIN,omega_a_tf)


	READ(IN,*,IOSTAT=in_stat) v_root_a0
	IF(IOError('read velocities of starting point of the first member',error)) RETURN
	WRITE(EIN,*)"Initial Velocity of the starting point of the first member"
    WRITE(EIN,*) '--------------------------------'
    CALL WriteVec(EIN,v_root_a0)

	READ(IN,*,IOSTAT=in_stat) v_root_a_tf
	IF(IOError('read time function for velocities of starting point of the first member',error)) RETURN
	WRITE(EIN,*)"Time Function for Initial Velocity of the starting point of the first member"
    WRITE(EIN,*) '--------------------------------'
    CALL WriteVec(EIN,v_root_a_tf)

	IF(analysis_flag==3) THEN
		READ(IN,*,IOSTAT=in_stat) nev
		IF(IOError('read number of frequencies/mode shapes',error)) RETURN
		WRITE(EIN,*) "Number of Frequencies/Mode Shapes= ", nev
	ENDIF
ENDIF


! Input and echo mesh control parameters
!----------------------------------------------------------------
READ(IN,*,IOSTAT=in_stat)  nkp, nmemb, ncond_pt,nmate, nframe,ncond_mb,ndistrfun,ntimefun, &
                        &   ncurv
IF(IOError('read mesh conontrol parameters',error)) RETURN

CALL TitlePrint(EIN, 'Mesh Control Parameters')
WRITE(EIN,*) "Number of Key Points             = ", nkp
WRITE(EIN,*) "Number of Members                = ", nmemb
WRITE(EIN,*) "Number of Point Conditions       = ", ncond_pt
WRITE(EIN,*) "Number of Cross-sections         = ", nmate
WRITE(EIN,*) "Number of Frames                 = ", nframe
WRITE(EIN,*) "Number of Distributed Load Cases = ", ncond_mb
WRITE(EIN,*) "Number of Time Functions         = ", ntimefun
WRITE(EIN,*) "Number of Distributed Functions  = ", ndistrfun
WRITE(EIN,*) "Number of Curvatures/Twist       = ", ncurv

! A small check for the control data
!---------------------------------
IF(nkp<=0 .OR. nmemb<=0 .OR. nmate<=0 .OR. ncond_pt<=0) THEN
   error='nkp, nmemb, nmat, ncond_pt must be greater than 0'
   RETURN
ENDIF


! Input and echo coordinates for key points
!=============================================================================================
ALLOCATE(coord(nkp,NDIM),STAT=allo_stat)
IF(MemoryError('coord',error)) GOTO 9999
coord=0.0D0

CALL TitlePrint(EIN, 'Key Point Coordinates')

DO i=1,nkp
   READ(IN,*,IOSTAT=in_stat)tmp_no,coord(tmp_no,:)
   IF(IOError('read key point coordinates',error)) GOTO 9999

   WRITE(EIN,*)"Point NO: ",tmp_no
   WRITE(EIN,*) '--------------------------------'
   CALL WriteVec(EIN,coord(tmp_no,:))
ENDDO


! Input and echo member properties: we need 7 numbers to describe a member including starting
! and ending points, structural/inertial property set # for the starting and ending points,
! frame # of the starting point, number of elements divided, geometry property set #
!=============================================================================================
ALLOCATE(member(nmemb,MEMB_CONST),STAT=allo_stat)
IF(MemoryError('member',error)) GOTO 9999
member=0

CALL TitlePrint(EIN, 'Member Definition')
WRITE(EIN,*) '      N1        N2     Sec#1  Sec#2  Frm#   #Elems   Geom#'
WRITE(EIN,*) '----------------------------------------------------------------'

DO i=1,nmemb
   READ(IN,*,IOSTAT=in_stat)tmp_no,member(tmp_no,:)
   IF(IOError('read member properties',error)) GOTO 9999
   WRITE(EIN,*)"Member NO: ",tmp_no
   WRITE(EIN,*) '--------------------------------'
   CALL WriteVec(EIN,member(tmp_no,:))
ENDDO

! Input and echo point conditions
!=============================================================================================
ALLOCATE(pt_condition(ncond_pt),STAT=allo_stat)
IF(MemoryError('pt_condition',error)) GOTO 9999
pt_condition=InitPI()

CALL TitlePrint(EIN, 'Prescribed Point Conditions')
DO i=1,ncond_pt
	pt_condition(i)=InputEchoPrescribedConditions(IN,EIN,error)
	IF(error/='')GOTO 9999
ENDDO


! Input and echo cross-sectional properties including flexibility matrix and mass matrix
!=============================================================================================
ALLOCATE(material(nmate,ndof_el-NSTRN,NSTRN),STAT=allo_stat)
IF(MemoryError('material',error)) GOTO 9999
material=0.0D0

CALL TitlePrint(EIN, 'Sectional Properties')

DO i=1,nmate

   READ(IN,*,IOSTAT=in_stat) tmp_no
   IF(IOError('read number of sectional properties',error)) GOTO 9999

   WRITE(EIN,*) 'Section No.         =', tmp_no
   WRITE(EIN,*) '--------------------------------'
   CALL TitlePrint(EIN, 'Sectional Flexibility Matrix')

   DO j=1,NSTRN
      READ(IN,*,IOSTAT=in_stat)material(tmp_no,j, :)
	  IF(IOError('read sectional flexibility matrix',error)) GOTO 9999
	  CALL WriteVec(EIN,material(tmp_no,j,:))
   ENDDO

   IF(analysis_flag/=0) THEN
   		CALL TitlePrint(EIN, 'Sectional Mass Matrix')

		DO j=NSTRN+1,NSTRN+NSTRN
			READ(IN,*,IOSTAT=in_stat)material(tmp_no,j, :)
			IF(IOError('read sectional mass matrix',error)) GOTO 9999
			CALL WriteVec(EIN,material(tmp_no,j,:))
		ENDDO
    ENDIF

ENDDO


! Input and echo frame
!=============================================================================================
IF(nframe>0) THEN
	ALLOCATE(frame(nframe,NDIM,NDIM),STAT=allo_stat)
	IF(MemoryError('frame',error)) GOTO 9999
	frame=0.0D0

	CALL TitlePrint(EIN, 'Member Frames')

	DO i=1,nframe

		READ(IN,*,IOSTAT=in_stat) tmp_no
		IF(IOError('read number of member frames',error)) GOTO 9999

		WRITE(EIN,*) 'Frame No.         =', tmp_no
		WRITE(EIN,*) '-----------------------------------'

		DO j=1,NDIM
			READ(IN,*,IOSTAT=in_stat)frame(tmp_no,j, :)
			IF(IOError('read member frames',error)) GOTO 9999
			CALL WriteVec(EIN,frame(tmp_no,j,:)) ! echo
		ENDDO

	ENDDO
ENDIF

! Input and echo distributed loads
!=============================================================================================
IF (ncond_mb>0) THEN
	ALLOCATE(mb_condition(ncond_mb),STAT=allo_stat)
	IF(MemoryError('mb_condition',error)) GOTO 9999
	mb_condition=InitPI()

	DO i=1,ncond_mb
		CALL TitlePrint(EIN, 'Distributed Loads')
		mb_condition(i)=InputEchoPrescribedConditions(IN,EIN,error)
		IF(error/='')GOTO 9999
	ENDDO
ENDIF

! Input and echo functions approximated using Chebychev polynomials up to the fifth order
!===========================================================================================
IF(ndistrfun>0) THEN
	ALLOCATE(distr_fun(ndistrfun,NSTRN),STAT=allo_stat)
	IF(MemoryError('distr_fun',error)) GOTO 9999
	distr_fun=0.0D0

    CALL TitlePrint(EIN, 'Distributed Functions Approximated Using Chebychev polynomials')

	DO i=1,ndistrfun
		READ(IN,*,IOSTAT=in_stat) tmp_no
	    IF(IOError('read function number for Chebychev polynomials',error)) GOTO 9999

		READ(IN,*,IOSTAT=in_stat)distr_fun(tmp_no, :)
		IF(IOError('read functions for distributed load',error)) GOTO 9999

		WRITE(EIN,*) 'Function No.=', tmp_no
		WRITE(EIN,*) '------------------------------------'
		CALL WriteVec(EIN,distr_fun(tmp_no,:))
    ENDDO

ENDIF

! Input and echo initial curvatures/twist
!=============================================================================================
IF (ncurv>0) THEN
	ALLOCATE(curvature(ncurv,NDIM),STAT=allo_stat)
	IF(MemoryError('curvature',error)) GOTO 9999
	curvature=0.0D0

    CALL TitlePrint(EIN, 'Initial Curvatures/Twist')

	DO i=1,ncurv

		READ(IN,*,IOSTAT=in_stat) tmp_no
		IF(IOError('read # for initial curvatures/twist',error)) GOTO 9999

	    READ(IN,*,IOSTAT=in_stat)curvature(tmp_no, :)
	    IF(IOError('read initial curvatures/twist',error)) GOTO 9999

		WRITE(EIN,*) 'Case No.=', tmp_no
		WRITE(EIN,*) '---------------------------------------'
		CALL WriteVec(EIN,curvature(tmp_no,:))
   ENDDO

ENDIF


! Input and echo time functions
!=============================================================================================
IF (ntimefun>0.OR.analysis_flag==2) THEN
	WRITE(EIN,*)

	READ(IN,*,IOSTAT=in_stat)simu_time
	IF(IOError('read simulation range',error)) GOTO 9999
	WRITE(EIN,'(A13,'//FMT_REAL//')') "Start Time =", simu_time(1)
    WRITE(EIN,'(A13,'//FMT_REAL//')') "End Time   =", simu_time(2)
ENDIF

IF (ntimefun>0) THEN
	ALLOCATE(time_function(ntimefun),STAT=allo_stat)
	IF(MemoryError('time_function',error)) GOTO 9999
	time_function=InitTF()

	DO i=1,ntimefun

		READ(IN,*,IOSTAT=in_stat) tmp_no
		IF(IOError('read number of time functions',error)) GOTO 9999
		IF(tmp_no<1) THEN
			error="time function number cannot be less than 1"
			GOTO 9999
		ENDIF
		WRITE(EIN,*) 'Time Function No.         =', tmp_no
		WRITE(EIN,*) '-----------------------------------'
		time_function(tmp_no)=InputEchoTimeFunctions(IN,EIN,error)
		IF(error/='')GOTO 9999

	ENDDO

ENDIF

! For time marching, we also need to read initial conditions, could be provided by
! users or automatically generatedly from steady state solution
!===============================================================
nelem=SUM(member(:,6))  ! total number of elements for the structure
IF(analysis_flag==2) THEN
	ALLOCATE(init_cond(nelem,12),STAT=allo_stat)  ! storing the initial coniditions
	IF(MemoryError('init_cond',error)) GOTO 9999

	! Create a file name for reading initial data
	!--------------------------------------------------
	init_name=TRIM(inp_name) // ".ini"
	IF(FileOpen(INIT, init_name, 'OLD', 'READ',error))	 RETURN

	WRITE(EIN,*)
	WRITE(EIN,*) 'Initial conditions'
	WRITE(EIN,*) '-----------------------------------'

	DO i=1,nelem
		READ(INIT,*,IOSTAT=in_stat)init_cond(i, 1:6)
		IF(IOError('read initial displacements/rotations',error)) GOTO 9999
		CALL WriteVec(EIN,init_cond(i, 1:NSTRN))
	ENDDO

	DO i=1,nelem
		READ(INIT,*,IOSTAT=in_stat)init_cond(i, 7:12)
		IF(IOError('read initial velocities',error)) GOTO 9999
        CALL WriteVec(EIN,init_cond(i, 7:12))
	ENDDO

	CLOSE(INIT)
ENDIF


! Creat arrays holding the outputs
!===============================================================
ALLOCATE(sol_pt(nstep,nkp,NDIM+NDOF_ND),STAT=allo_stat)  ! for each point there are 12 variables
IF(MemoryError('sol_pt',error)) GOTO 9999
ALLOCATE(sol_mb(nstep,nelem,NDIM+ndof_el),STAT=allo_stat) ! for each element, there are 18 variables.
IF(MemoryError('sol_mb',error)) GOTO 9999

IF(analysis_flag==3) THEN
	ALLOCATE(eigen_val(2,nev+1),STAT=allo_stat)
	IF(MemoryError('eigen_val',error)) GOTO 9999
	ALLOCATE(eigen_vec_pt(nev+1,nkp,NDIM+NDOF_ND),STAT=allo_stat)
	IF(MemoryError('eigen_vec_pt',error)) GOTO 9999
	ALLOCATE(eigen_vec_mb(nev+1,nelem,NDIM+ndof_el),STAT=allo_stat)
	IF(MemoryError('eigen_vec_mb',error)) GOTO 9999
ENDIF


WRITE(*,*) 'The inputs are echoed in ',TRIM(ech_name)
CLOSE(IN)

9999 IF(error/='')THEN
		IF(ALLOCATED(eigen_vec_mb))  DEALLOCATE(eigen_vec_mb)
		IF(ALLOCATED(eigen_vec_pt))  DEALLOCATE(eigen_vec_pt)
		IF(ALLOCATED(eigen_val)) 	 DEALLOCATE(eigen_val)
	    IF(ALLOCATED(sol_mb))		 DEALLOCATE(sol_mb)
    	IF(ALLOCATED(sol_pt))	     DEALLOCATE(sol_pt)
		IF(ALLOCATED(init_cond))	 DEALLOCATE(init_cond)
		IF(ALLOCATED(time_function)) DEALLOCATE(time_function)
		IF(ALLOCATED(curvature))	 DEALLOCATE(curvature)
		IF(ALLOCATED(distr_fun))     DEALLOCATE(distr_fun)
		IF(ALLOCATED(mb_condition))  DEALLOCATE(mb_condition)
		IF(ALLOCATED(frame)) 		 DEALLOCATE(frame)
    	IF(ALLOCATED(material))		 DEALLOCATE(material)
		IF(ALLOCATED(pt_condition))  DEALLOCATE(pt_condition)
		IF(ALLOCATED(member))		 DEALLOCATE(member)
		IF(ALLOCATED(coord)) 		 DEALLOCATE(coord)
     ENDIF

END SUBROUTINE Input
!*****************************************



!****************************************************************
!*                                                              *
!* To output the results for each point and each division       *
!* within each member                                           *
!*==============================================================*
!* x y z                                                        *
!* ux uy uz thetax thetay thetaz                                *
!* Fx Fy Fz Mx My Mz                                            *
!* x y z are the coordinates in the global coordinate system    *
!* for element the coordinate for the mid-point is used         *
!*															    *
!****************************************************************
SUBROUTINE Output

INTEGER::istep,ikp,imemb,j

INTEGER::ndiv ! number of divisions
INTEGER::div_no ! current division number


out_name=TRIM(inp_name) // ".out"
IF(FileOpen(OUT,  out_name,'REPLACE','WRITE',error)) RETURN

CALL TitlePrint(OUT, 'The Solution of Internal Variables')

DO istep=1,nstep

   IF(nstep/=1) WRITE(OUT,*)"Step #", istep

   DO ikp=1,nkp
	  WRITE(OUT,*)"Point #: ",ikp
	  WRITE(OUT,*) '--------------------------------'
      CALL WriteVec(OUT,sol_pt(istep,ikp,1:3))
      CALL WriteVec(OUT,sol_pt(istep,ikp,4:9))
      CALL WriteVec(OUT,sol_pt(istep,ikp,10:15))
	  WRITE(OUT,*)
   ENDDO

   div_no=0

   DO imemb=1,nmemb
		WRITE(OUT,*)"Member #: ",imemb
		WRITE(OUT,*) '--------------------------------'

		ndiv=member(imemb,6)

		DO j=1,ndiv
			div_no=div_no+1
			CALL WriteVec(OUT,sol_mb(istep,div_no,1:3))
			CALL WriteVec(OUT,sol_mb(istep,div_no,4:9))
			CALL WriteVec(OUT,sol_mb(istep,div_no,10:15))
			IF(ndof_el==18) CALL WriteVec(OUT,sol_mb(istep,div_no,16:21))
			WRITE(OUT,*)
		ENDDO
	ENDDO
	WRITE(OUT,*)
ENDDO

!DO istep=1,nstep
!    write(31,'(1x,6es25.14)')sol_pt(istep,1,4:9)
!    WRITE(32,'(1x,6es25.14)')sol_pt(istep,2,4:9)

    !write(30,'(1x,6es25.14)')sol_pt(istep,2,10:15)

!ENDDO

IF(analysis_flag==3) THEN

	DO istep=1,nev

		WRITE(OUT,*)"Eigenvalue #", istep
		CALL WriteVec(OUT,eigen_val(:,istep))

		DO ikp=1,nkp
			WRITE(OUT,*)"Point #: ",ikp
			WRITE(OUT,*) '--------------------------------'
			CALL WriteVec(OUT,eigen_vec_pt(istep,ikp,1:3))
			CALL WriteVec(OUT,eigen_vec_pt(istep,ikp,4:9))
			CALL WriteVec(OUT,eigen_vec_pt(istep,ikp,10:15))
			WRITE(OUT,*)
		ENDDO

		div_no=0

		DO imemb=1,nmemb
			WRITE(OUT,*)"Member #: ",imemb
			WRITE(OUT,*) '--------------------------------'

			ndiv=member(imemb,6)

			DO j=1,ndiv
				div_no=div_no+1
				CALL WriteVec(OUT,eigen_vec_mb(istep,div_no,1:3))
				CALL WriteVec(OUT,eigen_vec_mb(istep,div_no,4:9))
				CALL WriteVec(OUT,eigen_vec_mb(istep,div_no,10:15))
				IF(ndof_el==18) CALL WriteVec(OUT,eigen_vec_mb(istep,div_no,16:21))
				WRITE(OUT,*)
			ENDDO
		ENDDO
	ENDDO
ENDIF

WRITE(*,*) 'The results can be found in ',TRIM(out_name)
CLOSE(OUT)

! Construct given initial conditions
! Assume the initial displacement/rotations is taken from the last step of the
! steady state solution and the initial velocities are zero, how these values are given
! should be determined by the end user.
!---------------------------------------------------------------------
IF(analysis_flag==2) THEN
   init_name=TRIM(inp_name) // ".ini_out"

   IF(FileOpen(INIT, init_name, 'REPLACE','WRITE',error))	 RETURN

   !WRITE(INIT,'(1x,6ES25.15)')(sol_mb(nstep,j,4:9),j=1,div_no)
   !WRITE(INIT,'(1x,6ES25.15)')(sol_mb(nstep,j,10:15),j=1,div_no)

   WRITE(INIT,'(1x,6ES25.15)')(init_cond(j,1:6) ,j=1,div_no)
   WRITE(INIT,'(1x,6ES25.15)')(init_cond(j,7:12),j=1,div_no)

   !WRITE(INIT,'(1x,6ES25.15)') ((0.0D0,i=1,NSTRN),j=1,div_no)
   CLOSE(INIT)
ENDIF


END SUBROUTINE Output
!*****************************************


END MODULE IO
!=========================================================
