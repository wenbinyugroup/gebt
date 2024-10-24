
!***************************************************************
!*                                                             *
!* The file contains the analysis subroutine called by         *
!* the main program. This subroutine can be called directly by *
!* an outside environment where the arguments are defined and  *
!* passed to and from this subroutine                          *
!*                                                             *
!*                   Copyright (c) Wenbin Yu                   *
!*                  All rights reserved                        *
!*              Mechanical and Aerospace Engineering           *
!*                                                             *
!*                    Utah State University                    *
!*                                                             *
!***************************************************************
SUBROUTINE  Analysis(nkp,nelem,ndof_el,nmemb,ncond_pt,nmate, nframe,ndistrfun,ncurv,coord, &
                 &  member,pt_condition,material,niter,nstep,sol_pt,sol_mb,error,ncond_mb, &
                 &  ntimefun,frame,mb_condition,distr_fun,curvature,omega_a0,omega_a_tf, &
                 &  v_root_a0,v_root_a_tf,simu_time,time_function,analysis_flag,init_cond, &
                 &  nev,eigen_val,eigen_vec_pt,eigen_vec_mb)

   	
!DEC$ ATTRIBUTES DLLEXPORT ::  Analysis
USE InternalData
USE PreproModule
USE PrescribedCondition
USE TimeFunctionModule
USE Solve	
USE Eigen	
	
	
IMPLICIT NONE

INTEGER,INTENT(IN)         ::nkp ! number of key points
INTEGER,INTENT(IN)         ::nelem ! number of elements
INTEGER,INTENT(IN)         ::ndof_el            ! number of unknowns for each element
INTEGER,INTENT(IN)         ::nmemb ! number of beam members
INTEGER,INTENT(IN)         ::ncond_pt ! number of point conditions
INTEGER,INTENT(IN)         ::nmate ! number of different cross-sections
INTEGER,INTENT(IN)         ::nframe ! number of different cross-sectional frames
INTEGER,INTENT(IN)         ::ndistrfun ! number of distribution functions
INTEGER,INTENT(IN)         ::ncurv     ! number of sets of initial curvature/twists
REAL(DBL),INTENT(IN)       ::coord(nkp,NDIM)         ! coordinates for key points
INTEGER,INTENT(IN)         ::member(nmemb,MEMB_CONST)        ! member desriptions
TYPE(PrescriInf),INTENT(INOUT)::pt_condition(ncond_pt) ! prescribed concentrated infromation at points
REAL(DBL),INTENT(INOUT)       ::material(nmate,ndof_el-NSTRN,NSTRN)    ! cross-sectional properties
INTEGER,INTENT(IN)         ::niter              ! maximum number of iterations
INTEGER,INTENT(IN)         ::nstep              ! number of time steps
REAL(DBL),INTENT(OUT)      ::sol_pt(nstep,nkp,NDIM+NDOF_ND)      ! solutions for points
REAL(DBL),INTENT(OUT)      ::sol_mb(nstep,nelem,NDIM+ndof_el)      ! solutions for members
CHARACTER(*),INTENT(OUT)   ::error              ! error message
	
INTEGER,INTENT(IN)::ncond_mb,ntimefun	        ! number of prescribed distribution, number of time function                			                		
REAL(DBL),INTENT(IN)::frame(nframe,NDIM,NDIM)              ! the direction cosine matrix for each frame. The default frame is the global frame 
TYPE(PrescriInf),INTENT(INOUT)::mb_condition(ncond_mb) ! distributed load information
REAL(DBL),INTENT(IN)::distr_fun(ndistrfun,NSTRN)            ! distributed functions 
REAL(DBL),INTENT(IN)::curvature(ncurv,NDIM)            ! curvature information
REAL(DBL),INTENT(IN)::omega_a0(NDIM)                ! angular velocity of frame a
REAL(DBL),INTENT(IN)::v_root_a0(NDIM)               ! linear velocity of starting point of the first member
INTEGER,  INTENT(IN)::omega_a_tf(NDIM)             ! time function for angular velocity of frame a
INTEGER,  INTENT(IN)::v_root_a_tf(NDIM)            ! time function for linear velocity of starting point of the first member

REAL(DBL),INTENT(IN)::simu_time(2)              ! start and end time
TYPE(TimeFunction),INTENT(IN)::time_function(ntimefun) ! time functions
REAL(DBL),INTENT(INOUT)::init_cond(nelem,12)
INTEGER,INTENT(IN)     :: analysis_flag
INTEGER,INTENT(INOUT):: nev  ! upon return, nev is the number of converged eigen values might be the original nev given by the user in the inputs or nev+1 due to complex conjugate pairs
REAL(DBL),INTENT(OUT)::eigen_val(2,nev+1),eigen_vec_pt(nev+1,nkp,NDIM+NDOF_ND),eigen_vec_mb(nev+1,nelem,NDIM+ndof_el)  

! Variables used in this subroutine
!---------------------------------------------------------------------------------------	
INTEGER:: dof_con(nkp)                          ! the connecting condition for key point.
                                                ! 0, a boundary point; 1, a connection point.
REAL(DBL),ALLOCATABLE:: x(:)                    ! the solution vector of the linear system

REAL(DBL)            ::PHDot(nelem,6)

REAL(DBL),ALLOCATABLE:: eigen_vec(:,:)
  
    
INTEGER::i,j
REAL(DBL)::time_incr,time_current


REAL(DBL):: omega_a(NDIM),v_root_a(NDIM)

TYPE (MemberInf)::memb_info(SIZE(member,1))
REAL(DBL)::t0

! Create a file name for debugging data
!------------------------------------------
IF(DEBUG) THEN
  deb_name="gebt.deb" 
  IF(FileOpen(IOUT, deb_name, 'UNKNOWN', 'WRITE',error))	 RETURN 
ENDIF  

CALL Preprocess(nkp,nelem,ndof_el,member,material,frame,coord,curvature,&
                    & dof_con,memb_info,error)
IF(error/='') RETURN

ALLOCATE(x(nsize),STAT=allo_stat)
IF(MemoryError('x',error)) GOTO 9999
x=0.0D0


ALLOCATE(dof_all(nkp,NSTRN),STAT=allo_stat)
IF(MemoryError('dof_all',error)) GOTO 9999
dof_all=0

ALLOCATE(follower_all(nkp,NSTRN),STAT=allo_stat)
IF(MemoryError('cond_all',error)) GOTO 9999
follower_all=0

ALLOCATE(cond_all(nkp,NSTRN),STAT=allo_stat)
IF(MemoryError('cond_all',error)) GOTO 9999
cond_all=0.0D0

sol_pt=0.0D0;sol_mb=0.0D0 ! initialize

! Obtain the prescribed dof and follower condition for each point
!-----------------------------------------------------------------
CALL GetPrescribedDOF(nkp,pt_condition,dof_all,follower_all)

IF(analysis_flag==3)THEN

	! eigen vector for
	!-------------------------------------------------
	ALLOCATE(eigen_vec(nsize,nev+1),STAT=allo_stat)
	IF(MemoryError('eigen_vec',error)) GOTO 9999
	eigen_vec=0.0d0;eigen_vec_pt=0.0d0;eigen_vec_mb=0.0d0

ENDIF


! Initial angular velocity and linear velocities of the a frame
!-----------------------------------------------------------------------
!omega_a=omega_a0;v_root_a=v_root_a0
t0=0.0d0
omega_a =CurrentValues(omega_a0,omega_a_tf,time_function,t0)
v_root_a=CurrentValues(v_root_a0,v_root_a_tf,time_function,t0) 

IF(ntimefun>0.OR.analysis_flag==2)  time_incr=(simu_time(2)-simu_time(1))/nstep

init_flag=0  ! assign default value to init_flag

IF(analysis_flag==2) THEN

	two_divide_dt=2.0D0/time_incr !2/dt

	!Obtaining the complete set of initial conditions needed for 
	!time marching later
	!----------------------------------------
	init_flag=1
	!WRITE(*,*)"THE INITIAL STEP"
	IF(DEBUG) WRITE(IOUT,*)"THE INITIAL STEP"


	IF(niter==1) THEN  ! solve the linear system for initial conditions
		CALL LinearSolution(ndof_el,memb_info,v_root_a,omega_a,member,error,&
	    	& ncond_mb,mb_condition,distr_fun,(dof_con),x,init_cond)
	ELSE ! Solve the nonlinear system for initial conditions
	!---------------------------------------------------
		CALL NewtonRaphson(ndof_el,memb_info,v_root_a,omega_a,member,niter,error,&
	    	& ncond_mb,mb_condition,distr_fun,(dof_con),x,init_cond)
	ENDIF
	IF(error/='') GOTO 9999

	CALL ExtractElementValues(ndof_el,member,x,sol_mb(1,:,:))
	PHDot=sol_mb(1,:,4:9)   ! CTCabPDot, CTCabHDot
	sol_mb(1,:,4:9)=init_cond(:,1:6) ! insert the initial coniditions for ui, thetai into the solution.


	! set initial values as initial condition of x for time marching, replace Pdot, Hdot with u, theta
	!----------------------------------------------------------------------------------------------
	CALL InsertElementValues(ndof_el,member,x,init_cond)

	! obtain: 2/dt*value+\dot{value} at time t0
	!-------------------------------------------------------
	init_cond(:,1:6) =two_divide_dt*sol_mb(1,:,4:9)+init_cond(:,7:12) ! u, \theta
	init_cond(:,7:12)=two_divide_dt*CTCabPH(niter,member,memb_info,sol_mb(1,:,4:))+PHDot ! P, H

	init_flag=2

ENDIF


DO i=1,nstep
    ! evaluate the prescribed information for this step
	!---------------------------------------------------
	IF(ntimefun>0)THEN    
		time_current=time_incr*i+simu_time(1)
		
		CALL UpdatePI(pt_condition,time_function,time_current)
		IF(ncond_mb>0)	CALL UpdatePI(mb_condition,time_function,time_current)
	  
! updating the time varying angular/linear velocities of a frame
!------------------------------------------------------------------------
	    omega_a =CurrentValues(omega_a0,omega_a_tf,time_function,time_current)
	    v_root_a=CurrentValues(v_root_a0,v_root_a_tf,time_function,time_current)

   	ENDIF	


	!obtain the current prescribed value
	!note without time functions, the load steps are not meaningful
	!------------------------------------------------------------
    CALL GetPrescribedVal(nkp,pt_condition,cond_all)
	
	IF(nstep/=1) THEN
		WRITE(*,*) "STEP=",i
		IF(DEBUG) WRITE(IOUT,*) "STEP=",i
	ENDIF


    IF(niter==1) THEN  ! solve the linear system for ith step
      CALL LinearSolution(ndof_el,memb_info,v_root_a,omega_a,member,error,&
			& ncond_mb,mb_condition,distr_fun,(dof_con),x,init_cond)
	ELSE ! Solve the nonlinear system for ith step
	!---------------------------------------------------
		CALL NewtonRaphson(ndof_el,memb_info,v_root_a,omega_a,member,niter,error,&
			& ncond_mb,mb_condition,distr_fun,(dof_con),x,init_cond)
	ENDIF
	IF(error/='') GOTO 9999

    CALL ExtractSolution(ndof_el,member,coord,memb_info,x,dof_con,sol_pt(i,:,:),sol_mb(i,:,:))


!--------------------------------------------------------
!    CALCULATE DERIVATIVES FOR EACH STEP.
!    2/dt A(t+dt)+dot{A(t+dt)}= 2/dt A(t+dt)+2/dt*[A(t+dt)-A(t)]-\dot{A(t)}=4/dt*A(t+dt)-[2/dt*A(t)+\dot{A(t)}]
!------------------------------------------------
    IF(analysis_flag==2) THEN
		init_cond(:,1:6)=2.d0*two_divide_dt*sol_mb(i,:,4:9)-init_cond(:,1:6) 
		init_cond(:,7:12)=2.d0*two_divide_dt*CTCabPH(niter, member,memb_info ,sol_mb(i,:,4:))-init_cond(:,7:12)
    ENDIF
    
    IF(analysis_flag==3) THEN

	    IF(niter==1) x=0.0d0  ! eigenvalue analysis based on linear theory
	    CALL EigenSolve(ndof_el,memb_info,v_root_a,omega_a,member,pt_condition,niter,error,&
					& ncond_mb,mb_condition,distr_fun,(dof_con),x,&
					& nev,eigen_val,eigen_vec)
		IF(error/='') GOTO 9999
       DO j=1,nev
			CALL ExtractSolution(ndof_el,member,coord,memb_info,eigen_vec(:,j),dof_con,eigen_vec_pt(j,:,:),eigen_vec_mb(j,:,:))
        ENDDO
	ENDIF
ENDDO

DEALLOCATE(cond_all,follower_all,dof_all,x,index_kp,index_mb)

DO i=1,SIZE(member,1)
	DEALLOCATE(memb_info(i)%coordinate,memb_info(i)%triad,memb_info(i)%Le,memb_info(i)%mate)
ENDDO
IF(analysis_flag==3) DEALLOCATE(eigen_vec)

9999 IF(error/='') THEN
		IF(ALLOCATED(eigen_vec))DEALLOCATE(eigen_vec)
		IF(ALLOCATED(cond_all)) DEALLOCATE (cond_all)
 		IF(ALLOCATED(follower_all)) DEALLOCATE (follower_all)
 		IF(ALLOCATED(dof_all)) DEALLOCATE (dof_all)
 		IF(ALLOCATED(x)) DEALLOCATE (x)
	 ENDIF


END SUBROUTINE Analysis
!************************************************************
