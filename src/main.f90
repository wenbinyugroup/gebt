!***********************************************************
!*                          GEBT                           *
!*                                                         *
!*          --Geometrically Exact Beam Theory ---          *
!*                                                         *
!*=========================================================*
!*                                                         *
!*               Copyright (c) Wenbin Yu                   *
!*              All rights reserved                        *
!*          Mechanical and Aerospace Engineering           *
!*                                                         *
!*                Utah State University                    *
!*=========================================================*
!*                                                         *
!*         Created on      : June 26, 2008                 *
!*         Modified on     : July 15, 2009                 *
!*         Modified on     : July 6,  2010                 *
!*=========================================================*
!*                                                         *
!*                OBJECTIVE OF GEBT                        *
!*            =============================                *
!*                                                         * 
!* This code implements the mixed-formulation of intrinsic *
!* formulation of the geometrically exact beam theory      *
!* developed by Prof. Hodges. It is fully compatible with  *
!* VABS and it is designed to work together with VABS      *
!*                                                         * 
!*                      Alpha Version                      *
!*            =============================                *
!*                                                         * 
!* This version has the capability for linear static       *
!* analysis, which was finished during a summer faculty    *
!* fellowship at WPAFB                                     *
!*                                                         *
!*                      Beta Version                       *
!*            =============================                *
!*                                                         * 
!* This version uses the DLL concept. Seeks a steady state *  
!* solution of the linear theory.                          *
!*                                                         *
!*                                                         *
!*                      1.0                                *
!*            =============================                *
!*                                                         * 
!* This version enables the nonlinearity capability of GEBT*  
!*                                                         *
!*                      1.5                                *
!*            =============================                *
!*                                                         * 
!* This version enables the follower force capability      *  
!*                                                         *
!*                                                         *
!*                      2.0                                *
!*            =============================                *
!*                                                         * 
!* This version enables AD capability using DNAD           *  
!*                                                         *
!*                                                         *
!*                      3.0                                *
!*            =============================                *
!*                                                         * 
!* This version enables dynamic response                   *  
!*                                                         *
!*                                                         *
!*                      4.0                                *
!*            =============================                *
!*                                                         * 
!* New assembly without storing the NxN matrix             *
!*                                                         *
!*                                                         *
!*                                                         *
!*                      4.2                                *
!*            =============================                *
!*                                                         * 
!* Use MA28 for linea solver and Arpack for eigen solver   *
!* Use long names of input file, including spaces in the   *
!* path and file names  (7/18/2011)                        *
!* input stiffness matrix (08/07/2011)                     *
!* fix a bug related with input stiffness matrix (10/25/2011)*
!*                                                         *
!*                                                         *
!*                                                         *
!***********************************************************

PROGRAM GEBT

!-----------------------------------------------------------------
! This serves as a general-purpose driving program
! Goto is used here for the purpose to terminate the execution
! due to critical errors.
!------------------------------------------------------------------
USE CPUTime
USE IO
 
IMPLICIT NONE ! Cancelling the naming convention of Fortran

INTERFACE

SUBROUTINE Analysis(nkp,nelem,ndof_el,nmemb,ncond_pt,nmate, nframe,ndistrfun,ncurv,coord, &
                 &  member,pt_condition,material,niter,nstep,sol_pt,sol_mb,error,ncond_mb, &
                 &  ntimefun,frame,mb_condition,distr_fun,curvature,omega_a0,omega_a_tf, &
                 &  v_root_a0,v_root_a_tf,simu_time,time_function,analysis_flag,init_cond, &
                 &  nev,eigen_val,eigen_vec_pt,eigen_vec_mb)

   	
!DEC$ ATTRIBUTES DLLIMPORT ::  Analysis
USE	GlobalDataFun
USE PrescribedCondition
USE TimeFunctionModule
			
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

END SUBROUTINE Analysis

END INTERFACE



CALL TIC ! start the clock


CALL Input ! Read inputs for the beam analysis: in IO
IF(error/='')	 GOTO 9999
WRITE(*,*) 
WRITE(*,*) 'Finished reading inputs for the beam analysis.'


CALL Analysis(nkp,nelem,ndof_el,nmemb,ncond_pt,nmate, nframe,ndistrfun,ncurv,coord, &
          &  member,pt_condition,material,niter,nstep,sol_pt,sol_mb,error,ncond_mb, &
          & ntimefun,frame,mb_condition,distr_fun,curvature,omega_a0,omega_a_tf, &
          & v_root_a0,v_root_a_tf,simu_time,time_function,analysis_flag,init_cond, &
          &  nev,eigen_val,eigen_vec_pt,eigen_vec_mb)
   IF(error/='')	 GOTO 9999
WRITE(*,*) 'Finished the beam analysis.'

CALL Output ! Report results: in IO
IF(error/='')	 GOTO 9999
WRITE(*,*)
WRITE(*,*) 'Finished outputing results for the beam analysis.'


WRITE(*,*)
WRITE(*,*)'GEBT finished successfully'

! Determine the time spend by the system
!===========================================
WRITE(*,*)
WRITE(*,*)'Program runs  ', TOC(), 'seconds'

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

9999 IF(error/='') WRITE(*,*) error 
     CALL WriteError(EIN,error) ! Write error to the echo file: in GEBTIO
     CLOSE(EIN)

END PROGRAM GEBT
!**************************************************************


