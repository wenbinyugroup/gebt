!***********************************************************
!*               Copyright (c) Wenbin Yu                   *
!*              All rights reserved                        *
!*          Mechanical and Aerospace Engineering           *
!*                                                         *
!*                Utah State University                    *
!***********************************************************

!***************************************************************
!*                                                             *
!* This module contains information and calculation for an     *
!* element within a member 
!* Outputs: ElemEqn,ElemJacobian                               *
!***************************************************************

MODULE Element

USE InternalData
USE PrescribedCondition

IMPLICIT NONE

PRIVATE       ! So everything is private, except declared by PUBLIC

PUBLIC ElemEqn,ElemJacobian,ElemMass,ExtractElementProperties
PUBLIC load,exist_load,follower_load

! private data needed inside this module
!--------------------------------------------
REAL(DBL)::dL              ! length of the element
REAL(DBL)::Le              ! the ending arc length of the current element
REAL(DBL)::eCab(NDIM,NDIM) ! direction cosine matrix of the undeformed element
REAL(DBL)::eFlex(NSTRN,NSTRN) !flexibility matrix of the elment
REAL(DBL)::eMass(NSTRN,NSTRN) ! inverse of the mass matrix of the element
TYPE(DistriLoad)::load        ! distributed load
LOGICAL::exist_load,follower_load ! flags to indicate whether distributed load exist and whether they are follower forces

REAL(DBL)::Ui(NDIM),theta(NDIM),Fi(NDIM),Mi(NDIM),e1gamma(NDIM),kappa(NDIM),ePi(NDIM),Hi(NDIM),Vi(NDIM),OMEGAi(NDIM) ! e1gamma=e1+gamma
REAL(DBL)::ev_i(NDIM)    ! initial velocity of the mid point of the element
REAL(DBL)::eOmega_a(NDIM) ! initial angular velocity of the element

REAL(DBL)::eCT(NDIM,NDIM)         ! the transpose of the direction cosine matrix corresponding to elastic rotation
REAL(DBL)::eCTCAB(NDIM,NDIM)      ! eCT.Cab
REAL(DBL)::eCabhalfL(NDIM,NDIM)   ! eCab*dL/2
REAL(DBL)::eCTCabhalfL(NDIM,NDIM) ! eCTCab*dL/2

REAL(DBL)::UiDot(NDIM),ThetaDot(NDIM),CTCabPdot(NDIM),CTCabHdot(NDIM)

!=============================================

CONTAINS


!************************************************************
!*                                                          *
!* Evaluate nonlinear equations for each element, including *
!* the equation value due to x_0 and the term due to        *
!* distributed load. Note, it is put in the right hand side *
!*															*
!************************************************************ 
FUNCTION ElemEqn(ndof_el)

IMPLICIT NONE

INTEGER,INTENT(IN)  :: ndof_el

REAL(DBL)::ElemEqn(ndof_el+NDOF_ND) ! the functional value for each element
REAL(DBL):: tmpR(NDIM)

ElemEqn(1:3)=MATMUL(eCTCab,Fi)
ElemEqn(ndof_el+1:ndof_el+3)=-ElemEqn(1:3)

ElemEqn(4:6)= MATMUL(eCTCabhalfL,CrossProduct(e1gamma,Fi))
tmpR=MATMUL(eCTCab,Mi)
ElemEqn(ndof_el+4:ndof_el+6)=-tmpR+ElemEqn(4:6)
ElemEqn(4:6)=tmpR+ElemEqn(4:6)

ElemEqn(7:9)=MATMUL(eCTCabhalfL,e1gamma)-eCabhalfL(:,1)
ElemEqn(ndof_el+7:ndof_el+9)=Ui+ElemEqn(7:9)
ElemEqn(7:9)=-Ui+ElemEqn(7:9)

tmpR=MATMUL(eCabhalfL,kappa)
!ElemEqn(10:12)=tmpR+CrossProduct(theta*0.5d0,tmpR)+MATMUL(OuterProduct(theta,theta),tmpR*0.25d0)
!-----------------------------------------------------------------------------
! WM: Calculation of ElemEqn(10:12) has been changed for Wiener-Milenkovic parameters
ElemEqn(10:12)=tmpR-MATMUL(DOT_PRODUCT(theta,theta)*I3/16.0D0,tmpR) + &
               & CrossProduct(theta*0.5D0,tmpR) + &
               & MATMUL(OuterProduct(theta,theta)/8.0D0,tmpR)
! WM: End calculation of ElemEqn(10:12)
!-----------------------------------------------------------------------------
ElemEqn(ndof_el+10:ndof_el+12)=theta+ElemEqn(10:12)
ElemEqn(10:12)=-theta+ElemEqn(10:12)

IF(exist_load) THEN
	ElemEqn(1:6)=ElemEqn(1:6)+GetLoad(-1,dL,Le,eCT,load,follower_load)  ! add the contribution due to distributed load
	ElemEqn(ndof_el+1:ndof_el+6)=ElemEqn(ndof_el+1:ndof_el+6)+GetLoad(1,dL,Le,eCT,load,follower_load)  ! add the contribution due to distributed load
ENDIF

! for dynamic analysis
!------------------------------------------
IF(ndof_el==18) THEN

	! modify the right hand side for dynamic analysis including 
	! both steady state and initial conditions and time marching
	!--------------------------------------------------
	
	tmpR=CrossProduct(eOmega_a,MATMUL(eCTCabhalfL,ePi))
	IF(init_flag==2) tmpR=tmpR+two_divide_dt*MATMUL(eCTCabhalfL,ePi)	
	ElemEqn(1:3)=ElemEqn(1:3)-tmpR
    ElemEqn(ndof_el+1:ndof_el+3)=ElemEqn(ndof_el+1:ndof_el+3)-tmpR
    
	tmpR=CrossProduct(eOmega_a,MATMUL(eCTCabhalfL,Hi))+MATMUL(eCTCabhalfL,CrossProduct(Vi,ePi))
    IF(init_flag==2) tmpR=tmpR+two_divide_dt*MATMUL(eCTCabhalfL,Hi)
	ElemEqn(4:6)=ElemEqn(4:6)-tmpR
    ElemEqn(ndof_el+4:ndof_el+6)=ElemEqn(ndof_el+4:ndof_el+6)-tmpR

	! the functions of fpi and fhi
	!------------------------------------------
	ElemEqn(NDOF_ND+1:NDOF_ND+3)=ev_i+CrossProduct(eOmega_a,Ui)-MATMUL(eCTCab,Vi)
    IF(init_flag==2) ElemEqn(NDOF_ND+1:NDOF_ND+3)=ElemEqn(NDOF_ND+1:NDOF_ND+3)+two_divide_dt*Ui
	ElemEqn(NDOF_ND+4:NDOF_ND+6)=MATMUL(TRANSPOSE(eCTCab),eOmega_a)-OMEGAi
ENDIF

IF(init_flag/=0) THEN !needed for the initial step and time marching
	ElemEqn(1:3)=ElemEqn(1:3)-dL*0.5D0*CTCabPdot
	ElemEqn(ndof_el+1:ndof_el+3)=ElemEqn(ndof_el+1:ndof_el+3)-dL*0.5D0*CTCabPdot
    
	ElemEqn(4:6)=ElemEqn(4:6)-dL*0.5D0*CTCabHdot
	ElemEqn(ndof_el+4:ndof_el+6)=ElemEqn(ndof_el+4:ndof_el+6)-dL*0.5D0*CTCabHdot

    ElemEqn(NDOF_ND+1:NDOF_ND+3)=ElemEqn(NDOF_ND+1:NDOF_ND+3)+UiDot
!	ElemEqn(NDOF_ND+4:NDOF_ND+6)=ElemEqn(NDOF_ND+4:NDOF_ND+6)+ &
!	&   MATMUL(TRANSPOSE(eCab),(ThetaDot-CrossProduct(theta*.5D0,ThetaDot))/(1+DOT_PRODUCT(theta,theta)*0.25D0) )
!-----------------------------------------------------------------------------
! WM: Calculation of ElemEqn(NDOF_ND+4:NDOF_ND+6) has been changed for Wiener-Milenkovic parameters
	ElemEqn(NDOF_ND+4:NDOF_ND+6)=ElemEqn(NDOF_ND+4:NDOF_ND+6)+ &
	&	MATMUL(TRANSPOSE(eCab),(4.0D0*ThetaDot-0.25D0*DOT_PRODUCT(theta,theta)*ThetaDot- &
	&	2.0D0*MATMUL(Tilde(theta),ThetaDot) + 0.5D0*MATMUL(OuterProduct(theta,theta),ThetaDot)) &
	&   / (4-(2-DOT_PRODUCT(theta,theta)))**2)
! WM: End calculation of ElemEqn(NDOF_ND+4:NDOF_ND+6)
!-----------------------------------------------------------------------------
ENDIF
END FUNCTION ElemEqn
!***************************************************


!************************************************************
!*                                                          *
!*  Caculate the Jacobian matrix for each element           *
!*                                                          *
!*															*
!************************************************************ 
SUBROUTINE ElemJacobian(ndof_el,niter,elemJac)

IMPLICIT NONE

INTEGER,INTENT(IN)  :: ndof_el,niter

REAL(DBL),INTENT(OUT)::ElemJac(:,:) ! the coefficient matrix for each element

REAL(DBL)::ekttek(NDIM,NDIM,NDIM)   ! e_k.\theta^T+ \theta.e_k^T
REAL(DBL)::eCTtheta(NDIM,NDIM,NDIM)  !derivatives of C^T w.r.t theta, eCTtheta(k,:,:)= d_C^T/d_thetak
REAL(DBL)::temp31(NDIM),temp32(NDIM),temp33(NDIM,NDIM),temp333(NDIM,NDIM,NDIM) ! temporary arrays 
REAL(DBL)::Qpart1(NDIM,NDIM), Qpart2(NDIM,NDIM), Qpart3(NDIM,NDIM),Qpart4(NDIM,NDIM)

REAL(DBL)::tilde_omega(NDIM,NDIM),tt4,tt2

INTEGER:: i

ElemJac=0.0D0

CALL CT_THETA(theta,eCT,ekttek,eCTtheta)  ! obtain derivatives of C^T w.r.t theta

! d_fu/d_theta
!--------------------
temp31=MATMUL(eCab,Fi)
ElemJac(1:3,4:6)=-MATMUL3(eCTtheta,temp31)   ! coefficient of theta of Fu 
ElemJac(ndof_el+1:ndof_el+3,4:6)=-ElemJac(1:3,4:6)

! d_fu/d_F
!--------------------		
ElemJac(1:3,7:9)=-eCTCab  
ElemJac(ndof_el+1:ndof_el+3,7:9)=eCTCab

! d_fpsi/d_theta
!--------------------
temp31=MATMUL(eCab,Mi)
temp32=-MATMUL(eCab,CrossProduct(e1gamma,dL*0.5d0*Fi))
ElemJac(4:6,4:6)=MATMUL3(eCTtheta,temp32-temp31)  ! coefficient of theta of Fpsi1 
ElemJac(ndof_el+4:ndof_el+6,4:6)=MATMUL3(eCTtheta,temp32+temp31) 

! d_fpsi/d_F
!--------------------
ElemJac(4:6,7:9)=-MATMUL(eCTCabhalfL,Tilde(e1gamma)-MATMUL(Tilde(Fi),eFlex(1:3,1:3)))  
ElemJac(ndof_el+4:ndof_el+6,7:9)=ElemJac(4:6,7:9)

! d_fpsi/d_M
!--------------------
temp33=MATMUL(eCTCabhalfL,MATMUL(Tilde(Fi),eFlex(1:3,4:6)))	
ElemJac(4:6,10:12)=temp33-eCTCab
ElemJac(ndof_el+4:ndof_el+6,10:12)=temp33+eCTCab

!write(30,'(12e17.6)')(elemJacobian(i,:),i=4,6)
!write(30,'(12e17.6)')(elemJacobian(ndof_el+i,:),i=4,6)

! d_fF/d_u
!--------------------
ElemJac(7:9,1:3)=I3 
ElemJac(ndof_el+7:ndof_el+9,1:3)=-I3 

! d_fF/d_theta
!--------------------
temp31=-MATMUL(eCabhalfL,e1gamma)
ElemJac(7:9,4:6)=MATMUL3(eCTtheta,temp31)   
ElemJac(ndof_el+7:ndof_el+9,4:6)=ElemJac(7:9,4:6)


! d_fF/d_F,d_fF/d_M 
!--------------------
ElemJac(7:9, 7:12)=-MATMUL(eCTCabhalfL,eFlex(1:3,:)) 
ElemJac(ndof_el+7:ndof_el+9,7:12)=ElemJac(7:9, 7:12)

! d_fM/d_theta
!--------------------
DO i=1,NDIM
!	temp333(i,:,:)=Tilde(I3(i,:))*0.5d0+ekttek(i,:,:)*0.25d0  ! \tilde{ek}/2+ekttek/4
!-----------------------------------------------------------------------------
! WM: Calculation of temp333(:,:,:) has been changed for Wiener-Milenkovic parameters
	temp333(i,:,:)=-0.125D0*theta(i)*I3+Tilde(I3(i,:))*0.5d0+ekttek(i,:,:)*0.125d0
! WM: End calculation of temp333(:,:,:)
!-----------------------------------------------------------------------------
ENDDO

temp33=MATMUL3(temp333,MATMUL(eCabhalfL,kappa))
ElemJac(10:12,4:6)=I3-temp33
ElemJac(ndof_el+10:ndof_el+12,4:6)=-I3-temp33


! d_fM/d_F,d_fM/d_M 
!--------------------
!ElemJac(10:12,7:12)=-MATMUL(MATMUL(I3+Tilde(theta)*.5d0+OuterProduct(theta,theta)*.25d0,eCabhalfL),& 
!                        & eFlex(4:6,:)) 
!-----------------------------------------------------------------------------
! WM: Calculation of ElemJac(10:12,7:12) has been changed for Wiener-Milenkovic parameters
ElemJac(10:12,7:12)=-MATMUL(MATMUL((1-DOT_PRODUCT(theta,theta)/16.0D0)*I3+Tilde(theta)*.5d0+ &
				    & OuterProduct(theta,theta)*.125d0,eCabhalfL),& 
                        & eFlex(4:6,:)) 
! WM: End calculation of ElemJac(10:12,7:12)
!-----------------------------------------------------------------------------

ElemJac(ndof_el+10:ndof_el+12,7:12)=ElemJac(10:12,7:12)

! add the contribution due to follower distributed load:\diff C^T/\theta.f
!----------------------------------------------------------
IF(follower_load)THEN
	ElemJac(1:6,4:6)=ElemJac(1:6,4:6)-GetLoadJ(-1,dL,Le,eCTtheta,load)
	ElemJac(ndof_el+1:ndof_el+6,4:6)=ElemJac(ndof_el+1:ndof_el+6,4:6)-GetLoadJ(1,dL,Le,eCTtheta,load)  
ENDIF

! for dynamic analysis
!------------------------------------------
IF(ndof_el==18) THEN

	tilde_omega=Tilde(eOmega_a)
	IF(init_flag==2) tilde_omega=tilde_omega+two_divide_dt*I3

    ! additional terms for fui due to P
	!--------------------------------------------
	! coefficient of theta of Fu due to ePi
    temp31=MATMUL(eCabhalfL,ePi)
	temp33=MATMUL(tilde_omega,MATMUL3(eCTtheta,temp31))
	ElemJac(1:3,4:6)=ElemJac(1:3,4:6)+temp33   
	ElemJac(ndof_el+1:ndof_el+3,4:6)=ElemJac(ndof_el+1:ndof_el+3,4:6)+temp33   

	! coefficient of p of Fu
	ElemJac(1:3,13:15)=MATMUL(tilde_omega,eCTCabhalfL)    
	ElemJac(ndof_el+1:ndof_el+3,13:15)=ElemJac(1:3,13:15)

    ! additional terms for fpsi due to P,H
	!--------------------------------------------
	! coefficient of theta of fpsi due to P,H
    temp33=MATMUL3(eCTtheta,MATMUL(eCabHalfL,CrossProduct(Vi,ePi)))+ MATMUL(tilde_omega,MATMUL3(eCTtheta,MATMUL(eCabhalfL,Hi)))
	ElemJac(4:6,4:6)=ElemJac(4:6,4:6)+ temp33  
    ElemJac(ndof_el+4:ndof_el+6,4:6)=ElemJac(ndof_el+4:ndof_el+6,4:6)+ temp33  

	! coefficient of P,H of Fpsi
	ElemJac(4:6,13:15)=MATMUL(eCTCabhalfL,Tilde(Vi)-MATMUL(Tilde(ePi),eMass(1:3,1:3)))  
	ElemJac(4:6,16:18)=MATMUL(tilde_omega,eCTCabhalfL)-MATMUL(eCTCabhalfL,MATMUL(Tilde(ePi),eMass(1:3,4:6)))  
    ElemJac(ndof_el+4:ndof_el+6,13:18)=ElemJac(4:6,13:18)
	
	! the Jacobian of fPi
	!------------------------------------------
	ElemJac(NDOF_ND+1:NDOF_ND+3,1:3)=-tilde_omega
    ElemJac(NDOF_ND+1:NDOF_ND+3,4:6)=MATMUL3(eCTtheta,MATMUL(eCab,Vi))
    ElemJac(NDOF_ND+1:NDOF_ND+3,13:18)=MATMUL(eCTCab,eMass(1:3,:))
    
	! the Jacobian of fHi
	!------------------------------------------
	DO i=1,NDIM
		eCTtheta(i,:,:)=TRANSPOSE(eCTtheta(i,:,:))
	ENDDO
	ElemJac(NDOF_ND+4:NDOF_ND+6,4:6)=-MATMUL(TRANSPOSE(eCab),MATMUL3(eCTtheta,eOmega_a))
    ElemJac(NDOF_ND+4:NDOF_ND+6,13:18)=eMass(4:6,:)    

	IF(init_flag==2) THEN
	   IF(niter/=1)THEN  ! nonlinear solution
!			tt4=DOT_PRODUCT(theta,theta)*0.25D0
!-----------------------------------------------------------------------------
! WM: Calculation of temp33(:,:) has been changed for Wiener-Milenkovic parameters
    	    tt4=2.0D0+DOT_PRODUCT(theta,theta)*0.125D0
    	    tt2=DOT_PRODUCT(theta,theta)
			DO i=1,NDIM
!				temp333(i,:,:)=(Tilde(I3(i,:)*0.5d0*(1+tt4))+MATMUL(I3-Tilde(theta*0.5d0),ekttek(i,:,:)*0.25d0))/(1+tt4)**2
                Qpart1=0.0D0
                Qpart2=0.0D0
                Qpart3=0.0D0
                Qpart4=0.0D0
                
                Qpart1(1,1)=-2.0D0*theta(i)/tt4**3
                Qpart1(2,2)=Qpart1(1,1)
                Qpart1(3,3)=Qpart1(1,1)
                
                Qpart2(1,1)=theta(i)*tt2/(8*tt4**3) - theta(i)/(2*tt4**2)
                Qpart2(2,2)=Qpart2(1,1)
                Qpart2(3,3)=Qpart2(1,1)
                
                Qpart3(1,2)=-theta(i)*theta(3)/tt4**3
                Qpart3(1,3)=theta(i)*theta(2)/tt4**3
                Qpart3(2,1)=-Qpart3(1,2)
                Qpart3(2,3)=-theta(i)*theta(1)/tt4**3
                Qpart3(3,1)=-Qpart3(1,3)
                Qpart3(3,2)=-Qpart3(2,3)
                
                Qpart3 = Qpart3-Tilde(I3(i,:))*(2.0D0/tt4**2)
                
                Qpart4(1,1)=theta(1)**2
				Qpart4(1,2)=theta(1)*theta(2)
				Qpart4(1,3)=theta(1)*theta(3)
	 			Qpart4(2,1)=Qpart4(1,2)
				Qpart4(2,2)=theta(2)**2
				Qpart4(2,3)=theta(2)*theta(3)
				Qpart4(3,1)=Qpart4(1,3)
				Qpart4(3,2)=Qpart4(2,3)
				Qpart4(3,3)=theta(3)**2
				
				Qpart4=Qpart4*(-theta(i)/tt4**3)+2.0D0*ekttek(i,:,:)/tt4**2
				Qpart4=0.25D0*Qpart4
				
				temp333(i,:,:)=Qpart1+Qpart2+Qpart3+Qpart4
			ENDDO
			temp33=MATMUL3(temp333,ThetaDot)
!			temp33=temp33-(I3-Tilde(theta*.5D0))/(1+tt4)*two_divide_dt
            temp33=temp33-((4.0D0-0.25D0*tt2)*I3-2.0D0*Tilde(theta)+ &
                   & 0.5D0*OuterProduct(theta,theta))/tt4**2
! WM: End calculation of temp33(:,:)
!-----------------------------------------------------------------------------
        ELSE  ! linear solution
		    temp33=-I3*two_divide_dt
		ENDIF

		ElemJac(NDOF_ND+4:NDOF_ND+6,4:6)=ElemJac(NDOF_ND+4:NDOF_ND+6,4:6)+MATMUL(TRANSPOSE(eCab),temp33)
	ENDIF

ENDIF

IF(init_flag==1) THEN !needed for the initial step calculating the initial conditions
	ElemJac(:,1:6)=0.D0
	ElemJac(1:3,1:3)=dL*.5d0*I3
    ElemJac(4:6,4:6)=ElemJac(1:3,1:3)
    ElemJac(ndof_el+1:ndof_el+3,1:3)=ElemJac(1:3,1:3)
    ElemJac(ndof_el+4:ndof_el+6,4:6)=ElemJac(1:3,1:3)
ENDIF

!IF(DEBUG)THEN
!	WRITE(IOUT,*) 'ELEMENT COEFFICIENT MATRIX'
!	DO i=1, ndof_nd+ndof_el
!		WRITE(IOUT,'(1X,12ES20.12)')elemJac(i,:)
!	ENDDO
!ENDIF
END SUBROUTINE ElemJacobian
!***************************************************



!************************************************************
!*                                                          *
!*  Caculate the mass matrix for each element               *
!*                                                          *
!*															*
!************************************************************ 
SUBROUTINE ElemMass(elemM)

IMPLICIT NONE

REAL(DBL),INTENT(OUT)::ElemM(:,:) ! the mass matrix for each element

ElemM=0.0D0

! d_fu/d_theta_t
!--------------------
ElemM(1:3,4:6)=CT_THETA_T(theta,MATMUL(eCabhalfL,ePi))
ElemM(19:21,4:6)=ElemM(1:3,4:6)

! d_fu/d_P_t
!--------------------		
ElemM(1:3,13:15)=eCTCabhalfL  
ElemM(19:21,13:15)=eCTCabhalfL

! d_fpsi/d_theta_t
!--------------------
ElemM(4:6,4:6)=CT_THETA_T(theta,MATMUL(eCabhalfL,Hi)) 
ElemM(22:24,4:6)=ElemM(4:6,4:6)

! d_fpsi/d_H_t
!--------------------
ElemM(4:6,16:18)=eCTCabhalfL 
ElemM(22:24,16:18)=eCTCabhalfL 
	
! d_fPi/d_u_t
!------------------------------------------
ElemM(13:15,1:3)=-I3

    
! d_fHi/d_theta_t
!------------------------------------------
!ElemM(16:18,4:6)=MATMUL(TRANSPOSE(eCab),(Tilde(theta*.5D0)-I3)/(1+DOT_PRODUCT(theta,theta)*0.25D0))
!-----------------------------------------------------------------------------
! WM: Calculation of ElemM(16:18,4:6) has been changed for Wiener-Milenkovic parameters
ElemM(16:18,4:6)=-MATMUL(TRANSPOSE(eCab),((4.0D0-DOT_PRODUCT(theta,theta)*0.25D0*I3 - &
                 & 2.0D0*Tilde(theta)+0.5D0*OuterProduct(theta,theta))/(2.0D0+DOT_PRODUCT(theta,theta)*0.125D0)**2.0D0))
! WM: End calculation of ElemM(16:18,4:6)
!-----------------------------------------------------------------------------

END SUBROUTINE ElemMass
!***************************************************



!************************************************************
!*                                                          *                                      
!*  Extract element properties needed for element assembly  *
!************************************************************
SUBROUTINE  ExtractElementProperties(elem_no,memb_info_i,x_elem,v_root_a,omega_a,ndof_el,init_elem)  

INTEGER,INTENT(IN)::elem_no
TYPE (MemberInf),INTENT(IN)::memb_info_i

REAL(DBL),INTENT(IN)::x_elem(:) ! the 12 variables of the element u_i, \theta_i, F_i, M_i, for dynamic analysis, 6 more Pi, Hi
                                ! for initial step, x_elem contains, CTCabPdot, CTCabHdot, F_i, M_i, P_i, H_i
REAL(DBL),INTENT(IN)::v_root_a(:),omega_a(:)

REAL(DBL),INTENT(IN)::init_elem(:) ! initial step: the 12 initial values of the element u_i, \theta_i, \dot{u}_i, \dot{theta}_i 
                                   ! time marching: 2/dt ui+\dot{u}_i, 2/dt thetai+\dot{theta}_i, 2/dt CTCabP+\dot, 2/dt CTCabP+dot
                                   
!REAL(DBL)::t2
INTEGER,INTENT(IN)::ndof_el

!----------------------------------------------------
! Extract sectional properties for this element
!----------------------------------------------------
eFlex=memb_info_i%mate(elem_no,1:NSTRN,:) ! element flexibility matrix
IF(ndof_el==18) eMass=memb_info_i%mate(elem_no,NSTRN+1:,:)! element inverse mass matrix for all the analyses which are not static

!------------------------------------------------------------------

!----------------------------------------------------------------------------
! Extract b frame for this element
!----------------------------------------------------------------------------
eCab=memb_info_i%triad(elem_no,:,:)
dL=memb_info_i%dL
Le=memb_info_i%Le(elem_no) ! the ending x_1 of the element

!Angular velocity and linear velocity
!---------------------------------------
IF(ndof_el==18)	THEN
    eOmega_a=omega_a  ! angular velocity

	!linear velocity at the middle of the element
	!-----------------------------------------------
	ev_i=v_root_a+CrossProduct(eOmega_a,memb_info_i%coordinate(elem_no,:)-xyz_pt1)  
ENDIF

!----------------------------------------------------------------------------
! Extract internal variables for this element
!----------------------------------------------------------------------------
IF(init_flag==1) THEN  ! initial step
    Ui       = init_elem(1:3)          ! given initial displacements
	theta    = init_elem(4:6)          ! given initial rotations
	UiDot    = init_elem(7:9)          ! given initial linear velocities
	ThetaDot = init_elem(10:12)        ! given initial angular velocities
	CTCabPdot= x_elem(1:3)             ! calculated CTCabPdot from previous step
    CTCabHdot= x_elem(4:6)             ! calculated CTCabHdot from previous step
ELSE 
	Ui   =x_elem(1:3)
	theta=x_elem(4:6)
	IF(init_flag==2) THEN  ! time marching
		UiDot=-init_elem(1:3)
		ThetaDot=two_divide_dt*theta-init_elem(4:6)
		CTCabPdot=-init_elem(7:9)
		CTCabHdot=-init_elem(10:12)
	ENDIF
ENDIF

Fi   =x_elem(7:9)
Mi   =x_elem(10:12)

IF (ndof_el==18) THEN
	ePi   =x_elem(13:15)
	Hi    =x_elem(16:18)
ENDIF

!eFlex(4,4)=1/(1/eFlex(4,4)+(1/eFlex(5,5)+1/eFlex(6,6))*eFlex(1,1)*Fi(1))
!write(*,*)eflex(4,4)
e1gamma=MATMUL(eFlex(1:3,:),x_elem(7:12))+e1
kappa=MATMUL(eFlex(4:6,:),x_elem(7:12))

IF(ndof_el==18)THEN
	Vi=MATMUL(eMass(1:3,:),x_elem(13:18))
	OMEGAi=MATMUL(eMass(4:6,:),x_elem(13:18))
ENDIF

! Some numbers frequently needed for the element
!-----------------------------------------
eCabhalfL   =eCab*dL*0.5d0
eCT=DirCosineTRodrigues(theta)
eCTCab       =MATMUL(eCT,eCab)
eCTCabhalfL =eCTCab*dL*0.5d0


END SUBROUTINE ExtractElementProperties
!**********************************************************
   
END MODULE Element