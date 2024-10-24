C *******************************************************************
C COPYRIGHT (c) 1977 Hyprotech UK.
C All rights reserved.
C
C None of the comments in this Copyright notice between the lines
C of asterisks shall be removed or altered in any way.
C
C This Package is intended for compilation without modification,
C so most of the embedded comments have been removed.
C
C ALL USE IS SUBJECT TO LICENCE. For full details of an HSL ARCHIVE
C Licence, see http://hsl.rl.ac.uk/archive/cou.html
C
C Please note that for an HSL ARCHIVE Licence:
C
C 1. The Package must not be copied for use by any other person.
C    Supply of any part of the library by the Licensee to a third party
C    shall be subject to prior written agreement between AEA
C    Hyprotech UK Limited and the Licensee on suitable terms and
C    conditions, which will include financial conditions.
C 2. All information on the Package is provided to the Licensee on the
C    understanding that the details thereof are confidential.
C 3. All publications issued by the Licensee that include results obtained
C    with the help of one or more of the Packages shall acknowledge the
C    use of the Packages. The Licensee will notify the Numerical Analysis
C    Group at Rutherford Appleton Laboratory of any such publication.
C 4. The Packages may be modified by or on behalf of the Licensee
C    for such use in research applications but at no time shall such
C    Packages or modifications thereof become the property of the
C    Licensee. The Licensee shall make available free of charge to the
C    copyright holder for any purpose all information relating to
C    any modification.
C 5. Neither CCLRC nor Hyprotech UK Limited shall be liable for any
C    direct or consequential loss or damage whatsoever arising out of
C    the use of Packages by the Licensee.
C *******************************************************************
C
C######DATE   09 MAR 1989
      SUBROUTINE MC19AD(N,NA,A,IRN,ICN,R,C,W)
      INTEGER   N,NA,IRN(*),ICN(*)
      DOUBLE PRECISION A(*)
C      REAL R(N),C(N),W(N,5)
      DOUBLE PRECISION R(N),C(N),W(N,5)
      INTEGER LP,IFAIL
      COMMON/MC19BD/LP,IFAIL
      INTEGER I,I1,I2,ITER,J,K,L,MAXIT
C      REAL E,E1,EM,Q,Q1,QM,S,S1,SM,SMIN,U,V
	DOUBLE PRECISION E,E1,EM,Q,Q1,QM,S,S1,SM,SMIN,U,V
      EXTERNAL MC19CD
      INTRINSIC ALOG,DABS,FLOAT
      DATA MAXIT/100/,SMIN/0.1/
      IFAIL=1
      IF(N.LT.1)GO TO 230
      IFAIL=2
      IFAIL=0
      DO 5 I=1,N
      C(I)=0.
      R(I)=0.
5     CONTINUE
      DO 10 L=1,4
      DO 10 I=1,N
      W(I,L)=0.
10    CONTINUE
      IF(NA.LE.0)GO TO 250
      DO 30 K=1,NA
      U=DABS(A(K))
      IF(U.EQ.0.)GO TO 30
C      U=ALOG(U)
      U=LOG(U)
      I1=IRN(K)
      I2=ICN(K)
      IF(I1.GE.1 .AND. I1.LE.N .AND. I2.GE.1 .AND. I2.LE.N)GO TO 20
      IF(LP.GT.0)WRITE(LP,15)K,I1,I2
15    FORMAT(20H MC19 ERROR. ELEMENT,I5,10H IS IN ROW,I5,
     1 8H AND COL,I5)
      IFAIL=3
      GO TO 30
20    W(I1,1)=W(I1,1)+1.
      W(I2,2)=W(I2,2)+1.
      R(I1)=R(I1)+U
      W(I2,3)=W(I2,3)+U
   30 CONTINUE
      IF(IFAIL.EQ.3)GO TO 230
      DO 70 I=1,N
      IF(W(I,1).EQ.0.)W(I,1)=1.
      R(I)=R(I)/W(I,1)
      W(I,5)=R(I)
      IF(W(I,2).EQ.0.)W(I,2)=1.
      W(I,3)=W(I,3)/W(I,2)
70    CONTINUE
      SM=SMIN*FLOAT(NA)
      DO 80 K=1,NA
       IF(A(K).EQ.0.0 )GO TO 80
      I=IRN(K)
      J=ICN(K)
      R(I)=R(I)-W(J,3)/W(I,1)
   80 CONTINUE
      E=0.
      Q=1.
      S=0.
      DO 100 I=1,N
      S=S+W(I,1)*R(I)**2
100   CONTINUE
      IF(S.LE.SM)GO TO 186
      DO 185 ITER=1,MAXIT
      DO 130 K=1,NA
      IF(A(K).EQ.0.)GO TO 130
      I=ICN(K)
      J=IRN(K)
      C(I)=C(I)+R(J)
  130 CONTINUE
      S1=S
      S=0.
      DO 140 I=1,N
      V=-C(I)/Q
      C(I)=V/W(I,2)
      S=S+V*C(I)
140   CONTINUE
      E1=E
      E=Q*S/S1
      Q=1.-E
      IF(S.LE.SM)E=0.
      DO 150 I=1,N
      R(I)=R(I)*E*W(I,1)
150   CONTINUE
      IF(S.LE.SM)GO TO 190
      EM=E*E1
      DO 152 K=1,NA
      IF(A(K).EQ.0.0 ) GO TO 152
      I=IRN(K)
      J=ICN(K)
      R(I)=R(I)+C(J)
152   CONTINUE
      S1=S
      S=0.
      DO 155 I=1,N
      V=-R(I)/Q
      R(I)=V/W(I,1)
      S=S+V*R(I)
155   CONTINUE
      E1=E
      E=Q*S/S1
      Q1=Q
      Q=1.-E
      IF(S.LE.SM)Q=1.
      QM=Q*Q1
      DO 160 I=1,N
      W(I,4)=(EM*W(I,4)+C(I))/QM
      W(I,3)=W(I,3)+W(I,4)
160   CONTINUE
      IF(S.LE.SM)GO TO 186
      DO 180 I=1,N
      C(I)=C(I)*E*W(I,2)
180   CONTINUE
185   CONTINUE
186   DO 188 I=1,N
      R(I)=R(I)*W(I,1)
188   CONTINUE
190   DO 200 K=1,NA
      IF(A(K).EQ.0.0 )GO TO 200
      I=IRN(K)
      J=ICN(K)
      R(I)=R(I)+W(J,3)
  200 CONTINUE
      DO 220 I=1,N
      R(I)=R(I)/W(I,1)-W(I,5)
      C(I)=-W(I,3)
220   CONTINUE
      GO TO 250
230   IF(LP.GT.0)WRITE(LP,240)IFAIL
240   FORMAT(//13H ERROR RETURN,I2,10H FROM MC19)
250   RETURN
      END
      BLOCK DATA MC19CD
      INTEGER LP,IFAIL
      COMMON/MC19BD/LP,IFAIL
      DATA LP/6/
      END
