* *******************************************************************
* COPYRIGHT (c) 1977 Hyprotech UK
* All rights reserved.
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of an HSL ARCHIVE
* Licence, see http://hsl.rl.ac.uk/archive/cou.html
*
* Please note that for an HSL ARCHIVE Licence:
*
* 1. The Package must not be copied for use by any other person.
*    Supply of any part of the library by the Licensee to a third party
*    shall be subject to prior written agreement between AEA
*    Hyprotech UK Limited and the Licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither CCLRC nor Hyprotech UK Limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*
*######DATE 22 Feb 1993
C       Toolpack tool decs employed.
C       SAVE statements added.
C       MA28JD reference removed.
C       ZERO made PARAMETER.
C
C  EAT 21/6/93 EXTERNAL statement put in for block data so will work on VAXs.
C
C
      SUBROUTINE MA28AD(N,NZ,A,LICN,IRN,LIRN,ICN,U,IKEEP,IW,W,IFLAG)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION U
      INTEGER IFLAG,LICN,LIRN,N,NZ
      DOUBLE PRECISION A(LICN),W(N)
      INTEGER ICN(LICN),IKEEP(N,5),IRN(LIRN),IW(N,8)
      DOUBLE PRECISION UPRIV
      INTEGER I,I1,IEND,II,J,J1,J2,JAY,JJ,KNUM,LENGTH,MOVE,NEWJ1,NEWPOS
      EXTERNAL MA30AD,MC20AD,MC22AD,MC23AD,MC24AD
      EXTERNAL MA28JD
      INTRINSIC DABS,DMAX1,MAX0
      COMMON /MA28ED/LP,MP,LBLOCK,GROW
      COMMON /MA28FD/EPS,RMIN,RESID,IRNCP,ICNCP,MINIRN,MINICN,IRANK,
     +       ABORT1,ABORT2
      COMMON /MA28GD/IDISP(2)
      COMMON /MA28HD/TOL,THEMAX,BIG,DXMAX,ERRMAX,DRES,CGCE,NDROP,MAXIT,
     +       NOITER,NSRCH,ISTART,LBIG
      COMMON /MA30ED/NLP,ABORTA,ABORTB,ABORT3
      COMMON /MA30FD/MIRNCP,MICNCP,MIRANK,MIRN,MICN
      COMMON /MA30ID/TOL1,BIG1,NDROP1,NSRCH1,LBIG1
      COMMON /MC23BD/MLP,NUMNZ,NUM,LARGE,ABORT
      DOUBLE PRECISION BIG,BIG1,CGCE,DRES,DXMAX,EPS,ERRMAX,RESID,RMIN,
     +                 THEMAX,TOL,TOL1
      INTEGER ICNCP,IRANK,IRNCP,ISTART,LARGE,LP,MAXIT,MICN,MICNCP,
     +        MINICN,MINIRN,MIRANK,MIRN,MIRNCP,MLP,MP,NDROP,NDROP1,NLP,
     +        NOITER,NSRCH,NSRCH1,NUM,NUMNZ
      LOGICAL ABORT,ABORT1,ABORT2,ABORT3,ABORTA,ABORTB,GROW,LBIG,LBIG1,
     +        LBLOCK
      INTEGER IDISP
      SAVE /MA28ED/,/MA28FD/,/MA28GD/,/MA28HD/,/MA30ED/,/MA30FD/,
     +     /MA30ID/,/MC23BD/
      IFLAG = 0
      ABORTA = ABORT1
      ABORTB = ABORT2
      ABORT = ABORT1
      MLP = LP
      NLP = LP
      TOL1 = TOL
      LBIG1 = LBIG
      NSRCH1 = NSRCH
      UPRIV = U
      IF (N.GT.0) GO TO 10
      IFLAG = -8
      IF (LP.NE.0) WRITE (LP,FMT=99999) N
      GO TO 210
   10 IF (NZ.GT.0) GO TO 20
      IFLAG = -9
      IF (LP.NE.0) WRITE (LP,FMT=99998) NZ
      GO TO 210
   20 IF (LICN.GE.NZ) GO TO 30
      IFLAG = -10
      IF (LP.NE.0) WRITE (LP,FMT=99997) LICN
      GO TO 210
   30 IF (LIRN.GE.NZ) GO TO 40
      IFLAG = -11
      IF (LP.NE.0) WRITE (LP,FMT=99996) LIRN
      GO TO 210
   40 DO 50 I = 1,NZ
        IF (IRN(I).GT.0 .AND. IRN(I).LE.N .AND. ICN(I).GT.0 .AND.
     +      ICN(I).LE.N) GO TO 50
        IF (IFLAG.EQ.0 .AND. LP.NE.0) WRITE (LP,FMT=99995)
        IFLAG = -12
        IF (LP.NE.0) WRITE (LP,FMT=99994) I,A(I),IRN(I),ICN(I)
   50 CONTINUE
      IF (IFLAG.LT.0) GO TO 220
      CALL MC20AD(N,NZ,A,ICN,IW,IRN,0)
      DO 60 I = 1,N
        IKEEP(I,2) = 0
        IKEEP(I,1) = 0
   60 CONTINUE
      MOVE = 0
      THEMAX = ZERO
      J1 = IW(1,1)
      DO 130 I = 1,N
        IEND = NZ + 1
        IF (I.NE.N) IEND = IW(I+1,1)
        LENGTH = IEND - J1
        IF (LENGTH.EQ.0) GO TO 130
        J2 = IEND - 1
        NEWJ1 = J1 - MOVE
        DO 120 JJ = J1,J2
          J = ICN(JJ)
          THEMAX = DMAX1(THEMAX,DABS(A(JJ)))
          IF (IKEEP(J,2).EQ.I) GO TO 110
          IKEEP(J,2) = I
          IKEEP(J,3) = JJ - MOVE - NEWJ1
          IF (MOVE.EQ.0) GO TO 120
          NEWPOS = JJ - MOVE
          A(NEWPOS) = A(JJ)
          ICN(NEWPOS) = ICN(JJ)
          GO TO 120
  110     MOVE = MOVE + 1
          LENGTH = LENGTH - 1
          JAY = IKEEP(J,3) + NEWJ1
          IF (MP.NE.0) WRITE (MP,FMT=99993) I,J,A(JJ)
          A(JAY) = A(JAY) + A(JJ)
          THEMAX = DMAX1(THEMAX,DABS(A(JAY)))
  120   CONTINUE
        IKEEP(I,1) = LENGTH
        J1 = IEND
  130 CONTINUE
      KNUM = NZ - MOVE
      IF (.NOT.LBLOCK) GO TO 140
      CALL MC23AD(N,ICN,A,LICN,IKEEP,IDISP,IKEEP(1,2),IKEEP(1,3),
     +            IKEEP(1,5),IW(1,3),IW)
      IF (IDISP(1).GT.0) GO TO 170
      IFLAG = -7
      IF (IDISP(1).EQ.-1) IFLAG = -1
      IF (LP.NE.0) WRITE (LP,FMT=99992)
      GO TO 210
  140 DO 150 I = 1,KNUM
        II = KNUM - I + 1
        NEWPOS = LICN - I + 1
        ICN(NEWPOS) = ICN(II)
        A(NEWPOS) = A(II)
  150 CONTINUE
      IDISP(1) = 1
      IDISP(2) = LICN - KNUM + 1
      DO 160 I = 1,N
        IKEEP(I,2) = I
        IKEEP(I,3) = I
  160 CONTINUE
      IKEEP(1,5) = -1
  170 IF (LBIG) BIG1 = THEMAX
      IF (NSRCH.LE.N) GO TO 180
      CALL MA30AD(N,ICN,A,LICN,IKEEP,IKEEP(1,4),IDISP,IKEEP(1,2),
     +            IKEEP(1,3),IRN,LIRN,IW(1,2),IW(1,3),IW(1,4),IW(1,5),
     +            IW(1,6),IW(1,7),IW(1,8),IW,UPRIV,IFLAG)
      GO TO 190
  180 CALL MA30AD(N,ICN,A,LICN,IKEEP,IKEEP(1,4),IDISP,IKEEP(1,2),
     +            IKEEP(1,3),IRN,LIRN,IW(1,2),IW(1,3),IW(1,4),IW(1,5),
     +            IW,IW,IW(1,6),IW,UPRIV,IFLAG)
  190 MINIRN = MAX0(MIRN,NZ)
      MINICN = MAX0(MICN,NZ)
      IRNCP = MIRNCP
      ICNCP = MICNCP
      IRANK = MIRANK
      NDROP = NDROP1
      IF (LBIG) BIG = BIG1
      IF (IFLAG.GE.0) GO TO 200
      IF (LP.NE.0) WRITE (LP,FMT=99991)
      GO TO 210
  200 I1 = IDISP(1) - 1
      IF (I1.NE.0) CALL MC22AD(N,ICN,A,I1,IKEEP(1,5),IKEEP(1,2),
     +                         IKEEP(1,3),IW,IRN)
      I1 = IDISP(1)
      IEND = LICN - I1 + 1
      IF (GROW) CALL MC24AD(N,ICN,A(I1),IEND,IKEEP,IKEEP(1,4),W)
      IF (GROW) W(1) = W(1) + THEMAX
      IF (GROW .AND. N.GT.1) W(2) = THEMAX
      IF (IFLAG.GE.0 .AND. MOVE.NE.0) IFLAG = -14
      GO TO 220
  210 CONTINUE
  220 RETURN
99999 FORMAT (' ERROR RETURN FROM MA28A/AD BECAUSE ',
     +       'N OUT OF RANGE = ',I10)
99998 FORMAT (' ERROR RETURN FROM MA28A/AD BECAUSE ',
     +       'NZ NON POSITIVE = ',I10)
99997 FORMAT (' ERROR RETURN FROM MA28A/AD BECAUSE ',
     +       'LICN TOO SMALL = ',I10)
99996 FORMAT (' ERROR RETURN FROM MA28A/AD BECAUSE ',
     +       'LIRN TOO SMALL = ',I10)
99995 FORMAT (' ERROR RETURN FROM MA28A/AD BECAUSE INDICES FOUND OUT ',
     +       'OF RANGE')
99994 FORMAT (1X,I6,'TH ELEMENT WITH VALUE ',1P,D12.4,/,20X,
     +       ' IS OUT OF RANGE WITH INDICES ',I8,' ,',I8)
99993 FORMAT (' DUPLICATE ELEMENT IN POSITION ',I8,',',I8,
     +       ' WITH VALUE ',1P,D12.4)
99992 FORMAT (' ERROR RETURN FROM MA28A/AD BECAUSE ',
     +       'ERROR RETURN FROM MC23A/AD')
99991 FORMAT (' ERROR RETURN FROM MA28A/AD BECAUSE ',
     +       'ERROR RETURN FROM MA30A/AD')
      END
      SUBROUTINE MA28BD(N,NZ,A,LICN,IVECT,JVECT,ICN,IKEEP,IW,W,IFLAG)
      INTEGER IFLAG,LICN,N,NZ
      DOUBLE PRECISION A(LICN),W(N)
      INTEGER ICN(LICN),IKEEP(N,5),IVECT(NZ),IW(N,5),JVECT(NZ)
      INTEGER I1,IDUP,IEND
      EXTERNAL MA28DD,MA30BD,MC24AD
      COMMON /MA28ED/MP,LP,LBLOCK,GROW
      COMMON /MA28FD/EPS,RMIN,RESID,IRNCP,ICNCP,MINIRN,MINICN,IRANK,
     +       ABORT1,ABORT2
      COMMON /MA28GD/IDISP(2)
      COMMON /MA28HD/TOL,THEMAX,BIG,DXMAX,ERRMAX,DRES,CGCE,NDROP,MAXIT,
     +       NOITER,NSRCH,ISTART,LBIG
      COMMON /MA30ED/NLP,ABORTA,ABORTB,ABORT3
      COMMON /MA30GD/MEPS,MRMIN
      COMMON /MA30ID/TOL1,BIG1,NDROP1,NSRCH1,LBIG1
      DOUBLE PRECISION BIG,BIG1,CGCE,DRES,DXMAX,EPS,ERRMAX,MEPS,MRMIN,
     +                 RESID,RMIN,THEMAX,TOL,TOL1
      INTEGER ICNCP,IRANK,IRNCP,ISTART,LP,MAXIT,MINICN,MINIRN,MP,NDROP,
     +        NDROP1,NLP,NOITER,NSRCH,NSRCH1
      LOGICAL ABORT1,ABORT2,ABORT3,ABORTA,ABORTB,GROW,LBIG,LBIG1,LBLOCK
      INTEGER IDISP
      SAVE /MA28ED/,/MA28FD/,/MA28GD/,/MA28HD/,/MA30ED/,/MA30GD/,
     +     /MA30ID/
      IF (NDROP.EQ.0) GO TO 10
      IFLAG = -15
      WRITE (6,FMT=99999) IFLAG,NDROP
      GO TO 70
   10 IFLAG = 0
      MEPS = EPS
      NLP = LP
      IF (N.GT.0) GO TO 20
      IFLAG = -11
      IF (LP.NE.0) WRITE (LP,FMT=99998) N
      GO TO 60
   20 IF (NZ.GT.0) GO TO 30
      IFLAG = -10
      IF (LP.NE.0) WRITE (LP,FMT=99997) NZ
      GO TO 60
   30 IF (LICN.GE.NZ) GO TO 40
      IFLAG = -9
      IF (LP.NE.0) WRITE (LP,FMT=99996) LICN
      GO TO 60
   40 CALL MA28DD(N,A,LICN,IVECT,JVECT,NZ,ICN,IKEEP,IKEEP(1,4),
     +            IKEEP(1,5),IKEEP(1,2),IKEEP(1,3),IW(1,3),IW,W(1),
     +            IFLAG)
      THEMAX = W(1)
      IF (LBIG) BIG1 = THEMAX
      IDUP = 0
      IF (IFLAG.EQ. (N+1)) IDUP = 1
      IF (IFLAG.LT.0) THEN
        IF (LP.NE.0) WRITE (LP,FMT=99994)
        GO TO 60
      END IF
      CALL MA30BD(N,ICN,A,LICN,IKEEP,IKEEP(1,4),IDISP,IKEEP(1,2),
     +            IKEEP(1,3),W,IW,IFLAG)
      IF (LBIG) BIG1 = BIG
      RMIN = MRMIN
      IF (IFLAG.GE.0) GO TO 50
      IFLAG = -2
      IF (LP.NE.0) WRITE (LP,FMT=99995)
      GO TO 60
   50 I1 = IDISP(1)
      IEND = LICN - I1 + 1
      IF (GROW) CALL MC24AD(N,ICN,A(I1),IEND,IKEEP,IKEEP(1,4),W)
      IF (GROW) W(1) = W(1) + THEMAX
      IF (GROW .AND. N.GT.1) W(2) = THEMAX
      IF (IDUP.EQ.1 .AND. IFLAG.GE.0) IFLAG = -14
      GO TO 70
   60 CONTINUE
   70 RETURN
99999 FORMAT (' ERROR RETURN FROM MA28B/BD WITH IFLAG=',I4,/,I7,
     +       ' ENTRIES DROPPED FROM STRUCTURE BY MA28A/AD')
99998 FORMAT (' ERROR RETURN FROM MA28B/BD BECAUSE ',
     +       'N OUT OF RANGE = ',I10)
99997 FORMAT (' ERROR RETURN FROM MA28B/BD BECAUSE ',
     +       'NZ NON POSITIVE = ',I10)
99996 FORMAT (' ERROR RETURN FROM MA28B/BD BECAUSE ',
     +       'LICN TOO SMALL = ',I10)
99995 FORMAT (' ERROR RETURN FROM MA28B/BD BECAUSE ',
     +       'ERROR RETURN FROM MA30B/BD')
99994 FORMAT (' ERROR RETURN FROM MA28B/BD ')
      END
      SUBROUTINE MA28CD(N,A,LICN,ICN,IKEEP,RHS,W,MTYPE)
      INTEGER LICN,MTYPE,N
      DOUBLE PRECISION A(LICN),RHS(N),W(N)
      INTEGER ICN(LICN),IKEEP(N,5)
      EXTERNAL MA30CD
      COMMON /MA28FD/EPS,RMIN,RESID,IRNCP,ICNCP,MINIRN,MINICN,IRANK,
     +       ABORT1,ABORT2
      COMMON /MA28GD/IDISP(2)
      COMMON /MA30HD/MRESID
      DOUBLE PRECISION EPS,MRESID,RESID,RMIN
      INTEGER ICNCP,IRANK,IRNCP,MINICN,MINIRN
      LOGICAL ABORT1,ABORT2
      INTEGER IDISP
      SAVE /MA28FD/,/MA28GD/,/MA30HD/
	
      CALL MA30CD(N,ICN,A,LICN,IKEEP,IKEEP(1,4),IKEEP(1,5),IDISP,
     +            IKEEP(1,2),IKEEP(1,3),RHS,W,MTYPE)
      
	RESID = MRESID
	
      RETURN
      END
      SUBROUTINE MA28DD(N,A,LICN,IVECT,JVECT,NZ,ICN,LENR,LENRL,LENOFF,
     +                  IP,IQ,IW1,IW,W1,IFLAG)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION W1
      INTEGER IFLAG,LICN,N,NZ
      DOUBLE PRECISION A(LICN)
      INTEGER ICN(LICN),IP(N),IQ(N),IVECT(NZ),IW(N,2),IW1(N,3),
     +        JVECT(NZ),LENOFF(N),LENR(N),LENRL(N)
      DOUBLE PRECISION AA
      INTEGER I,IBLOCK,IDISP2,IDUMMY,II,INEW,IOLD,J1,J2,JCOMP,JDUMMY,JJ,
     +        JNEW,JOLD,MIDPT
      LOGICAL BLOCKL
      INTRINSIC DABS,DMAX1,IABS
      COMMON /MA28ED/LP,MP,LBLOCK,GROW
      COMMON /MA28GD/IDISP(2)
      INTEGER LP,MP
      LOGICAL GROW,LBLOCK
      INTEGER IDISP
      SAVE /MA28ED/,/MA28GD/
      BLOCKL = LENOFF(1) .GE. 0
      IBLOCK = 1
      IW(1,1) = 1
      IW(1,2) = IDISP(1)
      DO 10 I = 1,N
        IW1(I,3) = IBLOCK
        IF (IP(I).LT.0) IBLOCK = IBLOCK + 1
        II = IABS(IP(I)+0)
        IW1(II,1) = I
        JJ = IQ(I)
        JJ = IABS(JJ)
        IW1(JJ,2) = I
        IF (I.EQ.1) GO TO 10
        IF (BLOCKL) IW(I,1) = IW(I-1,1) + LENOFF(I-1)
        IW(I,2) = IW(I-1,2) + LENR(I-1)
   10 CONTINUE
      IDISP2 = IDISP(2)
      DO 170 I = 1,NZ
        IF (I.GT.IDISP2) GO TO 20
        IF (ICN(I).LT.0) GO TO 170
   20   IOLD = IVECT(I)
        JOLD = JVECT(I)
        AA = A(I)
        DO 140 IDUMMY = 1,NZ
          IF (IOLD.LE.N .AND. IOLD.GT.0 .AND. JOLD.LE.N .AND.
     +        JOLD.GT.0) GO TO 30
          IF (LP.NE.0) WRITE (LP,FMT=99999) I,A(I),IOLD,JOLD
          IFLAG = -12
          GO TO 180
   30     INEW = IW1(IOLD,1)
          JNEW = IW1(JOLD,2)
          IF (IW1(INEW,3)-IW1(JNEW,3)) 40,60,50
   40     IFLAG = -13
          IF (LP.NE.0) WRITE (LP,FMT=99998) IOLD,JOLD
          GO TO 180
   50     J1 = IW(INEW,1)
          J2 = J1 + LENOFF(INEW) - 1
          GO TO 110
   60     J1 = IW(INEW,2)
          IF (INEW.GT.JNEW) GO TO 70
          J2 = J1 + LENR(INEW) - 1
          J1 = J1 + LENRL(INEW)
          GO TO 110
   70     J2 = J1 + LENRL(INEW)
          DO 100 JDUMMY = 1,N
            MIDPT = (J1+J2)/2
            JCOMP = IABS(ICN(MIDPT)+0)
            IF (JNEW-JCOMP) 80,130,90
   80       J2 = MIDPT
            GO TO 100
   90       J1 = MIDPT
  100     CONTINUE
          IFLAG = -13
          IF (LP.NE.0) WRITE (LP,FMT=99997) IOLD,JOLD
          GO TO 180
  110     DO 120 MIDPT = J1,J2
            IF (IABS(ICN(MIDPT)+0).EQ.JNEW) GO TO 130
  120     CONTINUE
          IFLAG = -13
          IF (LP.NE.0) WRITE (LP,FMT=99997) IOLD,JOLD
          GO TO 180
  130     IF (ICN(MIDPT).LT.0) GO TO 160
          IF (MIDPT.GT.NZ .OR. MIDPT.LE.I) GO TO 150
          W1 = A(MIDPT)
          A(MIDPT) = AA
          AA = W1
          IOLD = IVECT(MIDPT)
          JOLD = JVECT(MIDPT)
          ICN(MIDPT) = -ICN(MIDPT)
  140   CONTINUE
  150   A(MIDPT) = AA
        ICN(MIDPT) = -ICN(MIDPT)
        GO TO 170
  160   A(MIDPT) = A(MIDPT) + AA
        IFLAG = N + 1
  170 CONTINUE
  180 W1 = ZERO
      DO 200 I = 1,IDISP2
        IF (ICN(I).LT.0) GO TO 190
        A(I) = ZERO
        GO TO 200
  190   ICN(I) = -ICN(I)
        W1 = DMAX1(W1,DABS(A(I)))
  200 CONTINUE
      RETURN
99999 FORMAT (' ELEMENT ',I6,' WITH VALUE ',1P,D12.4,/,10X,
     +       ' HAS INDICES ',I8,' ,',I8,' INDICES OUT OF RANGE')
99998 FORMAT (' NON-ZERO',I7,' ,',I6,' IN ZERO OFF-DIAGONAL BLOCK')
99997 FORMAT (' ELEMENT',I6,' ,',I6,' WAS NOT IN L/U PATTERN')
      END
      SUBROUTINE MA28ID(N,NZ,AORG,IRNORG,ICNORG,LICN,A,ICN,IKEEP,RHS,X,
     +                  R,W,MTYPE,PREC,IFLAG)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION PREC
      INTEGER IFLAG,LICN,MTYPE,N,NZ
      DOUBLE PRECISION A(LICN),AORG(NZ),R(N),RHS(N),W(N),X(N)
      INTEGER ICN(LICN),ICNORG(NZ),IKEEP(N,5),IRNORG(NZ)
      DOUBLE PRECISION CONVER,D,DD
      INTEGER I,ITERAT,NCOL,NROW
      EXTERNAL MA28CD
      INTRINSIC DABS,DMAX1
      COMMON /MA28ED/LP,MP,LBLOCK,GROW
      COMMON /MA28HD/TOL,THEMAX,BIG,DXMAX,ERRMAX,DRES,CGCE,NDROP,MAXIT,
     +       NOITER,NSRCH,ISTART,LBIG
      DOUBLE PRECISION BIG,CGCE,DRES,DXMAX,ERRMAX,THEMAX,TOL
      INTEGER ISTART,LP,MAXIT,MP,NDROP,NOITER,NSRCH
      LOGICAL GROW,LBIG,LBLOCK
      SAVE /MA28ED/,/MA28HD/
      NOITER = 0
      ERRMAX = ZERO
      DRES = ZERO
      IFLAG = 0
      IF (ISTART.EQ.1) GO TO 20
      DO 10 I = 1,N
        X(I) = RHS(I)
   10 CONTINUE
      CALL MA28CD(N,A,LICN,ICN,IKEEP,X,W,MTYPE)
   20 IF (MAXIT.EQ.0) GO TO 160
      DD = 0.0
      DO 30 I = 1,N
        DD = DMAX1(DD,DABS(X(I)))
   30 CONTINUE
      DXMAX = DD
      DO 120 ITERAT = 1,MAXIT
        D = DD
        DO 40 I = 1,N
          R(I) = RHS(I)
   40   CONTINUE
        IF (MTYPE.EQ.1) GO TO 60
        DO 50 I = 1,NZ
          NROW = IRNORG(I)
          NCOL = ICNORG(I)
          R(NCOL) = R(NCOL) - AORG(I)*X(NROW)
   50   CONTINUE
        GO TO 80
   60   DO 70 I = 1,NZ
          NROW = IRNORG(I)
          NCOL = ICNORG(I)
          R(NROW) = R(NROW) - AORG(I)*X(NCOL)
   70   CONTINUE
   80   DRES = 0.0
        DO 90 I = 1,N
          DRES = DMAX1(DRES,DABS(R(I)))
   90   CONTINUE
        IF (DRES.EQ.0.0) GO TO 150
        NOITER = NOITER + 1
        CALL MA28CD(N,A,LICN,ICN,IKEEP,R,W,MTYPE)
        DD = 0.0
        DO 100 I = 1,N
          DD = DMAX1(DD,DABS(R(I)))
  100   CONTINUE
        IF (DD.GT.D*CGCE .AND. ITERAT.GE.2) GO TO 130
        IF (DXMAX*10.0+DD.EQ.DXMAX*10.0) GO TO 140
        DXMAX = 0.0
        DO 110 I = 1,N
          X(I) = X(I) + R(I)
          DXMAX = DMAX1(DXMAX,DABS(X(I)))
  110   CONTINUE
        IF (DD.LT.PREC*DXMAX) GO TO 140
  120 CONTINUE
      IFLAG = -16
      WRITE (LP,FMT=99999) IFLAG,MAXIT
      GO TO 140
  130 IFLAG = -17
      CONVER = DD/D
      WRITE (LP,FMT=99998) IFLAG,CONVER,CGCE
  140 ERRMAX = DD
  150 CONTINUE
  160 RETURN
99999 FORMAT (' ERROR RETURN FROM MA28I/ID WITH IFLAG = ',I3,/,' MORE ',
     +       'THAN',I5,' ITERATIONS REQUIRED')
99998 FORMAT (' ERROR RETURN FROM MA28I WITH IFLAG = ',I3,/,' CONVERGE',
     +       'NCE RATE OF ',1P,E9.2,' TOO SLOW',/,
     +       ' MAXIMUM ACCEPTABLE RATE',' SET TO ',1P,E9.2)
      END
      BLOCK DATA MA28JD
      COMMON /MA28ED/LP,MP,LBLOCK,GROW
      COMMON /MA28FD/EPS,RMIN,RESID,IRNCP,ICNCP,MINIRN,MINICN,IRANK,
     +       ABORT1,ABORT2
      COMMON /MA28HD/TOL,THEMAX,BIG,DXMAX,ERRMAX,DRES,CGCE,NDROP,MAXIT,
     +       NOITER,NSRCH,ISTART,LBIG
      DOUBLE PRECISION BIG,CGCE,DRES,DXMAX,EPS,ERRMAX,RESID,RMIN,THEMAX,
     +                 TOL
      INTEGER ICNCP,IRANK,IRNCP,ISTART,LP,MAXIT,MINICN,MINIRN,MP,NDROP,
     +        NOITER,NSRCH
      LOGICAL ABORT1,ABORT2,GROW,LBIG,LBLOCK
      SAVE /MA28ED/,/MA28FD/,/MA28HD/
      DATA EPS/1.0D-4/,TOL/0.0D0/,CGCE/0.5D0/
      DATA MAXIT/16/
      DATA LP/6/,MP/6/,NSRCH/32768/,ISTART/0/
      DATA LBLOCK/.TRUE./,GROW/.TRUE./,LBIG/.FALSE./
      DATA ABORT1/.TRUE./,ABORT2/.TRUE./
      END
