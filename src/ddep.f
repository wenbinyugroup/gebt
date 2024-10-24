* *******************************************************************
* COPYRIGHT (c) 1967 Hyprotech UK
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
*######DATE 4 Oct 1992
C       Toolpack tool decs employed.
C       SAVE statement for COMMON FA01ED added.
C  EAT 21/6/93 EXTERNAL statement put in for block data on VAXs.
C
C
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
*######DATE 14 Jan 1993
C       Toolpack tool decs employed.
C       Reference MA30JD removed.
C       SAVE statements added.
C       ZERO, ONE and UMAX made PARAMETER.
C
C  EAT 21/6/93 EXTERNAL statement put in for block data on VAXs.
C   3/1/96. LENOFF made into an assumed-size array.
C
C
      SUBROUTINE MA30AD(NN,ICN,A,LICN,LENR,LENRL,IDISP,IP,IQ,IRN,LIRN,
     +                  LENC,IFIRST,LASTR,NEXTR,LASTC,NEXTC,IPTR,IPC,U,
     +                  IFLAG)
      DOUBLE PRECISION ZERO,UMAX
      PARAMETER (ZERO=0.0D0,UMAX=.999999999D0)
      DOUBLE PRECISION U
      INTEGER IFLAG,LICN,LIRN,NN
      DOUBLE PRECISION A(LICN)
      INTEGER ICN(LICN),IDISP(2),IFIRST(NN),IP(NN),IPC(NN),IPTR(NN),
     +        IQ(NN),IRN(LIRN),LASTC(NN),LASTR(NN),LENC(NN),LENR(NN),
     +        LENRL(NN),NEXTC(NN),NEXTR(NN)
      DOUBLE PRECISION AANEW,AMAX,ANEW,AU,PIVR,PIVRAT,SCALE
      INTEGER COLUPD,DISPC,I,I1,I2,IACTIV,IBEG,IDISPC,IDROP,IDUMMY,IEND,
     +        IFILL,IFIR,II,III,IJFIR,IJP1,IJPOS,ILAST,INDROW,IOP,IPIV,
     +        IPOS,IROWS,ISING,ISRCH,ISTART,ISW,ISW1,ITOP,J,J1,J2,JBEG,
     +        JCOST,JCOUNT,JDIFF,JDUMMY,JEND,JJ,JMORE,JNEW,JNPOS,JOLD,
     +        JPIV,JPOS,JROOM,JVAL,JZER,JZERO,K,KCOST,KDROP,L,LC,LENPIV,
     +        LL,LR,MOREI,MSRCH,N,NBLOCK,NC,NNM1,NR,NUM,NZ,NZ2,NZCOL,
     +        NZMIN,NZPC,NZROW,OLDEND,OLDPIV,PIVEND,PIVOT,PIVROW,ROWI
      EXTERNAL MA30DD
      EXTERNAL MA30JD
      INTRINSIC DABS,DMAX1,DMIN1,IABS,MAX0,MIN0
      COMMON /MA30ED/LP,ABORT1,ABORT2,ABORT3
      COMMON /MA30FD/IRNCP,ICNCP,IRANK,MINIRN,MINICN
      COMMON /MA30ID/TOL,BIG,NDROP,NSRCH,LBIG
      DOUBLE PRECISION BIG,TOL
      INTEGER ICNCP,IRANK,IRNCP,LP,MINICN,MINIRN,NDROP,NSRCH
      LOGICAL ABORT1,ABORT2,ABORT3,LBIG
      SAVE /MA30ED/,/MA30FD/,/MA30ID/
      MSRCH = NSRCH
      NDROP = 0
      MINIRN = 0
      MINICN = IDISP(1) - 1
      MOREI = 0
      IRANK = NN
      IRNCP = 0
      ICNCP = 0
      IFLAG = 0
      U = DMIN1(U,UMAX)
      U = DMAX1(U,ZERO)
      IBEG = IDISP(1)
      IACTIV = IDISP(2)
      NZROW = LICN - IACTIV + 1
      MINICN = NZROW + MINICN
      NUM = 1
      IPTR(1) = IACTIV
      IF (NN.EQ.1) GO TO 20
      NNM1 = NN - 1
      DO 10 I = 1,NNM1
        IF (IP(I).LT.0) NUM = NUM + 1
        IPTR(I+1) = IPTR(I) + LENR(I)
   10 CONTINUE
   20 ILAST = 0
      DO 1000 NBLOCK = 1,NUM
        ISTART = ILAST + 1
        DO 30 IROWS = ISTART,NN
          IF (IP(IROWS).LT.0) GO TO 40
   30   CONTINUE
        IROWS = NN
   40   ILAST = IROWS
        N = ILAST - ISTART + 1
        IF (N.NE.1) GO TO 90
        LENRL(ILAST) = 0
        ISING = ISTART
        IF (LENR(ILAST).NE.0) GO TO 50
        IRANK = IRANK - 1
        ISING = -ISING
        IF (IFLAG.NE.2 .AND. IFLAG.NE.-5) IFLAG = 1
        IF (.NOT.ABORT1) GO TO 80
        IDISP(2) = IACTIV
        IFLAG = -1
        IF (LP.NE.0) WRITE (LP,FMT=99999)
        GO TO 1120
   50   SCALE = DABS(A(IACTIV))
        IF (SCALE.EQ.ZERO) GO TO 60
        IF (LBIG) BIG = DMAX1(BIG,SCALE)
        GO TO 70
   60   ISING = -ISING
        IRANK = IRANK - 1
        IPTR(ILAST) = 0
        IF (IFLAG.NE.-5) IFLAG = 2
        IF (.NOT.ABORT2) GO TO 70
        IDISP(2) = IACTIV
        IFLAG = -2
        IF (LP.NE.0) WRITE (LP,FMT=99998)
        GO TO 1120
   70   A(IBEG) = A(IACTIV)
        ICN(IBEG) = ICN(IACTIV)
        IACTIV = IACTIV + 1
        IPTR(ISTART) = 0
        IBEG = IBEG + 1
        NZROW = NZROW - 1
   80   LASTR(ISTART) = ISTART
        IPC(ISTART) = -ISING
        GO TO 1000
   90   ITOP = LICN
        IF (ILAST.NE.NN) ITOP = IPTR(ILAST+1) - 1
        DO 100 I = ISTART,ILAST
          LENRL(I) = 0
          LENC(I) = 0
  100   CONTINUE
        IF (ITOP-IACTIV.LT.LIRN) GO TO 110
        MINIRN = ITOP - IACTIV + 1
        PIVOT = ISTART - 1
        GO TO 1100
  110   DO 120 II = IACTIV,ITOP
          I = ICN(II)
          LENC(I) = LENC(I) + 1
  120   CONTINUE
        IPC(ILAST) = LIRN + 1
        J1 = ISTART + 1
        DO 130 JJ = J1,ILAST
          J = ILAST - JJ + J1 - 1
          IPC(J) = IPC(J+1) - LENC(J+1)
  130   CONTINUE
        DO 150 INDROW = ISTART,ILAST
          J1 = IPTR(INDROW)
          J2 = J1 + LENR(INDROW) - 1
          IF (J1.GT.J2) GO TO 150
          DO 140 JJ = J1,J2
            J = ICN(JJ)
            IPOS = IPC(J) - 1
            IRN(IPOS) = INDROW
            IPC(J) = IPOS
  140     CONTINUE
  150   CONTINUE
        DISPC = IPC(ISTART)
        NZCOL = LIRN - DISPC + 1
        MINIRN = MAX0(NZCOL,MINIRN)
        NZMIN = 1
        DO 160 I = 1,N
          IFIRST(I) = 0
  160   CONTINUE
        DO 180 JJ = ISTART,ILAST
          J = ILAST - JJ + ISTART
          NZ = LENC(J)
          IF (NZ.NE.0) GO TO 170
          IPC(J) = 0
          GO TO 180
  170     IF (NSRCH.LE.NN) GO TO 180
          ISW = IFIRST(NZ)
          IFIRST(NZ) = -J
          LASTC(J) = 0
          NEXTC(J) = -ISW
          ISW1 = IABS(ISW)
          IF (ISW.NE.0) LASTC(ISW1) = J
  180   CONTINUE
        DO 210 II = ISTART,ILAST
          I = ILAST - II + ISTART
          NZ = LENR(I)
          IF (NZ.NE.0) GO TO 190
          IPTR(I) = 0
          LASTR(I) = 0
          GO TO 210
  190     ISW = IFIRST(NZ)
          IFIRST(NZ) = I
          IF (ISW.GT.0) GO TO 200
          NEXTR(I) = 0
          LASTR(I) = ISW
          GO TO 210
  200     NEXTR(I) = ISW
          LASTR(I) = LASTR(ISW)
          LASTR(ISW) = I
  210   CONTINUE
        DO 980 PIVOT = ISTART,ILAST
          NZ2 = NZMIN
          JCOST = HUGE(N)
          DO 340 L = 1,2
            PIVRAT = ZERO
            ISRCH = 1
            LL = L
            DO 330 NZ = NZ2,N
              IF (JCOST.LE. (NZ-1)**2) GO TO 420
              IJFIR = IFIRST(NZ)
              IF (IJFIR) 230,220,240
  220         IF (LL.EQ.1) NZMIN = NZ + 1
              GO TO 330
  230         LL = 2
              IJFIR = -IJFIR
              GO TO 290
  240         LL = 2
              DO 270 IDUMMY = 1,N
                IF (JCOST.LE. (NZ-1)**2) GO TO 420
                IF (ISRCH.GT.MSRCH) GO TO 420
                IF (IJFIR.EQ.0) GO TO 280
                I = IJFIR
                IJFIR = NEXTR(I)
                AMAX = ZERO
                J1 = IPTR(I) + LENRL(I)
                J2 = IPTR(I) + LENR(I) - 1
                DO 250 JJ = J1,J2
                  AMAX = DMAX1(AMAX,DABS(A(JJ)))
  250           CONTINUE
                AU = AMAX*U
                ISRCH = ISRCH + 1
                DO 260 JJ = J1,J2
                  IF (DABS(A(JJ)).LE.AU .AND. L.EQ.1) GO TO 260
                  J = ICN(JJ)
                  KCOST = (NZ-1)* (LENC(J)-1)
                  IF (KCOST.GT.JCOST) GO TO 260
                  PIVR = ZERO
                  IF (AMAX.NE.ZERO) PIVR = DABS(A(JJ))/AMAX
                  IF (KCOST.EQ.JCOST .AND. (PIVR.LE.PIVRAT.OR.
     +                NSRCH.GT.NN+1)) GO TO 260
                  JCOST = KCOST
                  IJPOS = JJ
                  IPIV = I
                  JPIV = J
                  IF (MSRCH.GT.NN+1 .AND. JCOST.LE. (NZ-1)**2) GO TO 420
                  PIVRAT = PIVR
  260           CONTINUE
  270         CONTINUE
  280         IJFIR = IFIRST(NZ)
              IJFIR = -LASTR(IJFIR)
  290         IF (JCOST.LE.NZ* (NZ-1)) GO TO 420
              IF (MSRCH.LE.NN) GO TO 330
              DO 320 IDUMMY = 1,N
                IF (IJFIR.EQ.0) GO TO 330
                J = IJFIR
                IJFIR = NEXTC(IJFIR)
                I1 = IPC(J)
                I2 = I1 + NZ - 1
                DO 310 II = I1,I2
                  I = IRN(II)
                  KCOST = (NZ-1)* (LENR(I)-LENRL(I)-1)
                  IF (KCOST.GE.JCOST) GO TO 310
                  J1 = IPTR(I) + LENRL(I)
                  J2 = IPTR(I) + LENR(I) - 1
                  AMAX = ZERO
                  DO 300 JJ = J1,J2
                    AMAX = DMAX1(AMAX,DABS(A(JJ)))
                    IF (ICN(JJ).EQ.J) JPOS = JJ
  300             CONTINUE
                  IF (DABS(A(JPOS)).LE.AMAX*U .AND. L.EQ.1) GO TO 310
                  JCOST = KCOST
                  IPIV = I
                  JPIV = J
                  IJPOS = JPOS
                  IF (AMAX.NE.ZERO) PIVRAT = DABS(A(JPOS))/AMAX
                  IF (JCOST.LE.NZ* (NZ-1)) GO TO 420
  310           CONTINUE
  320         CONTINUE
  330       CONTINUE
            MSRCH = N
            IRANK = IRANK - 1
  340     CONTINUE
          IF (IFLAG.NE.2 .AND. IFLAG.NE.-5) IFLAG = 1
          IRANK = IRANK - ILAST + PIVOT + 1
          IF (.NOT.ABORT1) GO TO 350
          IDISP(2) = IACTIV
          IFLAG = -1
          IF (LP.NE.0) WRITE (LP,FMT=99999)
          GO TO 1120
  350     K = PIVOT - 1
          DO 390 I = ISTART,ILAST
            IF (LASTR(I).NE.0) GO TO 390
            K = K + 1
            LASTR(I) = K
            IF (LENRL(I).EQ.0) GO TO 380
            MINICN = MAX0(MINICN,NZROW+IBEG-1+MOREI+LENRL(I))
            IF (IACTIV-IBEG.GE.LENRL(I)) GO TO 360
            CALL MA30DD(A,ICN,IPTR(ISTART),N,IACTIV,ITOP,.TRUE.)
            IF (IACTIV-IBEG.GE.LENRL(I)) GO TO 360
            MOREI = MOREI + IBEG - IDISP(1)
            IBEG = IDISP(1)
            IF (LP.NE.0) WRITE (LP,FMT=99997)
            IFLAG = -5
            IF (ABORT3) GO TO 1090
  360       J1 = IPTR(I)
            J2 = J1 + LENRL(I) - 1
            IPTR(I) = 0
            DO 370 JJ = J1,J2
              A(IBEG) = A(JJ)
              ICN(IBEG) = ICN(JJ)
              ICN(JJ) = 0
              IBEG = IBEG + 1
  370       CONTINUE
            NZROW = NZROW - LENRL(I)
  380       IF (K.EQ.ILAST) GO TO 400
  390     CONTINUE
  400     K = PIVOT - 1
          DO 410 I = ISTART,ILAST
            IF (IPC(I).NE.0) GO TO 410
            K = K + 1
            IPC(I) = K
            IF (K.EQ.ILAST) GO TO 990
  410     CONTINUE
  420     ISING = PIVOT
          IF (A(IJPOS).NE.ZERO) GO TO 430
          ISING = -ISING
          IF (IFLAG.NE.-5) IFLAG = 2
          IF (.NOT.ABORT2) GO TO 430
          IDISP(2) = IACTIV
          IFLAG = -2
          IF (LP.NE.0) WRITE (LP,FMT=99998)
          GO TO 1120
  430     OLDPIV = IPTR(IPIV) + LENRL(IPIV)
          OLDEND = IPTR(IPIV) + LENR(IPIV) - 1
          IF (NSRCH.LE.NN) GO TO 460
          COLUPD = NN + 1
          DO 450 JJ = OLDPIV,OLDEND
            J = ICN(JJ)
            LC = LASTC(J)
            NC = NEXTC(J)
            NEXTC(J) = -COLUPD
            IF (JJ.NE.IJPOS) COLUPD = J
            IF (NC.NE.0) LASTC(NC) = LC
            IF (LC.EQ.0) GO TO 440
            NEXTC(LC) = NC
            GO TO 450
  440       NZ = LENC(J)
            ISW = IFIRST(NZ)
            IF (ISW.GT.0) LASTR(ISW) = -NC
            IF (ISW.LT.0) IFIRST(NZ) = -NC
  450     CONTINUE
  460     I1 = IPC(JPIV)
          I2 = I1 + LENC(JPIV) - 1
          DO 480 II = I1,I2
            I = IRN(II)
            LR = LASTR(I)
            NR = NEXTR(I)
            IF (NR.NE.0) LASTR(NR) = LR
            IF (LR.LE.0) GO TO 470
            NEXTR(LR) = NR
            GO TO 480
  470       NZ = LENR(I) - LENRL(I)
            IF (NR.NE.0) IFIRST(NZ) = NR
            IF (NR.EQ.0) IFIRST(NZ) = LR
  480     CONTINUE
          IF (OLDPIV.EQ.IJPOS) GO TO 490
          AU = A(OLDPIV)
          A(OLDPIV) = A(IJPOS)
          A(IJPOS) = AU
          ICN(IJPOS) = ICN(OLDPIV)
          ICN(OLDPIV) = JPIV
  490     MINICN = MAX0(MINICN,NZROW+IBEG-1+MOREI+LENR(IPIV))
          IF (IACTIV-IBEG.GE.LENR(IPIV)) GO TO 500
          CALL MA30DD(A,ICN,IPTR(ISTART),N,IACTIV,ITOP,.TRUE.)
          OLDPIV = IPTR(IPIV) + LENRL(IPIV)
          OLDEND = IPTR(IPIV) + LENR(IPIV) - 1
          IF (IACTIV-IBEG.GE.LENR(IPIV)) GO TO 500
          MOREI = MOREI + IBEG - IDISP(1)
          IBEG = IDISP(1)
          IF (LP.NE.0) WRITE (LP,FMT=99997)
          IFLAG = -5
          IF (ABORT3) GO TO 1090
          IF (IACTIV-IBEG.GE.LENR(IPIV)) GO TO 500
          IFLAG = -4
          GO TO 1090
  500     IJPOS = 0
          J1 = IPTR(IPIV)
          DO 530 JJ = J1,OLDEND
            A(IBEG) = A(JJ)
            ICN(IBEG) = ICN(JJ)
            IF (IJPOS.NE.0) GO TO 510
            IF (ICN(JJ).EQ.JPIV) IJPOS = IBEG
            ICN(JJ) = 0
            GO TO 520
  510       K = IBEG - IJPOS
            J = ICN(JJ)
            ICN(JJ) = IQ(J)
            IQ(J) = -K
  520       IBEG = IBEG + 1
  530     CONTINUE
          IJP1 = IJPOS + 1
          PIVEND = IBEG - 1
          LENPIV = PIVEND - IJPOS
          NZROW = NZROW - LENRL(IPIV) - 1
          IPTR(IPIV) = OLDPIV + 1
          IF (LENPIV.EQ.0) IPTR(IPIV) = 0
          DO 560 JJ = IJPOS,PIVEND
            J = ICN(JJ)
            I1 = IPC(J)
            LENC(J) = LENC(J) - 1
            I2 = IPC(J) + LENC(J) - 1
            IF (I2.LT.I1) GO TO 550
            DO 540 II = I1,I2
              IF (IRN(II).NE.IPIV) GO TO 540
              IRN(II) = IRN(I2+1)
              GO TO 550
  540       CONTINUE
  550       IRN(I2+1) = 0
  560     CONTINUE
          NZCOL = NZCOL - LENPIV - 1
          NZPC = LENC(JPIV)
          IF (NZPC.EQ.0) GO TO 900
          DO 840 III = 1,NZPC
            II = IPC(JPIV) + III - 1
            I = IRN(II)
            IDROP = 0
            J1 = IPTR(I) + LENRL(I)
            IEND = IPTR(I) + LENR(I) - 1
            DO 570 JJ = J1,IEND
              IF (ICN(JJ).NE.JPIV) GO TO 570
              AU = ZERO
              IF (A(IJPOS).NE.ZERO) AU = -A(JJ)/A(IJPOS)
              IF (LBIG) BIG = DMAX1(BIG,DABS(AU))
              A(JJ) = A(J1)
              A(J1) = AU
              ICN(JJ) = ICN(J1)
              ICN(J1) = JPIV
              LENRL(I) = LENRL(I) + 1
              GO TO 580
  570       CONTINUE
  580       IF (LENPIV.EQ.0) GO TO 840
            ROWI = J1 + 1
            IOP = 0
            IF (ROWI.GT.IEND) GO TO 650
            DO 590 JJ = ROWI,IEND
              J = ICN(JJ)
              IF (IQ(J).GT.0) GO TO 590
              IOP = IOP + 1
              PIVROW = IJPOS - IQ(J)
              A(JJ) = A(JJ) + AU*A(PIVROW)
              IF (LBIG) BIG = DMAX1(DABS(A(JJ)),BIG)
              ICN(PIVROW) = -ICN(PIVROW)
              IF (DABS(A(JJ)).LT.TOL) IDROP = IDROP + 1
  590       CONTINUE
            IF (IDROP.EQ.0) GO TO 650
            JNEW = ROWI
            DO 630 JJ = ROWI,IEND
              IF (DABS(A(JJ)).LT.TOL) GO TO 600
              A(JNEW) = A(JJ)
              ICN(JNEW) = ICN(JJ)
              JNEW = JNEW + 1
              GO TO 630
  600         J = ICN(JJ)
              I1 = IPC(J)
              I2 = I1 + LENC(J) - 1
              DO 610 II = I1,I2
                IF (IRN(II).EQ.I) GO TO 620
  610         CONTINUE
  620         IRN(II) = IRN(I2)
              IRN(I2) = 0
              LENC(J) = LENC(J) - 1
              IF (NSRCH.LE.NN) GO TO 630
              IF (NEXTC(J).LT.0) GO TO 630
              LC = LASTC(J)
              NC = NEXTC(J)
              NEXTC(J) = -COLUPD
              COLUPD = J
              IF (NC.NE.0) LASTC(NC) = LC
              IF (LC.EQ.0) GO TO 622
              NEXTC(LC) = NC
              GO TO 630
  622         NZ = LENC(J) + 1
              ISW = IFIRST(NZ)
              IF (ISW.GT.0) LASTR(ISW) = -NC
              IF (ISW.LT.0) IFIRST(NZ) = -NC
  630       CONTINUE
            DO 640 JJ = JNEW,IEND
              ICN(JJ) = 0
  640       CONTINUE
            IDROP = IEND + 1 - JNEW
            IEND = JNEW - 1
            LENR(I) = LENR(I) - IDROP
            NZROW = NZROW - IDROP
            NZCOL = NZCOL - IDROP
            NDROP = NDROP + IDROP
  650       IFILL = LENPIV - IOP
            IF (IFILL.EQ.0) GO TO 750
            MINICN = MAX0(MINICN,MOREI+IBEG-1+NZROW+IFILL+LENR(I))
            DO 660 JDIFF = 1,IFILL
              JNPOS = IEND + JDIFF
              IF (JNPOS.GT.LICN) GO TO 670
              IF (ICN(JNPOS).NE.0) GO TO 670
  660       CONTINUE
            IEND = IEND + 1
            GO TO 750
  670       JMORE = IFILL - JDIFF + 1
            I1 = IPTR(I)
            DO 680 JDIFF = 1,JMORE
              JNPOS = I1 - JDIFF
              IF (JNPOS.LT.IACTIV) GO TO 690
              IF (ICN(JNPOS).NE.0) GO TO 700
  680       CONTINUE
  690       JNPOS = I1 - JMORE
            GO TO 710
  700       JNPOS = IACTIV - LENR(I) - IFILL
  710       IF (JNPOS.GE.IBEG) GO TO 730
            CALL MA30DD(A,ICN,IPTR(ISTART),N,IACTIV,ITOP,.TRUE.)
            I1 = IPTR(I)
            IEND = I1 + LENR(I) - 1
            JNPOS = IACTIV - LENR(I) - IFILL
            IF (JNPOS.GE.IBEG) GO TO 730
            MOREI = MOREI + IBEG - IDISP(1) - LENPIV - 1
            IF (LP.NE.0) WRITE (LP,FMT=99997)
            IFLAG = -5
            IF (ABORT3) GO TO 1090
            IBEG = IDISP(1)
            ICN(IBEG) = JPIV
            A(IBEG) = A(IJPOS)
            IJPOS = IBEG
            DO 720 JJ = IJP1,PIVEND
              IBEG = IBEG + 1
              A(IBEG) = A(JJ)
              ICN(IBEG) = ICN(JJ)
  720       CONTINUE
            IJP1 = IJPOS + 1
            PIVEND = IBEG
            IBEG = IBEG + 1
            IF (JNPOS.GE.IBEG) GO TO 730
            IFLAG = -4
            GO TO 1090
  730       IACTIV = MIN0(IACTIV,JNPOS)
            IPTR(I) = JNPOS
            DO 740 JJ = I1,IEND
              A(JNPOS) = A(JJ)
              ICN(JNPOS) = ICN(JJ)
              JNPOS = JNPOS + 1
              ICN(JJ) = 0
  740       CONTINUE
            IEND = JNPOS
  750       NZROW = NZROW + IFILL
            IDROP = 0
            DO 830 JJ = IJP1,PIVEND
              J = ICN(JJ)
              IF (J.LT.0) GO TO 820
              ANEW = AU*A(JJ)
              AANEW = DABS(ANEW)
              IF (AANEW.GE.TOL) GO TO 760
              IDROP = IDROP + 1
              NDROP = NDROP + 1
              NZROW = NZROW - 1
              MINICN = MINICN - 1
              IFILL = IFILL - 1
              GO TO 830
  760         IF (LBIG) BIG = DMAX1(AANEW,BIG)
              A(IEND) = ANEW
              ICN(IEND) = J
              IEND = IEND + 1
              MINIRN = MAX0(MINIRN,NZCOL+LENC(J)+1)
              JEND = IPC(J) + LENC(J)
              JROOM = NZPC - III + 1 + LENC(J)
              IF (JEND.GT.LIRN) GO TO 770
              IF (IRN(JEND).EQ.0) GO TO 810
  770         IF (JROOM.LT.DISPC) GO TO 780
              CALL MA30DD(A,IRN,IPC(ISTART),N,DISPC,LIRN,.FALSE.)
              IF (JROOM.LT.DISPC) GO TO 780
              JROOM = DISPC - 1
              IF (JROOM.GE.LENC(J)+1) GO TO 780
              GO TO 1100
  780         JBEG = IPC(J)
              JEND = IPC(J) + LENC(J) - 1
              JZERO = DISPC - 1
              DISPC = DISPC - JROOM
              IDISPC = DISPC
              DO 790 II = JBEG,JEND
                IRN(IDISPC) = IRN(II)
                IRN(II) = 0
                IDISPC = IDISPC + 1
  790         CONTINUE
              IPC(J) = DISPC
              JEND = IDISPC
              DO 800 II = JEND,JZERO
                IRN(II) = 0
  800         CONTINUE
  810         IRN(JEND) = I
              NZCOL = NZCOL + 1
              LENC(J) = LENC(J) + 1
              GO TO 830
  820         ICN(JJ) = -J
  830       CONTINUE
            IF (IDROP.EQ.0) GO TO 834
            DO 832 KDROP = 1,IDROP
              ICN(IEND) = 0
              IEND = IEND + 1
  832       CONTINUE
  834       LENR(I) = LENR(I) + IFILL
  840     CONTINUE
          I1 = IPC(JPIV)
          I2 = IPC(JPIV) + LENC(JPIV) - 1
          NZCOL = NZCOL - LENC(JPIV)
          DO 890 II = I1,I2
            I = IRN(II)
            IRN(II) = 0
            NZ = LENR(I) - LENRL(I)
            IF (NZ.NE.0) GO TO 850
            LASTR(I) = 0
            GO TO 890
  850       IFIR = IFIRST(NZ)
            IFIRST(NZ) = I
            IF (IFIR) 860,880,870
  860       LASTR(I) = IFIR
            NEXTR(I) = 0
            GO TO 890
  870       LASTR(I) = LASTR(IFIR)
            NEXTR(I) = IFIR
            LASTR(IFIR) = I
            GO TO 890
  880       LASTR(I) = 0
            NEXTR(I) = 0
            NZMIN = MIN0(NZMIN,NZ)
  890     CONTINUE
  900     IPC(JPIV) = -ISING
          LASTR(IPIV) = PIVOT
          IF (LENPIV.EQ.0) GO TO 980
          NZROW = NZROW - LENPIV
          JVAL = IJP1
          JZER = IPTR(IPIV)
          IPTR(IPIV) = 0
          DO 910 JCOUNT = 1,LENPIV
            J = ICN(JVAL)
            IQ(J) = ICN(JZER)
            ICN(JZER) = 0
            JVAL = JVAL + 1
            JZER = JZER + 1
  910     CONTINUE
          IF (NSRCH.GT.NN) GO TO 920
          DO 916 JJ = IJP1,PIVEND
            J = ICN(JJ)
            NZ = LENC(J)
            IF (NZ.NE.0) GO TO 914
            IPC(J) = 0
            GO TO 916
  914       NZMIN = MIN0(NZMIN,NZ)
  916     CONTINUE
          GO TO 980
  920     JJ = COLUPD
          DO 970 JDUMMY = 1,NN
            J = JJ
            IF (J.EQ.NN+1) GO TO 980
            JJ = -NEXTC(J)
            NZ = LENC(J)
            IF (NZ.NE.0) GO TO 924
            IPC(J) = 0
            GO TO 970
  924       IFIR = IFIRST(NZ)
            LASTC(J) = 0
            IF (IFIR) 930,940,950
  930       IFIRST(NZ) = -J
            IFIR = -IFIR
            LASTC(IFIR) = J
            NEXTC(J) = IFIR
            GO TO 970
  940       IFIRST(NZ) = -J
            NEXTC(J) = 0
            GO TO 960
  950       LC = -LASTR(IFIR)
            LASTR(IFIR) = -J
            NEXTC(J) = LC
            IF (LC.NE.0) LASTC(LC) = J
  960       NZMIN = MIN0(NZMIN,NZ)
  970     CONTINUE
  980   CONTINUE
  990   IF (ILAST.NE.NN) IACTIV = IPTR(ILAST+1)
 1000 CONTINUE
      IF (IRANK.EQ.NN) GO TO 1020
      DO 1010 I = 1,NN
        IF (IPC(I).LT.0) GO TO 1010
        ISING = IPC(I)
        IQ(ISING) = -IQ(ISING)
        IPC(I) = -ISING
 1010 CONTINUE
 1020 ISTART = IDISP(1)
      IEND = IBEG - 1
      IF (IEND.LT.ISTART) GO TO 1040
      DO 1030 JJ = ISTART,IEND
        JOLD = ICN(JJ)
        ICN(JJ) = -IPC(JOLD)
 1030 CONTINUE
 1040 DO 1050 II = 1,NN
        I = LASTR(II)
        NEXTR(I) = LENR(II)
        IPTR(I) = LENRL(II)
 1050 CONTINUE
      DO 1060 I = 1,NN
        LENRL(I) = IPTR(I)
        LENR(I) = NEXTR(I)
 1060 CONTINUE
      DO 1070 II = 1,NN
        I = LASTR(II)
        J = -IPC(II)
        NEXTR(I) = IABS(IP(II)+0)
        IPTR(J) = IABS(IQ(II)+0)
 1070 CONTINUE
      DO 1080 I = 1,NN
        IF (IP(I).LT.0) NEXTR(I) = -NEXTR(I)
        IP(I) = NEXTR(I)
        IF (IQ(I).LT.0) IPTR(I) = -IPTR(I)
        IQ(I) = IPTR(I)
 1080 CONTINUE
      IP(NN) = IABS(IP(NN)+0)
      IDISP(2) = IEND
      GO TO 1120
 1090 IDISP(2) = IACTIV
      IF (LP.EQ.0) GO TO 1120
      WRITE (LP,FMT=99996)
      GO TO 1110
 1100 IF (IFLAG.EQ.-5) IFLAG = -6
      IF (IFLAG.NE.-6) IFLAG = -3
      IDISP(2) = IACTIV
      IF (LP.EQ.0) GO TO 1120
      IF (IFLAG.EQ.-3) WRITE (LP,FMT=99995)
      IF (IFLAG.EQ.-6) WRITE (LP,FMT=99994)
 1110 PIVOT = PIVOT - ISTART + 1
      WRITE (LP,FMT=99993) PIVOT,NBLOCK,ISTART,ILAST
      IF (PIVOT.EQ.0) WRITE (LP,FMT=99992) MINIRN
 1120 RETURN
99999 FORMAT (' ERROR RETURN FROM MA30A/AD BECAUSE MATRIX IS STRUCTUR',
     +       'ALLY SINGULAR')
99998 FORMAT (' ERROR RETURN FROM MA30A/AD BECAUSE MATRIX IS NUMERICA',
     +       'LLY SINGULAR')
99997 FORMAT (' LU DECOMPOSITION DESTROYED TO CREATE MORE SPACE')
99996 FORMAT (' ERROR RETURN FROM MA30A/AD BECAUSE LICN NOT BIG ENOUG',
     +       'H')
99995 FORMAT (' ERROR RETURN FROM MA30A/AD BECAUSE LIRN NOT BIG ENOUG',
     +       'H')
99994 FORMAT (' ERROR RETURN FROM MA30A/AD LIRN AND LICN TOO SMALL')
99993 FORMAT (' AT STAGE ',I5,' IN BLOCK ',I5,' WITH FIRST ROW ',I5,
     +       ' AND LAST ROW ',I5)
99992 FORMAT (' TO CONTINUE SET LIRN TO AT LEAST ',I8)
      END
      SUBROUTINE MA30BD(N,ICN,A,LICN,LENR,LENRL,IDISP,IP,IQ,W,IW,IFLAG)
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      INTEGER IFLAG,LICN,N
      DOUBLE PRECISION A(LICN),W(N)
      INTEGER ICN(LICN),IDISP(2),IP(N),IQ(N),IW(N),LENR(N),LENRL(N)
      DOUBLE PRECISION AU,ROWMAX
      INTEGER I,IFIN,ILEND,IPIVJ,ISING,ISTART,J,JAY,JAYJAY,JFIN,JJ,
     +        PIVPOS
      LOGICAL STAB
      INTRINSIC DABS,DMAX1
      COMMON /MA30ED/LP,ABORT1,ABORT2,ABORT3
      COMMON /MA30GD/EPS,RMIN
      COMMON /MA30ID/TOL,BIG,NDROP,NSRCH,LBIG
      DOUBLE PRECISION BIG,EPS,RMIN,TOL
      INTEGER LP,NDROP,NSRCH
      LOGICAL ABORT1,ABORT2,ABORT3,LBIG
      SAVE /MA30ED/,/MA30GD/,/MA30ID/
      STAB = EPS .LE. ONE
      RMIN = EPS
      ISING = 0
      IFLAG = 0
      DO 10 I = 1,N
        W(I) = ZERO
   10 CONTINUE
      IW(1) = IDISP(1)
      IF (N.EQ.1) GO TO 25
      DO 20 I = 2,N
        IW(I) = IW(I-1) + LENR(I-1)
   20 CONTINUE
   25 DO 160 I = 1,N
        ISTART = IW(I)
        IFIN = ISTART + LENR(I) - 1
        ILEND = ISTART + LENRL(I) - 1
        IF (ISTART.GT.ILEND) GO TO 90
        DO 30 JJ = ISTART,IFIN
          J = ICN(JJ)
          W(J) = A(JJ)
   30   CONTINUE
        DO 70 JJ = ISTART,ILEND
          J = ICN(JJ)
          IPIVJ = IW(J) + LENRL(J)
          AU = -W(J)/A(IPIVJ)
          IF (LBIG) BIG = DMAX1(DABS(AU),BIG)
          W(J) = AU
          IPIVJ = IPIVJ + 1
          JFIN = IW(J) + LENR(J) - 1
          IF (IPIVJ.GT.JFIN) GO TO 70
          IF (LBIG) GO TO 50
          DO 40 JAYJAY = IPIVJ,JFIN
            JAY = ICN(JAYJAY)
            W(JAY) = W(JAY) + AU*A(JAYJAY)
   40     CONTINUE
          GO TO 70
   50     DO 60 JAYJAY = IPIVJ,JFIN
            JAY = ICN(JAYJAY)
            W(JAY) = W(JAY) + AU*A(JAYJAY)
            BIG = DMAX1(DABS(W(JAY)),BIG)
   60     CONTINUE
   70   CONTINUE
        DO 80 JJ = ISTART,IFIN
          J = ICN(JJ)
          A(JJ) = W(J)
          W(J) = ZERO
   80   CONTINUE
   90   PIVPOS = ILEND + 1
        IF (IQ(I).GT.0) GO TO 140
        IF (ISING.EQ.0) ISING = I
        IF (PIVPOS.GT.IFIN) GO TO 100
        IF (A(PIVPOS).NE.ZERO) GO TO 170
  100   IF (ISTART.GT.IFIN) GO TO 120
        DO 110 JJ = ISTART,IFIN
          IF (ICN(JJ).LT.ISING) GO TO 110
          IF (A(JJ).NE.ZERO) GO TO 170
  110   CONTINUE
  120   IF (PIVPOS.LE.IFIN) A(PIVPOS) = ONE
        IF (IP(I).GT.0 .AND. I.NE.N) GO TO 160
        DO 130 J = ISING,I
          IF ((LENR(J)-LENRL(J)).EQ.0) GO TO 130
          JJ = IW(J) + LENRL(J)
          A(JJ) = ZERO
  130   CONTINUE
        ISING = 0
        GO TO 160
  140   IF (PIVPOS.GT.IFIN) GO TO 170
        IF (A(PIVPOS).EQ.ZERO) GO TO 170
        IF (.NOT.STAB) GO TO 160
        ROWMAX = ZERO
        DO 150 JJ = PIVPOS,IFIN
          ROWMAX = DMAX1(ROWMAX,DABS(A(JJ)))
  150   CONTINUE
        IF (DABS(A(PIVPOS))/ROWMAX.GE.RMIN) GO TO 160
        IFLAG = I
        RMIN = DABS(A(PIVPOS))/ROWMAX
  160 CONTINUE
      GO TO 180
  170 IF (LP.NE.0) WRITE (LP,FMT=99999) I
      IFLAG = -I
  180 RETURN
99999 FORMAT (' ERROR RETURN FROM MA30B/BD SINGULARITY DETECTED IN RO',
     +       'W',I8)
      END
      SUBROUTINE MA30CD(N,ICN,A,LICN,LENR,LENRL,LENOFF,IDISP,IP,IQ,X,W,
     +                  MTYPE)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER LICN,MTYPE,N
      DOUBLE PRECISION A(LICN),W(N),X(N)
      INTEGER ICN(LICN),IDISP(2),IP(N),IQ(N),LENOFF(*),LENR(N),LENRL(N)
      DOUBLE PRECISION WI,WII
      INTEGER I,IB,IBACK,IBLEND,IBLOCK,IEND,IFIRST,II,III,ILAST,J,J1,J2,
     +        J3,JJ,JPIV,JPIVP1,K,LJ1,LJ2,LT,LTEND,NUMBLK
      LOGICAL NEG,NOBLOC
      INTRINSIC DABS,DMAX1,IABS
      COMMON /MA30HD/RESID
      DOUBLE PRECISION RESID
      SAVE /MA30HD/
      RESID = ZERO
      NOBLOC = LENOFF(1) .LT. 0
      IF (MTYPE.NE.1) GO TO 140
      NEG = .FALSE.
      IP(N) = -IP(N)
      DO 10 II = 1,N
        I = IP(II)
        I = IABS(I)
        W(II) = X(I)
   10 CONTINUE
      LT = 1
      IFIRST = 1
      IBLOCK = IDISP(1)
      DO 120 I = 1,N
        WI = W(I)
        IF (NOBLOC) GO TO 30
        IF (LENOFF(I).EQ.0) GO TO 30
        LTEND = LT + LENOFF(I) - 1
        DO 20 JJ = LT,LTEND
          J = ICN(JJ)
          WI = WI - A(JJ)*W(J)
   20   CONTINUE
        LT = LTEND + 1
   30   IF (IP(I).LT.0) NEG = .TRUE.
        IF (LENRL(I).EQ.0) GO TO 50
        IEND = IBLOCK + LENRL(I) - 1
        DO 40 JJ = IBLOCK,IEND
          J = ICN(JJ)
          WI = WI + A(JJ)*W(J)
   40   CONTINUE
   50   IBLOCK = IBLOCK + LENR(I)
        W(I) = WI
        IF (.NOT.NEG) GO TO 120
        J1 = IBLOCK
        IB = I
        IF (IQ(I).GT.0) GO TO 70
        DO 60 III = IFIRST,I
          IB = I - III + IFIRST
          IF (IQ(IB).GT.0) GO TO 70
          J1 = J1 - LENR(IB)
          RESID = DMAX1(RESID,DABS(W(IB)))
          W(IB) = ZERO
   60   CONTINUE
        GO TO 110
   70   DO 100 III = IFIRST,IB
          II = IB - III + IFIRST
          J2 = J1 - 1
          J1 = J1 - LENR(II)
          JPIV = J1 + LENRL(II)
          JPIVP1 = JPIV + 1
          IF (J2.LT.JPIVP1) GO TO 90
          WII = W(II)
          DO 80 JJ = JPIVP1,J2
            J = ICN(JJ)
            WII = WII - A(JJ)*W(J)
   80     CONTINUE
          W(II) = WII
   90     W(II) = W(II)/A(JPIV)
  100   CONTINUE
  110   IFIRST = I + 1
        NEG = .FALSE.
  120 CONTINUE
      DO 130 II = 1,N
        I = IQ(II)
        I = IABS(I)
        X(I) = W(II)
  130 CONTINUE
      IP(N) = -IP(N)
      GO TO 320
  140 DO 150 II = 1,N
        I = IQ(II)
        I = IABS(I)
        W(II) = X(I)
  150 CONTINUE
      LJ1 = IDISP(1)
      IBLOCK = IDISP(2) + 1
      ILAST = N
      IBLEND = IBLOCK
      DO 290 NUMBLK = 1,N
        IF (ILAST.EQ.0) GO TO 300
        IBLOCK = IBLOCK - LENR(ILAST)
        DO 160 K = 1,N
          II = ILAST - K
          IF (II.EQ.0) GO TO 170
          IF (IP(II).LT.0) GO TO 170
          IBLOCK = IBLOCK - LENR(II)
  160   CONTINUE
  170   IFIRST = II + 1
        J1 = IBLOCK
        DO 210 I = IFIRST,ILAST
          IF (W(I).EQ.ZERO) GO TO 200
          IF (IQ(I).LT.0) GO TO 220
          J2 = J1 + LENRL(I)
          WI = W(I)/A(J2)
          IF (LENR(I)-LENRL(I).EQ.1) GO TO 190
          J2 = J2 + 1
          J3 = J1 + LENR(I) - 1
          DO 180 JJ = J2,J3
            J = ICN(JJ)
            W(J) = W(J) - A(JJ)*WI
  180     CONTINUE
  190     W(I) = WI
  200     J1 = J1 + LENR(I)
  210   CONTINUE
        GO TO 240
  220   DO 230 II = I,ILAST
          RESID = DMAX1(RESID,DABS(W(II)))
          W(II) = ZERO
  230   CONTINUE
  240   J1 = IBLEND
        DO 280 IBACK = IFIRST,ILAST
          I = ILAST - IBACK + IFIRST
          J1 = J1 - LENR(I)
          IF (LENRL(I).EQ.0) GO TO 260
          J2 = J1 + LENRL(I) - 1
          DO 250 JJ = J1,J2
            J = ICN(JJ)
            W(J) = W(J) + A(JJ)*W(I)
  250     CONTINUE
  260     IF (NOBLOC) GO TO 280
          IF (LENOFF(I).EQ.0) GO TO 280
          LJ2 = LJ1 - 1
          LJ1 = LJ1 - LENOFF(I)
          DO 270 JJ = LJ1,LJ2
            J = ICN(JJ)
            W(J) = W(J) - A(JJ)*W(I)
  270     CONTINUE
  280   CONTINUE
        IBLEND = J1
        ILAST = IFIRST - 1
  290 CONTINUE
  300 DO 310 II = 1,N
        I = IP(II)
        I = IABS(I)
        X(I) = W(II)
  310 CONTINUE
  320 RETURN
      END
      SUBROUTINE MA30DD(A,ICN,IPTR,N,IACTIV,ITOP,REALS)
      INTEGER IACTIV,ITOP,N
      LOGICAL REALS
      DOUBLE PRECISION A(ITOP)
      INTEGER ICN(ITOP),IPTR(N)
      INTEGER J,JPOS,K,KL,KN
      COMMON /MA30FD/IRNCP,ICNCP,IRANK,MINIRN,MINICN
      INTEGER ICNCP,IRANK,IRNCP,MINICN,MINIRN
      SAVE /MA30FD/
      IF (REALS) ICNCP = ICNCP + 1
      IF (.NOT.REALS) IRNCP = IRNCP + 1
      DO 10 J = 1,N
        K = IPTR(J)
        IF (K.LT.IACTIV) GO TO 10
        IPTR(J) = ICN(K)
        ICN(K) = -J
   10 CONTINUE
      KN = ITOP + 1
      KL = ITOP - IACTIV + 1
      DO 30 K = 1,KL
        JPOS = ITOP - K + 1
        IF (ICN(JPOS).EQ.0) GO TO 30
        KN = KN - 1
        IF (REALS) A(KN) = A(JPOS)
        IF (ICN(JPOS).GE.0) GO TO 20
        J = -ICN(JPOS)
        ICN(JPOS) = IPTR(J)
        IPTR(J) = KN
   20   ICN(KN) = ICN(JPOS)
   30 CONTINUE
      IACTIV = KN
      RETURN
      END
      BLOCK DATA MA30JD
      COMMON /MA30ED/LP,ABORT1,ABORT2,ABORT3
      COMMON /MA30GD/EPS,RMIN
      COMMON /MA30ID/TOL,BIG,NDROP,NSRCH,LBIG
      DOUBLE PRECISION BIG,EPS,RMIN,TOL
      INTEGER LP,NDROP,NSRCH
      LOGICAL ABORT1,ABORT2,ABORT3,LBIG
      SAVE /MA30ED/,/MA30GD/,/MA30ID/
      DATA EPS/1.0D-4/,TOL/0.0D0/,BIG/0.0D0/
      DATA LP/6/,NSRCH/32768/
      DATA LBIG/.FALSE./
      DATA ABORT1/.TRUE./,ABORT2/.TRUE./,ABORT3/.FALSE./
      END
* COPYRIGHT (c) 1993 AEA Technology
*######DATE 21 Jan 1993
C       Toolpack tool decs employed.
C	Double version of MC13D (name change only)
C
      SUBROUTINE MC13DD(N,ICN,LICN,IP,LENR,IOR,IB,NUM,IW)
C     .. Scalar Arguments ..
      INTEGER LICN,N,NUM
C     ..
C     .. Array Arguments ..
      INTEGER IB(N),ICN(LICN),IOR(N),IP(N),IW(N,3),LENR(N)
C     ..
C     .. External Subroutines ..
      EXTERNAL MC13ED
C     ..
C     .. Executable Statements ..
      CALL MC13ED(N,ICN,LICN,IP,LENR,IOR,IB,NUM,IW(1,1),IW(1,2),IW(1,3))
      RETURN

      END
      SUBROUTINE MC13ED(N,ICN,LICN,IP,LENR,ARP,IB,NUM,LOWL,NUMB,PREV)
C
C ARP(I) IS ONE LESS THAN THE NUMBER OF UNSEARCHED EDGES LEAVING
C     NODE I.  AT THE END OF THE ALGORITHM IT IS SET TO A
C     PERMUTATION WHICH PUTS THE MATRIX IN BLOCK LOWER
C     TRIANGULAR FORM.
C IB(I) IS THE POSITION IN THE ORDERING OF THE START OF THE ITH
C     BLOCK.  IB(N+1-I) HOLDS THE NODE NUMBER OF THE ITH NODE
C     ON THE STACK.
C LOWL(I) IS THE SMALLEST STACK POSITION OF ANY NODE TO WHICH A PATH
C     FROM NODE I HAS BEEN FOUND.  IT IS SET TO N+1 WHEN NODE I
C     IS REMOVED FROM THE STACK.
C NUMB(I) IS THE POSITION OF NODE I IN THE STACK IF IT IS ON
C     IT, IS THE PERMUTED ORDER OF NODE I FOR THOSE NODES
C     WHOSE FINAL POSITION HAS BEEN FOUND AND IS OTHERWISE ZERO.
C PREV(I) IS THE NODE AT THE END OF THE PATH WHEN NODE I WAS
C     PLACED ON THE STACK.
C
C
C   ICNT IS THE NUMBER OF NODES WHOSE POSITIONS IN FINAL ORDERING HAVE
C     BEEN FOUND.
C     .. Scalar Arguments ..
      INTEGER LICN,N,NUM
C     ..
C     .. Array Arguments ..
      INTEGER ARP(N),IB(N),ICN(LICN),IP(N),LENR(N),LOWL(N),NUMB(N),
     +        PREV(N)
C     ..
C     .. Local Scalars ..
      INTEGER DUMMY,I,I1,I2,ICNT,II,ISN,IST,IST1,IV,IW,J,K,LCNT,NNM1,STP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN0
C     ..
C     .. Executable Statements ..
      ICNT = 0
C NUM IS THE NUMBER OF BLOCKS THAT HAVE BEEN FOUND.
      NUM = 0
      NNM1 = N + N - 1
C
C INITIALIZATION OF ARRAYS.
      DO 20 J = 1,N
        NUMB(J) = 0
        ARP(J) = LENR(J) - 1
   20 CONTINUE
C
C
      DO 120 ISN = 1,N
C LOOK FOR A STARTING NODE
        IF (NUMB(ISN).NE.0) GO TO 120
        IV = ISN
C IST IS THE NUMBER OF NODES ON THE STACK ... IT IS THE STACK POINTER.
        IST = 1
C PUT NODE IV AT BEGINNING OF STACK.
        LOWL(IV) = 1
        NUMB(IV) = 1
        IB(N) = IV
C
C THE BODY OF THIS LOOP PUTS A NEW NODE ON THE STACK OR BACKTRACKS.
        DO 110 DUMMY = 1,NNM1
          I1 = ARP(IV)
C HAVE ALL EDGES LEAVING NODE IV BEEN SEARCHED.
          IF (I1.LT.0) GO TO 60
          I2 = IP(IV) + LENR(IV) - 1
          I1 = I2 - I1
C
C LOOK AT EDGES LEAVING NODE IV UNTIL ONE ENTERS A NEW NODE OR
C     ALL EDGES ARE EXHAUSTED.
          DO 50 II = I1,I2
            IW = ICN(II)
C HAS NODE IW BEEN ON STACK ALREADY.
            IF (NUMB(IW).EQ.0) GO TO 100
C UPDATE VALUE OF LOWL(IV) IF NECESSARY.
   50     LOWL(IV) = MIN0(LOWL(IV),LOWL(IW))
C
C THERE ARE NO MORE EDGES LEAVING NODE IV.
          ARP(IV) = -1
C IS NODE IV THE ROOT OF A BLOCK.
   60     IF (LOWL(IV).LT.NUMB(IV)) GO TO 90
C
C ORDER NODES IN A BLOCK.
          NUM = NUM + 1
          IST1 = N + 1 - IST
          LCNT = ICNT + 1
C PEEL BLOCK OFF THE TOP OF THE STACK STARTING AT THE TOP AND
C     WORKING DOWN TO THE ROOT OF THE BLOCK.
          DO 70 STP = IST1,N
            IW = IB(STP)
            LOWL(IW) = N + 1
            ICNT = ICNT + 1
            NUMB(IW) = ICNT
            IF (IW.EQ.IV) GO TO 80
   70     CONTINUE
   80     IST = N - STP
          IB(NUM) = LCNT
C ARE THERE ANY NODES LEFT ON THE STACK.
          IF (IST.NE.0) GO TO 90
C HAVE ALL THE NODES BEEN ORDERED.
          IF (ICNT.LT.N) GO TO 120
          GO TO 130
C
C BACKTRACK TO PREVIOUS NODE ON PATH.
   90     IW = IV
          IV = PREV(IV)
C UPDATE VALUE OF LOWL(IV) IF NECESSARY.
          LOWL(IV) = MIN0(LOWL(IV),LOWL(IW))
          GO TO 110
C
C PUT NEW NODE ON THE STACK.
  100     ARP(IV) = I2 - II - 1
          PREV(IW) = IV
          IV = IW
          IST = IST + 1
          LOWL(IV) = IST
          NUMB(IV) = IST
          K = N + 1 - IST
          IB(K) = IV
  110   CONTINUE
C
  120 CONTINUE
C
C
C PUT PERMUTATION IN THE REQUIRED FORM.
  130 DO 140 I = 1,N
        II = NUMB(I)
  140 ARP(II) = I
      RETURN

      END
* *******************************************************************
* COPYRIGHT (c) 1975 Hyprotech UK
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
*######DATE 4 Oct 1992
C       Toolpack tool decs employed.
C       Array lengths given explicitly eg A(MAXA)
C
      SUBROUTINE MC20AD(NC,MAXA,A,INUM,JPTR,JNUM,JDISP)
      INTEGER JDISP,MAXA,NC
      DOUBLE PRECISION A(MAXA)
      INTEGER INUM(MAXA),JNUM(MAXA),JPTR(NC)
      DOUBLE PRECISION ACE,ACEP
      INTEGER I,ICE,ICEP,J,JA,JB,JCE,JCEP,K,KR,LOC,NULL
      NULL = -JDISP
C**      CLEAR JPTR
      DO 10 J = 1,NC
        JPTR(J) = 0
   10 CONTINUE
C**      COUNT THE NUMBER OF ELEMENTS IN EACH COLUMN.
      DO 20 K = 1,MAXA
        J = JNUM(K) + JDISP
        JPTR(J) = JPTR(J) + 1
   20 CONTINUE
C**      SET THE JPTR ARRAY
      K = 1
      DO 30 J = 1,NC
        KR = K + JPTR(J)
        JPTR(J) = K
        K = KR
   30 CONTINUE
C**      REORDER THE ELEMENTS INTO COLUMN ORDER.  THE ALGORITHM IS AN
      DO 50 I = 1,MAXA
        JCE = JNUM(I) + JDISP
        IF (JCE.EQ.0) GO TO 50
        ACE = A(I)
        ICE = INUM(I)
        JNUM(I) = NULL
        DO 40 J = 1,MAXA
          LOC = JPTR(JCE)
          JPTR(JCE) = JPTR(JCE) + 1
          ACEP = A(LOC)
          ICEP = INUM(LOC)
          JCEP = JNUM(LOC)
          A(LOC) = ACE
          INUM(LOC) = ICE
          JNUM(LOC) = NULL
          IF (JCEP.EQ.NULL) GO TO 50
          ACE = ACEP
          ICE = ICEP
          JCE = JCEP + JDISP
   40   CONTINUE
   50 CONTINUE
C**      RESET JPTR VECTOR.
      JA = 1
      DO 60 J = 1,NC
        JB = JPTR(J)
        JPTR(J) = JA
        JA = JB
   60 CONTINUE
      RETURN
      END
      SUBROUTINE MC20BD(NC,MAXA,A,INUM,JPTR)
      INTEGER MAXA,NC
      DOUBLE PRECISION A(MAXA)
      INTEGER INUM(MAXA),JPTR(NC)
      DOUBLE PRECISION ACE
      INTEGER ICE,IK,J,JJ,K,KDUMMY,KLO,KMAX,KOR
      INTRINSIC IABS
      KMAX = MAXA
      DO 50 JJ = 1,NC
        J = NC + 1 - JJ
        KLO = JPTR(J) + 1
        IF (KLO.GT.KMAX) GO TO 40
        KOR = KMAX
        DO 30 KDUMMY = KLO,KMAX
          ACE = A(KOR-1)
          ICE = INUM(KOR-1)
          DO 10 K = KOR,KMAX
            IK = INUM(K)
            IF (IABS(ICE).LE.IABS(IK)) GO TO 20
            INUM(K-1) = IK
            A(K-1) = A(K)
   10     CONTINUE
          K = KMAX + 1
   20     INUM(K-1) = ICE
          A(K-1) = ACE
          KOR = KOR - 1
   30   CONTINUE
   40   KMAX = KLO - 2
   50 CONTINUE
      RETURN
      END
* COPYRIGHT (c) 1992 AEA Technology
*######DATE 8 Oct 1992
C######8/10/92 Toolpack tool decs employed.
C######8/10/92 D version created by name change only.
      SUBROUTINE MC21AD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW)
C     .. Scalar Arguments ..
      INTEGER LICN,N,NUMNZ
C     ..
C     .. Array Arguments ..
      INTEGER ICN(LICN),IP(N),IPERM(N),IW(N,4),LENR(N)
C     ..
C     .. External Subroutines ..
      EXTERNAL MC21BD
C     ..
C     .. Executable Statements ..
      CALL MC21BD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW(1,1),IW(1,2),
     +            IW(1,3),IW(1,4))
      RETURN
C
      END
      SUBROUTINE MC21BD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,PR,ARP,CV,OUT)
C   PR(I) IS THE PREVIOUS ROW TO I IN THE DEPTH FIRST SEARCH.
C IT IS USED AS A WORK ARRAY IN THE SORTING ALGORITHM.
C   ELEMENTS (IPERM(I),I) I=1, ... N  ARE NON-ZERO AT THE END OF THE
C ALGORITHM UNLESS N ASSIGNMENTS HAVE NOT BEEN MADE.  IN WHICH CASE
C (IPERM(I),I) WILL BE ZERO FOR N-NUMNZ ENTRIES.
C   CV(I) IS THE MOST RECENT ROW EXTENSION AT WHICH COLUMN I
C WAS VISITED.
C   ARP(I) IS ONE LESS THAN THE NUMBER OF NON-ZEROS IN ROW I
C WHICH HAVE NOT BEEN SCANNED WHEN LOOKING FOR A CHEAP ASSIGNMENT.
C   OUT(I) IS ONE LESS THAN THE NUMBER OF NON-ZEROS IN ROW I
C WHICH HAVE NOT BEEN SCANNED DURING ONE PASS THROUGH THE MAIN LOOP.
C
C   INITIALIZATION OF ARRAYS.
C     .. Scalar Arguments ..
      INTEGER LICN,N,NUMNZ
C     ..
C     .. Array Arguments ..
      INTEGER ARP(N),CV(N),ICN(LICN),IP(N),IPERM(N),LENR(N),OUT(N),PR(N)
C     ..
C     .. Local Scalars ..
      INTEGER I,II,IN1,IN2,IOUTK,J,J1,JORD,K,KK
C     ..
C     .. Executable Statements ..
      DO 10 I = 1,N
        ARP(I) = LENR(I) - 1
        CV(I) = 0
        IPERM(I) = 0
   10 CONTINUE
      NUMNZ = 0
C
C
C   MAIN LOOP.
C   EACH PASS ROUND THIS LOOP EITHER RESULTS IN A NEW ASSIGNMENT
C OR GIVES A ROW WITH NO ASSIGNMENT.
      DO 100 JORD = 1,N
        J = JORD
        PR(J) = -1
        DO 70 K = 1,JORD
C LOOK FOR A CHEAP ASSIGNMENT
          IN1 = ARP(J)
          IF (IN1.LT.0) GO TO 30
          IN2 = IP(J) + LENR(J) - 1
          IN1 = IN2 - IN1
          DO 20 II = IN1,IN2
            I = ICN(II)
            IF (IPERM(I).EQ.0) GO TO 80
   20     CONTINUE
C   NO CHEAP ASSIGNMENT IN ROW.
          ARP(J) = -1
C   BEGIN LOOKING FOR ASSIGNMENT CHAIN STARTING WITH ROW J.
   30     CONTINUE
          OUT(J) = LENR(J) - 1
C INNER LOOP.  EXTENDS CHAIN BY ONE OR BACKTRACKS.
          DO 60 KK = 1,JORD
            IN1 = OUT(J)
            IF (IN1.LT.0) GO TO 50
            IN2 = IP(J) + LENR(J) - 1
            IN1 = IN2 - IN1
C FORWARD SCAN.
            DO 40 II = IN1,IN2
              I = ICN(II)
              IF (CV(I).EQ.JORD) GO TO 40
C   COLUMN I HAS NOT YET BEEN ACCESSED DURING THIS PASS.
              J1 = J
              J = IPERM(I)
              CV(I) = JORD
              PR(J) = J1
              OUT(J1) = IN2 - II - 1
              GO TO 70
C
   40       CONTINUE
C
C   BACKTRACKING STEP.
   50       CONTINUE
            J = PR(J)
            IF (J.EQ.-1) GO TO 100
   60     CONTINUE
C
   70   CONTINUE
C
C   NEW ASSIGNMENT IS MADE.
   80   CONTINUE
        IPERM(I) = J
        ARP(J) = IN2 - II - 1
        NUMNZ = NUMNZ + 1
        DO 90 K = 1,JORD
          J = PR(J)
          IF (J.EQ.-1) GO TO 100
          II = IP(J) + LENR(J) - OUT(J) - 2
          I = ICN(II)
          IPERM(I) = J
   90   CONTINUE
C
  100 CONTINUE
C
C   IF MATRIX IS STRUCTURALLY SINGULAR, WE NOW COMPLETE THE
C PERMUTATION IPERM.
      IF (NUMNZ.EQ.N) RETURN
      DO 110 I = 1,N
        ARP(I) = 0
  110 CONTINUE
      K = 0
      DO 130 I = 1,N
        IF (IPERM(I).NE.0) GO TO 120
        K = K + 1
        OUT(K) = I
        GO TO 130
C
  120   CONTINUE
        J = IPERM(I)
        ARP(J) = I
  130 CONTINUE
      K = 0
      DO 140 I = 1,N
        IF (ARP(I).NE.0) GO TO 140
        K = K + 1
        IOUTK = OUT(K)
        IPERM(IOUTK) = I
  140 CONTINUE
      RETURN
C
      END
* COPYRIGHT (c) 1993 AEA Technology
*######DATE 21 Jan 1993
C       Toolpack tool decs employed.
C
      SUBROUTINE MC22AD(N,ICN,A,NZ,LENROW,IP,IQ,IW,IW1)
C     .. Scalar Arguments ..
      INTEGER N,NZ
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(NZ)
      INTEGER ICN(NZ),IP(N),IQ(N),IW(N,2),IW1(NZ),LENROW(N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AVAL
      INTEGER I,ICHAIN,IOLD,IPOS,J,J2,JJ,JNUM,JVAL,LENGTH,NEWPOS
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC IABS
C     ..
C     .. Executable Statements ..
      IF (NZ.LE.0) GO TO 1000
      IF (N.LE.0) GO TO 1000
C SET START OF ROW I IN IW(I,1) AND LENROW(I) IN IW(I,2)
      IW(1,1) = 1
      IW(1,2) = LENROW(1)
      DO 10 I = 2,N
        IW(I,1) = IW(I-1,1) + LENROW(I-1)
   10 IW(I,2) = LENROW(I)
C PERMUTE LENROW ACCORDING TO IP.  SET OFF-SETS FOR NEW POSITION
C     OF ROW IOLD IN IW(IOLD,1) AND PUT OLD ROW INDICES IN IW1 IN
C     POSITIONS CORRESPONDING TO THE NEW POSITION OF THIS ROW IN A/ICN.
      JJ = 1
      DO 20 I = 1,N
        IOLD = IP(I)
        IOLD = IABS(IOLD)
        LENGTH = IW(IOLD,2)
        LENROW(I) = LENGTH
        IF (LENGTH.EQ.0) GO TO 20
        IW(IOLD,1) = IW(IOLD,1) - JJ
        J2 = JJ + LENGTH - 1
        DO 15 J = JJ,J2
   15   IW1(J) = IOLD
        JJ = J2 + 1
   20 CONTINUE
C SET INVERSE PERMUTATION TO IQ IN IW(.,2).
      DO 30 I = 1,N
        IOLD = IQ(I)
        IOLD = IABS(IOLD)
   30 IW(IOLD,2) = I
C PERMUTE A AND ICN IN PLACE, CHANGING TO NEW COLUMN NUMBERS.
C
C ***   MAIN LOOP   ***
C EACH PASS THROUGH THIS LOOP PLACES A CLOSED CHAIN OF COLUMN INDICES
C     IN THEIR NEW (AND FINAL) POSITIONS ... THIS IS RECORDED BY
C     SETTING THE IW1 ENTRY TO ZERO SO THAT ANY WHICH ARE SUBSEQUENTLY
C     ENCOUNTERED DURING THIS MAJOR SCAN CAN BE BYPASSED.
      DO 200 I = 1,NZ
        IOLD = IW1(I)
        IF (IOLD.EQ.0) GO TO 200
        IPOS = I
        JVAL = ICN(I)
C IF ROW IOLD IS IN SAME POSITIONS AFTER PERMUTATION GO TO 150.
        IF (IW(IOLD,1).EQ.0) GO TO 150
        AVAL = A(I)
C **  CHAIN LOOP  **
C EACH PASS THROUGH THIS LOOP PLACES ONE (PERMUTED) COLUMN INDEX
C     IN ITS FINAL POSITION  .. VIZ. IPOS.
        DO 100 ICHAIN = 1,NZ
C NEWPOS IS THE ORIGINAL POSITION IN A/ICN OF THE ELEMENT TO BE PLACED
C IN POSITION IPOS.  IT IS ALSO THE POSITION OF THE NEXT ELEMENT IN
C     THE CHAIN.
          NEWPOS = IPOS + IW(IOLD,1)
C IS CHAIN COMPLETE ?
          IF (NEWPOS.EQ.I) GO TO 130
          A(IPOS) = A(NEWPOS)
          JNUM = ICN(NEWPOS)
          ICN(IPOS) = IW(JNUM,2)
          IPOS = NEWPOS
          IOLD = IW1(IPOS)
          IW1(IPOS) = 0
C **  END OF CHAIN LOOP  **
  100   CONTINUE
  130   A(IPOS) = AVAL
  150   ICN(IPOS) = IW(JVAL,2)
C ***   END OF MAIN LOOP   ***
  200 CONTINUE
C
 1000 RETURN

      END
* COPYRIGHT (c) 1993 AEA Technology
*######DATE 21 Jan 1993
C       Toolpack tool decs employed.
C       SAVE statements added.
C       MC23CD reference removed.
C 12/12/94 Calls of MC13D and MC21A changed to MC13DD and MC21AD
C
C  EAT 21/6/93 EXTERNAL statement put in for block data on VAXs.
C
C
      SUBROUTINE MC23AD(N,ICN,A,LICN,LENR,IDISP,IP,IQ,LENOFF,IW,IW1)
C INPUT ... N,ICN .. A,ICN,LENR ....
C
C SET UP POINTERS IW(.,1) TO THE BEGINNING OF THE ROWS AND SET LENOFF
C     EQUAL TO LENR.
C     .. Scalar Arguments ..
      INTEGER LICN,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LICN)
      INTEGER ICN(LICN),IDISP(2),IP(N),IQ(N),IW(N,5),IW1(N,2),LENOFF(N),
     +        LENR(N)
C     ..
C     .. Local Scalars ..
      INTEGER I,I1,I2,IBEG,IBLOCK,IEND,II,ILEND,INEW,IOLD,IROWB,IROWE,J,
     +        JJ,JNEW,JNPOS,JOLD,K,LENI,NZ
C     ..
C     .. External Subroutines ..
      EXTERNAL MC13DD,MC21AD
C     ..
C     .. Data block external statement
      EXTERNAL MC23CD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX0,MIN0
C     ..
C     .. Common blocks ..
      COMMON /MC23BD/LP,NUMNZ,NUM,LARGE,ABORT
      INTEGER LARGE,LP,NUM,NUMNZ
      LOGICAL ABORT
C     ..
C     .. Save statement ..
      SAVE /MC23BD/
C     ..
C     .. Executable Statements ..
      IW1(1,1) = 1
      LENOFF(1) = LENR(1)
      IF (N.EQ.1) GO TO 20
      DO 10 I = 2,N
        LENOFF(I) = LENR(I)
   10 IW1(I,1) = IW1(I-1,1) + LENR(I-1)
C IDISP(1) POINTS TO THE FIRST POSITION IN A/ICN AFTER THE
C     OFF-DIAGONAL BLOCKS AND UNTREATED ROWS.
   20 IDISP(1) = IW1(N,1) + LENR(N)
C
C FIND ROW PERMUTATION IP TO MAKE DIAGONAL ZERO-FREE.
      CALL MC21AD(N,ICN,LICN,IW1,LENR,IP,NUMNZ,IW)
C
C POSSIBLE ERROR RETURN FOR STRUCTURALLY SINGULAR MATRICES.
      IF (NUMNZ.NE.N .AND. ABORT) GO TO 170
C
C IW1(.,2) AND LENR ARE PERMUTATIONS OF IW1(.,1) AND LENR/LENOFF
C     SUITABLE FOR ENTRY
C     TO MC13DD SINCE MATRIX WITH THESE ROW POINTER AND LENGTH ARRAYS
C     HAS MAXIMUM NUMBER OF NON-ZEROS ON THE DIAGONAL.
      DO 30 II = 1,N
        I = IP(II)
        IW1(II,2) = IW1(I,1)
   30 LENR(II) = LENOFF(I)
C
C FIND SYMMETRIC PERMUTATION IQ TO BLOCK LOWER TRIANGULAR FORM.
      CALL MC13DD(N,ICN,LICN,IW1(1,2),LENR,IQ,IW(1,4),NUM,IW)
C
      IF (NUM.NE.1) GO TO 60
C
C ACTION TAKEN IF MATRIX IS IRREDUCIBLE.
C WHOLE MATRIX IS JUST MOVED TO THE END OF THE STORAGE.
      DO 40 I = 1,N
        LENR(I) = LENOFF(I)
        IP(I) = I
   40 IQ(I) = I
      LENOFF(1) = -1
C IDISP(1) IS THE FIRST POSITION AFTER THE LAST ELEMENT IN THE
C     OFF-DIAGONAL BLOCKS AND UNTREATED ROWS.
      NZ = IDISP(1) - 1
      IDISP(1) = 1
C IDISP(2) IS THE POSITION IN A/ICN OF THE FIRST ELEMENT IN THE
C     DIAGONAL BLOCKS.
      IDISP(2) = LICN - NZ + 1
      LARGE = N
      IF (NZ.EQ.LICN) GO TO 230
      DO 50 K = 1,NZ
        J = NZ - K + 1
        JJ = LICN - K + 1
        A(JJ) = A(J)
   50 ICN(JJ) = ICN(J)
C 230 = RETURN
      GO TO 230
C
C DATA STRUCTURE REORDERED.
C
C FORM COMPOSITE ROW PERMUTATION ... IP(I) = IP(IQ(I)).
   60 DO 70 II = 1,N
        I = IQ(II)
   70 IW(II,1) = IP(I)
      DO 80 I = 1,N
   80 IP(I) = IW(I,1)
C
C RUN THROUGH BLOCKS IN REVERSE ORDER SEPARATING DIAGONAL BLOCKS
C     WHICH ARE MOVED TO THE END OF THE STORAGE.  ELEMENTS IN
C     OFF-DIAGONAL BLOCKS ARE LEFT IN PLACE UNLESS A COMPRESS IS
C     NECESSARY.
C
C IBEG INDICATES THE LOWEST VALUE OF J FOR WHICH ICN(J) HAS BEEN
C     SET TO ZERO WHEN ELEMENT IN POSITION J WAS MOVED TO THE
C     DIAGONAL BLOCK PART OF STORAGE.
      IBEG = LICN + 1
C IEND IS THE POSITION OF THE FIRST ELEMENT OF THOSE TREATED ROWS
C     WHICH ARE IN DIAGONAL BLOCKS.
      IEND = LICN + 1
C LARGE IS THE DIMENSION OF THE LARGEST BLOCK ENCOUNTERED SO FAR.
      LARGE = 0
C
C NUM IS THE NUMBER OF DIAGONAL BLOCKS.
      DO 150 K = 1,NUM
        IBLOCK = NUM - K + 1
C I1 IS FIRST ROW (IN PERMUTED FORM) OF BLOCK IBLOCK.
C I2 IS LAST ROW (IN PERMUTED FORM) OF BLOCK IBLOCK.
        I1 = IW(IBLOCK,4)
        I2 = N
        IF (K.NE.1) I2 = IW(IBLOCK+1,4) - 1
        LARGE = MAX0(LARGE,I2-I1+1)
C GO THROUGH THE ROWS OF BLOCK IBLOCK IN THE REVERSE ORDER.
        DO 140 II = I1,I2
          INEW = I2 - II + I1
C WE NOW DEAL WITH ROW INEW IN PERMUTED FORM (ROW IOLD IN ORIGINAL
C     MATRIX).
          IOLD = IP(INEW)
C IF THERE IS SPACE TO MOVE UP DIAGONAL BLOCK PORTION OF ROW GO TO 110
          IF (IEND-IDISP(1).GE.LENOFF(IOLD)) GO TO 110
C
C IN-LINE COMPRESS.
C MOVES SEPARATED OFF-DIAGONAL ELEMENTS AND UNTREATED ROWS TO
C     FRONT OF STORAGE.
          JNPOS = IBEG
          ILEND = IDISP(1) - 1
          IF (ILEND.LT.IBEG) GO TO 190
          DO 90 J = IBEG,ILEND
            IF (ICN(J).EQ.0) GO TO 90
            ICN(JNPOS) = ICN(J)
            A(JNPOS) = A(J)
            JNPOS = JNPOS + 1
   90     CONTINUE
          IDISP(1) = JNPOS
          IF (IEND-JNPOS.LT.LENOFF(IOLD)) GO TO 190
          IBEG = LICN + 1
C RESET POINTERS TO THE BEGINNING OF THE ROWS.
          DO 100 I = 2,N
  100     IW1(I,1) = IW1(I-1,1) + LENOFF(I-1)
C
C ROW IOLD IS NOW SPLIT INTO DIAG. AND OFF-DIAG. PARTS.
  110     IROWB = IW1(IOLD,1)
          LENI = 0
          IROWE = IROWB + LENOFF(IOLD) - 1
C BACKWARD SCAN OF WHOLE OF ROW IOLD (IN ORIGINAL MATRIX).
          IF (IROWE.LT.IROWB) GO TO 130
          DO 120 JJ = IROWB,IROWE
            J = IROWE - JJ + IROWB
            JOLD = ICN(J)
C IW(.,2) HOLDS THE INVERSE PERMUTATION TO IQ.
C     ..... IT WAS SET TO THIS IN MC13DD.
            JNEW = IW(JOLD,2)
C IF (JNEW.LT.I1) THEN ....
C ELEMENT IS IN OFF-DIAGONAL BLOCK AND SO IS LEFT IN SITU.
            IF (JNEW.LT.I1) GO TO 120
C ELEMENT IS IN DIAGONAL BLOCK AND IS MOVED TO THE END OF THE STORAGE.
            IEND = IEND - 1
            A(IEND) = A(J)
            ICN(IEND) = JNEW
            IBEG = MIN0(IBEG,J)
            ICN(J) = 0
            LENI = LENI + 1
  120     CONTINUE
C
          LENOFF(IOLD) = LENOFF(IOLD) - LENI
  130     LENR(INEW) = LENI
  140   CONTINUE
C
        IP(I2) = -IP(I2)
  150 CONTINUE
C RESETS IP(N) TO POSITIVE VALUE.
      IP(N) = -IP(N)
C IDISP(2) IS POSITION OF FIRST ELEMENT IN DIAGONAL BLOCKS.
      IDISP(2) = IEND
C
C THIS COMPRESS IS USED TO MOVE ALL OFF-DIAGONAL ELEMENTS TO THE
C     FRONT OF THE STORAGE.
      IF (IBEG.GT.LICN) GO TO 230
      JNPOS = IBEG
      ILEND = IDISP(1) - 1
      DO 160 J = IBEG,ILEND
        IF (ICN(J).EQ.0) GO TO 160
        ICN(JNPOS) = ICN(J)
        A(JNPOS) = A(J)
        JNPOS = JNPOS + 1
  160 CONTINUE
C IDISP(1) IS FIRST POSITION AFTER LAST ELEMENT OF OFF-DIAGONAL BLOCKS.
      IDISP(1) = JNPOS
      GO TO 230
C
C
C ERROR RETURN
  170 IF (LP.NE.0) WRITE (LP,FMT=180) NUMNZ

  180 FORMAT (/,' ERROR RETURN FROM MC23A  BECAUSE',/,10X,
     +       ' MATRIX IS STRUCTURALLY SINGULAR, RANK = ',I6)

      IDISP(1) = -1
      GO TO 230

  190 IF (LP.NE.0) WRITE (LP,FMT=200) N

  200 FORMAT (/,' ERROR RETURN FROM MC23A  BECAUSE',/,10X,
     +       ' LICN NOT BIG ENOUGH INCREASE BY ',I6)

      IDISP(1) = -2
C
  230 RETURN

      END
      BLOCK DATA MC23CD
C     .. Common blocks ..
      COMMON /MC23BD/LP,NUMNZ,NUM,LARGE,ABORT
      INTEGER LARGE,LP,NUM,NUMNZ
      LOGICAL ABORT
C     ..
C     .. Save statement ..
      SAVE /MC23BD/
C     ..
C     .. Data statements ..
      DATA LP/6/,ABORT/.FALSE./
C     ..
C     .. Executable Statements ..
      END
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
C       ZERO made PARAMETER.
C
      SUBROUTINE MC24AD(N,ICN,A,LICN,LENR,LENRL,W)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER LICN,N
      DOUBLE PRECISION A(LICN),W(N)
      INTEGER ICN(LICN),LENR(N),LENRL(N)
      DOUBLE PRECISION AMAXL,AMAXU,WROWL
      INTEGER I,J,J0,J1,J2,JJ
      INTRINSIC DABS,DMAX1
      AMAXL = ZERO
      DO 10 I = 1,N
   10 W(I) = ZERO
      J0 = 1
      DO 100 I = 1,N
        IF (LENR(I).EQ.0) GO TO 100
        J2 = J0 + LENR(I) - 1
        IF (LENRL(I).EQ.0) GO TO 50
        J1 = J0 + LENRL(I) - 1
        WROWL = ZERO
        DO 30 JJ = J0,J1
   30   WROWL = WROWL + DABS(A(JJ))
        AMAXL = DMAX1(AMAXL,WROWL)
        J0 = J1 + 1
   50   J0 = J0 + 1
        IF (J0.GT.J2) GO TO 90
        DO 80 JJ = J0,J2
          J = ICN(JJ)
   80   W(J) = DMAX1(DABS(A(JJ)),W(J))
   90   J0 = J2 + 1
  100 CONTINUE
      AMAXU = ZERO
      DO 200 I = 1,N
  200 AMAXU = DMAX1(AMAXU,W(I))
      W(1) = AMAXL*AMAXU
      RETURN
      END

       DOUBLE PRECISION FUNCTION FD15AD(T)
C----------------------------------------------------------------
C  Fortran 77 implementation of the Fortran 90 intrinsic
C    functions: EPSILON, TINY, HUGE and RADIX.  Note that
C    the RADIX result is returned as DOUBLE PRECISION.
C
C  The CHARACTER argument specifies the type of result:
C       
C   'E'  smallest positive real number: 1.0 + DC(1) > 1.0, i.e.
C          EPSILON(DOUBLE PRECISION)
C   'T'  smallest full precision positive real number, i.e.
C          TINY(DOUBLE PRECISION)
C   'H'  largest finite positive real number, i.e.
C          HUGE(DOUBLE PRECISION)
C   'R'  the base of the floating point arithematic, i.e.
C          RADIX(DOUBLE PRECISION)
C
C    any other value gives a result of zero.
C----------------------------------------------------------------
      CHARACTER T

      IF ( T.EQ.'E' ) THEN
         FD15AD = EPSILON(1.0D0)
      ELSE IF ( T.EQ.'T' ) THEN
         FD15AD = TINY(1.0D0)
      ELSE IF ( T.EQ.'H' ) THEN
         FD15AD = HUGE(1.0D0)
      ELSE IF ( T.EQ.'R' ) THEN
         FD15AD = DBLE(RADIX(1.0D0))
      ELSE
         FD15AD = 0.0D0
      ENDIF
      RETURN
      END
