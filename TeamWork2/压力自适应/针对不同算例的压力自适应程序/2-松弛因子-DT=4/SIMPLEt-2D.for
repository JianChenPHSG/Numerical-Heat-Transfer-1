Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  This computer program was copied from the graduate student course program
C  of the University of Minnesota of USA, and that of Xi'an Jiaotong
C  University of China. Part of it was re-formulated to meet the personal
C  computer environment, and for the transient simulations. The program is
C  used only for teaching purpose. No part of it may be published. You may
C  use it as a frame to re-develop your own code for your research purpose.
C  Tsinghua University, Beijing 100084, China, March, 2013
****************************************************************************
*----------------------------MAIN PROGRAM-----------------------------------
****************************************************************************
      LOGICAL LSOLVE,LPRINT,LBLK,LSTOP,LEND
      COMMON/CNTL/LSTOP,LEND
      COMMON /INDX/ITIME,MTIME,ITER
      COMMON/modl/REV_P,C_TIME
      REAL::T0,FINISH
****************************************************************************
      OPEN(08,FILE='OUTPUT.DAT')
	OPEN(UNIT=12,FILE='recyle.txt')
	REV_P=0.0
	DO i=1,11,1
*---------------------TIME----------------------------------------------------
      CALL CPU_TIME(T0)
*-------------------------------------------------------------------------------
	  ! use the cycle to obtain the constant relax data
**************************   revise  **********************************************
      CALL SETUP0
      CALL GRID
*-------------------------------------------------------------------------------
      C_TIME=0
	REV_P=REV_P+0.1
*----------------------------------------------------------------------------
      CALL START
      CALL SETUP1
      DO
	  DO
	   CALL DENSE
         CALL BOUND
	   CALL SETUP2
	   CALL OUTPUT
	   CALL ISSTOP
	   IF(LSTOP) EXIT
      END DO
	  CALL OUTPUT
        CALL TSTEP
	  CALL STEPUP
	  IF(LEND) THEN
*----------------------------------------------------------------------------------
         CALL CPU_TIME(FINISH)
         PRINT '("TIME=",f20.10,"seconds.")',FINISH-T0
		 WRITE(12,1) REV_P,C_TIME,FINISH-T0
1        FORMAT(3F20.10)  	 
*----------------------------------------------------------------------------------
	   EXIT
         ENDIF  
      END DO
      END DO
	CLOSE(08)
      END
*---------------------------------------------------------------------------
      SUBROUTINE DIFLOW
****************************************************************************
      COMMON/COEF/FLOW,DIFF,ACOF
****************************************************************************
      ACOF=DIFF
      IF(FLOW .EQ.0.0)RETURN
      TEMP=DIFF-ABS(FLOW)*0.1
      ACOF=0.
      IF(TEMP .LE. 0. ) RETURN
      TEMP=TEMP/DIFF
      ACOF=DIFF*TEMP**5
      RETURN
      END
*--------------------------------------------------------------------------
      SUBROUTINE SOLVE
****************************************************************************
	PARAMETER (NX=22,NY=22,NXY=22)
      CHARACTER (LEN=20) TITLE
      LOGICAL LSOLVE,LPRINT,LBLK,LSTOP,LEND
      COMMON F(NX,NY,13),FOLD(NX,NY,13),P(NX,NY),RHO(NX,NY),GAM(NX,NY),
     &  CON(NX,NY),AIP(NX,NY),AIM(NX,NY),AJP(NX,NY),AJM(NX,NY),
     &  AP(NX,NY),X(NXY),XU(NXY),XDIF(NXY),XCV(NXY),XCVS(NXY),
     &  Y(NXY),YV(NXY),YDIF(NXY),YCV(NXY),YCVS(NXY),
     &  YCVR(NXY),YCVRS(NXY),ARX(NXY),ARXJ(NXY),ARXJP(NXY),
     &  R(NXY),RMN(NXY),SX(NXY),SXMN(NXY),XCVI(NXY),XCVIP(NXY)
      COMMON DU(NX,NY),DV(NX,NY),FV(NXY),FVP(NXY),
     &  FX(NXY),FXM(NXY),FY(NXY),FYM(NXY),PT(NXY),QT(NXY)
      COMMON /INDX/NF,NFMAX,NP,NRHO,NGAM,L1,L2,L3,M1,M2,M3,
     &  IST,JST,ITER,LAST,TITLE(13),RELAX(13),TIME,DT,XL,YL,
     &  IPREF,JPREF,LSOLVE(10),LPRINT(13),LBLK(10),MODE,NTIMES(10),
     &  RHOCON,PRELAX(NX,NY)
      COMMON/CNTL/LSTOP,LEND
      COMMON ITIME,MTIME
      COMMON/modl/REV_P,C_TIME
****************************************************************************
      ISTF=IST-1
      JSTF=JST-1
      IT1=L2+IST
      IT2=L3+IST
      JT1=M2+JST
      JT2=M3+JST
****************************************************************************
      DO 999 NT=1,NTIMES(NF)
      DO 999 N=NF,NF
*---------------------------------------------------------------------------
      IF(.NOT. LBLK(NF)) GO TO 10
      PT(ISTF)=0.
      QT(ISTF)=0.
      DO 11 I=IST,L2
      BL=0.
      BLP=0.
      BLM=0.
      BLC=0.
      DO 12 J=JST,M2
      BL=BL+AP(I,J)
      IF(J .NE. M2) BL=BL-AJP(I,J)
      IF(J .NE. JST) BL=BL-AJM(I,J)
      BLP=BLP+AIP(I,J)
      BLM=BLM+AIM(I,J)
      BLC=BLC+CON(I,J)+AIP(I,J)*F(I+1,J,N)+AIM(I,J)*F(I-1,J,N)
     &   +AJP(I,J)*F(I,J+1,N)+AJM(I,J)*F(I,J-1,N)-AP(I,J)*F(I,J,N)
   12 CONTINUE
      DENOM=BL-PT(I-1)*BLM
      DENO=1.E15
      IF(ABS(DENOM/BL) .LT. 1.E-10) DENOM=1.E20*DENO
      PT(I)=BLP/DENOM
      QT(I)=(BLC+BLM*QT(I-1))/DENOM
   11 CONTINUE
      BL=0.
      DO 13 II=IST,L2
      I=IT1-II
      BL=BL*PT(I)+QT(I)
      DO 13 J=JST,M2
   13 F(I,J,N)=F(I,J,N)+BL
*---------------------------------------------------------------------------
      PT(JSTF)=0.
      QT(JSTF)=0.
      DO 21 J=JST,M2
      BL=0.
      BLP=0.
      BLM=0.
      BLC=0.
      DO 22 I=IST,L2
      BL=BL+AP(I,J)
      IF(I .NE. L2) BL=BL-AIP(I,J)
      IF(I .NE. IST) BL=BL-AIM(I,J)
      BLP=BLP+AJP(I,J)
      BLM=BLM+AJM(I,J)
      BLC=BLC+CON(I,J)+AIP(I,J)*F(I+1,J,N)+AIM(I,J)*F(I-1,J,N)
     &   +AJP(I,J)*F(I,J+1,N)+AJM(I,J)*F(I,J-1,N)-AP(I,J)*F(I,J,N)
   22 CONTINUE
      DENOM=BL-PT(J-1)*BLM
      IF (ABS(DENOM/BL) .LT. 1E-10) DENOM=1.E20*DENO
      PT(J)=BLP/DENOM
      QT(J)=(BLC+BLM*QT(J-1))/DENOM
   21 CONTINUE
      BL=0.
      DO 23 JJ=JST,M2
      J=JT1-JJ
      BL=BL*PT(J)+QT(J)
      DO 23 I=IST,L2
   23 F(I,J,N)=F(I,J,N)+BL
   10 CONTINUE
*-----------------------------------------------------------------------
      DO 90 J=JST,M2
      PT(ISTF)=0.
      QT(ISTF)=F(ISTF,J,N)
      DO 70 I=IST,L2
      DENOM=AP(I,J)-PT(I-1)*AIM(I,J)
      PT(I)=AIP(I,J)/DENOM
      TEMP=CON(I,J)+AJP(I,J)*F(I,J+1,N)+AJM(I,J)*F(I,J-1,N)
      QT(I)=(TEMP+AIM(I,J)*QT(I-1))/DENOM
   70 CONTINUE
      DO 80 II=IST,L2
      I=IT1-II
   80 F(I,J,N)=F(I+1,J,N)*PT(I)+QT(I)
   90 CONTINUE
*-----------------------------------------------------------------------
      DO 190 JJ=JST,M3
      J=JT2-JJ
      PT(ISTF)=0.
      QT(ISTF)=F(ISTF,J,N)
      DO 170 I=IST,L2
      DENOM=AP(I,J)-PT(I-1)*AIM(I,J)
      PT(I)=AIP(I,J)/DENOM
      TEMP=CON(I,J)+AJP(I,J)*F(I,J+1,N)+AJM(I,J)*F(I,J-1,N)
      QT(I)=(TEMP+AIM(I,J)*QT(I-1))/DENOM
  170 CONTINUE
      DO 180 II=IST,L2
      I=IT1-II
  180 F(I,J,N)=F(I+1,J,N)*PT(I)+QT(I)
  190 CONTINUE
*-----------------------------------------------------------------------
      DO 290 I=IST,L2
      PT(JSTF)=0.
      QT(JSTF)=F(I,JSTF,N)
      DO 270 J=JST,M2
      DENOM=AP(I,J)-PT(J-1)*AJM(I,J)
      PT(J)=AJP(I,J)/DENOM
      TEMP=CON(I,J)+AIP(I,J)*F(I+1,J,N)+AIM(I,J)*F(I-1,J,N)
      QT(J)=(TEMP+AJM(I,J)*QT(J-1))/DENOM
  270 CONTINUE
      DO 280 JJ=JST,M2
      J=JT1-JJ
  280 F(I,J,N)=F(I,J+1,N)*PT(J)+QT(J)
  290 CONTINUE
*-----------------------------------------------------------------------
      DO 390 II=IST,L3
      I=IT2-II
      PT(JSTF)=0.
      QT(JSTF)=F(I,JSTF,N)
      DO 370 J=JST,M2
      DENOM=AP(I,J)-PT(J-1)*AJM(I,J)
      PT(J)=AJP(I,J)/DENOM
      TEMP=CON(I,J)+AIP(I,J)*F(I+1,J,N)+AIM(I,J)*F(I-1,J,N)
      QT(J)=(TEMP+AJM(I,J)*QT(J-1))/DENOM
  370 CONTINUE
      DO 380 JJ=JST,M2
      J=JT1-JJ
  380 F(I,J,N)=F(I,J+1,N)*PT(J)+QT(J)
  390 CONTINUE
************************************************************************
  999 CONTINUE
      DO 400 J=2,M2
      DO 400 I=2,L2
      CON(I,J)=0.
      AP(I,J)=0.
  400 CONTINUE
      RETURN
      END
************************************************************************
      SUBROUTINE SETUP
************************************************************************
	PARAMETER (NX=22,NY=22,NXY=22)
      CHARACTER (LEN=20) TITLE
      LOGICAL LSOLVE,LPRINT,LBLK,LSTOP,LEND
      COMMON F(NX,NY,13),FOLD(NX,NY,13),P(NX,NY),RHO(NX,NY),GAM(NX,NY),
     &  CON(NX,NY),AIP(NX,NY),AIM(NX,NY),AJP(NX,NY),AJM(NX,NY),
     &  AP(NX,NY),X(NXY),XU(NXY),XDIF(NXY),XCV(NXY),XCVS(NXY),
     &  Y(NXY),YV(NXY),YDIF(NXY),YCV(NXY),YCVS(NXY),
     &  YCVR(NXY),YCVRS(NXY),ARX(NXY),ARXJ(NXY),ARXJP(NXY),
     &  R(NXY),RMN(NXY),SX(NXY),SXMN(NXY),XCVI(NXY),XCVIP(NXY)
      COMMON DU(NX,NY),DV(NX,NY),FV(NXY),FVP(NXY),
     &  FX(NXY),FXM(NXY),FY(NXY),FYM(NXY),PT(NXY),QT(NXY)
      COMMON/INDX/NF,NFMAX,NP,NRHO,NGAM,L1,L2,L3,M1,M2,M3,
     &  IST,JST,ITER,LAST,TITLE(13),RELAX(13),TIME,DT,XL,YL,
     &  IPREF,JPREF,LSOLVE(10),LPRINT(13),LBLK(10),MODE,NTIMES(10),
     &  RHOCON,PRELAX(NX,NY)
      COMMON/CNTL/LSTOP,LEND
      COMMON ITIME,MTIME   
      COMMON/modl/REV_P,C_TIME
      COMMON/SORC/SMAX,SSUM,RSMAX
      COMMON/COEF/FLOW,DIFF,ACOF
      DIMENSION U(NX,NY),V(NX,NY),PC(NX,NY)
      EQUIVALENCE (F(1,1,1),U(1,1)),(F(1,1,2),V(1,1)),(F(1,1,3),PC(1,1))
************************************************************************
    1 FORMAT(//15X,'COMPUTATION IN CARTISIAN COORDINATES')
    2 FORMAT(//15X,'COMPUTATION FOR AXISYMMETRICAL SITUATION')
    3 FORMAT(//15X,' COMPUTATION IN POLAR COORDINATES  ')
    4 FORMAT(1X,14X,40(1H*),//)
*-----------------------------------------------------------------------
      ENTRY SETUP0
      NFMAX=10
      NP=NFMAX+1
      NRHO=NFMAX+2
      NGAM=NFMAX+3
      LSTOP=.FALSE.
      DO 779 I=1,10
      LSOLVE(I)=.FALSE.
      LBLK(I)=.TRUE.
 779  NTIMES(I)=1
      DO 889 I=1,13
      LPRINT(I)=.FALSE.
 889  RELAX(I)=1.
      MODE=1
      LAST=5
      TIME=0.
      ITER=0
      DT=1.0E+10
      IPREF=1
      JPREF=1
      RHOCON=1
      RETURN
*-----------------------------------------------------------------------
      ENTRY SETUP1
      L2=L1-1
      L3=L2-1
      M2=M1-1
      M3=M2-1
      X(1)=XU(2)
      DO 5 I=2,L2
    5 X(I)=0.5*(XU(I+1)+XU(I))
      X(L1)=XU(L1)
      Y(1)=YV(2)
      DO 10 J=2,M2
   10 Y(J)=0.5*(YV(J+1)+YV(J))
      Y(M1)=YV(M1)
      DO 15 I=2,L1
   15 XDIF(I)=X(I)-X(I-1)
      DO 18 I=2,L2
   18 XCV(I)=XU(I+1)-XU(I)
      DO 20 I=3,L2
   20 XCVS(I)=XDIF(I)
      XCVS(3)=XCVS(3)+XDIF(2)
      XCVS(L2)=XCVS(L2)+XDIF(L1)
      DO 22 I=3,L3
      XCVI(I)=0.5*XCV(I)
   22 XCVIP(I)=XCVI(I)
      XCVIP(2)=XCV(2)
      XCVI(L2)=XCV(L2)
      DO 35 J=2,M1
   35 YDIF(J)=Y(J)-Y(J-1)
      DO 40 J=2,M2
   40 YCV(J)=YV(J+1)-YV(J)
      DO 45 J=3,M2
   45 YCVS(J)=YDIF(J)
      YCVS(3)=YCVS(3)+YDIF(2)
      YCVS(M2)=YCVS(M2)+YDIF(M1)
      IF (MODE .NE. 1) GO TO 55
      DO 52 J=1,M1
      RMN(J)=1.
   52 R(J)=1.
      GO TO 56
   55 DO 50 J=2,M1
   50 R(J)=R(J-1)+YDIF(J)
      RMN(2)=R(1)
      DO 60 J=3,M2
   60 RMN(J)=RMN(J-1)+YCV(J-1)
      RMN(M1)=R(M1)
   56 CONTINUE
      DO 57 J=1,M1
      SX(J)=1.
      SXMN(J)=1.
      IF(MODE .NE. 3) GO TO 57
      SX(J)=R(J)
      IF(J .NE. 1) SXMN(J)=RMN(J)
   57 CONTINUE
      DO 62 J=2,M2
      YCVR(J)=R(J)*YCV(J)
      ARX(J)=YCVR(J)
      IF (MODE .NE. 3) GO TO 62
      ARX(J)=YCV(J)
   62 CONTINUE
      DO 64 J=4,M3
   64 YCVRS(J)=0.5*(R(J)+R(J-1))*YDIF(J)
      YCVRS(3)=0.5*(R(3)+R(1))*YCVS(3)
      YCVRS(M2)=0.5*(R(M1)+R(M3))*YCVS(M2)
      IF(MODE .NE. 2) GO TO 67
      DO 65 J=3,M3
      ARXJ(J)=0.25*(1.+RMN(J)/R(J))*ARX(J)
   65 ARXJP(J)=ARX(J)-ARXJ(J)
      GO TO 68
   67 DO 66 J=3,M3
      ARXJ(J)=0.5*ARX(J)
   66 ARXJP(J)=ARXJ(J)
   68 ARXJP(2)=ARX(2)
      ARXJ(M2)=ARX(M2)
      DO 70 J=3,M3
      FV(J)=ARXJP(J)/ARX(J)
   70 FVP(J)=1.-FV(J)
      DO 85 I=3,L2
      FX(I)=0.5*XCV(I-1)/XDIF(I)
   85 FXM(I)=1.-FX(I)
      FX(2)=0.
      FXM(2)=1.
      FX(L1)=1.
      FXM(L1)=0.
      DO 90 J=3,M2
      FY(J)=0.5*YCV(J-1)/YDIF(J)
   90 FYM(J)=1.-FY(J)
      FY(2)=0.
      FYM(2)=1.
      FY(M1)=1.
      FYM(M1)=0.
*---CON,AP,U,V,RHO,PC AND P ARRAYS ARE INITIALIZED HERE----
      DO 95 J=1,M1
      DO 95 I=1,L1
      PC(I,J)=0.
      U(I,J)=0.
      V(I,J)=0.
      CON(I,J)=0.
      AP(I,J)=0.
      RHO(I,J)=RHOCON
      P(I,J)=0.
   95 CONTINUE
      IF(MODE .EQ. 1) WRITE(8,1)
      IF(MODE .EQ. 2) WRITE(8,2)
      IF(MODE .EQ. 3) WRITE(8,3)
      WRITE(8,4)
      RETURN
*----------------------------------------------------------------------
      ENTRY SETUP2
*---COEFFICIENTS FOR THE U EQUATION----
      NF=1
      IF(.NOT. LSOLVE(NF)) GO TO 100
      IST=3
      JST=2
      CALL GAMSOR
      REL=1.-RELAX(NF)
      DO 102 I=3,L2
      FL=XCVI(I)*V(I,2)*RHO(I,1)
      FLM=XCVIP(I-1)*V(I-1,2)*RHO(I-1,1)
      FLOW=R(1)*(FL+FLM)
      DIFF=R(1)*(XCVI(I)*GAM(I,1)+XCVIP(I-1)*GAM(I-1,1))/YDIF(2)
      CALL DIFLOW
  102 AJM(I,2)=ACOF+AMAX1(0.,FLOW)
      DO 103 J=2,M2
      FLOW=ARX(J)*U(2,J)*RHO(1,J)
      DIFF=ARX(J)*GAM(1,J)/(XCV(2)*SX(J))
      CALL DIFLOW
      AIM(3,J)=ACOF+AMAX1(0.,FLOW)
      DO 103 I=3,L2
      IF(I .EQ. L2) GO TO 104
      FL=U(I,J)*(FX(I)*RHO(I,J)+FXM(I)*RHO(I-1,J))
      FLP=U(I+1,J)*(FX(I+1)*RHO(I+1,J)+FXM(I+1)*RHO(I,J))
      FLOW=ARX(J)*0.5*(FL+FLP)
      DIFF=ARX(J)*GAM(I,J)/(XCV(I)*SX(J))
      GO TO 105
  104 FLOW=ARX(J)*U(L1,J)*RHO(L1,J)
      DIFF=ARX(J)*GAM(L1,J)/(XCV(L2)*SX(J))
  105 CALL DIFLOW
      AIM(I+1,J)=ACOF+AMAX1(0.,FLOW)
      AIP(I,J)=AIM(I+1,J)-FLOW
      IF (J .EQ. M2) GOTO 106
      FL=XCVI(I)*V(I,J+1)*(FY(J+1)*RHO(I,J+1)+FYM(J+1)*RHO(I,J))
      FLM=XCVIP(I-1)*V(I-1,J+1)*(FY(J+1)*RHO(I-1,J+1)+FYM(J+1)*
     &   RHO(I-1,J))
      GM=GAM(I,J)*GAM(I,J+1)/(YCV(J)*GAM(I,J+1)+YCV(J+1)*GAM(I,J)+
     &   1.0E-30)*XCVI(I)
      GMM=GAM(I-1,J)*GAM(I-1,J+1)/(YCV(J)*GAM(I-1,J+1)+YCV(J+1)*
     &   GAM(I-1,J)+1.E-30)*XCVIP(I-1)
      DIFF=RMN(J+1)*2.*(GM+GMM)
      GO TO 107
  106 FL=XCVI(I)*V(I,M1)*RHO(I,M1)
      FLM=XCVIP(I-1)*V(I-1,M1)*RHO(I-1,M1)
      DIFF=R(M1)*(XCVI(I)*GAM(I,M1)+XCVIP(I-1)*GAM(I-1,M1))/YDIF(M1)
  107 FLOW=RMN(J+1)*(FL+FLM)
      CALL DIFLOW
      AJM(I,J+1)=ACOF+AMAX1(0.,FLOW)
      AJP(I,J)=AJM(I,J+1)-FLOW
      VOL=YCVR(J)*XCVS(I)
      APT=(RHO(I,J)*XCVI(I)+RHO(I-1,J)*XCVIP(I-1))
     &   /(XCVS(I)*DT)
      AP(I,J)=AP(I,J)-APT
      CON(I,J)=CON(I,J)+APT*FOLD(I,J,NF)
      AP(I,J)=(-AP(I,J)*VOL+AIP(I,J)+AIM(I,J)+AJP(I,J)+AJM(I,J))
     &   /RELAX(NF)
      CON(I,J)=CON(I,J)*VOL+REL*AP(I,J)*U(I,J)
      DU(I,J)=VOL/(XDIF(I)*SX(J))
      CON(I,J)=CON(I,J)+DU(I,J)*(P(I-1,J)-P(I,J))
      DU(I,J)=DU(I,J)/AP(I,J)
  103 CONTINUE
      CALL SOLVE
  100 CONTINUE
*---COEFFICIENTS FOR THE V EQUATION----
      NF=2
      IF(.NOT. LSOLVE(NF)) GO TO 200
      IST=2
      JST=3
      CALL GAMSOR
      REL=1.-RELAX(NF)
      DO 202 I=2,L2
      AREA=R(1)*XCV(I)
      FLOW=AREA*V(I,2)*RHO(I,1)
      DIFF=AREA*GAM(I,1)/YCV(2)
      CALL DIFLOW
  202 AJM(I,3)=ACOF+AMAX1(0.,FLOW)
      DO 203 J=3,M2
      FL=ARXJ(J)*U(2,J)*RHO(1,J)
      FLM=ARXJP(J-1)*U(2,J-1)*RHO(1,J-1)
      FLOW=FL+FLM
      DIFF=(ARXJ(J)*GAM(1,J)+ARXJP(J-1)*GAM(1,J-1))/(XDIF(2)*SXMN(J))
      CALL DIFLOW
      AIM(2,J)=ACOF+AMAX1(0.,FLOW)
      DO 203 I=2,L2
      IF(I .EQ. L2)GO TO 204
      FL=ARXJ(J)*U(I+1,J)*(FX(I+1)*RHO(I+1,J)+FXM(I+1)*RHO(I,J))
      FLM=ARXJP(J-1)*U(I+1,J-1)*(FX(I+1)*RHO(I+1,J-1)+FXM(I+1)*
     &   RHO(I,J-1))
      GM=GAM(I,J)*GAM(I+1,J)/(XCV(I)*GAM(I+1,J)+XCV(I+1)*GAM(I,J)+
     &   1.E-30)*ARXJ(J)
      GMM=GAM(I,J-1)*GAM(I+1,J-1)/(XCV(I)*GAM(I+1,J-1)+XCV(I+1)*
     &   GAM(I,J-1)+1.0E-30)*ARXJP(J-1)
      DIFF=2.*(GM+GMM)/SXMN(J)
      GO TO 205
  204 FL=ARXJ(J)*U(L1,J)*RHO(L1,J)
      FLM=ARXJP(J-1)*U(L1,J-1)*RHO(L1,J-1)
      DIFF=(ARXJ(J)*GAM(L1,J)+ARXJP(J-1)*GAM(L1,J-1))/(XDIF(L1)*SXMN(J))
  205 FLOW=FL+FLM
      CALL DIFLOW
      AIM(I+1,J)=ACOF+AMAX1(0.,FLOW)
      AIP(I,J)=AIM(I+1,J)-FLOW
      IF(J .EQ. M2) GO TO 206
      AREA=R(J)*XCV(I)
      FL=V(I,J)*(FY(J)*RHO(I,J)+FYM(J)*RHO(I,J-1))*RMN(J)
      FLP=V(I,J+1)*(FY(J+1)*RHO(I,J+1)+FYM(J+1)*RHO(I,J))*RMN(J+1)
      FLOW=(FV(J)*FL+FVP(J)*FLP)*XCV(I)
      DIFF=AREA*GAM(I,J)/YCV(J)
      GO TO 207
  206 AREA=R(M1)*XCV(I)
      FLOW=AREA*V(I,M1)*RHO(I,M1)
      DIFF=AREA*GAM(I,M1)/YCV(M2)
  207 CALL DIFLOW
      AJM(I,J+1)=ACOF+AMAX1(0.,FLOW)
      AJP(I,J)=AJM(I,J+1)-FLOW
      VOL=YCVRS(J)*XCV(I)
      SXT=SX(J)
      IF(J .EQ. M2) SXT=SX(M1)
      SXB=SX(J-1)
      IF(J .EQ. 3) SXB=SX(1)
      APT=(ARXJ(J)*RHO(I,J)*0.5*(SXT+SXMN(J))+ARXJP(J-1)*RHO(I,J-1)*
     &   0.5*(SXB+SXMN(J)))/(YCVRS(J)*DT)
      AP(I,J)=AP(I,J)-APT
      CON(I,J)=CON(I,J)+APT*FOLD(I,J,NF)
      AP(I,J)=(-AP(I,J)*VOL+AIP(I,J)+AIM(I,J)+AJP(I,J)+AJM(I,J))
     &   /RELAX(NF)
      CON(I,J)=CON(I,J)*VOL+REL*AP(I,J)*V(I,J)
      DV(I,J)=VOL/YDIF(J)
      CON(I,J)=CON(I,J)+DV(I,J)*(P(I,J-1)-P(I,J))
      DV(I,J)=DV(I,J)/AP(I,J)
*---------------------------------------------------------
      PRELAX(I,J)=1-((AIP(I,J)+AIM(I,J)+AJP(I,J)+AJM(I,J))/AP(I,J))
*-----------------------------------------------------------
  203 CONTINUE
      CALL SOLVE
  200 CONTINUE
*---COEFIICIENTS FOR THE PRESSURE CORRECTION EQUATION----
      NF=3
      IF(.NOT. LSOLVE(NF)) GO TO 500
      IST=2
      JST=2
      CALL GAMSOR
      SMAX=0.
      SSUM=0.
      DO 410 J=2,M2
      DO 410 I=2,L2
      VOL=YCVR(J)*XCV(I)
  410 CON(I,J)=CON(I,J)*VOL
      DO 402 I=2,L2
      ARHO=R(1)*XCV(I)*RHO(I,1)
      CON(I,2)=CON(I,2)+ARHO*V(I,2)
  402 AJM(I,2)=0.
      DO 403 J=2,M2
      ARHO=ARX(J)*RHO(1,J)
      CON(2,J)=CON(2,J)+ARHO*U(2,J)
      AIM(2,J)=0.
      DO 403 I=2,L2
      IF(I .EQ. L2) GO TO 404
      ARHO=ARX(J)*(FX(I+1)*RHO(I+1,J)+FXM(I+1)*RHO(I,J))
      FLOW=ARHO*U(I+1,J)
      CON(I,J)=CON(I,J)-FLOW
      CON(I+1,J)=CON(I+1,J)+FLOW
      AIP(I,J)=ARHO*DU(I+1,J)
      AIM(I+1,J)=AIP(I,J)
      GO TO 405
  404 ARHO=ARX(J)*RHO(L1,J)
      CON(I,J)=CON(I,J)-ARHO*U(L1,J)
      AIP(I,J)=0.
  405 IF(J .EQ. M2) GO TO 406
      ARHO=RMN(J+1)*XCV(I)*(FY(J+1)*RHO(I,J+1)+FYM(J+1)*RHO(I,J))
      FLOW=ARHO*V(I,J+1)
      CON(I,J)=CON(I,J)-FLOW
      CON(I,J+1)=CON(I,J+1)+FLOW
      AJP(I,J)=ARHO*DV(I,J+1)
      AJM(I,J+1)=AJP(I,J)
      GO TO 407
  406 ARHO=RMN(M1)*XCV(I)*RHO(I,M1)
      CON(I,J)=CON(I,J)-ARHO*V(I,M1)
      AJP(I,J)=0.
  407 AP(I,J)=AIP(I,J)+AIM(I,J)+AJP(I,J)+AJM(I,J)
      PC(I,J)=0.
      SMAX=AMAX1(SMAX,ABS(CON(I,J)))
      SSUM=SSUM+CON(I,J)
  403 CONTINUE
      CALL SOLVE
*---COMEE HERE TO CORRECT THE PRESSURE AND VELOCITIES

*--------------------------------------------------------------------
       RELAX(NP)=REV_P
*      RELAX(NP)=MAXVAL(PRELAX)
*--------------------------------------------------------------------
      DO 501 J=2,M2
      DO 501 I=2,L2
      P(I,J)=P(I,J)+PC(I,J)*RELAX(NP)
      IF(I .NE. 2) U(I,J)=U(I,J)+DU(I,J)*(PC(I-1,J)-PC(I,J))
      IF(J .NE. 2) V(I,J)=V(I,J)+DV(I,J)*(PC(I,J-1)-PC(I,J))
  501 CONTINUE
  500 CONTINUE
*---COEFFICIENTS FOR OTHER EQUATIONS----
      IST=2
      JST=2
      DO 600 N=4,NFMAX
      NF=N
      IF(.NOT. LSOLVE(NF)) GO TO 600
      CALL GAMSOR
      REL=1.-RELAX(NF)
      DO 602 I=2,L2
      AREA=R(1)*XCV(I)
      FLOW=AREA*V(I,2)*RHO(I,1)
      DIFF=AREA*GAM(I,1)/YDIF(2)
      CALL DIFLOW
  602 AJM(I,2)=ACOF+AMAX1(0.,FLOW)
      DO 603 J=2,M2
      FLOW=ARX(J)*U(2,J)*RHO(1,J)
      DIFF=ARX(J)*GAM(1,J)/(XDIF(2)*SX(J))
      CALL DIFLOW
      AIM(2,J)=ACOF+AMAX1(0.,FLOW)
      DO 603 I=2,L2
      IF(I .EQ. L2) GO TO 604
      FLOW=ARX(J)*U(I+1,J)*(FX(I+1)*RHO(I+1,J)+FXM(I+1)*RHO(I,J))
      DIFF=ARX(J)*2.*GAM(I,J)*GAM(I+1,J)/((XCV(I)*GAM(I+1,J)+
     &   XCV(I+1)*GAM(I,J)+1.0E-30)*SX(J))
      GO TO 605
  604 FLOW=ARX(J)*U(L1,J)*RHO(L1,J)
      DIFF=ARX(J)*GAM(L1,J)/(XDIF(L1)*SX(J))
  605 CALL DIFLOW
      AIM(I+1,J)=ACOF+AMAX1(0.,FLOW)
      AIP(I,J)=AIM(I+1,J)-FLOW
      AREA=RMN(J+1)*XCV(I)
      IF(J .EQ. M2) GO TO 606
      FLOW=AREA*V(I,J+1)*(FY(J+1)*RHO(I,J+1)+FYM(J+1)*RHO(I,J))
      DIFF=AREA*2.*GAM(I,J)*GAM(I,J+1)/(YCV(J)*GAM(I,J+1)+
     &   YCV(J+1)*GAM(I,J)+1.0E-30)
      GO TO 607
  606 FLOW=AREA*V(I,M1)*RHO(I,M1)
      DIFF=AREA*GAM(I,M1)/YDIF(M1)
  607 CALL DIFLOW
      AJM(I,J+1)=ACOF+AMAX1(0.,FLOW)
      AJP(I,J)=AJM(I,J+1)-FLOW
      VOL=YCVR(J)*XCV(I)
      APT=RHO(I,J)/DT
      AP(I,J)=AP(I,J)-APT
      CON(I,J)=CON(I,J)+APT*FOLD(I,J,NF)
      AP(I,J)=(-AP(I,J)*VOL+AIP(I,J)+AIM(I,J)+AJP(I,J)+AJM(I,J))
     &   /RELAX(NF)
      CON(I,J)=CON(I,J)*VOL+REL*AP(I,J)*F(I,J,NF)
  603 CONTINUE
      CALL SOLVE
  600 CONTINUE
      ITER=ITER+1
*---------------------------------------------------------------------
      C_TIME=C_TIME+1
*------------------------------------------------------------------------
      IF(ITER .GE. LAST)THEN
       LSTOP=.TRUE.
       END IF
      RETURN
      END
*-----------------------------------------------------------------------
      SUBROUTINE SUPPLY
************************************************************************
	PARAMETER (NX=22,NY=22,NXY=22)
      CHARACTER (LEN=20) TITLE
      LOGICAL LSOLVE,LPRINT,LBLK,LSTOP,LEND
      COMMON F(NX,NY,13),FOLD(NX,NY,13),P(NX,NY),RHO(NX,NY),GAM(NX,NY),
     &  CON(NX,NY),AIP(NX,NY),AIM(NX,NY),AJP(NX,NY),AJM(NX,NY),
     &  AP(NX,NY),X(NXY),XU(NXY),XDIF(NXY),XCV(NXY),XCVS(NXY),
     &  Y(NXY),YV(NXY),YDIF(NXY),YCV(NXY),YCVS(NXY),
     &  YCVR(NXY),YCVRS(NXY),ARX(NXY),ARXJ(NXY),ARXJP(NXY),
     &  R(NXY),RMN(NXY),SX(NXY),SXMN(NXY),XCVI(NXY),XCVIP(NXY)
      COMMON DU(NX,NY),DV(NX,NY),FV(NXY),FVP(NXY),
     &  FX(NXY),FXM(NXY),FY(NXY),FYM(NXY),PT(NXY),QT(NXY)
      COMMON /INDX/NF,NFMAX,NP,NRHO,NGAM,L1,L2,L3,M1,M2,M3,
     &  IST,JST,ITER,LAST,TITLE(13),RELAX(13),TIME,DT,XL,YL,
     &  IPREF,JPREF,LSOLVE(10),LPRINT(13),LBLK(10),MODE,NTIMES(10),
     &  RHOCON,PRELAX(NX,NY)
      COMMON/CNTL/LSTOP,LEND
      COMMON ITIME,MTIME
      COMMON/modl/REV_P,C_TIME      
      DIMENSION U(NX,NY),V(NX,NY),PC(NX,NY)
      EQUIVALENCE (F(1,1,1),U(1,1)),(F(1,1,2),V(1,1)),(F(1,1,3),PC(1,1))
************************************************************************
   10 FORMAT(1X,26(1H*),3X,A10,3X,26(1H*))
   20 FORMAT(1X,4H I =,I6,6I9)
   30 FORMAT(1X,1HJ)
   40 FORMAT(1X,I2,3X,1P7E9.2)
   50 FORMAT(1X,1H )
   51 FORMAT(1X,' I =',2X,7(I4,5X))
   52 FORMAT(1X,' X =',1P7E9.2)
   53 FORMAT(1X,'TH =',1P7E9.2)
   54 FORMAT(1X,'J =',2X,7(I4,5X))
   55 FORMAT(1X,'Y =',1P7E9.2)
************************************************************************
      ENTRY UGRID
      XU(2)=0.
      DX=XL/FLOAT(L1-2)
      DO 1 I=3,L1
    1 XU(I)=XU(I-1)+DX
      YV(2)=0.
      DY=YL/FLOAT(M1-2)
      DO 2 J=3,M1
    2 YV(J)=YV(J-1)+DY
      RETURN
************************************************************************
      ENTRY PRINT_RESULT
      IF(.NOT. LPRINT(3)) GO TO 80
*---CALCULATE THE STREAM FUNTION----------------------------------------
      F(2,2,3)=0.
      DO 82 I=2,L1
      IF(I .NE. 2) F(I,2,3)=F(I-1,2,3)-RHO(I-1,1)*V(I-1,2)
     &   *R(1)*XCV(I-1)
      DO 82 J=3,M1
      RHOM=FX(I)*RHO(I,J-1)+FXM(I)*RHO(I-1,J-1)
   82 F(I,J,3)=F(I,J-1,3)+RHOM*U(I,J-1)*ARX(J-1)
   80 CONTINUE
*
      IF( .NOT. LPRINT(NP)) GO TO 90
*
*---CONSTRUCT BOUNDARY PRESSURES BY EXTRAPOLATION
      DO 91 J=2,M2
      P(1,J)=(P(2,J)*XCVS(3)-P(3,J)*XDIF(2))/XDIF(3)
   91 P(L1,J)=(P(L2,J)*XCVS(L2)-P(L3,J)*XDIF(L1))/XDIF(L2)
      DO 92 I=2,L2
      P(I,1)=(P(I,2)*YCVS(3)-P(I,3)*YDIF(2))/YDIF(3)
   92 P(I,M1)=(P(I,M2)*YCVS(M2)-P(I,M3)*YDIF(M1))/YDIF(M2)
      P(1,1)=P(2,1)+P(1,2)-P(2,2)
      P(L1,1)=P(L2,1)+P(L1,2)-P(L2,2)
      P(1,M1)=P(2,M1)+P(1,M2)-P(2,M2)
      P(L1,M1)=P(L2,M1)+P(L1,M2)-P(L2,M2)
      PREF=P(IPREF,JPREF)
      DO 93 J=1,M1
      DO 93 I=1,L1
   93 P(I,J)=P(I,J)-PREF
   90 CONTINUE
*
      WRITE (8,50)
      IEND=0
  301 IF(IEND .EQ. L1) GO TO 310
      IBEG=IEND+1
      IEND=IEND+7
      IEND=MIN0(IEND,L1)
      WRITE (8,50)
      WRITE(8,51) (I,I=IBEG,IEND)
      IF(MODE .EQ. 3) GO TO 302
      WRITE(8,52) (X(I),I=IBEG,IEND)
      GO TO 303
  302 WRITE (8,53) (X(I),I=IBEG,IEND)
  303 GO TO 301
  310 JEND=0
      WRITE(8,50)
  311 IF(JEND .EQ. M1) GO TO 320
      JBEG=JEND+1
      JEND=JEND+7
      JEND=MIN0(JEND,M1)
      WRITE(8,50)
      WRITE(8,54) (J,J=JBEG,JEND)
      WRITE(8,55) (Y(J),J=JBEG,JEND)
      GO TO 311
  320 CONTINUE
*
      DO 999 N=1,NGAM
      NF=N
      IF(.NOT. LPRINT(NF)) GO TO 999
      WRITE(8,50)
      WRITE(8,10) TITLE(NF)
      IFST=1
      JFST=1
      IF(NF .EQ. 1 .OR. NF .EQ. 3) IFST=2
      IF(NF .EQ. 2 .OR. NF .EQ. 3) JFST=2
      IBEG=IFST-7
  110 CONTINUE
      IBEG=IBEG+7
      IEND=IBEG+6
      IEND=MIN0(IEND,L1)
      WRITE(8,50)
      WRITE(8,20) (I,I=IBEG,IEND)
      WRITE(8,30)
      JFL=JFST+M1
      DO 115 JJ=JFST,M1
      J=JFL-JJ
      WRITE(8,40) J,(F(I,J,NF),I=IBEG,IEND)
  115 CONTINUE
      IF(IEND .LT. L1) GO TO 110
  999 CONTINUE
      RETURN
      END