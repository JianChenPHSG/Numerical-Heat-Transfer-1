CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE USER
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  This is the user part of the program. In this part, the users need
C  to specify all of the aspects related with the individual physical 
C  problems, including:
C  (1) the governing equations to be solved;
C  (2) the calculation domain, and the corresponding mesh generation;
C  (3) the physical properties of the materials;
C  (4) the boundary and/or initial conditions (if any);
C  (5) the source terms;
C  (6) the time interval if any, and the other parameters for the
C      iterations, i.e., the relaxation factors, the criteria for
C      ending the internal/external iterations.
***********************************************************************
      PARAMETER (NX=122,NY=22,NXY=122)
      CHARACTER (LEN=20) TITLE
      LOGICAL LSOLVE,LPRINT,LBLK,LSTOP,LEND
      COMMON F(NX,NY,13),FOLD(NX,NY,13),P(NX,NY),RHO(NX,NY),GAM(NX,NY),
     & CON(NX,NY),AIP(NX,NY),AIM(NX,NY),AJP(NX,NY),AJM(NX,NY),AP(NX,NY),
     & X(NXY),XU(NXY),XDIF(NXY),XCV(NXY),XCVS(NXY),
     & Y(NXY),YV(NXY),YDIF(NXY),YCV(NXY),YCVS(NXY),
     & YCVR(NXY),YCVRS(NXY),ARX(NXY),ARXJ(NXY),ARXJP(NXY),
     & R(NXY),RMN(NXY),SX(NXY),SXMN(NXY),XCVI(NXY),XCVIP(NXY)
      COMMON DU(NX,NY),DV(NX,NY),FV(NXY),FVP(NXY),
     & FX(NXY),FXM(NXY),FY(NXY),FYM(NXY),PT(NXY),QT(NXY)
      COMMON/INDX/NF,NFMAX,NP,NRHO,NGAM,L1,L2,L3,M1,M2,M3,
     &  IST,JST,ITER,LAST,TITLE(13),RELAX(13),TIME,DT,XL,YL,
     &  IPREF,JPREF,LSOLVE(10),LPRINT(13),LBLK(10),MODE,NTIMES(10),
     &  RHOCON,PRELAX(NX,NY)
      COMMON/CNTL/LSTOP,LEND
      COMMON ITIME
      COMMON/SORC/SMAX,SSUM,RSMAX
      COMMON/COEF/FLOW,DIFF,ACOF
      DIMENSION U(NX,NY),V(NX,NY),PC(NX,NY)
      EQUIVALENCE (F(1,1,1),U(1,1)),(F(1,1,2),V(1,1)),(F(1,1,3),PC(1,1))
      DIMENSION TH(NXY),THU(NXY),THDIF(NXY),THCV(NXY),THCVS(NXY)
      EQUIVALENCE(X,TH),(XU,THU),(XDIF,THDIF),(XCV,THCV),
     &  (XCVS,THCVS),(XL,THL)
***********************************************************************
*                                                                     *
*        Example-6: THERMAL PLASMA TUBE FLOW                          *
*                                                                     *
***********************************************************************
	DIMENSION T(NX,NY),AH(NX,NY),CP(NX,NY),VIS(NX,NY),AK(NX,NY)
	DIMENSION TDAT(120),TBRHO(120),TBVIS(120),TBK(120),TBCP(120),
     &          TBH(120),TBRAD(120),TBRHT(120)
	EQUIVALENCE (F(1,1,4),AH(1,1)),(F(1,1,5),T(1,1))
      
      ENTRY GRID
        XL=120.E-3
        YL=5.E-3
        L1=122
        M1=22
	  CALL UGRID
      RETURN
C----------------------------------------------------------------------
      ENTRY START
        MODE=2
        R(1)=0.
        ISTART=1
        
        LAST=10000
        RSMAX=1.E-8
	          
        TSTART=0.
	  DT=1.E10
	  MTIME=1
	  ITIME=1
	  TIME=TSTART
	  LEND=.FALSE.


        DO I=1,4
          LSOLVE(I)=.TRUE.
          LPRINT(I)=.TRUE.
          LBLK(I)=.FALSE.
        END DO
	  RELAX(1)=0.3
	  RELAX(2)=0.3	  
	  RELAX(4)=0.3
	  RELAX(NP)=0.3	  
	  RELAX(NRHO)=0.3
	  TITLE(1)='.VEL-U.'
	  TITLE(2)='.VEL-V.'
	  TITLE(3)='.STRF.'	  
	  TITLE(4)='.TEMP.'

	  OPEN(1,FILE='ARGON.DAT')
	    DO I=1,120
	      READ(1,100) TDAT(I),TBRHO(I),TBVIS(I),TBCP(I),
     &                  TBK(I),TBH(I),TBRAD(I)
	      TBRHT(I)=TBRHO(I)*TDAT(I)     
          END DO
	  CLOSE(1)
  100   FORMAT(7(E12.4,5X))  	  

   	  IF(ISTART.EQ.1) THEN
	    DO I=1,L1
	      DO J=1,M1
	        U(I,J)=0.
	        V(I,J)=0.
	        T(I,J)=500.
              AH(I,J)=PROP(T(I,J),TBH)
            END DO
	    END DO
	  ELSE
  	    OPEN(1,FILE='USER.DAT')
	      READ(1,*) ((U(I,J),I=1,L1),J=1,M1)
	      READ(1,*) ((V(I,J),I=1,L1),J=1,M1)
	      READ(1,*) ((P(I,J),I=1,L1),J=1,M1)
	      READ(1,*) ((T(I,J),I=1,L1),J=1,M1)
	      READ(1,*) ((AH(I,J),I=1,L1),J=1,M1)
	    CLOSE(1)
	  END IF  
	  DO K=1,NFMAX
	    DO I=1,L1
	      DO J=1,M1
	        FOLD(I,J,K)=F(I,J,K)
	      END DO
	    END DO
	  END DO
        TIME=TSTART+DT
      RETURN
C----------------------------------------------------------------------
      ENTRY DENSE
	  DO I=1,L1
	    DO J=1,M1
	      T(I,J)=ANPROP(AH(I,J),TBH,TDAT)
	      IF (T(I,J).GT.24100.) T(I,J)=24100.
	      IF (T(I,J).LT.300.) T(I,J)=300.
            DEN=PROP(T(I,J),TBRHT)/T(I,J)
	      IF (ITER.EQ.0) RHO(I,J)=DEN
	      IF (ITER.NE.0) RHO(I,J)=DEN*RELAX(NRHO)+RHO(I,J)*(1.-RELAX(NRHO))
            VIS(I,J)=PROP(T(I,J),TBVIS)
            AK(I,J)=PROP(T(I,J),TBK)
            CP(I,J)=PROP(T(I,J),TBCP)
          END DO
        END DO
      RETURN
C----------------------------------------------------------------------
      ENTRY BOUND
	  DO J=1,M1
	    U(2,J)=500.*(1.-(R(J)/YL)**2.)
	    V(1,J)=0.
	    T(1,J)=14700.*(1.-(R(J)/YL)**2.)+300.
	    AH(1,J)=PROP(T(1,J),TBH)       
          U(L1,J)=U(L2,J)
          V(L1,J)=V(L2,J)
          AH(L1,J)=AH(L2,J)
          T(L1,J)=ANPROP(AH(L1,J),TBH,TDAT)
	  END DO

	  DO I=1,L1
          U(I,1)=U(I,2)
          V(I,2)=0.
          AH(I,1)=AH(I,2)
          T(I,1)=ANPROP(AH(I,1),TBH,TDAT)
*******上边界无滑移条件	  
          U(I,M1)=0.
          V(I,M1)=0.
          T(I,M1)=300.
          AH(I,M1)=PROP(T(I,M1),TBH)
	  END DO
      RETURN
C----------------------------------------------------------------------
      ENTRY OUTPUT
        IF(ITER.GE.2000) THEN
          FLOWIN=0.
          FLOWOUT=0.
          DO J=2,M2
            FLOWIN=FLOWIN+RHO(1,J)*U(2,J)*ARX(J)
            FLOWOUT=FLOWOUT+RHO(L1,J)*U(L1,J)*ARX(J)          
          END DO          
          FACTOR=FLOWOUT/FLOWIN
          DO J=2,M2
            U(L1,J)=U(L1,J)/FACTOR
          END DO
        END IF
        
        IF(ITER.EQ.1) WRITE(8,400)
        IF(MOD(ITER,20).EQ.0) THEN
          WRITE(8,403) ITER,RELAX(NP),U(22,10),V(22,10),T(22,10),SMAX
          WRITE(*,403) ITER,RELAX(NP),U(22,10),V(22,10),T(22,10),SMAX
        END IF

  400   FORMAT(13X,'ITER',12X,'RELAX(NP)',12X,'U(62,10)',12X,'V(62,10)',
     &                      12X,'T(62,10)',12X,'SMAX')
  403   FORMAT(1X,I8,5F20.10)  
  
        IF(LSTOP.OR.(ITER.EQ.LAST)) CALL PRINT_RESULT        
      
        IF(LSTOP.OR.(ITER.EQ.LAST)) THEN
          FLOWIN=0.
          FLOWOUT=0.
          DO J=2,M2
            FLOWIN=FLOWIN+RHO(1,J)*U(2,J)*ARX(J)
            FLOWOUT=FLOWOUT+RHO(L1,J)*U(L1,J)*ARX(J)          
          END DO          
          WRITE(*,*) 'FLOWIN=',FLOWIN
          WRITE(*,*) 'FLOWOUT=',FLOWOUT
                  
          OPEN(1,FILE='STREAM_FUNCTION.DAT')
          WRITE(1,*) 'TITLE="SF"'
          WRITE(1,*) 'VARIABLES="X","Y","SF"'
          WRITE(1,*) 'ZONE I=',L1,',J=',2*M1-1
          DO J=M1,1,-1  
           DO I=1,L1
             WRITE(1,*) X(I)*1.E3,R(J)*1.E3,F(I,J,3)
            END DO
          END DO
          DO J=2,M1  
           DO I=1,L1
             WRITE(1,*) X(I)*1.E3,-R(J)*1.E3,F(I,J,3)
            END DO
          END DO
          CLOSE(1)

          OPEN(1,FILE='VEL_U.DAT')
          WRITE(1,*) 'TITLE="VEL_U"'
          WRITE(1,*) 'VARIABLES="X","Y","VEL_U"'
          WRITE(1,*) 'ZONE I=',L1,',J=',2*M1-1
          DO J=M1,1,-1  
           DO I=1,L1
             WRITE(1,*) X(I)*1.E3,R(J)*1.E3,U(I,J)
            END DO
          END DO
          DO J=2,M1  
           DO I=1,L1
             WRITE(1,*) X(I)*1.E3,-R(J)*1.E3,U(I,J)
            END DO
          END DO
          CLOSE(1)

          OPEN(1,FILE='TEMPERATURE.DAT')
          WRITE(1,*) 'TITLE="TEMP"'
          WRITE(1,*) 'VARIABLES="X","Y","TEMP"'
          WRITE(1,*) 'ZONE I=',L1,',J=',2*M1-1
          DO J=M1,1,-1  
           DO I=1,L1
             T(I,J)=ANPROP(AH(I,J),TBH,TDAT)
             WRITE(1,*) X(I)*1.E3,R(J)*1.E3,T(I,J)
            END DO
          END DO
          DO J=2,M1  
           DO I=1,L1
             T(I,J)=ANPROP(AH(I,J),TBH,TDAT)
             WRITE(1,*) X(I)*1.E3,-R(J)*1.E3,T(I,J)
            END DO
          END DO
          CLOSE(1)
          
	    OPEN(1,FILE='USEW.DAT')
	      WRITE(1,*) ((U(I,J),I=1,L1),J=1,M1)
	      WRITE(1,*) ((V(I,J),I=1,L1),J=1,M1)
	      WRITE(1,*) ((P(I,J),I=1,L1),J=1,M1)
	      WRITE(1,*) ((T(I,J),I=1,L1),J=1,M1)
	      WRITE(1,*) ((AH(I,J),I=1,L1),J=1,M1)
	    CLOSE(1)
        END IF          
      RETURN
C----------------------------------------------------------------------
*****等效扩散系数和源项
      ENTRY GAMSOR
        DO I=1,L1
	    DO J=1,M1
	      IF(NF.LE.3) THEN
	        GAM(I,J)=VIS(I,J)
	      ELSE
	        GAM(I,J)=AK(I,J)/CP(I,J)
            END IF
          END DO
        END DO
*********等效源项，认为边界上没有通量，等效到控制体中，所以对应的离散方程系数等于0        
	  IF(NF.EQ.1) THEN
          DO I=1,L1
            GAM(I,1)=0.
          END DO
          DO J=1,M1
            GAM(L1,J)=0.
          END DO
**********线性插值
          DO I=3,L2
            DO J=2,M2
              GAMP1=FX(I)*GAM(I,J+1)+FXM(I)*GAM(I-1,J+1)
              GAMP2=FX(I)*GAM(I,J)+FXM(I)*GAM(I-1,J)
              IF(J.EQ.M2) THEN
                GAMP=GAMP1
              ELSE
                GAMP=FY(J+1)*GAMP1+FYM(J+1)*GAMP2
              END IF
              GAMM1=FX(I)*GAM(I,J-1)+FXM(I)*GAM(I-1,J-1)
              IF(J.GT.2) THEN
                GAMM=FY(J)*GAMP2+FYM(J)*GAMM1
              ELSE
                GAMM=GAMM1
              END IF
              TEMP1=(GAMP*(V(I,J+1)-V(I-1,J+1))*RMN(J+1)-
     &             GAMM*(V(I,J)-V(I-1,J))*RMN(J))/YCV(J)/XCVS(I)/R(J)
              TEMP2=(GAM(I,J)*(U(I+1,J)-U(I,J))/XCV(I)-
     &             GAM(I-1,J)*(U(I,J)-U(I-1,J))/XCV(I-1))/XCVS(I)
	        CON(I,J)=TEMP1+TEMP2
	      END DO
	    END DO
	  END IF
	  
        IF(NF.EQ.2) THEN
          DO J=1,M1
           GAM(L1,J)=0.
*******只有出口可用等效源项法
          END DO

          DO I=2,L2
            DO J=3,M2
              TEMP1=(R(J)*GAM(I,J)*(V(I,J+1)-V(I,J))/YCV(J)-
     &               R(J-1)*GAM(I,J-1)*(V(I,J)-V(I,J-1))/YCV(J-1))
     &              /YCVS(J)/RMN(J)
              GAMP1=FY(J)*GAM(I+1,J)+FYM(J)*GAM(I+1,J-1)
              GAMP2=FY(J)*GAM(I,J)+FYM(J)*GAM(I,J-1)
              IF(I.EQ.L2) THEN
                GAMP=GAMP1
              ELSE
                GAMP=FX(I+1)*GAMP1+FXM(I+1)*GAMP2
              END IF
              GAMM1=FY(J)*GAM(I-1,J)+FYM(J)*GAM(I-1,J-1)
              IF(I.EQ.2) THEN
                GAMM=GAMM1
              ELSE 
                GAMM=FX(I)*GAMP2+FXM(I)*GAMM1
              END IF
              TEMP2=(GAMP*(U(I+1,J)-U(I+1,J-1))-GAMM*(U(I,J)-U(I,J-1)))
     &             /(XCV(I)*YCVS(J))
	        CON(I,J)=TEMP1+TEMP2
	        
              GAMP=FY(J)*GAM(I,J)+FYM(J)*GAM(I,J-1) 
              AP(I,J)=-2.*GAMP/(RMN(J)**2.)
            END DO
          END DO
        END IF
        
        IF(NF.EQ.4) THEN
          DO I=1,L1
            GAM(I,1)=0.
          END DO
          DO J=1,M1
            GAM(L1,J)=0.
          END DO

          DO I=2,L2
            DO J=2,M2
              T(I,J)=ANPROP(AH(I,J),TBH,TDAT)
              CON(I,J)=-PROP(T(I,J),TBRAD)
            END DO
          END DO
        END IF
      RETURN
C----------------------------------------------------------------------
      ENTRY TSTEP
        ITER=0
        LSTOP=.FALSE.
        IF (ITIME.GE.MTIME) LEND=.TRUE.
        TIME=TIME+DT
        ITIME=ITIME+1
      RETURN
C----------------------------------------------------------------------
      ENTRY STEPUP
        DO K=1,NFMAX
          DO J=1,M1
            DO I=1,L1
              FOLD(I,J,K)=F(I,J,K)
            END DO
          END DO
        END DO
      RETURN
C----------------------------------------------------------------------
      ENTRY ISSTOP
        IF(LSOLVE(1)) THEN
          IF (SMAX<RSMAX) LSTOP=.TRUE.
        END IF
      RETURN
C----------------------------------------------------------------------      
      END
C----------------------------------------------------------------------
	FUNCTION ANPROP(AH,VAL,TDAT)
	  DIMENSION VAL(120),TDAT(120)
        DO 10 I=1,120
          IF(AH.GT.VAL(I)) GOTO 10
          N=I-1
          IF(I.EQ.1) N=1
          GOTO 20
   10   CONTINUE      
   20   ANPROP=(AH-VAL(N))/(VAL(N+1)-VAL(N))*200.+TDAT(N)
	RETURN
	END
C----------------------------------------------------------------------
***正查表
	FUNCTION PROP(T,TABLE) 
	  DIMENSION TABLE(120)
	  IF(T.LT.300.) T=300.
	  IF(T.GT.24100.) T=24100.
	  J=IFIX((T-100.)/200.)
	  DELT=T-FLOAT(J)*200-100.
	  PROP=TABLE(J)+DELT/200.*(TABLE(J+1)-TABLE(J))
	RETURN 
	END