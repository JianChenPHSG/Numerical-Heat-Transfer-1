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
      PARAMETER (NX=26,NY=26,NXY=26)
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
     &  RHOCON
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
*            Example-1: Transient Navier-Stokes Equations             *
*                                                                     *
***********************************************************************
      ENTRY GRID
        XL=2.
        YL=2.
        L1=26
        M1=26
	  CALL UGRID
      RETURN
C----------------------------------------------------------------------
      ENTRY START
        MODE=1

        LAST=400
        RSMAX=1.E-7
	          
        TSTART=0.
	  DT=0.001
	  MTIME=100
	  ITIME=1
	  TIME=TSTART
	  LEND=.FALSE.

	  DO I=1,3
	    LSOLVE(I)=.TRUE.
          LPRINT(I)=.TRUE.
	  END DO
	  RELAX(1)=0.7
	  RELAX(2)=0.7
	  TITLE(1)='.U EQ.'
	  TITLE(2)='.V EQ.'

        PAI=3.1415926       
	  DO I=2,L1
	    DO J=1,M1
	      U(I,J)=-COS(PAI*XU(I))*SIN(PAI*Y(J))*EXP(-2*PAI*PAI*TIME)
	    END DO
	  END DO

	  DO I=1,L1
	    DO J=2,M1
	      V(I,J)=SIN(PAI*X(I))*COS(PAI*YV(J))*EXP(-2*PAI*PAI*TIME)
	    END DO
	  END DO

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
	  DO J=1,M1
	    DO I=1,L1
	      RHO(I,J)=1.
	    END DO
	  END DO
      RETURN
C----------------------------------------------------------------------
      ENTRY BOUND
	  PAI=3.1415926
	  DO J=1,M1
	    U(2,J)=-SIN(PAI*Y(J))*EXP(-2*PAI*PAI*TIME)
	    U(L1,J)=U(2,J)
	    V(1,J)=0
	    V(L1,J)=0
	  END DO

	  DO I=1,L1
	    U(I,1)=0
	    U(I,M1)=0
	    V(I,2)=SIN(PAI*X(I))*EXP(-2*PAI*PAI*TIME)
          V(I,M1)=V(I,2)
	  END DO
      RETURN
C----------------------------------------------------------------------
      ENTRY OUTPUT
	  PAI=3.1415926      
        TEMPU=-COS(PAI*XU(14))*SIN(PAI*Y(14))*EXP(-2.*PAI*PAI*TIME)
        TEMPD=ABS((TEMPU-U(14,14))/MAX(ABS(TEMPU),ABS(U(14,14))))
        IF(ITIME.EQ.1) WRITE(8,400)
        IF(ITIME.GE.1) GOTO 401

  400   FORMAT(13X,'TIME',14X,'U(14,14)',12X,'TU(14,14)',12X,'ERROR_U')
  401   WRITE(8,403) TIME,U(14,14),TEMPU,TEMPD
  403   FORMAT(1X,4F20.8)  
  
        IF(ITIME.EQ.MTIME) CALL PRINT_RESULT        
      
        IF(ITIME.EQ.50) THEN
          OPEN(1,FILE='VEL_U.DAT')
          WRITE(1,*) 'TITLE="VEL_U"'
          WRITE(1,*) 'VARIABLES="X","Y","U"'
          WRITE(1,*) 'ZONE I=',L1-1,',J=',M1
          DO J=1,M1  
           DO I=2,L1
             WRITE(1,*) XU(I),Y(J),U(I,J)
            END DO
          END DO
          CLOSE(1)
        
          OPEN(1,FILE='VEL_TU.DAT')
          WRITE(1,*) 'TITLE="VEL_TU"'
          WRITE(1,*) 'VARIABLES="X","Y","U"'
          WRITE(1,*) 'ZONE I=',L1-1,',J=',M1
          DO J=1,M1  
            DO I=2,L1
              TEMPU=-COS(PAI*XU(I))*SIN(PAI*Y(J))*EXP(-2.*PAI*PAI*TIME)
              WRITE(1,*) XU(I),Y(J),TEMPU
            END DO
          END DO
          CLOSE(1)
        END IF      
      RETURN
C----------------------------------------------------------------------
      ENTRY GAMSOR
***********************************************************************
C For the following source terms:
C d^2(u)/dx^2+(d/dy)(dv/dx) in u-Equation, and
C d^2(v)/dy^2+(d/dx)(du/dy) in v-Equation
C ISOURCE=1: The preceding source terms appear in u- or v-Equ.
C ISOURCE=2: The preceding source terms DO NOT appear in u- or v-Equ.
***********************************************************************
        ISOURCE=1

        DO I=1,L1
          DO J=1,M1
            GAM(I,J)=1.
          END DO
        END DO

        IF(ISOURCE.EQ.1) THEN	  
          IF (NF.EQ.1) THEN
            DO I=3,L2
              DO J=2,M2	
                CON(I,J)=((U(I+1,J)-U(I,J))/XCV(I)-
     &                    (U(I,J)-U(I-1,J))/XCV(I-1))/XDIF(I)
                CON(I,J)=CON(I,J)+((V(I,J+1)-V(I-1,J+1))/XDIF(I)-
     &                             (V(I,J)-V(I-1,J))/XDIF(I))/YCV(J)
              END DO
            END DO
          END IF
        
          IF (NF.EQ.2) THEN
            DO J=3,M2
              DO I=2,L2
                CON(I,J)=((V(I,J+1)-V(I,J))/YCV(J)-
     &                    (V(I,J)-V(I,J-1))/YCV(J-1))/YDIF(J)
                CON(I,J)=CON(I,J)+((U(I+1,J)-U(I+1,J-1))/YDIF(J)-
     &                             (U(I,J)-U(I,J-1))/YDIF(J))/XCV(I)
              END DO
            END DO
          END IF
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