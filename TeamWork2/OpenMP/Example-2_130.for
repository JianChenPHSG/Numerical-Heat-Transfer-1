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

      PARAMETER (NX=130,NY=130,NXY=130)
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
*            Example-2: Lid-Driven Cavity Flow                        *
*                                                                     *
***********************************************************************
      ENTRY GRID
        XL=1.
        YL=1.
        L1=130
        M1=130
	  CALL UGRID
      RETURN
C----------------------------------------------------------------------
      ENTRY START
        MODE=1

        LAST=200
        RSMAX=1E-9
	          
        TSTART=0.
	  DT=1.E10
	  MTIME=1
	  ITIME=1
	  TIME=TSTART
	  LEND=.FALSE.

	  DO I=1,3
	    LSOLVE(I)=.TRUE.
          LPRINT(I)=.TRUE.
	  END DO
	  RELAX(1)=0.5
	  RELAX(2)=0.5
	  TITLE(1)='.U EQ.'
	  TITLE(2)='.V EQ.'

        DO I=1,L1
	    DO J=1,M1
	      U(I,J)=0.
	      V(I,J)=0.
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
        UMAX=1.
	  DO J=1,M1
	    U(2,J)=0.
	    U(L1,J)=0.
	    V(1,J)=0.
	    V(L1,J)=0.
	  END DO

	  DO I=1,L1
	    U(I,1)=0.
	    U(I,M1)=UMAX
	    V(I,2)=0.
          V(I,M1)=0.
	  END DO
      RETURN
C----------------------------------------------------------------------
      ENTRY OUTPUT
        IF(ITER.EQ.1) WRITE(8,400)
        IF(ITER.GE.1) THEN
          WRITE(8,403) ITER,U(66,66),SMAX
          WRITE(*,403) ITER,U(66,66),SMAX          
        END IF

  400   FORMAT(13X,'ITER',14X,'U(66,66)',12X,'SMAX')
  403   FORMAT(1X,I8,2F20.10)  
  
        IF(LSTOP.OR.(ITER.EQ.LAST)) CALL PRINT_RESULT        
      
        IF(LSTOP.OR.(ITER.EQ.LAST)) THEN
          OPEN(1,FILE='VEL_U.DAT')
            DO J=1,M1  
            WRITE(1,*) U(66,J),Y(J)
          END DO
          CLOSE(1)

          OPEN(1,FILE='VEL_V.DAT')
          DO I=1,L1  
            WRITE(1,*) X(I),V(I,66)
          END DO
          CLOSE(1)
          
          OPEN(1,FILE='STREAM_FUNCTION.DAT')
          WRITE(1,*) 'TITLE="SF"'
          WRITE(1,*) 'VARIABLES="X","Y","SF"'
          WRITE(1,*) 'ZONE I=',L1,',J=',M1
          DO J=1,M1  
           DO I=1,L1
             WRITE(1,*) X(I),Y(J),F(I,J,3)
            END DO
          END DO
          CLOSE(1) 
          Pause       
        END IF         
      RETURN
C----------------------------------------------------------------------
      ENTRY GAMSOR
        RE=100.
        DO I=1,L1
          DO J=1,M1
            GAM(I,J)=1./RE
          END DO
        END DO
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
          IF (SMAX<RSMAX) THEN 
             LSTOP=.TRUE. 
             !CALL CPU_TIME(FINISH2)  
             !CALL SYSTEM_CLOCK(FINISH1) 
             !WRITE(* , *) "SYSTEM_TIME运算时间为(s):", FINISH1-START1
             !WRITE(* , *) "CPU_TIME运算时间为(s):", FINISH2-START2
             print *, 'hello'
          END IF
        END IF
      RETURN
C----------------------------------------------------------------------
      END