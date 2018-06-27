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
      PARAMETER (NX=1920,NY=128,NXY=1920)
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
      CHARACTER*11 FILENAME1
***********************************************************************
*                                                                     *
*            Example-3: FLOW OVER A BACKWARD-FACING STEP              *
*                                                                     *
***********************************************************************
      ENTRY GRID
        XL=30.
        YL=1.
        L1=1920
        M1=128
	  CALL UGRID
      RETURN
C----------------------------------------------------------------------
      ENTRY START
        MODE=1

        ISTART=1
        
        LAST=10000
        RSMAX=1.E-4
	          
        TSTART=0.
	  DT=1.E10
	  MTIME=1
	  ITIME=1
	  TIME=TSTART
	  LEND=.FALSE.

	  DO I=1,3
	    LSOLVE(I)=.TRUE.
          LPRINT(I)=.TRUE.
          LBLK(I)=.FALSE.
	  END DO
	  RELAX(1)=0.1
	  RELAX(2)=0.1
	  RELAX(11)=1.0	  
	  TITLE(1)='.U EQ.'
	  TITLE(2)='.V EQ.'

        IF(ISTART.EQ.1) THEN
          DO I=1,L1
	      DO J=1,M1
	        U(I,J)=0.
	        V(I,J)=0.
	      END DO
	    END DO
	  ELSE
	    OPEN(1,FILE='USER.DAT')
            DO I=1,L1
	        DO J=1,M1
	          READ(1,*) U(I,J),V(I,J)
	        END DO
	      END DO
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
	  DO J=1,M1
	    DO I=1,L1
	      RHO(I,J)=1.
	    END DO
	  END DO
      RETURN
C----------------------------------------------------------------------
      ENTRY BOUND
	  DO J=1,M1
          IF(J.LE.64) THEN
	      U(2,J)=0.
          ELSE
            U(2,J)=24.*(Y(J)-0.5)*(1.-Y(J))
          END IF
	    U(L1,J)=U(L2,J)
	    V(1,J)=0.
	    V(L1,J)=0.
	  END DO

	  DO I=1,L1
	    U(I,1)=0.
	    U(I,M1)=0.
	    V(I,2)=0.
          V(I,M1)=0.
	  END DO
      RETURN
C----------------------------------------------------------------------
      ENTRY OUTPUT
        DATA NSTART,NSTEP/20,5000/ 
        METHOD=2
        
        IF(METHOD.EQ.1) THEN
          FLOWIN=0.
          FLOWOUT=0.
          DO J=1,M1
            FLOWIN=FLOWIN+U(2,J)*ARX(J)
            FLOWOUT=FLOWOUT+U(L2,J)*ARX(J)
          END DO        
          FACTOR=FLOWOUT/FLOWIN
          DO J=1,M1
            U(L1,J)=U(L2,J)/FACTOR
          END DO
        END IF
         
        IF(LSTOP.OR.(ITER.EQ.LAST)) THEN
          WRITE(*,*) 'FLOWIN=',FLOWIN
          WRITE(*,*) 'FLOWOUT=',FLOWOUT
        END IF
                        
        IF(ITER.EQ.1) WRITE(8,400)
        IF(MOD(ITER,20).EQ.0) THEN
          WRITE(8,403) ITER,U(194,8),SMAX
          WRITE(*,403) ITER,U(194,8),SMAX          
        END IF

  400   FORMAT(13X,'ITER',14X,'U(194,8)',12X,'SMAX')
  403   FORMAT(1X,I8,2F20.10)  
  
        IF(MOD(ITER,5000).EQ.0) THEN
        CALL PRINT_RESULT
        NUM=ITER/NSTEP+NSTART
        FILENAME1(1:3)='SFC'
        FILENAME1(4:4)=CHAR(48+NUM/1000)
        NUM=NUM-NUM/1000*1000
        FILENAME1(5:5)=CHAR(48+NUM/100)
        NUM=NUM-NUM/100*100
        FILENAME1(6:6)=CHAR(48+NUM/10)
        NUM=NUM-NUM/10*10
        FILENAME1(7:7)=CHAR(48+NUM)
        FILENAME1(8:11)='.DAT'
        OPEN(1,FILE=FILENAME1)
          WRITE(1,*) 'TITLE="SF"'
          WRITE(1,*) 'VARIABLES="X","Y","SF"'
          WRITE(1,*) 'ZONE I=',L1,',J=',M1
          DO J=1,M1  
           DO I=1,L1
             WRITE(1,*) X(I),Y(J),F(I,J,3)
            END DO
          END DO
        CLOSE(1)
        END IF
        
        IF(LSTOP.OR.(ITER.EQ.LAST)) CALL PRINT_RESULT        
      
        IF(LSTOP.OR.(ITER.EQ.LAST)) THEN
          OPEN(1,FILE='VEL_U7.DAT')
            DO J=1,M1
            WRITE(1,*) u(450,J),Y(J)-0.5
          END DO
          CLOSE(1)

          OPEN(1,FILE='VEL_U15.DAT')
            DO J=1,M1
            WRITE(1,*) U(961,J),Y(J)-0.5
          END DO
          CLOSE(1)
          
          OPEN(1,FILE='STREAM_FUNCTION.DAT')
          WRITE(1,*) 'TITLE="SF"'
          WRITE(1,*) 'VARIABLES="X","Y","SF"'
          WRITE(1,*) 'ZONE I=',L1,',J=',M1
          DO J=1,M1  
           DO I=1,L1
             WRITE(1,*) X(I),(Y(J)-0.5),F(I,J,3)
            END DO
          END DO
          CLOSE(1)
          
	    OPEN(1,FILE='USEW.DAT')
            DO I=1,L1
	        DO J=1,M1
	          WRITE(1,*) U(I,J),V(I,J)
	        END DO
	      END DO
	    CLOSE(1)
        END IF          
      RETURN
C----------------------------------------------------------------------
      ENTRY GAMSOR
        RE=800.
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
          IF (SMAX<RSMAX) LSTOP=.TRUE.
        END IF
      RETURN
C----------------------------------------------------------------------
      END