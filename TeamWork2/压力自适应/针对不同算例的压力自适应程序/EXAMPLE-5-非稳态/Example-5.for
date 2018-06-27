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
      PARAMETER (NX=52,NY=52,NXY=52)
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
     &  IST,JST,ITER,C_TIME,LAST,TITLE(13),RELAX(13),TIME,DT,XL,YL,
     &  IPREF,JPREF,LSOLVE(10),LPRINT(13),LBLK(10),MODE,NTIMES(10),
     &  RHOCON,PRELAX(NX,NY)
      COMMON/CNTL/LSTOP,LEND
      COMMON ITIME,MTIME
      COMMON/SORC/SMAX,SSUM,RSMAX
      COMMON/COEF/FLOW,DIFF,ACOF
      DIMENSION U(NX,NY),V(NX,NY),PC(NX,NY)
      EQUIVALENCE (F(1,1,1),U(1,1)),(F(1,1,2),V(1,1)),(F(1,1,3),PC(1,1))
      DIMENSION TH(NXY),THU(NXY),THDIF(NXY),THCV(NXY),THCVS(NXY)
      EQUIVALENCE(X,TH),(XU,THU),(XDIF,THDIF),(XCV,THCV),
     &  (XCVS,THCVS),(XL,THL)
***********************************************************************
*                                                                     *
*        Example-5: IMPINGING FLOW WITH A ROTATING PLATE              *
*                                                                     *
***********************************************************************
      DIMENSION WR(NX,NY)
      EQUIVALENCE(F(1,1,4),WR(1,1))
      DATA UIN,OMEGA/1.,100./
      
      ENTRY GRID
        XL=1.E-2
        YL=1.E-2
        L1=52
        M1=52
	  CALL UGRID
      RETURN
C----------------------------------------------------------------------
      ENTRY START
        MODE=2
        R(1)=0.
        
        LAST=3000
        RSMAX=1E-9
	          
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
	  RELAX(1)=0.5
	  RELAX(2)=0.5  
	  RELAX(4)=0.5
	  RELAX(11)=0.8	  
	  TITLE(1)='.VEL-U.'
	  TITLE(2)='.VEL-V.'
	  TITLE(3)='.STRF.'	  
	  TITLE(4)='.VEL-WR.'

        DO I=1,L1
          DO J=1,M1
            U(I,J)=0.
            V(I,J)=0.
            WR(I,J)=0.
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
        DO I=1,L1
          DO J=1,M1
            RHO(I,J)=1.176
          END DO
        END DO
      RETURN
C----------------------------------------------------------------------
      ENTRY BOUND
	  DO J=1,M1
          IF(J.LE.6) THEN
	      U(2,J)=UIN
          ELSE
            U(2,J)=0.
          END IF
          V(1,J)=0.
          WR(1,J)=0.
          
          U(L1,J)=0.
          V(L1,J)=0.
          WR(L1,J)=OMEGA*(R(J)**2.)
	  END DO

	  DO I=1,L1
	    U(I,1)=U(I,2)
	    V(I,2)=0.
	    WR(I,1)=0.
	    
	    IF(I.LE.46) THEN
	      U(I,M1)=0.
	      V(I,M1)=0.
	      WR(I,M1)=0.
	    ELSE
	      U(I,M1)=0.
	      V(I,M1)=V(I,M2)*RMN(M2)/RMN(M1)
	      WR(I,M1)=WR(I,M2)
	    END IF
	  END DO
      RETURN
C----------------------------------------------------------------------
      ENTRY OUTPUT
        IF(ITER.GE.2000) THEN
        FLOWIN=0.
        FLOWOUT=0.
        DO J=1,6
          FLOWIN=FLOWIN+U(2,J)*ARX(J)
        END DO          
        DO I=47,L1
          FLOWOUT=FLOWOUT+V(I,M1)*XCV(I)*R(M1)
        END DO
        FACTOR=FLOWOUT/FLOWIN
        DO I=47,L1
          V(I,M1)=V(I,M1)/FACTOR
        END DO
        END IF
        
        IF((ITER.EQ.1).and.(ITIME.EQ.1)) WRITE(8,400)
*        IF(MOD(ITER,20).EQ.0) THEN
          WRITE(8,403) ITIME,ITER,C_TIME,RELAX(NP),U(25,25),V(25,25),
     &     WR(25,25),SMAX
          WRITE(*,403) ITIME,ITER,C_TIME,RELAX(NP),U(25,25),V(25,25),
     &     WR(25,25),SMAX          
*        END IF

  400   FORMAT('ITIME',13X,'ITER',12X,'C_TIME',12X,'RELAX(NP)',
     &         12X,'U(25,25)',12X,'V(25,25)',12X,'WR(25,25)',12X,'SMAX')
  403   FORMAT(1X,2I8,6F20.10)  
  
*        IF(LSTOP.OR.(ITER.EQ.LAST)) CALL PRINT_RESULT        
      
        IF(LSTOP.OR.(ITER.EQ.LAST)) THEN
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
        END IF          
      RETURN
C----------------------------------------------------------------------
      ENTRY GAMSOR
        
        DO I=1,L1
          DO J=1,M1
            GAM(I,J)=1.862E-5
          END DO
        END DO
        
        IF(NF.EQ.2) THEN
          DO J=3,M2
            DO I=2,L2
              RSWM=FY(J)*WR(I,J)+FYM(J)*WR(I,J-1)
              CON(I,J)=RHO(I,J)*(RSWM**2.)/(RMN(J)**3.)
              AP(I,J)=-2.*GAM(I,J)/(RMN(J)**2.)
            END DO
          END DO
        END IF
        
        IF(NF.EQ.4) THEN
          DO J=2,M2
            DO I=2,L2
              TEMP=2.*GAM(I,J)/(R(J)*YDIF(J))
              CON(I,J)=TEMP*WR(I,J-1)
              AP(I,J)=-TEMP
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