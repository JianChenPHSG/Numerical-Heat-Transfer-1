      SUBROUTINE MODREX
      !IMPLICIT NONE
      COMMON/modl/E,DE,E_1,E_2,ESUM,DSUM,REV_P,C_TIME
	
      REAL,DIMENSION(7,7)::A_P
	REAL,DIMENSION(7,7)::A_I
      REAL,DIMENSION(7,7)::A_D

	REAL,DIMENSION(21,1)::P_E
	REAL,DIMENSION(21,1)::P_DE
      REAL,DIMENSION(21,1)::P_P
	REAL,DIMENSION(21,1)::P_I
      REAL,DIMENSION(21,1)::P_D

      REAL,DIMENSION(7,1)::P_U
      REAL,DIMENSION(7,1)::P_DU
	INTEGER,DIMENSION(3,1)::B_U
      INTEGER,DIMENSION(3,1)::B_DU
	INTEGER::j=1,k=1,FLAG1,FLAG2
      
	REAL::NUM=0.0,DEN=0.0,KP=0.01,KI=0.04,KD=0.01,DKP=0.0,DKI=0.0,DKD=0.0 
	REAL::NB=-1.0,NM=-0.3,NS=-0.1,NO=0.0,PS=0.1,PM=0.2,PB=1.0
	REAL::X=0.0,Y=0.0,Z=0.0
	
  
	A_P(1,:)=(/PB,PB,PM,PB,PS,ZO,ZO/)
      A_P(2,:)=(/PB,PB,PM,PS,PS,ZO,NS/)
      A_P(3,:)=(/PM,PM,PM,PS,ZO,NS,NS/)
	A_P(4,:)=(/PM,PM,PS,ZO,NS,NM,NM/)
      A_P(5,:)=(/PS,PS,ZO,NS,NS,NM,NM/)
	A_P(6,:)=(/PS,ZO,NS,NM,NM,NM,NB/)
	A_P(7,:)=(/ZO,ZO,NM,NM,NM,NB,NB/)


	A_I(1,:)=(/NB,NB,NM,NM,NS,ZO,ZO/)
      A_I(2,:)=(/NB,NB,NM,NS,NS,ZO,NS/)
      A_I(3,:)=(/NB,NM,NS,NS,ZO,PS,PS/)
	A_I(4,:)=(/NM,NM,NS,ZO,PS,PM,PM/)
      A_I(5,:)=(/NM,NS,ZO,PS,PS,PB,PB/)
	A_I(6,:)=(/ZO,ZO,PS,PS,PM,PB,PB/)
	A_I(7,:)=(/ZO,ZO,PS,PM,PM,PB,PB/)


	A_D(1,:)=(/PS,NS,NB,NB,NB,NM,PS/)
      A_D(2,:)=(/PS,NS,NM,NM,NM,NS,ZO/)
      A_D(3,:)=(/ZO,NS,NS,NM,NS,NS,ZO/)
	A_D(4,:)=(/ZO,NS,NS,NS,NS,NS,ZO/)
      A_D(5,:)=(/ZO,ZO,ZO,ZO,ZO,ZO,ZO/)
	A_D(6,:)=(/PB,NS,PS,PS,PS,PS,PB/)
	A_D(7,:)=(/PB,PM,PM,PM,PS,PS,PB/)

	P_E(:,1)=(/NB,NB,NM,NB,NM,NS,NM,NS,NO,NS,NO,PS,NO,PS,PM,PS,PM,PB,PM,PB,
     &  PB/)

	P_DE(:,1)=(/NB,NB,NM,NB,NM,NS,NM,NS,NO,NS,NO,PS,NO,PS,PM,PS,PM,PB,PM,
     &  PB,PB/)

	P_P(:,1)=(/NB,NB,NM,NB,NM,NS,NM,NS,NO,NS,NO,PS,NO,PS,PM,PS,PM,PB,PM,PB,
     &  PB/)

	P_I(:,1)=(/NB,NB,NM,NB,NM,NS,NM,NS,NO,NS,NO,PS,NO,PS,PM,PS,PM,PB,PM,PB,
     &  PB/)
      P_D(:,1)=(/NB,NB,NM,NB,NM,NS,NM,NS,NO,NS,NO,PS,NO,PS,PM,PS,PM,PB,
     &  PM,PB,PB/)

	j=1
	k=1
      DO i=1,7,1
	   CALL TRIAN(E,P_E((i-1)*3+1,1),P_E((i-1)*3+2,1),P_E((i-1)*3+3,1),
     &  P_U(i,1)) 
	   CALL TRIAN(DE,P_DE((i-1)*3+1,1),P_DE((i-1)*3+2,1),P_DE((i-1)*3+3,1),
     &	  P_DU(i,1)) 
		 IF(P_U(i,1)/=0.0)	THEN
		    B_U(j,1)=i
			j=j+1
		ELSE 
		    FLAG1=i
		  END IF
		 IF(P_DU(i,1)/=0.0) THEN
		    B_DU(k,1)=i
			k=k+1
		 ELSE
		    FLAG2=i
		 END IF
	END DO  
   
      DO i=j,3,1
	   B_U(i,1)=FLAG1
	END DO

	DO i=k,3,1
	   B_DU(i,1)=FLAG2
	END DO

      DO i=1,3,1
	  DO j=1,3,1
		  NUM=NUM+P_U(B_U(i,1),1)*P_DU(B_DU(j,1),1)*A_P(B_U(i,1),B_DU(j,1))
		  DEN=DEN+P_U(B_U(i,1),1)*P_DU(B_DU(j,1),1)
        END DO
	END DO

	DKP=NUM/DEN
	IF(DKP>0.3)  THEN
         DKP=0.3
	ELSE IF(DKP<-0.3) THEN
	   DKP=-0.3
	END IF 
	KP=KP+DKP
	IF(KP<0.0)  THEN
	KP=0.0
	END IF 

	NUM=0.0
	DEN=0.0
      DO i=1,3,1
	  DO j=1,3,1
		  NUM=NUM+P_U(B_U(i,1),1)*P_DU(B_DU(j,1),1)*A_I(B_U(i,1),B_DU(j,1))
		  DEN=DEN+P_U(B_U(i,1),1)*P_DU(B_DU(j,1),1)
        END DO
	END DO
	DKI=NUM/DEN
	IF(DKI>0.9)  THEN
         DKI=0.9
	ELSE IF(DKI<-0.9) THEN
	   DKI=-0.9
	END IF 
	KI=KI+DKI
	IF(KI<0.0)  THEN
	KI=0.0
      END IF

	NUM=0.0
	DEN=0.0
	DO i=1,3,1
	  DO j=1,3,1
		  NUM=NUM+P_U(B_U(i,1),1)*P_DU(B_DU(j,1),1)*A_D(B_U(i,1),B_DU(j,1))
		  DEN=DEN+P_U(B_U(i,1),1)*P_DU(B_DU(j,1),1)

!		  PRINT *,'dassda' 
!	    PRINT *, NUM
!	    PRINT *, DEN
!	    PRINT *, 'DASDSA'

        END DO
	END DO
	DKD=NUM/DEN
	IF(DKD>0.6)  THEN
         DKD=0.6
	ELSE IF(DKD<-0.6) THEN
	   DKD=-0.6
	END IF 
	KD=KD+DKD
	IF(KD<0.0)  THEN
	KD=0.0
      END IF

	X=KP+KI+KD
	Y=-2*KD-KP
	Z=KD
      REV_P=(1+X*E+Y*E_1+Z*E_2)

	E_2=E_1
	E_1=E

!	WRITE(11,*) 'E'
!	WRITE(11,*)  E
!!	WRITE (11,*) 'DE'
!!	WRITE(11,*)  DE
!!	WRITE(11,*)  'E_1'
!!	WRITE(11,*)	E_1
!!	WRITE(11,*)	'E_2'
!!	WRITE(11,*)  E_2
	WRITE(11,*) 'REV_P'
      WRITE(11,*) REV_P

      END SUBROUTINE MODREX


      SUBROUTINE TRIAN(x,a,b,c,u)
	REAL::x,a,b,c,u
	IF((x .GE.a).AND.(x .LE.b)) THEN
		 u=(x-a)/(b-a)
	ELSEIF((x .GT.b).AND.(x .LT.c)) THEN
	   u=(c-x)/(c-b)
      ELSE
	   u=0
	END IF
	END SUBROUTINE TRIAN