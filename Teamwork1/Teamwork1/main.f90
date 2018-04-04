! -------------------------------------------------------
! Campare Discrete Format:Exponential format
! Group1
! 2018.4.30
! -------------------------------------------------------

PROGRAM Teamwork1
IMPLICIT NONE

REAL::Velo
INTEGER::N
WRITE(*,*) 'INPUT THE NUMBER OF CELLS'
READ(*,*)N
WRITE(*,*)'INPUT THE INLET VELOCITY'
READ(*,*)Velo

CALL Center_D(N,Velo)
CALL Exp_D(N,Velo)
Pause
END PROGRAM Teamwork1             
 
!Center Format
SUBROUTINE Center_D(N,Velo)
REAL::Velo
INTEGER::N
REAL::Len=1,Dens=1,Gama=0.1,D,F,Delta_x
REAL,ALLOCATABLE::A_e(:),A_p(:),A_w(:),X(:),Y(:),A(:),B(:),C(:),L(:),U(:),K(:)

ALLOCATE(A_e(N))
ALLOCATE(A_p(N))
ALLOCATE(A_w(N))
ALLOCATE(X(N+1))
ALLOCATE(Y(N))
ALLOCATE(A(N))
ALLOCATE(B(N))
ALLOCATE(C(N))
ALLOCATE(L(N))
ALLOCATE(U(N))
ALLOCATE(K(N))
Delta_x=Len/N
D=Gama/Delta_x
F=Dens*Velo 

DO i=2,N,1
!中心差分的系数
A_e(i)=D-F/2
A_w(i)=D+F/2
A_p(i)=A_e(i)+A_w(i)

A(i)=-A_w(i)
B(i)=A_p(i)
C(i)=-A_e(i)  
L(i)=B(i)
U(i)=C(i)
END DO

X(1)=1
X(N+1)=0
Y(2)=A_w(2)*X(1)
Y(N)=A_e(N)*X(N+1)

DO i=3,N-1,1
Y(i)=0
END DO

DO i=2,N-1,1
K(i)=A(i+1)/L(i)
L(i+1)=L(i+1)-U(i)*K(i)
Y(i+1)=Y(i+1)-Y(i)*K(i)
END DO

X(N)=Y(N)/L(N)

DO i=N-1,2,-1
X(i)=(Y(i)-U(i)*X(i+1))/L(i)
END DO

OPEN(UNIT=1,FILE='Center_D.txt')
DO i=1,N+1,1
WRITE(1,*) X(i)
END DO
RETURN
END SUBROUTINE 

! Exponential format
SUBROUTINE Exp_D(N,Velo)
REAL::Velo
INTEGER::N
REAL::Len=1,Dens=1,Gama=0.1,D,F,Delta_x
REAL,ALLOCATABLE::A_e(:),A_p(:),A_w(:),X(:),Y(:),A(:),B(:),C(:),L(:),U(:),K(:)

ALLOCATE(A_e(N))
ALLOCATE(A_p(N))
ALLOCATE(A_w(N))
ALLOCATE(X(N+1))
ALLOCATE(Y(N))
ALLOCATE(A(N))
ALLOCATE(B(N))
ALLOCATE(C(N))
ALLOCATE(L(N))
ALLOCATE(U(N))
ALLOCATE(K(N))
Delta_x=Len/N
D=Gama/Delta_x
F=Dens*Velo 

DO i=2,N,1
!指数格式的系数
A_e(i)=F/(EXP(F/D)-1)
A_w(i)=F*EXP(F/D)/(EXP(F/D)-1)
A_p(i)=A_e(i)+A_w(i)

!TDMA求解
A(i)=-A_w(i)
B(i)=A_p(i)
C(i)=-A_e(i)  
L(i)=B(i)
U(i)=C(i)
END DO

X(1)=1
X(N+1)=0
Y(2)=A_w(2)*X(1)
Y(N)=A_e(N)*X(N+1)

DO i=3,N-1,1
Y(i)=0
END DO

DO i=2,N-1,1
K(i)=A(i+1)/L(i)
L(i+1)=L(i+1)-U(i)*K(i)
Y(i+1)=Y(i+1)-Y(i)*K(i)
END DO

X(N)=Y(N)/L(N)

DO i=N-1,2,-1
X(i)=(Y(i)-U(i)*X(i+1))/L(i)
END DO

OPEN(UNIT=1,FILE='Exp_D.txt')
DO i=1,N+1,1
WRITE(1,*) X(i)
END DO
RETURN
END SUBROUTINE 


! 标准TDMA算法，还没有加入到里面，需要确定一下
!TRI-DIAGONAL SOLUTION ALGORITHM
!NOTE ARRAYS B AND D ARE DESTROYED IN THE PROCESS.
!FOR THE SOLUTION OF:     -C*T(I-1) + A*T(I) - B*T(I+1) = D
SUBROUTINE TDMA(A,B,C,D,N,X)
	IMPLICIT REAL*8 (A-H,O-Z)
	REAL*8 A(1),B(1),C(1),D(1),X(1)
	B(1)=B(1)/A(1)
	D(1)=D(1)/A(1)
	DO I = 2 , N
	    B(I)=B(I)/(A(I)-C(I)*B(I-1))
	    D(I)=(D(I)+C(I)*D(I-1))/(A(I)-C(I)*B(I-1))
     ENDDO
	X(N)=D(N)
	DO I = N-1 , 1 , -1                          !step size = -1
	    X(I)=B(I)*X(I+1)+D(I)
     ENDDO
	RETURN
	END

