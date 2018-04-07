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
CALL Upwind_D(N,Velo)
IF (Velo==2.5) THEN
CALL Mixed_D(N,Velo)
END IF
CALL Analytic_Sol(N,Velo)
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

!Upwind Format
SUBROUTINE Upwind_D(N,Velo)
REAL::Velo
INTEGER::N
REAL::Len=1,Dens=1,Gama=0.1,D,F,Delta_x,Zero
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
Zero=0

DO i=2,N,1
!上风格式的系数
A_e(i)=D+MAX(-F,Zero)
A_w(i)=D+MAX(F,Zero)
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

OPEN(UNIT=1,FILE='Upwind_D.txt')
DO i=1,N+1,1
WRITE(1,*) X(i)
END DO
RETURN
END SUBROUTINE 

!Mixed Format
SUBROUTINE Mixed_D(N,Velo)
REAL::Velo
INTEGER::N
REAL::Len=1,Dens=1,Gama=0.1,D,F,Delta_x,Zero
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
Zero=0

DO i=2,N,1
!混合格式的系数
A_e(i)=MAX(-F,D-F/2,Zero)
A_w(i)=MAX(F,D+F/2,Zero)
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

OPEN(UNIT=1,FILE='Mixed_D.txt')
DO i=1,N+1,1
WRITE(1,*) X(i)
END DO
RETURN
END SUBROUTINE 

!Analysic Solution
SUBROUTINE Analytic_Sol(N,Velo)
REAL::Velo
INTEGER::N
REAL::Len=1,Dens=1,Gama=0.1,Delta_x
REAL,ALLOCATABLE::X(:)
Delta_x=Len/N
ALLOCATE(X(N+2))

X(1)=1
X(N+2)=0;
DO i=2,N+1,1
X(i)=-1*(EXP(Dens*Velo*(i-1.5)*Delta_x/Gama)-1)/(EXP(Dens*Velo*Len/Gama)-1)+1
END DO
OPEN(UNIT=1,FILE='Analysic_Sol.txt')
DO i=1,N+2,1
WRITE(1,*) X(i)
END DO

RETURN
END SUBROUTINE 



!Center Format
!mayugao edit
SUBROUTINE Center_Format(N,Velo)
REAL::Velo
INTEGER::N
REAL::Len=1,Dens=1,Gama=0.1,D,F,Delta_x
REAL,ALLOCATABLE::A_e(:),A_p(:),A_w(:),X(:),A_c(:),d_(:),c_(:)


ALLOCATE(A_e(N))
ALLOCATE(A_p(N))
ALLOCATE(A_w(N))
ALLOCATE(A_c(N))
ALLOCATE(X(N))

ALLOCATE(d_(N))
ALLOCATE(c_(N))

Delta_x=Len/N
D=Gama/Delta_x
F=Dens*Velo 

!以下是马誉高推导的边界条件
A_e(1)=D-F/2
A_w(1)=0
A_p(1)=D+3*F/2
A_c(1)=(F+2*D)

A_e(N)=0
A_w(N)=D+F/2
A_p(N)=3*D+F/2
A_c(N)=0

!以下是冯夷宁组的边界条件
!A_e(1)=(D-F/2)
!A_w(1)=0
!A_p(1)=2*D
!A_c(1)=(D+F/2)

!A_e(N)=0
!A_w(N)=(D+F/2)
!A_p(N)=2*D
!A_c(N)=0

DO i=2,N-1,1
!中心差分的系数
A_e(i)=(D-F/2)
A_w(i)=(D+F/2)
A_p(i)=2*D
A_c(i)=0
END DO


call TDMA(A_p,A_e,A_w,A_c,N,X)
write(*,*) A_p,A_e,A_w,A_c,N,X

OPEN(UNIT=1,FILE='Center_D.txt')
DO i=1,N,1
WRITE(1,*) X(i)
END DO
RETURN
END SUBROUTINE 

! TDMA算法
! TRI-DIAGONAL SOLUTION ALGORITHM
! NOTE: ARRAYS B AND D ARE DESTROYED IN THE PROCESS.
! FOR THE SOLUTION OF:     -C*T(I-1) + A*T(I) - B*T(I+1) = D
! 改编自柯道友老师的讲解版本
! 可运行
SUBROUTINE TDMA(A_p,A_e,A_w,A_c,N,X)
    INTEGER::N
    REAL :: A_e(N),A_p(N),A_w(N),X(N),A_c(N)


	A_e(1)=(A_e(1))/(A_p(1))
	A_c(1)=A_c(1)/(A_p(1))
	DO I = 2 , N
	    A_e(I)=A_e(I)/((A_p(I))-A_w(I)*A_e(I-1))
	    A_c(I)=(A_c(I)+A_w(I)*A_c(I-1))/((A_p(I))-A_w(I)*A_e(I-1))
    ENDDO
	X(N)=A_c(N)
	DO I = N-1 , 1 , -1                          !step size = -1
	    X(I)=A_e(I)*X(I+1)+A_c(I)
    ENDDO
END SUBROUTINE


! A_p*T(I) = A_w*T(I-1)+A_e*T(I+1)+A_c
! TDMA求解方程
! 改编自冯夷宁版本
! 可运行，无差别
SUBROUTINE TDMA2(A_p,A_e,A_w,A_c,N,X)
    INTEGER::N
    REAL :: A_e(N),A_p(N),A_w(N),X(N),A_c(N),c_(N),d_(N)

    c_(1)=-A_e(1)/A_p(1)
    d_(1)=A_c(1)/A_p(1)
    DO I=2,(N-1)
        c_(I)=(-A_e(I))/(A_p(I)-c_(I-1)*(-A_w(I)))
        d_(I)=(A_c(I)-d_(I-1)*(-A_w(I)))/(A_p(I)-c_(I-1)*(-A_w(I)))
    ENDDO
        d_(N)=(A_c(N)-d_(N-1)*(-A_w(N)))/(A_p(N)-c_(N-1)*(-A_w(N-1)))
        X(N)=d_(N)
    DO I=1,(N-1)
        X(N-I)=d_(N-I)-c_(N-I)*X(N-I+1)
    ENDDO
END SUBROUTINE
