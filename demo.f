      PROGRAM DEMO
* MMAX et NMAX sont deux constantes. 
* Ce sont les dimensions maximales des matrices
      INTEGER MMAX, NMAX,I
      PARAMETER (MMAX = 50)
      PARAMETER (NMAX = 50)
* Déclaration de A
      DOUBLE PRECISION A(MMAX, NMAX),B(MMAX,MMAX),ATA(MMAX,MMAX),
     $W(MMAX),ATQ(MMAX),X(MMAX),Y(MMAX, NMAX)
     $,V(MMAX),XR(MMAX),LY(MMAX)
* Les vraies dimensions de A sont dans M et N (dim A = M x N)
      INTEGER M, N
*     
* Début des instructions
*
      WRITE (*,*) "Entrez vecteur X"
* La subroutine DREAD_MPL permet de lire une matrice au format MAPLE
      CALL DREAD_MPL (M, N, X, MMAX)
      WRITE (*,*) "Affichage"
      
      WRITE (*,*) "Entrez vecteur Y"
* La subroutine DREAD_MPL permet de lire une matrice au format MAPLE
      CALL DREAD_MPL (M,N,Y, NMAX)
      CALL DPRINT_MPL ('Y',M,1,Y,NMAX)
      CALL LOGIT(M,Y,MMAX,LY,MMAX,30.54D0)
      
      CALL CONSTRUCT_A (M,X,NMAX,A,MMAX)
      
      CALL DGEMM ('transpose','notranspose',2,M,M,1D0,A,MMAX,A,MMAX,0D0
     $,ATA,MMAX)
      CALL CHOLESKI(2,2,ATA,MMAX,B,MMAX)
      
      
      CALL DGEMM ('transpose','notranspose',2,M,M,1D0,A,MMAX,LY,MMAX,0D0
     $,ATQ,MMAX)
      
      
      CALL DESCENTE(2,B,NMAX,ATQ,W)
      
      CALL MONTEE (2,B,NMAX,W,XR)
      
      V = XR
      V(3) = 30.54D0
      
      CALL DPRINT_MPL ('X',M,1,X,NMAX)
      
      
      DO I=1,2
        WRITE (*,*) I
        CALL ITERATION(M,3,V,X,Y,MMAX)
        CALL DPRINT_MPL ('V1 ',3,1,V,NMAX)
      END DO
      
      
      END PROGRAM

      SUBROUTINE METTRE_A_ZERO (M, N, A, LDA)
* Paramètres formels
      INTEGER M, N, LDA
* Le paramètre LDA (Leading Dimension de A) est la dimension "informatique"
* de A. C'est le nombre de lignes du tableau A, tel qu'il a été déclaré.
      DOUBLE PRECISION A (LDA,*)
* Variables locales
      INTEGER I, J
      DO J = 1,N
         DO I = 1,M
            A(I,J) = 0D0
         END DO
      END DO
      END SUBROUTINE

      
      
      SUBROUTINE CHOLESKI (M,N,A,LDA,B,LDB)
        
        INTEGER M,N,LDA,LDB
        DOUBLE PRECISION A(LDA,*),B(LDB,*),SOM
        
        INTEGER I,J,K
        
        B(1,1)=SQRT(A(1,1))
        
        DO J = 1,N
           DO I=2,M
              IF (I.EQ.J) THEN
                 SOM = 0
                 DO K=1,I-1
                    SOM = SOM + B(I,K)**2
                 END DO
              B(J,J)=SQRT(A(J,J)-SOM)
              
              ELSE IF(J.GT.I) THEN
                 B(I,J)=0
              ELSE
                 SOM = 0
                 DO K = 1,I-1
                    SOM = SOM + B(I,K)*B(J,K)
                 END DO
                 B(I,J)=(A(I,J)-SOM)/B(J,J)
              END IF
           END DO
        END DO
              SOM = 0
              DO K=1,N-1
                 SOM = SOM + B(N,K)**2
              END DO

        B(N,N)=SQRT(A(N,N)-SOM)
        
      END SUBROUTINE
      
      SUBROUTINE DESCENTE (N,L,LDL,B,W)
        
        INTEGER N,LDL
        DOUBLE PRECISION L(LDL,*),B(LDL),W(LDL),C
        
        INTEGER I,J
        W(1)=B(1)/L(1,1)
        
        DO I = 2,N
           C = B(I)
           DO J=1,I-1
              C = C - L(I,J)*W(J)
           END DO
           W(I)=C/L(I,I)
        END DO
        
      END SUBROUTINE
      
      SUBROUTINE MONTEE (N,U,LDU,W,X)
        
        INTEGER N,LDU
        DOUBLE PRECISION U(LDU,*),W(LDU),X(LDU),C
        
        INTEGER I,J
        X(N)=W(N)/U(N,N)
        
        DO I = 2,N
           C = W(N-I+1)
           DO J=1,I-1
              C = C - U(N-J+1,N-I+1)*X(N-J+1)
           END DO
           X(N-I+1)=C/U(N-I+1,N-I+1)
        END DO
        
      END SUBROUTINE

      
      
      
      SUBROUTINE LOGIT (N,Y,LDY,R,LDR,KAPPA)
      
        INTEGER N,LDY,LDR
        DOUBLE PRECISION Y(LDY),R(LDR),KAPPA
        
        INTEGER I
        
        DO I = 1,N
            R(I)=LOG((Y(I)/KAPPA)/(1D0-(Y(I)/KAPPA)))
        END DO
        
      END SUBROUTINE
        
    
      SUBROUTINE CONSTRUCT_A (N,X,LDX,A,LDA)
      
        INTEGER N,LDX,LDA
        DOUBLE PRECISION X(LDX),A(LDA,*)
        
        INTEGER I
        
        DO I=1,N
           A(I,1)=X(I)
           A(I,2)=-1D0
        END DO
        
      END SUBROUTINE
      
      SUBROUTINE SIGMOIDE(ALPHA,KAPPA,RHO,X,S)
      
        DOUBLE PRECISION ALPHA, KAPPA, RHO, X, S
        
        S = (KAPPA/(1+EXP(ALPHA-RHO*X)))
      
      END SUBROUTINE
      
      
      SUBROUTINE RESIDU(N,ALPHA,KAPPA,RHO,R,LDR,X,LDX,Y,LDY)
      
        INTEGER N,LDR,LDX,LDY,I
        DOUBLE PRECISION R(LDR),ALPHA, KAPPA, RHO, X(LDX),Y(LDY), S
         
        DO I=1,N
            CALL SIGMOIDE(ALPHA,KAPPA,RHO,X(I),S)
            R(I)= S-Y(I)
        END DO
      
      END SUBROUTINE
           
      SUBROUTINE JACOBIENNE(N,ALPHA,KAPPA,RHO,J,LDJ,X,LDX)
      
        INTEGER N,LDJ,LDX
        DOUBLE PRECISION J(LDJ,*),ALPHA, KAPPA, RHO, X(LDX)
        
        INTEGER I
        
        DO I= 1,N
            J(I,1)= 1 / (1+EXP(ALPHA-RHO*X(I)))
            J(I,2)= -KAPPA *((EXP(ALPHA-RHO*X(I)))
     $/(1+EXP(ALPHA-RHO*X(I)))**2)
            J(I,3)= KAPPA*X(I)*((EXP(ALPHA-RHO*X(I)))
     $/(1+EXP(ALPHA-RHO*X(I)))**2)
        END DO
      
      END SUBROUTINE
        
      
      
      SUBROUTINE ITERATION(M,N,V,X,Y,MMAX)
      
        INTEGER M,N,MMAX
        DOUBLE PRECISION V(MMAX),X(MMAX),Y(MMAX),XR(MMAX),J(MMAX,MMAX),
     $R(MMAX,MMAX),P(MMAX,MMAX),ALPHA,KAPPA,RHO,ATA(MMAX,MMAX),
     $JTR(MMAX,MMAX),W(MMAX)
      
        ALPHA = V(2)
        RHO = V(1)
        KAPPA = V(3)
        WRITE(*,*) "ALPHA"
        WRITE(*,*) ALPHA
        WRITE(*,*) "RHO"
        WRITE(*,*) RHO
        WRITE(*,*) "KAPPA"
        WRITE(*,*) KAPPA
            
        CALL JACOBIENNE(M,ALPHA,KAPPA,RHO,J,MMAX,X,MMAX)
        
        CALL DPRINT_MPL ('J ',M,N,J,MMAX)
        
        CALL RESIDU(M,ALPHA,KAPPA,RHO,R,MMAX,X,MMAX,Y,MMAX)
        
        CALL DPRINT_MPL ('R ',M,1,R,MMAX)
        
        CALL DGEMM ('transpose','notranspose',N,M,M,1D0,J,MMAX,J,
     $MMAX,0D0,ATA,MMAX)
    
        CALL DPRINT_MPL ('JTJ ',N,N,ATA,MMAX)
        CALL CHOLESKI(N,N,ATA,MMAX,P,MMAX)
        CALL DPRINT_MPL ('P ',N,N,P,MMAX)
        CALL DGEMM ('transpose','notranspose',M,1,N,-1D0,J,MMAX,R,
     $MMAX,0D0,JTR,MMAX)
        CALL DPRINT_MPL ('-JTR ',N,1,JTR,MMAX)
        
        
        CALL DESCENTE(N,P,MMAX,JTR,W)
        
        CALL MONTEE (N,P,MMAX,W,XR)
        
        
        V(1) = V(1) + XR(1)
        V(2) = V(2) + XR(2)
        V(3) = V(3) + XR(3)
        
      
      END SUBROUTINE
