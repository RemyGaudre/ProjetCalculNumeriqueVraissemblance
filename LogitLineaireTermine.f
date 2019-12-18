      PROGRAM DEMO
* MMAX et NMAX sont deux constantes. 
* Ce sont les dimensions maximales des matrices
      INTEGER MMAX, NMAX
      PARAMETER (MMAX = 50)
      PARAMETER (NMAX = 50)
* Déclaration de A
      DOUBLE PRECISION A(MMAX, NMAX),B(MMAX,MMAX),ATA(MMAX,MMAX),
     $Q(MMAX),W(MMAX),ATQ(MMAX),X(MMAX),Y(MMAX, NMAX)
* Les vraies dimensions de A sont dans M et N (dim A = M x N)
      INTEGER M, N
*     
* Début des instructions
*
      WRITE (*,*) "Entrez vecteur X"
* La subroutine DREAD_MPL permet de lire une matrice au format MAPLE
      CALL DREAD_MPL (M, N, X, MMAX)
      WRITE (*,*) "Affichage"
      CALL DPRINT_MPL ('X', M, N, X, MMAX)
      
      WRITE (*,*) "Entrez vecteur Y"
* La subroutine DREAD_MPL permet de lire une matrice au format MAPLE
      CALL DREAD_MPL (M,N,Y, NMAX)
      CALL DPRINT_MPL ('Y', M,N, Y, NMAX)
      CALL LOGIT(M,Y,MMAX,Y,MMAX,30.54D0)
      CALL DPRINT_MPL ('LOGIT Y', M,N, Y, NMAX)
      CALL CONSTRUCT_A (M,X,NMAX,A,MMAX)
      CALL DPRINT_MPL ('A', M,2, A, MMAX)
      
      CALL DGEMM ('transpose','notranspose',2,M,M,1D0,A,MMAX,A,MMAX,0D0
     $,ATA,MMAX)
      CALL DPRINT_MPL ('ATA', 2,2, ATA, MMAX)
      CALL CHOLESKI(2,2,ATA,MMAX,B,MMAX)
      
      CALL DPRINT_MPL ('B', 2,2,B, MMAX)
      
      CALL DGEMM ('transpose','notranspose',2,M,M,1D0,A,MMAX,Y,MMAX,0D0
     $,ATQ,MMAX)
      CALL DPRINT_MPL ('ATQ', 2,1,ATQ, MMAX)
      
      
      CALL DESCENTE(2,B,NMAX,ATQ,W)
      
      CALL MONTEE (2,B,NMAX,W,X)
      
      CALL DPRINT_MPL ('La solution trouvee par la methode de CHOLESKI 
     $du système Ax = b est ',2,
     $1,X,NMAX)
     
        
*Creation de la matrice A de l'equation matricielle lineaire pour la resolution du logit
      
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
      
           
        
        
