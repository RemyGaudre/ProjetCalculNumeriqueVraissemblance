      PROGRAM resolutionlogit
* MMAX et NMAX sont deux constantes. 
* Ce sont les dimensions maximales des matrices
      INTEGER MMAX, NMAX
      PARAMETER (MMAX = 50)
      PARAMETER (NMAX = 10)
* Déclaration de A
      DOUBLE PRECISION A (MMAX, NMAX),B(NMAX,NMAX),ATA(NMAX,NMAX),
     $Q(NMAX),W(NMAX),ATQ(NMAX),X(NMAX)
* Les vraies dimensions de A sont dans M et N (dim A = M x N)
      INTEGER M, N
*     
* Début des instructions
*
      WRITE (*,*) "Entrez une matrice A"
* La subroutine DREAD_MPL permet de lire une matrice au format MAPLE
      CALL DREAD_MPL (M, N, A, MMAX)
      WRITE (*,*) "Affichage"
      CALL DPRINT_MPL ('A', M, N, A, MMAX)
      
      WRITE (*,*) "Entrez une matrice b"
* La subroutine DREAD_MPL permet de lire une matrice au format MAPLE
      CALL DREAD_MPL (M,1,Q, NMAX)
      CALL DPRINT_MPL ('b', M, 1, Q, MMAX)
* La subroutine DPRINT_MPL permet d'afficher une matrice

      CALL DGEMM ('transpose','notranspose',N,M,M,1D0,A,MMAX,A,MMAX,0D0
     $,ATA,NMAX)
      CALL CHOLESKI(N,N,ATA,NMAX,B,NMAX)
      
      
      
      CALL DGEMM ('transpose','notranspose',N,1,M,1D0,A,MMAX,Q,MMAX,0D0
     $,ATQ,NMAX)
      
      CALL DESCENTE(N,B,NMAX,ATQ,W)
      
      
      CALL MONTEE (N,B,NMAX,W,X)
      
      CALL DPRINT_MPL ('La solution trouvee par la methode de CHOLESKI 
     $du système Ax = b est ',N,
     $1,X,NMAX)
      
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

