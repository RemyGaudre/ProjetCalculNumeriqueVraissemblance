      PROGRAM DEMO
* MMAX et NMAX sont deux constantes. 
* Ce sont les dimensions maximales des matrices
      INTEGER MMAX, NMAX
      PARAMETER (MMAX = 50)
      PARAMETER (NMAX = 10)
* Déclaration de A
      DOUBLE PRECISION A (MMAX, NMAX)
* Les vraies dimensions de A sont dans M et N (dim A = M x N)
      INTEGER M, N
*
* Début des instructions
*
      WRITE (*,*) "Entrez une matrice"
* La subroutine DREAD_MPL permet de lire une matrice au format MAPLE
      CALL DREAD_MPL (M, N, A, MMAX)
* La subroutine DPRINT_MPL permet d'afficher une matrice
      WRITE (*,*) "Affichage"
      CALL DPRINT_MPL ('mat', M, N, A, MMAX)
      CALL METTRE_A_ZERO (M, N, A, MMAX)
      WRITE (*,*) "Nouvel affichage"
      CALL DPRINT_MPL ('mat', M, N, A, MMAX)
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

