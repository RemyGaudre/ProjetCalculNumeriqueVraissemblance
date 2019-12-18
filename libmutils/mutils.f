* F. Boulier. Janvier 2012. Modifié en juillet 2013.
*
* Utilitaires de manipulation de matrices
* Les identificateurs contiennent tous un '_', pour les distinguer
* des identificateurs des BLAS et des sous-programmes LAPACK.

* SUBROUTINE RAND_M (M, N, A, LDA, ISEED)
*
* Génération de matrice aléatoire.
*
* M       - INTEGER
*           entrée. Nombre de lignes de A
*
* N       - INTEGER
*           entrée. Nombre de colonnes de A
*
* A       - DOUBLE PRECISION array of DIMENSION (LDA,N)
*           sortie. Mathématiquement, A est une matrice M x N.
*           Reçoit des valeurs aléatoires.
*
* LDA     - INTEGER
*           entrée. Nombre de lignes « informatique » de la matrice qui 
*           contient A. Vaut M sauf si A est une sous-matrice d'une matrice.
*
* ISEED   - INTEGER array of DIMENSION (4)
*           entrée/sortie. ISEED(4) doit être impair.
*           Utilisé pour la génération de nombres aléatoires.
*

      SUBROUTINE RAND_M (M, N, A, LDA, ISEED)
      IMPLICIT NONE
      INTEGER M, N, LDA
      DOUBLE PRECISION A(LDA,N)
      INTEGER ISEED (4)
*
      CHARACTER DIST
      CHARACTER SYM
      DOUBLE PRECISION D (MIN(M,N))
*
      INTEGER MODE
      DOUBLE PRECISION DMAX, COND
      CHARACTER RSIGN
*
      CHARACTER GRADE
      DOUBLE PRECISION  DL(M)
      INTEGER MODEL
      DOUBLE PRECISION  CONDL
      DOUBLE PRECISION DR(N)
      INTEGER MODER
      DOUBLE PRECISION  CONDR
*
      CHARACTER PIVTNG
      INTEGER IPIVOT (MAX(M,N))
      DOUBLE PRECISION  SPARSE
      INTEGER KL, KU
      DOUBLE PRECISION  ANORM
      CHARACTER PACK
*      
      INTEGER IWORK (MAX(M,N))
      INTEGER INFO
*           'S' => UNIFORM( -1, 1 ) ( 'S' for symmetric )
      DIST = 'S'
*           If SYM='N', generated matrix is nonsymmetric.
      SYM = 'N'
*           MODE = 6 set D to random numbers from same distribution
*                    as the rest of the matrix.
      MODE = 6
      COND = 0.0D0
      DMAX = 0.0D0
*           'F' => diagonal unchanged
      RSIGN = 'F'
*           'N'  => no grading
      GRADE = 'N'
* The following are meaningless when GRADE = 'N'
      MODEL = 6
      CONDL = 0.0D0
      MODER = 6
      CONDR = 0.0D0
*           'N' or ' ' => none.
      PIVTNG = 'N'
*           Not sparsed
      SPARSE = 0.0D0
*           On entry specifies the lower bandwidth of the  matrix. 
      KL = M
*           On entry specifies the upper bandwidth of the  matrix. 
      KU = N
*           If ANORM is negative no scaling is done. Not modified.
      ANORM = -1.0D0
*           'N' => no packing
      PACK = 'N'
*
      CALL DLATMR (M, N, DIST, ISEED, SYM, D, MODE, COND, DMAX,
     $              RSIGN, GRADE, DL, MODEL, CONDL, DR, MODER,
     $              CONDR, PIVTNG, IPIVOT, KL, KU, SPARSE, ANORM,
     $              PACK, A, LDA, IWORK, INFO)
      IF (INFO .NE. 0) STOP 'Erreur dans DLATMR'
      END

* SUBROUTINE RAND_SYM (M, A, LDA, EIGENV, ISEED)
*
* Génération de matrice symétrique aléatoire avec des valeurs propres fixées.
* Calcule Q . EIGENV . Transpose (Q) avec Q orthogonale aléatoire.
*
* M       - INTEGER 
*           entrée. Nombre de lignes et de colonnes de A
*
* A       - DOUBLE PRECISION array of DIMENSION (LDA,M)
*           sortie. Mathématiquement, A est une matrice carrée M x M.
*           Reçoit des valeurs aléatoires
*
* LDA     - INTEGER
*           entrée. Nombre de lignes « informatique » de la matrice qui 
*           contient A. Vaut M sauf si A est une sous-matrice d'une matrice.
*
* EIGENV  - DOUBLE PRECISION array of DIMENSION (M)
*           entrée. Les valeurs propres désirées.
*
* ISEED   - INTEGER array of DIMENSION (4)
*           entrée/sortie. ISEED(4) doit être impair.
*           Utilisé pour la génération de nombres aléatoires.

      SUBROUTINE RAND_SYM (M, A, LDA, EIGENV, ISEED)
      IMPLICIT NONE
      INTEGER M, LDA
      DOUBLE PRECISION EIGENV(M)
      DOUBLE PRECISION A(LDA,M)
      INTEGER ISEED(4)
*
      DOUBLE PRECISION X(M,M), Y(M,M)
      DOUBLE PRECISION TAU(M)
      DOUBLE PRECISION WORK(M)
      INTEGER C, INFO
* Génération de X aléatoire
      CALL RAND_M (M, M, X, M, ISEED)
* Orthogonalisation de X
      CALL DGEQRF (M, M, X, M, TAU, WORK, M, INFO)
      IF (INFO .NE. 0) STOP 'Erreur dans DGEQRF'
      CALL DORGQR (M, M, M, X, M, TAU, WORK, M, INFO)
      IF (INFO .NE. 0) STOP 'Erreur dans DORGQR'
* Y = copie de X
      CALL DCOPY (M*M, X, 1, Y, 1)
* X = X . EIGEN
      DO C = 1,M
         CALL DSCAL (M, EIGENV(C), X(1,C), 1)
      END DO
* A = X . Y**T
      CALL DGEMM ('N', 'T', M, M, M, 1.0D0, X, M, Y, M, 0.0D0, A, LDA)
      END

* SUBROUTINE RAND_CM (M, A, LDA, KAPPA, ISEED)
*
* Génération de matrice aléatoire de condition approx. égale à KAPPA
*
* L'algorithme consiste à calculer A = U . D . Transpose (V) où :
* U et V sont deux matrices aléatoires orthogonalisées
* D = diag (1, KAPPA**(-1/(M-1)), ..., KAPPA**(M-2)/(M-1), 1)
* Voir [Afternotes on Num. Anal., page 129]
* Voir aussi [Matrix Anal. and Appl. Lin. Alg., page 321]
*
* M       - INTEGER
*           entrée. Nombre de lignes et de colonnes de A
*
* A       - DOUBLE PRECISION array of DIMENSION (LDA,M)
*           sortie. Mathématiquement, A est une matrice carrée M x M.
*           Reçoit des valeurs aléatoires
*
* LDA     - INTEGER
*           entrée. Nombre de lignes « informatique » de la matrice qui 
*           contient A. Vaut M sauf si A est une sous-matrice d'une matrice.
*
* KAPPA   - DOUBLE PRECISION
*           entrée. La condition désirée.
*
* ISEED   - INTEGER array of DIMENSION (4)
*           entrée/sortie. ISEED(4) doit être impair.
*           Utilisé pour la génération de nombres aléatoires.

      SUBROUTINE RAND_CM (M, A, LDA, KAPPA, ISEED)
      IMPLICIT NONE
      INTEGER M, LDA
      DOUBLE PRECISION A(LDA, M)
      DOUBLE PRECISION KAPPA
      INTEGER ISEED(4)
*
      DOUBLE PRECISION U(M,M)
      DOUBLE PRECISION V(M,M)
      DOUBLE PRECISION TAU(M)
      DOUBLE PRECISION WORK(M)
      INTEGER C, INFO
*
      CALL RAND_M (M, M, U, M, ISEED)
* Orthogonalisation de U
      CALL DGEQRF (M, M, U, M, TAU, WORK, M, INFO)
      IF (INFO .NE. 0) STOP 'Erreur dans DGEQRF'
      CALL DORGQR (M, M, M, U, M, TAU, WORK, M, INFO)
      IF (INFO .NE. 0) STOP 'Erreur dans DORGQR'
* Orthogonalisation de V
      CALL DGEQRF (M, M, V, M, TAU, WORK, M, INFO)
      IF (INFO .NE. 0) STOP 'Erreur dans DGEQRF'
      CALL DORGQR (M, M, M, V, M, TAU, WORK, M, INFO)
      IF (INFO .NE. 0) STOP 'Erreur dans DORGQR'
*
      DO C = 1,M
         CALL DSCAL (M, KAPPA**(-DBLE(C-1)/DBLE(M-1)), U(1,C), 1)
      END DO
      CALL DGEMM ('N', 'T', M, M, M, 1.0D0, U, M, V, M, 0.0D0, A, LDA)
      END

* SUBROUTINE PRINT_MMRK (M, N, A, LDA)
*
* Impression de la matrice A(M,N), suivant un format MatrixMarket
*
* M       - INTEGER
*           entrée. Nombre de lignes de A
*
* N       - INTEGER
*           entrée. Nombre de colonnes de A
* 
* A       - DOUBLE PRECISION array of DIMENSION (LDA,N)
*           entrée. Mathématiquement, A est une matrice M x N.
*           La matrice à afficher.
*
* LDA     - INTEGER
*           entrée. Nombre de lignes « informatique » de la matrice qui 
*           contient A. Vaut M sauf si A est une sous-matrice d'une matrice.

      SUBROUTINE PRINT_MMRK (M, N, A, LDA)
      IMPLICIT NONE
      INTEGER M, N, LDA
      DOUBLE PRECISION A (LDA,N)
      INTEGER R, C, NBNZRO
*
      NBNZRO = 0
      DO R=1,M
          DO C=1,N
             IF (A(R,C) .NE. 0.0D0) NBNZRO = NBNZRO + 1
          ENDDO
      ENDDO
      WRITE (*,*) M, N, NBNZRO
      DO R=1,M
          DO C=1,N
             IF (A(R,C) .NE. 0.0D0) WRITE (*,*) R, C, A(R,C)
          ENDDO
      ENDDO
      END

* SUBROUTINE READ_MMRK (M, N, A, LDA)
* 
* Lit une matrice suivant un format MatrixMarket.
* Description du format :
* La première ligne contient trois entiers : 
* - nombre de lignes M, nombre de colonnes N, nombre d'éléments non nuls NBNZ
* Les NBNZ lignes suivantes contiennent deux entiers et un flottant : 
* - indice de ligne entre 1 et M, de colonne entre 1 et N, valeur
* 
* Par exemple, la suite de lignes suivantes
*
* 2 2 3
* 1 1 3.0
* 2 1 -1.5
* 2 2 17
*
* Représente la matrice
*
*                                [3.0      0. ]
*                                [            ]
*                                [-1.5    17.0]
*
* M       - INTEGER
*           sortie. Le nombre de lignes de A
*
* N       - INTEGER
*           sortie. Le nombre de colonnes de A
*
* A       - DOUBLE PRECISION array of DIMENSION (LDA,N) si LDA > 0 ou bien
*           DOUBLE PRECISION array of DIMENSION (M,N)
*           sortie. Mathématiquement, A est une matrice M x N.
*           Reçoit la matrice lue.
*
* LDA     - INTEGER
*           entrée. Si LDA <= 0, ce paramètre est ignoré et M est utilisé 
*           à la place. Les colonnes de A sont alors rangées consécutivement 
*           en mémoire et les programmes qui utilisent A peuvent déclarer 
*           cette matrice A(M,N).
*
*           Si LDA > 0, nombre de lignes « informatique » de la matrice qui 
*           contient A. Doit alors être supérieur ou égal à M.

      SUBROUTINE READ_MMRK (M, N, A, LDA)
      IMPLICIT NONE
      DOUBLE PRECISION A(*)
      INTEGER M, N, LDA
      INTEGER NBNZRO, MM, NN
*
      READ (*,*) MM, NN, NBNZRO
      IF (LDA .LE. 0) THEN
         CALL READ_MMRK2 (MM, NN, A, MM, NBNZRO)
      ELSE IF (LDA .LT. MM) THEN
         STOP 'Erreur dans READ_MMRK'
      ELSE
         CALL READ_MMRK2 (MM, NN, A, LDA, NBNZRO)
      END IF
      IF (M .NE. MM) M = MM
      IF (N .NE. NN) N = NN
      END

* SUBROUTINE READ_MMRK2 (M, N, A, LDA, NBNZRO)
*
* Sous-programme de READ_MMRK. Lit NBNZRO lignes de trois nombres :
* indice de ligne, de colonne et valeur ; et les range dans A.
*
* M       - INTEGER
*           entrée. Nombre de lignes de A
*
* N       - INTEGER
*           entrée. Nombre de colonnes de A
* 
* A       - DOUBLE PRECISION array of DIMENSION (LDA,N)
*           sortie. Mathématiquement, A est une matrice M x N.
*           Reçoit la matrice lue.
*
* LDA     - INTEGER
*           entrée. Nombre de lignes « informatique » de la matrice qui 
*           contient A. Vaut M sauf si A est une sous-matrice d'une matrice.
*
* NBNZRO  - INTEGER
*           entrée. Nombre d'éléments non nuls de A

      SUBROUTINE READ_MMRK2 (M, N, A, LDA, NBNZRO)
      IMPLICIT NONE
      INTEGER M, N, LDA, NBNZRO
      DOUBLE PRECISION A(LDA,N), V
      INTEGER R, C, I
*
      DO C = 1,N
         DO R = 1,M
            A(R,C) = 0.0D0
         END DO
      END DO
      DO I = 1,NBNZRO
         READ (*,*) R, C, V
         IF (R .GT. M .OR. C .GT. N .OR. R .LE. 0 .OR. C .LE. 0) 
     $       STOP 'Erreur dans READ_MMRK2: indices incorrects'
         A(R,C) = V
      END DO
      END

* DOUBLE PRECISION FUNCTION COND1_M (M, A, LDA)
*
* Retourne une estimation du conditionnement de A pour la norme 1
*
* M       - INTEGER
*           entrée. Nombre de lignes et de colonnes de A
*
* A       - DOUBLE PRECISION array of DIMENSION (LDA,M)
*           entrée. Mathématiquement, A est une matrice carrée M x M
*           La matrice dont on cherche la condition
*
* LDA     - INTEGER
*           entrée. Nombre de lignes « informatique » de la matrice qui 
*           contient A. Vaut M sauf si A est une sous-matrice d'une matrice.

      DOUBLE PRECISION FUNCTION COND1_M (M, A, LDA)
      IMPLICIT NONE
      INTEGER M, LDA
      DOUBLE PRECISION A (LDA,M)
*
      DOUBLE PRECISION NRM1_M
      DOUBLE PRECISION ABAR (M,M)
      INTEGER IPIV (M)
      DOUBLE PRECISION WORK (4*M)
      DOUBLE PRECISION ANORM, RCOND
      INTEGER IWORK (M)
      INTEGER C, INFO
* ABAR = A
      DO C = 1,M
         CALL DCOPY (M, A(1,C), 1, ABAR(1,C), 1)
      END DO
* Décomposition LU de ABAR, nécessaire pour ce qui suit
      CALL DGETRF (M, M, ABAR, M, IPIV, INFO)
      IF (INFO .NE. 0) STOP 'Erreur dans DGETRF'
* ANORM = la norme 1 de A
*     CALL PRINT_MPL ('A', M, M, A, LDA)
      ANORM = NRM1_M (M, A, LDA)
*     WRITE (*,*) 'ANORM =', ANORM
* Calcul d'une estimation de l'inverse du conditionnement 
* à partir de la décomposition LU.
      CALL DGECON ('1', M, ABAR, M, ANORM, RCOND, WORK, IWORK, INFO)
      IF (INFO .NE. 0) STOP 'Erreur dans DGECON'
      COND1_M = 1.0D0 / RCOND
      END

* DOUBLE PRECISION FUNCTION NRM1_M (M, A, LDA)
*
* Retourne la norme 1 de A. Sous-fonction de COND1_M
*
* M       - INTEGER
*           entrée. Nombre de lignes et de colonnes de A
*
* A       - DOUBLE PRECISION array of DIMENSION (LDA,M)
*           entrée. Mathématiquement, A est une matrice carrée M x M
*           La matrice dont on cherche la norme 1
*
* LDA     - INTEGER
*           entrée. Nombre de lignes « informatique » de la matrice qui 
*           contient A. Vaut M sauf si A est une sous-matrice d'une matrice.

      DOUBLE PRECISION FUNCTION NRM1_M (M, A, LDA)
      IMPLICIT NONE
      INTEGER M, LDA
      DOUBLE PRECISION A (LDA,M)
*
      DOUBLE PRECISION ANORM, ASUM, DASUM
      INTEGER C
*
      ANORM = 0D0
      DO C = 1,M
         ASUM = DASUM (M, A(1,C), 1)
         IF (ASUM .GT. ANORM) ANORM = ASUM
      END DO
      NRM1_M = ANORM
      END

