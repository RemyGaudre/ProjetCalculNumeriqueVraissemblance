* F. Boulier. Janvier 2012. Modifié en juillet 2013. 
* Étendu aux types autres que double précision en 2016.
*
* Utilitaires permettant les E/S de matrices au format du paquetage
* MAPLE LinearAlgebra.
*
* Les identificateurs contiennent tous un '_', pour les distinguer
* des identificateurs des BLAS et des sous-programmes LAPACK.
*
* Les types autorisés sont : C, D, I, S

* SUBROUTINE CPRINT_MPL (IDENT, M, N, A, LDA)
*
* Impression de la matrice A suivant la syntaxe de MAPLE LinearAlgebra
* Les matrices 1x1 sont vues comme des nombres
* Les matrices dont une seule dimension vaut 1 sont vues comme des vecteurs
*
* IDENT   - CHARACTER array
*           entrée. Nom qu'on souhaite donner à la matrice à l'affichage
*
* M       - INTEGER
*           entrée. Nombre de lignes de A
*
* N       - INTEGER
*           entrée. Nombre de colonnes de A
* 
* A       - COMPLEX array of DIMENSION (LDA,N)
*           entrée. Mathématiquement, A est une matrice M x N.
*           La matrice à afficher.
*
* LDA     - INTEGER
*           entrée. Nombre de lignes « informatique » de la matrice qui 
*           contient A. Vaut M sauf si A est une sous-matrice d'une matrice.

      SUBROUTINE CPRINT_MPL (IDENT, M, N, A, LDA)
      IMPLICIT NONE
      CHARACTER*(*) IDENT
      INTEGER M, N, LDA
      COMPLEX A(LDA,*)
*
      INTEGER R, C
* ident :=
      WRITE (*,'(1X,A,'' :='')') IDENT
      IF (M .EQ. 1 .AND. N .EQ. 1) THEN
* matrices 1x1
*        WRITE (*,'(G11.5,SP,G11.5,"*I")') A(1,1)
         CALL PUT_COMPLEX (' ', A(1,1), .TRUE.)
         WRITE (*,*) ''
      ELSE IF (N .EQ. 1) THEN
* vecteurs colonne
*        WRITE (*,'(1X,A,G11.5,SP,G11.5,"*I",$)') '<', A(1,1)
         CALL PUT_COMPLEX ('<', A(1,1), .TRUE.)
         DO R = 2,M
*           WRITE (*,'(A,G11.5,SP,G11.5,"*I",$)') ',', A(R,1)
            CALL PUT_COMPLEX (',', A(R,1), .FALSE.)
         END DO
         WRITE (*,*) '>;'
      ELSE IF (M .EQ. 1) THEN
* vecteurs ligne
*        WRITE (*,'(1X,A,G11.5,SP,G11.5,"*I",$)') '<', A(1,1)
         CALL PUT_COMPLEX ('<', A(1,1), .TRUE.)
         DO C = 2,N
*           WRITE (*,'(A,G11.5,SP,G11.5,"*I",$)') '|', A(1,C)
            CALL PUT_COMPLEX ('|', A(1,C), .FALSE.)
         END DO
         WRITE (*,*) '>;'
      ELSE
* cas général
         DO R = 1,M
            IF (R .EQ. 1) THEN
*               WRITE (*,'(1X,A,G11.5,SP,G11.5,"*I",$)'), '<<', A(1,1)
                CALL PUT_COMPLEX ('<<', A(1,1), .TRUE.)
            ELSE
*               WRITE (*,'(1X,A,G11.5,SP,G11.5,"*I",$)'), ' <', A(R,1)
                CALL PUT_COMPLEX (' <', A(R,1), .TRUE.)
            END IF
            DO C = 2,N
*               WRITE (*,'(A,G11.5,SP,G11.5,"*I",$)') '|', A(R,C)
                CALL PUT_COMPLEX ('|', A(R,C), .FALSE.)
            END DO
            IF (R .EQ. M) THEN
                WRITE (*,*) '>>;'
            ELSE
                WRITE (*,*) '>,'
            END IF
         END DO
      END IF
      END SUBROUTINE

* SUBROUTINE PUT_COMPLEX (S, Z, FIRST)
*
* Sous-programme de CPRINT_MPL. Impression d'un nombre complexe.
*
* S     - CHARACTER array
*         entrée. Chaîne de caractères à imprimer devant le complexe
*
* Z     - COMPLEX
*         entrée. Le nombre complexe à imprimer
*
* FIRST - LOGICAL
*         entrée. Indique si l'impression se fait en première colonne.

      SUBROUTINE PUT_COMPLEX (S, Z, FIRST)
      CHARACTER*(*) S
      COMPLEX Z
      LOGICAL FIRST
*
      IF (IMAGPART(Z) .EQ. 0.) THEN
         IF (FIRST) THEN
            WRITE (*,'(1X,A,G11.5,11X,"  ",$)') S, REALPART(Z)
         ELSE
            WRITE (*,'(A,G11.5,11X,"  ",$)') S, REALPART(Z)
         END IF
      ELSE IF (REALPART(Z) .EQ. 0.) THEN
         IF (FIRST) THEN
            WRITE (*,'(1X,A,11X,G11.5,"*I",$)') S, IMAGPART(Z)
         ELSE
            WRITE (*,'(A,11X,G11.5,"*I",$)') S, IMAGPART(Z)
         END IF
      ELSE
         IF (FIRST) THEN
            WRITE (*,'(1X,A,G11.5,SP,G11.5,"*I",$)') S, Z
         ELSE
            WRITE (*,'(A,G11.5,SP,G11.5,"*I",$)') S, Z
         END IF
      END IF
      END SUBROUTINE

* SUBROUTINE DPRINT_MPL (IDENT, M, N, A, LDA)
*
* Impression de la matrice A suivant la syntaxe de MAPLE LinearAlgebra
* Les matrices 1x1 sont vues comme des nombres
* Les matrices dont une seule dimension vaut 1 sont vues comme des vecteurs
*
* IDENT   - CHARACTER array
*           entrée. Nom qu'on souhaite donner à la matrice à l'affichage
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

      SUBROUTINE DPRINT_MPL (IDENT, M, N, A, LDA)
      IMPLICIT NONE
      CHARACTER*(*) IDENT
      INTEGER M, N, LDA
      DOUBLE PRECISION A(LDA,*)
*
      INTEGER R, C
* ident :=
      WRITE (*,'(1X,A,'' :='')') IDENT
      IF (M .EQ. 1 .AND. N .EQ. 1) THEN
* matrices 1x1
         WRITE (*,*) A(1,1)
      ELSE IF (N .EQ. 1) THEN
* vecteurs colonne
         WRITE (*,'(1X,A,G16.10,$)') '<', A(1,1)
         DO R = 2,M
            WRITE (*,'(A,G16.10,$)') ',', A(R,1)
         END DO
         WRITE (*,*) '>;'
      ELSE IF (M .EQ. 1) THEN
* vecteurs ligne
         WRITE (*,'(1X,A,G16.10,$)') '<', A(1,1)
         DO C = 2,N
            WRITE (*,'(A,G16.10,$)') '|', A(1,C)
         END DO
         WRITE (*,*) '>;'
      ELSE
* cas général
         DO R = 1,M
            IF (R .EQ. 1) THEN
                WRITE (*,'(1X,A,G16.10,$)') '<<', A(1,1)
            ELSE
                WRITE (*,'(1X,A,G16.10,$)') ' <', A(R,1)
            END IF
            DO C = 2,N
                WRITE (*,'(A,G16.10,$)') '|', A(R,C)
            END DO
            IF (R .EQ. M) THEN
                WRITE (*,*) '>>;'
            ELSE
                WRITE (*,*) '>,'
            END IF
         END DO
      END IF
      END SUBROUTINE

* Pour des raisons de compatibilité avec les vieilles versions.

      SUBROUTINE PRINT_MPL (IDENT, M, N, A, LDA)
      IMPLICIT NONE
      CHARACTER*(*) IDENT
      INTEGER M, N, LDA
      DOUBLE PRECISION A(LDA,*)
      CALL DPRINT_MPL (IDENT, M, N, A, LDA)
      END SUBROUTINE

* SUBROUTINE IPRINT_MPL (IDENT, M, N, A, LDA)
*
* Impression de la matrice A suivant la syntaxe de MAPLE LinearAlgebra
* Les matrices 1x1 sont vues comme des nombres
* Les matrices dont une seule dimension vaut 1 sont vues comme des vecteurs
*
* IDENT   - CHARACTER array
*           entrée. Nom qu'on souhaite donner à la matrice à l'affichage
*
* M       - INTEGER
*           entrée. Nombre de lignes de A
*
* N       - INTEGER
*           entrée. Nombre de colonnes de A
* 
* A       - INTEGER array of DIMENSION (LDA,N)
*           entrée. Mathématiquement, A est une matrice M x N.
*           La matrice à afficher.
*
* LDA     - INTEGER
*           entrée. Nombre de lignes « informatique » de la matrice qui 
*           contient A. Vaut M sauf si A est une sous-matrice d'une matrice.

      SUBROUTINE IPRINT_MPL (IDENT, M, N, A, LDA)
      IMPLICIT NONE
      CHARACTER*(*) IDENT
      INTEGER M, N, LDA
      INTEGER A(LDA,*)
*
      INTEGER R, C
* ident :=
      WRITE (*,'(1X,A,'' :='')') IDENT
      IF (M .EQ. 1 .AND. N .EQ. 1) THEN
* matrices 1x1
         WRITE (*,*) A(1,1)
      ELSE IF (N .EQ. 1) THEN
* vecteurs colonne
         WRITE (*,'(1X,A,I12,$)') '<', A(1,1)
         DO R = 2,M
            WRITE (*,'(A,I12,$)') ',', A(R,1)
         END DO
         WRITE (*,*) '>;'
      ELSE IF (M .EQ. 1) THEN
* vecteurs ligne
         WRITE (*,'(1X,A,I12,$)') '<', A(1,1)
         DO C = 2,N
            WRITE (*,'(A,I12,$)') '|', A(1,C)
         END DO
         WRITE (*,*) '>;'
      ELSE
* cas général
         DO R = 1,M
            IF (R .EQ. 1) THEN
                WRITE (*,'(1X,A,I12,$)') '<<', A(1,1)
            ELSE
                WRITE (*,'(1X,A,I12,$)') ' <', A(R,1)
            END IF
            DO C = 2,N
                WRITE (*,'(A,I12,$)') '|', A(R,C)
            END DO
            IF (R .EQ. M) THEN
                WRITE (*,*) '>>;'
            ELSE
                WRITE (*,*) '>,'
            END IF
         END DO
      END IF
      END SUBROUTINE

* SUBROUTINE SPRINT_MPL (IDENT, M, N, A, LDA)
*
* Impression de la matrice A suivant la syntaxe de MAPLE LinearAlgebra
* Les matrices 1x1 sont vues comme des nombres
* Les matrices dont une seule dimension vaut 1 sont vues comme des vecteurs
*
* IDENT   - CHARACTER array
*           entrée. Nom qu'on souhaite donner à la matrice à l'affichage
*
* M       - INTEGER
*           entrée. Nombre de lignes de A
*
* N       - INTEGER
*           entrée. Nombre de colonnes de A
* 
* A       - REAL array of DIMENSION (LDA,N)
*           entrée. Mathématiquement, A est une matrice M x N.
*           La matrice à afficher.
*
* LDA     - INTEGER
*           entrée. Nombre de lignes « informatique » de la matrice qui 
*           contient A. Vaut M sauf si A est une sous-matrice d'une matrice.

      SUBROUTINE SPRINT_MPL (IDENT, M, N, A, LDA)
      IMPLICIT NONE
      CHARACTER*(*) IDENT
      INTEGER M, N, LDA
      REAL A(LDA,*)
*
      INTEGER R, C
* ident :=
      WRITE (*,'(1X,A,'' :='')') IDENT
      IF (M .EQ. 1 .AND. N .EQ. 1) THEN
* matrices 1x1
         WRITE (*,*) A(1,1)
      ELSE IF (N .EQ. 1) THEN
* vecteurs colonne
         WRITE (*,'(1X,A,G11.5,$)') '<', A(1,1)
         DO R = 2,M
            WRITE (*,'(A,G11.5,$)') ',', A(R,1)
         END DO
         WRITE (*,*) '>;'
      ELSE IF (M .EQ. 1) THEN
* vecteurs ligne
         WRITE (*,'(1X,A,G11.5,$)') '<', A(1,1)
         DO C = 2,N
            WRITE (*,'(A,G11.5,$)') '|', A(1,C)
         END DO
         WRITE (*,*) '>;'
      ELSE
* cas général
         DO R = 1,M
            IF (R .EQ. 1) THEN
                WRITE (*,'(1X,A,G11.5,$)') '<<', A(1,1)
            ELSE
                WRITE (*,'(1X,A,G11.5,$)') ' <', A(R,1)
            END IF
            DO C = 2,N
                WRITE (*,'(A,G11.5,$)') '|', A(R,C)
            END DO
            IF (R .EQ. M) THEN
                WRITE (*,*) '>>;'
            ELSE
                WRITE (*,*) '>,'
            END IF
         END DO
      END IF
      END SUBROUTINE

***********************************************************************
* PARSER
***********************************************************************

* SUBROUTINE CREAD_MPL (M, N, A, LDA)
*
* Lit une matrice ou un vecteur au format du paquetage LinearAlgebra de MAPLE
* Note : si on connaît à l'appel l'une des dimensions de la matrice ou du 
* vecteur, il est possible de passer une constante comme paramètre effectif.
* Note : les matrices dont une dimension vaut 1 sont assimilées à des vecteurs.
*
* Les lignes commençant par '#' sont ignorées
*
* M       - INTEGER
*           entrée/sortie. Reçoit le nombre de lignes.
*
* N       - INTEGER
*           entrée/sortie. Reçoit le nombre de colonnes.
*
* A       - COMPLEX array of DIMENSION (LDA,*)
*           sortie. Reçoit la matrice lue.
*
* LDA     - INTEGER
*           entrée. Nombre de lignes « informatique » de la matrice qui 
*           contient A.

      SUBROUTINE CREAD_MPL (M, N, A, LDA)
      INTEGER M, N, LDA
      COMPLEX A(LDA,*)
*
      CHARACTER*65536 BUFFER
      INTEGER LEVEL, OLDK, I, J, K, L, NEXT_TOKEN
      LOGICAL FIRST, MATRIX
*
* Première boucle : on copie le flux d'entrée dans une chaîne : BUFFER(1:L)
* On cherche la plus petite chaîne bien parenthésée encadrée par "<" et ">"
*
      FIRST = .TRUE.
* FIRST = .TRUE. si on n'a pas encore rencontré le premier "<"
      LEVEL = 0
* LEVEL = niveau d'imbrication des "<" et ">"
      K = 0
      L = 0
* L = longueur de BUFFER
* K = l'indice dans BUFFER auquel l'analyse lexicale est parvenue
      DO WHILE (FIRST .OR. LEVEL .GE. 1)
         K = NEXT_TOKEN (K+1, L, BUFFER)
         IF (K .EQ. 0) THEN
*
* Plus de tokens : lecture d'une nouvelle ligne
*
            K = L
* Saut des lignes commençant par #
            BUFFER(L+1:L+1) = "#"
            DO WHILE (BUFFER(L+1:L+1) .EQ. "#")
               READ (*,FMT='(A)') BUFFER(L+1:)
            END DO
* Calcul de la « vraie » longueur
            L = LEN(BUFFER)
            DO WHILE (BUFFER (L:L) .EQ. " ")
               L = L-1
            END DO
* Insertion d'espaces dont l'absence gêne le parseur de flottants
            I = K
            DO WHILE (I .LT. L)
               IF (I .GT. 1 .AND. BUFFER (I-1:I-1) .NE. " " .AND. 
     $            (BUFFER (I:I) .EQ. "|" .OR. 
     $             BUFFER (I:I) .EQ. ">")) THEN
                  DO J = L,I,-1
                     BUFFER (J+1:J+1) = BUFFER (J:J)
                  END DO
                  BUFFER (I:I) = " "
                  L = L+1
               END IF
               I = I+1
            END DO
*
* Fin du pré-traitement d'une nouvelle ligne
*
         ELSE IF (BUFFER (K:K) .EQ. "<") THEN
            FIRST = .FALSE.
            LEVEL = LEVEL+1
         ELSE IF (BUFFER (K:K) .EQ. ">") THEN
            LEVEL = LEVEL-1
         END IF
      END DO
*
* Seconde boucle : remplissage de A
*
* La matrice est stockée sous forme texte dans BUFFER (1:L)
* K = l'indice dans BUFFER auquel l'analyse lexicale est parvenue
*
      K = NEXT_TOKEN (1, L, BUFFER)
      IF (BUFFER (K:K) .NE. "<") 
     $   STOP "READ_MPL: token initial '<' attendu"
      J = NEXT_TOKEN (K+1, L, BUFFER)
      MATRIX = BUFFER (J:J) .EQ. "<"
* MATRIX = .TRUE. si on lit une matrice, .FALSE. si on lit un vecteur
      LEVEL = 1
* LEVEL = niveau d'imbrication des "<" et ">"
      J = 1
      I = 1
* I et J indices de ligne et de colonne dans A
      FIRST = .TRUE.
* FIRST = .TRUE. si on n'a pas encore fait un tour de la boucle extérieure
* Utile dans le cas de vecteurs, pour lesquels on ne fera qu'un tour 
      DO WHILE (FIRST .OR. (MATRIX .AND. LEVEL .GE. 1))
         IF ((.NOT. FIRST) .OR. MATRIX) THEN
            K = NEXT_TOKEN (K+1, L, BUFFER)
         END IF
         FIRST = .FALSE.
         IF (BUFFER (K:K) .EQ. ">") THEN
            LEVEL = LEVEL-1
         ELSE IF (BUFFER (K:K) .EQ. "<") THEN
            LEVEL = LEVEL+1
            CALL GET_COMPLEX (BUFFER(K+1:), A(I,J))
            DO WHILE (LEVEL .EQ. 2)
               K = NEXT_TOKEN (K+1, L, BUFFER)
               IF (BUFFER (K:K) .EQ. ",") THEN
                  I = I+1
                  CALL GET_COMPLEX (BUFFER(K+1:), A(I,J))
               ELSE IF (BUFFER (K:K) .EQ. "|") THEN
                  J = J+1
                  CALL GET_COMPLEX (BUFFER(K+1:), A(I,J))
               ELSE IF (BUFFER (K:K) .EQ. ">") THEN
                  LEVEL = LEVEL-1
               ELSE
                  STOP 'READ_MPL: erreur de syntaxe'
               END IF
            END DO
            IF (MATRIX) THEN
               OLDK = K
               K = NEXT_TOKEN (K+1, L, BUFFER)
               IF (BUFFER (K:K) .EQ. "|") THEN
                  I = 1
                  J = J+1
               ELSE IF (BUFFER (K:K) .EQ. ",") THEN
                  I = I+1
                  J = 1
               ELSE
                  K = OLDK
               END IF
            END IF
         END IF
      END DO
* Permet de passer des constantes comme paramètres effectifs
      IF (M .NE. I) M = I
      IF (N .NE. J) N = J
      END 

* SUBROUTINE GET_COMPLEX (T, Z)
*
* Sous-programme de CREAD_MPL. 
* Lit dans la chaîne T un nombre complexe écrit suivant une syntaxe
* naturelle et l'affecte à Z.
*
* Exemples: I, -I, 3*I+4 -7.2E4+5E2*I, 345, 0, -1
*
* T     - CHARACTER array
*         entrée. Chaîne de caractères terminée par un caractère
*         différent de "01234567890.DE*I"
*
* Z     - COMPLEXE
*         sortie. Reçoit le nombre complexe lu

      SUBROUTINE GET_COMPLEX (T, Z)
      CHARACTER*(*) T
      COMPLEX Z
*
      INTEGER I, STX, STY
      REAL X, Y
*
      I = 1
      CALL GET_REAL (T, I, X, STX)
      IF (STX .EQ. 0) THEN
         STOP "GET_COMPLEX: nombre complexe attendu"
      END IF
      CALL GET_REAL (T, I, Y, STY)
      IF (STX .EQ. 1) THEN
         IF (STY .EQ. 0) THEN
            Z = CMPLX(X, 0.)
         ELSE IF (STY .EQ. 1) THEN
            STOP "GET_COMPLEX: partie imaginaire attendue"
         ELSE
            Z = CMPLX(X, Y)
         END IF
      ELSE
         IF (STY .EQ. 0) THEN
            Z = CMPLX(0., X)
         ELSE IF (STY .EQ. 1) THEN
            Z = CMPLX(Y, X)
         ELSE
            STOP "GET_COMPLEX: partie reelle attendue"
         END IF
      END IF
      END SUBROUTINE

* SUBROUTINE GET_REAL (T, I, X, ST)
*
* Sous-programme de GET_COMPLEX. Lit dans la chaîne T, à partir de
* l'indice i, la partie réelle ou la partie imaginaire d'un nombre
* complexe et l'affecte à X.
*
* T     - CHARACTER array
*         entrée. Chaîne de caractères
*
* I     - INTEGER
*         entrée/sortie. En entrée, l'indice à partir duquel lire le
*         nombre. En sortie, l'indice du premier caractère qui suit
*         le nombre lu.
*
* X     - REAL
*         sortie. La partie du réelle du nombre lu.
*
* ST    - INTEGER
*         sortie. Code de retour
*         0 : aucun nombre n'a été lu
*         1 : un réel a été lu
*         2 : un nombre imaginaire a t lu

      SUBROUTINE GET_REAL (T, I, X, ST)
      CHARACTER*(*) T
      INTEGER I, ST
      REAL X
*
      CHARACTER*128 U
      LOGICAL LOOP, IMAG, STARTED
      INTEGER J
*
      ST = 1
* La chaîne U reçoit le nombre à lire. J est un indice dans U.
      J = 1
      DO WHILE (T(I:I) .EQ. " ")
         I = I+1
      END DO
      IF (T(I:I) .EQ. "+" .OR. T(I:I) .EQ. "-") THEN
         U(J:J) = T(I:I)
         J = J+1
         I = I+1
      END IF
* STARTED indique si des chiffres ont commencé à être lus.
* IMAG indique si "I" a été lu.
      LOOP = .TRUE.
      STARTED = .FALSE.
      IMAG = .FALSE.
      DO WHILE (LOOP)
         DO WHILE (T(I:I) .EQ. " ")
            I = I+1
         END DO
         IF (T(I:I) .EQ. "D" .OR. T(I:I) .EQ. "E") THEN
            STARTED = .TRUE.
            U(J:J) = T(I:I)
            J = J+1
            I = I+1
            DO WHILE (T(I:I) .EQ. " ")
               I = I+1
            END DO
            IF (T(I:I) .EQ. "+" .OR. T(I:I) .EQ. "-") THEN
               U(J:J) = T(I:I)
               J = J+1
               I = I+1
            END IF
         ELSE IF (T(I:I) .EQ. "I") THEN
            IF (STARTED .OR. IMAG) THEN
               STOP 'GET_COMPLEX: partie imaginaire mal formulee'
            END IF
            I = I+1
            DO WHILE (T(I:I) .EQ. " ")
               I = I+1
            END DO
            IF (T(I:I) .EQ. "*") THEN
               I = I+1
            ELSE
               U(J:J) = "1"
               J = J+1
            END IF
            IMAG = .TRUE.
            ST = 2
         ELSE IF (T(I:I) .EQ. "*") THEN
            IF (IMAG .OR. .NOT. STARTED) THEN
               STOP 'GET_COMPLEX: partie imaginaire mal formulee'
            END IF
            I = I+1
            DO WHILE (T(I:I) .EQ. " ")
               I = I+1
            END DO
            IF (T(I:I) .EQ. "I") THEN
               I = I+1
            ELSE
               STOP 'GET_COMPLEX: partie imaginaire mal formulee'
            END IF
            IMAG = .TRUE.
            ST = 2
         ELSE IF (INDEX ("0123456789.", T(I:I)) .NE. 0) THEN
            STARTED = .TRUE.
            U(J:J) = T(I:I)
            J = J+1
            I = I+1
         ELSE
            LOOP = .FALSE.
         END IF
      END DO
      IF (J .EQ. 1) THEN
         ST = 0
      ELSE
         READ (U(1:J-1),*) X
      END IF
      END SUBROUTINE

* INTEGER FUNCTION NEXT_TOKEN (K, L, BUFFER)
*
* Recherche du prochain token la sous-chaîne de caractères BUFFER(K:L).
* Retourne l'indice du prochain token trouvé dans l'intervalle K:L, 0 
* si pas trouvé. 
*
* K       - INTEGER
*           entrée. Indice de début de la sous-chaîne
*
* L       - INTEGER
*           entrée. Indice de fin de la sous-chaîne
*
* BUFFER  - STRING
*           entrée. La chaîne contenant la sous-chaîne

      INTEGER FUNCTION NEXT_TOKEN (K, L, BUFFER)
      INTEGER K, L
      CHARACTER*(*) BUFFER
*
      CHARACTER*(*) TOKENS
      PARAMETER (TOKENS = "<,>|")
      INTEGER I, P, R
      R = 0
      DO I = 1,LEN(TOKENS)
         P = INDEX (BUFFER (K:L), TOKENS (I:I))
         IF (P .NE. 0 .AND. (R .EQ. 0 .OR. P .LT. R)) R = P
      END DO
      IF (R .EQ. 0) THEN
         NEXT_TOKEN = R
      ELSE
         NEXT_TOKEN = R+K-1
      END IF
      END

* SUBROUTINE DREAD_MPL (M, N, A, LDA)
*
* Lit une matrice ou un vecteur au format du paquetage LinearAlgebra de MAPLE
* Note : si on connaît à l'appel l'une des dimensions de la matrice ou du 
* vecteur, il est possible de passer une constante comme paramètre effectif.
* Note : les matrices dont une dimension vaut 1 sont assimilées à des vecteurs.
*
* Les lignes commençant par '#' sont ignorées
*
* M       - INTEGER
*           entrée/sortie. Reçoit le nombre de lignes.
*
* N       - INTEGER
*           entrée/sortie. Reçoit le nombre de colonnes.
*
* A       - DOUBLE PRECISION array of DIMENSION (LDA,*)
*           sortie. Reçoit la matrice lue.
*
* LDA     - INTEGER
*           entrée. Nombre de lignes « informatique » de la matrice qui 
*           contient A.

      SUBROUTINE DREAD_MPL (M, N, A, LDA)
      INTEGER M, N, LDA
      DOUBLE PRECISION A(LDA,*)
*
      CHARACTER*65536 BUFFER
      INTEGER LEVEL, OLDK, I, J, K, L, NEXT_TOKEN
      LOGICAL FIRST, MATRIX
*
* Première boucle : on copie le flux d'entrée dans une chaîne : BUFFER(1:L)
* On cherche la plus petite chaîne bien parenthésée encadrée par "<" et ">"
*
      FIRST = .TRUE.
* FIRST = .TRUE. si on n'a pas encore rencontré le premier "<"
      LEVEL = 0
* LEVEL = niveau d'imbrication des "<" et ">"
      K = 0
      L = 0
* L = longueur de BUFFER
* K = l'indice dans BUFFER auquel l'analyse lexicale est parvenue
      DO WHILE (FIRST .OR. LEVEL .GE. 1)
         K = NEXT_TOKEN (K+1, L, BUFFER)
         IF (K .EQ. 0) THEN
*
* Plus de tokens : lecture d'une nouvelle ligne
*
            K = L
* Saut des lignes commençant par #
            BUFFER(L+1:L+1) = "#"
            DO WHILE (BUFFER(L+1:L+1) .EQ. "#")
               READ (*,FMT='(A)') BUFFER(L+1:)
            END DO
* Calcul de la « vraie » longueur
            L = LEN(BUFFER)
            DO WHILE (BUFFER (L:L) .EQ. " ")
               L = L-1
            END DO
* Insertion d'espaces dont l'absence gêne le parseur de flottants
            I = K
            DO WHILE (I .LT. L)
               IF (I .GT. 1 .AND. BUFFER (I-1:I-1) .NE. " " .AND. 
     $            (BUFFER (I:I) .EQ. "|" .OR. 
     $             BUFFER (I:I) .EQ. ">")) THEN
                  DO J = L,I,-1
                     BUFFER (J+1:J+1) = BUFFER (J:J)
                  END DO
                  BUFFER (I:I) = " "
                  L = L+1
               END IF
               I = I+1
            END DO
*
* Fin du pré-traitement d'une nouvelle ligne
*
         ELSE IF (BUFFER (K:K) .EQ. "<") THEN
            FIRST = .FALSE.
            LEVEL = LEVEL+1
         ELSE IF (BUFFER (K:K) .EQ. ">") THEN
            LEVEL = LEVEL-1
         END IF
      END DO
*
* Seconde boucle : remplissage de A
*
* La matrice est stockée sous forme texte dans BUFFER (1:L)
* K = l'indice dans BUFFER auquel l'analyse lexicale est parvenue
*
      K = NEXT_TOKEN (1, L, BUFFER)
      IF (BUFFER (K:K) .NE. "<") 
     $   STOP "READ_MPL: token initial '<' attendu"
      J = NEXT_TOKEN (K+1, L, BUFFER)
      MATRIX = BUFFER (J:J) .EQ. "<"
* MATRIX = .TRUE. si on lit une matrice, .FALSE. si on lit un vecteur
      LEVEL = 1
* LEVEL = niveau d'imbrication des "<" et ">"
      J = 1
      I = 1
* I et J indices de ligne et de colonne dans A
      FIRST = .TRUE.
* FIRST = .TRUE. si on n'a pas encore fait un tour de la boucle extérieure
* Utile dans le cas de vecteurs, pour lesquels on ne fera qu'un tour 
      DO WHILE (FIRST .OR. (MATRIX .AND. LEVEL .GE. 1))
         IF ((.NOT. FIRST) .OR. MATRIX) THEN
            K = NEXT_TOKEN (K+1, L, BUFFER)
         END IF
         FIRST = .FALSE.
         IF (BUFFER (K:K) .EQ. ">") THEN
            LEVEL = LEVEL-1
         ELSE IF (BUFFER (K:K) .EQ. "<") THEN
            LEVEL = LEVEL+1
            READ (BUFFER(K+1:),*) A(I,J)
            DO WHILE (LEVEL .EQ. 2)
               K = NEXT_TOKEN (K+1, L, BUFFER)
               IF (BUFFER (K:K) .EQ. ",") THEN
                  I = I+1
                  READ (BUFFER(K+1:),*) A(I,J)
               ELSE IF (BUFFER (K:K) .EQ. "|") THEN
                  J = J+1
                  READ (BUFFER(K+1:),*) A(I,J)
               ELSE IF (BUFFER (K:K) .EQ. ">") THEN
                  LEVEL = LEVEL-1
               ELSE
                  STOP 'READ_MPL: erreur de syntaxe'
               END IF
            END DO
            IF (MATRIX) THEN
               OLDK = K
               K = NEXT_TOKEN (K+1, L, BUFFER)
               IF (BUFFER (K:K) .EQ. "|") THEN
                  I = 1
                  J = J+1
               ELSE IF (BUFFER (K:K) .EQ. ",") THEN
                  I = I+1
                  J = 1
               ELSE
                  K = OLDK
               END IF
            END IF
         END IF
      END DO
* Permet de passer des constantes comme paramètres effectifs
      IF (M .NE. I) M = I
      IF (N .NE. J) N = J
      END 

* Pour des raisons de compatibilité avec les vieilles versions

      SUBROUTINE READ_MPL (M, N, A, LDA)
      INTEGER M, N, LDA
      DOUBLE PRECISION A(LDA,*)
      CALL DREAD_MPL (M, N, A, LDA)
      END SUBROUTINE

* SUBROUTINE SREAD_MPL (M, N, A, LDA)
*
* Lit une matrice ou un vecteur au format du paquetage LinearAlgebra de MAPLE
* Note : si on connaît à l'appel l'une des dimensions de la matrice ou du 
* vecteur, il est possible de passer une constante comme paramètre effectif.
* Note : les matrices dont une dimension vaut 1 sont assimilées à des vecteurs.
*
* Les lignes commençant par '#' sont ignorées
*
* M       - INTEGER
*           entrée/sortie. Reçoit le nombre de lignes.
*
* N       - INTEGER
*           entrée/sortie. Reçoit le nombre de colonnes.
*
* A       - REAL array of DIMENSION (LDA,*)
*           sortie. Reçoit la matrice lue.
*
* LDA     - INTEGER
*           entrée. Nombre de lignes « informatique » de la matrice qui 
*           contient A.

      SUBROUTINE SREAD_MPL (M, N, A, LDA)
      INTEGER M, N, LDA
      REAL A(LDA,*)
*
      CHARACTER*65536 BUFFER
      INTEGER LEVEL, OLDK, I, J, K, L, NEXT_TOKEN
      LOGICAL FIRST, MATRIX
*
* Première boucle : on copie le flux d'entrée dans une chaîne : BUFFER(1:L)
* On cherche la plus petite chaîne bien parenthésée encadrée par "<" et ">"
*
      FIRST = .TRUE.
* FIRST = .TRUE. si on n'a pas encore rencontré le premier "<"
      LEVEL = 0
* LEVEL = niveau d'imbrication des "<" et ">"
      K = 0
      L = 0
* L = longueur de BUFFER
* K = l'indice dans BUFFER auquel l'analyse lexicale est parvenue
      DO WHILE (FIRST .OR. LEVEL .GE. 1)
         K = NEXT_TOKEN (K+1, L, BUFFER)
         IF (K .EQ. 0) THEN
*
* Plus de tokens : lecture d'une nouvelle ligne
*
            K = L
* Saut des lignes commençant par #
            BUFFER(L+1:L+1) = "#"
            DO WHILE (BUFFER(L+1:L+1) .EQ. "#")
               READ (*,FMT='(A)') BUFFER(L+1:)
            END DO
* Calcul de la « vraie » longueur
            L = LEN(BUFFER)
            DO WHILE (BUFFER (L:L) .EQ. " ")
               L = L-1
            END DO
* Insertion d'espaces dont l'absence gêne le parseur de flottants
            I = K
            DO WHILE (I .LT. L)
               IF (I .GT. 1 .AND. BUFFER (I-1:I-1) .NE. " " .AND. 
     $            (BUFFER (I:I) .EQ. "|" .OR. 
     $             BUFFER (I:I) .EQ. ">")) THEN
                  DO J = L,I,-1
                     BUFFER (J+1:J+1) = BUFFER (J:J)
                  END DO
                  BUFFER (I:I) = " "
                  L = L+1
               END IF
               I = I+1
            END DO
*
* Fin du pré-traitement d'une nouvelle ligne
*
         ELSE IF (BUFFER (K:K) .EQ. "<") THEN
            FIRST = .FALSE.
            LEVEL = LEVEL+1
         ELSE IF (BUFFER (K:K) .EQ. ">") THEN
            LEVEL = LEVEL-1
         END IF
      END DO
*
* Seconde boucle : remplissage de A
*
* La matrice est stockée sous forme texte dans BUFFER (1:L)
* K = l'indice dans BUFFER auquel l'analyse lexicale est parvenue
*
      K = NEXT_TOKEN (1, L, BUFFER)
      IF (BUFFER (K:K) .NE. "<") 
     $   STOP "READ_MPL: token initial '<' attendu"
      J = NEXT_TOKEN (K+1, L, BUFFER)
      MATRIX = BUFFER (J:J) .EQ. "<"
* MATRIX = .TRUE. si on lit une matrice, .FALSE. si on lit un vecteur
      LEVEL = 1
* LEVEL = niveau d'imbrication des "<" et ">"
      J = 1
      I = 1
* I et J indices de ligne et de colonne dans A
      FIRST = .TRUE.
* FIRST = .TRUE. si on n'a pas encore fait un tour de la boucle extérieure
* Utile dans le cas de vecteurs, pour lesquels on ne fera qu'un tour 
      DO WHILE (FIRST .OR. (MATRIX .AND. LEVEL .GE. 1))
         IF ((.NOT. FIRST) .OR. MATRIX) THEN
            K = NEXT_TOKEN (K+1, L, BUFFER)
         END IF
         FIRST = .FALSE.
         IF (BUFFER (K:K) .EQ. ">") THEN
            LEVEL = LEVEL-1
         ELSE IF (BUFFER (K:K) .EQ. "<") THEN
            LEVEL = LEVEL+1
            READ (BUFFER(K+1:),*) A(I,J)
            DO WHILE (LEVEL .EQ. 2)
               K = NEXT_TOKEN (K+1, L, BUFFER)
               IF (BUFFER (K:K) .EQ. ",") THEN
                  I = I+1
                  READ (BUFFER(K+1:),*) A(I,J)
               ELSE IF (BUFFER (K:K) .EQ. "|") THEN
                  J = J+1
                  READ (BUFFER(K+1:),*) A(I,J)
               ELSE IF (BUFFER (K:K) .EQ. ">") THEN
                  LEVEL = LEVEL-1
               ELSE
                  STOP 'READ_MPL: erreur de syntaxe'
               END IF
            END DO
            IF (MATRIX) THEN
               OLDK = K
               K = NEXT_TOKEN (K+1, L, BUFFER)
               IF (BUFFER (K:K) .EQ. "|") THEN
                  I = 1
                  J = J+1
               ELSE IF (BUFFER (K:K) .EQ. ",") THEN
                  I = I+1
                  J = 1
               ELSE
                  K = OLDK
               END IF
            END IF
         END IF
      END DO
* Permet de passer des constantes comme paramètres effectifs
      IF (M .NE. I) M = I
      IF (N .NE. J) N = J
      END 

* SUBROUTINE IREAD_MPL (M, N, A, LDA)
*
* Lit une matrice ou un vecteur au format du paquetage LinearAlgebra de MAPLE
* Note : si on connaît à l'appel l'une des dimensions de la matrice ou du 
* vecteur, il est possible de passer une constante comme paramètre effectif.
* Note : les matrices dont une dimension vaut 1 sont assimilées à des vecteurs.
*
* Les lignes commençant par '#' sont ignorées
*
* M       - INTEGER
*           entrée/sortie. Reçoit le nombre de lignes.
*
* N       - INTEGER
*           entrée/sortie. Reçoit le nombre de colonnes.
*
* A       - INTEGER array of DIMENSION (LDA,*)
*           sortie. Reçoit la matrice lue.
*
* LDA     - INTEGER
*           entrée. Nombre de lignes « informatique » de la matrice qui 
*           contient A.

      SUBROUTINE IREAD_MPL (M, N, A, LDA)
      INTEGER M, N, LDA
      INTEGER A(LDA,*)
*
      CHARACTER*65536 BUFFER
      INTEGER LEVEL, OLDK, I, J, K, L, NEXT_TOKEN
      LOGICAL FIRST, MATRIX
*
* Première boucle : on copie le flux d'entrée dans une chaîne : BUFFER(1:L)
* On cherche la plus petite chaîne bien parenthésée encadrée par "<" et ">"
*
      FIRST = .TRUE.
* FIRST = .TRUE. si on n'a pas encore rencontré le premier "<"
      LEVEL = 0
* LEVEL = niveau d'imbrication des "<" et ">"
      K = 0
      L = 0
* L = longueur de BUFFER
* K = l'indice dans BUFFER auquel l'analyse lexicale est parvenue
      DO WHILE (FIRST .OR. LEVEL .GE. 1)
         K = NEXT_TOKEN (K+1, L, BUFFER)
         IF (K .EQ. 0) THEN
*
* Plus de tokens : lecture d'une nouvelle ligne
*
            K = L
* Saut des lignes commençant par #
            BUFFER(L+1:L+1) = "#"
            DO WHILE (BUFFER(L+1:L+1) .EQ. "#")
               READ (*,FMT='(A)') BUFFER(L+1:)
            END DO
* Calcul de la « vraie » longueur
            L = LEN(BUFFER)
            DO WHILE (BUFFER (L:L) .EQ. " ")
               L = L-1
            END DO
* Insertion d'espaces dont l'absence gêne le parseur de flottants
            I = K
            DO WHILE (I .LT. L)
               IF (I .GT. 1 .AND. BUFFER (I-1:I-1) .NE. " " .AND. 
     $            (BUFFER (I:I) .EQ. "|" .OR. 
     $             BUFFER (I:I) .EQ. ">")) THEN
                  DO J = L,I,-1
                     BUFFER (J+1:J+1) = BUFFER (J:J)
                  END DO
                  BUFFER (I:I) = " "
                  L = L+1
               END IF
               I = I+1
            END DO
*
* Fin du pré-traitement d'une nouvelle ligne
*
         ELSE IF (BUFFER (K:K) .EQ. "<") THEN
            FIRST = .FALSE.
            LEVEL = LEVEL+1
         ELSE IF (BUFFER (K:K) .EQ. ">") THEN
            LEVEL = LEVEL-1
         END IF
      END DO
*
* Seconde boucle : remplissage de A
*
* La matrice est stockée sous forme texte dans BUFFER (1:L)
* K = l'indice dans BUFFER auquel l'analyse lexicale est parvenue
*
      K = NEXT_TOKEN (1, L, BUFFER)
      IF (BUFFER (K:K) .NE. "<") 
     $   STOP "READ_MPL: token initial '<' attendu"
      J = NEXT_TOKEN (K+1, L, BUFFER)
      MATRIX = BUFFER (J:J) .EQ. "<"
* MATRIX = .TRUE. si on lit une matrice, .FALSE. si on lit un vecteur
      LEVEL = 1
* LEVEL = niveau d'imbrication des "<" et ">"
      J = 1
      I = 1
* I et J indices de ligne et de colonne dans A
      FIRST = .TRUE.
* FIRST = .TRUE. si on n'a pas encore fait un tour de la boucle extérieure
* Utile dans le cas de vecteurs, pour lesquels on ne fera qu'un tour 
      DO WHILE (FIRST .OR. (MATRIX .AND. LEVEL .GE. 1))
         IF ((.NOT. FIRST) .OR. MATRIX) THEN
            K = NEXT_TOKEN (K+1, L, BUFFER)
         END IF
         FIRST = .FALSE.
         IF (BUFFER (K:K) .EQ. ">") THEN
            LEVEL = LEVEL-1
         ELSE IF (BUFFER (K:K) .EQ. "<") THEN
            LEVEL = LEVEL+1
            READ (BUFFER(K+1:),*) A(I,J)
            DO WHILE (LEVEL .EQ. 2)
               K = NEXT_TOKEN (K+1, L, BUFFER)
               IF (BUFFER (K:K) .EQ. ",") THEN
                  I = I+1
                  READ (BUFFER(K+1:),*) A(I,J)
               ELSE IF (BUFFER (K:K) .EQ. "|") THEN
                  J = J+1
                  READ (BUFFER(K+1:),*) A(I,J)
               ELSE IF (BUFFER (K:K) .EQ. ">") THEN
                  LEVEL = LEVEL-1
               ELSE
                  STOP 'READ_MPL: erreur de syntaxe'
               END IF
            END DO
            IF (MATRIX) THEN
               OLDK = K
               K = NEXT_TOKEN (K+1, L, BUFFER)
               IF (BUFFER (K:K) .EQ. "|") THEN
                  I = 1
                  J = J+1
               ELSE IF (BUFFER (K:K) .EQ. ",") THEN
                  I = I+1
                  J = 1
               ELSE
                  K = OLDK
               END IF
            END IF
         END IF
      END DO
* Permet de passer des constantes comme paramètres effectifs
      IF (M .NE. I) M = I
      IF (N .NE. J) N = J
      END 

