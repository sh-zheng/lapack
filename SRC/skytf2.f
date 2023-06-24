*> \brief \b SKYTF2 computes the factorization of a real skew-symmetric matrix, using the Bunch partial pivoting method (unblocked algorithm).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download SKYTF2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/skytf2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/skytf2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/skytf2.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE SKYTF2( UPLO, N, A, LDA, IPIV, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       REAL               A( LDA, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SKYTF2 computes the factorization of a real skew-symmetric matrix A using
*> the Bunch block diagonal pivoting method:
*>
*>    A = U*D*(-U)**T  or  A = L*D*(-L)**T
*>
*> where U (or L) is a product of permutation and unit upper (lower)
*> triangular matrices, U**T is the transpose of U, and D is skew-symmetric
*> and block diagonal with 2-by-2 diagonal blocks. But if A is singular,
*> there are some 1-by-1 diagonal blocks, which are all 0.
*>
*> This is the unblocked version of the algorithm, calling Level 2 BLAS.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the upper or lower triangular part of the
*>          skew-symmetric matrix A is stored:
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>          On entry, the skew-symmetric matrix A.  If UPLO = 'U', the
*>          strictly upper triangular part of A contains the upper
*>          triangular part of the matrix A, and the leading N-by-N lower
*>          triangular part of A is not referenced.  If UPLO = 'L', the
*>          strictly lower triangular part of A contains the lower
*>          triangular part of the matrix A, and the leading N-by-N upper
*>          triangular part of A is not referenced.
*>
*>          On exit, the block diagonal matrix D and the multipliers used
*>          to obtain the factor U or L (see below for further details).
*>
*>          The block diagonal matrix D is always 2-by-2 unless A is
*>          singular. In this case the fore (if UPLO = 'U') or the
*>          back (if UPLO = 'L') blocks are all 0.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          Details of the interchanges of D.
*>
*>          The elements of array IPIV are combined in pair,
*>          starting from the last (if UPLO = 'U') element,
*>          or the first (if UPLO = 'L') element, and the second
*>          (if UPLO = 'U') or the first (if UPLO = 'L') element in
*>          the pair always keeps the value 0.  If the order of A
*>          is odd, the first (if UPLO = 'U') or the last (if UPLO = 'L')
*>          element of IPIV is 0, which is the only element not in pair.
*>          
*>          In the following, we only use the first (if UPLO = 'L') or
*>          the second (if UPLO = 'U') element k in the pair value 
*>          to refer to IPIV(k).
*>
*>          IPIV(k) as an INTEGER
*>          = 0: there was no interchange.
*>          > 0: rows and columns k-1 and IPIV(k) were interchanged, if
*>               UPLO = 'U', and rows and columns k+1 and IPIV(k) were
*>               interchanged, if UPLO = 'L'.
*>          < 0: rows and columns k and k-1 were interchanged,
*>               then rows and columns k-1 and -IPIV(k) were interchanged, if
*>               UPLO = 'U', and rows and columns k and k+1 were interchanged,
*>               then rows and columns k+1 and -IPIV(k) were interchanged, if
*>               UPLO = 'L'.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
*>          > 0: if INFO = i, D(i-1,i) (if UPLO = 'U') or D(i+1,i) (if UPLO = 'L')
*>               is exactly zero.  The factorization has been completed,
*>               but the block diagonal matrix D is exactly singular,
*>               so the solution could not be computed.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup realKYcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  If UPLO = 'U', then A = U*D*(-U)**T, where
*>     U = P(n)*U(n)* ... *P(k)U(k)* ...,
*>  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
*>  1 in steps of 2, and D is a block diagonal matrix with 2-by-2
*>  diagonal blocks D(k).  P(k) is a permutation matrix as defined by
*>  IPIV(k), and U(k) is a unit upper triangular matrix, such that if
*>  the diagonal block D(k) is of order 2, then
*>
*>             (   I    v    0   )   k-s
*>     U(k) =  (   0    I    0   )   s
*>             (   0    0    I   )   n-k
*>                k-s   s   n-k
*>
*>  The upper triangle of D(k) overwrites A(k-1,k), and v overwrites
*>  A(1:k-2,k-1:k).
*>
*>  If UPLO = 'L', then A = L*D*(-L)**T, where
*>     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
*>  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
*>  n in steps of 2, and D is a block diagonal matrix with 2-by-2
*>  diagonal blocks D(k).  P(k) is a permutation matrix as defined by
*>  IPIV(k), and L(k) is a unit lower triangular matrix, such that if
*>  the diagonal block D(k) is of order 2, then
*>
*>             (   I    0     0   )  k-1
*>     L(k) =  (   0    I     0   )  s
*>             (   0    v     I   )  n-k-s+1
*>                k-1   s  n-k-s+1
*>
*>  The lower triangle of D(k) overwrites A(k+1,k), and v overwrites
*>  A(k+2:n,k:k+1).
*>
*>  If n is odd or INFO > 0, A is singular. In this situation, D may have
*>  1-by-1 diagonal blocks(the values are all 0).
*> \endverbatim
*
*> \par Contributors:
*  ==================
*>
*> \verbatim
*>
*>  09-29-06 - patch from
*>    Bobby Cheng, MathWorks
*>
*>    Replace l.204 and l.372
*>         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
*>    by
*>         IF( (MAX( ABSAKK, COLMAX ).EQ.ZERO) .OR. SISNAN(ABSAKK) ) THEN
*>
*>  01-01-96 - Based on modifications by
*>    J. Lewis, Boeing Computer Services Company
*>    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA
*>  1-96 - Based on modifications by J. Lewis, Boeing Computer Services
*>         Company
*>
*> \endverbatim
*
*  =====================================================================
      SUBROUTINE SKYTF2( UPLO, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      REAL               EIGHT, SEVTEN
      PARAMETER          ( EIGHT = 8.0E+0, SEVTEN = 17.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, IMAX1, IMAX2, J, JMAX, K, KK, KP, KSTEP
     $                   JP, JJ
      REAL               ABSAKP1K, ALPHA, COLMAX1, COLMAX2, D11, D12,
     $                   D21, D22, R1, ROWMAX, T, WK, WKM1, WKP1
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, SISNAN
      INTEGER            ISAMAX
      EXTERNAL           LSAME, ISAMAX, SISNAN
*     ..
*     .. External Subroutines ..
      EXTERNAL           SSCAL, SSWAP, SSYR, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SKYTF2', -INFO )
         RETURN
      END IF

      IF( UPPER ) THEN
*
*        Factorize A as U*D*(-U)**T using the upper triangle of A
*        K is the main loop index, decreasing from N to 1 in steps
*        of 2
*
         K = N
   10    CONTINUE
*
*        If K <= 1, exit from loop
*
         IF( K.EQ.1 ) THEN
            IF( INFO.EQ.0 )
     $         INFO = K
            KP = 0
            IPIV( K ) = KP
            GO TO 100
         END IF

         IF( K.LT.1 )
     $      GO TO 100
         KSTEP = 2
*
*        Determine rows and columns to be interchanged
*
         ABSAKP1K = ABS( A( K-1, K ) )
*
*        IMAX1 is the row-index of the absolute value largest element in
*        row 1 to K-2, column K.
*        IMAX2 is the row-index of the absolute value largest element in
*        row 1 to K-2 column K-1.
*        COLMAX1 and COLMAX2 are their absolute values.
*
         IMAX1 = ISAMAX( K-2, A( 1, K ), 1 )
         COLMAX1 = ABS( A( IMAX1, K ) )
         IMAX2 = ISAMAX( K-2, A( 1, K-1 ), 1 )
         COLMAX2 = ABS( A( IMAX2, K-1 ) )
*
         IF( MAX(MAX( ABSAKP1K, COLMAX1 ), COLMAX2).EQ.ZERO ) THEN
*
*           Column K and K+1 is zero or underflow: set INFO and continue
*
            IF( INFO.EQ.0 )
     $         INFO = K
            KP = 0
            IPIV( K ) = KP
         ELSE
            IF( ABSAKP1K.GE.MAX( COLMAX1, COLMAX2 ) ) THEN
*
*              No interchange
*
               KP = 0
               IPIV( K ) = KP
            ELSE
               IF( COLMAX1.GE.COLMAX2 ) THEN

*
*                 Absolute value largest element is in column K
*                 Interchange rows and columns K-1 and IMAX1
*                  
                  KP = IMAX1
                  IPIV( K ) = KP

                  CALL SSWAP( K-IMAX1-2, A( IMAX1, IMAX1+1 ), LDA,
     $                     A( IMAX1+1, K-1 ), 1 )

                  CALL SSCAL( K-IMAX1-2, -ONE, A( IMAX1, IMAX1+1 ),
     $                     LDA )

                  CALL SSCAL( K-IMAX1-2, -ONE, A( IMAX1+1, K-1 ),
     $                     1 )
             
                  CALL SSWAP( IMAX1-1, A( 1, IMAX1 ), 1,
     $                     A( 1, K-1 ), 1 ) 

                  A( IMAX1, K-1 ) = -A( IMAX1, K-1 )

*                 
*                 Interchange rows K-1 and IMAX1 in last K-1 columns of A
*     
                  CALL SSWAP( N-K+1, A( K-1, K ), LDA,
     $                        A( IMAX1, K ), LDA )
               ELSE
*
*                 Absolute value largest element is in column K-1
*                 Interchange rows and columns K and K-1, then Interchange K-1 and IMAX2
*              
                  KP = -IMAX2
                  IPIV( K ) = KP

                  CALL SSWAP( K-2, A( 1, K ), 1, A( 1, K-1 ),
     $                     1 )
             
                  A( K-1, K ) = -A( K-1, K )

                  CALL SSWAP( K-IMAX2-2, A( IMAX2, IMAX2+1 ), LDA,
     $                     A( IMAX2+1, K-1 ), 1 )

                  CALL SSCAL( K-IMAX2-2, -ONE, A( IMAX2, IMAX2+1 ),
     $                     LDA )

                  CALL SSCAL( K-IMAX2-2, -ONE, A( IMAX2+1, K-1 ),
     $                     1 )
             
                  CALL SSWAP( IMAX2-1, A( 1, IMAX2 ), 1,
     $                     A( 1, K-1 ), 1 ) 

                  A( IMAX2, K-1 ) = -A( IMAX2, K-1 )

*                 
*                 Interchange rows K and K-1, then K-1 and IMAX2 in last K+1 columns of A
*
                  IF( K.LT.N ) THEN
                     CALL SSWAP( N-K, A( K, K+1 ), LDA, A( K-1, K+1 ),
     $                           LDA )
                  END IF

                  CALL SSWAP( N-K+1, A( K-1, K ), LDA,
     $                        A( IMAX2, K ), LDA )

               END IF
            END IF
*
*           Update the lower triangle of A11 (= A(1:k-2,1:k-2))
*
            D21 = ONE/A( K-1, K )

            DO 20 J = 1, K-2
*
               WK = -A( J, K-1 )*D21
               WKM1 = A( J, K )*D21
*
               DO 30 I = J+1, K-2
                  A( J, I ) = A( J, I ) + A( I, K )*WK +
     $               A( I, K-1 )*WKM1
30             CONTINUE

20          CONTINUE

*
*           Update C*S^-1
*
            DO 80 J = 1, K-2
               T = A( J, K-1 )
               A( J, K-1 ) = A( J, K )*D21
               A( J, K ) = -T*D21
80          CONTINUE
         END IF
*
*        Decrease K and return to the start of the main loop
*
         K = K - KSTEP
         GO TO 10
*
100      CONTINUE
*
*
*        Put U12 in standard form by partially undoing the interchanges
*        of rows in columns 1:k-1 looping backwards from k-1 to 1
*
         J = N - K - 1
110      CONTINUE
*
*           Undo the interchanges (if any) of rows JJ and JP at each
*           step J
*
*           (Here, J is a diagonal index)

            IF( J.GT.1 ) THEN
               JJ = N-J+1
               JP = IPIV( N-J+1 )

               IF( JP.LT.0 ) THEN
                  JP = -JP
*                 (Here, J is a diagonal index)
                  CALL SSWAP( J-1, A( JP, N-J+2 ), LDA,
     $                  A( JJ-1, N-J+2 ), LDA )
                  CALL SSWAP( J-1, A( JJ-1, N-J+2 ), LDA,
     $                  A( JJ, N-J+2 ), LDA )
               ELSEIF( JP.GT.0 ) THEN
                  CALL SSWAP( J-1, A( JP, N-J+2 ), LDA,
     $                  A( JJ-1, N-J+2 ), LDA )
               END IF
               
            END IF
*        (NOTE: Here, J is used to determine row length. Length J
*        of the rows to swap back doesn't include diagonal element)

         J = J - 2
         IF( J.GT.1 )
     $      GO TO 110
     
*
      ELSE
*
*        Factorize A as L*D*(-L)**T using the lower triangle of A
*        K is the main loop index, increasing from 1 to N in steps
*        of 2
*
         K = 1
   40    CONTINUE
*
*        If K >= N, exit from loop
*
         IF( K.EQ.N ) THEN
            IF( INFO.EQ.0 )
     $         INFO = K
            KP = 0
            IPIV( K ) = KP
            GO TO 70
         END IF

         IF( K.GT.N )
     $      GO TO 70
         KSTEP = 2
*
*        Determine rows and columns to be interchanged
*
         ABSAKP1K = ABS( A( K+1, K ) )
*
*        IMAX1 is the row-index of the absolute value largest element in
*        row K+2 to N, column K.
*        IMAX2 is the row-index of the absolute value largest element in
*        row K+2 to N, column K+1.
*        COLMAX1 and COLMAX2 are their absolute values.
*
         IMAX1 = K+1 + ISAMAX( N-K-1, A( K+2, K ), 1 )
         COLMAX1 = ABS( A( IMAX1, K ) )
         IMAX2 = K+1 + ISAMAX( N-K-1, A( K+2, K+1 ), 1 )
         COLMAX2 = ABS( A( IMAX2, K+1 ) )
*
         IF( MAX(MAX( ABSAKP1K, COLMAX1 ), COLMAX2).EQ.ZERO ) THEN
*
*           Column K and K+1 is zero or underflow: set INFO and continue
*
            IF( INFO.EQ.0 )
     $         INFO = K
            KP = 0
            IPIV( K ) = KP

         ELSE
            IF( ABSAKP1K.GE.MAX( COLMAX1, COLMAX2 ) ) THEN
*
*              no interchange
*
               KP = 0
               IPIV( K ) = KP
            
            ELSE
               IF( COLMAX1.GE.COLMAX2 ) THEN
*
*              Absolute value largest element is in column K
*              Interchange rows and columns K+1 and IMAX1
*                  
                  KP = IMAX1
                  IPIV( K ) = KP

                  CALL SSWAP( IMAX1-K-2, A( IMAX1, K+2 ), LDA,
     $                     A( K+2, K+1 ), 1 )

                  CALL SSCAL( IMAX1-K-2, -ONE, A( IMAX1, K+2 ),
     $                     LDA )

                  CALL SSCAL( IMAX1-K-2, -ONE, A( K+2, K+1 ),
     $                     1 )
             
                  CALL SSWAP( N-IMAX1, A( IMAX1+1, IMAX1 ), 1,
     $                     A( IMAX1+1, K+1 ), 1 ) 

                  A( IMAX1, K+1 ) = -A( IMAX1, K+1 )

*                 
*                 Interchange rows K+1 and IMAX1 in first K-1 columns of A
*     
                  CALL SSWAP( K, A( K+1, 1 ), LDA, A( IMAX1, 1 ),
     $                     LDA )

               ELSE
*
*                 Absolute value largest element is in column K+1
*                 Interchange rows and columns K and K+1, then Interchange K+1 and IMAX2
*
                  KP = -IMAX2
                  IPIV( K ) = KP

                  CALL SSWAP( N-K-1, A( K+2, K ), 1, A( K+2, K+1 ),
     $                     1 )
             
                  A( K+1, K ) = -A( K+1, K )

                  CALL SSWAP( IMAX2-K-2, A( IMAX2, K+2 ), LDA,
     $                     A( K+2, K+1 ), 1 )

                  CALL SSCAL( IMAX2-K-2, -ONE, A( IMAX2, K+2 ),
     $                     LDA )

                  CALL SSCAL( IMAX2-K-2, -ONE, A( K+2, K+1 ),
     $                     1 )
             
                  CALL SSWAP( N-IMAX2, A( IMAX2+1, IMAX2 ), 1,
     $                     A( IMAX2+1, K+1 ), 1 ) 

                  A( IMAX2, K+1 ) = -A( IMAX2, K+1 )

*                 
*                 Interchange rows K and K+1, then K+1 and IMAX2 in first K-1 columns of A
*
                  CALL SSWAP( K-1, A( K, 1 ), LDA, A( K+1, 1 ),
     $                     LDA )

                  CALL SSWAP( K, A( K+1, 1 ), LDA, A( IMAX2, 1 ),
     $                     LDA )


               END If
            END If

*
*           Update the lower triangle of A22 (= A(k+2:n,k+2:n))
*
            D21 = ONE/A( K+1, K )

            DO 60 J = K+2, N
*
               WK = -A( J, K+1 )*D21
               WKP1 = A( J, K )*D21
*
               DO 50 I = K+2, J-1
                  A( J, I ) = A( J, I ) + A( I, K )*WK +
     $               A( I, K+1 )*WKP1
50             CONTINUE

60          CONTINUE

*
*           Update C*S^-1
*
            DO 90 J = K+2, N
               T = A( J, K )
               A( J, K ) = -A( J, K+1 )*D21
               A( J, K+1 ) = T*D21
90          CONTINUE
         END IF

         K = K + KSTEP
         GO TO 40
*
   70    CONTINUE
*
*
*        Put L21 in standard form by partially undoing the interchanges
*        of rows in columns 1:k-1 looping backwards from k-1 to 1
*
         J = K - 2
120      CONTINUE
*
*           Undo the interchanges (if any) of rows JJ and JP at each
*           step J
*
*           (Here, J is a diagonal index)

            IF( J.GT.1 ) THEN
               JJ = J
               JP = IPIV( J )

               IF( JP.LT.0 ) THEN
                  JP = -JP
*                 (Here, J is a diagonal index)
                  CALL SSWAP( J-1, A( JP, 1 ), LDA, A( JJ+1, 1 ), LDA )
                  CALL SSWAP( J-1, A( JJ+1, 1 ), LDA, A( JJ, 1 ), LDA )
               ELSEIF( JP.GT.0 ) THEN
                  CALL SSWAP( J-1, A( JP, 1 ), LDA, A( JJ+1, 1 ), LDA )
               END IF
               
            END IF
*        (NOTE: Here, J is used to determine row length. Length J
*        of the rows to swap back doesn't include diagonal element)

         J = J - 2
         IF( J.GT.1 )
     $      GO TO 120

*
      END IF

      RETURN
*
*     End of SKYTF2
*
      END
