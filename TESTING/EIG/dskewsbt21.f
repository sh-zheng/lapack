*> \brief \b DSKEWSBT21
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DSKEWSBT21( UPLO, N, KA, KS, A, LDA, E, U, LDU, WORK,
*                          RESULT )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            KA, KS, LDA, LDU, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), E( * ), RESULT( 2 ),
*      $                   U( LDU, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DSKEWSBT21  generally checks a decomposition of the form
*>
*>         A = U S U**T
*>
*> where **T means transpose, A is skew-symmetric banded, U is
*> orthogonal, and S is diagonal (if KS=0) or skew-symmetric
*> tridiagonal (if KS=1).
*>
*> Specifically:
*>
*>         RESULT(1) = | A - U S U**T | / ( |A| n ulp ) and
*>         RESULT(2) = | I - U U**T | / ( n ulp )
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER
*>          If UPLO='U', the upper triangle of A and V will be used and
*>          the (strictly) lower triangle will not be referenced.
*>          If UPLO='L', the lower triangle of A and V will be used and
*>          the (strictly) upper triangle will not be referenced.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The size of the matrix.  If it is zero, DSKEWSBT21 does nothing.
*>          It must be at least zero.
*> \endverbatim
*>
*> \param[in] KA
*> \verbatim
*>          KA is INTEGER
*>          The bandwidth of the matrix A.  It must be at least zero.  If
*>          it is larger than N-1, then max( 0, N-1 ) will be used.
*> \endverbatim
*>
*> \param[in] KS
*> \verbatim
*>          KS is INTEGER
*>          The bandwidth of the matrix S.  It may only be zero or one.
*>          If zero, then S is diagonal, and E is not referenced.  If
*>          one, then S is skew-symmetric tri-diagonal.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA, N)
*>          The original (unfactored) matrix.  It is assumed to be
*>          skew-symmetric, and only the upper (UPLO='U') or only the lower
*>          (UPLO='L') will be referenced.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of A.  It must be at least 1
*>          and at least min( KA, N-1 ).
*> \endverbatim
*>
*> \param[in] E
*> \verbatim
*>          E is DOUBLE PRECISION array, dimension (N-1)
*>          The off-diagonal of the (skew-symmetric tri-) diagonal matrix S.
*>          If UPLO = 'U', E(1) is the (1,2), E(2) is the (2,3) element, etc.
*>          If UPLO = 'L', E(1) is the (2,1), E(2) is the (3,2) element, etc.
*>          Not referenced if KS=0.
*> \endverbatim
*>
*> \param[in] U
*> \verbatim
*>          U is DOUBLE PRECISION array, dimension (LDU, N)
*>          The orthogonal matrix in the decomposition, expressed as a
*>          dense matrix (i.e., not as a product of Householder
*>          transformations, Givens transformations, etc.)
*> \endverbatim
*>
*> \param[in] LDU
*> \verbatim
*>          LDU is INTEGER
*>          The leading dimension of U.  LDU must be at least N and
*>          at least 1.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (N**2+N)
*> \endverbatim
*>
*> \param[out] RESULT
*> \verbatim
*>          RESULT is DOUBLE PRECISION array, dimension (2)
*>          The values computed by the two tests described above.  The
*>          values are currently limited to 1/ulp, to avoid overflow.
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
*> \ingroup double_eig
*
*  =====================================================================
      SUBROUTINE DSKEWSBT21( UPLO, N, KA, KS, A, LDA, E, U, LDU, WORK,
     $                   RESULT )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            KA, KS, LDA, LDU, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), E( * ), RESULT( 2 ),
     $                   U( LDU, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER
      CHARACTER          CUPLO
      INTEGER            IKA, I, J, K, JC, JR, LW
      DOUBLE PRECISION   ANORM, ULP, UNFL, WNORM, ALPHA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANGE, DLANSKEWSB
      EXTERNAL           LSAME, DLAMCH, DLANGE, DLANSKEWSB
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Constants
*
      RESULT( 1 ) = ZERO
      RESULT( 2 ) = ZERO
      IF( N.LE.0 )
     $   RETURN
*
      IKA = MAX( 0, MIN( N-1, KA ) )
      LW = ( N*( N+1 ) ) / 2
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         LOWER = .FALSE.
         CUPLO = 'U'
      ELSE
         LOWER = .TRUE.
         CUPLO = 'L'
      END IF
*
      UNFL = DLAMCH( 'Safe minimum' )
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
*
*     Some Error Checks
*
*     Do Test 1
*
*     Norm of A:
*
      ANORM = MAX( DLANSKEWSB( '1', CUPLO, N, IKA, A, LDA, WORK ),
     $             UNFL )
*
*     Compute error matrix:    Error = A - U S U**T
*
*     Copy A from SB to SP storage format.
*
      J = 0
      DO 50 JC = 1, N
         IF( LOWER ) THEN
            J = J + 1
            WORK( J ) = ZERO
            DO 10 JR = 2, MIN( IKA+1, N+1-JC )
               J = J + 1
               WORK( J ) = A( JR, JC )
   10       CONTINUE
            DO 20 JR = IKA + 2, N + 1 - JC
               J = J + 1
               WORK( J ) = ZERO
   20       CONTINUE
         ELSE
            DO 30 JR = IKA + 2, JC
               J = J + 1
               WORK( J ) = ZERO
   30       CONTINUE
            DO 40 JR = MIN( IKA, JC-1 ), 0, -1
               J = J + 1
               IF( JR .EQ. 0 ) THEN
                  WORK( J ) = ZERO
               ELSE
                  WORK( J ) = A( IKA+1-JR, JC )
               END IF
   40       CONTINUE
         END IF
   50 CONTINUE
*
      IF( N.GT.1 .AND. KS.EQ.1 ) THEN
         DO 70 J = 1, N - 1
            IF( LOWER ) THEN
               ALPHA = E( J )
            ELSE
               ALPHA = -E( J )
            END IF
            K = 0
            IF( LOWER ) THEN
               DO 75 JC = 1, N
                  DO 74 JR = JC, N
                     K = K + 1
                     WORK( K ) = WORK( K ) +
     $                  ALPHA*( U( JR, J )*U( JC, J+1 ) -
     $                          U( JR, J+1 )*U( JC, J ) )
   74             CONTINUE
   75          CONTINUE
            ELSE
               DO 77 JC = 1, N
                  DO 76 JR = 1, JC
                     K = K + 1
                     WORK( K ) = WORK( K ) +
     $                  ALPHA*( U( JR, J )*U( JC, J+1 ) -
     $                          U( JR, J+1 )*U( JC, J ) )
   76             CONTINUE
   77          CONTINUE
            END IF
   70    CONTINUE
      END IF
*
*     Compute WNORM = skew-symmetric 1-norm of error matrix
*
      DO 80 I = 1, N
         WORK( LW+I ) = ZERO
   80 CONTINUE
      K = 0
      IF( LOWER ) THEN
         DO 82 JC = 1, N
            DO 81 JR = JC, N
               K = K + 1
               IF( JR.NE.JC ) THEN
                  ALPHA = ABS( WORK( K ) )
                  WORK( LW+JC ) = WORK( LW+JC ) + ALPHA
                  WORK( LW+JR ) = WORK( LW+JR ) + ALPHA
               END IF
   81       CONTINUE
   82    CONTINUE
      ELSE
         DO 84 JC = 1, N
            DO 83 JR = 1, JC
               K = K + 1
               IF( JR.NE.JC ) THEN
                  ALPHA = ABS( WORK( K ) )
                  WORK( LW+JC ) = WORK( LW+JC ) + ALPHA
                  WORK( LW+JR ) = WORK( LW+JR ) + ALPHA
               END IF
   83       CONTINUE
   84    CONTINUE
      END IF
      WNORM = ZERO
      DO 85 I = 1, N
         WNORM = MAX( WNORM, WORK( LW+I ) )
   85 CONTINUE
*
      IF( ANORM.GT.WNORM ) THEN
         RESULT( 1 ) = ( WNORM / ANORM ) / ( N*ULP )
      ELSE
         IF( ANORM.LT.ONE ) THEN
            RESULT( 1 ) = ( MIN( WNORM, N*ANORM ) / ANORM ) / ( N*ULP )
         ELSE
            RESULT( 1 ) = MIN( WNORM / ANORM, DBLE( N ) ) / ( N*ULP )
         END IF
      END IF
*
*     Do Test 2
*
*     Compute  U U**T - I
*
      CALL DGEMM( 'N', 'C', N, N, N, ONE, U, LDU, U, LDU, ZERO, WORK,
     $            N )
*
      DO 90 J = 1, N
         WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - ONE
   90 CONTINUE
*
      RESULT( 2 ) = MIN( DLANGE( '1', N, N, WORK, N, WORK( N**2+1 ) ),
     $              DBLE( N ) ) / ( N*ULP )
*
      RETURN
*
*     End of DSKEWSBT21
*
      END
