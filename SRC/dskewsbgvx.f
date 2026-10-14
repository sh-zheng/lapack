*> \brief \b DSKEWSBGVX
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download DSKEWSBGVX + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dskewsbgvx.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dskewsbgvx.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dskewsbgvx.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE DSKEWSBGVX( JOBZ, RANGE, UPLO, N, KA, KB, AB, LDAB, BB,
*                          LDBB, Q, LDQ, VL, VU, IL, IU, ABSTOL, M, W, Z,
*                          LDZ, WORK, IWORK, IFAIL, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBZ, RANGE, UPLO
*       INTEGER            IL, INFO, IU, KA, KB, LDAB, LDBB, LDQ, LDZ, M,
*      $                   N, K
*       DOUBLE PRECISION   ABSTOL, VL, VU
*       ..
*       .. Array Arguments ..
*       INTEGER            IFAIL( * ), IWORK( * )
*       DOUBLE PRECISION   AB( LDAB, * ), BB( LDBB, * ), Q( LDQ, * ),
*      $                   W( * ), WORK( * ), Z( LDZ, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DSKEWSBGVX computes selected eigenvalues, and optionally, eigenvectors
*> of a real generalized skew-symmetric-definite banded eigenproblem, of
*> the form A*x=(lambda)*B*x.  Here A and B are assumed to be banded. A is
*> skew-symmetric, B is symmetric and positive definite.  Eigenvalues and
*> eigenvectors can be selected by specifying either all eigenvalues,
*> a range of values or a range of indices for the desired eigenvalues.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOBZ
*> \verbatim
*>          JOBZ is CHARACTER*1
*>          = 'N':  Compute eigenvalues only;
*>          = 'V':  Compute eigenvalues and eigenvectors.
*> \endverbatim
*>
*> \param[in] RANGE
*> \verbatim
*>          RANGE is CHARACTER*1
*>          = 'A': all eigenvalues will be found.
*>          = 'V': all eigenvalues in the half-open interval (VL,VU]
*>                 will be found.
*>          = 'I': the IL-th through IU-th eigenvalues will be found.
*> \endverbatim
*>
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          = 'U':  Upper triangles of A and B are stored;
*>          = 'L':  Lower triangles of A and B are stored.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices A and B.  N >= 0.
*> \endverbatim
*>
*> \param[in] KA
*> \verbatim
*>          KA is INTEGER
*>          The number of superdiagonals of the matrix A if UPLO = 'U',
*>          or the number of subdiagonals if UPLO = 'L'.  KA >= 1.
*> \endverbatim
*>
*> \param[in] KB
*> \verbatim
*>          KB is INTEGER
*>          The number of superdiagonals of the matrix B if UPLO = 'U',
*>          or the number of subdiagonals if UPLO = 'L'.  KB >= 0.
*> \endverbatim
*>
*> \param[in,out] AB
*> \verbatim
*>          AB is DOUBLE PRECISION array, dimension (LDAB, N)
*>          On entry, the upper or lower triangle of the skew-symmetric band
*>          matrix A, stored in the first ka rows (upper triangle case) or
*>          2:ka+1 rows (lower triangle case) of the array.  The j-th column
*>          of A is stored in the j-th column of the array AB as follows:
*>          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<j;
*>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<i<=min(n,j+ka).
*>
*>          On exit, the contents of AB are destroyed.
*> \endverbatim
*>
*> \param[in] LDAB
*> \verbatim
*>          LDAB is INTEGER
*>          The leading dimension of the array AB.  LDAB >= KA+1.
*> \endverbatim
*>
*> \param[in,out] BB
*> \verbatim
*>          BB is DOUBLE PRECISION array, dimension (LDBB, N)
*>          On entry, the upper or lower triangle of the symmetric band
*>          matrix B, stored in the first kb+1 rows of the array.  The
*>          j-th column of B is stored in the j-th column of the array BB
*>          as follows:
*>          if UPLO = 'U', BB(ka+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<j;
*>          if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<i<=min(n,j+kb).
*>
*>          On exit, the factor S from the split Cholesky factorization
*>          B = S**T*S, as returned by DPBSTF.
*> \endverbatim
*>
*> \param[in] LDBB
*> \verbatim
*>          LDBB is INTEGER
*>          The leading dimension of the array BB.  LDBB >= KB+1.
*> \endverbatim
*>
*> \param[out] Q
*> \verbatim
*>          Q is DOUBLE PRECISION array, dimension (LDQ, N)
*>          If JOBZ = 'V', the n-by-n matrix used in the reduction of
*>          A*x = (lambda)*B*x to standard form, i.e. C*x = (lambda)*x,
*>          and consequently C to tridiagonal form.
*>          If JOBZ = 'N', the array Q is not referenced.
*> \endverbatim
*>
*> \param[in] LDQ
*> \verbatim
*>          LDQ is INTEGER
*>          The leading dimension of the array Q.  If JOBZ = 'N',
*>          LDQ >= 1. If JOBZ = 'V', LDQ >= max(1,N).
*> \endverbatim
*>
*> \param[in] VL
*> \verbatim
*>          VL is DOUBLE PRECISION
*>
*>          If RANGE='V', the lower bound of the interval to
*>          be searched for eigenvalues.  Eigenvalues with positive
*>          imaginary part less than or equal to VL, or greater
*>          than VU, will not be returned.  VL < VU. If VL < 0,
*>          zero eigenvalues will be included.
*>          Not referenced if RANGE = 'A' or 'I'.
*> \endverbatim
*>
*> \param[in] VU
*> \verbatim
*>          VU is DOUBLE PRECISION
*>          If RANGE='V', the upper bound of the interval to
*>          be searched for eigenvalues.  Eigenvalues with positive
*>          imaginary part less than or equal to VL, or greater
*>          than VU, will not be returned.  VL < VU, VU >= 0.
*>          Not referenced if RANGE = 'A' or 'I'.
*> \endverbatim
*>
*> \param[in] IL
*> \verbatim
*>          IL is INTEGER
*>          If RANGE='I', the index of eigenvalue with the
*>          smallest positive imaginary part to be returned.
*>          1 <= IL <= IU <= ceil(N/2), if N > 0; IL = 1 and
*>          IU = 0 if N = 0.
*>          Not referenced if RANGE = 'A' or 'V'.
*> \endverbatim
*>
*> \param[in] IU
*> \verbatim
*>          IU is INTEGER
*>           If RANGE='I', the index of eigenvalue with the
*>          largest positive imaginary part to be returned.
*>          1 <= IL <= IU <= ceil(N/2), if N > 0; IL = 1 and
*>          IU = 0 if N = 0.
*>          Not referenced if RANGE = 'A' or 'V'.
*> \endverbatim
*>
*> \param[in] ABSTOL
*> \verbatim
*>          ABSTOL is DOUBLE PRECISION
*>          The absolute error tolerance for the eigenvalues.
*>          An approximate eigenvalue is accepted as converged
*>          when it is determined to lie in an interval [a,b]
*>          of width less than or equal to
*>
*>                  ABSTOL + EPS *   max( |a|,|b| ) ,
*>
*>          where EPS is the machine precision.  If ABSTOL is less than
*>          or equal to zero, then  EPS*|T|  will be used in its place,
*>          where |T| is the 1-norm of the tridiagonal matrix obtained
*>          by reducing A to tridiagonal form.
*>
*>          Eigenvalues will be computed most accurately when ABSTOL is
*>          set to twice the underflow threshold 2*SLAMCH('S'), not zero.
*>          If this routine returns with INFO>0, indicating that some
*>          eigenvectors did not converge, try setting ABSTOL to
*>          2*SLAMCH('S').
*> \endverbatim
*>
*> \param[out] M
*> \verbatim
*>          M is INTEGER
*>          The total number of eigenvalues found.  0 <= M <= N.
*>          If RANGE = 'A', M = N. If RANGE = 'I', when N is
*>          odd and IL = 1, M = 2 * (IU-IL) + 1, otherwise M = 
*>          2 * (IU-IL+1).
*> \endverbatim
*>
*> \param[out] W
*> \verbatim
*>          W is DOUBLE PRECISION array, dimension (N)
*>          On normal exit, the first M elements contain the selected
*>          eigenvalues in block descending order. Refer to SSKEWSYGV. 
*> \endverbatim
*>
*> \param[out] Z
*> \verbatim
*>          Z is DOUBLE PRECISION array, dimension (LDZ, N)
*>          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
*>          contain the block form orthonormal eigenvectors of the
*>          matrix A corresponding to the selected eigenvalues, with the
*>          i-th column of Z holding the block form eigenvector associated
*>          with W(i). The block form eigenvectors (v1, v2) corresponding
*>          to a conjugate purely imaginary eigenvalue pair (ri, -ri) means
*>          A*v1 = r*v2 and A*v2 = -r*v1.  The eigenvectors are
*>          normalized so Z**T*B*Z = I.
*>          If JOBZ = 'N', then Z is not referenced.
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of the array Z.  LDZ >= 1, and if
*>          JOBZ = 'V', LDZ >= max(1,N).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (8*N)
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (5*N)
*> \endverbatim
*>
*> \param[out] IFAIL
*> \verbatim
*>          IFAIL is INTEGER array, dimension (M)
*>          If JOBZ = 'V', then if INFO = 0, the first M elements of
*>          IFAIL are zero.  If INFO > 0, then IFAIL contains the
*>          indices of the eigenvalues that failed to converge.
*>          If JOBZ = 'N', then IFAIL is not referenced.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          <= N:    if INFO = i, then i eigenvectors failed to converge
*>                   in DSKEWSTEIN; their indices are stored in IFAIL.
*>          N+1..2N: if INFO = N + i, for 1 <= i <= N, then DPBSTF
*>                   returned an error code; the leading principal minor
*>                   of order i of B is not positive. No eigenvalues or
*>                   eigenvectors were computed.
*>          > 2N:    if INFO = 2*N + i, then DSKEWSTEBZ returned
*>                   INFO = i; see DSKEWSTEBZ for details.
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
*> \ingroup hbgvx
*
*> \par Contributors:
*  ==================
*>
*>     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA
*
*  =====================================================================
      SUBROUTINE DSKEWSBGVX( JOBZ, RANGE, UPLO, N, KA, KB, AB, LDAB, BB,
     $                   LDBB, Q, LDQ, VL, VU, IL, IU, ABSTOL, M, W, Z,
     $                   LDZ, WORK, IWORK, IFAIL, INFO )
      IMPLICIT NONE
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE, UPLO
      INTEGER            IL, INFO, IU, KA, KB, LDAB, LDBB, LDQ, LDZ, M,
     $                   N
      DOUBLE PRECISION   ABSTOL, VL, VU
*     ..
*     .. Array Arguments ..
      INTEGER            IFAIL( * ), IWORK( * )
      DOUBLE PRECISION   AB( LDAB, * ), BB( LDBB, * ), Q( LDQ, * ),
     $                   W( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALLEIG, INDEIG, TEST, UPPER, VALEIG, WANTZ
      CHARACTER          ORDER, VECT
      INTEGER            I, IINFO, INDE, INDEE, INDISP, K,
     $                   INDIWO, INDWRK, ITMP1, J, JJ, NSPLIT
      DOUBLE PRECISION   TMP1
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMV, DLACPY, DPBSTF, DSKEWSBGST,
     $                   DSKEWSBTRD, DSCAL,
     $                   DSKEWSTEBZ, DSKEWSTEIN, DSKEWSTEQR, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
      UPPER = LSAME( UPLO, 'U' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )
*
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( UPPER .OR. LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( KA.LT.1 ) THEN
         INFO = -5
      ELSE IF( KB.LT.0 .OR. KB.GT.KA ) THEN
         INFO = -6
      ELSE IF( LDAB.LT.KA+1 ) THEN
         INFO = -8
      ELSE IF( LDBB.LT.KB+1 ) THEN
         INFO = -10
      ELSE IF( LDQ.LT.1 .OR. ( WANTZ .AND. LDQ.LT.N ) ) THEN
         INFO = -12
      ELSE
         IF( VALEIG ) THEN
            IF( N.GT.0 .AND. VU.LE.VL )
     $         INFO = -14
         ELSE IF( INDEIG ) THEN
            IF( IL.LT.1 .OR. IL.GT.MAX( 1, N ) ) THEN
               INFO = -15
            ELSE IF ( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) THEN
               INFO = -16
            END IF
         END IF
      END IF
      IF( INFO.EQ.0) THEN
         IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
            INFO = -21
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSKEWSBGVX', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      M = 0
      IF( N.EQ.0 )
     $   RETURN
*
*     Form a split Cholesky factorization of B.
*
      CALL DPBSTF( UPLO, N, KB, BB, LDBB, INFO )
      IF( INFO.NE.0 ) THEN
         INFO = N + INFO
         RETURN
      END IF
*
*     Transform problem to standard eigenvalue problem.
*
      CALL DSKEWSBGST( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Q,
     $             LDQ, WORK, IINFO )
*
*     Reduce skew-symmetric band matrix to tridiagonal form.
*
      INDE = 1
      INDWRK = INDE + N
      IF( WANTZ ) THEN
         VECT = 'U'
      ELSE
         VECT = 'N'
      END IF
      CALL DSKEWSBTRD( VECT, UPLO, N, KA, AB, LDAB, WORK( INDE ),
     $             Q, LDQ, WORK( INDWRK ), IINFO )
      IF(UPPER)
     $   CALL DSCAL(N-1, -ONE, WORK( INDE ), 1)
*
*     If all eigenvalues are desired and ABSTOL is less than or equal
*     to zero, then CALL DSKEWSTEQR.  If this fails for some
*     eigenvalue, then try DSKEWSTEBZ.
*
      TEST = .FALSE.
      IF( INDEIG ) THEN
         IF( IL.EQ.1 .AND. IU.EQ.(N+1)/2 ) THEN
            TEST = .TRUE.
         END IF
      END IF
*      IF( ( ALLEIG .OR. TEST ) .AND. ( ABSTOL.LE.ZERO ) ) THEN
      IF( .FALSE. ) THEN
         INDEE = INDWRK + 2*N
         CALL DCOPY( N-1, WORK( INDE ), 1, WORK( INDEE ), 1 )
         CALL DLACPY( 'A', N, N, Q, LDQ, Z, LDZ )
         CALL DSKEWSTEQR( JOBZ, N, W, WORK( INDEE ), Z, LDZ,
     $                WORK( INDWRK ), INFO )
         CALL DCOPY( N-1, WORK( INDEE ), 1, W, 1 )
         W( N ) = ZERO
         IF( INFO.EQ.0 ) THEN
            DO 10 I = 1, N
               IFAIL( I ) = 0
   10         CONTINUE
         END IF
         IF( INFO.EQ.0 ) THEN
            M = N
            GO TO 30
         END IF
         INFO = 0
      END IF
*
*     Otherwise, call DSKEWSTEBZ and, if eigenvectors are desired,
*     call DSKEWSTEIN.
*
      IF( WANTZ ) THEN
         ORDER = 'B'
      ELSE
         ORDER = 'E'
      END IF
      INDISP = 1 + N
      INDIWO = INDISP + N
      CALL DSKEWSTEBZ( RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL,
     $             WORK( INDE ), M, NSPLIT, W,
     $             IWORK( 1 ), IWORK( INDISP ), WORK( INDWRK ),
     $             IWORK( INDIWO ), IINFO )
      IF( IINFO.NE.0 ) THEN
         INFO = 2*N + IINFO
         IF( IINFO.NE.1 )
     $      GO TO 30
      END IF
*
      IF( WANTZ ) THEN
         CALL DSKEWSTEIN( N, WORK( INDE ), M, W,
     $                IWORK( 1 ), IWORK( INDISP ), Z, LDZ,
     $                WORK( INDWRK ), IWORK( INDIWO ), IFAIL, IINFO )
         IF( IINFO.NE.0 .AND. INFO.EQ.0 )
     $      INFO = IINFO
*
*        Apply transformation matrix used in reduction to tridiagonal
*        form to eigenvectors returned by DSKEWSTEIN.
*
         DO 20 J = 1, M
            CALL DCOPY( N, Z( 1, J ), 1, WORK( 1 ), 1 )
            CALL DGEMV( 'N', N, N, ONE, Q, LDQ, WORK, 1, ZERO,
     $                  Z( 1, J ), 1 )
   20    CONTINUE
      END IF
*
   30 CONTINUE
*
*     If eigenvalues are not in order, then sort them, along with
*     eigenvectors.
*
      IF( WANTZ ) THEN
         JJ = 1
         DO WHILE(JJ.LT.(M-1))
            IF(W(JJ).EQ.ZERO) THEN
               DO K = JJ+1,M-1,2
                  IF(W(K).EQ.ZERO) THEN
                     DO I = JJ, K-2
                        W(I) = W(I+1)
                        CALL DSWAP( N, Z( 1, I ), 1, Z( 1, I+1 ),
     $                              1 )
                        IF( INFO.NE.0 ) THEN
                           ITMP1 = IFAIL( I )
                           IFAIL( I ) = IFAIL( I+1 )
                           IFAIL( I+1 ) = ITMP1
                        END IF
                     END DO
                     W(K-1) = ZERO
                     JJ = K+1
                     EXIT
                  ELSEIF(MOD(M,2).EQ.1 .AND. K.EQ.(M-1)) THEN
                     DO I = JJ, K-1
                        W(I) = W(I+1)
                        CALL DSWAP( N, Z( 1, I ), 1, Z( 1, I+1 ),
     $                              1 )
                        IF( INFO.NE.0 ) THEN
                           ITMP1 = IFAIL( I )
                           IFAIL( I ) = IFAIL( I+1 )
                           IFAIL( I+1 ) = ITMP1
                        END IF
                     END DO
                     CALL DSWAP( N, Z( 1, K ), 1, Z( 1, K+1 ),
     $                           1 )
                     IF( INFO.NE.0 ) THEN
                        ITMP1 = IFAIL( K )
                        IFAIL( K ) = IFAIL( K+1 )
                        IFAIL( K+1 ) = ITMP1
                     END IF
                     W(K) = ZERO
                     JJ = K+1
                     EXIT
                  ELSEIF(MOD(M,2).EQ.0 .AND. K.EQ.(M-2)) THEN
                     DO I = JJ, K-1
                        W(I) = W(I+1)
                        CALL DSWAP( N, Z( 1, I ), 1, Z( 1, I+1 ),
     $                              1 )
                        IF( INFO.NE.0 ) THEN
                           ITMP1 = IFAIL( I )
                           IFAIL( I ) = IFAIL( I+1 )
                           IFAIL( I+1 ) = ITMP1
                        END IF
                     END DO
                     CALL DSWAP( N, Z( 1, K ), 1, Z( 1, K+1 ), 1 )
                     IF( INFO.NE.0 ) THEN
                        ITMP1 = IFAIL( K )
                        IFAIL( K ) = IFAIL( K+1 )
                        IFAIL( K+1 ) = ITMP1
                     END IF
                     W(K) = ZERO
                     JJ = K+1
                     EXIT
                  END IF
               END DO
               IF (JJ.LT.(M-1)) THEN
                  CYCLE
               END IF
            END IF
            JJ = JJ+2
         END DO
*
         DO 60 JJ = 1, M-1, 2
            I = JJ
            TMP1 = ABS(W(JJ))
            DO 50 K = JJ+2, M-1, 2
               IF(ABS(W(K)).GT.TMP1) THEN
                  I = K
                  TMP1 = ABS(W(K))
               END IF
   50       CONTINUE
            IF(I.NE.JJ) THEN
               CALL DSWAP( 1, W( I ), 1, W( JJ ), 1 )
               CALL DSWAP( N, Z( 1, I ), 1, Z( 1, JJ ), 1 )
               CALL DSWAP( N, Z( 1, I+1 ), 1, Z( 1, JJ+1 ), 1 )
               IF( INFO.NE.0 ) THEN
                  ITMP1 = IFAIL( I )
                  IFAIL( I ) = IFAIL( JJ )
                  IFAIL( JJ ) = ITMP1
                  ITMP1 = IFAIL( I+1 )
                  IFAIL( I+1 ) = IFAIL( JJ+1 )
                  IFAIL( JJ+1 ) = ITMP1
               END IF
            END IF
   60    CONTINUE
      END IF
*
      RETURN
*
*     End of DSKEWSBGVX
*
      END
