*> \brief \b SSKEWSBGV
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download SSKEWSBGV + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sskewsbgv.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sskewsbgv.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sskewsbgv.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE SSKEWSBGV( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, W, Z,
*                         LDZ, WORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBZ, UPLO
*       INTEGER            INFO, KA, KB, LDAB, LDBB, LDZ, N
*       ..
*       .. Array Arguments ..
*       REAL               AB( LDAB, * ), BB( LDBB, * ), W( * ),
*      $                   WORK( * ), Z( LDZ, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SSKEWSBGV computes all the eigenvalues, and optionally, the eigenvectors
*> of a real generalized skew-symmetric-definite banded eigenproblem, of
*> the form A*x=(lambda)*B*x. Here A and B are assumed to be banded. A is
*> skew-symmetric, B is symmetric and positive definite.
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
*>          or the number of subdiagonals if UPLO = 'L'. KA >= 1.
*> \endverbatim
*>
*> \param[in] KB
*> \verbatim
*>          KB is INTEGER
*>          The number of superdiagonals of the matrix B if UPLO = 'U',
*>          or the number of subdiagonals if UPLO = 'L'. KB >= 0.
*> \endverbatim
*>
*> \param[in,out] AB
*> \verbatim
*>          AB is REAL array, dimension (LDAB, N)
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
*>          BB is REAL array, dimension (LDBB, N)
*>          On entry, the upper or lower triangle of the symmetric band
*>          matrix B, stored in the first kb+1 rows of the array.  The
*>          j-th column of B is stored in the j-th column of the array BB
*>          as follows:
*>          if UPLO = 'U', BB(kb+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<=j;
*>          if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb).
*>
*>          On exit, the factor S from the split Cholesky factorization
*>          B = S**T*S, as returned by SPBSTF.
*> \endverbatim
*>
*> \param[in] LDBB
*> \verbatim
*>          LDBB is INTEGER
*>          The leading dimension of the array BB.  LDBB >= KB+1.
*> \endverbatim
*>
*> \param[out] W
*> \verbatim
*>          W is REAL array, dimension (N)
*>          If INFO = 0, the (N-1) lower subdiagonal elements of the
*>          block diagonal matrix at front, and zero at last.
*>		    The matrix consists of 2-by-2 skew-symmetric blocks, and zeros.
*>          The values in W, which represent blocks, are always
*>          positive, and sorted in descending order.
*>          The eigenvalues of each blocks can be evaluated directly.
*> \endverbatim
*>
*> \param[out] Z
*> \verbatim
*>          Z is REAL array, dimension (LDZ, N)
*>          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z,
*>          which leads to the block diagonal form in W. The matrix are
*>          normalized so that Z**T*B*Z = I. The eigenvectors of the
*>          matrix can be evaluated directly.
*>          If JOBZ = 'N', then Z is not referenced.
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of the array Z.  LDZ >= 1, and if
*>          JOBZ = 'V', LDZ >= N.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (3*N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, and i is:
*>             <= N:  the algorithm failed to converge:
*>                    i off-diagonal elements of an intermediate
*>                    tridiagonal form did not converge to zero;
*>             > N:   if INFO = N + i, for 1 <= i <= N, then SPBSTF
*>                    returned INFO = i: B is not positive definite.
*>                    The factorization of B could not be completed and
*>                    no eigenvalues or eigenvectors were computed.
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
*> \ingroup hbgv
*
*  =====================================================================
      SUBROUTINE SSKEWSBGV( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB,
     $                  LDBB, W, Z,
     $                  LDZ, WORK, INFO )
      IMPLICIT NONE
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, KA, KB, LDAB, LDBB, LDZ, N
*     ..
*     .. Array Arguments ..
      REAL               AB( LDAB, * ), BB( LDBB, * ), W( * ),
     $                   WORK( * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
*     .. Local Scalars ..
      LOGICAL            UPPER, WANTZ
      CHARACTER          VECT
      INTEGER            IINFO, INDE, INDWRK
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           SPBSTF, SSKEWSBGST, SSKEWSBTRD, SSKEWSTEQR,
     $                   SCOPY, SSCAL, XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
      UPPER = LSAME( UPLO, 'U' )
*
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( UPPER .OR. LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( KA.LT.1 ) THEN
         INFO = -4
      ELSE IF( KB.LT.0 .OR. KB.GT.KA ) THEN
         INFO = -5
      ELSE IF( LDAB.LT.KA+1 ) THEN
         INFO = -7
      ELSE IF( LDBB.LT.KB+1 ) THEN
         INFO = -9
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSKEWSBGV', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Form a split Cholesky factorization of B.
*
      CALL SPBSTF( UPLO, N, KB, BB, LDBB, INFO )
      IF( INFO.NE.0 ) THEN
         INFO = N + INFO
         RETURN
      END IF
*
*     Transform problem to standard eigenvalue problem.
*
      INDE = 1
      INDWRK = INDE + N
      CALL SSKEWSBGST( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Z,
     $             LDZ, WORK( INDWRK ), IINFO )
*
*     Reduce to tridiagonal form.
*
      IF( WANTZ ) THEN
         VECT = 'U'
      ELSE
         VECT = 'N'
      END IF
      CALL SSKEWSBTRD( VECT, UPLO, N, KA, AB, LDAB, WORK( INDE ), Z,
     $             LDZ,
     $             WORK( INDWRK ), IINFO )
*
*     Call SSKEWSTEQR.
*
      IF(UPPER)
     $   CALL SSCAL( N-1, -ONE, WORK( INDE ), 1 )
      CALL SCOPY( N-1, WORK( INDE ), 1, W, 1 )
      CALL SSKEWSTEQR( JOBZ, N, W, Z, LDZ,
     $             WORK( INDWRK ),
     $             INFO )
      W(N) = ZERO
      RETURN
*
*     End of SSKEWSBGV
*
      END
