*> \brief <b> SKYEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for KY matrices</b>
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download SKYEV + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/skyev.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/skyev.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/skyev.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE SKYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBZ, UPLO
*       INTEGER            INFO, LDA, LWORK, N
*       ..
*       .. Array Arguments ..
*       REAL               A( LDA, * ), W( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SKYEV computes all eigenvalues and, optionally, eigenvectors of a
*> real skew-symmetric matrix A.
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
*>          = 'U':  Upper triangle of A is stored;
*>          = 'L':  Lower triangle of A is stored.
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
*>          A is REAL array, dimension (LDA, N)
*>          On entry, the skew-symmetric matrix A.  If UPLO = 'U', the
*>          strictly N-by-N upper triangular part of A contains the
*>          upper triangular part of the matrix A.  If UPLO = 'L',
*>          the strictly N-by-N lower triangular part of A contains
*>          the lower triangular part of the matrix A.
*>          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
*>          orthonormal eigenvectors of the matrix A.
*>          The eigenvectors of skew-symmetric matrix are conjugate
*>          complex (if the corresponding eigenvalues are conjugate pure
*>          imaginary) or real (if the corresponding eigenvalues are 0).
*>          A stores half of the complex (corresponding to the eigenvalues
*>          with positive imaginary part) and all the real eigenvectors.
*>          The complex eigenvector needs twice the storage than the real
*>          eigenvector. They are stored in the same position as the
*>          eigenvalues in W.
*>          If JOBZ = 'N', then on exit the strictly lower triangle
*>          (if UPLO='L') or the upper triangle (if UPLO='U') of A,
*>          is destroyed.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] W
*> \verbatim
*>          W is REAL array, dimension (N)
*>          If INFO = 0, the eigenvalues.
*>          If the matrix is not singular, the elements in W are always
*>          organized in pairs, and the positive imaginary part of the
*>          eigenvalues are stored in the first element in pairs, in
*>          ascending order. The conjugate negative imaginary part are
*>          stored in the second element in pairs.
*>          If the matrix is singular, the 0 eigenvalues are stored
*>          firstly in W, each takes one element, and the other
*>          eigenvalues are stored in pairs, as mentioned above.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The length of the array WORK.  If JOBZ = 'N',
*>          LWORK >= max(1,3*N). For optimal efficiency, LWORK >= (NB+2)*N.
*>          If JOBZ = 'V', LWORK >= max(1,N**2+3*N). For optimal efficiency,
*>          LWORK >= N**2+(NB+2)*N. where NB is the blocksize for SKYTRD
*>          returned by ILAENV.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, the algorithm failed to converge; i
*>                off-diagonal elements of an intermediate tridiagonal
*>                form did not converge to zero.
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
*> \ingroup realKYeigen
*
*  =====================================================================
      SUBROUTINE SKYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), W( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER, LQUERY, WANTZ
      INTEGER            IINFO, IMAX, INDE, INDTAU, INDWRK, ISCALE,
     $                   LLWORK, LWKOPT, NB, MADDI, LMWORK, I, J,
     $                   RECZ, RECP, LTWORK, NBT, CRES
      REAL               ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA,
     $                   SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      REAL               SLAMCH, SLANSY
      EXTERNAL           ILAENV, LSAME, SLAMCH, SLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLASCL, SORGTR, SSCAL, SSTEQR, SSTERF, SKYTRD,
     $                   SCOPY, SAXPY, SGEMM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
      LQUERY = ( LWORK.EQ.-1 )
*
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
*
      MADDI = 0
      IF( WANTZ )
     $   MADDI = N**2

      IF( INFO.EQ.0 ) THEN
         NB = ILAENV( 1, 'SKYTRD', UPLO, N, -1, -1, -1 )
         LWKOPT = MAX( 1, MADDI+( NB+2 )*N )
         WORK( 1 ) = LWKOPT
*
         IF( LWORK.LT.MAX( 1, MADDI+3*N ) .AND. .NOT.LQUERY )
     $      INFO = -8
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SKYEV ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         RETURN
      END IF
*
      IF( N.EQ.1 ) THEN
         W( 1 ) = ZERO
         WORK( 1 ) = 3
         IF( WANTZ )
     $      A( 1, 1 ) = ONE
         RETURN
      END IF
*
*     Get machine constants.
*
      SAFMIN = SLAMCH( 'Safe minimum' )
      EPS = SLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = SQRT( BIGNUM )
*
*     Scale matrix to allowable range, if necessary.
*
      ANRM = SLANSY( 'M', UPLO, N, A, LDA, WORK )
      ISCALE = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 )
     $   CALL SLASCL( UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO )
*
*     Call SKYTRD to reduce skew-symmetric matrix to tridiagonal form.
*
      LMWORK = 1
      INDE = LMWORK + MADDI
      INDTAU = INDE + N
      INDWRK = INDTAU + N
      LLWORK = LWORK - INDWRK + 1
      LTWORK = LWORK - INDE + 1
      NBT = LTWORK / (2*N)
      CALL SKYTRD( UPLO, N, A, LDA, W, WORK( INDE ), WORK( INDTAU ),
     $             WORK( INDWRK ), LLWORK, IINFO )
*
*     For eigenvalues only, call SSTERF.  For eigenvectors, first call
*     SORGTR to generate the orthogonal matrix, then call SSTEQR.
*
      IF( .NOT.WANTZ ) THEN
*
*        Convert W into symmetric tridiagonal form, using the upper off-diagonal elements.
*        
         IF( LOWER ) THEN
            CALL SSCAL( N-1, -ONE, WORK( INDE ), 1)
         END IF
*
         CALL SSTERF( N, W, WORK( INDE ), INFO )
*
*        Restore W into skew-symmetric form.
*   
         J = 1
         IF( MOD(N,2) .NE. 0 ) THEN
            W( J ) = ZERO
            J = J+1
         END IF
         DO 20 I = (N+1)/2+1, N
            IF( W( I ) .LE. ZERO ) THEN
               W( I ) = ZERO
            END IF
            W( J ) = W( I )
            W( J+1 ) = -W( J )
            J = J+2
   20    CONTINUE
      ELSE
         CALL SORGTR( UPLO, N, A, LDA, WORK( INDTAU ), WORK( INDWRK ),
     $                LLWORK, IINFO )
*
*        Rearrange A.
*  
         DO 30 I = 1, (N-1)/2+1
            CALL SCOPY( N, A(1, 2*(I-1)+1), 1, WORK(LMWORK+(I-1)*N), 1 )
            IF( MOD(I+1, 2).EQ.1 )
     $         CALL SSCAL( N, -ONE, WORK(LMWORK+(I-1)*N), 1)
   30    CONTINUE
*
         DO 40 I = (N-1)/2+2, N
            CALL SCOPY( N, A(1, 2*(I-((N-1)/2+1))), 1,
     $           WORK(LMWORK+(I-1)*N), 1 )
            IF( MOD(I-((N-1)/2), 2).EQ.1 )
     $         CALL SSCAL( N, -ONE, WORK(LMWORK+(I-1)*N), 1)
   40    CONTINUE
*
         CALL SCOPY( N**2, WORK(LMWORK), 1, A, 1 )
*
*        Convert W into symmetric tridiagonal form, using the upper off-diagonal elements.
*        
         IF( LOWER ) THEN
            CALL SSCAL( N-1, -ONE, WORK( INDE ), 1)
         END IF
*
         CALL SSTEQR( 'I', N, W, WORK( INDE ), WORK( LMWORK ), LDA,
     $               WORK( INDTAU ), INFO )
*
*        Restore W into skew-symmetric form and evaluate corresponding eigenvectors.
*  
         J = 1
         RECZ = N/2+1
         RECP = N/2+1
         IF( MOD(N,2) .NE. 0 ) THEN
            W( J ) = ZERO
            J = J+1
            RECZ = N/2+1
            RECP = N/2+2
         END IF
         DO 50 I = (N+1)/2+1, N
            IF( W( I ) .LE. ZERO ) THEN
               W( I ) = ZERO
               RECZ = RECZ - 1
               RECP = RECP + 1
            END IF
            W( J ) = W( I )
            W( J+1 ) = -W( J )
            J = J+2
   50    CONTINUE
*
         DO 60 I = RECZ, RECP-1
            CALL SCOPY( (N-1)/2+1, WORK(LMWORK+(I-1)*N), 2,
     $                 WORK(INDE), 1 )
            CALL SCOPY( N-((N-1)/2+1), WORK(LMWORK+(I-1)*N+1), 2,
     $                 WORK(INDE+(N-1)/2+1), 1 )
            CALL SCOPY( N, WORK(INDE), 1, WORK(LMWORK+(I-RECZ)*N), 1 )
   60    CONTINUE
*
         DO 70 I = RECP, N
            CALL SCOPY( (N-1)/2+1, WORK(LMWORK+(I-1)*N), 2,
     $                 WORK(INDE), 1 )
            CALL SCOPY( N-((N-1)/2+1), WORK(LMWORK+(I-1)*N+1), 2,
     $                 WORK(INDE+(N-1)/2+1), 1 )
            CALL SCOPY( N, WORK(INDE), 1, WORK(LMWORK+(I-1)*N), 1 )
   70    CONTINUE
*
         CRES = RECP-RECZ
         J = 0
         DO 80 I = 1, CRES/NBT
            CALL SGEMM( 'N', 'N', N, NBT, (N-1)/2+1, ONE, A(1, 1),
     $                 LDA, WORK(LMWORK+J*NBT*N), N,
     $                 ZERO, WORK(INDE), N )
            CALL SGEMM( 'N', 'N', N, NBT, N-((N-1)/2+1), ONE,
     $                 A(1, (N-1)/2+2), LDA,
     $                 WORK(LMWORK+J*NBT*N+(N-1)/2+1),
     $                 N, ZERO, WORK(INDE+NBT*N), N )
            CALL SCOPY( NBT*N, WORK(INDE), 1,
     $                 WORK(LMWORK+J*NBT*N), 1 )
            CALL SAXPY( NBT*N, ONE, WORK(INDE+NBT*N), 1,
     $                 WORK(LMWORK+J*NBT*N), 1 )
            J = J+1
   80    CONTINUE
         IF (CRES-J*NBT > 0) THEN
            CALL SGEMM( 'N', 'N', N, CRES-J*NBT, (N-1)/2+1, ONE,
     $                 A(1, 1), LDA, WORK(LMWORK+J*NBT*N),
     $                 N, ZERO, WORK(INDE), N )
            CALL SGEMM( 'N', 'N', N, CRES-J*NBT, N-((N-1)/2+1),
     $                 ONE, A(1, (N-1)/2+2), LDA,
     $                 WORK(LMWORK+J*NBT*N+(N-1)/2+1),
     $                 N, ZERO, WORK(INDE+NBT*N), N )
            CALL SCOPY( (CRES-J*NBT)*N, WORK(INDE), 1,
     $                 WORK( LMWORK+J*NBT*N ), 1 )
            CALL SAXPY( (CRES-J*NBT)*N, ONE, WORK(INDE+NBT*N),
     $                 1, WORK( LMWORK+J*NBT*N ), 1 )
         END IF
*
         CRES = MIN(RECZ-1, N-RECP+1)
         J = 0
         DO 90 I = 1, CRES/NBT
            CALL SGEMM( 'N', 'N', N, NBT, (N-1)/2+1, ONE, A(1, 1),
     $                 LDA, WORK(LMWORK+(RECP-1+J*NBT)*N),
     $                 N, ZERO, WORK(INDE), N )
            CALL SGEMM( 'N', 'N', N, NBT, N-((N-1)/2+1),
     $                 ONE, A(1, (N-1)/2+2), LDA, 
     $                 WORK(LMWORK+(RECP-1+J*NBT)*N+(N-1)/2+1),
     $                 N, ZERO, WORK(INDE+NBT*N), N )
            CALL SCOPY( NBT*N, WORK(INDE), 1,
     $                 WORK(LMWORK+(N-2*CRES+J*2*NBT)*N), 2)
            CALL SCOPY( NBT*N, WORK(INDE+NBT*N), 1,
     $                 WORK(LMWORK+(N-2*CRES+J*2*NBT)*N+1), 2)
            J = J+1
   90    CONTINUE
         IF (CRES-J*NBT > 0) THEN
            CALL SGEMM( 'N', 'N', N, CRES-J*NBT, (N-1)/2+1, ONE,
     $                 A(1, 1), LDA, WORK(LMWORK+(RECP-1+J*NBT)*N),
     $                 N, ZERO, WORK(INDE), N )
            CALL SGEMM( 'N', 'N', N, CRES-J*NBT, N-((N-1)/2+1),
     $                 ONE, A(1, (N-1)/2+2), LDA,
     $                 WORK(LMWORK+(RECP-1+J*NBT)*N+(N-1)/2+1),
     $                 N, ZERO, WORK(INDE+NBT*N), N )
            CALL SCOPY( (CRES-J*NBT)*N, WORK(INDE), 1,
     $                 WORK(LMWORK+(N-2*CRES+J*2*NBT)*N), 2 )
            CALL SCOPY( (CRES-J*NBT)*N, WORK(INDE+NBT*N), 1,
     $                 WORK(LMWORK+(N-2*CRES+J*2*NBT)*N+1), 2 )
         END IF
*
         CALL SCOPY( N**2, WORK(LMWORK), 1, A, 1 )
      END IF
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = N
         ELSE
            IMAX = INFO - 1
         END IF
         CALL SSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF
*
*     Set WORK(1) to optimal workspace size.
*
      WORK( 1 ) = LWKOPT
*
      RETURN
*
*     End of SKYEV
*
      END
