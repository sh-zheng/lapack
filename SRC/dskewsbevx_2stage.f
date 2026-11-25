*> \brief <b> DSKEWSBEVX_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b>
*
*  @generated from dsbevx_2stage.f, fortran d -> s, Sat Nov  5 23:58:06 2016
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download DSKEWSBEVX_2STAGE + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dskewsbevx_2stage.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dskewsbevx_2stage.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dskewsbevx_2stage.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE DSKEWSBEVX_2STAGE( JOBZ, RANGE, UPLO, N, KD, AB, LDAB, Q,
*                                 LDQ, VL, VU, IL, IU, ABSTOL, M, W, Z,
*                                 LDZ, WORK, LWORK, IWORK, IFAIL, INFO )
*
*       IMPLICIT NONE
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBZ, RANGE, UPLO
*       INTEGER            IL, INFO, IU, KD, LDAB, LDQ, LDZ, M, N, LWORK
*       DOUBLE PRECISION               ABSTOL, VL, VU
*       ..
*       .. Array Arguments ..
*       INTEGER            IFAIL( * ), IWORK( * )
*       DOUBLE PRECISION   AB( LDAB, * ), Q( LDQ, * ), W( * ), WORK( * ),
*      $                   Z( LDZ, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DSKEWSBEVX_2STAGE computes selected eigenvalues and, optionally, eigenvectors
*> of a real skew-symmetric band matrix A using the 2stage technique for
*> the reduction to tridiagonal. Eigenvalues and eigenvectors can
*> be selected by specifying either a range of values or a range of
*> indices for the desired eigenvalues.
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
*>                  Not available in this release.
*> \endverbatim
*>
*> \param[in] RANGE
*> \verbatim
*>          RANGE is CHARACTER*1
*>          = 'A': all eigenvalues will be found;
*>          = 'V': all eigenvalues in the half-open interval (VL,VU]
*>                 will be found;
*>          = 'I': the IL-th through IU-th eigenvalues will be found.
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
*> \param[in] KD
*> \verbatim
*>          KD is INTEGER
*>          The number of superdiagonals of the matrix A if UPLO = 'U',
*>          or the number of subdiagonals if UPLO = 'L'.  KD >= 1.
*> \endverbatim
*>
*> \param[in,out] AB
*> \verbatim
*>          AB is DOUBLE PRECISION array, dimension (LDAB, N)
*>          On entry, the upper or lower triangle of the skew-symmetric band
*>          matrix A, stored in the first KD+1 rows of the array.  The
*>          j-th column of A is stored in the j-th column of the array AB
*>          as follows:
*>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<j;
*>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<i<=min(n,j+kd).
*>
*>          On exit, AB is overwritten by values generated during the
*>          reduction to tridiagonal form.  If UPLO = 'U', the first
*>          superdiagonal and the diagonal of the tridiagonal matrix T
*>          are returned in rows KD and KD+1 of AB, and if UPLO = 'L',
*>          the diagonal and first subdiagonal of T are returned in the
*>          first two rows of AB.
*> \endverbatim
*>
*> \param[in] LDAB
*> \verbatim
*>          LDAB is INTEGER
*>          The leading dimension of the array AB.  LDAB >= KD + 1.
*> \endverbatim
*>
*> \param[out] Q
*> \verbatim
*>          Q is DOUBLE PRECISION array, dimension (LDQ, N)
*>          If JOBZ = 'V', the N-by-N orthogonal matrix used in the
*>                         reduction to tridiagonal form.
*>          If JOBZ = 'N', the array Q is not referenced.
*> \endverbatim
*>
*> \param[in] LDQ
*> \verbatim
*>          LDQ is INTEGER
*>          The leading dimension of the array Q.  If JOBZ = 'V', then
*>          LDQ >= max(1,N).
*> \endverbatim
*>
*> \param[in] VL
*> \verbatim
*>          VL is DOUBLE PRECISION
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
*>          If RANGE='I', the index of eigenvalue with the
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
*>          by reducing AB to tridiagonal form.
*>
*>          Eigenvalues will be computed most accurately when ABSTOL is
*>          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
*>          If this routine returns with INFO>0, indicating that some
*>          eigenvectors did not converge, try setting ABSTOL to
*>          2*DLAMCH('S').
*>
*>          See "Computing Small Singular Values of Bidiagonal Matrices
*>          with Guaranteed High Relative Accuracy," by Demmel and
*>          Kahan, LAPACK Working Note #3.
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
*>          eigenvalues in block descending order. Refer to SKYEV. 
*> \endverbatim
*>
*> \param[out] Z
*> \verbatim
*>          Z is DOUBLE PRECISION array, dimension (LDZ, max(1,M))
*>          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
*>          contain the block form orthonormal eigenvectors of the
*>          matrix A corresponding to the selected eigenvalues, with the
*>          i-th column of Z holding the block form eigenvector associated
*>          with W(i). The block form eigenvectors (v1, v2) corresponding
*>          to a conjugate purely imaginary eigenvalue pair (ri, -ri) means
*>          A*v1 = r*v2 and A*v2 = -r*v1.
*>          If an eigenvector fails to converge, then that column of Z
*>          contains the latest approximation to the eigenvector, and the
*>          index of the eigenvector is returned in IFAIL.
*>          If JOBZ = 'N', then Z is not referenced.
*>          Note: the user must ensure that at least max(1,M) columns are
*>          supplied in the array Z; if RANGE = 'V', the exact value of M
*>          is not known in advance and an upper bound must be used.
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
*>          WORK is DOUBLE PRECISION array, dimension (LWORK)
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The length of the array WORK. LWORK >= 1, when N <= 1;
*>          otherwise
*>          If JOBZ = 'N' and N > 1, LWORK must be queried.
*>                                   LWORK = MAX(1, 7*N, dimension) where
*>                                   dimension = (2KD+1)*N + KD*NTHREADS + 2*N
*>                                   where KD is the size of the band.
*>                                   NTHREADS is the number of threads used when
*>                                   openMP compilation is enabled, otherwise =1.
*>          If JOBZ = 'V' and N > 1, LWORK must be queried. Not yet available
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (5*N)
*> \endverbatim
*>
*> \param[out] IFAIL
*> \verbatim
*>          IFAIL is INTEGER array, dimension (N)
*>          If JOBZ = 'V', then if INFO = 0, the first M elements of
*>          IFAIL are zero.  If INFO > 0, then IFAIL contains the
*>          indices of the eigenvectors that failed to converge.
*>          If JOBZ = 'N', then IFAIL is not referenced.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit.
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>          > 0:  if INFO = i, and i is:
*>                <= N: then i eigenvectors failed to converge in
*>                     DSKEWSTEIN; their indices are stored in IFAIL.
*>                > N: DSKEWSTEBZ returned INFO = INFO - N;
*>                     see DSKEWSTEBZ for details.
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
*> \ingroup hbevx_2stage
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  All details about the 2stage techniques are available in:
*>
*>  Azzam Haidar, Hatem Ltaief, and Jack Dongarra.
*>  Parallel reduction to condensed forms for symmetric eigenvalue problems
*>  using aggregated fine-grained and memory-aware kernels. In Proceedings
*>  of 2011 International Conference for High Performance Computing,
*>  Networking, Storage and Analysis (SC '11), New York, NY, USA,
*>  Article 8 , 11 pages.
*>  http://doi.acm.org/10.1145/2063384.2063394
*>
*>  A. Haidar, J. Kurzak, P. Luszczek, 2013.
*>  An improved parallel singular value algorithm and its implementation
*>  for multicore hardware, In Proceedings of 2013 International Conference
*>  for High Performance Computing, Networking, Storage and Analysis (SC '13).
*>  Denver, Colorado, USA, 2013.
*>  Article 90, 12 pages.
*>  http://doi.acm.org/10.1145/2503210.2503292
*>
*>  A. Haidar, R. Solca, S. Tomov, T. Schulthess and J. Dongarra.
*>  A novel hybrid CPU-GPU generalized eigensolver for electronic structure
*>  calculations based on fine-grained memory aware tasks.
*>  International Journal of High Performance Computing Applications.
*>  Volume 28 Issue 2, Pages 196-209, May 2014.
*>  http://hpc.sagepub.com/content/28/2/196
*>
*> \endverbatim
*
*  =====================================================================
      SUBROUTINE DSKEWSBEVX_2STAGE( JOBZ, RANGE, UPLO, N, KD, AB, LDAB,
     $                          Q,
     $                          LDQ, VL, VU, IL, IU, ABSTOL, M, W, Z,
     $                          LDZ, WORK, LWORK, IWORK, IFAIL, INFO )
*
      IMPLICIT NONE
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE, UPLO
      INTEGER            IL, INFO, IU, KD, LDAB, LDQ, LDZ, M, N, LWORK
      DOUBLE PRECISION   ABSTOL, VL, VU
*     ..
*     .. Array Arguments ..
      INTEGER            IFAIL( * ), IWORK( * )
      DOUBLE PRECISION   AB( LDAB, * ), Q( LDQ, * ), W( * ), WORK( * ),
     $                   Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALLEIG, INDEIG, LOWER, TEST, VALEIG, WANTZ,
     $                   LQUERY
      CHARACTER          ORDER
      INTEGER            I, IINFO, IMAX, INDE, INDEE, INDIBL,
     $                   INDISP, INDIWO, INDWRK, ISCALE, ITMP1, J, JJ,
     $                   LLWORK, LWMIN, LHTRD, LWTRD, IB, INDHOUS,
     $                   NSPLIT, K
      DOUBLE PRECISION   ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN,
     $                   SIGMA, SMLNUM, TMP1, VLL, VUU
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV2STAGE
      DOUBLE PRECISION   DLAMCH, DLANSKEWSB
      EXTERNAL           LSAME, DLAMCH, DLANSKEWSB, ILAENV2STAGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMV, DLACPY,
     $                   DSCAL,
     $                   DSKEWSTEBZ, DSKEWSTEIN, DSKEWSTEQR, DSWAP,
     $                   DSKEWSYTRD_SB2ST, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )
      LOWER = LSAME( UPLO, 'L' )
      LQUERY = ( LWORK.EQ.-1 )
*
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( KD.LT.1 ) THEN
         INFO = -5
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -7
      ELSE IF( WANTZ .AND. LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE
         IF( VALEIG ) THEN
            IF( N.GT.0 .AND. VU.LE.VL .OR. VU.LT.ZERO ) THEN
               INFO = -11
            END IF
         ELSE IF( INDEIG ) THEN
            IF( IL.LT.1 .OR. IL.GT.MAX( 1, (N+1)/2 ) ) THEN
               INFO = -12
            ELSE IF( IU.LT.MIN( (N+1)/2, IL ) .OR. IU.GT.(N+1)/2 )
     $      THEN
               INFO = -13
            END IF
         END IF
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) )
     $      INFO = -18
      END IF
*
      IF( INFO.EQ.0 ) THEN
         IF( N.LE.1 ) THEN
            LWMIN = 1
            WORK( 1 ) = LWMIN
         ELSE
            IB    = ILAENV2STAGE( 2, 'SSYTRD_SB2ST', JOBZ,
     $                            N, KD, -1, -1 )
            LHTRD = ILAENV2STAGE( 3, 'SSYTRD_SB2ST', JOBZ,
     $                            N, KD, IB, -1 )
            LWTRD = ILAENV2STAGE( 4, 'SSYTRD_SB2ST', JOBZ,
     $                            N, KD, IB, -1 )
            LWMIN = 2*N + LHTRD + LWTRD
            WORK( 1 )  = LWMIN
         ENDIF
*
         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY )
     $      INFO = -20
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSKEWSBEVX_2STAGE ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      M = 0
      IF( N.EQ.0 )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         IF( ALLEIG .OR. INDEIG ) THEN
            M = 1
            W( 1 ) = ZERO
         ELSE
            IF( VL.LT.ZERO ) THEN
               M = 1
               W( 1 ) = ZERO
            END IF
         END IF
         IF( WANTZ )
     $      Z( 1, 1 ) = ONE
         RETURN
      END IF
*
*     Get machine constants.
*
      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS    = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN   = SQRT( SMLNUM )
      RMAX   = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
*
*     Scale matrix to allowable range, if necessary.
*
      ISCALE = 0
      ABSTLL = ABSTOL
      IF( VALEIG ) THEN
         VLL = VL
         VUU = VU
      ELSE
         VLL = ZERO
         VUU = ZERO
      END IF
      ANRM = DLANSKEWSB( 'M', UPLO, N, KD, AB, LDAB, WORK )
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 ) THEN
         IF( LOWER ) THEN
            DO 10 J = 1, N-1
               CALL DSCAL( MIN( KD, N-J ), SIGMA, AB( 2, J ), 1 )
   10       CONTINUE
         ELSE
            DO 20 J = 2, N
               CALL DSCAL( MIN( KD, J-1 ), SIGMA,
     $                     AB( MAX(KD+2-J, 1), J ), 1 )
   20       CONTINUE
         END IF
         IF( ABSTOL.GT.0 )
     $      ABSTLL = ABSTOL*SIGMA
         IF( VALEIG ) THEN
            VLL = VL*SIGMA
            VUU = VU*SIGMA
         END IF
      END IF
*
*     Call DSKEWSYTRD_SB2ST to reduce skew-symmetric band matrix to tridiagonal form.
*
      INDE    = 1 + N
      INDHOUS = INDE + N
      INDWRK  = INDHOUS + LHTRD
      LLWORK  = LWORK - INDWRK + 1
*
      CALL DSKEWSYTRD_SB2ST( "N", JOBZ, UPLO, N, KD, AB, LDAB,
     $                    WORK( INDE ), WORK( INDHOUS ), LHTRD,
     $                    WORK( INDWRK ), LLWORK, IINFO )
      IF(.NOT.LOWER)
     $   CALL DSCAL(N-1, -ONE, WORK( INDE ), 1)
*
*     If all eigenvalues are desired and ABSTOL is less than or equal
*     to zero, then call SSTERF or DSKEWSTEQR.  If this fails for some
*     eigenvalue, then try DSKEWSTEBZ.
*
      TEST = .FALSE.
      IF (INDEIG) THEN
         IF (IL.EQ.1 .AND. IU.EQ.N) THEN
            TEST = .TRUE.
         END IF
      END IF
*      IF ((ALLEIG .OR. TEST) .AND. (ABSTOL.LE.ZERO)) THEN
      IF( .FALSE. ) THEN
         INDEE = INDWRK + 2*N
         IF( WANTZ ) THEN
            CALL DLACPY( 'A', N, N, Q, LDQ, Z, LDZ )
         END IF
         CALL DCOPY( N-1, WORK( INDE ), 1, WORK( INDEE ), 1 )
         CALL DSKEWSTEQR( JOBZ, N, WORK( INDEE ), Z, LDZ,
     $                WORK( INDWRK ), INFO )
         CALL DCOPY( N-1, WORK( INDEE ), 1, W, 1 )
         W( N ) = ZERO
         IF( INFO.EQ.0 ) THEN
            DO I = 1, N
               IFAIL( I ) = 0
            END DO
         END IF
         IF( INFO.EQ.0 ) THEN
            M = N
            GO TO 30
         END IF
         INFO = 0
      END IF
*
*     Otherwise, call DSKEWSTEBZ and, if eigenvectors are desired, DSKEWSTEIN.
*
      IF( WANTZ ) THEN
         ORDER = 'B'
      ELSE
         ORDER = 'E'
      END IF
      INDIBL = 1
      INDISP = INDIBL + N
      INDIWO = INDISP + N
      CALL DSKEWSTEBZ( RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL,
     $             WORK( INDE ), M, NSPLIT, W,
     $             IWORK( INDIBL ), IWORK( INDISP ), WORK( INDWRK ),
     $             IWORK( INDIWO ), IINFO )
      IF( IINFO.NE.0 ) THEN
         INFO = N + IINFO
         IF( IINFO.NE.1 )
     $      GO TO 30
      END IF
*
      IF( WANTZ ) THEN
         CALL DSKEWSTEIN( N, WORK( INDE ), M, W,
     $                IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ,
     $                WORK( INDWRK ), IWORK( INDIWO ), IFAIL, IINFO )
         IF( IINFO.NE.0 .AND. INFO.EQ.0 )
     $      INFO = IINFO
*
*        Apply orthogonal matrix used in reduction to tridiagonal
*        form to eigenvectors returned by DSKEWSTEIN.
*
         DO J = 1, M
            CALL DCOPY( N, Z( 1, J ), 1, WORK( 1 ), 1 )
            CALL DGEMV( 'N', N, N, ONE, Q, LDQ, WORK, 1, ZERO,
     $                  Z( 1, J ), 1 )
         END DO
      END IF
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
   30 CONTINUE
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = M
         ELSE
            IMAX = INFO - 1
         END IF
         CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF
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
*     Set WORK(1) to optimal workspace size.
*
      WORK( 1 ) = LWMIN
*
      RETURN
*
*     End of DSKEWSBEVX_2STAGE
*
      END
