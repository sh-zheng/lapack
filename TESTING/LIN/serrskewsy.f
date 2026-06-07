*> \brief \b SERRSKEWSY
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE SERRSKEWSY( PATH, NUNIT )
*
*       .. Scalar Arguments ..
*       CHARACTER*3        PATH
*       INTEGER            NUNIT
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SERRSKEWSY tests the error exits for the REAL routines
*> for skew-symmetric indefinite matrices.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] PATH
*> \verbatim
*>          PATH is CHARACTER*3
*>          The LAPACK path name for the routines to be tested.
*> \endverbatim
*>
*> \param[in] NUNIT
*> \verbatim
*>          NUNIT is INTEGER
*>          The unit number for output.
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
*> \ingroup single_lin
*
*  =====================================================================
      SUBROUTINE SERRSKEWSY( PATH, NUNIT )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER*3        PATH
      INTEGER            NUNIT
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NMAX
      PARAMETER          ( NMAX = 4 )
*     ..
*     .. Local Scalars ..
      CHARACTER*2        C2
      INTEGER            I, INFO, J
      REAL               ANRM, RCOND
*     ..
*     .. Local Arrays ..
      INTEGER            IP( NMAX ), IW( NMAX )
      REAL               A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ),
     $                   E( NMAX ), R1( NMAX ), R2( NMAX ), W( 3*NMAX ),
     $                   X( NMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      EXTERNAL           LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, SSKEWSYTRI2X, SSKEWSYTF2,
     $                   SSKEWSYTRF, SSKEWSYTRF_AA, SSKEWSYTRI, SSKEWSYTRS,
     $                   SSKEWSYTRS_AA, SSKEWSYTRI2
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      CHARACTER*32       SRNAMT
      INTEGER            INFOT, NOUT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL
*     ..
*     .. Executable Statements ..
*
      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
      C2 = PATH( 2: 3 )
*
*     Set the variables to innocuous values.
*
      DO 20 J = 1, NMAX
         DO 10 I = 1, NMAX
            A( I, J ) = 1. / REAL( I+J )
            AF( I, J ) = 1. / REAL( I+J )
   10    CONTINUE
         B( J ) = 0.E+0
         E( J ) = 0.E+0
         R1( J ) = 0.E+0
         R2( J ) = 0.E+0
         W( J ) = 0.E+0
         X( J ) = 0.E+0
         IP( J ) = J
         IW( J ) = J
   20 CONTINUE
      ANRM = 1.0
      RCOND = 1.0
      OK = .TRUE.
*
      IF( LSAMEN( 2, C2, 'KY' ) ) THEN
*
*        Test error exits of the routines that use factorization
*        of a skew-symmetric indefinite matrix with patrial
*        (Bunch-Kaufman) pivoting.
*
*        SSKEWSYTRF
*
         SRNAMT = 'SSKEWSYTRF'
         INFOT = 1
         CALL SSKEWSYTRF( '/', 0, A, 1, IP, W, 1, INFO )
         CALL CHKXER( 'SSKEWSYTRF', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SSKEWSYTRF( 'U', -1, A, 1, IP, W, 1, INFO )
         CALL CHKXER( 'SSKEWSYTRF', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL SSKEWSYTRF( 'U', 2, A, 1, IP, W, 4, INFO )
         CALL CHKXER( 'SSKEWSYTRF', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL SSKEWSYTRF( 'U', 0, A, 1, IP, W, 0, INFO )
         CALL CHKXER( 'SSKEWSYTRF', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL SSKEWSYTRF( 'U', 0, A, 1, IP, W, -2, INFO )
         CALL CHKXER( 'SSKEWSYTRF', INFOT, NOUT, LERR, OK )
*
*        SSKEWSYTF2
*
         SRNAMT = 'SSKEWSYTF2'
         INFOT = 1
         CALL SSKEWSYTF2( '/', 0, A, 1, IP, INFO )
         CALL CHKXER( 'SSKEWSYTF2', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SSKEWSYTF2( 'U', -1, A, 1, IP, INFO )
         CALL CHKXER( 'SSKEWSYTF2', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL SSKEWSYTF2( 'U', 2, A, 1, IP, INFO )
         CALL CHKXER( 'SSKEWSYTF2', INFOT, NOUT, LERR, OK )
*
*        SSKEWSYTRI
*
         SRNAMT = 'SSKEWSYTRI'
         INFOT = 1
         CALL SSKEWSYTRI( '/', 0, A, 1, IP, W, INFO )
         CALL CHKXER( 'SSKEWSYTRI', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SSKEWSYTRI( 'U', -1, A, 1, IP, W, INFO )
         CALL CHKXER( 'SSKEWSYTRI', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL SSKEWSYTRI( 'U', 2, A, 1, IP, W, INFO )
         CALL CHKXER( 'SSKEWSYTRI', INFOT, NOUT, LERR, OK )
*
*        SSKEWSYTRI2
*
         SRNAMT = 'SSKEWSYTRI2'
         INFOT = 1
         CALL SSKEWSYTRI2( '/', 0, A, 1, IP, W, IW(1), INFO )
         CALL CHKXER( 'SSKEWSYTRI2', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SSKEWSYTRI2( 'U', -1, A, 1, IP, W, IW(1), INFO )
         CALL CHKXER( 'SSKEWSYTRI2', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL SSKEWSYTRI2( 'U', 2, A, 1, IP, W, IW(1), INFO )
         CALL CHKXER( 'SSKEWSYTRI2', INFOT, NOUT, LERR, OK )
*
*        SSKEWSYTRI2X
*
         SRNAMT = 'SSKEWSYTRI2X'
         INFOT = 1
         CALL SSKEWSYTRI2X( '/', 0, A, 1, IP, W, 1, INFO )
         CALL CHKXER( 'SSKEWSYTRI2X', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SSKEWSYTRI2X( 'U', -1, A, 1, IP, W, 1, INFO )
         CALL CHKXER( 'SSKEWSYTRI2X', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL SSKEWSYTRI2X( 'U', 2, A, 1, IP, W, 1, INFO )
         CALL CHKXER( 'SSKEWSYTRI2X', INFOT, NOUT, LERR, OK )
*
*        SSKEWSYTRS
*
         SRNAMT = 'SSKEWSYTRS'
         INFOT = 1
         CALL SSKEWSYTRS( '/', 0, 0, A, 1, IP, B, 1, INFO )
         CALL CHKXER( 'SSKEWSYTRS', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SSKEWSYTRS( 'U', -1, 0, A, 1, IP, B, 1, INFO )
         CALL CHKXER( 'SSKEWSYTRS', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL SSKEWSYTRS( 'U', 0, -1, A, 1, IP, B, 1, INFO )
         CALL CHKXER( 'SSKEWSYTRS', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL SSKEWSYTRS( 'U', 2, 1, A, 1, IP, B, 2, INFO )
         CALL CHKXER( 'SSKEWSYTRS', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL SSKEWSYTRS( 'U', 2, 1, A, 2, IP, B, 1, INFO )
         CALL CHKXER( 'SSKEWSYTRS', INFOT, NOUT, LERR, OK )
*
      ELSE IF( LSAMEN( 2, C2, 'KA' ) ) THEN
*
*        Test error exits of the routines that use factorization
*        of a skew-symmetric indefinite matrix with Aasen's algorithm.
*
*        SSKEWSYTRF_AA
*
         SRNAMT = 'SSKEWSYTRF_AA'
         INFOT = 1
         CALL SSKEWSYTRF_AA( '/', 0, A, 1, IP, W, 1, INFO )
         CALL CHKXER( 'SSKEWSYTRF_AA', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SSKEWSYTRF_AA( 'U', -1, A, 1, IP, W, 1, INFO )
         CALL CHKXER( 'SSKEWSYTRF_AA', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL SSKEWSYTRF_AA( 'U', 2, A, 1, IP, W, 4, INFO )
         CALL CHKXER( 'SSKEWSYTRF_AA', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL SSKEWSYTRF_AA( 'U', 0, A, 1, IP, W, 0, INFO )
         CALL CHKXER( 'SSKEWSYTRF_AA', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL SSKEWSYTRF_AA( 'U', 0, A, 1, IP, W, -2, INFO )
         CALL CHKXER( 'SSKEWSYTRF_AA', INFOT, NOUT, LERR, OK )
*
*        SSKEWSYTRS_AA
*
         SRNAMT = 'SSKEWSYTRS_AA'
         INFOT = 1
         CALL SSKEWSYTRS_AA( '/', 0, 0, A, 1, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'SSKEWSYTRS_AA', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SSKEWSYTRS_AA( 'U', -1, 0, A, 1, IP, B, 1, W, 1,
     $                       INFO )
         CALL CHKXER( 'SSKEWSYTRS_AA', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL SSKEWSYTRS_AA( 'U', 0, -1, A, 1, IP, B, 1, W, 1,
     $                       INFO )
         CALL CHKXER( 'SSKEWSYTRS_AA', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL SSKEWSYTRS_AA( 'U', 2, 1, A, 1, IP, B, 2, W, 1, INFO )
         CALL CHKXER( 'SSKEWSYTRS_AA', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL SSKEWSYTRS_AA( 'U', 2, 1, A, 2, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'SSKEWSYTRS_AA', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL SSKEWSYTRS_AA( 'U', 0, 1, A, 2, IP, B, 1, W, 0, INFO )
         CALL CHKXER( 'SSKEWSYTRS_AA', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL SSKEWSYTRS_AA( 'U', 0, 1, A, 2, IP, B, 1, W, -2,
     $                       INFO )
         CALL CHKXER( 'SSKEWSYTRS_AA', INFOT, NOUT, LERR, OK )
      END IF
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
*     End of SERRSKEWSY
*
      END
