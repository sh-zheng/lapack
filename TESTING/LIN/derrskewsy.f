*> \brief \b DERRSKEWSY
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DERRSKEWSY( PATH, NUNIT )
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
*> DERRSKEWSY tests the error exits for the DOUBLE PRECISION routines
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
*> \ingroup double_lin
*
*  =====================================================================
      SUBROUTINE DERRSKEWSY( PATH, NUNIT )
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
      DOUBLE PRECISION   ANRM, RCOND
*     ..
*     .. Local Arrays ..
      INTEGER            IP( NMAX ), IW( NMAX )
      DOUBLE PRECISION   A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ),
     $                   E( NMAX ), R1( NMAX ), R2( NMAX ), W( 3*NMAX ),
     $                   X( NMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      EXTERNAL           LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, DSKEWSYTRI2X, DSKEWSYTF2,
     $                   DSKEWSYTRF, DSKEWSYTRF_AA, DSKEWSYTRI, DSKEWSYTRS,
     $                   DSKEWSYTRS_AA, DSKEWSYTRI2
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
      INTRINSIC          DBLE
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
            A( I, J ) = 1. / DBLE( I+J )
            AF( I, J ) = 1. / DBLE( I+J )
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
*        DSKEWSYTRF
*
         SRNAMT = 'DSKEWSYTRF'
         INFOT = 1
         CALL DSKEWSYTRF( '/', 0, A, 1, IP, W, 1, INFO )
         CALL CHKXER( 'DSKEWSYTRF', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DSKEWSYTRF( 'U', -1, A, 1, IP, W, 1, INFO )
         CALL CHKXER( 'DSKEWSYTRF', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DSKEWSYTRF( 'U', 2, A, 1, IP, W, 4, INFO )
         CALL CHKXER( 'DSKEWSYTRF', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DSKEWSYTRF( 'U', 0, A, 1, IP, W, 0, INFO )
         CALL CHKXER( 'DSKEWSYTRF', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DSKEWSYTRF( 'U', 0, A, 1, IP, W, -2, INFO )
         CALL CHKXER( 'DSKEWSYTRF', INFOT, NOUT, LERR, OK )
*
*        DSKEWSYTF2
*
         SRNAMT = 'DSKEWSYTF2'
         INFOT = 1
         CALL DSKEWSYTF2( '/', 0, A, 1, IP, INFO )
         CALL CHKXER( 'DSKEWSYTF2', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DSKEWSYTF2( 'U', -1, A, 1, IP, INFO )
         CALL CHKXER( 'DSKEWSYTF2', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DSKEWSYTF2( 'U', 2, A, 1, IP, INFO )
         CALL CHKXER( 'DSKEWSYTF2', INFOT, NOUT, LERR, OK )
*
*        DSKEWSYTRI
*
         SRNAMT = 'DSKEWSYTRI'
         INFOT = 1
         CALL DSKEWSYTRI( '/', 0, A, 1, IP, W, INFO )
         CALL CHKXER( 'DSKEWSYTRI', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DSKEWSYTRI( 'U', -1, A, 1, IP, W, INFO )
         CALL CHKXER( 'DSKEWSYTRI', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DSKEWSYTRI( 'U', 2, A, 1, IP, W, INFO )
         CALL CHKXER( 'DSKEWSYTRI', INFOT, NOUT, LERR, OK )
*
*        DSKEWSYTRI2
*
         SRNAMT = 'DSKEWSYTRI2'
         INFOT = 1
         CALL DSKEWSYTRI2( '/', 0, A, 1, IP, W, IW(1), INFO )
         CALL CHKXER( 'DSKEWSYTRI2', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DSKEWSYTRI2( 'U', -1, A, 1, IP, W, IW(1), INFO )
         CALL CHKXER( 'DSKEWSYTRI2', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DSKEWSYTRI2( 'U', 2, A, 1, IP, W, IW(1), INFO )
         CALL CHKXER( 'DSKEWSYTRI2', INFOT, NOUT, LERR, OK )
*
*        DSKEWSYTRI2X
*
         SRNAMT = 'DSKEWSYTRI2X'
         INFOT = 1
         CALL DSKEWSYTRI2X( '/', 0, A, 1, IP, W, 1, INFO )
         CALL CHKXER( 'DSKEWSYTRI2X', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DSKEWSYTRI2X( 'U', -1, A, 1, IP, W, 1, INFO )
         CALL CHKXER( 'DSKEWSYTRI2X', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DSKEWSYTRI2X( 'U', 2, A, 1, IP, W, 1, INFO )
         CALL CHKXER( 'DSKEWSYTRI2X', INFOT, NOUT, LERR, OK )
*
*        DSKEWSYTRS
*
         SRNAMT = 'DSKEWSYTRS'
         INFOT = 1
         CALL DSKEWSYTRS( '/', 0, 0, A, 1, IP, B, 1, INFO )
         CALL CHKXER( 'DSKEWSYTRS', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DSKEWSYTRS( 'U', -1, 0, A, 1, IP, B, 1, INFO )
         CALL CHKXER( 'DSKEWSYTRS', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DSKEWSYTRS( 'U', 0, -1, A, 1, IP, B, 1, INFO )
         CALL CHKXER( 'DSKEWSYTRS', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DSKEWSYTRS( 'U', 2, 1, A, 1, IP, B, 2, INFO )
         CALL CHKXER( 'DSKEWSYTRS', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DSKEWSYTRS( 'U', 2, 1, A, 2, IP, B, 1, INFO )
         CALL CHKXER( 'DSKEWSYTRS', INFOT, NOUT, LERR, OK )
*
      ELSE IF( LSAMEN( 2, C2, 'KA' ) ) THEN
*
*        Test error exits of the routines that use factorization
*        of a skew-symmetric indefinite matrix with Aasen's algorithm.
*
*        DSKEWSYTRF_AA
*
         SRNAMT = 'DSKEWSYTRF_AA'
         INFOT = 1
         CALL DSKEWSYTRF_AA( '/', 0, A, 1, IP, W, 1, INFO )
         CALL CHKXER( 'DSKEWSYTRF_AA', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DSKEWSYTRF_AA( 'U', -1, A, 1, IP, W, 1, INFO )
         CALL CHKXER( 'DSKEWSYTRF_AA', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DSKEWSYTRF_AA( 'U', 2, A, 1, IP, W, 4, INFO )
         CALL CHKXER( 'DSKEWSYTRF_AA', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DSKEWSYTRF_AA( 'U', 0, A, 1, IP, W, 0, INFO )
         CALL CHKXER( 'DSKEWSYTRF_AA', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DSKEWSYTRF_AA( 'U', 0, A, 1, IP, W, -2, INFO )
         CALL CHKXER( 'DSKEWSYTRF_AA', INFOT, NOUT, LERR, OK )
*
*        DSKEWSYTRS_AA
*
         SRNAMT = 'DSKEWSYTRS_AA'
         INFOT = 1
         CALL DSKEWSYTRS_AA( '/', 0, 0, A, 1, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'DSKEWSYTRS_AA', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DSKEWSYTRS_AA( 'U', -1, 0, A, 1, IP, B, 1, W, 1,
     $                       INFO )
         CALL CHKXER( 'DSKEWSYTRS_AA', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DSKEWSYTRS_AA( 'U', 0, -1, A, 1, IP, B, 1, W, 1,
     $                       INFO )
         CALL CHKXER( 'DSKEWSYTRS_AA', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DSKEWSYTRS_AA( 'U', 2, 1, A, 1, IP, B, 2, W, 1, INFO )
         CALL CHKXER( 'DSKEWSYTRS_AA', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DSKEWSYTRS_AA( 'U', 2, 1, A, 2, IP, B, 1, W, 1, INFO )
         CALL CHKXER( 'DSKEWSYTRS_AA', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL DSKEWSYTRS_AA( 'U', 0, 1, A, 2, IP, B, 1, W, 0, INFO )
         CALL CHKXER( 'DSKEWSYTRS_AA', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL DSKEWSYTRS_AA( 'U', 0, 1, A, 2, IP, B, 1, W, -2,
     $                       INFO )
         CALL CHKXER( 'DSKEWSYTRS_AA', INFOT, NOUT, LERR, OK )
      END IF
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
*     End of DERRSKEWSY
*
      END
