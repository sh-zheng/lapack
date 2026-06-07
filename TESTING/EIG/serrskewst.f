*> \brief \b SERRSKEWST
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE SERRSKEWST( PATH, NUNIT )
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
*> SERRSKEWST tests the error exits for SSKEWSYTRD, SSKEWSTEQR and SSKEWSYEV.
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
*> \ingroup single_eig
*
*  =====================================================================
      SUBROUTINE SERRSKEWST( PATH, NUNIT )
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
*     NMAX has to be at least 3 or LIW may be too small
*     .. Parameters ..
      INTEGER            NMAX, LIW, LW
      PARAMETER          ( NMAX = 3, LIW = 12*NMAX, LW = 20*NMAX )
*     ..
*     .. Local Scalars ..
      CHARACTER*2        C2
      INTEGER            I, INFO, J, M, N, NSPLIT, NT
*     ..
*     .. Local Arrays ..
      INTEGER            I1( NMAX ), I2( NMAX ), I3( NMAX ), IW( LIW )
      REAL               A( NMAX, NMAX ), C( NMAX, NMAX ), D( NMAX ),
     $                   E( NMAX ), Q( NMAX, NMAX ), R( NMAX ),
     $                   TAU( NMAX ), W( LW ), X( NMAX ),
     $                   Z( NMAX, NMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      EXTERNAL           LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           CHKXER, SSKEWSTEQR, SSKEWSYEV, SSKEWSTEV,
     $                   SSKEWSYTRD, SSKEWSTEBZ, SSKEWSTEIN,
     $                   SSKEWSYEVX, SSKEWSTEVX, SSKEWSTEDC,
     $                   SSKEWSTEVD, SSKEWSYEVD
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
   10    CONTINUE
   20 CONTINUE
      DO 30 J = 1, NMAX
         D( J ) = REAL( J )
         E( J ) = 0.0
         I1( J ) = J
         I2( J ) = J
         TAU( J ) = 1.
   30 CONTINUE
      OK = .TRUE.
      NT = 0
*
*     Test error exits for the KT path.
*
      IF( LSAMEN( 2, C2, 'KT' ) ) THEN
*
*        SSKEWSYTRD
*
         SRNAMT = 'SSKEWSYTRD'
         INFOT = 1
         CALL SSKEWSYTRD( '/', 0, A, 1, E, TAU, W, 1, INFO )
         CALL CHKXER( 'SSKEWSYTRD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SSKEWSYTRD( 'U', -1, A, 1, E, TAU, W, 1, INFO )
         CALL CHKXER( 'SSKEWSYTRD', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL SSKEWSYTRD( 'U', 2, A, 1, E, TAU, W, 1, INFO )
         CALL CHKXER( 'SSKEWSYTRD', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL SSKEWSYTRD( 'U', 0, A, 1, E, TAU, W, 0, INFO )
         CALL CHKXER( 'SSKEWSYTRD', INFOT, NOUT, LERR, OK )
         NT = NT + 4
*
*        SSKEWSTEBZ
*
         SRNAMT = 'SSKEWSTEBZ'
         INFOT = 1
         CALL SSKEWSTEBZ( '/', 'E', 0, 0.0, 1.0, 1, 0, 0.0, E, M,
     $                NSPLIT, X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'SSKEWSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SSKEWSTEBZ( 'A', '/', 0, 0.0, 0.0, 0, 0, 0.0, E, M,
     $                NSPLIT, X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'SSKEWSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL SSKEWSTEBZ( 'A', 'E', -1, 0.0, 0.0, 0, 0, 0.0, E, M,
     $                NSPLIT, X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'SSKEWSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL SSKEWSTEBZ( 'V', 'E', 0, 0.0, 0.0, 0, 0, 0.0, E, M,
     $                NSPLIT, X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'SSKEWSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL SSKEWSTEBZ( 'I', 'E', 0, 0.0, 0.0, 0, 0, 0.0, E, M,
     $                NSPLIT, X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'SSKEWSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL SSKEWSTEBZ( 'I', 'E', 1, 0.0, 0.0, 2, 1, 0.0, E, M,
     $                NSPLIT, X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'SSKEWSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL SSKEWSTEBZ( 'I', 'E', 1, 0.0, 0.0, 1, 0, 0.0, E, M,
     $                NSPLIT, X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'SSKEWSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL SSKEWSTEBZ( 'I', 'E', 1, 0.0, 0.0, 1, 2, 0.0, E, M,
     $                NSPLIT, X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'SSKEWSTEBZ', INFOT, NOUT, LERR, OK )
         NT = NT + 8
*
*        SSKEWSTEIN
*
         SRNAMT = 'SSKEWSTEIN'
         INFOT = 1
         CALL SSKEWSTEIN( -1, E, 0, X, I1, I2, Z, 1, W, IW, I3,
     $                    INFO )
         CALL CHKXER( 'SSKEWSTEIN', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL SSKEWSTEIN( 0, E, -1, X, I1, I2, Z, 1, W, IW, I3,
     $                    INFO )
         CALL CHKXER( 'SSKEWSTEIN', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL SSKEWSTEIN( 0, E, 1, X, I1, I2, Z, 1, W, IW, I3,
     $                    INFO )
         CALL CHKXER( 'SSKEWSTEIN', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL SSKEWSTEIN( 2, E, 0, X, I1, I2, Z, 1, W, IW, I3,
     $                    INFO )
         CALL CHKXER( 'SSKEWSTEIN', INFOT, NOUT, LERR, OK )
         NT = NT + 4
*
*        SSKEWSTEQR
*
         SRNAMT = 'SSKEWSTEQR'
         INFOT = 1
         CALL SSKEWSTEQR( '/', 0, E, Z, 1, W, INFO )
         CALL CHKXER( 'SSKEWSTEQR', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SSKEWSTEQR( 'N', -1, E, Z, 1, W, INFO )
         CALL CHKXER( 'SSKEWSTEQR', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL SSKEWSTEQR( 'V', 2, E, Z, 1, W, INFO )
         CALL CHKXER( 'SSKEWSTEQR', INFOT, NOUT, LERR, OK )
         NT = NT + 3
*
*        SSKEWSYEV
*
         SRNAMT = 'SSKEWSYEV '
         INFOT = 1
         CALL SSKEWSYEV( '/', 'U', 0, A, 1, X, W, 1, INFO )
         CALL CHKXER( 'SSKEWSYEV ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SSKEWSYEV( 'N', '/', 0, A, 1, X, W, 1, INFO )
         CALL CHKXER( 'SSKEWSYEV ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL SSKEWSYEV( 'N', 'U', -1, A, 1, X, W, 1, INFO )
         CALL CHKXER( 'SSKEWSYEV ', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL SSKEWSYEV( 'N', 'U', 2, A, 1, X, W, 3, INFO )
         CALL CHKXER( 'SSKEWSYEV ', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL SSKEWSYEV( 'N', 'U', 2, A, 2, X, W, 2, INFO )
         CALL CHKXER( 'SSKEWSYEV ', INFOT, NOUT, LERR, OK )
         NT = NT + 5
*
*        SSKEWSTEV
*
         SRNAMT = 'SSKEWSTEV '
         INFOT = 1
         CALL SSKEWSTEV( '/', 0, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'SSKEWSTEV ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SSKEWSTEV( 'N', -1, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'SSKEWSTEV ', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL SSKEWSTEV( 'V', 2, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'SSKEWSTEV ', INFOT, NOUT, LERR, OK )
         NT = NT + 3
*
*        SSKEWSYEVX
*
         SRNAMT = 'SSKEWSYEVX'
         INFOT = 1
         CALL SSKEWSYEVX( '/', 'A', 'U', 0, A, 1, 0.0, 0.0, 0, 0,
     $                0.0, M, X, Z, 1, W, 1, IW, I3, INFO )
         CALL CHKXER( 'SSKEWSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SSKEWSYEVX( 'N', '/', 'U', 0, A, 1, 0.0, 1.0, 1, 0,
     $                0.0, M, X, Z, 1, W, 1, IW, I3, INFO )
         CALL CHKXER( 'SSKEWSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL SSKEWSYEVX( 'N', 'A', '/', 0, A, 1, 0.0, 0.0, 0, 0,
     $                0.0, M, X, Z, 1, W, 1, IW, I3, INFO )
         INFOT = 4
         CALL SSKEWSYEVX( 'N', 'A', 'U', -1, A, 1, 0.0, 0.0, 0,
     $                0, 0.0, M, X, Z, 1, W, 1, IW, I3, INFO )
         CALL CHKXER( 'SSKEWSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL SSKEWSYEVX( 'N', 'A', 'U', 2, A, 1, 0.0, 0.0, 0, 0,
     $                0.0, M, X, Z, 1, W, 16, IW, I3, INFO )
         CALL CHKXER( 'SSKEWSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL SSKEWSYEVX( 'N', 'V', 'U', 1, A, 1, 0.0, 0.0, 0, 0,
     $                0.0, M, X, Z, 1, W, 8, IW, I3, INFO )
         CALL CHKXER( 'SSKEWSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL SSKEWSYEVX( 'N', 'V', 'U', 1, A, 1, -2.0, -1.0, 0, 0,
     $                0.0, M, X, Z, 1, W, 8, IW, I3, INFO )
         CALL CHKXER( 'SSKEWSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL SSKEWSYEVX( 'N', 'I', 'U', 1, A, 1, 0.0, 0.0, 0, 0,
     $                0.0, M, X, Z, 1, W, 8, IW, I3, INFO )
         CALL CHKXER( 'SSKEWSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL SSKEWSYEVX( 'N', 'I', 'U', 1, A, 1, 0.0, 0.0, 2, 1,
     $                0.0, M, X, Z, 1, W, 8, IW, I3, INFO )
         CALL CHKXER( 'SSKEWSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL SSKEWSYEVX( 'N', 'I', 'U', 3, A, 3, 0.0, 0.0, 2, 1,
     $                0.0, M, X, Z, 1, W, 16, IW, I3, INFO )
         CALL CHKXER( 'SSKEWSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL SSKEWSYEVX( 'N', 'I', 'U', 2, A, 2, 0.0, 0.0, 1, 2,
     $                0.0, M, X, Z, 1, W, 8, IW, I3, INFO )
         CALL CHKXER( 'SSKEWSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 15
         CALL SSKEWSYEVX( 'V', 'A', 'U', 2, A, 2, 0.0, 0.0, 0, 0,
     $                0.0, M, X, Z, 1, W, 16, IW, I3, INFO )
         CALL CHKXER( 'SSKEWSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 17
         CALL SSKEWSYEVX( 'V', 'A', 'U', 1, A, 1, 0.0, 0.0, 0, 0,
     $                0.0, M, X, Z, 1, W, 0, IW, I3, INFO )
         CALL CHKXER( 'SSKEWSYEVX', INFOT, NOUT, LERR, OK )
         NT = NT + 13
*
*        SSKEWSTEVX
*
         SRNAMT = 'SSKEWSTEVX'
         INFOT = 1
         CALL SSKEWSTEVX( '/', 'A', 0, E, 0.0, 0.0, 0, 0, 0.0,
     $                M, X, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'SSKEWSTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SSKEWSTEVX( 'N', '/', 0, E, 0.0, 1.0, 1, 0, 0.0,
     $                M, X, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'SSKEWSTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL SSKEWSTEVX( 'N', 'A', -1, E, 0.0, 0.0, 0, 0, 0.0,
     $                M, X, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'SSKEWSTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL SSKEWSTEVX( 'N', 'V', 1, E, 0.0, 0.0, 0, 0, 0.0,
     $                M, X, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'SSKEWSTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL SSKEWSTEVX( 'N', 'V', 1, E, -2.0, -1.0, 0, 0, 0.0,
     $                M, X, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'SSKEWSTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL SSKEWSTEVX( 'N', 'I', 1, E, 0.0, 0.0, 0, 0, 0.0,
     $                M, X, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'SSKEWSTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL SSKEWSTEVX( 'N', 'I', 1, E, 0.0, 0.0, 2, 1, 0.0,
     $                M, X, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'SSKEWSTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL SSKEWSTEVX( 'N', 'I', 3, E, 0.0, 0.0, 2, 1, 0.0,
     $                M, X, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'SSKEWSTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL SSKEWSTEVX( 'N', 'I', 2, E, 0.0, 0.0, 1, 2, 0.0,
     $                M, X, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'SSKEWSTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL SSKEWSTEVX( 'V', 'A', 2, E, 0.0, 0.0, 0, 0, 0.0,
     $                M, X, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'SSKEWSTEVX', INFOT, NOUT, LERR, OK )
         NT = NT + 10
*
*        SSKEWSTEDC
*
         SRNAMT = 'SSKEWSTEDC'
         INFOT = 1
         CALL SSKEWSTEDC( '/', 0, E, Z, 1, W, 1, IW, 1, INFO )
         CALL CHKXER( 'SSKEWSTEDC', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SSKEWSTEDC( 'N', -1, E, Z, 1, W, 1, IW, 1, INFO )
         CALL CHKXER( 'SSKEWSTEDC', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL SSKEWSTEDC( 'V', 2, E, Z, 1, W, 23, IW, 28, INFO )
         CALL CHKXER( 'SSKEWSTEDC', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL SSKEWSTEDC( 'N', 1, E, Z, 1, W, 0, IW, 1, INFO )
         CALL CHKXER( 'SSKEWSTEDC', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL SSKEWSTEDC( 'I', 2, E, Z, 2, W, 0, IW, 12, INFO )
         CALL CHKXER( 'SSKEWSTEDC', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL SSKEWSTEDC( 'V', 2, E, Z, 2, W, 0, IW, 28, INFO )
         CALL CHKXER( 'SSKEWSTEDC', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL SSKEWSTEDC( 'N', 1, E, Z, 1, W, 1, IW, 0, INFO )
         CALL CHKXER( 'SSKEWSTEDC', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL SSKEWSTEDC( 'I', 2, E, Z, 2, W, 19, IW, 0, INFO )
         CALL CHKXER( 'SSKEWSTEDC', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL SSKEWSTEDC( 'V', 2, E, Z, 2, W, 23, IW, 0, INFO )
         CALL CHKXER( 'SSKEWSTEDC', INFOT, NOUT, LERR, OK )
         NT = NT + 9
*
*        SSKEWSTEVD
*
         SRNAMT = 'SSKEWSTEVD'
         INFOT = 1
         CALL SSKEWSTEVD( '/', 0, D, E, Z, 1, W, 1, IW, 1, INFO )
         CALL CHKXER( 'SSKEWSTEVD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SSKEWSTEVD( 'N', -1, D, E, Z, 1, W, 1, IW, 1, INFO )
         CALL CHKXER( 'SSKEWSTEVD', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL SSKEWSTEVD( 'V', 2, D, E, Z, 1, W, 19, IW, 12, INFO )
         CALL CHKXER( 'SSKEWSTEVD', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL SSKEWSTEVD( 'N', 1, D, E, Z, 1, W, 0, IW, 1, INFO )
         CALL CHKXER( 'SSKEWSTEVD', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL SSKEWSTEVD( 'V', 2, D, E, Z, 2, W, 4, IW, 12, INFO )
         CALL CHKXER( 'SSKEWSTEVD', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL SSKEWSTEVD( 'N', 0, D, E, Z, 1, W, 1, IW, 0, INFO )
         CALL CHKXER( 'SSKEWSTEVD', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL SSKEWSTEVD( 'V', 2, D, E, Z, 2, W, 19, IW, 7, INFO )
         CALL CHKXER( 'SSKEWSTEVD', INFOT, NOUT, LERR, OK )
         NT = NT + 7
*
*        SSKEWSYEVD
*
         SRNAMT = 'SSKEWSYEVD'
         INFOT = 1
         CALL SSKEWSYEVD( '/', 'U', 0, A, 1, X, W, 1, IW, 1, INFO )
         CALL CHKXER( 'SSKEWSYEVD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SSKEWSYEVD( 'N', '/', 0, A, 1, X, W, 1, IW, 1, INFO )
         CALL CHKXER( 'SSKEWSYEVD', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL SSKEWSYEVD( 'N', 'U', -1, A, 1, X, W, 1, IW, 1, INFO )
         CALL CHKXER( 'SSKEWSYEVD', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL SSKEWSYEVD( 'N', 'U', 2, A, 1, X, W, 3, IW, 1, INFO )
         CALL CHKXER( 'SSKEWSYEVD', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL SSKEWSYEVD( 'N', 'U', 1, A, 1, X, W, 0, IW, 1, INFO )
         CALL CHKXER( 'SSKEWSYEVD', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL SSKEWSYEVD( 'N', 'U', 2, A, 2, X, W, 0, IW, 1, INFO )
         CALL CHKXER( 'SSKEWSYEVD', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL SSKEWSYEVD( 'V', 'U', 2, A, 2, X, W, 10, IW, 12,
     $                    INFO )
         CALL CHKXER( 'SSKEWSYEVD', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL SSKEWSYEVD( 'N', 'U', 1, A, 1, X, W, 1, IW, 0, INFO )
         CALL CHKXER( 'SSKEWSYEVD', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL SSKEWSYEVD( 'N', 'U', 2, A, 2, X, W, 5, IW, 0, INFO )
         CALL CHKXER( 'SSKEWSYEVD', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL SSKEWSYEVD( 'V', 'U', 2, A, 2, X, W, 27, IW, 7, INFO )
         CALL CHKXER( 'SSKEWSYEVD', INFOT, NOUT, LERR, OK )
         NT = NT + 10
      END IF
*
*     Print a summary line.
*
      IF( OK ) THEN
         WRITE( NOUT, FMT = 9999 )PATH, NT
      ELSE
         WRITE( NOUT, FMT = 9998 )PATH
      END IF
*
 9999 FORMAT( 1X, A3, ' routines passed the tests of the error exits',
     $      ' (', I3, ' tests done)' )
 9998 FORMAT( ' *** ', A3, ' routines failed the tests of the error ',
     $      'exits ***' )
*
      RETURN
*
*     End of SERRSKEWST
*
      END
