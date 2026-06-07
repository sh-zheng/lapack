*> \brief \b DERRSKEWST
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DERRSKEWST( PATH, NUNIT )
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
*> DERRSKEWST tests the error exits for DSKEWSYTRD, DSKEWSTEQR and DSKEWSYEV.
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
*> \ingroup double_eig
*
*  =====================================================================
      SUBROUTINE DERRSKEWST( PATH, NUNIT )
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
      DOUBLE PRECISION   A( NMAX, NMAX ), C( NMAX, NMAX ), D( NMAX ),
     $                   E( NMAX ), Q( NMAX, NMAX ), R( NMAX ),
     $                   TAU( NMAX ), W( LW ), X( NMAX ),
     $                   Z( NMAX, NMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      EXTERNAL           LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           CHKXER, DSKEWSTEQR, DSKEWSYEV, DSKEWSTEV,
     $                   DSKEWSYTRD, DSKEWSTEBZ, DSKEWSTEIN,
     $                   DSKEWSYEVX, DSKEWSTEVX, DSKEWSTEDC,
     $                   DSKEWSTEVD, DSKEWSYEVD
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
   10    CONTINUE
   20 CONTINUE
      DO 30 J = 1, NMAX
         D( J ) = DBLE( J )
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
*        DSKEWSYTRD
*
         SRNAMT = 'DSKEWSYTRD'
         INFOT = 1
         CALL DSKEWSYTRD( '/', 0, A, 1, E, TAU, W, 1, INFO )
         CALL CHKXER( 'DSKEWSYTRD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DSKEWSYTRD( 'U', -1, A, 1, E, TAU, W, 1, INFO )
         CALL CHKXER( 'DSKEWSYTRD', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DSKEWSYTRD( 'U', 2, A, 1, E, TAU, W, 1, INFO )
         CALL CHKXER( 'DSKEWSYTRD', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DSKEWSYTRD( 'U', 0, A, 1, E, TAU, W, 0, INFO )
         CALL CHKXER( 'DSKEWSYTRD', INFOT, NOUT, LERR, OK )
         NT = NT + 4
*
*        DSKEWSTEBZ
*
         SRNAMT = 'DSKEWSTEBZ'
         INFOT = 1
         CALL DSKEWSTEBZ( '/', 'E', 0, 0.0, 1.0, 1, 0, 0.0, E, M,
     $                NSPLIT, X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'DSKEWSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DSKEWSTEBZ( 'A', '/', 0, 0.0, 0.0, 0, 0, 0.0, E, M,
     $                NSPLIT, X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'DSKEWSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DSKEWSTEBZ( 'A', 'E', -1, 0.0, 0.0, 0, 0, 0.0, E, M,
     $                NSPLIT, X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'DSKEWSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DSKEWSTEBZ( 'V', 'E', 0, 0.0, 0.0, 0, 0, 0.0, E, M,
     $                NSPLIT, X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'DSKEWSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DSKEWSTEBZ( 'I', 'E', 0, 0.0, 0.0, 0, 0, 0.0, E, M,
     $                NSPLIT, X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'DSKEWSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DSKEWSTEBZ( 'I', 'E', 1, 0.0, 0.0, 2, 1, 0.0, E, M,
     $                NSPLIT, X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'DSKEWSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DSKEWSTEBZ( 'I', 'E', 1, 0.0, 0.0, 1, 0, 0.0, E, M,
     $                NSPLIT, X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'DSKEWSTEBZ', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DSKEWSTEBZ( 'I', 'E', 1, 0.0, 0.0, 1, 2, 0.0, E, M,
     $                NSPLIT, X, I1, I2, W, IW, INFO )
         CALL CHKXER( 'DSKEWSTEBZ', INFOT, NOUT, LERR, OK )
         NT = NT + 8
*
*        DSKEWSTEIN
*
         SRNAMT = 'DSKEWSTEIN'
         INFOT = 1
         CALL DSKEWSTEIN( -1, E, 0, X, I1, I2, Z, 1, W, IW, I3,
     $                    INFO )
         CALL CHKXER( 'DSKEWSTEIN', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DSKEWSTEIN( 0, E, -1, X, I1, I2, Z, 1, W, IW, I3,
     $                    INFO )
         CALL CHKXER( 'DSKEWSTEIN', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DSKEWSTEIN( 0, E, 1, X, I1, I2, Z, 1, W, IW, I3,
     $                    INFO )
         CALL CHKXER( 'DSKEWSTEIN', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DSKEWSTEIN( 2, E, 0, X, I1, I2, Z, 1, W, IW, I3,
     $                    INFO )
         CALL CHKXER( 'DSKEWSTEIN', INFOT, NOUT, LERR, OK )
         NT = NT + 4
*
*        DSKEWSTEQR
*
         SRNAMT = 'DSKEWSTEQR'
         INFOT = 1
         CALL DSKEWSTEQR( '/', 0, E, Z, 1, W, INFO )
         CALL CHKXER( 'DSKEWSTEQR', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DSKEWSTEQR( 'N', -1, E, Z, 1, W, INFO )
         CALL CHKXER( 'DSKEWSTEQR', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DSKEWSTEQR( 'V', 2, E, Z, 1, W, INFO )
         CALL CHKXER( 'DSKEWSTEQR', INFOT, NOUT, LERR, OK )
         NT = NT + 3
*
*        DSKEWSYEV
*
         SRNAMT = 'DSKEWSYEV '
         INFOT = 1
         CALL DSKEWSYEV( '/', 'U', 0, A, 1, X, W, 1, INFO )
         CALL CHKXER( 'DSKEWSYEV ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DSKEWSYEV( 'N', '/', 0, A, 1, X, W, 1, INFO )
         CALL CHKXER( 'DSKEWSYEV ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DSKEWSYEV( 'N', 'U', -1, A, 1, X, W, 1, INFO )
         CALL CHKXER( 'DSKEWSYEV ', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DSKEWSYEV( 'N', 'U', 2, A, 1, X, W, 3, INFO )
         CALL CHKXER( 'DSKEWSYEV ', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DSKEWSYEV( 'N', 'U', 2, A, 2, X, W, 2, INFO )
         CALL CHKXER( 'DSKEWSYEV ', INFOT, NOUT, LERR, OK )
         NT = NT + 5
*
*        DSKEWSTEV
*
         SRNAMT = 'DSKEWSTEV '
         INFOT = 1
         CALL DSKEWSTEV( '/', 0, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'DSKEWSTEV ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DSKEWSTEV( 'N', -1, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'DSKEWSTEV ', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DSKEWSTEV( 'V', 2, D, E, Z, 1, W, INFO )
         CALL CHKXER( 'DSKEWSTEV ', INFOT, NOUT, LERR, OK )
         NT = NT + 3
*
*        DSKEWSYEVX
*
         SRNAMT = 'DSKEWSYEVX'
         INFOT = 1
         CALL DSKEWSYEVX( '/', 'A', 'U', 0, A, 1, 0.0, 0.0, 0, 0,
     $                0.0, M, X, Z, 1, W, 1, IW, I3, INFO )
         CALL CHKXER( 'DSKEWSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DSKEWSYEVX( 'N', '/', 'U', 0, A, 1, 0.0, 1.0, 1, 0,
     $                0.0, M, X, Z, 1, W, 1, IW, I3, INFO )
         CALL CHKXER( 'DSKEWSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DSKEWSYEVX( 'N', 'A', '/', 0, A, 1, 0.0, 0.0, 0, 0,
     $                0.0, M, X, Z, 1, W, 1, IW, I3, INFO )
         INFOT = 4
         CALL DSKEWSYEVX( 'N', 'A', 'U', -1, A, 1, 0.0, 0.0, 0,
     $                0, 0.0, M, X, Z, 1, W, 1, IW, I3, INFO )
         CALL CHKXER( 'DSKEWSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DSKEWSYEVX( 'N', 'A', 'U', 2, A, 1, 0.0, 0.0, 0, 0,
     $                0.0, M, X, Z, 1, W, 16, IW, I3, INFO )
         CALL CHKXER( 'DSKEWSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DSKEWSYEVX( 'N', 'V', 'U', 1, A, 1, 0.0, 0.0, 0, 0,
     $                0.0, M, X, Z, 1, W, 8, IW, I3, INFO )
         CALL CHKXER( 'DSKEWSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DSKEWSYEVX( 'N', 'V', 'U', 1, A, 1, -2.0, -1.0, 0, 0,
     $                0.0, M, X, Z, 1, W, 8, IW, I3, INFO )
         CALL CHKXER( 'DSKEWSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL DSKEWSYEVX( 'N', 'I', 'U', 1, A, 1, 0.0, 0.0, 0, 0,
     $                0.0, M, X, Z, 1, W, 8, IW, I3, INFO )
         CALL CHKXER( 'DSKEWSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL DSKEWSYEVX( 'N', 'I', 'U', 1, A, 1, 0.0, 0.0, 2, 1,
     $                0.0, M, X, Z, 1, W, 8, IW, I3, INFO )
         CALL CHKXER( 'DSKEWSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL DSKEWSYEVX( 'N', 'I', 'U', 3, A, 3, 0.0, 0.0, 2, 1,
     $                0.0, M, X, Z, 1, W, 16, IW, I3, INFO )
         CALL CHKXER( 'DSKEWSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL DSKEWSYEVX( 'N', 'I', 'U', 2, A, 2, 0.0, 0.0, 1, 2,
     $                0.0, M, X, Z, 1, W, 8, IW, I3, INFO )
         CALL CHKXER( 'DSKEWSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 15
         CALL DSKEWSYEVX( 'V', 'A', 'U', 2, A, 2, 0.0, 0.0, 0, 0,
     $                0.0, M, X, Z, 1, W, 16, IW, I3, INFO )
         CALL CHKXER( 'DSKEWSYEVX', INFOT, NOUT, LERR, OK )
         INFOT = 17
         CALL DSKEWSYEVX( 'V', 'A', 'U', 1, A, 1, 0.0, 0.0, 0, 0,
     $                0.0, M, X, Z, 1, W, 0, IW, I3, INFO )
         CALL CHKXER( 'DSKEWSYEVX', INFOT, NOUT, LERR, OK )
         NT = NT + 13
*
*        DSKEWSTEVX
*
         SRNAMT = 'DSKEWSTEVX'
         INFOT = 1
         CALL DSKEWSTEVX( '/', 'A', 0, E, 0.0, 0.0, 0, 0, 0.0,
     $                M, X, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'DSKEWSTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DSKEWSTEVX( 'N', '/', 0, E, 0.0, 1.0, 1, 0, 0.0,
     $                M, X, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'DSKEWSTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DSKEWSTEVX( 'N', 'A', -1, E, 0.0, 0.0, 0, 0, 0.0,
     $                M, X, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'DSKEWSTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DSKEWSTEVX( 'N', 'V', 1, E, 0.0, 0.0, 0, 0, 0.0,
     $                M, X, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'DSKEWSTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DSKEWSTEVX( 'N', 'V', 1, E, -2.0, -1.0, 0, 0, 0.0,
     $                M, X, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'DSKEWSTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DSKEWSTEVX( 'N', 'I', 1, E, 0.0, 0.0, 0, 0, 0.0,
     $                M, X, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'DSKEWSTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DSKEWSTEVX( 'N', 'I', 1, E, 0.0, 0.0, 2, 1, 0.0,
     $                M, X, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'DSKEWSTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DSKEWSTEVX( 'N', 'I', 3, E, 0.0, 0.0, 2, 1, 0.0,
     $                M, X, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'DSKEWSTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DSKEWSTEVX( 'N', 'I', 2, E, 0.0, 0.0, 1, 2, 0.0,
     $                M, X, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'DSKEWSTEVX', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL DSKEWSTEVX( 'V', 'A', 2, E, 0.0, 0.0, 0, 0, 0.0,
     $                M, X, Z, 1, W, IW, I3, INFO )
         CALL CHKXER( 'DSKEWSTEVX', INFOT, NOUT, LERR, OK )
         NT = NT + 10
*
*        DSKEWSTEDC
*
         SRNAMT = 'DSKEWSTEDC'
         INFOT = 1
         CALL DSKEWSTEDC( '/', 0, E, Z, 1, W, 1, IW, 1, INFO )
         CALL CHKXER( 'DSKEWSTEDC', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DSKEWSTEDC( 'N', -1, E, Z, 1, W, 1, IW, 1, INFO )
         CALL CHKXER( 'DSKEWSTEDC', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DSKEWSTEDC( 'V', 2, E, Z, 1, W, 23, IW, 28, INFO )
         CALL CHKXER( 'DSKEWSTEDC', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DSKEWSTEDC( 'N', 1, E, Z, 1, W, 0, IW, 1, INFO )
         CALL CHKXER( 'DSKEWSTEDC', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DSKEWSTEDC( 'I', 2, E, Z, 2, W, 0, IW, 12, INFO )
         CALL CHKXER( 'DSKEWSTEDC', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DSKEWSTEDC( 'V', 2, E, Z, 2, W, 0, IW, 28, INFO )
         CALL CHKXER( 'DSKEWSTEDC', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL DSKEWSTEDC( 'N', 1, E, Z, 1, W, 1, IW, 0, INFO )
         CALL CHKXER( 'DSKEWSTEDC', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL DSKEWSTEDC( 'I', 2, E, Z, 2, W, 19, IW, 0, INFO )
         CALL CHKXER( 'DSKEWSTEDC', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL DSKEWSTEDC( 'V', 2, E, Z, 2, W, 23, IW, 0, INFO )
         CALL CHKXER( 'DSKEWSTEDC', INFOT, NOUT, LERR, OK )
         NT = NT + 9
*
*        DSKEWSTEVD
*
         SRNAMT = 'DSKEWSTEVD'
         INFOT = 1
         CALL DSKEWSTEVD( '/', 0, D, E, Z, 1, W, 1, IW, 1, INFO )
         CALL CHKXER( 'DSKEWSTEVD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DSKEWSTEVD( 'N', -1, D, E, Z, 1, W, 1, IW, 1, INFO )
         CALL CHKXER( 'DSKEWSTEVD', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DSKEWSTEVD( 'V', 2, D, E, Z, 1, W, 19, IW, 12, INFO )
         CALL CHKXER( 'DSKEWSTEVD', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DSKEWSTEVD( 'N', 1, D, E, Z, 1, W, 0, IW, 1, INFO )
         CALL CHKXER( 'DSKEWSTEVD', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DSKEWSTEVD( 'V', 2, D, E, Z, 2, W, 4, IW, 12, INFO )
         CALL CHKXER( 'DSKEWSTEVD', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL DSKEWSTEVD( 'N', 0, D, E, Z, 1, W, 1, IW, 0, INFO )
         CALL CHKXER( 'DSKEWSTEVD', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL DSKEWSTEVD( 'V', 2, D, E, Z, 2, W, 19, IW, 7, INFO )
         CALL CHKXER( 'DSKEWSTEVD', INFOT, NOUT, LERR, OK )
         NT = NT + 7
*
*        DSKEWSYEVD
*
         SRNAMT = 'DSKEWSYEVD'
         INFOT = 1
         CALL DSKEWSYEVD( '/', 'U', 0, A, 1, X, W, 1, IW, 1, INFO )
         CALL CHKXER( 'DSKEWSYEVD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DSKEWSYEVD( 'N', '/', 0, A, 1, X, W, 1, IW, 1, INFO )
         CALL CHKXER( 'DSKEWSYEVD', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DSKEWSYEVD( 'N', 'U', -1, A, 1, X, W, 1, IW, 1, INFO )
         CALL CHKXER( 'DSKEWSYEVD', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DSKEWSYEVD( 'N', 'U', 2, A, 1, X, W, 3, IW, 1, INFO )
         CALL CHKXER( 'DSKEWSYEVD', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DSKEWSYEVD( 'N', 'U', 1, A, 1, X, W, 0, IW, 1, INFO )
         CALL CHKXER( 'DSKEWSYEVD', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DSKEWSYEVD( 'N', 'U', 2, A, 2, X, W, 0, IW, 1, INFO )
         CALL CHKXER( 'DSKEWSYEVD', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DSKEWSYEVD( 'V', 'U', 2, A, 2, X, W, 10, IW, 12,
     $                    INFO )
         CALL CHKXER( 'DSKEWSYEVD', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL DSKEWSYEVD( 'N', 'U', 1, A, 1, X, W, 1, IW, 0, INFO )
         CALL CHKXER( 'DSKEWSYEVD', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL DSKEWSYEVD( 'N', 'U', 2, A, 2, X, W, 5, IW, 0, INFO )
         CALL CHKXER( 'DSKEWSYEVD', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL DSKEWSYEVD( 'V', 'U', 2, A, 2, X, W, 27, IW, 7, INFO )
         CALL CHKXER( 'DSKEWSYEVD', INFOT, NOUT, LERR, OK )
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
*     End of DERRSKEWST
*
      END
