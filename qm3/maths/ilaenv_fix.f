C ---------------------------------------------------------------------
C     stupid fix for C binding...
C     don't know why character varaibles NAME and OPTS are segmenting
C     when called from C
C
      INTEGER          FUNCTION XILAENV( N1 )
      INTEGER N1
      XILAENV = ILAENV( 1, "DGETRI", "", N1, -1, -1, -1 )
      END
C ---------------------------------------------------------------------
