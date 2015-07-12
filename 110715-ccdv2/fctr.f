      SUBROUTINE FACTR(N,FCTR)
C
C COMPUTES FACTORIAL OF N
C
      IMPLICIT REAL*8(A-H,O-Z)
      IF(N.EQ.0)GO TO 11
      FCTR=1.
      DO 10 I=1,N
   10 FCTR=FLOAT(I)*FCTR
      GO TO 12
   11 FCTR=1.
   12 RETURN
      END
