      SUBROUTINE KIN(OMEGA,FATHER,GROUP,MOTHER,MAXPEO,MOMEGA,NPTOT)
C
C     THIS SUBROUTINE COMPUTES TWICE THE KINSHIP MATRIX FOR A PEDIGREE
C     AND STORES IT IN THE LOWER TRIANGLE OF OMEGA.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION OMEGA(MOMEGA,MOMEGA)
      INTEGER FATHER(MAXPEO),GROUP(MAXPEO),MOTHER(MAXPEO),COUNT
C
C     INITIALIZE OMEGA.  THE DIAGONAL ENTRIES OF OMEGA WILL BE USED TO KEEP
C     TRACK OF WHO IS CURRENTLY INCLUDED IN THE RECURSIVE CALCULATIONS.
C
      DO 10 I=1,NPTOT
      DO 10 J=1,I
 10   OMEGA(I,J)=0.0D0
      COUNT=0
C
C     FOR EACH FOUNDER I SET OMEGA(I,I) TO 1.
C
      DO 20 I=1,NPTOT
      IF (FATHER(I).EQ.0) THEN
      COUNT=COUNT+1
      OMEGA(I,I)=1.0D0
      END IF
 20   CONTINUE
C
C     FILL IN OMEGA FOR THE NON-FOUNDERS.
C
 70   DO 30 I=1,NPTOT
      IF (OMEGA(I,I).NE.0.0D0) GO TO 30
      I1=FATHER(I)
      I2=MOTHER(I)
      IF (OMEGA(I1,I1).EQ.0.0D0.OR.OMEGA(I2,I2).EQ.0.0D0) GO TO 30
C
C     THE KINSHIP COEFFICIENTS FOR THE PARENTS OF I ARE ALREADY
C     KNOWN.  CALCULATE I'S KINSHIP COEFFICIENTS AGAINST ALL
C     PREVIOUSLY CONSIDERED PEOPLE AND AGAINST HIMSELF.
C
      DO 40 J=1,NPTOT
      IF (OMEGA(J,J).NE.0.0D0) THEN
      OMEGA(MAX(I,J),MIN(I,J))=0.5D0*(OMEGA(MAX(I1,J),MIN(I1,J))
     1+OMEGA(MAX(I2,J),MIN(I2,J)))
      END IF
 40   CONTINUE
      COUNT=COUNT+1
      OMEGA(I,I)=1.0D0+0.5D0*OMEGA(MAX(I1,I2),MIN(I1,I2))
C
      IGRP=IABS(GROUP(I))/100000
      IF (IGRP.NE.0) THEN
C
C     PICK UP THE MZ TWINS OF I.
C
      DO 50 J=1,NPTOT
      IF (J.NE.I.AND.IABS(GROUP(J))/100000.EQ.IGRP) THEN
      COUNT=COUNT+1
      OMEGA(J,J)=OMEGA(I,I)
      OMEGA(MAX(I,J),MIN(I,J))=OMEGA(I,I)
      DO 60 K=1,NPTOT
      IF (OMEGA(K,K).NE.0.0D0) THEN
      OMEGA(MAX(J,K),MIN(J,K))=OMEGA(MAX(I,K),MIN(I,K))
      END IF
 60   CONTINUE
      END IF
 50   CONTINUE
      END IF
 30   CONTINUE
C
C     CHECK TO SEE IF ANOTHER PASS THROUGH THE PEDIGREE IS NEEDED.
C
      IF (COUNT.LT.NPTOT) GO TO 70
      END
