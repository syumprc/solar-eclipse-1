      SUBROUTINE PEDTST(UNIFRM,DEGREE,PEDCUT,TPAR,NPED,UNIT3,NORMAL)
C
C     THIS SUBROUTINE FLAGS OUTLIER PEDEIGREES BASED ON THE PEDIGREE
C     QUADRATIC FORMS.  THE PVALUES OF THESE QUADRATIC FORMS ARE THEN
C     SUBJECTED TO VARIOUS EMPIRICAL DISTRIBUTION FUNCTION TESTS.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION UNIFRM(NPED)
      INTEGER DEGREE(NPED),PED,UNIT3,LID
      LOGICAL NORMAL
      CHARACTER*18 FIRST_ID
C
C     OUTPUT TWO HEADERS.
C
      WRITE(UNIT3,10)
 10   FORMAT(//,'                POSSIBLE OUTLIER PEDIGREES'
     1//' FOR THE NON-PROBANDS IN A PEDIGREE, IT IS POSSIBLE TO COMPUTE'
     2/' A QUADRATIC FORM WHICH UNDER A MULTIVARIATE NORMAL MODEL'
     3/' FOLLOWS AN APPROXIMATE CHI-SQUARE DISTRIBUTION.  THE DEGREES'
     3/' OF FREEDOM N EQUALS THE NUMBER OF NON-PROBANDS WITH ALL'
     4/' VARIABLES MEASURED TIMES THE NUMBER OF TRAITS.  UNDER A'
     5/' MULTIVARIATE T MODEL THIS QUADRATIC FORM DIVIDED BY N FOLLOWS'
     6/' AN APPROXIMATE F DISTRIBUTION WITH DEGREES OF FREEDOM N AND'
     7/' V, WHERE V IS THE CONVERGED T PARAMETER.  THOSE PEDIGREES'
     8/' WITH PVALUES FOR THE QUADRATIC FORM LESS THAN THE CUTOFF ARE'
     9/' LISTED BELOW.  THE VALUE OF N IS LISTED UNDER THE HEADING'
     1/' ''DEGREES OF FREEDOM''.  PEDIGREES INDENTIFIED BY INDIV. ID.')
C
      WRITE(UNIT3,20)
 20   FORMAT(/' PEDIGREE  DEG. FREE.  QUADRATIC FORM      PVALUE'/)
C
C     FIND THE OUTLIER PEDIGREES AND OUTPUT THEM.
C
      DO 30 PED=1,NPED
      QDFORM=UNIFRM(PED)
C
C     FOR THE NORMAL MODEL USE A CHI-SQUARE DISTRIBUTION AND FOR THE
C     T MODEL USE AN F DISTRIBUTION.  THE T PARAMETER IS PASSED VIA
C     TPAR.
C
      IF (NORMAL) THEN
      UNIFRM(PED)=CHISQ(QDFORM,DEGREE(PED))
      ELSE
      DFN=DBLE(DEGREE(PED))
      POINT=QDFORM/DFN
      UNIFRM(PED)=1.D0-FDIST(POINT,DFN,TPAR)
      END IF
      IF (UNIFRM(PED).LE.PEDCUT) THEN
      CALL FIRSTID (PED, FIRST_ID)
      LID = LNBLNK (FIRST_ID)
      IF (QDFORM.LT.100.) THEN
      WRITE(UNIT3,40) FIRST_ID(1:LID),DEGREE(PED),QDFORM,UNIFRM(PED)
 40   FORMAT(1X,A,6X,I6,9X,F8.4,8X,F7.5)
      ELSE
      WRITE(UNIT3,50) FIRST_ID(1:LID),DEGREE(PED),QDFORM,UNIFRM(PED)
 50   FORMAT(1X,A,6X,I6,7X,G10.4,8X,F7.5)
      END IF
      END IF
 30   CONTINUE
C
C     OUTPUT A HEADER FOR THE EDF TESTS.  THEN CALL SUBROUTINE EDFTST.
C
      WRITE(UNIT3,60)
 60   FORMAT(//' IN THE FOLLOWING STATISTICAL TESTS OF THE OVERALL FIT'
     1/' OF THE MODEL, THE PVALUES ASSOCIATED WITH THESE QUADRATIC'
     2/' FORMS ARE VIEWED AS RANDOM VARIABLES OR DEVIATES.')
      CALL EDFTST(UNIFRM,NPED,UNIT3)
      END
