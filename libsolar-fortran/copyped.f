C
C File:       copyped.f
C 
      SUBROUTINE COPYPED (UNIT1,UNIT2,IARRAY,RARRAY,LENI,LENR,NPED
     1,NVAR)
C
C Copy pedigree data into scratch file
C (Code from the now obsolete INPUT.F, after call to preped.)
C  Note: Problem definitions were not written to unit 1.)
C
      INTEGER UNIT1, UNIT2
      INTEGER LENI
      INTEGER*8 LENR
      INTEGER IARRAY(LENI)
      DOUBLE PRECISION RARRAY(LENR)
      INTEGER NPED, NVAR
      INTEGER NPBAND, NPEO, NPTOT
      SAVE INTEGER IPED, MAXPED
C
C     RETRIEVE THE PEDIGREE DATA FROM UNIT1 AND WRITE IT ON UNIT2.
C     This has changed.  Unit 2 is now always the scratch ped data file.
C     It is not closed and reopened.  It is opened in OPENOUT.
C
      CALL RESCRATCH(UNIT1)
C
C     NEXT READ AND REWRITE THE PEDIGREE DATA.  THE ORDER OF READING
C     ARRAYS IS FATHER, GROUP, MOTHER, AND VAR.
C
      MAXPED=NPED
      DO 50 IPED=1,MAXPED
C        PRINT *,"Copying pedigree",IPED
         CALL ISCRAT (NPBAND,1,UNIT1,.FALSE.)
         CALL ISCRAT (NPBAND,1,UNIT2,.TRUE.)

         CALL ISCRAT (NPEO,1,UNIT1,.FALSE.)
         CALL ISCRAT (NPEO,1,UNIT2,.TRUE.)

         CALL ISCRAT (NPTOT,1,UNIT1,.FALSE.)
         CALL ISCRAT (NPTOT,1,UNIT2,.TRUE.)

C        PRINT *, "NPTOT appears to be ",NPTOT

         CALL ISCRAT(IARRAY,NPTOT,UNIT1,.FALSE.)
         CALL ISCRAT(IARRAY,NPTOT,UNIT2,.TRUE.)

         CALL ISCRAT(IARRAY,NPTOT,UNIT1,.FALSE.)
         CALL ISCRAT(IARRAY,NPTOT,UNIT2,.TRUE.)

         CALL ISCRAT(IARRAY,NPTOT,UNIT1,.FALSE.)
         CALL ISCRAT(IARRAY,NPTOT,UNIT2,.TRUE.)

         MVAR=NPTOT*NVAR
         CALL RSCRAT(RARRAY,MVAR,UNIT1,.FALSE.)
         CALL RSCRAT(RARRAY,MVAR,UNIT2,.TRUE.)
 50   CONTINUE
C
C     READ VMEAN AND VVAR INTO THE FIRST TWO SECTIONS OF RARRAY.
C     THESE WILL BE PASSED TO THE REST OF THE PROGRAM AND REAPPEAR
C     WITH THEIR PROPER NAMES IN SUBROUTINE SEARCH.
C
      CALL RSCRAT (RARRAY, NVAR, UNIT1, .FALSE.)
      CALL RSCRAT (RARRAY(NVAR+1), NVAR, UNIT1, .FALSE.)
      END