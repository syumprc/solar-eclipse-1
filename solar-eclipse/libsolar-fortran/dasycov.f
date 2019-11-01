      SUBROUTINE DASYCOV(CNSTR,CVALUE,DF,HESS,PAR,TABLE,WORK1             
     1,WORK2,NSWEEP,PNAME,CNORM,SMALL,TOL,IOUNIT,MAXPAR,MAXTAB           
     2,MCNSTR,NCNSTR,NPAR,INASYCV)
C                                                                        
C     THIS SUBROUTINE COMPUTES THE ASYMPTOTIC COVARIANCE MATRIX          
C     FOR THE PARAMETER ESTIMATES.                                       
C                                                                        
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                 
      DOUBLE PRECISION CNSTR(MCNSTR,MAXPAR),CVALUE(MCNSTR),DF(MAXPAR)    
     1,HESS(MAXPAR,MAXPAR),PAR(MAXPAR),TABLE(MAXTAB,MAXTAB)              
     2,WORK1(MAXPAR),WORK2(MAXTAB)                                       
      INTEGER NSWEEP(MAXPAR)                                             
      CHARACTER*8 PNAME(MAXPAR)                                          
      LOGICAL INASYCV
C                                                                        
      NTIMES=1                                                           
C                                                                        
C     SET UP THE TABLEAU.  IF THE UPPER LEFT BLOCK IS NOT POSITIVE       
C     DEFINITE, THEN INCREASE THE CONTRIBUTION OF THE CONSTRAINTS.       
C     THIS WILL EVENTUALLY WORK DUE TO THEOREMS 6.1(CHAPS.1 AND 2)       
C     IN: M.R.HESTENES(1981). 'OPTIMIZATION THEORY: THE FINITE           
C     DIMENSIONAL CASE'  KRIEGER PUBLISHING. HUNTINGTON,NEW,YORK.        
C                                                                        
 130  CALL SETTAB(CNSTR,CVALUE,DF,HESS,PAR,TABLE,WORK2,CNORM             
     1,MAXPAR,MAXTAB,MCNSTR,NCNSTR,NPAR,NTAB)                            
C                                                                        
C     STORE THE DIAGONAL ELEMENTS OF TABLE IN WORK1 FOR CHECKING         
C     TOLERANCE.                                                         
C                                                                        
      DO 10 I=1,NPAR                                                     
      IF (NSWEEP(I).EQ.1.AND.TABLE(I,I).LE.0.0D0) GO TO 20               
 10   WORK1(I)=TABLE(I,I)                                                
C                                                                        
C     SWEEP ONLY ON THOSE PARAMETERS WHICH ARE NOT ON A BOUNDARY         
C     AND PASS THE TOLERANCE TEST.                                       
C                                                                        
      DO 30 I=1,NPAR                                                     
      IF (NSWEEP(I).NE.1) GO TO 30                                       
      IF (TABLE(I,I)/WORK1(I).LT.TOL) GO TO 20                           
      CALL SWEEP(TABLE,WORK2,I,MAXTAB,NTAB,.FALSE.)                      
 30   CONTINUE                                                           
C                                                                        
C     NEXT SWEEP ON THE CONSTRAINTS.  THE NEGATIVE OF THE                
C     ASYMPTOTIC COVARIANCE MATRIX SHOULD APPEAR IN THE                  
C     UPPER LEFT BLOCK OF THE TABLEAU.                                   
C                                                                        
      DO 40 I=1,NCNSTR                                                   
      K=I+NPAR                                                           
 40   IF (TABLE(K,K).LT.0.0D0)                                           
     :CALL SWEEP(TABLE,WORK2,K,MAXTAB,NTAB,.FALSE.)                      
C                                                                        
C     OUTPUT THE RESULTS.                                                
C                                                                        
      DO 50 J=1,NPAR                                                     
      IF (ABS(TABLE(J,J)).LT.SMALL) NSWEEP(J)=0                          
      IF (NSWEEP(J).EQ.1) THEN                                           
      IF (TABLE(J,J).GT.0.0D0) GO TO 20                                  
      WORK2(J)=SQRT(-TABLE(J,J))                                         
      ELSE                                                               
      WORK2(J)=0.0D0                                                     
      END IF                                                             
      DO 50 I=1,J                                                        
      IF (NSWEEP(I).EQ.1.AND.NSWEEP(J).EQ.1) THEN                        
      TABLE(I,J)=-TABLE(I,J)/(WORK2(I)*WORK2(J))                         
      ELSE                                                               
      TABLE(I,J)=0.0D0                                                   
      END IF                                                             
 50   CONTINUE                                                           
      IF (IOUNIT.GE.0) WRITE(IOUNIT,60)                                  
 60   FORMAT(//T18,'Asymptotic Standard Errors of the Parameters'/)
      IF (IOUNIT.GE.0) CALL WCNAMES (IOUNIT)
      IF (IOUNIT.GE.0) WRITE(IOUNIT,80) (WORK2(I),I=1,NPAR)              
 80   FORMAT(100(/6(1X,D11.4),:))                                        
      IF (IOUNIT.GE.0) WRITE(IOUNIT,90)                                  
 90   FORMAT(//T16,'Asymptotic Correlation Matrix of the Parameters'/)
      IF (IOUNIT.GE.0) CALL WCNAMES (IOUNIT)
      DO 100 J=1,NPAR                                                    
 100  IF (IOUNIT.GE.0) WRITE(IOUNIT,80) (TABLE(I,J),I=1,J)               
C                                                                        
C     COPY THE STANDARD ERRORS AND CORRELATIONS INTO HESS.               
C                                                                        
      DO 110 J=1,NPAR                                                    
      DO 120 I=1,J-1                                                     
      HESS(I,J)=TABLE(I,J)                                               
 120  HESS(J,I)=HESS(I,J)                                                
 110  HESS(J,J)=WORK2(J)                                                 
      INASYCV=.TRUE.
      RETURN                                                             
C                                                                        
C     INCREASE THE CONTRIBUTION OF THE CONSTRAINTS BY DECREASING         
C     CNORM.  RESET THE TABLEAU AND TRY AGAIN.                           
C                                                                        
 20   IF (NTIMES.LT.6) THEN                                              
      NTIMES=NTIMES+1                                                    
      CNORM=1.0D-1*CNORM                                                 
      GO TO 130                                                          
      END IF                                                             
      DO 140 J=1,NPAR                                                    
      DO 140 I=1,NPAR                                                    
 140  HESS(I,J)=0.0D0                                                    
      IF (IOUNIT.GE.0) WRITE(IOUNIT,150)                                 
 150  FORMAT(/T9,'Standard errors could not be computed')
      WRITE (4,150)
      WRITE (6,150)

      INASYCV=.FALSE.
      END                                                                
