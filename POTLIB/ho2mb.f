C                                                                               
      SUBROUTINE PREPOT                                                         
C                                                                               
C   System:    HO2                                                              
C   Reference: C. F. Melius and R. J. Blint                                     
C              Chem. Phys. Lett. 64, 183 (1979).                                
C                                                                               
C   PREPOT must be called once before any calls to POT.                         
C   The potential parameters are included in DATA statements.                   
C   and in the block data subprogram PTPACM.                                    
C   Coordinates, potential energy, and derivatives are passed                   
C   The potential energy in the three asymptotic valleys are                    
C   stored in the common block ASYCM:                                           
C                  COMMON /ASYCM/ EASYAB, EASYBC, EASYAC                        
C   The potential energy in the AB valley, EASYAB, is equal to the potential    
C   energy of the O "infinitely" far from the HO diatomic, with the             
C   HO diatomic at its equilibrium configuration.  Similarly, the terms         
C   EASYBC and EASYAC represent the O2 and the HO asymptotic valleys,           
C   respectively.                                                               
C   The other potential parameters are passed through the common blocks         
C   SATOCM and LEPSCM.                                                          
C   All the information passed through the common blocks PT1CM, ASYCM,          
C   SATOCM, and LEPSCM are in hartree atomic units.                             
C                                                                               
C        This potential is written such that:                                   
C                       R(1) = R(H-O)                                           
C                       R(2) = R(O-O)                                           
C                       R(3) = R(O-H)                                           
C   The classical potential energy is set equal to zero for the H               
C   infinitely far from the equilibrium O2 diatomic.                            
C                                                                               
C   The flags that indicate what calculations should be carried out in          
C   the potential routine are passed through the common block PT2CM:            
C   where:                                                                      
C        NASURF - which electronic states availalble                            
C                 (1,1) = 1 as only gs state available                          
C        NDER  = 0 => no derivatives should be calculated                       
C        NDER  = 1 => calculate first derivatives                               
C        NFLAG - these integer values can be used to flag options               
C                within the potential;                                          
C                                                                               
C                                                                               
C   Potential parameters' default settings                                      
C                  Variable            Default value                            
C                  NDER                1                                        
C                  NFLAG(18)           6                                        
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                                  
C         GENERAL INFORMATION                                                   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                                  
C                                                                               
C          THE POTENTIAL CONSISTS OF THREE MORSE TERMS FOR                      
C          THE THREE DIATOMICS, HO,OO, AND HO RESPECTIVELY                      
C          PLUS TWO INTERACTION TERMS.                                          
C                                                                               
C          EACH POTENTIAL TERM WILL BE CALCULATED INDIVIDUALLY                  
C          TO SIMPLIFY READING THE CODE.                                        
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                           
C          DEFINE VARIABLES, ARRAYS, AND COMMON BLOCK                           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                           
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
      COMMON /PT3CM/ EZERO(ISURF+1)                                             
C                                                                               
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),                            
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
C                                                                               
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
      COMMON /ASYCM/ EASYAB,EASYBC,EASYAC                                       
C                                                                               
      REAL*8 NUM,DENOM,GAMMA,JUNK,EXPO,DEXPO,MEW                                
C                                                                               
C                                                                               
C                                                                               
         COMMON /POTCM/  C(5,3,3),DE(3),GAMMA(3),REQ(3),                        
     +                   ZEROAD,ALPHA,BETA                                      
         COMMON /PRECM/  DCSTH(3),DV(3),JUNK(4),A(5),                           
     +                   DA(5,3),EXPO(3),DEXPO(3),DMEW(3),                      
     +                   B(5,3),DB(5,3,3),DPROD1(3),                            
     +                   DPROD2(3),DSEC(3),DTHIRD(3),                           
     +                   DPROD3(3),ZKCAL                                        
C                                                                               
      IF(NATOMS.GT.25) THEN                                                     
         WRITE(NFLAG(18),1111)                                                  
 1111    FORMAT(2X,'STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS')            
         STOP                                                                   
      END IF                                                                    
C                                                                               
C                                                                               
C   Initialize the potential values in the three asymptotic valleys.            
C                                                                               
      EASYAB = DE(1)                                                            
      EASYBC = DE(2)                                                            
      EASYAC = DE(3)                                                            
C                                                                               
C   Echo potential parameters to unit IPRT.                                     
C                                                                               
      WRITE (NFLAG(18), 1000)                                                   
      WRITE (NFLAG(18),1050) (DE(IT),IT=1,3)                                    
      WRITE (NFLAG(18),1150) (GAMMA(IT),IT=1,3)                                 
      WRITE (NFLAG(18),1200) (REQ(IT),IT=1,3)                                   
C                                                                               
 1000 FORMAT (/, 2X, T5, 'PREPOT has been called for H + OO',                   
     *        /, 2X, T5, 'Potential energy function by Melius ',                
     *                   'and Blint',                                           
     *        //, 2X, T5, 'Potential energy surface parameters:',               
     *        /, 2X, T10, 'Bond', T47, 'H-O', T58, 'O-O', T69, 'O-H')           
 1050 FORMAT(2X,T10,'Dissociation energies (hartrees):',                        
     *       T44, F10.5, T55, F10.5, T66, F10.5)                                
 1150 FORMAT(2X,T10,'Morse betas (reciprocal bohr):',                           
     *       T44, F10.5, T55, F10.5, T66, F10.5)                                
 1200 FORMAT(2X,T10,'Equilibrium bond lengths (bohr):',                         
     *       T44, F10.5, T55, F10.5, T66, F10.5)                                
C                                                                               
C   Convert ZEROAD to eV and kcal/mol and echo to unit IPRT.                    
C                                                                               
      ZEV = ZEROAD*27.21106D0                                                   
      ZKCAL = ZEROAD*627.5095D0                                                 
      WRITE (NFLAG(18),1300) ZEROAD,ZEV,ZKCAL                                   
C                                                                               
 1300 FORMAT(/,2X,T10,'Constant added to the energy:',                          
     *       /,T38, 1PE20.12,' hartree atomic units',                           
     *       /,T38, 1PE20.12,' eV',/,T38,1PE20.12,' kcal/mol')                  
C                                                                               
CCCCCCCCCCCC                                                                    
C                                                                               
C                                                                               
      EZERO(1)=DE(2)                                                            
C                                                                               
       DO I=1,5                                                                 
          REF(I) = ' '                                                          
       END DO                                                                   
C                                                                               
       REF(1)='C. F. Melius and R. J. Blint'                                    
       REF(2)='Chem. Phys. Lett. 64, 183(1979)'                                 
C                                                                               
      INDEXES(1) = 1                                                            
      INDEXES(2) = 8                                                            
      INDEXES(3) = 8                                                            
C                                                                               
C                                                                               
C                                                                               
      IRCTNT=2                                                                  
C                                                                               
      CALL POTINFO                                                              
C                                                                               
      CALL ANCVRT                                                               
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
CCCCCCCCCCCC                                                                    
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                 
C        MOVE ONTO POT SUBROUTINE                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                 
C                                                                               
      SUBROUTINE POT                                                            
C                                                                               
C   The potential energy in the AB valley, EASYAB, is equal to the potential    
C   energy of the O "infinitely" far from the HO diatomic, with the             
C   HO diatomic at its equilibrium configuration.  Similarly, the terms         
C   EASYBC and EASYAC represent the O2 and the HO asymptotic valleys,           
C   respectively.                                                               
C                                                                               
C                                                                               
C        This potential is written such that:                                   
C                       R(1) = R(H-O)                                           
C                       R(2) = R(O-O)                                           
C                       R(3) = R(O-H)                                           
C   The classical potential energy is set equal to zero for the H               
C   infinitely far from the equilibrium O2 diatomic.                            
C                                                                               
C      ENTRY POT                                                                
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
      COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                            
      COMMON /PT3CM/ EZERO(ISURF+1)                                             
      COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                         
      COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                         
C                                                                               
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),                            
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),                                    
     +               PENGYIJ(JSURF),                                            
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),                   
     +               DIJCART(NATOM,3,JSURF)                                     
C                                                                               
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
      COMMON /ASYCM/ EASYAB,EASYBC,EASYAC                                       
C                                                                               
      REAL*8 NUM,DENOM,GAMMA,JUNK,EXPO,DEXPO,MEW                                
C                                                                               
C                                                                               
C                                                                               
         COMMON /POTCM/  C(5,3,3),DE(3),GAMMA(3),REQ(3),                        
     +                   ZEROAD,ALPHA,BETA                                      
         COMMON /PRECM/  DCSTH(3),DV(3),JUNK(4),A(5),                           
     +                   DA(5,3),EXPO(3),DEXPO(3),DMEW(3),                      
     +                   B(5,3),DB(5,3,3),DPROD1(3),                            
     +                   DPROD2(3),DSEC(3),DTHIRD(3),                           
     +                   DPROD3(3),ZKCAL                                        
C                                                                               
      CALL CARTOU                                                               
      CALL CARTTOR                                                              
C                                                                               
C   Check the values of NASURF and NDER for validity.                           
C                                                                               
      IF (NASURF(1,1) .EQ. 0) THEN                                              
         WRITE(NFLAG(18), 900) NASURF(1,1)                                      
         STOP                                                                   
      ENDIF                                                                     
         IF (NDER .GT. 1) THEN                                                  
             WRITE (NFLAG(18), 910) NDER                                        
             STOP 'SURF 2'                                                      
         ENDIF                                                                  
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                                
C         SET RAB,RBC,AND RAC                                                   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                                   
C                                                                               
      RBC = R(2)                                                                
C                                                                               
C           RAB IS DEFINED AS THE SMALLER OF R(1) AND R(3)                      
C           RAC IS DEFINED AS THE LARGER OF R(1) AND R(3)                       
C                                                                               
      IF (R(1).GT.R(3)) GO TO 40                                                
C                                                                               
      RAB = R(1)                                                                
      RAC = R(3)                                                                
      GO TO 50                                                                  
C                                                                               
   40 RAB = R(3)                                                                
      RAC = R(1)                                                                
   50 CONTINUE                                                                  
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                        
C          CALCULATE USEFUL EXPONENTIALS                                        
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                        
C                                                                               
      EXPO(1) = EXP(-ALPHA*RAB)                                                 
      EXPO(2) = EXP(-BETA*RBC)                                                  
      EXPO(3) = EXP(-ALPHA*RAC)                                                 
C                                                                               
C          AND THEIR DERIVATIVES                                                
C                                                                               
         IF (NDER .EQ. 1) THEN                                                  
             DEXPO(1) = -ALPHA*EXPO(1)                                          
             DEXPO(2) = -BETA*EXPO(2)                                           
             DEXPO(3) = -ALPHA*EXPO(3)                                          
         ENDIF                                                                  
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                               
C          CALCULATE VAB AND DV(1)                                              
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                               
C                                                                               
C          THE POTENTIAL CALCULATED HERE AND THE POTENTIALS VAC AND VBC         
C          CORRESPOND TO THE VOH POTENTIAL DESCRIBED IN EQUATION 3.             
C          DV(1) IS THE DERIVATIVE OF THIS TERM WITH RESPECT TO RAB.            
C                                                                               
      JUNK(1) = EXP(GAMMA(1)*(REQ(1)-RAB))                                      
      VAB = DE(1)*(JUNK(1)*JUNK(1)-2.0D0*JUNK(1))                               
      IF (NDER .EQ. 1)                                                          
     *    DV(1) = 2.0D0*DE(1)*(1.0D0-JUNK(1))*GAMMA(1)*JUNK(1)                  
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                               
C          CALCULATE VAC AND DV(3)                                              
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                               
C                                                                               
      JUNK(1) = EXP(GAMMA(3)*(REQ(3)-RAC))                                      
      VAC = DE(3)*(JUNK(1)*JUNK(1)-2.0D0*JUNK(1))                               
      IF (NDER .EQ. 1)                                                          
     *    DV(3) = 2.0D0*DE(3)*(1.0D0-JUNK(1))*GAMMA(3)*JUNK(1)                  
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                               
C          CALCULATE VBC AND DV(2)                                              
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                               
C                                                                               
      JUNK(1) = EXP(GAMMA(2)*(REQ(2)-RBC))                                      
      VBC = DE(2)*(JUNK(1)*JUNK(1)-2.0D0*JUNK(1))                               
      IF (NDER .EQ. 1)                                                          
     *    DV(2) = -2.0D0*DE(2)*(JUNK(1)-1.0D0)*GAMMA(2)*JUNK(1)                 
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                    
C          CALCULATE CSTH AND DERIVATIVES                                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                    
C                                                                               
C          TH IS THE BEND ANGLE OF THE THREE ATOM SYSTEM.                       
C          IT IS NEVER USED EXPLICITELY WITHIN THE PROGRAM                      
C          HOWEVER COSINE(TH) IS.  COSINE(TH) IS CALCULATED                     
C          FROM THE LAW OF COSINES)                                             
C                                                                               
      NUM = RBC*RBC+RAB*RAB-RAC*RAC                                             
      DENOM = 2.0D0*RAB*RBC                                                     
C                                                                               
      CSTH = NUM/DENOM                                                          
C                                                                               
      IF (NDER .EQ. 1) THEN                                                     
          DCSTH(1) = (2.0D0*RAB/DENOM)-(2.0D0*RBC)*CSTH/DENOM                   
          DCSTH(2) = (2.0D0*RBC/DENOM)-(2.0D0*RAB)*CSTH/DENOM                   
          DCSTH(3) = -(2.0D0*RAC/DENOM)                                         
      ENDIF                                                                     
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                
C          CALCULATE THE B IJ 'S--SEE EQUATION 5                                
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                
C                                                                               
      DO 80 I = 1, 5                                                            
         DO 70 J = 1, 3                                                         
            B(I,J) = C(I,J,1)*(1.0D0+C(I,J,2)*CSTH*(1.0D0+C(I,J,3)*             
     *               CSTH))                                                     
C                                                                               
C               CALCULATE THE DERIVATIVES TOO.                                  
C                                                                               
            IF (NDER .EQ. 1) THEN                                               
                DO 60 IT = 1, 3                                                 
                      DB(I,J,IT) = C(I,J,1)*C(I,J,2)*(CSTH*C(I,J,3)*            
     *                             DCSTH(IT)+DCSTH(IT)*                         
     *                             (1.0D0+C(I,J,3)*CSTH))                       
   60           CONTINUE                                                        
            ENDIF                                                               
   70    CONTINUE                                                               
   80 CONTINUE                                                                  
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                
C          NEXT CALCULATE MEW--SEE EQUATION FOUR                                
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                
C                                                                               
      MEW = EXPO(2)-0.02538745816335D0                                          
      IF (NDER .EQ. 1) THEN                                                     
          DMEW(2) = -BETA*EXPO(2)                                               
          DMEW(1) = 0.0D0                                                       
          DMEW(3) = 0.0D0                                                       
      ENDIF                                                                     
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                 
C          NOW CALCULATE A I 'S AND DERIVATIVES--SEE EQUATION 4                 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                 
C                                                                               
      DO 100 IT = 1, 5                                                          
         JUNK(1) = 1.0D0+B(IT,3)*MEW                                            
         JUNK(2) = MEW*B(IT,2)                                                  
         JUNK(3) = 1.0D0+JUNK(2)*JUNK(1)                                        
C                                                                               
         A(IT) = B(IT,1)*JUNK(3)                                                
         IF (NDER .EQ. 1) THEN                                                  
             DO 90 IT1 = 1, 3                                                   
                   DA(IT,IT1) = DB(IT,1,IT1)*JUNK(3)+B(IT,1)*                   
     *                          ((DB(IT,2,IT1)*MEW+B(IT,2)*DMEW(IT1))*          
     *                          JUNK(1)+JUNK(2)*(DB(IT,3,IT1)*MEW+              
     *                          B(IT,3)*DMEW(IT1)))                             
   90        CONTINUE                                                           
         ENDIF                                                                  
  100 CONTINUE                                                                  
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                          
C          CALCULATE THE SECOND TERM IN EQUATION THREE                          
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                          
C                                                                               
      PROD1 = EXPO(1)*EXPO(3)*A(1)                                              
C                                                                               
      PROD2 = 1.0D0+(A(3)*RAB*RAC)/(RAB+RAC)                                    
C                                                                               
      IF (NDER .EQ. 1) THEN                                                     
          DPROD1(1) = DEXPO(1)*EXPO(3)*A(1)+EXPO(1)*EXPO(3)*DA(1,1)             
          DPROD1(2) = EXPO(1)*EXPO(3)*DA(1,2)                                   
          DPROD1(3) = DEXPO(3)*EXPO(1)*A(1)+EXPO(1)*EXPO(3)*DA(1,3)             
C                                                                               
          DPROD2(2) = (DA(3,2)*RAB*RAC)/(RAB+RAC)                               
          DPROD2(1) = ((DA(3,1)*RAB*RAC+A(3)*RAC)*(RAB+RAC)-                    
     *                 (A(3)*RAB*RAC))/((RAB+RAC)*(RAB+RAC))                    
          DPROD2(3) = ((DA(3,3)*RAB*RAC+A(3)*RAB)*(RAB+RAC)-                    
     *                 (A(3)*RAB*RAC))/((RAB+RAC)*(RAB+RAC))                    
      ENDIF                                                                     
C                                                                               
C          CALCULATE PROD3                                                      
C                                                                               
      PROD3 = 1.0D0+A(2)*(RAB+RAC)*PROD2                                        
C                                                                               
      IF (NDER .EQ. 1) THEN                                                     
          DO 120 IT = 1, 3                                                      
                 DPROD3(IT) = DA(2,IT)*(RAB+RAC)*PROD2+A(2)*                    
     *                        (RAB+RAC)*DPROD2(IT)                              
                 IF (IT.EQ.2) GO TO 110                                         
                 DPROD3(IT) = DPROD3(IT)+A(2)*PROD2                             
  110            CONTINUE                                                       
  120     CONTINUE                                                              
      ENDIF                                                                     
C                                                                               
      SECOND = PROD1*PROD3                                                      
C                                                                               
      IF (NDER .EQ. 1) THEN                                                     
          DO 130 IT = 1, 3                                                      
                 DSEC(IT) = PROD1*DPROD3(IT)+DPROD1(IT)*PROD3                   
  130     CONTINUE                                                              
      ENDIF                                                                     
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                             
C          CALCULATE THE THIRD TERM                                             
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                             
C                                                                               
      JUNK(1) = EXPO(1)+EXPO(3)                                                 
      JUNK(2) = (1.0D0+A(5)*(RAB+RAC))                                          
C                                                                               
      THIRD = JUNK(1)*EXPO(2)*A(4)*JUNK(2)                                      
C                                                                               
      IF (NDER .EQ. 1) THEN                                                     
          DO 140 IT1 = 1, 2                                                     
                 IT = 1                                                         
                 IF (IT1.EQ.2) IT = 3                                           
                 DTHIRD(IT) = DEXPO(IT)*EXPO(2)*A(4)*JUNK(2)+                   
     *                        JUNK(1)*EXPO(2)*DA(4,IT)*JUNK(2)+                 
     *                        JUNK(1)*EXPO(2)*A(4)*(DA(5,IT)*                   
     *                        (RAB+RAC)+A(5))                                   
  140     CONTINUE                                                              
          DTHIRD(2) = JUNK(1)*(DEXPO(2)*A(4)*JUNK(2)+EXPO(2)*(DA(4,2)*          
     *                JUNK(2)+A(4)*(DA(5,2)*(RAB+RAC))))                        
      ENDIF                                                                     
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                                
C          FINALLY SUM IT ALL UP                                                
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                                
C                                                                               
      ENGYGS = VAB+VAC+VBC+SECOND+THIRD+ EZERO(1)                               
      ENGYGS = ENGYGS+ZEROAD                                                    
C                                                                               
      IF (NDER .EQ. 1) THEN                                                     
          DO 150 IT = 1, 3                                                      
                 DEGSDR(IT) = DV(IT)+DSEC(IT)+DTHIRD(IT)                        
  150     CONTINUE                                                              
C                                                                               
          IF (R(1).LT.R(3)) GO TO 160                                           
          JUNK(1) = DEGSDR(3)                                                   
          DEGSDR(3) = DEGSDR(1)                                                 
          DEGSDR(1) = JUNK(1)                                                   
  160     CONTINUE                                                              
      ENDIF                                                                     
C                                                                               
CCCCCCCCCCCCCCCC                                                                
C                                                                               
      CALL EUNITZERO                                                            
      IF(NDER.NE.0) THEN                                                        
         CALL RTOCART                                                           
         CALL DEDCOU                                                            
      ENDIF                                                                     
C                                                                               
      RETURN                                                                    
C                                                                               
CCCCCCCCCCCCCCCC                                                                
C                                                                               
1000  FORMAT (/, 2X, T5, 'PREPOT has been called for H + OO',                   
     *        /, 2X, T5, 'Potential energy function by Melius ',                
     *                   'and Blint',                                           
     *        //, 2X, T5, 'Potential energy surface parameters:',               
     *        /, 2X, T10, 'Bond', T47, 'H-O', T58, 'O-O', T69, 'O-H')           
 1050 FORMAT(2X,T10,'Dissociation energies (hartrees):',                        
     *       T44, F10.5, T55, F10.5, T66, F10.5)                                
 1150 FORMAT(2X,T10,'Morse betas (reciprocal bohr):',                           
     *       T44, F10.5, T55, F10.5, T66, F10.5)                                
 1200 FORMAT(2X,T10,'Equilibrium bond lengths (bohr):',                         
     *       T44, F10.5, T55, F10.5, T66, F10.5)                                
 1300 FORMAT(/,2X,T10,'Constant added to the energy:',                          
     *       /,T38, 1PE20.12,' hartree atomic units',                           
     *       /,T38, 1PE20.12,' eV',/,T38,1PE20.12,' kcal/mol')                  
 900  FORMAT(/,2X,T5,13HNASURF(1,1) =,I5,                                       
     *       /,2X,T5,24HThis value is unallowed.                                
     *       /,2X,T5,31HOnly gs surface=>NASURF(1,1)=1 )                        
910   FORMAT(/, 2X,'POT has been called with NDER = ',I5,                       
     *       /, 2X, 'This value of NDER is not allowed in this ',               
     *              'version of the potential.')                                
C                                                                               
      END                                                                       
C                                                                               
C*****                                                                          
C                                                                               
         BLOCK DATA PTPACM                                                      
         IMPLICIT REAL*8 (A-H,O-Z)                                              
C                                                                               
      CHARACTER*75 REF(5)                                                       
      REAL*8 JUNK
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
      COMMON /PT3CM/ EZERO(ISURF+1)                                             
C                                                                               
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),                            
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
C                                                                               
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
      COMMON /ASYCM/ EASYAB,EASYBC,EASYAC                                       
C                                                                               
         COMMON /POTCM/  C(5,3,3),DE(3),GAMMA(3),REQ(3),                        
     +                   ZEROAD,ALPHA,BETA                                      
         COMMON /PRECM/  DCSTH(3),DV(3),JUNK(4),A(5),                           
     +                   DA(5,3),EXPO(3),DEXPO(3),DMEW(3),                      
     +                   B(5,3),DB(5,3,3),DPROD1(3),                            
     +                   DPROD2(3),DSEC(3),DTHIRD(3),                           
     +                   DPROD3(3),ZKCAL                                        
C                                                                               
C   Initialize the flags and the I/O unit numbers for the potential             
C                                                                               
      DATA NASURF /1,35*0/                                                      
      DATA NDER /0/                                                             
         DATA NFLAG /1,1,15*0,6,0,0/                                            
C                                                                               
      DATA ANUZERO /0.0D0/                                                      
      DATA ICARTR,MSURF,MDER/3,0,1/                                             
      DATA NULBL /25*0/                                                         
      DATA NATOMS /3/                                                           
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                             
C          Initialize the constants                                             
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                             
C                                                                               
      DATA C / 77.45D0,    -0.4071D0,    -0.508D0,                              
     *         0.1489D0,    1.013D0,    -10.495D0,                              
     *         9.1484D0,   10.273D0,      0.0D0,                                
     +         0.0D0,      -3.050D0,    -33.78D0,                               
     +       -22.951D0,     0.0D0,        0.0D0,                                
     +        -1.3699D0,    0.4101D0,     0.3411D0,                             
     *        16.906D0,     0.0005637D0, -0.6359D0,                             
     *        -0.06253D0,  -0.1225D0,     0.0D0,                                
     +         0.0D0,      -0.00906D0,    0.4435D0,                             
     +         0.7439D0,    0.0D0,        0.0D0,                                
     *        -0.3498D0,    0.6617D0,     0.8677D0,                             
     +         1.3858D0,  233.36D0,      -4.756D0,                              
     +        -1.626D0,     0.6337D0,     0.0D0,                                
     +         0.0D0,    -476.32D0,      -1.2104D0,                             
     +        -1.2762D0,    0.0D0,        0.0D0 /                               
C                                                                               
C   The array DE contains the dissociation energies in hartree atomic units.    
C                                                                               
      DATA DE / 0.1559D0, 0.1779D0, 0.1559D0/                                   
C                                                                               
C   The array GAMMA contains the Morse Betas in reciprocal bohrs.               
C                                                                               
      DATA GAMMA / 1.2670D0, 1.4694D0, 1.2670D0/                                
C                                                                               
C   The array REQ contains the equilibrium bond lengths in bohr.                
C                                                                               
      DATA REQ / 1.8460D0, 2.3158D0, 1.8460D0/                                  
C                                                                               
C   The variable ZEROAD contains the constant                                   
C   added to the energy in hartree atomic units.                                
C                                                                               
      DATA ZEROAD /-0.022000295721D0/                                           
C                                                                               
C          SET ALPHA, UNITS: RECIPROCAL BOHR                                    
C                                                                               
      DATA ALPHA /0.9172D0/                                                     
C                                                                               
C          SET BETA, UNITS: RECIPROCAL BOHR                                     
C                                                                               
      DATA BETA /1.4694D0/                                                      
C                                                                               
         END                                                                    
C                                                                               
C*****                                                                          
