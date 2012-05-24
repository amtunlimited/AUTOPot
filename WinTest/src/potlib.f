C                                                                               
      SUBROUTINE PREPOT                                                         
C                                                                               
C   System:          H3                                                         
C   Common name:     TK                                                         
C   Reference:       D. G. Truhlar and A. Kuppermann                            
C                    J. Chem. Phys. 56, 2232 (1984).                            
C Number of bodies: 3
C Interface: potlib2001
C Number of electronic states: 1
C Number of derivatives: 1
C                                                                               
C   PREPOT must be called once before any calls to POT.                         
C   The potential parameters are included in DATA statements.                   
C   and in the block data subprogram PTPACM.                                    
C   Coordinates, potential energy, and derivatives are passed                   
C                                                                               
C   COORDINATES:                                                                
C                                                                               
C            R(1) = RAB                                                         
C            R(2) = RBC                                                         
C            R(3) = RAC                                                         
C                                                                               
C   The potential energy in the three asymptotic valleys are                    
C   stored in the common block ASYCM:                                           
C                  COMMON /ASYCM/ EASYAB, EASYBC, EASYAC                        
C   The potential energy in the AB valley, EASYAB, is equal to the potential    
C   energy of the H "infinitely" far from the H2 diatomic, with the             
C   H2 diatomic at its equilibrium configuration.  For this potential energy    
C   surface the potential energy in the three asymptotic valleys are equivalent.
C   All the information passed through the common blocks PT1CM and ASYCM        
C   is in Hartree atomic units.                                                 
C                                                                               
C   The potential energy is defined to be equal to zero when one H is           
C   "infinitely" far from the H2 diatomic, and the H2 diatomic bond             
C   length is equal to the H2 equilibrium bond length.                          
C                                                                               
C   The flags that indicate what calculations should be carried out in          
C   the potential routine are passed through the common block PT2CM:            
C   where:                                                                      
C        NASURF - which electronic states availalble                            
C                 (1,1) = 1 as only gs state available                          
C        NDER  - the order of the derivatives of the energy with respect to     
C                the coordinates to be calculates.                              
C        NDER  = 1 => calculate first derivatives                               
C                This is the only option supported in this version              
C                of the potential.                                              
C        NFLAG - these integer values can be used to flag options               
C                within the potential;                                          
C                                                                               
C                                                                               
C   Potential parameters' default settings                                      
C                  Variable            Default value                            
C                  NDER                1                                        
C                  NFLAG(18)           6                                        
C                                                                               
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
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
      DOUBLE PRECISION NAB,NBC,NAB1,NBC1,JUNK                                   
C                                                                               
C   Conversion constants                                                        
C   BOHR (Angstrom/bohr) converts Angstroms to Bohr                             
C   HARTRE (kcal/Hartree) converts from Hartree to kcal/mol                     
C                                                                               
         PARAMETER (BOHR = 0.52917706D0)                                        
         PARAMETER (HARTRE = 627.5095D0)                                        
C                                                                               
         COMMON /PRECM/  DE,BETA,A,GAMMA,RBB,RHH,                               
     +                   NAB,NBC,NAB1,NBC1,                                     
     +                   DBCOOR(2),JUNK(4),REQ(3)                               
         COMMON /PRE2CM/ ALFD,A2,SN,D,ALFN,XN,RNMSN,                            
     +                   RN,RNMXN,ALFDIF,R2T,L                                  
C                                                                               
      IF(NATOMS.GT.25) THEN                                                     
         WRITE(NFLAG(18),1111)                                                  
 1111    FORMAT(2X,'STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS')            
         STOP                                                                   
      END IF                                                                    
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                      
C      SET COLLINEAR POTENTIAL PARAMETERS                                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                      
C                                                                               
      CALL PREPOT2                                                              
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                       
C   Echo the bending potential parameters                                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                       
C                                                                               
C   DE    = dissociation energy in kcal/mol                                     
C   BETA  = Morse beta in reciprocal Angstroms                                  
C   A     = Pauling parameter in Angstroms                                      
C   GAMMA = gamma for the bend correction                                       
C           gamma is dimensionless and ranges from 0.5 to 0.65                  
C   RBB   = parameter used in the scaling calculation for the                   
C           bending correction                                                  
C                                                                               
      WRITE (NFLAG(18),600) DE, DE, DE, (REQ(I),I=1,3),BETA, BETA, BETA         
      WRITE (NFLAG(18), 610) A, GAMMA, RBB                                      
C                                                                               
600   FORMAT (/,2X,T5,'PREPOT has been called for the H3 ',                     
     *                'potential energy surface TK',                            
     *       //,2X,T5,'Potential energy surface parameters:',                   
     *        /,2X,T5,'Bond', T47, 'H-H', T58, 'H-H', T69, 'H-H',               
     *        /,2X,T5,'Dissociation energies (kcal/mol):',                      
     *        T44, F10.5, T55, F10.5, T66, F10.5,                               
     *        /,2X,T5,'Equilibrium bond lengths (Angstroms):',                  
     *        T44, F10.5, T55, F10.5, T66, F10.5,                               
     *        /,2X,T5,'Morse beta parameters (Angstroms**-1):',                 
     *        T44, F10.5, T55, F10.5, T66, F10.5)                               
610   FORMAT(/,2X,T5,'Pauling parameter (Angstroms):',T44,F10.5,                
     *       /,2X,T5,'Gamma for the bending correction:',T44,F10.5,             
     *       /,2X,T5,'RBB (Angstroms):',T44,F10.5)                              
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                                 
C       CONVERT TO ATOMIC UNITS                                                 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                                 
C                                                                               
      DO 50 IT = 1,3                                                            
      REQ(IT) = REQ(IT)/BOHR                                                    
   50 CONTINUE                                                                  
      BETA = BETA * BOHR                                                        
      DE = DE/HARTRE                                                            
      A = A/BOHR                                                                
      RHH = 1.4007977D0                                                         
      RBB = RBB /BOHR                                                           
C                                                                               
CCCCCCCCCCC                                                                     
C                                                                               
C***********************                                                        
C                      *                                                        
C   The vale of D is   *                                                        
C   set in PREPOT2     *                                                        
C                      *                                                        
C***********************                                                        
       EZERO(1)=D                                                               
C                                                                               
       DO I=1,5                                                                 
          REF(I) = ' '                                                          
       END DO                                                                   
C   Reference:                                                                  
C                                                                               
C                                                                               
       REF(1)='D. G. Truhlar and A. Kuppermann'                                 
       REF(2)='J. Chem. Phys. 56, 2232(1984)'                                   
C                                                                               
      INDEXES(1) = 1                                                            
      INDEXES(2) = 1                                                            
      INDEXES(3) = 1                                                            
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
CCCCCCCCCCCC                                                                    
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                 
C        MOVE ONTO THE POT SUBROUTINE                                           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                 
C                                                                               
      SUBROUTINE POT                                                            
C                                                                               
C                                                                               
C   COORDINATES:                                                                
C                                                                               
C            R(1) = RAB                                                         
C            R(2) = RBC                                                         
C            R(3) = RAC                                                         
C                                                                               
C   The potential energy in the AB valley, EASYAB, is equal to the potential    
C   energy of the H "infinitely" far from the H2 diatomic, with the             
C   H2 diatomic at its equilibrium configuration.  For this potential energy    
C   surface the potential energy in the three asymptotic valleys are equivalent.
C   All the information passed through the common blocks PT1CM and ASYCM        
C   is in Hartree atomic units.                                                 
C                                                                               
C   The potential energy is defined to be equal to zero when one H is           
C   "infinitely" far from the H2 diatomic, and the H2 diatomic bond             
C   length is equal to the H2 equilibrium bond length.                          
C                                                                               
CCCCCCCCCCCCCCCC                                                                
C      ENTRY POT                                                                
CCCCCCCCCCCCCCCC                                                                
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                    
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
         DOUBLE PRECISION NAB,NBC,NAB1,NBC1,JUNK                                
C                                                                               
C   Conversion constants                                                        
C   BOHR (Angstrom/bohr) converts Angstroms to Bohr                             
C   HARTRE (kcal/Hartree) converts from Hartree to kcal/mol                     
C                                                                               
         PARAMETER (BOHR = 0.52917706D0)                                        
         PARAMETER (HARTRE = 627.5095D0)                                        
C                                                                               
         COMMON /PRECM/  DE,BETA,A,GAMMA,RBB,RHH,                               
     +                   NAB,NBC,NAB1,NBC1,                                     
     +                   DBCOOR(2),JUNK(4),REQ(3)                               
         COMMON /PRE2CM/ ALFD,A2,SN,D,ALFN,XN,RNMSN,                            
     +                   RN,RNMXN,ALFDIF,R2T,L                                  
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
             STOP 'POT 2'                                                       
         ENDIF                                                                  
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                    
C        FIRST CALCULATE SOME USEFUL NUMBERS                                    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                    
C                                                                               
C        NAB AND NBC ARE THE BOND ORDERS OF THE DIATOMICS                       
C        AB AND BC RESPECTIVELY. (EQUATION 2.2)                                 
C                                                                               
      NAB = EXP((REQ(1) - R(1))/A)                                              
C                                                                               
      IF (NAB.LT.1.D-15) NAB = 1.D-15                                           
C                                                                               
      NBC = EXP((REQ(2) - R(2))/A)                                              
C                                                                               
      IF (NBC.LT.1.D-15) NBC = 1.D-15                                           
C                                                                               
C       CALCULATE BOND ORDERS ALONG THE REACTION COORDINATE                     
C       SEE EQUATION 2.7A, 2.7B                                                 
C                                                                               
      C = 1.D0 - NAB/NBC                                                        
      IF (ABS(C).GT.1.D-14) GOTO 90                                             
      NAB1=.5D0                                                                 
      NBC1=.5D0                                                                 
      GO TO 95                                                                  
   90 C = C/NAB                                                                 
      NAB1 = (2.D0 + C - SQRT(4.D0+C*C))/(2.D0*C)                               
      NBC1 = 1.D0 - NAB1                                                        
C                                                                               
C          ROB IS USED IN THE BENDING POTENTIAL CALCULATIONS                    
C                                                                               
C          FIRST TRAP OUT ANY ZERO ARGUEMENTS                                   
C                                                                               
   95 IF (NAB1*NBC1.GT.0.0D0) GO TO 103                                         
           ROB = 1.D0                                                           
      GO TO 107                                                                 
C                                                                               
  103 STUFF = A * LOG(NAB1*NBC1)                                                
      ROB = (RBB - STUFF)/(RHH - STUFF)                                         
C                                                                               
  107 CONTINUE                                                                  
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                           
C        CALCULATE BENDING CORRECTION                                           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                           
C                                                                               
C        CALCULATE V(R(3))                                                      
C                                                                               
      JUNK(3) = EXP(-BETA*(R(3)-REQ(3))/ROB)                                    
      VAC = GAMMA * DE * JUNK(3)* (1.D0 + 0.5D0*JUNK(3))                        
C                                                                               
C        CALCULATE V(R(1)+R(2))                                                 
C                                                                               
      JUNK(4) = EXP(BETA*(REQ(3)-R(2)-R(1))/ROB)                                
      VABBC= GAMMA * DE * JUNK(4) * (1.D0 + 0.5D0*JUNK(4))                      
C                                                                               
C        HERE IS THE BENDING CORRECTION                                         
C                                                                               
      BCOOR = VAC - VABBC                                                       
C                                                                               
C          WHILE IT'S CONVENIENT, CALCULATE THE DERIVATIVE                      
C          OF BCOOR WITH RESPECT TO R(1)  R(2) AND R(3).                        
C                                                                               
C        CALCULATE DAC (THE DERIVATIVE WITH RESPECT TO R(3))                    
C                                                                               
      DAC = -GAMMA*DE*BETA*JUNK(3)*(1.D0+JUNK(3))/ROB                           
      DBCOOR(1) = GAMMA * DE * BETA * JUNK(4) * (1.D0 + JUNK(4))/ROB            
      DBCOOR(2)= DBCOOR(1)                                                      
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                     
C        NOW CALCULATE THE COLLINEAR ENERGY                                     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                     
C                                                                               
C          STORE R(3) IN RACTMP THEN SET R(3) TO THE COLLINEAR                  
C          GEOMETRY AND CALL POT                                                
C                                                                               
      RACTMP = R(3)                                                             
      R(3) = R(2) + R(1)                                                        
C                                                                               
      CALL POT2                                                                 
C                                                                               
      R(3) = RACTMP                                                             
C                                                                               
      ENGYGS = ENGYGS + BCOOR                                                   
C                                                                               
      CALL EUNITZERO                                                            
      IF(NDER.NE.0) THEN                                                        
         DEGSDR(1) = DEGSDR(1) + DEGSDR(3)                                      
         DEGSDR(2) = DEGSDR(2) + DEGSDR(3)                                      
         DEGSDR(3) = DAC                                                        
         DO 300 IT = 1,2                                                        
              DEGSDR(IT) = DBCOOR(IT) + DEGSDR(IT)                              
  300    CONTINUE                                                               
         CALL RTOCART                                                           
         CALL DEDCOU                                                            
      ENDIF                                                                     
C                                                                               
600   FORMAT (/,2X,T5,'PREPOT has been called for the H3 ',                     
     *                'potential energy surface TK',                            
     *       //,2X,T5,'Potential energy surface parameters:',                   
     *        /,2X,T5,'Bond', T47, 'H-H', T58, 'H-H', T69, 'H-H',               
     *        /,2X,T5,'Dissociation energies (kcal/mol):',                      
     *        T44, F10.5, T55, F10.5, T66, F10.5,                               
     *        /,2X,T5,'Equilibrium bond lengths (Angstroms):',                  
     *        T44, F10.5, T55, F10.5, T66, F10.5,                               
     *        /,2X,T5,'Morse beta parameters (Angstroms**-1):',                 
     *        T44, F10.5, T55, F10.5, T66, F10.5)                               
610   FORMAT(/,2X,T5,'Pauling parameter (Angstroms):',T44,F10.5,                
     *       /,2X,T5,'Gamma for the bending correction:',T44,F10.5,             
     *       /,2X,T5,'RBB (Angstroms):',T44,F10.5)                              
 900  FORMAT(/,2X,T5,13HNASURF(1,1) =,I5,                                       
     *       /,2X,T5,24HThis value is unallowed.                                
     *       /,2X,T5,31HOnly gs surface=>NASURF(1,1)=1 )                        
910   FORMAT(/, 2X,'POT has been called with NDER = ',I5,                       
     *       /, 2X, 'This value of NDER is not supported in this ',             
     *              'version of the potential.')                                
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE PREPOT2                                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
C                                                                               
      DOUBLE PRECISION NAB,NBC,NAB1,NBC1,JUNK                                   
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
      PARAMETER (HARTRE = 627.5095D0)                                           
C                                                                               
      LOGICAL RABL,RBCL                                                         
C                                                                               
C         COORDINATES,POTENTIAL ENERGY, AND DERIVATIVES ARE                     
C                                                                               
      COMMON /PRECM/  DE,BETA,A,GAMMA,RBB,RHH,                                  
     +                NAB,NBC,NAB1,NBC1,                                        
     +                DBCOOR(2),JUNK(4),REQ(3)                                  
      COMMON /PRE2CM/ ALFD,A2,SN,D,ALFN,XN,RNMSN,                               
     +                RN,RNMXN,ALFDIF,R2T,L                                     
C                                                                               
C  THIS SUBPROGRAM SHOULD BE CALLED                                             
C  ONCE BEFORE THE FIRST CALL TO POT                                            
C                                                                               
      RN=RNMSN+SN                                                               
      RNMXN=RN-XN                                                               
      ALFDIF=ALFD-ALFN                                                          
      R2T=1.4142136D0*RNMSN                                                     
      D = D /HARTRE                                                             
C                                                                               
C    Set the values of the classical energy in the three asymptotic valleys     
C                                                                               
             EASYAB = D                                                         
             EASYBC = EASYAB                                                    
             EASYAC = EASYAB                                                    
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
C     ENTRY POT2                                                                
C                                                                               
      SUBROUTINE POT2                                                           
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
C                                                                               
      DOUBLE PRECISION NAB,NBC,NAB1,NBC1,JUNK                                   
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
      PARAMETER (HARTRE = 627.5095D0)                                           
C                                                                               
      LOGICAL RABL,RBCL                                                         
C                                                                               
C         COORDINATES,POTENTIAL ENERGY, AND DERIVATIVES ARE                     
C         PASSED THROUGH COMMON /PT1CM/ R(3), ENGYGS, DEGSDR(3)                 
C                                                                               
      COMMON /PRECM/  DE,BETA,A,GAMMA,RBB,RHH,                                  
     +                NAB,NBC,NAB1,NBC1,                                        
     +                DBCOOR(2),JUNK(4),REQ(3)                                  
      COMMON /PRE2CM/ ALFD,A2,SN,D,ALFN,XN,RNMSN,                               
     +                RN,RNMXN,ALFDIF,R2T,L                                     
C                                                                               
C THIS VERSION WORKS ONLY FOR L .LT. 4                                          
C                                                                               
      RABL=R(1).GT.RN                                                           
      RBCL=R(2).GT.RN                                                           
      IF (RABL.OR.RBCL) GO TO 12                                                
      RNMX=RN-R(1)                                                              
      RNMY = RN - R(2)                                                          
      THETA= ATAN(RNMY/RNMX)                                                    
      ST = SIN(THETA)                                                           
      CT = COS(THETA)                                                           
      STT = 2.D0*ST*CT                                                          
      CTT = CT*CT - ST*ST                                                       
      LM4=L-4                                                                   
    3 STTL=STT**L                                                               
      STT4=STTL*STT**(-LM4)                                                     
    6 GAM=D*(1.D0-A2*STTL)                                                      
      ALF=ALFN+ALFDIF*STT4                                                      
      BRACK=R2T-RNMXN                                                           
      RCRC=2.D0*RN*RN+R(1)*R(1)+R(2)*R(2)-2.D0*RN*(R(1)+R(2))                   
      SQUIG= BRACK*STT4+RNMXN-SQRT(RCRC)                                        
C                                                                               
C ENERGY AND DERIVATIVE INFORMATION                                             
C                                                                               
      T1 = STTL*CTT*3.D0*A2*D                                                   
       DGAMX = -T1/RNMX                                                         
       DGAMY = T1/RNMY                                                          
       T2 = STT4*CTT*4.D0                                                       
       T1 = T2*ALFDIF                                                           
       DALFX = T1/RNMX                                                          
       DALFY = -T1/RNMY                                                         
       T1 = T2*BRACK                                                            
      DSQIGX = CT + T1/RNMX                                                     
       DSQIGY = ST - T1/RNMY                                                    
      GO TO 13                                                                  
   12 GAM=D                                                                     
      ALF=ALFN                                                                  
      IF(RABL) SQUIG = R(2)-XN                                                  
      IF(RBCL) SQUIG = R(1) - XN                                                
      IF(RABL.AND.RBCL) SQUIG = SQUIG + R(2) - RN                               
   13 CONTINUE                                                                  
      T1 = EXP(-ALF*SQUIG)                                                      
      QUAN = 1.D0-T1                                                            
      T2 = QUAN*QUAN - 1.D0                                                     
      ENGYGS = GAM*T2 + EZERO(1)                                                
      IF(NDER.EQ.0) RETURN                                                      
C                                                                               
C EXCLUSIVELY DERIVATIVES                                                       
C                                                                               
      IF(RABL.OR.RBCL) GO TO 20                                                 
      DEGSDR(1) = DGAMX*T2                                                      
      DEGSDR(2) = DGAMY * T2                                                    
      T2 = T1*QUAN*2.D0*GAM                                                     
      DEGSDR(1) = DEGSDR(1) + T2*(ALF*DSQIGX + SQUIG*DALFX)                     
      DEGSDR(2) = DEGSDR(2) + T2*(ALF*DSQIGY + SQUIG*DALFY)                     
      DEGSDR(3) = 0.D0                                                          
      RETURN                                                                    
   20 T2 = 2.D0*ALF*GAM*T1*QUAN                                                 
      DEGSDR(1) = 0.D0                                                          
      DEGSDR(2) = 0.D0                                                          
      IF(RABL) DEGSDR(2) = T2                                                   
      IF(RBCL) DEGSDR(1) = T2                                                   
      DEGSDR(3) = 0.D0                                                          
      RETURN                                                                    
      END                                                                       
C                                                                               
C*****                                                                          
C                                                                               
         BLOCK DATA PTPACM                                                      
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                    
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
         DOUBLE PRECISION NAB,NBC,NAB1,NBC1,JUNK                                
         COMMON /PRECM/  DE,BETA,A,GAMMA,RBB,RHH,                               
     +                   NAB,NBC,NAB1,NBC1,                                     
     +                   DBCOOR(2),JUNK(4),REQ(3)                               
         COMMON /PRE2CM/ ALFD,A2,SN,D,ALFN,XN,RNMSN,                            
     +                   RN,RNMXN,ALFDIF,R2T,L                                  
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
C   POTENTIAL ENERGY PARAMETERS                                                 
C                                                                               
C   DE    = dissociation energy in kcal/mol                                     
C   BETA  = Morse beta in reciprocal Angstroms                                  
C   A     = Pauling parameter in Angstroms                                      
C   GAMMA = gamma for the bend correction                                       
C           gamma is dimensionless and ranges from 0.5 to 0.65                  
C   RBB   = parameter used in the scaling calculation for the                   
C           bending correction                                                  
C                                                                               
         DATA DE /109.47D0/                                                     
         DATA BETA /1.97354/                                                    
         DATA REQ /0.741287D0, 0.741287D0, 0.741287D0/                          
         DATA A /0.26D0/                                                        
         DATA GAMMA /0.5D0/                                                     
         DATA RBB /0.741270D0/                                                  
C                                                                               
C        DATA HAS ENERGY IN UNITS OF KCAL/MOLE  IN DATA STATEMENT               
C        THIS IN CONVERTED TO HARTREES.  LENGTH IS IN BOHR.                     
C                                                                               
C   DATA FOR POT2                                                               
C                                                                               
      DATA ALFD,A2,SN,D,ALFN,XN                                                 
     +     /0.70152D0,0.08939D0,1.765D0,109.47D0,1.04435D0,1.40083D0/           
      DATA L,RNMSN                                                              
     +     /3,1.85630D0/                                                        
C                                                                               
         END                                                                    
C                                                                               
C*****                                                                          
