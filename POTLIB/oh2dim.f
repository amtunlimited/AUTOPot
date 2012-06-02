C                                                                               
      SUBROUTINE PREPOT                                                         
C                                                                               
C   System:          OH2                                                        
C   Functional form: Anti-Morse Bends plus Bowman's collinear spline            
C   Common name:     OH2DIM                                                     
C   Reference:       A. Whitlock, J. T. Muckerman, E. R. Fisher,
C                    Report from Research Institute for Engineering Siences,
C                    Wayne State University, 1976.
C                    A. F. Wagner, G. C. Schatz, and J. M. Bowman,
C                    J. Chem. Phys., Vol. 74, 4960(1981);
C                    J. M. Bowman and K. T. Lee in "Potential Energy Surfaces
C                    and Dynamics Calculations" D. G. Truhlar, Ed., Plenum, 1981,
C                    p. 359.
C                                                                               
C   PREPOT must be called once before any calls to POT.                         
C   The control flags for the potential, such as IPRT, are                      
C   initialized in the block data subprogram PTPACM.                            
C   Coordinates, potential energy, and derivatives are passed                   
C   The potential energy in the three asymptotic valleys are                    
C   stored in the common block ASYCM:                                           
C                  COMMON /ASYCM/ EASYAB, EASYBC, EASYAC                        
C   The potential energy in the AB valley, EASYAB, is equal to the potential    
C   energy of the H "infinitely" far from the OH diatomic, with the             
C   OH diatomic at its equilibrium configuration.  Similarly, the terms         
C   EASYBC and EASYAC represent the H2 and the OH asymptotic valleys,           
C   respectively.                                                               
C   All the information passed through the common blocks PT1CM and ASYCM        
C   is in Hartree atomic units.                                                 
C                                                                               
C   This potential is written such that:                                        
C                  R(1) = R(O-H)                                                
C                  R(2) = R(H-H)                                                
C                  R(3) = R(H-O)                                                
C   The zero of energy is defined at O "infinitely" far from the H2 diatomic.   
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
C          GENERAL INFORMATION                                                  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                                  
C                                                                               
C   POTA.FOR HAS ANTI-MORSE BEND PLUS BOWMAN'S COLLINEAR SPLINE ROUTINE.        
C   ANTI-MORSE BENDING POTENTIAL IS ADDED TO COLLINEAR POTENTIAL.               
C   THIS SUBROUTINE CAN BE USED TO ADD A BENDING CORRECTION TO                  
C   ANY COLLINEAR POTENTIAL SUBROUTINE.                                         
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCC                                                     
C          DEFINE VARIABLES                                                     
CCCCCCCCCCCCCCCCCCCCCCCCCCC                                                     
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      REAL*8 NAB,NAB1,NBC,NBC1,JUNK,HRTREE,BOHR                                 
      DIMENSION JUNK(4)                                                         
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
C Conversion constant from Angstroms to Bohr (Angstrom/Bohr)                    
C                                                                               
         PARAMETER (BOHR = .52917706D0)                                         
C                                                                               
C Conversion constant from Hartree to kcal (kcal/Hartree)                       
C                                                                               
         PARAMETER (HRTREE = 627.5095D0)                                        
C                                                                               
         COMMON /NSTEP/  STEP                                                   
         COMMON /POTCM/  XBC(3),XAB(3),CBC(3),ABC(3),CAB(3),                    
     +                   AAB(3),PHI(40),BSP(40,3),                              
     +                   DE,BETA,REQ(3),A,GAMMA,RBB,RHH,                        
     +                   R1S,R2S,VSP,DBCOOR(2),NPHI                             
C                                                                               
      IF(NATOMS.GT.25) THEN                                                     
         WRITE(NFLAG(18),1111)                                                  
 1111    FORMAT(2X,'STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS')            
         STOP                                                                   
      END IF                                                                    
C                                                                               
      WRITE(NFLAG(18), 600)                                                     
600   FORMAT (/,2X,T5,'PREPOT2 has been called for the OH2 ',                   
     *                'potential energy surface DIM',                           
     *        /,2X,T5,'Parameters for the potential energy surface',/)          
C                                                                               
      CALL PREPOT2                                                              
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                        
C      BENDING POTENTIAL PARAMETERS                                             
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                        
C                                                                               
      WRITE(NFLAG(18),3)                                                        
 3    FORMAT(/,2X,T5,'Parameters for the bending potential ',                   
     +       'initialized in BLOCK DATA')                                       
C                                                                               
C   Echo the potential parameters to the file linked to UNIT NFLAG(18)          
C                                                                               
      WRITE (NFLAG(18),25) (REQ(IT), IT = 1,3)                                  
      WRITE(NFLAG(18),10) DE                                                    
      WRITE(NFLAG(18),21) BETA                                                  
      WRITE (NFLAG(18),30) A                                                    
      WRITE (NFLAG(18),31) GAMMA                                                
      WRITE(NFLAG(18),35) RBB                                                   
C                                                                               
25    FORMAT(2X,T10,'Equilibrium bond distances (Angstroms):',                  
     *       /,2X,T15,3(1PE20.10,1X))                                           
10    FORMAT(2X,T10,'Dissociation energies (kcal/mol):',                        
     *          T50,1PE20.10)                                                   
21    FORMAT(2X,T10,'Morse Beta parameters (Angstroms**-1):',                   
     *          T50,1PE20.10)                                                   
35    FORMAT(2X,T10,'RBB (Angstroms):',                                         
     *          T50,1PE20.10)                                                   
30    FORMAT(2X,T10,'Pauling parameter (Angstroms):',                           
     *          T50,1PE20.10)                                                   
31    FORMAT(2X,T10,'Gamma (unitless):',                                        
     *          T50,1PE20.10)                                                   
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                                 
C       CONVERT TO ATOMIC UNITS                                                 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                                 
C                                                                               
      DO 50 IT = 1,3                                                            
      REQ(IT) = REQ(IT)/BOHR                                                    
   50 CONTINUE                                                                  
      BETA = BETA * BOHR                                                        
      DE = DE/HRTREE                                                            
      A = A/BOHR                                                                
      RHH = 1.4007977D0                                                         
      RBB = RBB /BOHR                                                           
C                                                                               
CCCCCCCCCCCC                                                                    
C                                                                               
C                                                                               
      EZERO(1)=VSP                                                              
C                                                                               
       DO I=1,5                                                                 
          REF(I) = ' '                                                          
       END DO                                                                   
C                                                                               
       REF(1)='A. Whitlock, J. T. Muckerman, E. R. Fisher'
       REF(2)='Report Res. Inst. for Eng. Sciences., Wayne State University, 1976.'
       REF(3)='A.F. Wagner, G.C. Schatz, J.M. Bowman, J. Chem. Phys. 74, 4960(1981).
       REF(4)='J.M. Bowman and K.T. Lee in "Potential Energy Surfaces and Dynamics' 
       REF(5)='Calculations" D. G. Truhlar, Ed., Plenum, 1981, p. 359.'
C                                                                               
      INDEXES(1) = 8                                                            
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
C                                                                               
CCCCCCCCCCCC                                                                    
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                 
C        MOVE TO POT SUBROUTINE                                                 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                 
C                                                                               
      SUBROUTINE POT                                                            
C                                                                               
C   The potential energy in the AB valley, EASYAB, is equal to the potential    
C   energy of the H "infinitely" far from the OH diatomic, with the             
C   OH diatomic at its equilibrium configuration.  Similarly, the terms         
C   EASYBC and EASYAC represent the H2 and the OH asymptotic valleys,           
C   respectively.                                                               
C                                                                               
C   This potential is written such that:                                        
C                  R(1) = R(O-H)                                                
C                  R(2) = R(H-H)                                                
C                  R(3) = R(H-O)                                                
C   The zero of energy is defined at O "infinitely" far from the H2 diatomic.   
C                                                                               
CCCCCCCCCCCCCCCC                                                                
C      ENTRY POT                                                         8/17R79
CCCCCCCCCCCCCCCC                                                                
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      REAL*8 NAB,NAB1,NBC,NBC1,JUNK,HRTREE,BOHR                                 
      DIMENSION JUNK(4)                                                         
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
C Conversion constant from Angstroms to Bohr (Angstrom/Bohr)                    
C                                                                               
         PARAMETER (BOHR = .52917706D0)                                         
C                                                                               
C Conversion constant from Hartree to kcal (kcal/Hartree)                       
C                                                                               
         PARAMETER (HRTREE = 627.5095D0)                                        
C                                                                               
         COMMON /NSTEP/  STEP                                                   
         COMMON /POTCM/  XBC(3),XAB(3),CBC(3),ABC(3),CAB(3),                    
     +                   AAB(3),PHI(40),BSP(40,3),                              
     +                   DE,BETA,REQ(3),A,GAMMA,RBB,RHH,                        
     +                   R1S,R2S,VSP,DBCOOR(2),NPHI                             
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
      VAC = GAMMA *1.0D-14 * DE * JUNK(3)* (1.D0 + 0.5D0*JUNK(3))               
C                                                                               
C        CALCULATE V(R(1)+R(2))                                                 
C                                                                               
      JUNK(4) = EXP(BETA*(REQ(3)-R(2)-R(1))/ROB)                                
      VABBC= GAMMA *1.0D-14 * DE * JUNK(4) * (1.D0 + 0.5D0*JUNK(4))             
C                                                                               
C        HERE IS THE BENDING CORRECTION                                         
C                                                                               
      BCOOR = (VAC - VABBC)*1.0D14                                              
C                                                                               
C          WHILE IT'S CONVENIENT, CALCULATE THE DERIVATIVE                      
C          OF BCOOR WITH RESPECT TO R(1), R(2), AND R(3).                       
C                                                                               
C        CALCULATE DAC (THE DERIVATIVE WITH RESPECT TO R(3))                    
C                                                                               
      IF (NDER .EQ. 1) THEN                                                     
          DAC = -GAMMA*DE*BETA*JUNK(3)*(1.D0+JUNK(3))/ROB                       
          DBCOOR(1) = GAMMA * DE * BETA * JUNK(4) *                             
     *                (1.D0 + JUNK(4))/ROB                                      
          DBCOOR(2)= DBCOOR(1)                                                  
      ENDIF                                                                     
C                                                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                     
C        NOW CALCULATE THE COLLINEAR ENERGY                                     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                     
C                                                                               
C          STORE R(3) IN RACTMP THEN SET R(3) TO THE COLLINEAR                  
C          GEOMETRY AND CALL POT2                                               
C                                                                               
      RACTMP = R(3)                                                             
      R(3) = R(2) + R(1)                                                        
C                                                                               
      CALL POT2                                                         8/17R79 
C                                                                               
      IF (NDER .EQ. 1) THEN                                                     
          DEGSDR(1) = DEGSDR(1) + DEGSDR(3)                                     
          DEGSDR(2) = DEGSDR(2) + DEGSDR(3)                                     
          DEGSDR(3) = DAC                                                       
      ENDIF                                                                     
      R(3) = RACTMP                                                             
C                                                                               
      ENGYGS = ENGYGS + BCOOR                                                   
C                                                                               
      IF (NDER .EQ. 1) THEN                                                     
          DO 300 IT = 1,2                                                       
                 DEGSDR(IT) = DBCOOR(IT) + DEGSDR(IT)                           
300       CONTINUE                                                              
      ENDIF                                                                     
C                                                                               
 900  FORMAT(/,2X,T5,13HNASURF(1,1) =,I5,                                       
     *       /,2X,T5,24HThis value is unallowed.                                
     *       /,2X,T5,31HOnly gs surface=>NASURF(1,1)=1 )                        
910   FORMAT(/, 2X,'POT has been called with NDER = ',I5,                       
     *       /, 2X, 'This value of NDER is not allowed in this ',               
     *              'version of the potential.')                                
C                                                                               
      CALL EUNITZERO                                                            
      IF(NDER.NE.0) THEN                                                        
         CALL RTOCART                                                           
         CALL DEDCOU                                                            
      ENDIF                                                                     
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE PREPOT2                                                        
C                                                                               
C   POTENTIAL FUNCTION FOR RATE1D SERIES                                        
C   RMO-SPLINE TYPE FIT WITH BEBO TYPE CONTINUATION INTO THE ASYMBTOTIC REGION  
C   POTENTIAL ROUTINE INTERACTS WITH REST OF PROGRAM WITH UNITS OF EV/BOHR      
C   POTENTIAL ROUTINE ASSUMES SPLINE INFO IS IN UNITS OF KCAL/ANGSTROM          
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      REAL*8 NAB,NAB1,NBC,NBC1,JUNK,HRTREE,BOHR                                 
      DIMENSION JUNK(4)                                                         
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
C Conversion constant from Angstroms to Bohr (Angstrom/Bohr)                    
C                                                                               
      PARAMETER (BOHR = .52917706D0)                                            
C                                                                               
C Conversion constant from Hartree to kcal (kcal/Hartree)                       
C                                                                               
      PARAMETER (HRTREE = 627.5095D0)                                           
C                                                                               
      PARAMETER (DETRAD=.0174532D0)                                             
C                                                                               
      COMMON /NSTEP/  STEP                                                      
      COMMON /POTCM/  XBC(3),XAB(3),CBC(3),ABC(3),CAB(3),                       
     +                AAB(3),PHI(40),BSP(40,3),                                 
     +                DE,BETA,REQ(3),A,GAMMA,RBB,RHH,                           
     +                R1S,R2S,VSP,DBCOOR(2),NPHI                                
C                                                                               
      DIMENSION       V(4),CSP(48,3),SCRAP(48)                                  
C                                                                               
      CHARACTER*80 TITLE                                                        
C                                                                               
      TITLE=' RMOFIT to DIM surface with BEBO continuations '                   
C                                                                               
      WRITE(NFLAG(18),845)TITLE                                                 
      WRITE(NFLAG(18),125)STEP                                                  
      WRITE(NFLAG(18),111)XAB(1), XBC(1), XAB(1),                               
     *               XAB(2), XBC(2), XAB(2),                                    
     *               XAB(3), XBC(3), XAB(3)                                     
      WRITE(NFLAG(18),105)R1S,R2S,VSP                                           
      WRITE(NFLAG(18),113)(PHI(I),(BSP(I,J),J = 1,3),I = 1,NPHI)                
      WRITE(NFLAG(18),117)(CBC(I),ABC(I),I = 1,3),                              
     +                    (CAB(I),AAB(I),I = 1,3)                               
C                                                                               
C   Initialize the energy in the asymptotic valleys                             
C                                                                               
      EASYAB = XAB(1)/HRTREE                                                    
      EASYBC = XBC(1)/HRTREE                                                    
      EASYAC = EASYAB                                                           
C                                                                               
2     FORMAT(F12.10)                                                            
105   FORMAT(/,2X,T5,'R1S, R2S, VSP = ',3(F10.5,1X))                            
111   FORMAT (/,2X,T5,'Bond', T47, 'O-H', T58, 'H-H', T69, 'H-O',               
     *        /,2X,T5,'Dissociation energies (kcal/mol):',                      
     *        T44, F10.5, T55, F10.5, T66, F10.5,                               
     *        /,2X,T5,'Equilibrium bond lengths (Angstroms):',                  
     *        T44, F10.5, T55, F10.5, T66, F10.5,                               
     *        /,2X,T5,'Morse beta parameters (Angstroms**-1):',                 
     *        T44, F10.5, T55, F10.5, T66, F10.5)                               
113   FORMAT(/,2X,T24,'Phi',T39,'De',T59,'Re',T73,'Beta',/,                     
     *      (2X,T20,F10.5,T35,F10.5,T55,F10.5,T70,F10.5))                       
117   FORMAT(/,2X,T5,'For the asymptote: MO parameter value = ',                
     *               '(asymp. value)+C*EXP(-A*(RI-RIS))',                       
     *       /,2X,T10,'CBC(I); I=1,3: ',3(F15.5, 1X),                           
     *       /,2X,T10,'ABC(I); I=1,3: ',3(F15.5, 1X),                           
     *       /,2X,T10,'CAB(I); I=1,3: ',3(F15.5, 1X),                           
     *       /,2X,T10,'AAB(I); I=1,3: ',3(F15.5, 1X))                           
125   FORMAT(2X,T5,'Step size for the numerical derivative ',                   
     *               'calculation = ',1PE20.10)                                 
401   FORMAT(16I5)                                                              
403   FORMAT(3E20.7)                                                            
407   FORMAT(4E20.7)                                                            
843   FORMAT(A80)                                                               
845   FORMAT(2X,T5,//'Title card in potential input data file: ',//,A80)        
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE POT2                                                           
C                                                                               
C COCK SPLINES                                                                  
C                                                                               
C      ENTRY POT2                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
      REAL*8 NAB,NAB1,NBC,NBC1,JUNK,HRTREE,BOHR                                 
      DIMENSION JUNK(4)                                                         
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
      PARAMETER (HRTREE = 627.5095D0)                                           
      PARAMETER (DETRAD=.0174532D0)                                             
C                                                                               
      COMMON /NSTEP/  STEP                                                      
      COMMON /POTCM/  XBC(3),XAB(3),CBC(3),ABC(3),CAB(3),                       
     +                AAB(3),PHI(40),BSP(40,3),                                 
     +                DE,BETA,REQ(3),A,GAMMA,RBB,RHH,                           
     +                R1S,R2S,VSP,DBCOOR(2),NPHI                                
C                                                                               
      DIMENSION V(4),CSP(48,3),SCRAP(48)                                        
C                                                                               
      IC=1                                                                      
   10 GO TO (81,82,83,84),IC                                                    
  81  R1=R(1)                                                                   
      R2=R(2)                                                                   
      R3=R(3)                                                                   
      GO TO 85                                                                  
   82 R1=R(1) + STEP                                                            
      R2=R(2)                                                                   
      R3=R(3)                                                                   
      GO TO 85                                                                  
   83 R1=R(1)                                                                   
      R2=R(2) + STEP                                                            
      R3=R(3)                                                                   
      GO TO 85                                                                  
   84 R1=R(1)                                                                   
      R2=R(2)                                                                   
      R3=R(3) + STEP                                                            
 85   R1=R1*0.529177D0                                                          
      R2=R2*0.529177D0                                                          
      R3=R3*0.529177D0                                                          
      DR1=R1S - R1                                                              
      DR2=R2S - R2                                                              
      ANGLE=ATAN(DR1/DR2)                                                       
      ANGLE=ANGLE/DETRAD                                                        
      DO 500 K=1,3                                                              
      PPP=SPLINE(1,NPHI,PHI,BSP(1,K),CSP(1,K),SCRAP,ANGLE)                      
 500  CONTINUE                                                                  
C                                                                               
C SETUP ANY CONSTANTS                                                           
C                                                                               
      DBC = XBC(1)                                                              
      DERG = XBC(1)/23.061D0                                                    
      RERG = (R2S-XBC(2))/.529177D0                                             
      BETRG = XBC(3)*.529177D0                                                  
      DEPR = XAB(1)/23.061D0                                                    
      REPR = (R1S-XAB(2))/.529177D0                                             
      BETPR = XAB(3)*.529177D0                                                  
C                                                                               
C ENTRY INTO POTENTIAL CALCULATOR                                               
C                                                                               
      IF(R1.GT.R1S) GO TO 100                                                   
      IF(R2.GT.R2S) GO TO 102                                                   
C                                                                               
C INTERACTION REGION OF POLAR COORDINATES                                       
C                                                                               
      DR1=R1S-R1                                                                
      DR2=R2S-R2                                                                
      ANGLE=ATAN(DR1/DR2)                                                       
      ANGLE=ANGLE/DETRAD                                                        
      RR=SQRT(DR1*DR1+DR2*DR2)                                                  
      DIS=SPLINE(2,NPHI,PHI,BSP(1,1),CSP(1,1),SCRAP,ANGLE)                      
      REQ2=SPLINE(2,NPHI,PHI,BSP(1,2),CSP(1,2),SCRAP,ANGLE)                     
      BET=SPLINE(2,NPHI,PHI,BSP(1,3),CSP(1,3),SCRAP,ANGLE)                      
C      D(1)=SPLINE(3,NPHI,PHI,BSP(1,1),CSP(1,1),SCRAP,ANGLE)                    
C      D(2)=SPLINE(3,NPHI,PHI,BSP(1,2),CSP(1,2),SCRAP,ANGLE)                    
C30      D(1)=D(1)/627.5095D0*0.529177D0                                        
C      D(2)=D(2)/627.5095D0*0.529177D0                                          
C                                                                               
C 30   X = DIS*((1.D0-EXP(BET*(RR-REQ2)))**2.0D0 -1.D0) + VSP                   
 30   X = DIS*((1.D0-EXP(BET*(RR-REQ2)))**2.0D0 -1.D0) +                        
     +    EZERO(1) - (ANUZERO*627.5095D0)                                       
C       1                                                                       
      V(IC) = X                                                                 
      IC=IC+1                                                                   
      IF(IC.LT.5) GO TO 10                                                      
      DO 333 I=1,3                                                              
 333  V(I)=V(I)/627.5095D0                                                      
      IF (NDER .EQ. 1) THEN                                                     
          DEGSDR(1)=(V(2)-V(1))/STEP                                            
          DEGSDR(2)=(V(3)-V(1))/STEP                                            
          DEGSDR(3)=0.0D0                                                       
      ENDIF                                                                     
      ENGYGS=V(1)                                                               
   93 RETURN                                                                    
C                                                                               
C LARGE R1 ASYMPTOTIC REGION                                                    
C                                                                               
  100 IF(R2.GT.R2S) GO TO 104                                                   
      RR = R2S-R2                                                               
      SEP = R1-R1S                                                              
      DIS = XBC(1) + CBC(1)*EXP(-MIN(ABC(1)*SEP,25.D0))                         
      REQ2 = XBC(2) + CBC(2)*EXP(-MIN(ABC(2)*SEP,25.D0))                        
      BET = XBC(3) + CBC(3)*EXP(-MIN(ABC(3)*SEP,25.D0))                         
C      IF(R1.GT.10.D0)D(1)=1.0D1                                                
C      IF(R1.GT.10.D0)D(2)=1.0D-8                                               
C      IF(R1.GT.10.D0)GO TO 30                                                  
C      D(1)=-CBC(1)*ABC(1)*((1.D0-EXP(-BET*(RR-REQ2)))**2-1.D0)-                
C     +     DIS*2.D0*(1.D0-EXP(-BET*(RR-REQ2)))*EXP(-BET*(RR-REQ2))*            
C     +     ((-CBC(3)*ABC(3)*EXP(-ABC(3)*SEP)*(RR-REQ2) -                       
C     +     CBC(2)*ABC(2)*EXP(ABC(2)*SEP)*BET))*BET*(RR-REQ2)                   
C      D(2)=2.D0*DIS*(1.D0-EXP(-BET*(RR-REQ2)))*BET*(RR-REQ2)*                  
C     +     EXP(-BET*(RR-REQ2))*BET                                             
       GO TO 30                                                                 
C                                                                               
C LARGE R2 ASYMPTOTIC REGION                                                    
C                                                                               
  102 CONTINUE                                                                  
      RR = R1S-R1                                                               
      SEP = R2-R2S                                                              
      DIS = XAB(1) + CAB(1)*EXP(-MIN(AAB(1)*SEP,25.D0))                         
      REQ2 = XAB(2) + CAB(2)*EXP(-MIN(AAB(2)*SEP,25.D0))                        
      BET = XAB(3) + CAB(3)*EXP(-MIN(AAB(3)*SEP,25.D0))                         
C      IF(R2.GT.10.D0)D(2)=1.0D1                                                
C      IF(R2.GT.10.D0)D(1)=1.0D-8                                               
C      IF(R2.GT.10.D0)GO TO 30                                                  
C      D(1)=-CAB(1)*AAB(1)*((1.D0-EXP(-BET*(RR-REQ2)))**2-1.D0)-                
C     +     DIS*2.D0*(1.D0-EXP(-BET*(RR-REQ2)))*EXP(-BET*(RR-REQ2))*            
C     +     ((-CAB(3)*AAB(3)*EXP(-AAC(3)*SEP)*(RR-REQ2) -                       
C     +     CAB(2)*AAB(2)*EXP(AAB(2)*SEP)*BET))*BET*(RR-REQ2)                   
C      D(2)=2.D0*DIS*(1.D0-EXP(-BET*(RR-REQ2)))*BET*(RR-REQ2)*                  
C     +     EXP(-BET*(RR-REQ2))*BET                                             
       GO TO 30                                                                 
C                                                                               
C LARGE R1,R2 3 BODY BREAK-UP REGION                                            
C                                                                               
  104 CONTINUE                                                                  
      J=J+1                                                                     
      ENGYGS=XBC(1)/627.5095D0                                                  
      IF (NDER .EQ. 1) THEN                                                     
          DEGSDR(1)=1.0D-9                                                      
          DEGSDR(2)=1.0D-9                                                      
          DEGSDR(3)=0.0D0                                                       
      ENDIF                                                                     
      GO TO 93                                                                  
      END                                                                       
C                                                                               
      FUNCTION SPLINE(ISW,NN,X,Y,C,D,XPT)                                       
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
C      DIMENSION X(*),Y(*),C(*),D(*)                                            
      DIMENSION X(NN+2),Y(NN+2),C(NN+2),D(NN+2)                                 
C                                                                               
C                                                                               
C  THIS IS A SUBROUTINE FOR FITTING DATA WITH A CUBIC SPLINE                    
C  POLYNOMIAL AND EVALUATING THAT POLYNOMIAL AT A GIVEN POINT                   
C  OR ITS DERIVATIVE AT A GIVEN POINT                                           
C                                                                               
C  CALLING SEQUENCE .......                                                     
C                                                                               
C     ISW ... CONTROL OPTION                                                    
C         ISW=1  IF A CUBIC SPLINE IS TO BE FITTED TO THE SET OF KNOTS          
C                DEFINED BY THE ARRAYS X AND Y.  THE SPLINE COEFFICIENTS        
C                ARE STORED IN THE ARRAY C.                                     
C         ISW=2  IF THE SPLINE DEFINED BY THE COEFFICIENT ARRAY 'C' IS          
C                TO BE EVALUATED (INTERPOLATED) AT THE POINT DEFINED BY         
C                THE PARAMETER 'XPT'.                                           
C         ISW=3  AS IN ISW=2, ONLY THE DERIVATIVE/3.D0 IS ALSO CALCULATED AT XPT
C         ISW=4  THE DERIVATIVE CALCULATED BY THE LAST USE OF SPLINE WITH ISW=3 
C                IS RETURNED.                                                   
C                                                                               
C     NN ... THE NUMBER OF KNOTS (DATA POINTS) TO WHICH THE SPLINE IS TO        
C            BE FITTED                                                          
C                                                                               
C     X,Y ... THE ARRAYS DEFINING THE KNOTS.  THE X-VALUES MUST BE IN           
C             INCREASING ORDER.  THE ARRAYS MUST BE DIMENSIONED AT LEAST        
C             NN.                                                               
C                                                                               
C     C ... THE ARRAY THAT CONTAINS THE CUBIC SPLINE COEFFICIENTS.              
C           MUST BE DIMENSIONED AT LEAST NN+2 .                                 
C                                                                               
C     D ... A WORK SPACE.  MUST BE DIMENSIONED AT LEAST NN+2 .                  
C                                                                               
C     XPT ... THE POINT AT WHICH THE INTERPOLATION IS DESIRED (IF ISW IS        
C              SET TO 2).  THE VALUE OF SPLINE IS SET TO THE                    
C              INTERPOLATED VALUE.                                              
C                                                                               
C *****  USER NOTES  *****                                                      
C                                                                               
C     INTERPOLATION INVOLVES AT LEAST TWO STEPS .......                         
C                                                                               
C       A.  CALL SPLINE WITH THE KNOTS.  THIS SETS UP THE                       
C           COEFFICIENT ARRAY C.                                                
C           EG.  DUMY=SPLINE(1,NN,X,Y,C,D,XPT)                                  
C                                                                               
C       B.  CALL SPLINE WITH THE ARRAY C WHICH WAS DEFINED BY THE               
C           PREVIOUS CALL AND WILL BE USED TO FIND THE VALUE AT THE             
C           POINT 'XPT' .                                                       
C           EG.   VALUE=SPLINE(2,NN,X,Y,C,D,XPT)                                
C                                                                               
C     STEP 'A' NEED BE EXECUTED ONLY ONCE FOR A GIVEN SET OF KNOTS.             
C     STEP B MAY BE EXECUTED AS MANY TIMES AS NECESSARY.                        
C                                                                               
C     Output from this subprogram is written to unit 6.                         
C                                                                               
2     N=NN                                                                      
      NP1=N+1                                                                   
      NP2=N+2                                                                   
      Z=XPT                                                                     
24    GO TO (4,5,6,7),ISW                                                       
4     C(1)=Y(1)                                                                 
      D(1)=1.0D0                                                                
      C(NP1)=0.0D0                                                              
      D(NP1)=0.0D0                                                              
      C(NP2)=0.0D0                                                              
      D(NP2)=0.0D0                                                              
      DO 41 I=2,N                                                               
      C(I)=Y(I)-Y(1)                                                            
41    D(I)=X(I)-X(1)                                                            
      DO 410 I=3,NP2                                                            
      IF(D(I-1).NE.0)GO TO 43                                                   
      WRITE(NFLAG(18),1001)                                                     
      STOP 'SPLINE 1'                                                           
43    PIVOT=1.0D0/D(I-1)                                                        
      IF(I.GE.NP2)GO TO 45                                                      
      SUPD=X(I-1)-X(I-2)                                                        
      IF(SUPD.GE.0)GO TO 44                                                     
      WRITE(NFLAG(18),1000)                                                     
      STOP 'SPLINE 2'                                                           
44    SUPD=SUPD*SUPD*SUPD                                                       
      GO TO 46                                                                  
45    SUPD=1.0D0                                                                
46    DFACT=SUPD*PIVOT                                                          
      CFACT=C(I-1)*PIVOT                                                        
      IF(I.GT.N)GO TO 48                                                        
      DO 47 J=I,N                                                               
      V=X(J)-X(I-2)                                                             
      C(J)=C(J)-D(J)*CFACT                                                      
47    D(J)=V*V*V-D(J)*DFACT                                                     
48    CONTINUE                                                                  
      IF(I.GE.NP2)GO TO 49                                                      
      C(NP1)=C(NP1)-D(NP1)*CFACT                                                
      D(NP1)=1.0D0-D(NP1)*DFACT                                                 
49    C(NP2)=C(NP2)-D(NP2)*CFACT                                                
410   D(NP2)=X(I-2)-D(NP2)*DFACT                                                
      DO 411 I=1,N                                                              
      J=NP2-I                                                                   
      IF(J.NE.NP1)GO TO 413                                                     
      V=1.0D0                                                                   
      GO TO 414                                                                 
413   V=X(J)-X(J-1)                                                             
      V=V*V*V                                                                   
414   IF(D(J+1).NE.0)GO TO 415                                                  
      WRITE(NFLAG(18),1001)                                                     
      STOP 'SPLINE 3'                                                           
415   C(J+1)=C(J+1)/D(J+1)                                                      
411   C(J)=C(J)-C(J+1)*V                                                        
      IF(D(2).NE.0)GO TO 416                                                    
      WRITE(NFLAG(18),1001)                                                     
      STOP 'SPLINE 4'                                                           
416   C(2)=C(2)/D(2)                                                            
      RETURN                                                                    
5     SPLINE=C(1)+C(2)*(Z-X(1))                                                 
      DO 51 I=1,N                                                               
      V=Z-X(I)                                                                  
      IF(V.GT.0)GO TO 51                                                        
      RETURN                                                                    
51    SPLINE=SPLINE+C(I+2)*V*V*V                                                
      RETURN                                                                    
    6 CONTINUE                                                                  
      SPLINE=C(1)+C(2)*(Z-X(1))                                                 
      DERIV = C(2)/3.D0                                                         
      DO 53 I = 1,N                                                             
      V=Z-X(I)                                                                  
      IF(V.LE.0) RETURN                                                         
      V2 = V*V                                                                  
      SPLINE = SPLINE + C(I+2)*V2*V                                             
      DERIV = DERIV + C(I+2)*V2                                                 
   53 CONTINUE                                                                  
      RETURN                                                                    
    7 CONTINUE                                                                  
      SPLINE = 3.D0*DERIV                                                       
      RETURN                                                                    
1000  FORMAT(1X,5X,'***** ERROR IN SPLINE ... UNORDERED X-VALUES                
     1*****')                                                                   
1001  FORMAT(1X,5X,' ***** ERROR IN SPLINE ... DIVIDE FAULT *****')             
      END                                                                       
C                                                                               
C*****                                                                          
C                                                                               
         BLOCK DATA PTPACM                                                      
         IMPLICIT REAL*8 (A-H,O-Z)                                              
C                                                                               
C         REAL*8 NAB,NAB1,NBC,NBC1,JUNK,HRTREE,BOHR                             
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
         COMMON /NSTEP/  STEP                                                   
         COMMON /POTCM/  XBC(3),XAB(3),CBC(3),ABC(3),CAB(3),                    
     +                   AAB(3),PHI(40),BSP(40,3),                              
     +                   DE,BETA,REQ(3),A,GAMMA,RBB,RHH,                        
     +                   R1S,R2S,VSP,DBCOOR(2),NPHI                             
C                                                                               
C   Initialize the flags and the I/O unit numbers for the potential             
C                                                                               
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
         DATA STEP /0.000000008d0/                                              
         DATA NPHI /17/                                                         
         DATA R1S /0.2212079d+01/                                               
         DATA R2S /0.2215309d+01/                                               
         DATA VSP /0.1094890d+03/                                               
         DATA XBC / 0.1094890d+03, 0.1473409d+01, 0.1941998d+01/                
         DATA XAB / 0.1056400d+03, 0.1241479d+01, 0.2304798d+01/                
         DATA CBC /-0.4815674d-01, 0.3814697d-05, 0.4415512d-03/                
         DATA ABC / 0.7873216d+01, 0.7873216d+01, 0.7873216d+01/                
         DATA CAB /-0.1482804d+01,-0.2084732d-02,-0.3727913d-02/                
         DATA AAB / 0.1932070d+01, 0.1932070d+01, 0.1932070d+01/                
         DATA (PHI(I), I=1,17)                                                  
     +            / 0.0000000d+00, 0.2000000d+01,                               
     +              0.4000000d+01, 0.6000000d+01,                               
     +              0.1800000d+02, 0.2400000d+02,                               
     +              0.3000000d+02, 0.3400000d+02,                               
     +              0.3800000d+02, 0.4400000d+02,                               
     +              0.4800000d+02, 0.5400000d+02,                               
     +              0.6000000d+02, 0.7600000d+02,                               
     +              0.8600000d+02, 0.8800000d+02,                               
     +              0.9000000d+02/                                              
                                                                                
         DATA (BSP(I,1),I=1,17)                                                 
     +            / 0.1094409d+03, 0.1094212d+03,                               
     +              0.1094004d+03, 0.1093775d+03,                               
     +              0.1088412d+03, 0.1077417d+03,                               
     +              0.1055117d+03, 0.1030257d+03,                               
     +              0.9966266d+02, 0.9630693d+02,                               
     +              0.9637480d+02, 0.9811765d+02,                               
     +              0.1000138d+03, 0.1029822d+03,                               
     +              0.1038929d+03, 0.1040309d+03,                               
     +              0.1041572d+03/                                              
         DATA (BSP(I,2),I=1,17)                                                 
     +            / 0.1473413d+01, 0.1474257d+01,                               
     +              0.1476847d+01, 0.1481213d+01,                               
     +              0.1545280d+01, 0.1599182d+01,                               
     +              0.1655956d+01, 0.1681611d+01,                               
     +              0.1680401d+01, 0.1626349d+01,                               
     +              0.1572269d+01, 0.1486160d+01,                               
     +              0.1408802d+01, 0.1273705d+01,                               
     +              0.1241724d+01, 0.1239825d+01,                               
     +              0.1239394d+01/                                              
         DATA (BSP(I,3),I=1,17)                                                 
     +        / 0.1942440d+01, 0.1940545d+01,                                   
     +          0.1936212d+01, 0.1929558d+01,                                   
     +          0.1825750d+01, 0.1719964d+01,                                   
     +          0.1603891d+01, 0.1538857d+01,                                   
     +          0.1490752d+01, 0.1525683d+01,                                   
     +          0.1620584d+01, 0.1788983d+01,                                   
     +          0.1942842d+01, 0.2222638d+01,                                   
     +          0.2294156d+01, 0.2299114d+01,                                   
     +          0.2301070d+01/                                                  
C                                                                               
C       DISSOCIATION ENERGY IN KCAL/MOLE                                        
C                                                                               
         DATA DE /106.56d0/                                                     
C                                                                               
C       MORSE BETAS--RECIPROCAL ANGSTROMS                                       
C                                                                               
         DATA BETA /2.07942d0/                                                  
C                                                                               
C       THE EQUILIBRIUM DISTANCES IN ANGSTROMS                                  
C                                                                               
         DATA REQ /0.96966d0,0.74144d0,0.96966d0/                               
C                                                                               
C       PAULING PARAMETER IN ANGSTROMS                                          
C                                                                               
         DATA A /0.26d0/                                                        
C                                                                               
C       GAMMA FOR THE BEND CORRECTION                                           
C       GAMMA IS DIMENSIONLESS AND RANGES FROM .5 TO .65                        
C                                                                               
         DATA GAMMA / 0.4131d0/                                                 
C                                                                               
C       RBB IN ANGSTROMS--THIS NUMBER IS USED IN A SCALING                      
C       CALCULATION FOR THE BENDING CORRECTION                                  
C                                                                               
         DATA RBB / 0.74144d0/                                                  
C                                                                               
         END                                                                    
C                                                                               
C*****                                                                          
