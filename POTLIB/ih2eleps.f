C                                                                               
      SUBROUTINE PREPOT                                                         
C                                                                               
C   System:           IH2                                                       
C   Functional form:  Extended-LEPS (London-Erying-Polyani-Sato)                
C                                                                               
C   References for the potential parameters and the potential functional form:  
C   Potential Parameters: D. S. Perry, J. C. Polanyi, and C. W. Wilson, Jr.     
C                         Chem. Phys. 3, 317 (1974)                             
C   Functional form:      P. J. Kuntz, E. M. Nemth, J. C. Polanyi,              
C                         S. D. Rosner, and C. E. Young                         
C                         J. Chem. Phys. 44, 1168 (1966)                        
C                                                                               
C   Note: The parameters for this potential energy surface differ               
C         from those in C. A. Parr, Ph. D. Thesis, Calif. Inst. of              
C         Tech., Pasadena, 1969.                                                
C                                                                               
C   PREPOT must be called once before any calls to POT.                         
C   The potential parameters are included in the block data subprogram PTPACM.  
C   Coordinates, potential energy, and derivatives are passed                   
C   The potential energy in the three asymptotic valleys are                    
C   stored in the common block ASYCM:                                           
C                  COMMON /ASYCM/ EASYAB, EASYBC, EASYAC                        
C   The potential energy in the AB valley, EASYAB, is equal to the potential    
C   energy of the H "infinitely" far from the IH diatomic, with the             
C   IH diatomic at its equilibrium configuration.  Similarly, the terms         
C   EASYBC and EASYAC represent the H2 and the IH asymptotic valleys,           
C   respectively.                                                               
C   The other potential parameters are passed through the common blocks         
C   SATOCM and LEPSCM.                                                          
C   All the information passed through the common blocks PT1CM, ASYCM,          
C   SATOCM, and LEPSCM are in hartree atomic units.                             
C                                                                               
C        This potential is written such that:                                   
C                       R(1) = R(I-H)                                           
C                       R(2) = R(H-H)                                           
C                       R(3) = R(H-I)                                           
C   The classical potential energy is set equal to zero for the I               
C   infinitely far from the equilibrium H2 diatomic.                            
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
         COMMON /SATOCM/ DE(3), RE(3), BETA(3), Z(3)                            
         COMMON /LEPSCM/ ZPO(3), OP3Z(3), ZP3(3), TZP3(3), TOP3Z(3),            
     *                   DO4Z(3), B(3)                                          
         PARAMETER (CKCAL = 627.5095D0)                                         
         PARAMETER (CANGS =   0.529177106D0)                                    
C                                                                               
C   Echo the potential parameters                                               
C                                                                               
      IF(NATOMS.GT.25) THEN                                                     
         WRITE(NFLAG(18),1111)                                                  
 1111    FORMAT(2X,'STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS')            
         STOP                                                                   
      END IF                                                                    
C                                                                               
         WRITE (NFLAG(18), 100) DE, RE, BETA, Z                                 
C                                                                               
100   FORMAT (/, 2X, T5, 'PREPOT has been called for IH2 ',                     
     *                   'potential - eLEPS functional form',                   
     *        //, 2X, T5, 'Potential energy surface parameters:',               
     *        /, 2X, T5, 'Bond', T47, 'I-H', T58, 'H-H', T69, 'H-I',            
     *        /, 2X, T5, 'Dissociation energies (kcal/mol):',                   
     *        T44, F10.5, T55, F10.5, T66, F10.5,                               
     *        /, 2X, T5, 'Equilibrium bond lengths (Angstroms):',               
     *        T44, F10.5, T55, F10.5, T66, F10.5,                               
     *        /, 2X, T5, 'Morse beta parameters (Angstroms**-1):',              
     *        T44, F10.5, T55, F10.5, T66, F10.5,                               
     *        /, 2X, T5, 'Sato parameters:',                                    
     *        T44, F10.5, T55, F10.5, T66, F10.5,/)                             
C                                                                               
      DO  10 I = 1,3                                                            
C                                                                               
C   Convert to atomic units                                                     
C                                                                               
             DE(I)  = DE(I)/CKCAL                                               
             RE(I)   = RE(I)/CANGS                                              
             BETA(I) = BETA(I)*CANGS                                            
C                                                                               
C   Compute useful constants                                                    
C                                                                               
             ZPO(I)   = 1.0D0 + Z(I)                                            
             OP3Z(I)  = 1.0D0 + 3.0D0*Z(I)                                      
             TOP3Z(I) = 2.0D0*OP3Z(I)                                           
             ZP3(I)   = Z(I) + 3.0D0                                            
             TZP3(I)  = 2.0D0*ZP3(I)                                            
             DO4Z(I)  = DE(I)/4.0D0/ZPO(I)                                      
             B(I)     = BETA(I)*DO4Z(I)*2.0D0                                   
10    CONTINUE                                                                  
C                                                                               
C    Set the values of the classical energy in the three asymptotic valleys     
C                                                                               
             EASYAB = DE(1)                                                     
             EASYBC = DE(2)                                                     
             EASYAC = DE(3)                                                     
C                                                                               
C                                                                               
      EZERO(1)=DE(2)                                                            
C                                                                               
       DO I=1,5                                                                 
          REF(I) = ' '                                                          
       END DO                                                                   
C                                                                               
       REF(1)='D. S. Perry, J. C. Polanyi, C. W. Wilson, Jr.,'                  
       REF(2)='Chem. Phys. 3, 317(1974) [PES Parameters]'                       
       REF(3)='P. J. Kuntz, E. M. Nemeth, J. C. Polanyi,'                       
       REF(4)='S. D. Rosner, C. E. Young,'                                      
       REF(5)='J. Chem. Phys. 44, 1168(1966) [Functional Form]'                 
C                                                                               
      INDEXES(1) = 53                                                           
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
C*****                                                                          
C                                                                               
         SUBROUTINE POT                                                         
C                                                                               
C   System:          ABC                                                        
C   Functional form: Extended LEPS (London-Erying-Polyani-Sato)                 
C   References:      P. J. Kuntz, E. M. Nemth, J. C. Polanyi, S. D. Rosner,     
C                    and C. E. Young                                            
C                    J. Chem. Phys. 44, 1168 (1966)                             
C                                                                               
C   The potential parameters must be passed through the common blocks           
C   PT1CM, ASYCM, PT2CM, SATOCM, and LEPSCM.                                    
C   All information passed through the common blocks PT1CM, ASYCM,              
C   SATOCM, and LEPSCM must be in Hartree atomic units.                         
C                                                                               
C        For the reaction: A + BC -> AB + C we write:                           
C                          R(1) = R(A-B)                                        
C                          R(2) = R(B-C)                                        
C                          R(3) = R(C-A)                                        
C                                                                               
C   NOTE: The potential energy at the reactant asymptote, that is at            
C         A infinitely far from the BC diatomic, BC diatomic at its             
C         equilibrium configuration, is set equal to zero.                      
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
         COMMON /SATOCM/ DE(3), RE(3), BETA(3), Z(3)                            
         COMMON /LEPSCM/ ZPO(3), OP3Z(3), ZP3(3), TZP3(3), TOP3Z(3),            
     *                   DO4Z(3), B(3)                                          
         DIMENSION X(3), COUL(3), EXCH(3)                                       
         PARAMETER (R2 = 1.41421356D0)                                          
C                                                                               
      CALL CARTOU                                                               
      CALL CARTTOR                                                              
C                                                                               
C   Initialize the variable used in the calculation of the energy.              
C                                                                               
         ENGYGS = 0.D0                                                          
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
C   Compute the energy.                                                         
C                                                                               
         DO 10 I = 1,3                                                          
               X(I)    = EXP(-BETA(I)*(R(I)-RE(I)))                             
               COUL(I) = DO4Z(I)*(ZP3(I)*X(I)-TOP3Z(I))*X(I)                    
               EXCH(I) = DO4Z(I)*(OP3Z(I)*X(I)-TZP3(I))*X(I)                    
               ENGYGS  = ENGYGS + COUL(I)                                       
10       CONTINUE                                                               
C                                                                               
         RAD = SQRT((EXCH(1)-EXCH(2))**2 + (EXCH(2)-EXCH(3))**2 +               
     *              (EXCH(3)-EXCH(1))**2)                                       
C                                                                               
         ENGYGS = ENGYGS - RAD/R2 + EZERO(1)                                    
C                                                                               
C   Compute the derivatives of the energy with respect to the internal          
C   coordinates.                                                                
C                                                                               
         IF (NDER .EQ. 1) THEN                                                  
             S = EXCH(1) + EXCH(2) + EXCH(3)                                    
             DO 20 I = 1,3                                                      
                   DEGSDR(I) = B(I)*X(I)*((3.0D0*EXCH(I)-S)/R2*                 
     *                       (OP3Z(I)*X(I)-ZP3(I))/RAD-                         
     *                       ZP3(I)*X(I)+OP3Z(I))                               
20           CONTINUE                                                           
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
         COMMON /SATOCM/ DE(3), RE(3), BETA(3), Z(3)                            
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
C   Initialize the potential parameters; the energy parameters                  
C   are in kcal/mol, and the lengths are in Angstroms.                          
C                                                                               
         DATA DE/ 73.78D0, 109.43D0, 73.78D0/                                   
         DATA RE/ 1.604D0, 0.742D0, 1.604D0/                                    
         DATA BETA/ 1.75D0, 1.942D0, 1.75D0/                                    
         DATA Z/ 0.09D0, 0.23D0, 0.09D0/                                        
C                                                                               
         END                                                                    
C                                                                               
C*****                                                                          
