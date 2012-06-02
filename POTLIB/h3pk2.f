                                                                                
      SUBROUTINE PREPOT                                                         
C                                                                               
C   System:          H3                                                         
C   Common name:     PK2                                                        
C   Reference:       R. N. Porter and M. Karplus                                
C                    J. Chem. Phys. 40, 1105 (1964).                            
C                                                                               
C   PREPOT must be called once before any calls to POT.                         
C   The potential parameters are included in the block data subprogram PTPACM.  
C   Coordinates, potential energy, and derivatives are passed                   
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
C        NDER  - order of the derivatives of the energy with respect to         
C                the coordinates to be computed                                 
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
         COMMON /PARMCM/ D1,D3,RE,ALPH,BETA,EPS,DEL,FKAP,GAM                    
C                                                                               
C   Conversion constants                                                        
C   CEV (eV) converts from Hartree to kcal/mol                                  
C                                                                               
      PARAMETER (CEV = 27.21161D0)                                              
      PARAMETER (HARTRE = 627.5095D0)                                           
C                                                                               
C      DIMENSION EX(3),S(3),SS(3),E1(3),E3(3),SUME(3),                          
C     +          DIFE(3),Q(3),FJ(3),DE1(3),DE3(3),DQ(3),                        
C     +          DJ(3),DS(3),RR(3),DELJ(3)                                      
C                                                                               
C   CONVERT TO ATOMIC UNITS                                                     
C                                                                               
      IF(NATOMS.GT.25) THEN                                                     
         WRITE(NFLAG(18),1111)                                                  
 1111    FORMAT(2X,'STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS')            
         STOP                                                                   
      END IF                                                                    
C                                                                               
      D1 = D1/CEV                                                               
      D3 = D3/CEV                                                               
      EPS = EPS/CEV                                                             
      DEL = DEL/CEV                                                             
C                                                                               
C   echo the potential parameters                                               
C                                                                               
      WRITE (NFLAG(18), 600) D1,D3,EPS,DEL,RE,ALPH,BETA,FKAP,GAM                
C                                                                               
600   FORMAT (/,2X,T5,'PREPOT has been called for the H3 ',                     
     *                'potential energy surface PK2',                           
     *       //,2X,T5,'Potential energy surface parameters:',                   
     *        /,2X,T15,'D1',  T25,F15.8,T50,'D3',  T60,F15.8,                   
     *        /,2X,T15,'EPS', T25,F15.8,T50,'DEL', T60,F15.8,                   
     *        /,2X,T15,'RE',  T25,F15.8,T50,'ALPH',T60,F15.8,                   
     *        /,2X,T15,'BETA',T25,F15.8,T50,'FKAP',T60,F15.8,                   
     *        /,2X,T15,'GAM', T25,F15.8)                                        
C                                                                               
C    Set the values of the classical energy in the three asymptotic valleys     
C                                                                               
             EASYAB = D1                                                        
             EASYBC = EASYAB                                                    
             EASYAC = EASYAB                                                    
C                                                                               
      EZERO(1)=D1                                                               
C                                                                               
       DO I=1,5                                                                 
          REF(I) = ' '                                                          
       END DO                                                                   
C   Reference:                                                                  
C                                                                               
C                                                                               
       REF(1)='R. N. Porter and M. Karplus'                                     
       REF(2)='J. Chem. Phys. 40, 1105(1964)'                                   
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
C                                                                               
      SUBROUTINE POT                                                            
C                                                                               
C   The potential energy in the AB valley, EASYAB, is equal to the potential    
C   energy of the H "infinitely" far from the H2 diatomic, with the             
C   H2 diatomic at its equilibrium configuration.  For this potential energy    
C   surface the potential energy in the three asymptotic valleys are equivalent.
C                                                                               
C   The potential energy is defined to be equal to zero when one H is           
C   "infinitely" far from the H2 diatomic, and the H2 diatomic bond             
C   length is equal to the H2 equilibrium bond length.                          
C                                                                               
C      ENTRY POT                                                                
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
         COMMON /PARMCM/ D1,D3,RE,ALPH,BETA,EPS,DEL,FKAP,GAM                    
C                                                                               
C   Conversion constants                                                        
C   CEV (eV) converts from Hartree to kcal/mol                                  
C                                                                               
      PARAMETER (CEV = 27.21161D0)                                              
      PARAMETER (HARTRE = 627.5095D0)                                           
C                                                                               
      DIMENSION EX(3),S(3),SS(3),E1(3),E3(3),SUME(3),                           
     +          DIFE(3),Q(3),FJ(3),DE1(3),DE3(3),DQ(3),                         
     +          DJ(3),DS(3),RR(3),DELJ(3)                                       
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
      DO 12 I = 1,3                                                             
      GAMR = GAM*R(I)                                                           
      EXPG = EXP(-GAMR)                                                         
      ZETA = 1.0D0 + FKAP*EXPG                                                  
      RHO = R(I)*ZETA                                                           
      RR(I) = 1.0D0/R(I)                                                        
      EX(I) = DEL*EXP(-R(I)-R(I))                                               
      EXPR = EXP(-RHO)                                                          
      S(I) = (1.0D0+RHO+(RHO*RHO)/3.0D0)*EXPR                                   
      SS(I) = S(I)*S(I)                                                         
      REL = RE-R(I)                                                             
      EXA = EXP(ALPH*REL)                                                       
      EXB = EXP(BETA*REL)                                                       
      E1(I) = D1*(1.0D0-EXA)**2-D1                                              
      E3(I) = D3*(1.0D0+EXB)**2-D3                                              
      SUME(I) = E1(I)+E3(I)                                                     
      DIFE(I) = E1(I)-E3(I)                                                     
      Q(I) = 0.5D0*(SUME(I)+SS(I)*DIFE(I))                                      
      FJ(I) = 0.5D0*(DIFE(I)+SS(I)*SUME(I))                                     
      DE1(I) = 2.0D0*D1*ALPH*EXA*(1.0D0-EXA)                                    
      DE3(I) = 2.0D0*D3*BETA*EXB*(-1.0D0-EXB)                                   
      DS(I) = -RHO*(1.0D0+RHO)*EXPR*(ZETA+GAMR*(1.0D0-ZETA))/3.0D0              
      DJ(I) = S(I)*SUME(I)*DS(I)+0.5D0*((1.0D0+SS(I))*DE1(I)-(1.0D0-SS(I        
     V))*                                                                       
     1   DE3(I))                                                                
      DELJ(I) = (1.0D0+RR(I))*EX(I)                                             
   12 DQ(I) = S(I)*DIFE(I)*DS(I)+0.5D0*((1.0D0+SS(I))*DE1(I)+(1.0D0-SS(I        
     V))*                                                                       
     1   DE3(I))                                                                
      DJ(1) = DJ(1) + 2.0D0*S(1)*DS(1)*(DELJ(2)+DELJ(3))                        
      DJ(2) = DJ(2) + 2.0D0*S(2)*DS(2)*(DELJ(1)+DELJ(3))                        
      DJ(3) = DJ(3) + 2.0D0*S(3)*DS(3)*(DELJ(1)+DELJ(2))                        
      FJ(1) = FJ(1) + SS(1)*(DELJ(2)+DELJ(3))                                   
      FJ(2) = FJ(2) + SS(2)*(DELJ(1)+DELJ(3))                                   
      FJ(3) = FJ(3) + SS(3)*(DELJ(1)+DELJ(2))                                   
      S12 = SS(1)-SS(2)                                                         
      S23 = SS(2)-SS(3)                                                         
      S13 = S12+S23                                                             
      FJ12 = FJ(1)-FJ(2)                                                        
      FJ23 = FJ(2)-FJ(3)                                                        
      FJ13 = FJ12+FJ23                                                          
      DJ12 = -SS(2)*EX(1)*(2.0D0+2.0D0*RR(1)+RR(1)*RR(1))                       
      DJ13 = -SS(3)*EX(1)*(2.0D0+2.0D0*RR(1)+RR(1)*RR(1))                       
      DJ21 = -SS(1)*EX(2)*(2.0D0+2.0D0*RR(2)+RR(2)*RR(2))                       
      DJ23 = -SS(3)*EX(2)*(2.0D0+2.0D0*RR(2)+RR(2)*RR(2))                       
      DJ31 = -SS(1)*EX(3)*(2.0D0+2.0D0*RR(3)+RR(3)*RR(3))                       
      DJ32 = -SS(2)*EX(3)*(2.0D0+2.0D0*RR(3)+RR(3)*RR(3))                       
      DJ112 = DJ(1)-DJ12                                                        
      DJ113 = DJ(1)-DJ13                                                        
      DJ212 = DJ21-DJ(2)                                                        
      DJ223 = DJ(2)-DJ23                                                        
      DJ313 = DJ31-DJ(3)                                                        
      DJ323 = DJ32-DJ(3)                                                        
      DJ123 = DJ12-DJ13                                                         
      DJ213 = DJ21-DJ23                                                         
      DJ312 = DJ31-DJ32                                                         
      S123 = S(1)*S(2)*S(3)                                                     
      QT = Q(1)+Q(2)+Q(3)-EPS*S123                                              
      CS123 = 1.0D0-S123                                                        
      A = CS123*CS123-0.5D0*(S12*S12+S23*S23+S13*S13)                           
      B = -QT*CS123+0.5D0*(FJ12*S12+FJ23*S23+FJ13*S13)                          
      C = QT*QT-0.5D0*(FJ12*FJ12+FJ23*FJ23+FJ13*FJ13)                           
      A1 = -(DS(1)+DS(1))*(CS123*S(2)*S(3)+S(1)*(S12+S13))                      
      A2 = -(DS(2)+DS(2))*(CS123*S(1)*S(3)+S(2)*(S23-S12))                      
      A3 = -(DS(3)+DS(3))*(CS123*S(1)*S(2)-S(3)*(S13+S23))                      
      F = CS123*EPS+QT                                                          
      B1 = -CS123*DQ(1)+DS(1)*(S(2)*S(3)*F+S(1)*(FJ12+FJ13))                    
     1       +0.5D0*(S12*DJ112+S23*DJ123+S13*DJ113)                             
      B2 = -CS123*DQ(2)+DS(2)*(S(1)*S(3)*F+S(2)*(FJ23-FJ12))                    
     1       +0.5D0*(S12*DJ212+S23*DJ223+S13*DJ213)                             
      B3 = -CS123*DQ(3)+DS(3)*(S(1)*S(2)*F-S(3)*(FJ13+FJ23))                    
     1       +0.5D0*(S12*DJ312+S23*DJ323+S13*DJ313)                             
      F = QT+QT                                                                 
      C1 = F*(DQ(1)-EPS*S(2)*S(3)*DS(1))-(FJ12*DJ112+FJ23*DJ123                 
     1     +FJ13*DJ113)                                                         
      C2 = F*(DQ(2)-EPS*S(1)*S(3)*DS(2))-(FJ12*DJ212+FJ23*DJ223                 
     1     +FJ13*DJ213)                                                         
      C3 = F*(DQ(3)-EPS*S(2)*S(1)*DS(3))-(FJ12*DJ312+FJ23*DJ323                 
     1     +FJ13*DJ313)                                                         
      BAC = B*B-A*C                                                             
      RBAC = SQRT(BAC)                                                          
      ENGYGS0 = -B - RBAC                                                       
      IF(ENGYGS0.NE.0.D0) ENGYGS0 = ENGYGS0/A                                   
      ENGYGS = ENGYGS0 + EZERO(1)                                               
C                                                                               
      CALL EUNITZERO                                                            
      IF(NDER.NE.0) THEN                                                        
         T = 2.0D0*RBAC                                                    02JUL
         DEGSDR(1) = ENGYGS0*ENGYGS0*A1+(ENGYGS0+ENGYGS0)*B1+C1                 
         IF(DEGSDR(1).NE.0.D0) DEGSDR(1) = DEGSDR(1)/T                          
         DEGSDR(2) = ENGYGS0*ENGYGS0*A2+(ENGYGS0+ENGYGS0)*B2+C2                 
         IF(DEGSDR(2).NE.0.D0) DEGSDR(2) = DEGSDR(2)/T                          
         DEGSDR(3) = ENGYGS0*ENGYGS0*A3+(ENGYGS0+ENGYGS0)*B3+C3                 
         IF(DEGSDR(3).NE.0.D0) DEGSDR(3) = DEGSDR(3)/T                          
         CALL RTOCART                                                           
         CALL DEDCOU                                                            
      ENDIF                                                                     
C                                                                               
600   FORMAT (/,2X,T5,'PREPOT has been called for the H3 ',                     
     *                'potential energy surface PK2',                           
     *       //,2X,T5,'Potential energy surface parameters:',                   
     *        /,2X,T15,'D1',  T25,F15.8,T50,'D3',  T60,F15.8,                   
     *        /,2X,T15,'EPS', T25,F15.8,T50,'DEL', T60,F15.8,                   
     *        /,2X,T15,'RE',  T25,F15.8,T50,'ALPH',T60,F15.8,                   
     *        /,2X,T15,'BETA',T25,F15.8,T50,'FKAP',T60,F15.8,                   
     *        /,2X,T15,'GAM', T25,F15.8)                                        
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
         COMMON /PARMCM/ D1,D3,RE,ALPH,BETA,EPS,DEL,FKAP,GAM                    
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
C   Initialize the parameters for the potential                                 
C                                                                               
      DATA D1,D3,RE,ALPH,BETA,EPS,DEL,FKAP,GAM                                  
     +     / 4.7466D0,1.9668D0,1.40083D0,1.04435D0,                             
     +       1.000122D0,-17.5D0,28.2D0,.60D0,0.65D0/                            
         END                                                                    
C                                                                               
C*****                                                                          
