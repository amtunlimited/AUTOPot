C                                                                               
      SUBROUTINE PREPOT                                                         
C                                                                               
C   System:          H3                                                         
C   Functional Form: Double many-body expansion                                 
C   Common Name:     DMBE                                                       
C   Reference:       A. J. C. Varandas, F. B. Brown, C. A. Mead,                
C                    D. G. Truhlar, and N. C. Blais                             
C                    J. Chem. Phys. 86, 6258 (1987).                            
C                                                                               
C   Calling Sequence:                                                           
C      PREPOT - initializes the potential's variables and                       
C               must be called once before any calls to POT                     
C      POT    - driver for the evaluation of the energy and the derivatives     
C               of the energy with respect to the coordinates for a given       
C               geometric configuration                                         
C                                                                               
C   Units:                                                                      
C      energies    - hartrees                                                   
C      coordinates - bohrs                                                      
C      derivatives - hartrees/bohr                                              
C                                                                               
C   Surfaces:                                                                   
C      ground electronic state                                                  
C                                                                               
C   Zero of energy:                                                             
C      The classical potential energy is set equal to zero for the H            
C      infinitely far from the H2 diatomic and R(H2) set equal to the           
C      H2 equilibrium diatomic value.                                           
C                                                                               
C   Parameters:                                                                 
C      Set in the BLOCK DATA subprogram PTPACM                                  
C                                                                               
C   Coordinates:                                                                
C      Internal, Definition: R(1) = R( first H -- second H)                     
C                            R(2) = R(second H -- third H )                     
C                            R(3) = R( third H -- first H )                     
C                                                                               
C   Common Blocks (used between the calling program and this potential):        
C        passes the coordinates, ground state electronic energy, and            
C        derivatives of the ground electronic state energy with respect         
C        to the coordinates.                                                    
C        passes the control flags where                                         
C              = 1, first excited electronic state energy and derivatives       
C              = 2, ground and first excited electronic state energies and      
C                   derivatives                                                 
C        NDER  - not used                                                       
C        NFLAG - Control flags                                                  
C                                                                               
C              derivatives for the excited electronic states                    
C                                                                               
C                                                                               
C              derivatives for the nonadiabatic coupling between ground and     
C              excited states.                                                  
C                                                                               
C   Default Parameter Values:                                                   
C      Variable      Default value                                              
C      NFLAG(18)        6                                                       
C                                                                               
C***********************************************************************        
C     Calculates double many-body expansion of the H3 potential.                
C        ENGYGS  - energy of surface 1 (Ground State).        [ARRAY]           
C        DEGSDR  - derivatives for surface 1 (Ground State).  [ARRAY]           
C        ENGYES  - energy of surface 2 (Excited State).       [ARRAY]           
C        DEESDR  - derivatives for surface 2 (Excited State). [ARRAY]           
C        ENGYIJ  - nonadiabatic coupling energy.              [ARRAY]           
C        CAPF   - nonadiabatic coupling.  [SUBROUTINE TO CALCULATE ENGYIJ]      
C***********************************************************************        
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
      COMMON /PT3CM/  EZERO(ISURF+1)                                            
C                                                                               
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),                            
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
C                                                                               
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
C************************************************************************       
C                                                                       *       
C  Common block for H3 potential parameters, set in BLOCK DATA PTPACM   *       
C                                                                       *       
C************************************************************************       
C                                                                               
      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,                                 
     +                AL0,AL1,AL2,AL3,AZ2,                                      
     +                BETA1,BETA2,BETA3,BET0,BET1,                              
     +                CD0,CD1,CHH(3),CK0(3),CK1(3),                             
     +                DAMPA(3),DAMPB(3),                                        
     +                HFD,HFA1,HFA2,HFA3,HFGAM,                                 
     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,             
     +                RHOL,SQRT3,XPAR(15)                                       
C                                                                               
C******************************************************                         
C                                                     *                         
C  Common block for coordinates (computed in H3COOR)  *                         
C                                                     *                         
C******************************************************                         
C                                                                               
      COMMON /H3COCM/ PER,PER2,R12,R22,R32,                                     
     +                QCOORD,DQ1,DQ2,DQ3,                                       
     +                RHO,DRHO1,DRHO2,DRHO3,                                    
     +                S,S2,DS1,DS2,DS3,                                         
     +                CPHI3,DCPHI1,DCPHI2,DCPHI3                                
C                                                                               
C**************************************************                             
C                                                 *                             
C  Common block for 2-body correlation energies   *                             
C  and damping factors (computed in H3COR2)       *                             
C                                                 *                             
C**************************************************                             
C                                                                               
      COMMON /H3CRCM/ CORRS(3),DCORRS(3),                                       
     +                CORRT(3),DCORRT(3),                                       
     +                DAMP(3,3),DDAMP(3,3)                                      
C                                                                               
C**************************************************                             
C                                                 *                             
C  units conversion for anomalous adiabatic       *                             
C  coupling pengyij                               *                             
C                                                 *                             
C**************************************************                             
C                                                                               
      COMMON /IJCONVERT/CNVRTQ,CNVRTS,CNVRTCHI,ICOUPLE                          
C                                                                               
C     DIMENSION DLEP(3), DLEP2(3), DV2(3), DV3(3), DVA(3), DC3(3)               
C     DIMENSION R(3)                                                            
C                                                                               
      IF(NATOMS.GT.25) THEN                                                     
         WRITE(NFLAG(18),1111)                                                  
 1111    FORMAT(2X,'STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS')            
         STOP                                                                   
      END IF                                                                    
C                                                                               
      WRITE (NFLAG(18), 600)                                                    
600   FORMAT(/,1X,'*****','Potential Energy Surface',1X,'*****',/,              
     *       /,1X,T5,'H3 DMBE potential energy surface')                        
      IF (NASURF(1,1) .GT. 0) WRITE (NFLAG(18), 610)                            
610   FORMAT(/,1X,T5,'Ground electronic state surface')                         
      IF (NASURF(2,2) .GT. 0) WRITE (NFLAG(18), 620)                            
620   FORMAT(/,1X,T5,'First excited electronic state surface')                  
C                                                                               
C NON STANDARD COUPLING BETWEEN SURFACES                                        
C                                                                               
      ICOUPLE = 0                                                               
      IF (NFLAG(4).NE.0) THEN                                                   
         ICOUPLE = 1                                                            
         WRITE (NFLAG(18), 630)                                                 
630      FORMAT(/,1X,T5,'Ground and first excited coupling surface')            
         WRITE (NFLAG(18),*)                                                    
     x      '      This coupling is not in standard diabatic form'              
         WRITE (NFLAG(18),*)                                                    
     x      '      It is composed of 3 surfaces:'                               
         WRITE (NFLAG(18),*)                                                    
     x      '         <psi|dpse/dQ>,<psi|dpse/ds>,<psi|dpse/dchi>'              
         WRITE (NFLAG(18),*)                                                    
     x      '      It is activated by NFLAG(4), not NASURF'                     
         WRITE (NFLAG(18),*)                                                    
     x      '      No derivatives of these surfaces are coded'                  
         WRITE (NFLAG(18),*)                                                    
     x      '      To avoid interface problems with utility programs,',         
     x      ' nasurf12,21 set to 0'                                             
         NASURF(1,2) = 0                                                        
         NASURF(2,1) = 0                                                        
         CNVRTD = 1.D0                                                          
         IF(NFLAG(1).EQ.2) CNVRTD = .529177249D0                                
         CNVRTE = 1.D0                                                          
         IF(NFLAG(2).EQ.2) CNVRTE = 1000.D0                                     
         IF(NFLAG(2).EQ.3) CNVRTE = 27.211395D0                                 
         IF(NFLAG(2).EQ.4) CNVRTE = 627.5095D0                                  
         IF(NFLAG(2).EQ.5) CNVRTE = 219474.6D0                                  
         IF(NFLAG(2).EQ.6) CNVRTE = 2625.500D0                                  
         CNVRTQ = CNVRTE/CNVRTD                                                 
         CNVRTS = CNVRTE/CNVRTD                                                 
         CNVRTCHI = CNVRTE                                                      
      ENDIF                                                                     
      WRITE (NFLAG(18), 640)                                                    
640   FORMAT(//,1X,'*****')                                                     
C                                                                               
      NEXP = 4                                                                  
      DO 10 ID = 1,3                                                            
         NEXP = NEXP + 2                                                        
         DAMPA(ID) = ALPH0/DBLE(NEXP)**ALPH1                                    
         DAMPB(ID) = BET0*EXP(-BET1*DBLE(NEXP))                                 
   10 CONTINUE                                                                  
C                                                                               
C                                                                               
       DO I=1,5                                                                 
          REF(I) = ' '                                                          
       END DO                                                                   
C                                                                               
       REF(1)='A. J. C. Varandas, F. B. Brown, C. A. Mead,'                     
       REF(2)='D. G. Truhlar, N. C. Blais,'                                     
       REF(3)='J. Chem. Phys. 86, 6258(1987)'                                   
C                                                                               
      INDEXES(1) = 1                                                            
      INDEXES(2) = 1                                                            
      INDEXES(3) = 1                                                            
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
C   Surfaces:                                                                   
C      ground electronic state                                                  
C                                                                               
C   Zero of energy:                                                             
C      The classical potential energy is set equal to zero for the H            
C      infinitely far from the H2 diatomic and R(H2) set equal to the           
C      H2 equilibrium diatomic value.                                           
C                                                                               
C   Coordinates:                                                                
C      Internal, Definition: R(1) = R( first H -- second H)                     
C                            R(2) = R(second H -- third H )                     
C                            R(3) = R( third H -- first H )                     
C                                                                               
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (NATOM = 25)                                                    
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
      PARAMETER (PI = 3.141592653589793D0)                                      
C                                                                               
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)                         
      COMMON /PT3CM/  EZERO(ISURF+1)                                            
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                        
      COMMON /PT5CM/  ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                        
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
C************************************************************************       
C                                                                       *       
C  Common block for H3 potential parameters, set in BLOCK DATA PTPACM   *       
C                                                                       *       
C************************************************************************       
C                                                                               
      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,                                 
     +                AL0,AL1,AL2,AL3,AZ2,                                      
     +                BETA1,BETA2,BETA3,BET0,BET1,                              
     +                CD0,CD1,CHH(3),CK0(3),CK1(3),                             
     +                DAMPA(3),DAMPB(3),                                        
     +                HFD,HFA1,HFA2,HFA3,HFGAM,                                 
     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,             
     +                RHOL,SQRT3,XPAR(15)                                       
C                                                                               
C******************************************************                         
C                                                     *                         
C  Common block for coordinates (computed in H3COOR)  *                         
C                                                     *                         
C******************************************************                         
C                                                                               
      COMMON /H3COCM/ PER,PER2,R12,R22,R32,                                     
     +                QCOORD,DQ1,DQ2,DQ3,                                       
     +                RHO,DRHO1,DRHO2,DRHO3,                                    
     +                S,S2,DS1,DS2,DS3,                                         
     +                CPHI3,DCPHI1,DCPHI2,DCPHI3                                
C                                                                               
C**************************************************                             
C                                                 *                             
C  Common block for 2-body correlation energies   *                             
C  and damping factors (computed in H3COR2)       *                             
C                                                 *                             
C**************************************************                             
C                                                                               
      COMMON /H3CRCM/ CORRS(3),DCORRS(3),                                       
     +                CORRT(3),DCORRT(3),                                       
     +                DAMP(3,3),DDAMP(3,3)                                      
C                                                                               
C**************************************************                             
C                                                 *                             
C  units conversion for anomalous adiabatic       *                             
C  coupling pengyij                               *                             
C                                                 *                             
C**************************************************                             
C                                                                               
      COMMON /IJCONVERT/CNVRTQ,CNVRTS,CNVRTCHI,ICOUPLE                          
C                                                                               
      DIMENSION DLEP(3),DLEP2(3),                                               
     +          DV2(3),DV3(3),DVA(3),DC3(3)                                     
C     DIMENSION R(3)                                                            
C                                                                               
      CALL CARTOU                                                               
      CALL CARTTOR                                                              
C                                                                               
C Set up symmetry coordinates and their derivatives                             
C                                                                               
      CALL H3COOR (R(1), R(2), R(3))                                            
C                                                                               
C Set up 2-body correlation energies, dispersion damping terms, and             
C    their derivatives                                                          
C                                                                               
      CALL H3COR2                                                               
C                                                                               
C Get Leps potentials and derivatives                                           
C                                                                               
      CALL H3LEPS (VLEPS, VLEPS2, DLEP, DLEP2)                                  
C                                                                               
C Get VA term and its derivatives                                               
C                                                                               
      CALL H3VA (R(1), R(2), R(3), VA, DVA)                                     
C                                                                               
C Get VII term and its derivatives                                              
C                                                                               
      CALL H3VII (R(1), R(2), R(3), VII, DV2)                                   
C                                                                               
C Get VIII term and its derivatives                                             
C                                                                               
      CALL H3VIII (VIII, DV3)                                                   
C                                                                               
C Get 3-body correlation and its derivative                                     
C                                                                               
      CALL H3COR3 (CE3, DC3)                                                    
C                                                                               
C 2-body correlation term                                                       
C                                                                               
      CE2 = CORRS(1) + CORRS(2) + CORRS(3)                                      
C                                                                               
C ground state                                                                  
C                                                                               
      IF (NASURF(1,1) .GT. 0) THEN                                              
          VC = VA + VII + VIII + CE2 + CE3                                      
          ENGYGS = VLEPS + VC + EZERO(1)                                        
          IF(NDER.EQ.1) THEN                                                    
             DO 20 I = 1,3                                                      
                T = DVA(I) + DV2(I) + DV3(I) + DCORRS(I) + DC3(I)               
                DEGSDR(I) =( DLEP(I) + T)                                       
20           CONTINUE                                                           
          ENDIF                                                                 
      ENDIF                                                                     
C                                                                               
C first electronic state                                                        
C                                                                               
      IF (NASURF(2,2) .GT. 0) THEN                                              
          INDES = 1                                                             
          VC = VA + VII + VIII + CE2 + CE3                                      
          ENGYES(INDES) = VLEPS2 + VC + EZERO(2)                                
          IF(NDER.EQ.1) THEN                                                    
             DO 30 I = 1,3                                                      
                T = DVA(I) + DV2(I) + DV3(I) + DCORRS(I) + DC3(I)               
                DEESDR(I,INDES) =( DLEP2(I) + T)                                
30           CONTINUE                                                           
          ENDIF                                                                 
      ENDIF                                                                     
C                                                                               
C nonadiabatic coupling in an anomalous form where pengyij is treated as:       
C -  a vector with <psi|dpsi/dQ>,<psi|dpsi/ds>,<psi|dpsi/dchi>                  
C -  already transformed to the user's units                                    
C    (thus no RTOCART, DEDCOU calls via setting is nasurf(i,j))                 
C -  having no coded derivatives (no matter what nder is set to )               
C                                                                               
      IF (ICOUPLE .GT. 0) THEN                                                  
         CALL CAPF(R(1), R(2), R(3))                                            
         PENGYIJ(1) = ENGYIJ(1) * CNVRTQ                                        
         PENGYIJ(2) = ENGYIJ(2) * CNVRTS                                        
         PENGYIJ(3) = ENGYIJ(3) * CNVRTCHI                                      
      ENDIF                                                                     
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
      BLOCK DATA PTPACM                                                         
C                                                                               
C***********************************************************************        
C     Data for the DMBE H3 surface.                                             
C***********************************************************************        
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
      COMMON /PT3CM/  EZERO(ISURF+1)                                            
C                                                                               
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),                            
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
C                                                                               
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
C************************************************************************       
C                                                                       *       
C  Common block for H3 potential parameters, set in BLOCK DATA PTPACM   *       
C                                                                       *       
C************************************************************************       
C                                                                               
      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,                                 
     +                AL0,AL1,AL2,AL3,AZ2,                                      
     +                BETA1,BETA2,BETA3,BET0,BET1,                              
     +                CD0,CD1,CHH(3),CK0(3),CK1(3),                             
     +                DAMPA(3),DAMPB(3),                                        
     +                HFD,HFA1,HFA2,HFA3,HFGAM,                                 
     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,             
     +                RHOL,SQRT3,XPAR(15)                                       
C                                                                               
C******************************************************                         
C                                                     *                         
C  Common block for coordinates (computed in H3COOR)  *                         
C                                                     *                         
C******************************************************                         
C                                                                               
      COMMON /H3COCM/ PER,PER2,R12,R22,R32,                                     
     +                QCOORD,DQ1,DQ2,DQ3,                                       
     +                RHO,DRHO1,DRHO2,DRHO3,                                    
     +                S,S2,DS1,DS2,DS3,                                         
     +                CPHI3,DCPHI1,DCPHI2,DCPHI3                                
C                                                                               
C**************************************************                             
C                                                 *                             
C  Common block for 2-body correlation energies   *                             
C  and damping factors (computed in H3COR2)       *                             
C                                                 *                             
C**************************************************                             
C                                                                               
      COMMON /H3CRCM/ CORRS(3),DCORRS(3),                                       
     +                CORRT(3),DCORRT(3),                                       
     +                DAMP(3,3),DDAMP(3,3)                                      
C                                                                               
      DATA EZERO(1) /0.174474112D0/                                             
      DATA EZERO(2) /0.174474112D0/                                             
                                                                                
      DATA NASURF /1,35*0/                                                      
      DATA NDER /0/                                                             
      DATA NFLAG /1,1,15*0,6,0,0/                                               
C                                                                               
      DATA ANUZERO /0.0D0/                                                      
      DATA ICARTR,MSURF,MDER/3,1,1/                                             
      DATA NULBL /25*0/                                                         
      DATA NATOMS /3/                                                           
C                                                                               
      DATA RHOL /2.4848D0/                                                      
      DATA ALPH0, ALPH1, BET0, BET1                                             
     +     /25.9528D0,1.1868D0,15.7381D0,0.09729D0/                             
      DATA ALPHA5, BETA1, BETA2, BETA3                                          
     +     /8.2433033D-03,0.53302897D0,0.39156612D-01,0.69996945D0/             
      DATA ALPH2 /4.735364D-1/                                                  
      DATA AL0, AL1, AL2, AL3                                                   
     +     /0.45024472D+01,-0.62467617D+01,0.40966542D+01,                      
     +      0.21813012D+01/                                                     
      DATA AZ2 /4.453649D-4/                                                    
      DATA CD0, CD1 /6.333404D-03,-1.726839D-03/                                
      DATA CHH /  6.499027D0,      1.243991D+02,  3285.8D0/                     
      DATA CK0 / -1.372843D-01,   -1.638459D-01, -1.973814D-01/                 
      DATA CK1 /  1.011204D0,      9.988099D-01,  9.399411D-01/                 
      DATA HFD, HFA1, HFA2, HFA3, HFGAM                                         
     +     /0.15796326D0,2.1977034D0,1.2932502D0,0.64375666D0,                  
     +      2.1835071D0/                                                        
      DATA H2RM, H2R0, H2RMT /1.401D0, 6.928203D0, 7.82D0/                      
      DATA H2TA, H2TB1, H2TB2, H2TB3, H2TB4                                     
     +     /0.448467D0,-0.056687D0,0.246831D0,-0.018419D0,                      
     +      0.000598D0/                                                         
C      DATA PI / 3.14159265358979D0/                                            
      DATA SQRT3 /1.73205080756887D0/                                           
      DATA XPAR /-0.9286579D-02,  0.2811592D-03, -0.4665659D-05,                
     +            0.2069800D-07,  0.2903613D+02, -0.2934824D+01,                
     +            0.7181886D0,   -0.3753218D0,   -0.1114538D+01,                
     +            0.2134221D+01, -0.4164343D0,    0.2022584D0,                  
     +           -0.4662687D-01, -0.4818623D+02,  0.2988468D0/                  
      END                                                                       
C                                                                               
      SUBROUTINE H3COOR (R1, R2, R3)                                            
C                                                                               
C**********************************************************************         
C   Calculates D3H symmetry coordinates and derivatives                         
C**********************************************************************         
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
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)                         
      COMMON /PT3CM/  EZERO(ISURF+1)                                            
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                        
      COMMON /PT5CM/  ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)                       
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
C************************************************************************       
C                                                                       *       
C  Common block for H3 potential parameters, set in BLOCK DATA PTPACM   *       
C                                                                       *       
C************************************************************************       
C                                                                               
      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,                                 
     +                AL0,AL1,AL2,AL3,AZ2,                                      
     +                BETA1,BETA2,BETA3,BET0,BET1,                              
     +                CD0,CD1,CHH(3),CK0(3),CK1(3),                             
     +                DAMPA(3),DAMPB(3),                                        
     +                HFD,HFA1,HFA2,HFA3,HFGAM,                                 
     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,             
     +                RHOL,SQRT3,XPAR(15)                                       
C                                                                               
C******************************************************                         
C                                                     *                         
C  Common block for coordinates (computed in H3COOR)  *                         
C                                                     *                         
C******************************************************                         
C                                                                               
      COMMON /H3COCM/ PER,PER2,R12,R22,R32,                                     
     +                QCOORD,DQ1,DQ2,DQ3,                                       
     +                RHO,DRHO1,DRHO2,DRHO3,                                    
     +                S,S2,DS1,DS2,DS3,                                         
     +                CPHI3,DCPHI1,DCPHI2,DCPHI3                                
C                                                                               
C**************************************************                             
C                                                 *                             
C  Common block for 2-body correlation energies   *                             
C  and damping factors (computed in H3COR2)       *                             
C                                                 *                             
C**************************************************                             
C                                                                               
      COMMON /H3CRCM/ CORRS(3),DCORRS(3),                                       
     +                CORRT(3),DCORRT(3),                                       
     +                DAMP(3,3),DDAMP(3,3)                                      
C                                                                               
      COMMON/COUPCM/PHI,DPHI(3),SPHI3,DSPHI3(3)                                 
C                                                                               
      PER = R1+R2+R3                                                            
      PER2 = PER*PER                                                            
C                                                                               
C   QCOORD and its derivatives                                                  
C                                                                               
      R12 = R1*R1                                                               
      R22 = R2*R2                                                               
      R32 = R3*R3                                                               
      QCOORD = R12+R22+R32                                                      
      DQ1 = 2.0D0*R1                                                            
      DQ2 = 2.0D0*R2                                                            
      DQ3 = 2.0D0*R3                                                            
C                                                                               
C   RHO and its derivatives                                                     
C                                                                               
      RHO = SQRT(QCOORD/3.0D0)                                                  
      T = 1.0D0/(6.0D0*RHO)                                                     
      DRHO1 = T*DQ1                                                             
      DRHO2 = T*DQ2                                                             
      DRHO3 = T*DQ3                                                             
C                                                                               
C   S, CPHI3 (cos(phi3)), and their derivatives                                 
C                                                                               
      GAMMA = 2.0D0*R12-R22-R32                                                 
      GAM2 = GAMMA*GAMMA                                                        
      DGM1 = 2.0D0*DQ1                                                          
      DGM2 = -DQ2                                                               
      DGM3 = -DQ3                                                               
      BETA = SQRT3*(R22-R32)                                                    
      BET2 = BETA*BETA                                                          
      DBT1 = 0.0D0                                                              
      DBT2 = SQRT3*DQ2                                                          
      DBT3 = -SQRT3*DQ3                                                         
      T12 = BET2 + GAM2                                                         
      T1 = SQRT(T12)                                                            
      S = T1/QCOORD                                                             
      S2 = S*S                                                                  
      IF (S .EQ. 0.0D0) THEN                                                    
         DS1 = 0.0D0                                                            
         DS2 = 0.0D0                                                            
         DS3 = 0.0D0                                                            
C                                                                               
C   For S=0, CPHI3 and its derivative should not be                             
C   used anywhere but set to zero anyway.                                       
C                                                                               
         CPHI3 = 0.0D0                                                          
         DCPHI1 = 0.0D0                                                         
         DCPHI2 = 0.0D0                                                         
         DCPHI3 = 0.0D0                                                         
      ELSE                                                                      
         DS1 = S*((BETA*DBT1+GAMMA*DGM1)/T12 - DQ1/QCOORD)                      
         DS2 = S*((BETA*DBT2+GAMMA*DGM2)/T12 - DQ2/QCOORD)                      
         DS3 = S*((BETA*DBT3+GAMMA*DGM3)/T12 - DQ3/QCOORD)                      
         T2 = 1.0D0/(T1*T12)                                                    
         CPHI3 = GAMMA*(3.0D0*BET2 - GAM2)*T2                                   
         T3 = 3.0D0*BETA*(3.0D0*GAM2 - BET2)*T2/T12                             
         DCPHI1 = T3*(GAMMA*DBT1 - BETA*DGM1)                                   
         DCPHI2 = T3*(GAMMA*DBT2 - BETA*DGM2)                                   
         DCPHI3 = T3*(GAMMA*DBT3 - BETA*DGM3)                                   
      END IF                                                                    
      IF (NASURF(2,2).EQ.0.AND.NFLAG(4).EQ.0) GO TO 10                          
      ARG=GAMMA/T1                                                              
      IF(ARG.LE.-1.D0)ARG=-1.D0                                                 
      IF(ARG.GE.1.D0)ARG=1.D0                                                   
      THETA=ACOS(ARG)                                                           
      PHI=PI-THETA                                                              
      PHI3=3.0D0*PHI                                                            
      SPHI3=SIN(PHI3)                                                           
      RTAN=CPHI3/SPHI3                                                          
      DSPHI3(1)=-RTAN*DCPHI1                                                    
      DSPHI3(2)=-RTAN*DCPHI2                                                    
      DSPHI3(3)=-RTAN*DCPHI3                                                    
      DPHI(1)=-DCPHI1/(SPHI3*3.D0)                                              
      DPHI(2)=-DCPHI2/(SPHI3*3.D0)                                              
      DPHI(3)=-DCPHI3/(SPHI3*3.D0)                                              
10    CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE H3COR2                                                         
C                                                                               
C***********************************************************************        
C     Calculates 2-body correlation energies for singlet and triplet            
C     states of H2, the damping factors for the dispersion terms, and           
C     their 1st derivatives                                                     
C***********************************************************************        
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
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)                         
      COMMON /PT3CM/  EZERO(ISURF+1)                                            
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                        
      COMMON /PT5CM/  ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)                       
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
C      DIMENSION RR(3)                                                          
C                                                                               
C************************************************************************       
C                                                                       *       
C  Common block for H3 potential parameters, set in BLOCK DATA PTPACM   *       
C                                                                       *       
C************************************************************************       
C                                                                               
      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,                                 
     +                AL0,AL1,AL2,AL3,AZ2,                                      
     +                BETA1,BETA2,BETA3,BET0,BET1,                              
     +                CD0,CD1,CHH(3),CK0(3),CK1(3),                             
     +                DAMPA(3),DAMPB(3),                                        
     +                HFD,HFA1,HFA2,HFA3,HFGAM,                                 
     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,             
     +                RHOL,SQRT3,XPAR(15)                                       
C                                                                               
C******************************************************                         
C                                                     *                         
C  Common block for coordinates (computed in H3COR2)  *                         
C                                                     *                         
C******************************************************                         
C                                                                               
      COMMON /H3COCM/ PER,PER2,R12,R22,R32,                                     
     +                QCOORD,DQ1,DQ2,DQ3,                                       
     +                RHO,DRHO1,DRHO2,DRHO3,                                    
     +                S,S2,DS1,DS2,DS3,                                         
     +                CPHI3,DCPHI1,DCPHI2,DCPHI3                                
C                                                                               
C**************************************************                             
C                                                 *                             
C  Common block for 2-body correlation energies   *                             
C  and damping factors (computed in H3COR2)       *                             
C                                                 *                             
C**************************************************                             
C                                                                               
      COMMON /H3CRCM/ CORRS(3),DCORRS(3),                                       
     +                CORRT(3),DCORRT(3),                                       
     +                DAMP(3,3),DDAMP(3,3)                                      
C                                                                               
C                                                                               
C   Loop over three coordinates                                                 
C                                                                               
      DO 10 J = 1,3                                                             
         RSAVE = R(J)                                                           
         CORRS(J) = 0.0D0                                                       
         DCORRS(J) = 0.0D0                                                      
         CORRT(J) = 0.0D0                                                       
         DCORRT(J) = 0.0D0                                                      
         T1 = 1.0D0/RSAVE                                                       
         T = T1*T1                                                              
         T2 = T*T                                                               
         T1 = T2*T1                                                             
         NEXP = 4                                                               
C                                                                               
C      Loop over terms in dispersion expansion                                  
C                                                                               
         DO 1 ID = 1,3                                                          
            NEXP = NEXP + 2                                                     
            FNEXP = DBLE(NEXP)                                                  
            T2 = T2*T                                                           
            T1 = T1*T                                                           
C                                                                               
C     singlet                                                                   
C                                                                               
C        Calculate damping function for the nth dispersion coefficient          
C        and its 1st derivative                                                 
C                                                                               
            DENOM = H2RM + 2.5D0*H2R0                                           
            X = 2.0D0*RSAVE/DENOM                                               
            TT = DAMPB(ID)*X                                                    
            TT1 = DAMPA(ID) + TT                                                
            POL = X*TT1                                                         
            DPOL = TT1 + TT                                                     
            TT = EXP(-POL)                                                      
            TT1 = 1.0D0 - TT                                                    
            TT2 = TT1**(NEXP-1)                                                 
            D = TT1*TT2                                                         
            DD = 2.0D0*FNEXP*DPOL*TT*TT2/DENOM                                  
C                                                                               
C     store damping factors (including CHH coeff and 1/R**NEXP)                 
C     FOR later use in computing the 3-body correlation energy.                 
C                                                                               
            DAMP(ID,J) = CHH(ID)*D*T2                                           
            DDAMP(ID,J) = CHH(ID)*(DD*T2 - FNEXP*D*T1)                          
            CORRS(J) = CORRS(J) - DAMP(ID,J)                                    
            DCORRS(J) = DCORRS(J) - DDAMP(ID,J)                                 
C                                                                               
C     triplet                                                                   
C        Calculate damping function for the nth dispersion                      
C        coefficient and its 1st derivative                                     
C                                                                               
            DENOM = H2RMT + 2.5D0*H2R0                                          
            X = 2.0D0*RSAVE/DENOM                                               
            TT = DAMPB(ID)*X                                                    
            TT1 = DAMPA(ID) + TT                                                
            POL = X*TT1                                                         
            DPOL = TT1 + TT                                                     
            TT = EXP(-POL)                                                      
            TT1 = 1.0D0 - TT                                                    
            TT2 = TT1**(NEXP-1)                                                 
            D = TT1*TT2                                                         
            DD = 2.0D0*FNEXP*DPOL*TT*TT2/DENOM                                  
            CORRT(J) = CORRT(J) - CHH(ID)*D*T2                                  
            DCORRT(J) = DCORRT(J) - CHH(ID)*(DD*T2 - FNEXP*D*T1)                
    1    CONTINUE                                                               
   10 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE H3COR3 (CE3, DC3)                                              
C                                                                               
C***********************************************************************        
C     Calculates 3-body correlation energy and its 1st derivatives              
C***********************************************************************        
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
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)                         
      COMMON /PT3CM/  EZERO(ISURF+1)                                            
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                        
      COMMON /PT5CM/  ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)                       
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
      DIMENSION DC3(3), G(3), GD(3), H(3), HD(3)                                
C                                                                               
C************************************************************************       
C                                                                       *       
C  Common block for H3 potential parameters, set in BLOCK DATA PTPACM   *       
C                                                                       *       
C************************************************************************       
C                                                                               
      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,                                 
     +                AL0,AL1,AL2,AL3,AZ2,                                      
     +                BETA1,BETA2,BETA3,BET0,BET1,                              
     +                CD0,CD1,CHH(3),CK0(3),CK1(3),                             
     +                DAMPA(3),DAMPB(3),                                        
     +                HFD,HFA1,HFA2,HFA3,HFGAM,                                 
     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,             
     +                RHOL,SQRT3,XPAR(15)                                       
C                                                                               
C******************************************************                         
C                                                     *                         
C  Common block for coordinates (computed in H3COOR)  *                         
C                                                     *                         
C******************************************************                         
C                                                                               
      COMMON /H3COCM/ PER,PER2,R12,R22,R32,                                     
     +                QCOORD,DQ1,DQ2,DQ3,                                       
     +                RHO,DRHO1,DRHO2,DRHO3,                                    
     +                S,S2,DS1,DS2,DS3,                                         
     +                CPHI3,DCPHI1,DCPHI2,DCPHI3                                
C                                                                               
C**************************************************                             
C                                                 *                             
C  Common block for 2-body correlation energies   *                             
C  and damping factors (computed in H3COR2)       *                             
C                                                 *                             
C**************************************************                             
C                                                                               
      COMMON /H3CRCM/ CORRS(3),DCORRS(3),                                       
     +                CORRT(3),DCORRT(3),                                       
     +                DAMP(3,3),DDAMP(3,3)                                      
C                                                                               
      CE3 = 0.0D0                                                               
      DC3(1) = 0.0D0                                                            
      DC3(2) = 0.0D0                                                            
      DC3(3) = 0.0D0                                                            
C                                                                               
C   Loop over terms in dispersion expansion; dispersion damping                 
C   FACTORS are computed in H3COR2 and passed through COMMON /H3CRCM/.          
C                                                                               
      DO 30 ID = 1,3                                                            
C                                                                               
C   Loop over 3 coordinates                                                     
C                                                                               
         DO 10 I = 1,3                                                          
            RR = R(I)                                                           
C                                                                               
C     Calculate G function and its 1st derivative                               
C                                                                               
            T = CK0(ID)*EXP(-CK1(ID)*(RR-H2RM))                                 
            G(I) = 1.0D0 + T                                                    
            GD(I) = -CK1(ID)*T                                                  
C                                                                               
C     Calculate H function and its 1st derivative.                              
C                                                                               
            T = CK1(ID)*RR                                                      
            SGNT = 1.0D0                                                        
            IF (T .LT. 0.0D0) THEN                                              
               T = -T                                                           
               SGNT = -1.0D0                                                    
            END IF                                                              
            T = EXP(-T)                                                         
            T2 = T*T                                                            
            T1 = 1.0D0/(1.0D0+T2)                                               
            HYSEC = 2.0D0*T*T1                                                  
            HYTAN = SGNT*(1.0D0-T2)*T1                                          
            T1 = HYTAN**5                                                       
            H(I) = HYTAN*T1                                                     
            HD(I) = 6.0D0*CK1(ID)*T1*HYSEC*HYSEC                                
   10    CONTINUE                                                               
         DO 20 I = 1,3                                                          
            I2 = MOD(I,3) + 1                                                   
            I3 = MOD(I+1,3) + 1                                                 
            T = 1.0D0 - 0.5D0*(G(I2)*H(I3) + G(I3)*H(I2))                       
            CE3 = CE3 + T*DAMP(ID,I)                                            
            DC3(I) = DC3(I) + T*DDAMP(ID,I)                                     
            DC3(I2) = DC3(I2) -                                                 
     *         0.5D0*(GD(I2)*H(I3)+G(I3)*HD(I2))*DAMP(ID,I)                     
            DC3(I3) = DC3(I3) -                                                 
     *         0.5D0*(G(I2)*HD(I3)+GD(I3)*H(I2))*DAMP(ID,I)                     
   20    CONTINUE                                                               
   30 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE H3LEPS (VLEPS, VLEPS2, DLEP, DLEP2)                            
C                                                                               
C***********************************************************************        
C     Calculates 3-body extended-Hartree-Fock energy defined by a               
C     LEPS-type function.                                                       
C***********************************************************************        
C                                                                               
C   VLEPS  IS LEPS LOWER SURFACE.                                               
C   VLEPS2 IS LEPS UPPER SURFACE.                                               
C   DLEP(3) ARE LEPS LOWER SURFACE DERIVATIVES                                  
C   DLEP2(3) "   "   UPPER    "       "                                         
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
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)                         
      COMMON /PT3CM/  EZERO(ISURF+1)                                            
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                        
      COMMON /PT5CM/  ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)                       
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
C************************************************************************       
C                                                                       *       
C  Common block for H3 potential parameters, set in BLOCK DATA PTPACM   *       
C                                                                       *       
C************************************************************************       
C                                                                               
      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,                                 
     +                AL0,AL1,AL2,AL3,AZ2,                                      
     +                BETA1,BETA2,BETA3,BET0,BET1,                              
     +                CD0,CD1,CHH(3),CK0(3),CK1(3),                             
     +                DAMPA(3),DAMPB(3),                                        
     +                HFD,HFA1,HFA2,HFA3,HFGAM,                                 
     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,             
     +                RHOL,SQRT3,XPAR(15)                                       
C                                                                               
C******************************************************                         
C                                                     *                         
C  Common block for coordinates (computed in H3COOR)  *                         
C                                                     *                         
C******************************************************                         
C                                                                               
      COMMON /H3COCM/ PER,PER2,R12,R22,R32,                                     
     +                QCOORD,DQ1,DQ2,DQ3,                                       
     +                RHO,DRHO1,DRHO2,DRHO3,                                    
     +                S,S2,DS1,DS2,DS3,                                         
     +                CPHI3,DCPHI1,DCPHI2,DCPHI3                                
C                                                                               
C**************************************************                             
C                                                 *                             
C  Common block for 2-body correlation energies   *                             
C  and damping factors (computed in H3COR2)       *                             
C                                                 *                             
C**************************************************                             
C                                                                               
      COMMON /H3CRCM/ CORRS(3),DCORRS(3),                                       
     +                CORRT(3),DCORRT(3),                                       
     +                DAMP(3,3),DDAMP(3,3)                                      
C                                                                               
      DIMENSION DLEP(3), DLEP2(3)                                               
      DIMENSION DF(3),DQ(3),XJ(3),                                              
     +          DXJ(3,3)                                                        
C                                                                               
C  Calculate switching term for the Hartree-Fock component of the               
C  diatomic triplet state function.  The notation of Thompson et al.            
C  (J. Chem. Phys., 82,5597,1985; and references therein) is used               
C  throughout.                                                                  
C                                                                               
      T1 = -AZ2*QCOORD                                                          
      IF (S .EQ. 0.0D0) THEN                                                    
         F = EXP(T1*QCOORD)                                                     
         T1 = 2.0D0*T1*F                                                        
         DF(1) = T1*DQ1                                                         
         DF(2) = T1*DQ2                                                         
         DF(3) = T1*DQ3                                                         
      ELSE                                                                      
         T2 = 1.0D0 + S*S2*CPHI3                                                
         F = EXP(T1*QCOORD*T2)                                                  
         T1 = T1*F                                                              
         TDQ = 2.0D0*T2                                                         
         T2 = S2*QCOORD                                                         
         TDS = 3.0D0*CPHI3*T2                                                   
         TDC = S*T2                                                             
         DF(1) = T1*(TDQ*DQ1 + TDS*DS1 + TDC*DCPHI1)                            
         DF(2) = T1*(TDQ*DQ2 + TDS*DS2 + TDC*DCPHI2)                            
         DF(3) = T1*(TDQ*DQ3 + TDS*DS3 + TDC*DCPHI3)                            
      END IF                                                                    
      OMF = 1.0D0 - F                                                           
C                                                                               
      QSUM = 0.0D0                                                              
      DO 10 I = 1,3                                                             
         DQ(I) = 0.0D0                                                          
   10 CONTINUE                                                                  
C                                                                               
C  Loop over 3 coordinates                                                      
C                                                                               
      DO 50 I = 1,3                                                             
         RR = R(I)                                                              
         RRI = 1.0D0/RR                                                         
C                                                                               
C     Calculate extended-Hartree-Fock curve for ground singlet state of         
C     H2 [A.J.C. Varandas & J.D. Silva, J. Chem. Soc. Faraday II                
C     (submitted)] and 1st derivative of 2-body extended-Hartree-Fock           
C     curve for the ground-singlet state of H2.                                 
C                                                                               
         DR = RR - H2RM                                                         
         T = EXP(-HFGAM*DR)                                                     
         ES = -HFD*(1.0D0 + DR*(HFA1 + DR*(HFA2 + DR*HFA3)))*T                  
         DESDR = -HFGAM*ES - HFD*(HFA1 + DR*(2.0D0*HFA2 +                       
     *            3.0D0*DR*HFA3))*T                                             
C                                                                               
C     Compute the HFACE (Hartree-Fock-approximate correlation energy)           
C     potential for the lowest triplet state of H2 [A.J.C. Varandas &           
C     J. Brandao Mol. Phys.,45,1982,857] without the 2-body correlation         
C     energy.                                                                   
C                                                                               
         T = RR*(H2TB1 + RR*(H2TB2 + RR*(H2TB3 + RR*H2TB4)))                    
         AT = H2TA*EXP(-T)*RRI                                                  
         T = H2TB1 + RR*(2.0D0*H2TB2 + RR*(3.0D0*H2TB3 +                        
     *      RR*4.0D0*H2TB4))                                                    
         DAT = -AT*(T + RRI)                                                    
C                                                                               
C     Add in triplet and subtract singlet 2-body correlation terms              
C                                                                               
         AT = AT + CORRT(I) - CORRS(I)                                          
         DAT = DAT + DCORRT(I) - DCORRS(I)                                      
C                                                                               
C     Calculate effective diatomic triplet state curve.                         
C                                                                               
         T = EXP(-AL3*RR)                                                       
         WE= (AL0 + RR*(AL1 + RR*AL2))*T                                        
         DWE = (AL1 + 2.0D0*AL2*RR)*T - AL3*WE                                  
C                                                                               
C     Triplet energy                                                            
C                                                                               
         ET = F*WE + OMF*AT                                                     
C                                                                               
C     Contribution to coulomb term                                              
C                                                                               
         QSUM = QSUM + 0.5D0*(ES + ET)                                          
C                                                                               
C     Contribution to exchange term                                             
C                                                                               
         XJ(I) = 0.5D0*(ES - ET)                                                
C                                                                               
C     Derivatives of coulomb and exchange parts                                 
C                                                                               
         FTERM = WE - AT                                                        
         DO 20 J = 1,3                                                          
            T = 0.5D0*FTERM*DF(J)                                               
            DQ(J) = DQ(J) + T                                                   
            DXJ(I,J) = -T                                                       
   20    CONTINUE                                                               
         T = F*DWE + OMF*DAT                                                    
         DQ(I) = DQ(I) + 0.5D0*(DESDR + T)                                      
         DXJ(I,I) = DXJ(I,I) + 0.5D0*(DESDR - T)                                
   50 CONTINUE                                                                  
C                                                                               
      F1 = (XJ(1) - XJ(2))**2                                                   
      F2 = (XJ(2) - XJ(3))**2                                                   
      F3 = (XJ(3) - XJ(1))**2                                                   
      EXCH = SQRT(0.5D0*(F1 + F2 + F3))                                         
      VLEPS = QSUM - EXCH                                                       
      VLEPS2 = QSUM + EXCH                                                      
      XTERM1 = 2.0D0*XJ(1) - XJ(2) - XJ(3)                                      
      XTERM2 = 2.0D0*XJ(2) - XJ(1) - XJ(3)                                      
      XTERM3 = 2.0D0*XJ(3) - XJ(1) - XJ(2)                                      
      DO 80 I = 1,3                                                             
         TEXCH = EXCH                                                           
         IF(TEXCH.LE.0.0D0)TEXCH = 1.0D0                                        
         XTERM = 0.5D0/TEXCH*(XTERM1*DXJ(1,I)  + XTERM2*DXJ(2,I)                
     *      + XTERM3*DXJ(3,I))                                                  
         DLEP(I) = DQ(I) - XTERM                                                
         DLEP2(I) = DQ(I) + XTERM                                               
   80 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE H3VA (R1, R2, R3, VA, DVA)                                     
C                                                                               
C***********************************************************************        
C     Calculates Va correction energy and its 1st derivative.                   
C***********************************************************************        
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
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)                         
      COMMON /PT3CM/  EZERO(ISURF+1)                                            
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                        
      COMMON /PT5CM/  ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)                       
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
      DIMENSION DVA(3)                                                          
C                                                                               
C************************************************************************       
C                                                                       *       
C  Common block for H3 potential parameters, set in BLOCK DATA PTPACM   *       
C                                                                       *       
C************************************************************************       
C                                                                               
      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,                                 
     +                AL0,AL1,AL2,AL3,AZ2,                                      
     +                BETA1,BETA2,BETA3,BET0,BET1,                              
     +                CD0,CD1,CHH(3),CK0(3),CK1(3),                             
     +                DAMPA(3),DAMPB(3),                                        
     +                HFD,HFA1,HFA2,HFA3,HFGAM,                                 
     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,             
     +                RHOL,SQRT3,XPAR(15)                                       
C                                                                               
C******************************************************                         
C                                                     *                         
C  Common block for coordinates (computed in H3COOR)  *                         
C                                                     *                         
C******************************************************                         
C                                                                               
      COMMON /H3COCM/ PER,PER2,R12,R22,R32,                                     
     +                QCOORD,DQ1,DQ2,DQ3,                                       
     +                RHO,DRHO1,DRHO2,DRHO3,                                    
     +                S,S2,DS1,DS2,DS3,                                         
     +                CPHI3,DCPHI1,DCPHI2,DCPHI3                                
C                                                                               
C**************************************************                             
C                                                 *                             
C  Common block for 2-body correlation energies   *                             
C  and damping factors (computed in H3COR2)       *                             
C                                                 *                             
C**************************************************                             
C                                                                               
C      COMMON /H3CRCM/ CORRS(3),DCORRS(3),                                      
C     +                CORRT(3),DCORRT(3),                                      
C     +                DAMP(3,3),DDAMP(3,3)                                     
C                                                                               
      T = ALPHA5*PER2                                                           
      EXPVA = EXP(-T*PER)                                                       
      IF (EXPVA .EQ. 0.0D0) THEN                                                
         VA = 0.0D0                                                             
         DVA(1) = 0.0D0                                                         
         DVA(2) = 0.0D0                                                         
         DVA(3) = 0.0D0                                                         
      ELSE                                                                      
         VA = 0.0D0                                                             
         DV = 0.0D0                                                             
         T3 = (R1-R2)*(R2-R3)*(R3-R1)                                           
         T1 = T3*T3                                                             
         T2 = 1.0D0                                                             
         T4 = 2.0D0                                                             
         DO 1 J=1,4                                                             
            T2 = T2*T1                                                          
            VA = VA + XPAR(J)*T2                                                
            DV = DV + T4*XPAR(J)*T3                                             
            T3 = T3*T1                                                          
            T4 = T4 + 2.0D0                                                     
    1    CONTINUE                                                               
         VA = VA*EXPVA                                                          
         T1 = DV*EXPVA                                                          
         T2 = 3.0D0*T*VA                                                        
         DVA(1) = T1*(R2-R3)*(R3+R2-2.0D0*R1) - T2                              
         DVA(2) = T1*(R3-R1)*(R1+R3-2.0D0*R2) - T2                              
         DVA(3) = T1*(R1-R2)*(R1+R2-2.0D0*R3) - T2                              
      END IF                                                                    
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE H3VII (R1, R2, R3, E, DV2)                                     
C                                                                               
C***********************************************************************        
C     Calculates VII correction energy and its 1st derivative.                  
C***********************************************************************        
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
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)                         
      COMMON /PT3CM/  EZERO(ISURF+1)                                            
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                        
      COMMON /PT5CM/  ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)                       
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
      DIMENSION DV2(3)                                                          
C                                                                               
C************************************************************************       
C                                                                       *       
C  Common block for H3 potential parameters, set in BLOCK DATA PTPACM   *       
C                                                                       *       
C************************************************************************       
C                                                                               
      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,                                 
     +                AL0,AL1,AL2,AL3,AZ2,                                      
     +                BETA1,BETA2,BETA3,BET0,BET1,                              
     +                CD0,CD1,CHH(3),CK0(3),CK1(3),                             
     +                DAMPA(3),DAMPB(3),                                        
     +                HFD,HFA1,HFA2,HFA3,HFGAM,                                 
     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,             
     +                RHOL,SQRT3,XPAR(15)                                       
C                                                                               
C******************************************************                         
C                                                     *                         
C  Common block for coordinates (computed in H3COOR)  *                         
C                                                     *                         
C******************************************************                         
C                                                                               
      COMMON /H3COCM/ PER,PER2,R12,R22,R32,                                     
     +                QCOORD,DQ1,DQ2,DQ3,                                       
     +                RHO,DRHO1,DRHO2,DRHO3,                                    
     +                S,S2,DS1,DS2,DS3,                                         
     +                CPHI3,DCPHI1,DCPHI2,DCPHI3                                
C                                                                               
C**************************************************                             
C                                                 *                             
C  Common block for 2-body correlation energies   *                             
C  and damping factors (computed in H3COR2)       *                             
C                                                 *                             
C**************************************************                             
C                                                                               
C      COMMON /H3CRCM/ CORRS(3),DCORRS(3),                                      
C     +                CORRT(3),DCORRT(3),                                      
C     +                DAMP(3,3),DDAMP(3,3)                                     
C                                                                               
C   Compute B1 function.                                                        
C                                                                               
      R1I = 1.0D0/R1                                                            
      R2I = 1.0D0/R2                                                            
      R3I = 1.0D0/R3                                                            
      COS1 = 0.5D0*R2I*R3I*(R12-R22-R32)                                        
      COS2 = 0.5D0*R1I*R3I*(R22-R12-R32)                                        
      COS3 = 0.5D0*R1I*R2I*(R32-R12-R22)                                        
      WB = 1.0D0+COS1+COS2+COS3                                                 
      WB2 = WB*WB                                                               
C                                                                               
C   WB derivatives                                                              
C                                                                               
      WB1P = (R1*R3I - 1.0D0)*R2I - R3I - (COS2+COS3)*R1I                       
      WB2P = (R2*R3I - 1.0D0)*R1I - R3I - (COS1+COS3)*R2I                       
      WB3P = (R3*R2I - 1.0D0)*R1I - R2I - (COS1+COS2)*R3I                       
C                                                                               
C   EB1 term                                                                    
C                                                                               
      EXP1 = EXP(-BETA1*PER)                                                    
      EXP3 = EXP(-BETA3*PER)                                                    
      EB1T = (XPAR(5) + XPAR(6)*PER)*EXP1                                       
      EB3T = (XPAR(14) + XPAR(15)*PER2)*EXP3                                    
      EB1 = WB*(EB1T + EB3T)                                                    
C                                                                               
C   EB1 derivatives                                                             
C                                                                               
      EB1PR = WB*(-BETA1*EB1T - BETA3*EB3T + XPAR(6)*EXP1                       
     *   + 2.0D0*PER*XPAR(15)*EXP3)                                             
      EB1PWB = EB1T + EB3T                                                      
      EB1P1 = EB1PWB*WB1P + EB1PR                                               
      EB1P2 = EB1PWB*WB2P + EB1PR                                               
      EB1P3 = EB1PWB*WB3P + EB1PR                                               
C                                                                               
C   EB2 term                                                                    
C                                                                               
      T1 = BETA2*PER                                                            
      EXP2 = EXP(-T1*PER)                                                       
      EB2 = WB2*(XPAR(7) + WB*(XPAR(8) + WB*XPAR(9)))*EXP2                      
C                                                                               
C   EB2 derivatives                                                             
C                                                                               
      EB2PWB = WB*(2.0D0*XPAR(7) + WB*(3.0D0*XPAR(8) + WB*4.0D0*XPAR(9))        
     V)*EXP2                                                                    
      EB2PR = -2.0D0*T1*EB2                                                     
      EB2P1 = EB2PWB*WB1P + EB2PR                                               
      EB2P2 = EB2PWB*WB2P + EB2PR                                               
      EB2P3 = EB2PWB*WB3P + EB2PR                                               
C                                                                               
C   EB4 term                                                                    
C      EB4A                                                                     
C                                                                               
      T2 = XPAR(10)*EXP1                                                        
      T3 = WB*XPAR(11)*EXP2                                                     
      EB4A = WB*(T2 + T3)                                                       
C                                                                               
C      EB4A derivatives                                                         
C                                                                               
      EB4APW = T2 + 2.0D0*T3                                                    
      EB4APR = -WB*(BETA1*T2 + 2.0D0*T1*T3)                                     
      EB4AP1 = EB4APW*WB1P + EB4APR                                             
      EB4AP2 = EB4APW*WB2P + EB4APR                                             
      EB4AP3 = EB4APW*WB3P + EB4APR                                             
C                                                                               
C      EB4B                                                                     
C                                                                               
      T2 = XPAR(12)*EXP1                                                        
      T3 = XPAR(13)*EXP2                                                        
      EB4B = WB*(T2 + T3)                                                       
C                                                                               
C      EB4B derivatives                                                         
C                                                                               
      EB4BPW = T2 + T3                                                          
      EB4BPR = -WB*(BETA1*T2 + 2.0D0*T1*T3)                                     
      EB4BP1 = EB4BPW*WB1P + EB4BPR                                             
      EB4BP2 = EB4BPW*WB2P + EB4BPR                                             
      EB4BP3 = EB4BPW*WB3P + EB4BPR                                             
      DR12 = R1 - R2                                                            
      DR23 = R2 - R3                                                            
      DR31 = R3 - R1                                                            
      EQ = DR12*DR12 + DR23*DR23 + DR31*DR31                                    
      EQP1 = 2.0D0*(DR12 - DR31)                                                
      EQP2 = 2.0D0*(-DR12 + DR23)                                               
      EQP3 = 2.0D0*(-DR23 + DR31)                                               
      RI = R1I + R2I + R3I                                                      
      EB4 = EB4A*RI + EB4B*EQ                                                   
C                                                                               
C   EB4 derivatives                                                             
C                                                                               
      EB4P1 = EB4AP1*RI - EB4A/R12 + EB4BP1*EQ + EB4B*EQP1                      
      EB4P2 = EB4AP2*RI - EB4A/R22 + EB4BP2*EQ + EB4B*EQP2                      
      EB4P3 = EB4AP3*RI - EB4A/R32 + EB4BP3*EQ + EB4B*EQP3                      
      E = EB1 + EB2 + EB4                                                       
      DV2(1) = EB1P1 + EB2P1 + EB4P1                                            
      DV2(2) = EB1P2 + EB2P2 + EB4P2                                            
      DV2(3) = EB1P3 + EB2P3 + EB4P3                                            
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE H3VIII (VIII, DV3)                                             
C                                                                               
C**********************************************************************         
C     Calculates VIII correction energy and it 1st derivatives                  
C***********************************************************************        
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
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)                         
      COMMON /PT3CM/  EZERO(ISURF+1)                                            
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                        
      COMMON /PT5CM/  ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)                       
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
      DIMENSION DV3(3)                                                          
C                                                                               
C************************************************************************       
C                                                                       *       
C  Common block for H3 potential parameters, set in BLOCK DATA PTPACM   *       
C                                                                       *       
C************************************************************************       
C                                                                               
      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,                                 
     +                AL0,AL1,AL2,AL3,AZ2,                                      
     +                BETA1,BETA2,BETA3,BET0,BET1,                              
     +                CD0,CD1,CHH(3),CK0(3),CK1(3),                             
     +                DAMPA(3),DAMPB(3),                                        
     +                HFD,HFA1,HFA2,HFA3,HFGAM,                                 
     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,             
     +                RHOL,SQRT3,XPAR(15)                                       
C                                                                               
C******************************************************                         
C                                                     *                         
C  Common block for coordinates (computed in H3COOR)  *                         
C                                                     *                         
C******************************************************                         
C                                                                               
      COMMON /H3COCM/ PER,PER2,R12,R22,R32,                                     
     +                QCOORD,DQ1,DQ2,DQ3,                                       
     +                RHO,DRHO1,DRHO2,DRHO3,                                    
     +                S,S2,DS1,DS2,DS3,                                         
     +                CPHI3,DCPHI1,DCPHI2,DCPHI3                                
C                                                                               
C**************************************************                             
C                                                 *                             
C  Common block for 2-body correlation energies   *                             
C  and damping factors (computed in H3COR2)       *                             
C                                                 *                             
C**************************************************                             
C                                                                               
C      COMMON /H3CRCM/ CORRS(3),DCORRS(3),                                      
C     +                CORRT(3),DCORRT(3),                                      
C     +                DAMP(3,3),DDAMP(3,3)                                     
C                                                                               
      IF (S .EQ. 0.0D0) THEN                                                    
         VIII = 0.0D0                                                           
         DV3(1)=0.0D0                                                           
         DV3(2)=0.0D0                                                           
         DV3(3)=0.0D0                                                           
      ELSE                                                                      
         T1 = RHO - RHOL                                                        
         T5 = ALPH2*T1                                                          
         EXPV3 = EXP(-T5*T1)                                                    
         IF (EXPV3 .EQ. 0.0D0) THEN                                             
            VIII = 0.0D0                                                        
            DV3(1) = 0.0D0                                                      
            DV3(2) = 0.0D0                                                      
            DV3(3) = 0.0D0                                                      
         ELSE                                                                   
            T1 = S2*CPHI3                                                       
            T2 = 1.0D0 + S*T1                                                   
            T3 = S2*T2                                                          
            T4 = CD0 + CD1*RHO                                                  
            VIII = T3*T4*EXPV3                                                  
            TDS   = (2.0D0*S*T2 + 3.0D0*S2*T1)*T4                               
            TDCPH = S2*S*S2*T4                                                  
            TDRHO = T3*(CD1 - 2.0D0*T5*T4)                                      
            DV3(1) = (TDS*DS1 + TDCPH*DCPHI1 + TDRHO*DRHO1)*EXPV3               
            DV3(2) = (TDS*DS2 + TDCPH*DCPHI2 + TDRHO*DRHO2)*EXPV3               
            DV3(3) = (TDS*DS3 + TDCPH*DCPHI3 + TDRHO*DRHO3)*EXPV3               
         END IF                                                                 
      END IF                                                                    
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE CAPF(R1,R2,R3)                                                 
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
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)                         
      COMMON /PT3CM/  EZERO(ISURF+1)                                            
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                        
      COMMON /PT5CM/  ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)                       
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
C************************************************************************       
C                                                                       *       
C  Common block for H3 potential parameters, set in BLOCK DATA PTPACM   *       
C                                                                       *       
C************************************************************************       
C                                                                               
      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,                                 
     +                AL0,AL1,AL2,AL3,AZ2,                                      
     +                BETA1,BETA2,BETA3,BET0,BET1,                              
     +                CD0,CD1,CHH(3),CK0(3),CK1(3),                             
     +                DAMPA(3),DAMPB(3),                                        
     +                HFD,HFA1,HFA2,HFA3,HFGAM,                                 
     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,             
     +                RHOL,SQRT3,XPAR(15)                                       
C                                                                               
C******************************************************                         
C                                                     *                         
C  Common block for coordinates (computed in H3COOR)  *                         
C                                                     *                         
C******************************************************                         
C                                                                               
      COMMON /H3COCM/ PER,PER2,R12,R22,R32,                                     
     +                QCOORD,DQ1,DQ2,DQ3,                                       
     +                RHO,DRHO1,DRHO2,DRHO3,                                    
     +                S,S2,DS1,DS2,DS3,                                         
     +                CPHI3,DCPHI1,DCPHI2,DCPHI3                                
C                                                                               
C**************************************************                             
C                                                 *                             
C  Common block for 2-body correlation energies   *                             
C  and damping factors (computed in H3COR2)       *                             
C                                                 *                             
C**************************************************                             
C                                                                               
C      COMMON /H3CRCM/ CORRS(3),DCORRS(3),                                      
C     +                CORRT(3),DCORRT(3),                                      
C     +                DAMP(3,3),DDAMP(3,3)                                     
C                                                                               
      COMMON/COUPCM/PHI,DPHI(3),SPHI3,DSPHI3(3)                                 
      DIMENSION DX(3), DCAPK(4),DFZ(3),                                         
     +          DFU(3),DGZ(3),F(3,JSURF),RCAPF(3),                              
     +          DS(3), DCPHI(3)                                                 
      RCAPF(1)=R1                                                               
      RCAPF(2)=R2                                                               
      RCAPF(3)=R3                                                               
      CALL H3COOR(RCAPF(1),RCAPF(2),RCAPF(3))                                   
C                                                                               
      DS(1) = DS1                                                               
      DS(2) = DS2                                                               
      DS(3) = DS3                                                               
      DCPHI(1) = DCPHI1                                                         
      DCPHI(2) = DCPHI2                                                         
      DCPHI(3) = DCPHI3                                                         
C                                                                               
      X=SQRT(QCOORD/3.0D0)                                                      
      XI=1.0D0/X                                                                
      XI3=X/3.0D0                                                               
      X2=X*X                                                                    
      CALL KAPPA(X,CAPK,DCAPK)                                                  
      FZ=0.375D0*X*DCAPK(1)                                                     
      FU=0.03515625D0*(X*DCAPK(1)-X2*DCAPK(2)+X2*XI3*DCAPK(3))                  
      GZ=0.046875D0*(-X*DCAPK(1)+X2*DCAPK(2))                                   
      DENOM=FZ+GZ*CPHI3*S+FU*S2                                                 
      ARG=GZ*SPHI3*S/DENOM                                                      
C     FCT=PHI-ATAN(ARG)                                                         
      DO 10 I=1,3                                                               
         DX(I)=RCAPF(I)*XI                                                      
         DFZ(I)=0.125D0*DX(I)*(DCAPK(1)+X*DCAPK(2))                             
         DFU(I)=0.01171875D0*DX(I)*(DCAPK(1)-X*DCAPK(2)+X2*XI3*DCAPK(4))        
         DGZ(I)=0.015625D0*DX(I)*(-DCAPK(1)+X*DCAPK(2)+X2*DCAPK(3))             
         BRAK=S*DGZ(I)+GZ*DS(I)                                                 
         DARG=-ARG*(DFZ(I)+CPHI3*BRAK+S*(GZ*DCPHI(I)+S*DFU(I)+                  
     +         2.0D0*FU*DS(I)))                                                 
         DARG=(DARG+GZ*S*DSPHI3(I)+SPHI3*BRAK)/DENOM                            
         ENGYIJ(I)=DPHI(I)-1.0D0/(1.0D0+ARG**2)*DARG                            
   10 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE H2COR2 (RR)                                                    
C                                                                               
C*****************************************************************              
C     CALCULATES 1-4 DERIVATIVES OF 2-BODY CORRELATION TERMS FOR                
C     BOTH SINGLET AND TRIPLET STATES AND RETURNS THEIR DIFFERENCE              
C******************************************************************             
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
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)                         
      COMMON /PT3CM/  EZERO(ISURF+1)                                            
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                        
      COMMON /PT5CM/  ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)                       
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
C************************************************************************       
C                                                                       *       
C  Common block for H3 potential parameters, set in BLOCK DATA PTPACM   *       
C                                                                       *       
C************************************************************************       
C                                                                               
      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,                                 
     +                AL0,AL1,AL2,AL3,AZ2,                                      
     +                BETA1,BETA2,BETA3,BET0,BET1,                              
     +                CD0,CD1,CHH(3),CK0(3),CK1(3),                             
     +                DAMPA(3),DAMPB(3),                                        
     +                HFD,HFA1,HFA2,HFA3,HFGAM,                                 
     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,             
     +                RHOL,SQRT3,XPAR(15)                                       
C                                                                               
C******************************************************                         
C                                                     *                         
C  Common block for coordinates (computed in H3COOR)  *                         
C                                                     *                         
C******************************************************                         
C                                                                               
      COMMON /H3COCM/ PER,PER2,R12,R22,R32,                                     
     +                QCOORD,DQ1,DQ2,DQ3,                                       
     +                RHO,DRHO1,DRHO2,DRHO3,                                    
     +                S,S2,DS1,DS2,DS3,                                         
     +                CPHI3,DCPHI1,DCPHI2,DCPHI3                                
C                                                                               
C**************************************************                             
C                                                 *                             
C  Common block for 2-body correlation energies   *                             
C  and damping factors (computed in H3COR2)       *                             
C                                                 *                             
C**************************************************                             
C                                                                               
C      COMMON /H3CRCM/ CORRS(3),DCORRS(3),                                      
C     +                CORRT(3),DCORRT(3),                                      
C     +                DAMP(3,3),DDAMP(3,3)                                     
C                                                                               
      COMMON/DMPCM/DPOL,Y,TT,TT1,PARB                                           
      COMMON/COR2CM/CORRS,CORRT,DCOR2(4)                                        
         RSAVE = RR                                                             
         CORRS = 0.0D0                                                          
         DCORRS = 0.0D0                                                         
         CORRT = 0.0D0                                                          
         DCORRT = 0.0D0                                                         
      DO 11 J=1,4                                                               
   11 DCOR2(J)=0.0D0                                                            
         T1 = 1.0D0/RSAVE                                                       
         T = T1*T1                                                              
         T2 = T*T                                                               
         T1 = T2*T1                                                             
         NEXP = 4                                                               
C      Loop over terms in dispersion expansion                                  
         DO 1 ID = 1,3                                                          
            NEXP = NEXP + 2                                                     
            FNEXP = DBLE(NEXP)                                                  
            T2 = T2*T                                                           
            T1 = T1*T                                                           
C                                                                               
C******SINGLET                                                                  
C                                                                               
            DENOM = H2RM + 2.5D0*H2R0                                           
            Y=2.0D0/DENOM                                                       
            X=RSAVE*Y                                                           
            PARB=DAMPB(ID)                                                      
            TT = DAMPB(ID)*X                                                    
            TT1 = DAMPA(ID) + TT                                                
            POL = X*TT1                                                         
            DPOL = TT1 + TT                                                     
            TT = EXP(-POL)                                                      
            TT1 = 1.0D0 - TT                                                    
            TT2 = TT1**(NEXP-1)                                                 
            D = TT1*TT2                                                         
            DD = 2.0D0*FNEXP*DPOL*TT*TT2/DENOM                                  
            CORRS = CORRS - CHH(ID)*D*T2                                        
            DCORRS = DCORRS - CHH(ID)*(DD*T2-FNEXP*D*T1)                        
            CALL DAMP(RSAVE,NEXP,D2,D3,D4)                                      
            TEM2=D2                                                             
            TEM3=D3                                                             
            TEM4=D4                                                             
C                                                                               
C*****TRIPLET                                                                   
C                                                                               
            DENOM = H2RMT + 2.5D0*H2R0                                          
            Y=2.0D0/DENOM                                                       
            X=RSAVE*Y                                                           
            TT = DAMPB(ID)*X                                                    
            TT1 = DAMPA(ID) + TT                                                
            POL = X*TT1                                                         
            DPOL = TT1 + TT                                                     
            TT = EXP(-POL)                                                      
            TT1 = 1.0D0 - TT                                                    
            TT2 = TT1**(NEXP-1)                                                 
            D = TT1*TT2                                                         
            DD = 2.0D0*FNEXP*DPOL*TT*TT2/DENOM                                  
            CORRT = CORRT - CHH(ID)*D*T2                                        
            DCORRT = DCORRT - CHH(ID)*(DD*T2 - FNEXP*D*T1)                      
            CALL DAMP(RSAVE,NEXP,D2,D3,D4)                                      
            DCOR2(2)=DCOR2(2)+CHH(ID)*(D2-TEM2)                                 
            DCOR2(3)=DCOR2(3)+CHH(ID)*(D3-TEM3)                                 
            DCOR2(4)=DCOR2(4)+CHH(ID)*(D4-TEM4)                                 
    3       CONTINUE                                                            
    1    CONTINUE                                                               
      DCOR2(1)=DCORRT-DCORRS                                                    
   10 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE DAMP(RMP,N,D2,D3,D4)                                           
C                                                                               
C**********************************************************************         
C      CALCULATES 2ND,3RD AND 4TH DERIVATIVES OF 2-BODY CORRELATION             
C      ENERGIES FOR SINGLET AND TRIPLET STATES FOR USE IN NON-ADIABATIC         
C      COUPLINGS OF H3 SURFACES 1 AND 2. SENT TO SUBROUTINE H3COR2              
C**********************************************************************         
C                                                                               
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                        
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)                         
      COMMON /PT3CM/  EZERO(ISURF+1)                                            
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                        
      COMMON /PT5CM/  ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)                       
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
C************************************************************************       
C                                                                       *       
C  Common block for H3 potential parameters, set in BLOCK DATA PTPACM   *       
C                                                                       *       
C************************************************************************       
C                                                                               
C      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,                                
C     +                AL0,AL1,AL2,AL3,AZ2,                                     
C     +                BETA1,BETA2,BETA3,BET0,BET1,                             
C     +                CD0,CD1,CHH(3),CK0(3),CK1(3),                            
C     +                DAMPA(3),DAMPB(3),                                       
C     +                HFD,HFA1,HFA2,HFA3,HFGAM,                                
C     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,            
C     +                RHOL,SQRT3,XPAR(15)                                      
C                                                                               
C******************************************************                         
C                                                     *                         
C  Common block for coordinates (computed in H3COOR)  *                         
C                                                     *                         
C******************************************************                         
C                                                                               
C      COMMON /H3COCM/ PER,PER2,R12,R22,R32,                                    
C     +                QCOORD,DQ1,DQ2,DQ3,                                      
C     +                RHO,DRHO1,DRHO2,DRHO3,                                   
C     +                S,S2,DS1,DS2,DS3,                                        
C     +                CPHI3,DCPHI1,DCPHI2,DCPHI3                               
C                                                                               
C**************************************************                             
C                                                 *                             
C  Common block for 2-body correlation energies   *                             
C  and damping factors (computed in H3COR2)       *                             
C                                                 *                             
C**************************************************                             
C                                                                               
C      COMMON /H3CRCM/ CORRS(3),DCORRS(3),                                      
C     +                CORRT(3),DCORRT(3),                                      
C     +                DAMP(3,3),DDAMP(3,3)                                     
C                                                                               
      COMMON/DMPCM/DPOL,Y,TEXP,B,DAMPB                                          
      DPOL2=DPOL**2                                                             
      D1B=Y*DPOL*TEXP                                                           
      D2B=Y**2*(2.0D0*DAMPB-DPOL2)*TEXP                                         
      D3B=Y**3*DPOL*(DPOL2-6.0D0*DAMPB)*TEXP                                    
      D4B=Y**4*(-12.D0*DAMPB*DAMPB+DPOL2*(12.D0*DAMPB-DPOL2))*TEXP              
      T=DBLE(N)                                                                 
      TM=DBLE(N-1)                                                              
      TM2=DBLE(N-2)                                                             
      TM3=DBLE(N-3)                                                             
      TP=DBLE(N+1)                                                              
      TP2=DBLE(N+2)                                                             
      TP3=DBLE(N+3)                                                             
      TH=3.0D0                                                                  
      TW=2.0D0                                                                  
      FR=4.0D0                                                                  
      RR=1.0D0/RMP                                                              
      RRN=RR**N                                                                 
      BM4=B**(N-4)                                                              
      BM3=B*BM4                                                                 
      BM2=B*BM3                                                                 
      D2=T*BM2*(-TP*(B*RR)**2+TW*T*B*RR*D1B-TM*D1B**2-B*D2B)*RRN                
      D3=T*BM3*(TP*TP2*(B*RR)**3-TH*T*TP*(B*RR)**2*D1B+TH*TM*T*B*D1B**2*        
     1RR-TM2*TM*D1B**3+TH*T*RR*B**2*D2B-TH*TM*B*D1B*D2B-B**2*D3B)*RRN           
      D4=T*BM4*(-TP*TP2*TP3*(B*RR)**4+FR*T*TP*TP2*(B*RR)**3*D1B-TW*TH*TM        
     1*T*TP*(B*RR*D1B)**2+FR*TM2*TM*T*B*RR*D1B**3-TW*TH*T*TP*B*(B*RR)**         
     22*D2B+TH*FR*TM*T*B*B*RR*D1B*D2B-TW*TH*TM2*TM*B*D1B*D1B*D2B+FR*T*          
     3RR*B**3*D3B-TM3*TM2*TM*D1B**4-TH*TM*(B*D2B)**2-FR*TM*B*B*D1B              
     4*D3B-B**3*D4B)*RRN                                                        
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE KAPPA(QLC,CAPK,DCAPK)                                          
C                                                                               
C*****************************************************************              
C     CALCULATES FUNCTION CAPK AND ITS FIRST ROUR DERIVATIVES                   
C*****************************************************************              
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
      COMMON /PT1CM/  R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)                         
      COMMON /PT3CM/  EZERO(ISURF+1)                                            
      COMMON /PT4CM/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                        
      COMMON /PT5CM/  ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)                       
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
C************************************************************************       
C                                                                       *       
C  Common block for H3 potential parameters, set in BLOCK DATA PTPACM   *       
C                                                                       *       
C************************************************************************       
C                                                                               
      COMMON /H3DMCM/ ALPH2,ALPHA5,ALPH0,ALPH1,                                 
     +                AL0,AL1,AL2,AL3,AZ2,                                      
     +                BETA1,BETA2,BETA3,BET0,BET1,                              
     +                CD0,CD1,CHH(3),CK0(3),CK1(3),                             
     +                DAMPA(3),DAMPB(3),                                        
     +                HFD,HFA1,HFA2,HFA3,HFGAM,                                 
     +                H2RM,H2RMT,H2R0,H2TA,H2TB1,H2TB2,H2TB3,H2TB4,             
     +                RHOL,SQRT3,XPAR(15)                                       
C                                                                               
C******************************************************                         
C                                                     *                         
C  Common block for coordinates (computed in H3COOR)  *                         
C                                                     *                         
C******************************************************                         
C                                                                               
      COMMON /H3COCM/ PER,PER2,R12,R22,R32,                                     
     +                QCOORD,DQ1,DQ2,DQ3,                                       
     +                RHO,DRHO1,DRHO2,DRHO3,                                    
     +                S,S2,DS1,DS2,DS3,                                         
     +                CPHI3,DCPHI1,DCPHI2,DCPHI3                                
C                                                                               
C**************************************************                             
C                                                 *                             
C  Common block for 2-body correlation energies   *                             
C  and damping factors (computed in H3COR2)       *                             
C                                                 *                             
C**************************************************                             
C                                                                               
C      COMMON /H3CRCM/ CORRS(3),DCORRS(3),                                      
C     +                CORRT(3),DCORRT(3),                                      
C     +                DAMP(3,3),DDAMP(3,3)                                     
C                                                                               
C                                                                               
C                                                                               
      COMMON/COR2CM/CORRS,CORRT,DCOR2(4)                                        
      DIMENSION DCAPK(4)                                                        
      CALL H2COR2(QLC)                                                          
      T1 = -AZ2*QCOORD                                                          
      FZERO=EXP(T1*QCOORD)                                                      
      RR = QLC                                                                  
      RRI = 1.0D0/RR                                                            
C                                                                               
C*****CALCULATE 1-4 DERIVATIVES OF VEHF                                         
C                                                                               
      DR = RR - H2RM                                                            
      T = EXP(-HFGAM*DR)                                                        
      ES = -HFD*(1.0D0 + DR*(HFA1 + DR*(HFA2 + DR*HFA3)))*T                     
      DESDR = -HFGAM*ES - HFD*(HFA1 + DR*(2.0D0*HFA2 +                          
     *      3.0D0*DR*HFA3))*T                                                   
      P1=HFA1+2.0D0*HFA2*DR+3.0D0*HFA3*DR*DR                                    
      P3=6.0D0*HFA3                                                             
      P2=2.0D0*HFA2 +P3*DR                                                      
      D1VHF=DESDR                                                               
      D2VHF=-HFGAM*D1VHF-HFD*T*(P2-HFGAM*P1)                                    
      D3VHF=-HFGAM*D2VHF-HFD*T*(P3-2.D0*HFGAM*P2+HFGAM*HFGAM*P1)                
      D4VHF=-HFGAM*D3VHF+HFD*T*(3.D0*HFGAM*P3-3.D0*HFGAM*HFGAM*P2+              
     1 HFGAM*HFGAM*HFGAM*P1)                                                    
      DCAPK(1)=D1VHF                                                            
      DCAPK(2)=D2VHF                                                            
      DCAPK(3)=D3VHF                                                            
      DCAPK(4)=D4VHF                                                            
C                                                                               
C*****CALCULATE 1-4 DERIVATIVES OF UEHF                                         
C                                                                               
         T = RR*(H2TB1 + RR*(H2TB2 + RR*(H2TB3 + RR*H2TB4)))                    
         AT = H2TA*EXP(-T)*RRI                                                  
         T = H2TB1 + RR*(2.0D0*H2TB2 + RR*(3.0D0*H2TB3 +                        
     *      RR*4.0D0*H2TB4))+RRI                                                
         DAT = -AT*T                                                            
         Q1=-RRI**2+2.D0*H2TB2+6.D0*H2TB3*RR+12.D0*H2TB4*RR**2                  
         Q2=2.D0*RRI**3+6.D0*H2TB3+24.D0*H2TB4*RR                               
         Q3=-6.D0*RRI**4+24.D0*H2TB4                                            
         D1UEHF=DAT                                                             
         D2UEHF=-T*D1UEHF-AT*Q1                                                 
         D3UEHF=-T*D2UEHF-2.D0*Q1*D1UEHF-AT*Q2                                  
         D4UEHF=-T*D3UEHF-3.D0*Q1*D2UEHF-3.D0*Q2*D1UEHF-AT*Q3                   
C                                                                               
C     Add in triplet and subtract singlet 2-body correlation terms              
C                                                                               
         AT = AT + CORRT - CORRS                                                
C                                                                               
C*****INCLUDE CORRELATION PART OF CAPK,(USUBC-VSUBC)*(1-FZERO)                  
C                                                                               
      OMFZ=1.0D0-FZERO                                                          
      DCAPK(1)=DCAPK(1)-(D1UEHF+DCOR2(1))*OMFZ                                  
      DCAPK(2)=DCAPK(2)-(D2UEHF+DCOR2(2))*OMFZ                                  
      DCAPK(3)=DCAPK(3)-(D3UEHF+DCOR2(3))*OMFZ                                  
      DCAPK(4)=DCAPK(4)-(D4UEHF+DCOR2(4))*OMFZ                                  
C                                                                               
C*****CALCULATE 1-4 DERIVATIVES OF W-EFFECTIVE                                  
C                                                                               
         T = EXP(-AL3*RR)                                                       
         WE= (AL0 + RR*(AL1 + RR*AL2))*T                                        
         DWE = (AL1 + 2.0D0*AL2*RR)*T - AL3*WE                                  
      O2=2.0D0*AL2                                                              
      O1=AL1+O2*RR                                                              
      D1W=DWE                                                                   
      D2W=-AL3*D1W+T*(O2-AL3*O1)                                                
      D3W=-AL3*D2W+T*(AL3*O1-2.0D0*O2)*AL3                                      
      D4W=-AL3*D3W+T*(3.0D0*O2-AL3*O1)*AL3**2                                   
      DCAPK(1)=DCAPK(1)-FZERO*D1W                                               
      DCAPK(2)=DCAPK(2)-FZERO*D2W                                               
      DCAPK(3)=DCAPK(3)-FZERO*D3W                                               
      DCAPK(4)=DCAPK(4)-FZERO*D4W                                               
      CAPK=ES-FZERO*WE-OMFZ*AT                                                  
      RETURN                                                                    
      END                                                                       
C                                                                               
