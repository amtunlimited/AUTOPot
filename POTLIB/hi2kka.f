      subroutine prepot                                                         
c                                                                               
c System:          HI2                                                          
c Functional form: eLEPS (extended London-Eyring-Polanyi-Sato)                  
c Common name:     KK                                                           
c Reference:       J. A. Kaye and A. Kuppermann, Chem. Phys. Lett.              
c                  77, 573 (1981).                                              
c                                                                               
c PREPOT must be called once before any calls to POT are made.  It should       
c never be called more than once.                                               
c The potential parameters are included in DATA statements.                     
c Coordinates, potential energy, derivatives, and other information             
c is passed through three common blocks:                                        
c           COMMON /PT1CM/ R(3), ENGYGS, DEGSDR(3)                              
c           COMMON /ASYCM/ EASYAB, EASYBC, EASYAC                               
c The common block PT1CM contains the internal coordinates, the potential       
c energy, and the derivatives with respect to internal coordinates.  The        
c common block PT2CM contains the a several variables which control             
c various aspects of the calculation as outlined below.  The last common        
c block, ASYCM, contains information about the differences in the               
c dissociation energies which is useful for certain calculations.               
c                                                                               
c Internuclear distances should be expressed in bohr and energies, etc. are     
c returned in Hartree atomic units.                                             
c                                                                               
c The coordinates are defined as follows:                                       
c           R(1) = R(I-H)                                                       
c           R(2) = R(H-I')                                                      
c           R(3) = R(I-I')                                                      
c                                                                               
c The zero of energy for this potential occurs when the I atom is               
c "infinitely" far from the HI molecule, and the HI molecule is at              
c its equilibrium separation.                                                   
c                                                                               
c The array DEGSDR contains the analytic first derivatives of the potential     
c energy with respect to the internal coordinates if the variable NDER is       
c equal to 1 (see not below).  The derivatives are defined as follows:          
c           DEGSDR(1) = dV(R)/dR(1)                                             
c           DEGSDR(2) = dV(R)/dR(2)                                             
c           DEGSDR(3) = dV(R)/dR(3)                                             
c where V(R) = V { R(1), R(2), R(3) }.                                          
c                                                                               
c The variable NDER controls whether derivatives are calculated.                
c If NDER = 1 first derivatives are calculated, otherwise no derivatives        
c are calculated.  Note:  the default value of NDER is 1.                       
c                                                                               
c The array NFLAG is unused in this potential surface.                          
c                                                                               
c The common block ASYCM is defined as follows:                                 
c           EASYAB = DE(I-H)                                                    
c           EASYBC = DE(H-I')                                                   
c           EASYAC = DE(I-I')                                                   
c where DE is the dissociation energy in atomic units.                          
c                                                                               
C   Potential parameters' default settings                                      
C                  Variable            Default value                            
C                  NDER                1                                        
C                  NFLAG(18)           6                                        
C                                                                               
      implicit real*8 (a-h,o-z)                                                 
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
      PARAMETER (cevau = 27.2113961d0)                                          
C                                                                               
C                                                                               
      COMMON /PRECM/ de(3),re(3),b(3),s(3),hde(3),                              
     +               hpf(3),pfs(3),pft(3),pf                                    
C                                                                               
C      save de,re,b,s,cevau                                                     
c                                                                               
c Write header potential parameters to UNIT NFLAG(18).                          
c                                                                               
      IF(NATOMS.GT.25) THEN                                                     
         WRITE(NFLAG(18),1111)                                                  
 1111    FORMAT(2X,'STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS')            
         STOP                                                                   
      END IF                                                                    
C                                                                               
      write(NFLAG(18),*)                                                        
      write(NFLAG(18),*) 'PREPOT has been called for the HI2 KK '               
      write(NFLAG(18),*) 'potential energy surface.  Scalar version '           
      write(NFLAG(18),*) '(optimized) including analytic first '                
      WRITE(NFLAG(18),*) 'derivatives in internal coordinates.'                 
      write(NFLAG(18),*) 'Last modified 7 February 1997 (TCA).'                 
      write(NFLAG(18),1000) de(1),de(2),de(3),re(1),re(2),re(3),b(1),           
     +                      b(2),b(3),s(1),s(2),s(3),cevau                      
c                                                                               
c Calculate first derivatives (default).  If no first derivatives are           
c desired, nder should be set equal to zero after prepot is called.             
c                                                                               
c      nder=1                                                                   
c                                                                               
c Convert to atomic units and calculate some useful constants.                  
c                                                                               
      do 10 i=1,3                                                               
        de(i)=de(i)/cevau                                                       
        pf=0.5d0*de(i)*((1.d0-s(i))/(1.d0+s(i)))                                
        hde(i)=0.5d0*de(i)                                                      
        hpf(i)=0.5d0*pf                                                         
        pfs(i)=-de(i)*b(i)                                                      
        pft(i)=-pf*b(i)                                                         
10    continue                                                                  
                                                                                
      easyab=de(1)                                                              
      easybc=de(2)                                                              
      easyac=de(3)                                                              
                                                                                
 1000 format(/,10x,'H-I  ',7x,'I-I'' ',7x,'H-I'' ',                             
     &       /,1x,'De',5x,3(f8.4,4x),'eV',                                      
     &       /,1x,'Re',5x,3(f8.4,4x),'bohr',                                    
     &       /,1x,'B ',5x,3(f8.4,4x),'bohr**-1',                                
     &       /,1x,'S ',5x,3(f8.4,4x),                                           
     &      //,1x,'1.00 hartree = ',f11.7,' eV',/)                              
                                                                                
C                                                                               
      EZERO(1)=DE(2)                                                            
C                                                                               
       DO I=1,5                                                                 
          REF(I) = ' '                                                          
       END DO                                                                   
C                                                                               
       REF(1)='J. A. Kaye and A. Kuppermann,'                                   
       REF(2)='Chem. Phys. Lett. 77, 573(1981)'                                 
C                                                                               
      INDEXES(1) = 53                                                           
      INDEXES(2) = 1                                                            
      INDEXES(3) = 53                                                           
C                                                                               
C                                                                               
C                                                                               
      IRCTNT=2                                                                  
C                                                                               
      CALL POTINFO                                                              
C                                                                               
      CALL ANCVRT                                                               
C                                                                               
      return                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE POT                                                            
C                                                                               
c                                                                               
c The coordinates are defined as follows:                                       
c           R(1) = R(I-H)                                                       
c           R(2) = R(H-I')                                                      
c           R(3) = R(I-I')                                                      
c                                                                               
c The zero of energy for this potential occurs when the I atom is               
c "infinitely" far from the HI molecule, and the HI molecule is at              
c its equilibrium separation.                                                   
c                                                                               
C      entry pot                                                                
      implicit real*8 (a-h,o-z)                                                 
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
      PARAMETER (cevau = 27.2113961d0)                                          
C                                                                               
C                                                                               
      COMMON /PRECM/ de(3),re(3),b(3),s(3),hde(3),                              
     +               hpf(3),pfs(3),pft(3),pf                                    
C                                                                               
C      save de,re,b,s,cevau                                                     
C                                                                               
      CALL CARTOU                                                               
      CALL CARTTOR                                                              
c                                                                               
c compute coulomb and exchange terms                                            
c                                                                               
        y1=dexp(b(1)*(re(1)-r(1)))                                              
        hs1=hde(1)*(y1-2.d0)*y1                                                 
        ht1=hpf(1)*(y1+2.d0)*y1                                                 
        wq1=hs1+ht1                                                             
        wj1=hs1-ht1                                                             
        y2=dexp(b(2)*(re(2)-r(2)))                                              
        hs2=hde(2)*(y2-2.d0)*y2                                                 
        ht2=hpf(2)*(y2+2.d0)*y2                                                 
        wq2=hs2+ht2                                                             
        wj2=hs2-ht2                                                             
        y3=dexp(b(3)*(re(3)-r(3)))                                              
        hs3=hde(3)*(y3-2.d0)*y3                                                 
        ht3=hpf(3)*(y3+2.d0)*y3                                                 
        wq3=hs3+ht3                                                             
        wj3=hs3-ht3                                                             
c                                                                               
c compute potential energy                                                      
c                                                                               
        rad=dsqrt(wj1*(wj1-wj2)+wj2*(wj2-wj3)+wj3*(wj3-wj1))                    
        ENGYGS=EZERO(1)  +wq1+wq2+wq3-rad                                       
c                                                                               
c  compute first derivatives if necessary                                       
c                                                                               
        if (nder .eq. 1) then                                                   
          radi=1.d0/rad                                                         
          hradi=0.5d0*radi                                                      
          hds1dr1=pfs(1)*(y1-1.d0)*y1                                           
          hds2dr2=pfs(2)*(y2-1.d0)*y2                                           
          hds3dr3=pfs(3)*(y3-1.d0)*y3                                           
          hdt1dr1=pft(1)*(y1+1.d0)*y1                                           
          hdt2dr2=pft(2)*(y2+1.d0)*y2                                           
          hdt3dr3=pft(3)*(y3+1.d0)*y3                                           
          dq1dr1=hds1dr1+hdt1dr1                                                
          dq2dr2=hds2dr2+hdt2dr2                                                
          dq3dr3=hds3dr3+hdt3dr3                                                
          dj1dr1=hds1dr1-hdt1dr1                                                
          dj2dr2=hds2dr2-hdt2dr2                                                
          dj3dr3=hds3dr3-hdt3dr3                                                
          DEGSDR(1)=dq1dr1-hradi*dj1dr1*(2.d0*wj1-wj2-wj3)                      
          DEGSDR(2)=dq2dr2-hradi*dj2dr2*(2.d0*wj2-wj1-wj3)                      
          DEGSDR(3)=dq3dr3-hradi*dj3dr3*(2.d0*wj3-wj1-wj2)                      
        end if                                                                  
                                                                                
      CALL EUNITZERO                                                            
      IF(NDER.NE.0) THEN                                                        
         CALL RTOCART                                                           
         CALL DEDCOU                                                            
      ENDIF                                                                     
C                                                                               
      return                                                                    
      end                                                                       
                                                                                
      BLOCK DATA PTPACM                                                         
C                                                                               
      implicit real*8 (a-h,o-z)                                                 
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
      COMMON /PRECM/ de(3),re(3),b(3),s(3),hde(3),                              
     +               hpf(3),pfs(3),pft(3),pf                                    
c                                                                               
c Constants used by the code.                                                   
c de is the dissociation energy (eV),                                           
c re is the equilibrium bond length (bohr),                                     
c b is the Morse beta parameter (bohr**-1),                                     
c s is the Sato parameter (unitless) for the appropriate diatom                 
c   (1 = I-H, 2 = H-I', 3 = I-I').                                              
c The other constants may be found in the literature reference.                 
c cevau is the conversion factor between eV and hartrees.                       
c                                                                               
      data de     / 3.3303d0, 3.3303d0, 1.5567d0 /                              
      data re     / 3.0236d0, 3.0236d0, 5.0457d0 /                              
      data b      / 0.9260d0, 0.9260d0, 0.9843d0 /                              
      data s      / 0.2000d0, 0.2000d0, 0.1250d0 /                              
c                                                                               
c                                                                               
      DATA NASURF /1,35*0 /                                                     
      DATA NDER /0/                                                             
      data NFLAG /1,1,15*0,6,0,0/                                               
      DATA ANUZERO /0.0D0/                                                      
                                                                                
      DATA ICARTR,MSURF,MDER/3,0,1/                                             
      DATA NULBL /25*0/                                                         
                                                                                
      DATA NATOMS /3/                                                           
                                                                                
      END                                                                       
