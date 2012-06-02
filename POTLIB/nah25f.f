C                                                                               
      subroutine prepot                                                         
C                                                                               
C                                                                               
C   System:          NaH2                                                       
C   Common name:     5F                                                         
C   Number of electronic surfaces: 2                                            
C   Number of derivatives: 0                                                    
C                                                                               
C   Reference:       P. Halvick and D. G. Truhlar, J. Chem. Phys. 96,           
C                    2895 (1992); E 100, 4718 (1994)                            
C                                                                               
C   Calling Sequence:                                                           
C                                                                               
C      PREPOT - initializes the potential's variables and                       
C               must be called once before any calls to POT                     
C      POT    - driver for the evaluation of the energy                         
C                                                                               
C   Units:                                                                      
C                                                                               
C      energies    - hartrees                                                   
C      coordinates - bohr                                                       
C      derivatives - hartrees/bohr                                              
C                                                                               
C   Surfaces:                                                                   
C                                                                               
C      ground and first excited electronic states, diabatic representation      
C                                                                               
C   Zero of energy:                                                             
C                                                                               
C      The zero of energy occurs on the lower adiabatic surface when            
C      the Na atom is "infinitely" far from the H2 diatom and the               
C      H2 diatom is at its equilibrium geometry.  On this surface this occurs   
C      at R(H-H)=1.402721754031838 bohr with an energy of 1.842768357662727E-4  
C      hartrees (=5.014430 meV)                                                 
C                                                                               
      implicit double precision (a-h,o-z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
         COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                         
         COMMON /PT3CM/ EZERO(ISURF+1)                                          
         COMMON /ASYCM/ EASYBC, EASYAC, EASYBA                                  
         COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                      
         COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                      
C                                                                               
      COMMON/INFOCM/CARTNU(NATOM,3),INDEXES(NATOM),                             
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
         COMMON /PRE11CM/ dt2,ret2,betat2,                                      
     +                    y90c1,y90c2,y90c3,y90c4,                              
     +                     y0c1, y0c2, y0c3, y0c4,                              
     +                    ads1,ars1,abs1,                                       
     +                    cds1,crs1,cbs1,                                       
     +                    swa,swr0,swa2,swr20,                                  
     +                    c2a,c2b,c2c,deh2                                      
C                                                                               
         COMMON /PRE22CM/ c1s1,c2s1,c3s1,c4s1,c5s1,c6s1,                        
     +                    deh222,una2p,                                         
     +                    c2a22,c2b22,c2c22,cprime22,                           
     +                    y0d1,y0r1,y0b1,y90d1,y90r1,y90b1,                     
     +                    y0d2,y0r2,y0b2,y90d2,y90r2,y90b2,                     
     +                    y0swa,y0swr0,y90swa,y90swr0,                          
     +                    a0g1,b0g1,r0g1,a0g2,b0g2,r0g2,                        
     +                    a0g3,b0g3,r0g3,                                       
     +                    a90g1,b90g1,r90g1,a90g2,b90g2,r90g2,                  
     +                    a90g3,b90g3,r90g3                                     
C                                                                               
      IF(NATOMS.GT.25) THEN                                                     
         WRITE(NFLAG(18),1111)                                                  
 1111    FORMAT(2X,'STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS')            
         STOP                                                                   
      END IF                                                                    
C                                                                               
        call prepot11                                                           
        call prepot12                                                           
        call prepot22                                                           
C                                                                               
      EZERO(1)=DEH2                                                             
      EZERO(2)=DEH222                                                           
C                                                                               
       DO I=1,5                                                                 
          REF(I) = ' '                                                          
       END DO                                                                   
C                                                                               
       REF(1)='P. Halvick, D. G. Truhlar,'                                      
       REF(2)='J. Chem. Phys. 96, 2895(1992); E 100, 4718(1994)'                
C                                                                               
      INDEXES(1) = 11                                                           
      INDEXES(2) = 1                                                            
      INDEXES(3) = 1                                                            
C                                                                               
      IRCTNT=2                                                                  
C                                                                               
      CALL POTINFO                                                              
C                                                                               
      CALL ANCVRT                                                               
C                                                                               
      return                                                                    
C                                                                               
      end                                                                       
C                                                                               
      SUBROUTINE POT                                                            
C                                                                               
      implicit double precision (a-h,o-z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
         COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                         
         COMMON /PT3CM/ EZERO(ISURF+1)                                          
         COMMON /ASYCM/ EASYBC, EASYAC, EASYBA                                  
         COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                      
         COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                      
C                                                                               
      COMMON/INFOCM/CARTNU(NATOM,3),INDEXES(NATOM),                             
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
      CALL CARTOU                                                               
      CALL CARTTOR                                                              
C                                                                               
        if (nasurf(1,1).gt. 0) call pot11                                       
        if (nasurf(2,2).gt. 0) call pot22                                       
        if (nasurf(1,2)+nasurf(2,1).gt. 0) call pot12                           
C                                                                               
      CALL EUNITZERO                                                            
      IF(NDER.NE.0) THEN                                                        
         CALL RTOCART                                                           
         CALL DEDCOU                                                            
      ENDIF                                                                     
C                                                                               
      return                                                                    
C                                                                               
      end                                                                       
c                                                                               
c  U11 surface - lower diabatic surface                                         
c                                                                               
      subroutine prepot11                                                       
      implicit double precision (a-h,o-z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
         PARAMETER (CANGAU = 0.52917706D0)                                      
         PARAMETER (CHARKC = 627.5095D0)                                        
         PARAMETER (cevau=1.d0/27.211611d0)                                     
C                                                                               
         COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                         
         COMMON /PT3CM/ EZERO(ISURF+1)                                          
         COMMON /ASYCM/ EASYBC, EASYAC, EASYBA                                  
         COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                      
         COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                      
C                                                                               
      COMMON/INFOCM/CARTNU(NATOM,3),INDEXES(NATOM),                             
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
         COMMON /PRE11CM/ dt2,ret2,betat2,                                      
     +                    y90c1,y90c2,y90c3,y90c4,                              
     +                     y0c1, y0c2, y0c3, y0c4,                              
     +                    ads1,ars1,abs1,                                       
     +                    cds1,crs1,cbs1,                                       
     +                    swa,swr0,swa2,swr20,                                  
     +                    c2a,c2b,c2c,deh2                                      
C                                                                               
      rac2 = 1.d0/dsqrt(2.d0)                                                   
C                                                                               
      WRITE(NFLAG(18),100) dt2,ret2,betat2,y90c1,y90c2,y90c3,y90c4,             
     +             y0c1,y0c2,y0c3,y0c4,                                         
     +             ads1,ars1,abs1,cds1,crs1,cbs1,                               
     +             swa,swr0,swa2,swr20,c2a,c2b,c2c                              
C                                                                               
100   format(/,24('='),' NaH2 - U11 potential ',24('='),                        
     &  /'   dt2 =',1pe14.6,2x,'  ret2 =',e14.6,2x,'betat2 =',e14.6,            
     &  /' y90c1 =',e14.6,2x,' y90c2 =',e14.6,2x,' y90c3 =',e14.6,              
     &  /' y90c4 =',e14.6,2x,'  y0c1 =',e14.6,2x,'  y0c2 =',e14.6,              
     &  /'  y0c3 =',e14.6,2x,'  y0c4 =',e14.6,2x,'  ads1 =',e14.6,              
     &  /'  ars1 =',e14.6,2x,'  abs1 =',e14.6,2x,'  cds1 =',e14.6,              
     &  /'  crs1 =',e14.6,2x,'  cbs1 =',e14.6,2x,'   swa =',e14.6,              
     &  /'  swr0 =',e14.6,2x,'  swa2 =',e14.6,2x,' swr20 =',e14.6,              
     &  /'   c2a =',e14.6,2x,'   c2b =',e14.6,2x,'   c2c =',e14.6,              
     &  /70('='))                                                               
C                                                                               
      return                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE POT11                                                          
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
         PARAMETER (CANGAU = 0.52917706D0)                                      
         PARAMETER (CHARKC = 627.5095D0)                                        
         PARAMETER (cevau=1.d0/27.211611d0)                                     
C                                                                               
         COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                         
         COMMON /PT3CM/ EZERO(ISURF+1)                                          
         COMMON /ASYCM/ EASYBC, EASYAC, EASYBA                                  
         COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                      
         COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                      
C                                                                               
      COMMON/INFOCM/CARTNU(NATOM,3),INDEXES(NATOM),                             
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
         COMMON /PRE11CM/ dt2,ret2,betat2,                                      
     +                    y90c1,y90c2,y90c3,y90c4,                              
     +                     y0c1, y0c2, y0c3, y0c4,                              
     +                    ads1,ars1,abs1,                                       
     +                    cds1,crs1,cbs1,                                       
     +                    swa,swr0,swa2,swr20,                                  
     +                    c2a,c2b,c2c,deh2                                      
C                                                                               
      rac2 = 1.d0/dsqrt(2.d0)                                                   
C                                                                               
      r1s=r(1)*r(1)                                                             
      r2s=r(2)*r(2)                                                             
      r3s=r(3)*r(3)                                                             
      rnag=dsqrt(0.5d0*dabs(r1s+r3s-0.5d0*r2s))                                 
      if (rnag.lt.1.d-10) stop 'RNAG < 1.D-10'                                  
C                                                                               
      if (r(2).ne.0.d0) then                                                    
         cstsq=(0.5d0*(r1s-r3s)/(r(2)*rnag))**2                                 
      else                                                                      
         cstsq=0.d0                                                             
      end if                                                                    
C                                                                               
      c1=y90c1+(y0c1-y90c1)*cstsq                                               
      c2=y90c2+(y0c2-y90c2)*cstsq                                               
      c3=y90c3+(y0c3-y90c3)*cstsq                                               
      c4=y90c4+(y0c4-y90c4)*cstsq                                               
C                                                                               
      swt=0.5d0*(1.d0-dtanh(swa*(r(1)-swr0)))*                                  
     +    0.5d0*(1.d0-dtanh(swa*(r(3)-swr0)))*                                  
     +    0.5d0*(1.d0+dtanh(swa2*(r(2)-swr20)))                                 
C                                                                               
      axs1=dexp(-abs1*(r(1)-ars1))                                              
      as1=ads1*axs1*(axs1-2.d0)                                                 
      cxs1=dexp(-cbs1*(r(1)-crs1))                                              
      cs1=cds1*cxs1*(cxs1-2.d0)                                                 
      s1=as1+(cs1-as1)*swt                                                      
C                                                                               
      s2=(149.480d0-59.6557d0*r(2))*dexp(-4.13792d0*r(2))+                      
     +   r2s*(-23.7299d0+3.91747d0*r(2))*dexp(-1.41350d0*r(2))                  
C                                                                               
      axs3=dexp(-abs1*(r(3)-ars1))                                              
      as3=ads1*axs3*(axs3-2.d0)                                                 
      cxs3=dexp(-cbs1*(r(3)-crs1))                                              
      cs3=cds1*cxs3*(cxs3-2.d0)                                                 
      s3=as3+(cs3-as3)*swt                                                      
C                                                                               
      t1=c1*dexp(-c2*r(1)**6)+c3*dexp(-c4*r(1))                                 
C                                                                               
      xt2=dexp(-betat2*(r(2)-ret2))                                             
      t2=dt2*xt2*(xt2-2.d0)                                                     
C                                                                               
      t3=c1*dexp(-c2*r(3)**6)+c3*dexp(-c4*r(3))                                 
C                                                                               
      coul1=0.5d0*(s1+t1)                                                       
      coul2=0.5d0*(s2+t2)                                                       
      coul3=0.5d0*(s3+t3)                                                       
C                                                                               
      exch1=0.5d0*(s1-t1)                                                       
      exch2=0.5d0*(s2-t2)                                                       
      exch3=0.5d0*(s3-t3)                                                       
C                                                                               
      w=(exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2                      
      cplg2=c2a*dexp(-c2b*w-c2c*(r(1)+r(2)+r(3)))                               
C                                                                               
      ENGYGS = (coul1+coul2+coul3+ EZERO(1)  -                                  
     +          rac2*dsqrt(w+cplg2*cplg2))*cevau                                
C                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c  U12 surface - coupling surface                                               
c                                                                               
      subroutine prepot12                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
         PARAMETER (CANGAU = 0.52917706D0)                                      
         PARAMETER (CHARKC = 627.5095D0)                                        
         PARAMETER (cevau=1.d0/27.211611d0)                                     
C                                                                               
         COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                         
         COMMON /PT3CM/ EZERO(ISURF+1)                                          
         COMMON /ASYCM/ EASYBC, EASYAC, EASYBA                                  
         COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                      
         COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                      
C                                                                               
      COMMON/INFOCM/CARTNU(NATOM,3),INDEXES(NATOM),                             
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
         COMMON /PRE11CM/ dt2,ret2,betat2,                                      
     +                    y90c1,y90c2,y90c3,y90c4,                              
     +                     y0c1, y0c2, y0c3, y0c4,                              
     +                    ads1,ars1,abs1,                                       
     +                    cds1,crs1,cbs1,                                       
     +                    swa,swr0,swa2,swr20,                                  
     +                    c2a,c2b,c2c,deh2                                      
C                                                                               
         COMMON /PRE12CM/ eps,ampn,ampp,ampe,ampr,                              
     +                        erln,erlp,erle,erlr,                              
     +                        ersn,ersp,erse,ersr,                              
     +                        rmxn,rmxp,rmxe,rmxr                               
C                                                                               
      rac2 = 1.d0/dsqrt(2.d0)                                                   
C                                                                               
      WRITE(NFLAG(18),100) ampn,ampp,ampe,ampr,erln,erlp,erle,erlr,             
     +             ersn,ersp,erse,ersr,rmxn,rmxp,rmxe,rmxr               1G25T93
C                                                                               
100   format(/,28('='),' NaH2 - U12 coupling ',29('='),                         
     &  /'ampn=',1pe13.5,'  ampp=',e13.5,'  ampe=',e13.5,'  ampr=',e13.5,       
     &  /'erln=',e13.5,'  erlp=',e13.5,'  erle=',e13.5,'  erlr=',e13.5,         
     &  /'ersn=',e13.5,'  ersp=',e13.5,'  erse=',e13.5,'  ersr=',e13.5,         
     &  /'rmxn=',e13.5,'  rmxp=',e13.5,'  rmxe=',e13.5,'  rmxr=',e13.5,         
     &  /78('='))                                                        1G25T93
C                                                                               
      return                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE POT12                                                          
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
         PARAMETER (CANGAU = 0.52917706D0)                                      
         PARAMETER (CHARKC = 627.5095D0)                                        
         PARAMETER (cevau=1.d0/27.211611d0)                                     
C                                                                               
         COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                         
         COMMON /PT3CM/ EZERO(ISURF+1)                                          
         COMMON /ASYCM/ EASYBC, EASYAC, EASYBA                                  
         COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                      
         COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                      
C                                                                               
      COMMON/INFOCM/CARTNU(NATOM,3),INDEXES(NATOM),                             
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
         COMMON /PRE11CM/ dt2,ret2,betat2,                                      
     +                    y90c1,y90c2,y90c3,y90c4,                              
     +                     y0c1, y0c2, y0c3, y0c4,                              
     +                    ads1,ars1,abs1,                                       
     +                    cds1,crs1,cbs1,                                       
     +                    swa,swr0,swa2,swr20,                                  
     +                    c2a,c2b,c2c,deh2                                      
C                                                                               
         COMMON /PRE12CM/ eps,ampn,ampp,ampe,ampr,                              
     +                        erln,erlp,erle,erlr,                              
     +                        ersn,ersp,erse,ersr,                              
     +                        rmxn,rmxp,rmxe,rmxr                               
C                                                                               
      rac2 = 1.d0/dsqrt(2.d0)                                                   
C                                                                               
      amp=ampn+(ampp-ampn)*0.5d0*(dtanh(ampe*(r(2)-ampr))+1.d0)                 
      erl=erln+(erlp-erln)*0.5d0*(dtanh(erle*(r(2)-erlr))+1.d0)                 
      ers=ersn+(ersp-ersn)*0.5d0*(dtanh(erse*(r(2)-ersr))+1.d0)                 
      rmx=rmxn+(rmxp-rmxn)*0.5d0*(dtanh(rmxe*(r(2)-rmxr))+1.d0)                 
      f1=amp*dexp(-(erl+ers/(eps+r(1)**3))*(r(1)-rmx)**2)                       
      f3=amp*dexp(-(erl+ers/(eps+r(3)**3))*(r(3)-rmx)**2)                       
      ang=((r(1)-r(3))/r(2))**2                                                 
C                                                                               
      ENGYIJ(1)=ang*(f1+f3)*cevau                                               
C                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c  U22 surface - upper diabatic surface                                         
c                                                                               
      subroutine prepot22                                                       
      implicit double precision (a-h,o-z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
         PARAMETER (CANGAU = 0.52917706D0)                                      
         PARAMETER (CHARKC = 627.5095D0)                                        
         PARAMETER (cevau=1.d0/27.211611d0)                                     
C                                                                               
         COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                         
         COMMON /PT3CM/ EZERO(ISURF+1)                                          
         COMMON /ASYCM/ EASYBC, EASYAC, EASYBA                                  
         COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                      
         COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                      
C                                                                               
      COMMON/INFOCM/CARTNU(NATOM,3),INDEXES(NATOM),                             
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
         COMMON /PRE11CM/ dt2,ret2,betat2,                                      
     +                    y90c1,y90c2,y90c3,y90c4,                              
     +                     y0c1, y0c2, y0c3, y0c4,                              
     +                    ads1,ars1,abs1,                                       
     +                    cds1,crs1,cbs1,                                       
     +                    swa,swr0,swa2,swr20,                                  
     +                    c2a,c2b,c2c,deh2                                      
C                                                                               
         COMMON /PRE12CM/ eps,ampn,ampp,ampe,ampr,                              
     +                        erln,erlp,erle,erlr,                              
     +                        ersn,ersp,erse,ersr,                              
     +                        rmxn,rmxp,rmxe,rmxr                               
C                                                                               
         COMMON /PRE22CM/ c1s1,c2s1,c3s1,c4s1,c5s1,c6s1,                        
     +                    deh222,una2p,                                         
     +                    c2a22,c2b22,c2c22,cprime22,                           
     +                    y0d1,y0r1,y0b1,y90d1,y90r1,y90b1,                     
     +                    y0d2,y0r2,y0b2,y90d2,y90r2,y90b2,                     
     +                    y0swa,y0swr0,y90swa,y90swr0,                          
     +                    a0g1,b0g1,r0g1,a0g2,b0g2,r0g2,                        
     +                    a0g3,b0g3,r0g3,                                       
     +                    a90g1,b90g1,r90g1,a90g2,b90g2,r90g2,                  
     +                    a90g3,b90g3,r90g3                                     
C                                                                               
      rac2 = 1.d0/dsqrt(2.d0)                                                   
C                                                                               
      WRITE(NFLAG(18),100) y90d1,  y90r1,y90b1,y90d2,y90r2,y90b2,               
     +            y90swa,y90swr0,                                               
     +              y0d1,   y0r1, y0b1, y0d2, y0r2, y0b2,                       
     +             y0swa, y0swr0,                                               
     +              c1s1,   c2s1, c3s1, c4s1, c5s1, c6s1,                       
     +             c2a22,  c2b22,c2c22, a0g1, b0g1, r0g1,                       
     +              a0g2,   b0g2, r0g2, a0g3, b0g3, r0g3,                       
     +             a90g1,  b90g1,r90g1,a90g2,b90g2,r90g2,                       
     +             a90g3,  b90g3,r90g3                                          
C                                                                               
100   format(/,24('='),' NaH2 - U22 potential ',25('='),                        
     &    /'  y90d1 =',1pe14.6,'    y90r1 =',e14.6,'  y90b1 =',e14.6,           
     &    /'  y90d2 =',e14.6,  '    y90r2 =',e14.6,'  y90b2 =',e14.6,           
     &    /' y90swa =',e14.6,  '  y90swr0 =',e14.6,                             
     &    /'   y0d1 =',1pe14.6,'     y0r1 =',e14.6,'   y0b1 =',e14.6,           
     &    /'   y0d2 =',e14.6,  '     y0r2 =',e14.6,'   y0b2 =',e14.6,           
     &    /'  y0swa =',e14.6,  '   y0swr0 =',e14.6,                             
     &    /'   c1s1 =',1pe14.6,'     c2s1 =',e14.6,'   c3s1 =',e14.6,           
     &    /'   c4s1 =',1pe14.6,'     c5s1 =',e14.6,'   c6s1 =',e14.6,           
     &    /'    c2a =',1pe14.6,'      c2b =',e14.6,'    c2c =',e14.6,           
     &    /'   a0g1 =',1pe14.6,'     b0g1 =',e14.6,'   r0g1 =',e14.6,           
     &    /'   a0g2 =',1pe14.6,'     b0g2 =',e14.6,'   r0g2 =',e14.6,           
     &    /'   a0g3 =',1pe14.6,'     b0g3 =',e14.6,'   r0g3 =',e14.6,           
     &    /'  a90g1 =',1pe14.6,'    b90g1 =',e14.6,'  r90g1 =',e14.6,           
     &    /'  a90g2 =',1pe14.6,'    b90g2 =',e14.6,'  r90g2 =',e14.6,           
     &    /'  a90g3 =',1pe14.6,'    b90g3 =',e14.6,'  r90g3 =',e14.6,           
     &    /71('='),/)                                                           
C                                                                               
      return                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE POT22                                                          
      implicit double precision (a-h,o-z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
         PARAMETER (CANGAU = 0.52917706D0)                                      
         PARAMETER (CHARKC = 627.5095D0)                                        
         PARAMETER (cevau=1.d0/27.211611d0)                                     
C                                                                               
         COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                         
         COMMON /PT3CM/ EZERO(ISURF+1)                                          
         COMMON /ASYCM/ EASYBC, EASYAC, EASYBA                                  
         COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                      
         COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                      
C                                                                               
      COMMON/INFOCM/CARTNU(NATOM,3),INDEXES(NATOM),                             
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
         COMMON /PRE11CM/ dt2,ret2,betat2,                                      
     +                    y90c1,y90c2,y90c3,y90c4,                              
     +                     y0c1, y0c2, y0c3, y0c4,                              
     +                    ads1,ars1,abs1,                                       
     +                    cds1,crs1,cbs1,                                       
     +                    swa,swr0,swa2,swr20,                                  
     +                    c2a,c2b,c2c,deh2                                      
C                                                                               
         COMMON /PRE12CM/ eps,ampn,ampp,ampe,ampr,                              
     +                        erln,erlp,erle,erlr,                              
     +                        ersn,ersp,erse,ersr,                              
     +                        rmxn,rmxp,rmxe,rmxr                               
C                                                                               
         COMMON /PRE22CM/ c1s1,c2s1,c3s1,c4s1,c5s1,c6s1,                        
     +                    deh222,una2p,                                         
     +                    c2a22,c2b22,c2c22,cprime22,                           
     +                    y0d1,y0r1,y0b1,y90d1,y90r1,y90b1,                     
     +                    y0d2,y0r2,y0b2,y90d2,y90r2,y90b2,                     
     +                    y0swa,y0swr0,y90swa,y90swr0,                          
     +                    a0g1,b0g1,r0g1,a0g2,b0g2,r0g2,                        
     +                    a0g3,b0g3,r0g3,                                       
     +                    a90g1,b90g1,r90g1,a90g2,b90g2,r90g2,                  
     +                    a90g3,b90g3,r90g3                                     
C                                                                               
      rac2 = 1.d0/dsqrt(2.d0)                                                   
C                                                                               
      r1s=r(1)*r(1)                                                             
      r2s=r(2)*r(2)                                                             
      r3s=r(3)*r(3)                                                             
      rnag=dsqrt(0.5d0*dabs(r1s+r3s-0.5d0*r2s))                                 
      if (rnag.lt.1.d-10) stop 'RNAG < 1.D-10'                                  
      cstsq=(0.5d0*(r1s-r3s)/(r(2)*rnag))**2                                    
C                                                                               
      swa22=y90swa+(y0swa-y90swa)*cstsq                                         
      swr022=y90swr0+(y0swr0-y90swr0)*cstsq                                     
C                                                                               
      s1=(c1s1+c2s1*r(1))*dexp(c3s1*r(1))+                                      
     +   r1s*(c4s1+c5s1*r(1))*dexp(c6s1*r(1))                                   
C                                                                               
      x90t2=dexp(-y90b2*(r(2)-y90r2))                                           
      t902a=y90d2*x90t2*(x90t2-2.d0)                                            
      x0t2=dexp(-y0b2*(r(2)-y0r2))                                              
      t02a=y0d2*x0t2*(x0t2-2.d0)                                                
      t2a=t902a+(t02a-t902a)*cstsq                                              
      swt=0.5d0*(1.d0-dtanh( swa22*(0.5d0*(r(1)+r(3))-swr022)))                 
      s2a=(149.480d0-59.6557d0*r(2))*dexp(-4.13792d0*r(2))+r2s*                 
     +    (-23.7299d0+3.91747d0*r(2))*dexp(-1.41350d0*r(2))+una2p               
      s2b=144.893d0*dexp(-3.85716d0*r(2))+r(2)*(37.5919d0-                      
     +    4.32985d0*r(2)-0.003807d0*r2s)*dexp(-1.52496d0*r(2))                  
      t2=swt*(t2a-s2b)+s2b                                                      
      s2=0.5d0*(s2a+t2-dtanh((s2a-t2)/cprime22)*(s2a-t2))                       
C                                                                               
      s3=(c1s1+c2s1*r(3))*dexp(c3s1*r(3))+                                      
     +   r3s*(c4s1+c5s1*r(3))*dexp(c6s1*r(3))                                   
C                                                                               
      x90t1=dexp(-y90b1*(r(1)-y90r1))                                           
      t901=y90d1*x90t1*(x90t1+2.d0)                                             
      x0t1=dexp(-y0b1*(r(1)-y0r1))                                              
      t01=y0d1*x0t1*(x0t1+2.d0)                                                 
      t1=t901+(t01-t901)*cstsq                                                  
C                                                                               
      x90t3=dexp(-y90b1*(r(3)-y90r1))                                           
      t903=y90d1*x90t3*(x90t3+2.d0)                                             
      x0t3=dexp(-y0b1*(r(3)-y0r1))                                              
      t03=y0d1*x0t3*(x0t3+2.d0)                                                 
      t3=t903+(t03-t903)*cstsq                                                  
C                                                                               
      coul1=0.5d0*(s1+t1)                                                       
      coul2=0.5d0*(s2+t2)                                                       
      coul3=0.5d0*(s3+t3)                                                       
C                                                                               
      exch1=0.5d0*(s1-t1)                                                       
      exch2=0.5d0*(s2-t2)                                                       
      exch3=0.5d0*(s3-t3)                                                       
C                                                                               
      w=(exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2                      
      cplg2=c2a22*dexp(-c2b22*w-c2c22*(r(1)+r(2)+r(3)))                         
C                                                                               
      E=coul1+coul2+coul3+ EZERO(2) - rac2*dsqrt(w+cplg2*cplg2)                 
C                                                                               
      g1 = a0g1*dexp(-b0g1*(rnag-r0g1)**2)                                      
      g2 = a0g2*dexp(-b0g2*(rnag-r0g2)**2)                                      
      g3 = a0g3*dexp(-b0g3*(rnag-r0g3)**2)                                      
      g4 = a90g1*dexp(-b90g1*(rnag-r90g1)**2)                                   
      g5 = a90g2*dexp(-b90g2*(rnag-r90g2)**2)                                   
      g6 = a90g3*dexp(-b90g3*(rnag-r90g3)**2)                                   
C                                                                               
      vdisp = g4+g5+g6+(g1+g2+g3-g4-g5-g6)*cstsq                                
      vdisp = vdisp*dexp(-0.5d0*(r(2)/2.d0)**7)                                 
C                                                                               
      ENGYES(1)=(E+vdisp)*cevau                                                 
C                                                                               
      return                                                                    
      end                                                                       
C                                                                               
      BLOCK DATA PTPACM                                                         
      implicit double precision (a-h,o-z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER (N3ATOM = 75)                                                   
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
         PARAMETER (CANGAU = 0.52917706D0)                                      
         PARAMETER (CHARKC = 627.5095D0)                                        
         PARAMETER (cevau=1.d0/27.211611d0)                                     
C                                                                               
         COMMON /PT1CM/ R(N3ATOM),ENGYGS,DEGSDR(N3ATOM)                         
         COMMON /PT3CM/ EZERO(ISURF+1)                                          
         COMMON /ASYCM/ EASYBC, EASYAC, EASYBA                                  
         COMMON /PT4CM/ ENGYES(ISURF),DEESDR(N3ATOM,ISURF)                      
         COMMON /PT5CM/ ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)                      
C                                                                               
      COMMON/INFOCM/CARTNU(NATOM,3),INDEXES(NATOM),                             
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
         COMMON /PRE11CM/ dt2,ret2,betat2,                                      
     +                    y90c1,y90c2,y90c3,y90c4,                              
     +                     y0c1, y0c2, y0c3, y0c4,                              
     +                    ads1,ars1,abs1,                                       
     +                    cds1,crs1,cbs1,                                       
     +                    swa,swr0,swa2,swr20,                                  
     +                    c2a,c2b,c2c,deh2                                      
C                                                                               
         COMMON /PRE12CM/ eps,ampn,ampp,ampe,ampr,                              
     +                        erln,erlp,erle,erlr,                              
     +                        ersn,ersp,erse,ersr,                              
     +                        rmxn,rmxp,rmxe,rmxr                               
C                                                                               
         COMMON /PRE22CM/ c1s1,c2s1,c3s1,c4s1,c5s1,c6s1,                        
     +                    deh222,una2p,                                         
     +                    c2a22,c2b22,c2c22,cprime22,                           
     +                    y0d1,y0r1,y0b1,y90d1,y90r1,y90b1,                     
     +                    y0d2,y0r2,y0b2,y90d2,y90r2,y90b2,                     
     +                    y0swa,y0swr0,y90swa,y90swr0,                          
     +                    a0g1,b0g1,r0g1,a0g2,b0g2,r0g2,                        
     +                    a0g3,b0g3,r0g3,                                       
     +                    a90g1,b90g1,r90g1,a90g2,b90g2,r90g2,                  
     +                    a90g3,b90g3,r90g3                                     
C                                                                               
      DATA NASURF /1,35*0/                                                      
      DATA NDER /0/                                                             
      DATA NFLAG /1,1,15*0,6,0,0/                                               
C                                                                               
      DATA ANUZERO /0.0D0/                                                      
      DATA ICARTR,MSURF,MDER/3,1,0/                                             
      DATA NULBL /25*0/                                                         
      DATA NATOMS /3/                                                           
C                                                                               
      DATA dt2    /2.885656d-2/                                                 
      DATA ret2   /1.237624d0/                                                  
      DATA betat2 /1.038036d+1/                                                 
      DATA y90c1  /6.279847d0/                                                  
      DATA y90c2  /3.849003d-3/                                                 
      DATA y90c3  /5.546257d0/                                                  
      DATA y90c4  /5.573967d-1/                                                 
      DATA y0c1   /1.060311d+1/                                                 
      DATA y0c2   /6.058605d-3/                                                 
      DATA y0c3   /4.573378d0/                                                  
      DATA y0c4   /5.182708d-1/                                                 
      DATA ads1   /0.08d0/                                                      
      DATA ars1   /5.015440d0/                                                  
      DATA abs1   /7.536929d-1/                                                 
      DATA cds1   /1.323042d0/                                                  
      DATA crs1   /4.00d0/                                                      
      DATA cbs1   /5.4070d-1/                                                   
      DATA swa    /0.7d0/                                                       
      DATA swr0   /6.5d0/                                                       
      DATA swa2   /1.1d0/                                                       
      DATA swr20  /1.9d0/                                                       
      DATA c2a    /1.5d0/                                                       
      DATA c2b    /1.0d0/                                                       
      DATA c2c    /0.15d0/                                                      
      DATA deh2   /4.7469d0/                                                    
C                                                                               
      DATA eps/1.d-4/                                                           
C                                                                               
      DATA ampn  / 2.35737419d-1/                                               
      DATA ampp  / 1.51495000d+0/                                               
      DATA ampe  / 5.02447648d-1/                                               
      DATA ampr  / 4.05991868d+0/                                               
      DATA erln  / 2.30065420d+2/                                               
      DATA erlp  / 7.10000000d-2/                                               
      DATA erle  / 3.46985415d+0/                                               
      DATA erlr  / 8.95373168d-1/                                               
      DATA ersn  / 2.00000000d+0/                                               
      DATA ersp  / 1.16830000d+1/                                               
      DATA erse  / 4.54595148d+1/                                               
      DATA ersr  /-2.38192975d+0/                                               
      DATA rmxn  / 2.33484458d+0/                                               
      DATA rmxp  / 2.62222000d+0/                                               
      DATA rmxe  / 2.04523701d+1/                                               
      DATA rmxr  / 1.98941949d+0/                                               
C                                                                               
      DATA c1s1  / 2.614452d+2/                                                 
      DATA c2s1  /-3.155166d-4/                                                 
      DATA c3s1  /-1.785467d+0/                                                 
      DATA c4s1  /-3.548074d+0/                                                 
      DATA c5s1  / 3.356766d-1/                                                 
      DATA c6s1  /-7.621511d-1/                                                 
C                                                                               
      DATA deh222    /  4.7469d0/                                               
C                                                                               
      DATA una2p     /  2.1037d0/                                               
C                                                                               
      DATA c2a22     /  0.8d0/                                                  
      DATA c2b22     /  0.8d0/                                                  
      DATA c2c22     /  0.15d0/                                                 
C                                                                               
      DATA cprime22  /  1.0d-1/                                                 
C                                                                               
      DATA y0d1      /  2.417800d-1/                                            
      DATA y0r1      /  3.708055d0/                                             
      DATA y0b1      /  6.960127d-1/                                            
C                                                                               
      DATA y90d1     /  1.372065d-2/                                            
      DATA y90r1     /  7.233552d0/                                             
      DATA y90b1     /  5.127734d-1/                                            
C                                                                               
      DATA y0d2      /  8.870172d-1/                                            
      DATA y0r2      /  2.300265d0/                                             
      DATA y0b2      /  8.497374d-1/                                            
C                                                                               
      DATA y90d2     /  1.106187d0/                                             
      DATA y90r2     /  2.002516d0/                                             
      DATA y90b2     /  9.514887d-1/                                            
C                                                                               
      DATA y0swa     /  0.614393d0/                                             
      DATA y0swr0    /  9.006323d0/                                             
      DATA y90swa    /  0.872116d0/                                             
      DATA y90swr0   /  8.194467d0/                                             
C                                                                               
      DATA a0g1      / -8.210245d-2/                                            
      DATA b0g1      /  4.658698d0/                                             
      DATA r0g1      /  4.552860d0/                                             
C                                                                               
      DATA a0g2      /  3.546427d-2/                                            
      DATA b0g2      /  6.456787d-1/                                            
      DATA r0g2      /  6.332225d0/                                             
C                                                                               
      DATA a0g3      / -1.562327d-2/                                            
      DATA b0g3      /  1.170537d-1/                                            
      DATA r0g3      /  1.097310d+1/                                            
C                                                                               
      DATA a90g1     / -9.929297d-2/                                            
      DATA b90g1     /  4.367220d0/                                             
      DATA r90g1     /  4.402447d0/                                             
C                                                                               
      DATA a90g2     /  4.696181d-2/                                            
      DATA b90g2     /  5.470459d-1/                                            
      DATA r90g2     /  6.289526d0/                                             
C                                                                               
      DATA a90g3     / -1.219900d-2/                                            
      DATA b90g3     /  1.002291d-2/                                            
      DATA r90g3     /  1.118725d+1/                                            
C                                                                               
      END                                                                       
