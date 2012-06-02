      subroutine prepot                                                         
C                                                                               
C   System:          NaFH                                                       
C   Common name:     N3                                                         
C   Number of derivatives: 0                                                    
C   Number of electronic surfaces: 2                                            
C                                                                               
C   Protocol:                                                                   
C                                                                               
C      PREPOT - initializes the potential's variables                           
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
         COMMON /PRE11CM/ scale,cprime,dt2,ret2,betat2,                         
     +                    c1,c2,c3,c4,                                          
     +                    ads1,ars1,abs1,                                       
     +                    cds1,crs1,cbs1,                                       
     +                    swa,swr0,swa2,swr20,                                  
     +                    c2a,c2b,c2c,deh2,cfm                                  
C                                                                               
         COMMON /PRE22CM/ c1s1,c2s1,c3s1,c4s1,c5s1,c6s1,                        
     +                    dnaht1,dnaht2,dnaht3,dnaht4,                          
     +                    dnaht5,dnaht6,dnaht7,                                 
     +                    dnaht11,dnaht61,deh222,una2p,                         
     +                    c2a22,c2b22,c2c22,cprime22,                           
     +                    y0d1,y0r1,y0b1,y90d1,y90r1,y90b1,                     
     +                    y0d2,y0r2,y0b2,y90d2,y90r2,y90b2,                     
     +                    swa22,swr022,                                         
     +                    a0g1,b0g1,r0g1,a0g2,b0g2,r0g2,                        
     +                    a0g3,b0g3,r0g3,                                       
     +                    a90g1,b90g1,r90g1,a90g2,b90g2,r90g2,                  
     +                    a90g3,b90g3,r90g3,                                    
     +                    h2t1,h2t2,h2t3,h2t4,h2t5,h2t6,h2t7                    
C                                                                               
      IF(NATOMS.GT.25) THEN                                                     
         WRITE(NFLAG(18),1111)                                                  
 1111    FORMAT(2X,'STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS')            
         STOP                                                                   
      END IF                                                                    
C                                                                               
        WRITE(NFLAG(18),*) 'PREPOT called for MHF N3 potential '                
        write(NFLAG(18),*) 'energy surface.'                                    
        WRITE(NFLAG(18),*) 'Development version.'                               
        WRITE(NFLAG(18),*) 'Last modified 7 October 1997.'                      
        WRITE(NFLAG(18),*)                                                      
                                                                                
        call prepot11                                                           
        call prepot12                                                           
        call prepot22                                                           
C                                                                               
C                                                                               
      EZERO(1)=DEH2                                                             
      EZERO(2)=DEH222                                                           
C                                                                               
       DO I=1,5                                                                 
          REF(I) = ' '                                                          
       END DO                                                                   
C                                                                               
       REF(1)='No Reference'                                                    
C                                                                               
      INDEXES(1) = 11                                                           
      INDEXES(2) = 9                                                            
      INDEXES(3) = 1                                                            
C                                                                               
      IRCTNT=2                                                                  
C                                                                               
      CALL POTINFO                                                              
C                                                                               
      CALL ANCVRT                                                               
C                                                                               
      return                                                                    
                                                                                
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
      if (NASURF(1,1).NE.0) call pot11                                          
      if (NASURF(2,2).NE.0) call pot22                                          
      if (NASURF(1,2)+NASURF(2,1).NE.0) call pot12                              
                                                                                
      CALL EUNITZERO                                                            
      IF(NDER.NE.0) THEN                                                        
         CALL RTOCART                                                           
         CALL DEDCOU                                                            
      ENDIF                                                                     
C                                                                               
      return                                                                    
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
         PARAMETER (cevau=1.d0/27.2113961d0)                                    
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
         COMMON /PRE11CM/ scale,cprime,dt2,ret2,betat2,                         
     +                    c1,c2,c3,c4,                                          
     +                    ads1,ars1,abs1,                                       
     +                    cds1,crs1,cbs1,                                       
     +                    swa,swr0,swa2,swr20,                                  
     +                    c2a,c2b,c2c,deh2,cfm                                  
C                                                                               
      rac2 = 1.d0/dsqrt(2.d0)                                                   
C                                                                               
      WRITE(NFLAG(18),100) dt2,ret2,betat2,c1,c2,c3,c4,                         
     +             ads1,ars1,abs1,cds1,crs1,cbs1,                               
     +             swa,swr0,swa2,swr20,c2a,c2b,c2c                              
C                                                                               
100   format(/24('='),' MHF - U11 potential ',24('='),                          
     &  /'   dt2 =',1pe14.6,2x,'  ret2 =',e14.6,2x,'betat2 =',e14.6,            
     &  /' c1 =',e14.6,2x,' c2 =',e14.6,2x,' c3 =',e14.6,                       
     &  /' c4 =',e14.6,2x,                                                      
     &  /'  ads1 =',e14.6,                                                      
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
         PARAMETER (cevau=1.d0/27.2113961d0)                                    
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
         COMMON /PRE11CM/ scale,cprime,dt2,ret2,betat2,                         
     +                    c1,c2,c3,c4,                                          
     +                    ads1,ars1,abs1,                                       
     +                    cds1,crs1,cbs1,                                       
     +                    swa,swr0,swa2,swr20,                                  
     +                    c2a,c2b,c2c,deh2,cfm                                  
C                                                                               
      rac2 = 1.d0/dsqrt(2.d0)                                                   
C                                                                               
C***********************************************************************        
C                                                                               
      r2s=r(2)*r(2)                                                             
      rr1 = r(1)*scale                                                          
      rr3 = r(3)*scale                                                          
      rr1s = rr1**2                                                             
      rr3s = rr3**2                                                             
C                                                                               
C*****SWITCHING FUNCTION FOR THE SINGLETS BELOW, AND ITS DERIVATIVES:           
C                                                                               
      hlp1 = 1.d0-dtanh(swa*(rr1-swr0))                                         
      hlp2 = 1.d0+dtanh(swa2*(r(2)-swr20))                                      
      hlp3 = 1.d0-dtanh(swa*(rr3-swr0))                                         
      swt=0.125d0*hlp1*hlp2*hlp3                                                
C                                                                               
C*****SINGLETS AND THEIR DERIVATIVES:                                           
C*****MH SINGLET:                                                               
C                                                                               
      axs1=dexp(-0.5d0*abs1*(rr1-ars1))                                         
      as1=0.5d0*ads1*axs1*(axs1-2.d0)                                           
                                                                                
      cxs1=dexp(-2.80d0*abs1*(rr1-3.500d0))                                     
      cs1=0.35d0*cxs1*(cxs1+2.d0)                                               
                                                                                
      s1=as1+(cs1-as1)*swt                                                      
C                                                                               
C*****H2 SINGLET, AND ITS DERIVATIVES:  (ORIGINAL)                              
C                                                                               
      s2a=(109.384d0+463.033d0*r(2))*dexp(-6.240d0*r(2))-                       
     +     r2s*(13.654d0+20.717d0*r(2))*dexp(-2.032d0*r(2))+                    
     +     1.34357d0                                                            
C                                                                               
      hlp1 = dexp(-8.d0*(r(2)-1.411655d0))                                      
      hlp2 = 0.2d0*dexp(-2.032d0*r(2))                                          
      s2h = hlp1 - hlp2                                                         
C                                                                               
      hlp=dtanh((s2a-s2h)/cprime)                                               
C                                                                               
c      hlp1=(s2a-s2h)/(cprime*dcosh((s2a-s2h)/cprime)**2)                       
C                                                                               
      s2a=0.5d0*(s2a+s2h-hlp*(s2a-s2h))                                         
C                                                                               
      a = 1.00075d0                                                             
      dn = 3.39687d0                                                            
C                                                                               
      hlp = dexp(-a*(r(2)-1.42d0))                                              
      s2i = dn*hlp*(hlp-2.d0)                                                   
                                                                                
      hlp = 0.5d0*(r(1)+r(3)-8.d0)                                              
      hlp = 0.5d0*(1.d0-dtanh(hlp))                                             
      s2 = s2a+(s2i-s2a)*hlp                                                    
C                                                                               
C*****MF SINGLET:                                                               
C                                                                               
      axs3=dexp(-1.3d0*abs1*(rr3-ars1))                                         
      as3=ads1*axs3*(axs3-2.d0)                                                 
                                                                                
      cxs3=dexp(-cbs1*(rr3-crs1))                                               
      cs3=cds1*cxs3*(cxs3-2.d0)                                                 
                                                                                
      s3=as3+(cs3-as3)*swt                                                      
C                                                                               
C*****TRIPLETS AND THEIR FIRST DERIVATIVES:                                     
C     MH TRIPLET                                                                
C                                                                               
      hlp1=dexp(-1.10d0*c2*(rr1-0.6d0))                                         
      hlp2=dexp(-1.65d0*c4*(rr1-0.6d0))                                         
      t1=c1*hlp1+c3*hlp2                                                        
C                                                                               
C     HH TRIPLET                                                                
C                                                                               
      xt2=dexp(-betat2*(r(2)-ret2))                                             
      t2=dt2*xt2*(xt2-2.d0)                                                     
C                                                                               
C     MF TRIPLET                                                                
C                                                                               
      hlp1=dexp(-c2*(rr3-0.55d0))                                               
      hlp2=dexp(-2.5d0*c4*(rr3-0.55d0))                                         
      t3=c1*hlp1+c3*hlp2                                                        
C                                                                               
C*****COULOMB INTEGRALS AND THEIR FIRST DERIVATIVES:                            
C                                                                               
      coul1=0.5d0*(s1+t1)                                                       
      coul2=0.5d0*(s2+t2)                                                       
      coul3=0.5d0*(s3+t3)                                                       
C                                                                               
C*****EXCHANGE INTEGRALS AND THEIR FIRST DERIVATIVES:                           
C                                                                               
      exch1=0.5d0*(s1-t1)                                                       
      exch2=0.5d0*(s2-t2)                                                       
      exch3=0.5d0*(s3-t3)                                                       
C                                                                               
C***********************************************************************        
C                                                                               
      w=(exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2                      
      cplg2=c2a*dexp(-c2b*w-c2c*(rr1+r(2)+rr3))                                 
      rootw = rac2*dsqrt(w+cplg2*cplg2)                                         
      ENGYGS=((coul1+coul2+coul3+ EZERO(1)  -rootw)*                            
     +            cevau-cfm)                                                    
C                                                                               
C***********************************************************************        
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
         PARAMETER (cevau=1.d0/27.2113961d0)                                    
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
         COMMON /PRE11CM/ scale,cprime,dt2,ret2,betat2,                         
     +                    c1,c2,c3,c4,                                          
     +                    ads1,ars1,abs1,                                       
     +                    cds1,crs1,cbs1,                                       
     +                    swa,swr0,swa2,swr20,                                  
     +                    c2a,c2b,c2c,deh2,cfm                                  
C                                                                               
         COMMON /PRE12CM/ ahf,anaf,anah,r21,r23,r13,                            
     +                    ccMF(4),ccMH(4)                                       
C                                                                               
      rac2 = 1.d0/dsqrt(2.d0)                                                   
C                                                                               
      WRITE(NFLAG(18),100)                                                      
      WRITE(NFLAG(18),101)(ccMF(i),i=1,4)                                       
      WRITE(NFLAG(18),102)(ccMH(i),i=1,4)                                       
      WRITE(NFLAG(18),106)ahf                                                   
      WRITE(NFLAG(18),107)anaf                                                  
      WRITE(NFLAG(18),108)anah                                                  
      WRITE(NFLAG(18),109)r21                                                   
      WRITE(NFLAG(18),110)r23                                                   
      WRITE(NFLAG(18),111)r13                                                   
C                                                                               
  100 format(/25('='),' MFH - U12 coupling ',25('='))                           
  101 format(2x,'ccMF = ',5x,4e14.6/10x,2e14.6)                                 
  102 format(2x,'ccMH = ',4x,4e14.6/11x,2e14.6)                                 
  106 format(2x,'ahf = ',7x,3e14.6)                                             
  107 format(2x,'anaf = ',6x,3e14.6)                                            
  108 format(2x,'anah = ',6x,3e14.6)                                            
  109 format(2x,'r21 = ',7x,3e14.6)                                             
  110 format(2x,'r23 = ',7x,3e14.6)                                             
  111 format(2x,'r13 = ',7x,3e14.6,/,71('='),/)                                 
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
         PARAMETER (cevau=1.d0/27.2113961d0)                                    
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
         COMMON /PRE11CM/ scale,cprime,dt2,ret2,betat2,                         
     +                    c1,c2,c3,c4,                                          
     +                    ads1,ars1,abs1,                                       
     +                    cds1,crs1,cbs1,                                       
     +                    swa,swr0,swa2,swr20,                                  
     +                    c2a,c2b,c2c,deh2,cfm                                  
C                                                                               
         COMMON /PRE12CM/ ahf,anaf,anah,r21,r23,r13,                            
     +                    ccMF(4),ccMH(4)                                       
C                                                                               
      rac2 = 1.d0/dsqrt(2.d0)                                                   
C                                                                               
         hlp = ccMF(3)/(1.d0+ccMF(4)*r(3))                                      
         unaf = dexp(-ccMF(1)-ccMF(2)*r(3)- hlp*r(3))                           
         unaf = unaf*27.2113961d0                                               
         hlp = ccMH(3)/(1.d0+ccMH(4)*r(1))                                      
         unah = dexp(-ccMH(1)-ccMH(2)*r(1)- hlp*r(1))                           
C                                                                               
C-----------------------------------------------------------------------        
C                                                                               
         hlp = anah*r(1)-anaf*r(3)-r13                                          
         sw1 = 0.5d0*(1.d0+dtanh(hlp))                                          
         hlp = ahf*r(2)-anah*r(1)-r21                                           
         sw21 = 0.5d0*(1.d0+dtanh(hlp))                                         
         hlp = ahf*r(2)-anaf*r(3)-r23                                           
         sw23 = 0.5d0*(1.d0+dtanh(hlp))                                         
                                                                                
         E = (unaf*sw1*sw23+unah*(1.d0-sw1)*sw21)*cevau                         
C                                                                               
C     MODIFICATION: coupling cut-off at large H-H separations                   
C                                                                               
         gmc = 0.3d0                                                            
         amc = 0.40d0                                                           
         rmc = 2.7d0                                                            
         hlp= gmc*(1.d0-dtanh(amc*(r(2)-rmc)))                                  
         ENGYIJ(1) = E*hlp                                                      
C                                                                               
C***********************************************************************        
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
         PARAMETER (cevau=1.d0/27.2113961d0)                                    
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
         COMMON /PRE11CM/ scale,cprime,dt2,ret2,betat2,                         
     +                    c1,c2,c3,c4,                                          
     +                    ads1,ars1,abs1,                                       
     +                    cds1,crs1,cbs1,                                       
     +                    swa,swr0,swa2,swr20,                                  
     +                    c2a,c2b,c2c,deh2,cfm                                  
C                                                                               
         COMMON /PRE12CM/ ahf,anaf,anah,r21,r23,r13,                            
     +                    ccMF(4),ccMH(4)                                       
C                                                                               
         COMMON /PRE22CM/ c1s1,c2s1,c3s1,c4s1,c5s1,c6s1,                        
     +                    dnaht1,dnaht2,dnaht3,dnaht4,                          
     +                    dnaht5,dnaht6,dnaht7,                                 
     +                    dnaht11,dnaht61,deh222,una2p,                         
     +                    c2a22,c2b22,c2c22,cprime22,                           
     +                    y0d1,y0r1,y0b1,y90d1,y90r1,y90b1,                     
     +                    y0d2,y0r2,y0b2,y90d2,y90r2,y90b2,                     
     +                    swa22,swr022,                                         
     +                    a0g1,b0g1,r0g1,a0g2,b0g2,r0g2,                        
     +                    a0g3,b0g3,r0g3,                                       
     +                    a90g1,b90g1,r90g1,a90g2,b90g2,r90g2,                  
     +                    a90g3,b90g3,r90g3,                                    
     +                    h2t1,h2t2,h2t3,h2t4,h2t5,h2t6,h2t7                    
C                                                                               
C      write(6,100)  y90d1,y90r1,y90b1,y90d2,y90r2,y90b2,y90swa,y90swr0,        
C     +               y0d1, y0r1, y0b1, y0d2, y0r2, y0b2, y0swa, y0swr0,        
C     +               c1s1, c2s1, c3s1, c4s1, c5s1, c6s1,                       
C     +              c2a22,c2b22,c2c22,                                         
C     +               a0g1, b0g1, r0g1, a0g2, b0g2, r0g2, a0g3, b0g3, r0g3,     
C     +              a90g1,b90g1,r90g1,a90g2,b90g2,r90g2,a90g3,b90g3,r90g3      
C                                                                               
      rac2 = 1.d0/dsqrt(2.d0)                                                   
C                                                                               
      WRITE(NFLAG(18),100)  y90d1,y90r1,y90b1,y90d2,y90r2,y90b2,swa22,          
     +                      swr022,y0d1,y0r1,y0b1,y0d2,y0r2,y0b2,swa22,         
     +                      swr022,c1s1, c2s1, c3s1, c4s1, c5s1, c6s1,          
     +                      c2a22,c2b22,c2c22,                                  
     +                      a0g1, b0g1, r0g1, a0g2, b0g2, r0g2,                 
     +                      a0g3, b0g3, r0g3,                                   
     +                      a90g1,b90g1,r90g1,a90g2,b90g2,r90g2,                
     +                      a90g3,b90g3,r90g3                                   
C                                                                               
100   format(/24('='),' MHF - U22 potential ',25('='),                          
     &    /'  y90d1 =',1pe14.6,'    y90r1 =',e14.6,'  y90b1 =',e14.6,           
     &    /'  y90d2 =',e14.6,  '    y90r2 =',e14.6,'  y90b2 =',e14.6,           
     &    /'  swa22 =',e14.6,  '   swr022 =',e14.6,                             
     &    /'   y0d1 =',1pe14.6,'     y0r1 =',e14.6,'   y0b1 =',e14.6,           
     &    /'   y0d2 =',e14.6,  '     y0r2 =',e14.6,'   y0b2 =',e14.6,           
     &    /'  swa22 =',e14.6,  '   swr022 =',e14.6,                             
     &    /'   c1s1 =',1pe14.6,'     c2s1 =',e14.6,'   c3s1 =',e14.6,           
     &    /'   c4s1 =',1pe14.6,'     c5s1 =',e14.6,'   c6s1 =',e14.6,           
     &    /'    c2a =',1pe14.6,'      c2b =',e14.6,'    c2c =',e14.6,           
     &    /'   a0g1 =',1pe14.6,'     b0g1 =',e14.6,'   r0g1 =',e14.6,           
     &    /'   a0g2 =',1pe14.6,'     b0g2 =',e14.6,'   r0g2 =',e14.6,           
     &    /'   a0g3 =',1pe14.6,'     b0g3 =',e14.6,'   r0g3 =',e14.6,           
     &    /'  a90g1 =',1pe14.6,'    b90g1 =',e14.6,'  r90g1 =',e14.6,           
     &    /'  a90g2 =',1pe14.6,'    b90g2 =',e14.6,'  r90g2 =',e14.6,           
     &    /'  a90g3 =',1pe14.6,'    b90g3 =',e14.6,'  r90g3 =',e14.6,           
     &    /71('='))                                                             
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
         PARAMETER (cevau=1.d0/27.2113961d0)                                    
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
         COMMON /PRE11CM/ scale,cprime,dt2,ret2,betat2,                         
     +                    c1,c2,c3,c4,                                          
     +                    ads1,ars1,abs1,                                       
     +                    cds1,crs1,cbs1,                                       
     +                    swa,swr0,swa2,swr20,                                  
     +                    c2a,c2b,c2c,deh2,cfm                                  
C                                                                               
         COMMON /PRE12CM/ ahf,anaf,anah,r21,r23,r13,                            
     +                    ccMF(4),ccMH(4)                                       
C                                                                               
         COMMON /PRE22CM/ c1s1,c2s1,c3s1,c4s1,c5s1,c6s1,                        
     +                    dnaht1,dnaht2,dnaht3,dnaht4,                          
     +                    dnaht5,dnaht6,dnaht7,                                 
     +                    dnaht11,dnaht61,deh222,una2p,                         
     +                    c2a22,c2b22,c2c22,cprime22,                           
     +                    y0d1,y0r1,y0b1,y90d1,y90r1,y90b1,                     
     +                    y0d2,y0r2,y0b2,y90d2,y90r2,y90b2,                     
     +                    swa22,swr022,                                         
     +                    a0g1,b0g1,r0g1,a0g2,b0g2,r0g2,                        
     +                    a0g3,b0g3,r0g3,                                       
     +                    a90g1,b90g1,r90g1,a90g2,b90g2,r90g2,                  
     +                    a90g3,b90g3,r90g3,                                    
     +                    h2t1,h2t2,h2t3,h2t4,h2t5,h2t6,h2t7                    
C                                                                               
      rac2 = 1.d0/dsqrt(2.d0)                                                   
C                                                                               
C***********************************************************************        
C                                                                               
      r1s=r(1)*r(1)                                                             
      r2s=r(2)*r(2)                                                             
      r3s=r(3)*r(3)                                                             
      rr1 = r(1)*scale                                                          
      rr3 = r(3)*scale                                                          
      rr1s = rr1**2                                                             
      rr3s = rr3**2                                                             
C                                                                               
C*****SWITCHING FUNCTION FOR THE H2 TRIPLET, AND ITS DERIVATIVES:               
C                                                                               
      swt=0.5d0*(1.d0-dtanh( swa22*(0.5d0*(rr1+rr3)-swr022)))                   
C                                                                               
C*****SINGLETS AND THEIR DERIVATIVES:                                           
C*****MH SINGLET:                                                               
C                                                                               
      raaa = 2.700d0                                                            
      aaa = 2.40d0                                                              
      aaa1 = 2.60d0                                                             
      daa = 1.80d0                                                              
      hlpNAH1 = dexp(-aaa1*(r(1)-raaa))                                         
      hlpNAH11 = dexp(-1.65d0*aaa*(r(1)-raaa))                                  
      s1 = 0.30d0*daa*(hlpNAH11-2.d0*hlpNAH1)                                   
C                                                                               
C*****H2 SINGLET                                                                
C                                                                               
      s2a=(109.384d0+463.033d0*r(2))*dexp(-6.240d0*r(2))-                       
     +    r2s*(13.654d0+20.717d0*r(2))*dexp(-2.032d0*r(2))+                     
     +    una2p                                                                 
C                                                                               
      b = 1.4d0                                                                 
      a = 2.d0*b                                                                
      r00 = 1.411655d0                                                          
      dn = -4.74044d0/(1.d0-a/b)                                                
      hlp1 = dexp(-a*(r(2)-r00))                                                
      hlp2 = dexp(-b*(r(2)-r00))                                                
      s2i = dn*(hlp1-(a/b)*hlp2)+una2p                                          
C                                                                               
      hlp = 0.5d0*(r(1)+r(3)-7.0d0)                                             
      hlp = 0.5d0*(1.d0-dtanh(hlp))                                             
      s2a = s2a+(s2i-s2a)*hlp                                                   
C                                                                               
C*****MF SINGLET:                                                               
C                                                                               
      hlpNAH3 = dexp(-aaa*(r(3)-raaa))                                          
      hlpNAH33 = dexp(-1.65d0*aaa*(r(3)-raaa))                                  
      s3 = daa*(hlpNAH33-2.d0*hlpNAH3)                                          
C                                                                               
C*****TRIPLETS AND THEIR DERIVATIVES:                                           
C*****H2 TRIPLET (interaction region):                                          
C                                                                               
      t2a=0.5d0*(dexp(-h2t2*(r(2)-h2t1))+(h2t3+h2t5*r(2)+                       
     +    h2t6*r2s+h2t7*r2s*r(2))*dexp(-h2t4*r(2)))                             
C                                                                               
C*****H2 TRIPLET (asymptotic)                                                   
C                                                                               
      s2b=147.214d0*dexp(-3.915d0*r(2))+r(2)*(41.275d0+10.505d0*                
     +    r(2)+4.408d0*r2s)*dexp(-2.032d0*r(2))                                 
C                                                                               
C*****H2 TRIPLET: TOTAL (switched between the asymptotic and the ineraction one)
C                                                                               
      t2=swt*(t2a-s2b)+s2b                                                      
C                                                                               
C*****TRIPLETS AND THEIR DERIVATIVES:  (!!!!!MH TRIPLET)                        
C                                                                               
      t1=3.d0*(dexp(-dnaht2*(rr1-dnaht11))+                                     
     +   (dnaht7+dnaht3*rr1+dnaht4*rr1s+dnaht5*rr1*rr1s)*                       
     +   dexp(-dnaht61*rr1))                                                    
C                                                                               
C   !!!!!MF TRIPLET: (should be the same as the first one, t1)                  
C                                                                               
      t3=3.d0*(dexp(-dnaht2*(rr3-dnaht1))+                                      
     +   (dnaht7+dnaht3*rr3+dnaht4*rr3s+dnaht5*rr3*rr3s)*                       
     +   dexp(-dnaht6*rr3))                                                     
C                                                                               
C*****TOTAL H2 SINGLET:                                                         
C     !!!!!H2 SINGLET: total (switched between the true singlet and a triplet   
C                        in the asymptotic region at large HH distances)        
C                                                                               
      hlp=dtanh((s2a-t2)/cprime22)                                              
      s2=0.5d0*(s2a+t2-hlp*(s2a-t2))                                            
C                                                                               
C*****COULOMB INTEGRALS AND THEIR FIRST DERIVATIVES:                            
C                                                                               
      coul1=0.5d0*(s1+t1)                                                       
      coul2=0.5d0*(s2+t2)                                                       
      coul3=0.5d0*(s3+t3)                                                       
C                                                                               
C*****EXCHANGE INTEGRALS AND THEIR FIRST DERIVATIVES:                           
C                                                                               
      exch1=0.5d0*(s1-t1)                                                       
      exch2=0.5d0*(s2-t2)                                                       
      exch3=0.5d0*(s3-t3)                                                       
C                                                                               
C***********************************************************************        
C                                                                               
      w=(exch1-exch2)**2+(exch2-exch3)**2+(exch3-exch1)**2                      
      cplg2=c2a22*dexp(-c2b22*w-c2c22*(rr1+r(2)+rr3))                           
      rootw = rac2*dsqrt(w+cplg2*cplg2)                                         
      E = coul1+coul2+coul3+ EZERO(2)  -rootw                                   
C                                                                               
c***********************************************************************        
C                                                                               
      rf1 = 3.00d0                                                              
      af1 = 0.60d0                                                              
      fff = 1.d0 - 0.5d0*(1.d0+dtanh(af1*(r(3)-rf1)))                           
      ENGYES(1)=((E-1.35d0)/(1.00d0+1.80d0*fff))*cevau                          
C                                                                               
C***********************************************************************        
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
         PARAMETER (cevau=1.d0/27.2113961d0)                                    
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
         COMMON /PRE11CM/ scale,cprime,dt2,ret2,betat2,                         
     +                    c1,c2,c3,c4,                                          
     +                    ads1,ars1,abs1,                                       
     +                    cds1,crs1,cbs1,                                       
     +                    swa,swr0,swa2,swr20,                                  
     +                    c2a,c2b,c2c,deh2,cfm                                  
C                                                                               
         COMMON /PRE12CM/ ahf,anaf,anah,r21,r23,r13,                            
     +                    ccMF(4),ccMH(4)                                       
C                                                                               
C                                                                               
         COMMON /PRE22CM/ c1s1,c2s1,c3s1,c4s1,c5s1,c6s1,                        
     +                    dnaht1,dnaht2,dnaht3,dnaht4,                          
     +                    dnaht5,dnaht6,dnaht7,                                 
     +                    dnaht11,dnaht61,deh222,una2p,                         
     +                    c2a22,c2b22,c2c22,cprime22,                           
     +                    y0d1,y0r1,y0b1,y90d1,y90r1,y90b1,                     
     +                    y0d2,y0r2,y0b2,y90d2,y90r2,y90b2,                     
     +                    swa22,swr022,                                         
     +                    a0g1,b0g1,r0g1,a0g2,b0g2,r0g2,                        
     +                    a0g3,b0g3,r0g3,                                       
     +                    a90g1,b90g1,r90g1,a90g2,b90g2,r90g2,                  
     +                    a90g3,b90g3,r90g3,                                    
     +                    h2t1,h2t2,h2t3,h2t4,h2t5,h2t6,h2t7                    
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
      DATA scale  /1.5d0/                                                       
C                                                                               
      DATA cprime /4.0d-1/                                                      
      DATA dt2    /2.885656d-2/                                                 
      DATA ret2   /1.237624d0/                                                  
      DATA betat2 /1.038036d+1/                                                 
      DATA c1     /1255969.4d0/                                                 
      DATA c2     /4.7d0/                                                       
      DATA c3     /5.546257d0/                                                  
      DATA c4     /5.573967d-1/                                                 
      DATA ads1   /0.08d0/                                                      
      DATA ars1   /4.0d0/                                                       
      DATA abs1   /1.55d0/                                                      
      DATA cds1   /2.40d0/                                                      
      DATA crs1   /4.15d0/                                                      
      DATA cbs1   /1.20d0/                                                      
      DATA swa    /0.7d0/                                                       
      DATA swr0   /6.5d0/                                                       
      DATA swa2   /1.1d0/                                                       
      DATA swr20  /1.9d0/                                                       
      DATA c2a    /1.5d0/                                                       
      DATA c2b    /1.0d0/                                                       
      DATA c2c    /0.15d0/                                                      
      DATA deh2   /3.39687d0/                                                   
      DATA cfm    /0.0d0/                                                       
C                                                                               
      DATA ccMF  /3.9960D0, 5.4730D-01, -1.1428D0, 2.0440D-01/                  
      DATA ccMH  /1.00d0,   0.80d0,     -2.67d0,   0.456d0/                     
C                                                                               
      DATA ahf  /  3.0000d0/                                                    
      DATA anaf / 18.0000d-1/                                                   
      DATA anah / 16.0000d-1/                                                   
      DATA r21  /  0.7000d0/                                                    
      DATA r23  / -3.3006d-1/                                                   
      DATA r13  /  1.2431d0/                                                    
C                                                                               
      DATA c1s1  / 203.128d0/                                                   
      DATA c2s1  /   1.6167d0/                                                  
      DATA c3s1  /-100.124d0/                                                   
      DATA c4s1  /  68.078d0/                                                   
      DATA c5s1  / -15.832d0/                                                   
      DATA c6s1  /   1.303d0/                                                   
C                                                                               
      DATA dnaht1    /  3.35000d0/                                              
      DATA dnaht2    / 10.000d0/                                                
      DATA dnaht3    / 54.6393d0/                                               
      DATA dnaht4    /-17.63343d0/                                              
      DATA dnaht5    /  3.0085d0/                                               
      DATA dnaht6    /  1.450d0/                                                
      DATA dnaht7    / -2.87977d0/                                              
      DATA dnaht11   /  2.85000d0/                                              
      DATA dnaht61   /  1.300d0/                                                
C                                                                               
      DATA deh222    /  4.7404422d0/                                            
C                                                                               
      DATA una2p     /  2.1037d0/                                               
C                                                                               
      DATA c2a22     /  0.8d0/                                                  
      DATA c2b22     /  0.8d0/                                                  
      DATA c2c22     /  0.15d0/                                                 
C                                                                               
      DATA cprime22  /  8.0d-1/                                                 
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
      DATA swa22     /  0.78495d0/                                              
      DATA swr022    /  7.47654d0/                                              
                                                                                
      DATA a0g1      / -8.210245d-2/                                            
      DATA b0g1      /  4.658698d0/                                             
      DATA r0g1      /  4.552860d0/                                             
      DATA a0g2      /  3.546427d-2/                                            
      DATA b0g2      /  6.456787d-1/                                            
      DATA r0g2      /  6.332225d0/                                             
      DATA a0g3      / -1.562327d-2/                                            
      DATA b0g3      /  1.170537d-1/                                            
      DATA r0g3      /  1.097310d+1/                                            
C                                                                               
      DATA a90g1     / -9.929297d-2/                                            
      DATA b90g1     /  4.367220d0/                                             
      DATA r90g1     /  4.402447d0/                                             
      DATA a90g2     /  4.696181d-2/                                            
      DATA b90g2     /  5.470459d-1/                                            
      DATA r90g2     /  6.289526d0/                                             
      DATA a90g3     / -1.219900d-2/                                            
      DATA b90g3     /  1.002291d-2/                                            
      DATA r90g3     /  1.118725d+1/                                            
                                                                                
      DATA h2t1      /  1.5d0/                                                  
      DATA h2t2      /  4.47009d0/                                              
      DATA h2t3      /  4.48534d0/                                              
      DATA h2t4      /  1.77908d0/                                              
      DATA h2t5      / -6.46276d0/                                              
      DATA h2t6      /  3.15396d0/                                              
      DATA h2t7      / -1.59912d0/                                              
C                                                                               
      END                                                                       
                                                                                
                                                                                
