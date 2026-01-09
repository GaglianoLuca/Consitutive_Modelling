C     Umat coded by William Fuentes and Merita Tafili 
C     Modified Cam Clay, Implicit integration
C     EDITED 07/01/2026 by Luca Gagliano, CohesiveCamClay (Gaume, 2018)
C     Inifinitesimal strains
      SUBROUTINE MCC(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      use tensor_tools
      
      implicit none

      CHARACTER*80 CMNAME
	  
      integer :: NTENS,NSTATEV,NPROPS,NDI,NSHR,NOEL,
     &           NPT,LAYER,KSPT,KSTEP,KINC
	 
      real(8) :: SSE,SPD,SCD,RPL,DRPLDT,DTIME,TEMP,DTEMP,
     &           PNEWDT,CELENT
	 
      real(8) :: STRESS(NTENS),STATEV(NSTATEV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

C ---------------------------------------------------------------     
        real(8) :: T(3,3), DEPS(3,3), JAC(3,3,3,3),
     1 evoid, epsP(3,3), pc, lambda, kappa, M
     2 ,nu, Kbulk, E, G, ELMOD(3,3,3,3), Idev(3,3,3,3), TrialSTRESS(3,3)
     3 , TrialNormS, TrialDevSTRESS(3,3), p , q, Depsp2(3,3)
     4 , TrialNormSt,Trialp,Trialq, TrialF, DeltaPhi, theta
     5 ,TolG, TolF, Gpc, denom, Gpcderiv, Fyield, dfdp
     6 ,dfdq, dfdpc ,dpdDPhi, dqdDPhi, dpcdDPhi ,derivF
     7 ,normaldev(3,3) ,Depsp(3,3), NormDepsp, Depsp1(3,3)
     8 , caa,ca1,ca2,ca3,ca4,ca5,ca6,cbb,cb1,cb2,tsi, NormSt, term1,pcn
        
      real(8) :: a1,a2,a3,a4,a5,a6,a7,a8,a9,Tf1
      real(8) :: t11,t12,t13,tmin,tmax,tdelta,ttest
      real(8) :: b1,b2,b3,cc1,cc2,cc3,beta, Msq,pc11
      real(8) :: p1,g1,g2,g3,g4,g5,pc1,pc2,q1,q2,dpde3(3,3)
      real(8) :: dqde(3,3), c(3,3,3,3), dpde(3,3), df1(3,3),dg1
      real(8) :: da1,db1,dc1,dd1,de1,de2,d1(3,3),d2(3,3),d3(3,3)
      real(8) :: d4(3,3), d5, d6,d7, dDeltaPhide(3,3),dpde2(3,3)
      real(8) :: dqde1, dqde2(3,3),dnde(3,3,3,3),dpde1(3,3)
       
      real(8) :: d31, d32(3,3),d41, d42(3,3),dnde1,dnde2(3,3,3,3)
      real(8) :: dnde3(3,3,3,3),ct1(3,3,3,3),ct2(3,3,3,3),ct3(3,3,3,3)
      real(8) :: a00,a01,a02,a03,a04,a05,a06,a07,a08,a09,a10,a11,
     1 a12,a13,a14,a15,a16,a17,a18, a19,
     2 a20,a21,a22(3,3,3,3),a23(3,3,3,3),a24(3,3,3,3),a25(3,3,3,3),
     3 a26(3,3,3,3),a27(3,3,3,3),a28(3,3,3,3),a29(3,3,3,3),a30(3,3,3,3)
     
  
      integer maxiter, kcounter, j
C     
      CALL Initial(STRESS,T, DSTRAN, DEPS, NTENS,NDI, NSHR)       
C ---------------------------------------------------------------       
C     State variables 
C     Void ratio as state variable (scalar)
      evoid=STATEV(1)
C     Plastic strains as state variable (tensor R3xR3)
      CALL Vectortomatrix(STATEV(2:7),epsP) 
C     Overconsolidated mean stress state variable
      pc=abs(STATEV(8))
C ---------------------------------------------------------------      
C     Defining material parameters 
      lambda=PROPS(1) !compresion index
      kappa=PROPS(2) !swelling index
      M=PROPS(3) !critic state line slope
      nu=PROPS(4) !poisson modulus    
C ---------------------------------------------------------------
C     Mean stress
      p=-Tr(T)/3.0d0   
      if (abs(p)<1.0d-3) then 
      p=1.0d-3
      endif
      if (abs(p)>abs(pc)) then 
      pc=p
      endif
C --------------------------------------------------------------- 
C     ELASTIC MODULUS 
C     Material constants         
      Kbulk=abs((1.0d0+evoid)*p)/(kappa) !bulk modulus
      if (abs(Kbulk)<1.0d0) then 
      Kbulk=1.0d0
      endif 
      if (nu==0.5d0) then
      nu=0.4999999999999999d0
      endif
      E=Kbulk*3.0d0*(1.0d0-2.0d0*nu)!young modulus
      G=E/(2.0d0*(1.0d0+nu))!shear modulus
      Msq=M*M
      
      !***********************************************************************
      !***********************************************************************
      !***********************************************************************
      !***********************************************************************
      !EDIT HERE FOR COHESION PARAMETER
      !***********************************************************************
      !***********************************************************************
      !***********************************************************************
      !***********************************************************************
      beta=0.1d0 !cohesion parameter
      !0 = no cohesion -->Modified Cam Clay
      !0.05 snow slab
      !0.2 weak layer snow 
      CALL ELMOD1(E, nu, ELMOD) 
C ---------------------------------------------------------------  
C ---------------------------------------------------------------  
C     Trial Elastic trial step
C     Deviator fourth unit tensor
      CALL Idev1(Idev) 
C     Trial stress 
      TrialSTRESS=(ELMOD.xx.DEPS)+T
C     Trial deviator stress  
      TrialDevSTRESS=Idev.xx.TrialSTRESS   
C     Norm of trial stress 
      TrialNormS=norm(TrialSTRESS)  
      TrialNormSt= norm(TrialDevSTRESS)    
      call pq(TrialSTRESS,Trialp,Trialq)    
C     Trial yield function F
      Tf1=Trialq*Trialq
      TrialF=Tf1*(1.0d0+2.0d0*beta)+Msq*(Trialp+beta*pc)*(Trialp-pc)
C ---------------------------------------------------------------       
      if (TrialF<0.0d0) then
C     Elastic Steps, Trial step is ok
      T=TrialSTRESS
      JAC=ELMOD
C     Void ratio 
      evoid=evoid+(1.0d0+evoid)*(Tr(Deps))
C --------------------------------------------------------------- 
C     Saving State variables 
C     Alpha isotropic as state variable (scalar)
      STATEV(1)=evoid
C     Plastic strains as state variable (tensor R3xR3)
      CALL Matrixtovector(STATEV(2:7),epsp)
      STATEV(8)=pc
C ---------------------------------------------------------------       
      else ! Plastic corrector step

      DeltaPhi=0.0d0     !initializing consistency parameter
      evoid=evoid+(1.0d0+evoid)*(tr(Deps))
      theta=(1.0d0+evoid)/(lambda-kappa) !state variable (for softening hardening law)
      pcn=pc !preconsolidation pressure at the beginning of the iterations 
       
C ---------------------------------------------------------------
C     Defining tolerances for G and F and maximum iterations
      TolG=1.0E-7*pcn           ! Tolerance for preconsolidation
      TolF=1.0E-7*TrialNormS    ! Tolerance for yield surface
      maxiter=50                ! maximum iterations
C ---------------------------------------------------------------
C     Cycle for iterating with DeltaPhi 
      do kcounter=1,maxiter
C ---------------------------------------------------------------
C     Cycle for iterating with pc   
      do j=0,maxiter  
      
C     Target function
      g1=TrialP+Kbulk*DeltaPhi*Msq*pc*(1.0d0-beta)
      g2=1.0d0+Kbulk*DeltaPhi*2.0d0*Msq
      g3=theta*DeltaPhi*(2.0d0*Msq*(g1/g2)+beta*pc*Msq-Msq*pc)
      Gpc=pcn*exp(g3)-pc
      denom=2.0d0*Kbulk*DeltaPhi !constant 
C     Derivative of G(abs(pc)) with respect to time 
      g4=theta*DeltaPhi*2.0d0*Msq*Kbulk*DeltaPhi*Msq*(1.0d0-beta)
      g5=g4/g2
      Gpcderiv=Gpc*(g5+beta*Msq-Msq)-1.0d0
C     Next value for abs(pc)     
      pc=pc-Gpc/Gpcderiv
      if (abs(Gpc)<TolG) then
      exit
      endif
      enddo
C --------------------------------------------------------------- 
C --------------------------------------------------------------- 
C     Current q 
      q=Trialq-3.0d0*G*DeltaPhi*(2.0d0*q+4.0d0*beta*q)
C     Current p      
      p1=Trialp+Kbulk*DeltaPhi*Msq*pc*(1.d0-beta)
      p=p1/(1.0d0+2.0d0*DeltaPhi*Kbulk*Msq)
C     Current yield function F     
      Fyield=q*q*(1.0d0+2.0d0*beta)+Msq*(p+beta*pc)*(p-pc)      
C     Calculating derivatives of F with respect to p,q,pc
      dfdp=2.0d0*Msq*p+beta*pc*Msq-Msq*pc
      dfdq=2.0d0*q+4.0d0*beta*q
      dfdpc=Msq*(-p+beta*p-2.0d0*beta*pc)
C     Calculating derivatives of p ,q,and pc with respect to DeltaPhi (consistency parameter)
      a1=1+2*Kbulk*DeltaPhi*Msq !blu
      a2=-Kbulk*Msq*(2*p+beta*pc-pc) !a
      a3=theta*DeltaPhi*Msq*(2*p+beta*pc-pc) !b
      a4=Kbulk*DeltaPhi*Msq*pcn*exp(a3)*theta*Msq*DeltaPhi
      a5=1-pcn*exp(a3)*Msq*theta*DeltaPhi*(beta-1.0d0)     
      dpdDPhi=(a2/a1-(a4/(a5*a1)))/(1+a4/(a5*a1))
      q1=-3.0d0*G*(2.0d0*q+4.0d0*beta*q)
      dqdDPhi=q1/(1.0d0+6.d0*G*DeltaPhi*(1.0d0+2.0d0*beta))
      pc11=pcn*exp(a3)*(theta*Msq*(2.0d0*p+beta*pc-pc))
      pc1=pcn*exp(a3)*(theta*Msq*DeltaPhi*2.0d0*dpdDPhi)
      pc2=1.0d0-pcn*exp(a3)*Msq*theta*DeltaPhi*(beta-1.0d0)
      dpcdDPhi=(pc1+pc11)/pc2 
C     Calculating derivative of F with respect to DeltaPhi
      derivF=dfdp*dpdDPhi+dfdq*dqdDPhi+dfdpc*dpcdDPhi   
C     Calculating new DeltaPhi for next iteration
      DeltaPhi=DeltaPhi-Fyield/derivF  
      if (abs(Fyield)<TolF) then
      exit         
      endif  
      enddo 
	  
C     End of iterations 
C ---------------------------------------------------------------
C     Increment of plastic strains DepsP
C ---------------------------------------------------------------   
C     Delta plastic strains
      if  (TrialNormSt/=0.0d0) then
C     Flow rule direction
      normaldev=TrialDevSTRESS/TrialNormSt    
      Depsp1=(1.0d0/3.0d0)*(-dfdp)*delta
      Depsp2=sqrt(3.0d0/2.0d0)*(dfdq)*normaldev
      Depsp=DeltaPhi*(Depsp1+Depsp2)
      NormDepsp=norm(Depsp)
      else
      Depsp=0.0d0
      endif
C ---------------------------------------------------------------
C     Stress
      T=TrialSTRESS-(ELMOD.xx.Depsp)
C     Plastic strain 
      epsP=epsP+Depsp
C     Void ratio 
c      evoid=evoid+(1.0d0+evoid)*Tr(Deps)
      call pq(T,p,q) 
C --------------------------------------------------------------- 
C     Saving State variables 
C     Alpha isotropic as state variable (scalar)
      STATEV(1)=evoid
C     Plastic strains as state variable (tensor R3xR3)
      CALL Matrixtovector(STATEV(2:7),epsp)
      if (p.gt.0.0d0) then
           t11=-beta*Msq
           t12=Msq*beta*p-Msq*p !B
           t13=p*p*Msq+(q*q*(1.0d0+2.0d0*beta)) !C
           tdelta=t12*t12-4.0d0*t11*t13
           tmin=(-t12-sqrt(tdelta))/(2*t11)
           tmax=(-t12+sqrt(tdelta))/(2*t11)  
	    
           if (tmin>tmax) then
           tmax=tmin
           endif
            term1=tmax
           
           if (beta==0) then
                term1=(q/M)**2.0d0/p+p
           endif
           
	  endif
	  pc=term1

	  STATEV(8)=pc
  
     
C ---------------------------------------------------------------
      normaldev=-normaldev
C ---------------------------------------------------------------  
C     Consistent elastoplastic modulus
      if  (NormDepsp/=0) then       
      a00=1.0d0+Kbulk*DeltaPhi*2.0d0*Msq !A
      a01=Kbulk*DeltaPhi*Msq*(beta-1.0d0) !B
      a02=Kbulk*(2.0d0*Msq*p+Msq*beta*pc-Msq*pc) !C
      a03=theta*DeltaPhi*(2.0d0*Msq*p-Msq*pc+Msq*beta*pc) !D
      a04=1.0d0-pcn*exp(a03)*(Msq*beta-Msq)*theta*DeltaPhi !E
      a05=pcn*exp(a03)*theta*DeltaPhi*2.0d0*Msq !F
      a06=pcn*exp(a03)*(2.0d0*Msq*p-Msq*pc+Msq*beta*pc) !H
      a07=1.0d0+(a01*a05)/(a00*a04) !J
      a08=2.0d0*G*sqrt(3.0d0/2.0d0) !N (3x3)
      a09=3.0d0*G*(2.0d0*q*(1.0d0+2.0d0*beta)) !L
      a10=1.0d0+3.0d0*G*DeltaPhi*(2.0d0*(1.0d0+2.0d0*beta)) !M
      a11=2.0d0*q*(1.0d0+2.0d0*beta) !R
      a12=Msq*(2.0d0*p-pc+beta*pc) !S
      a13=Msq*(-p+beta*p-2.0d0*beta*pc) !T
      a14=-a11*a09/a10-a12*(a01*a06+a02*a04)/(a00*a04*a07) ! 1 Q
      a15=a14+(a13*a06*a04*a07-a05*a13*(a01*a06+a02*a04)/(a04*a04*a07))!Q
      a16=-a11*a08/(a10*a15) !3x3 1 deDeltaPhide
      a17=-(a12*Kbulk)/(a00*a07*a15) !3x3 2 deDeltaPhide
      a18=a05*Kbulk*a13/(a04*a00*a07*a15) !3x3 3 deDeltaPhide
      a19=a16+a17+a18 !3x3 deDeltaPhide
      a20=a08/a10-(a09/a10)*a19 !dqde 
      a21=Kbulk/(a00*a07)-(a01*a06+a02*a04)/(a00*a04*a07)*a19 !dpde
 
      dnde1=(2.0d0*G)/TrialNormSt
      dnde2=(1.0d0/3.0d0)*(delta.out.delta)
      dnde3=normaldev.out.normaldev
      dnde=dnde1*(Idelta-dnde2-dnde3)
      a22=(Kbulk/(a00*a07))*(delta.out.delta) !1 C
      a23=-(a01*a06+a02*a04)/(a00*a04*a07)*a16*(delta.out.normaldev) !2 C
      a24=-(a01*a06+a02*a04)/(a00*a04*a07)*a17*(delta.out.delta) !3 C
      a25=-(a01*a06+a02*a04)/(a00*a04*a07)*a18*(delta.out.delta) !4 C
      a26=sqrt(2.0d0/3.0d0)*q*dnde  !5 C
      a27=sqrt(2.0d0/3.0d0)*(a08/a10)*(normaldev.out.normaldev) !6 C
      a28=sqrt(2.0d0/3.0d0)*(-(a09/a10))*a16*(normaldev.out.delta) !7 C
      a29=sqrt(2.0d0/3.0d0)*(-(a09/a10))*a17*(normaldev.out.delta) !8 C
      a30=sqrt(2.0d0/3.0d0)*(-(a09/a10))*a18*(normaldev.out.delta) !9 C
       
      c=a22+a23+a24+a25+a26+a27+a28+a29+a30
      NormSt=q*sqrt(2.0d0/3.0d0)  
      tsi=NormSt/TrialNormSt !flow rule

C     Consistent elastoplastic modulus    
      JAC=c
      
       else
      JAC=ELMOD
       endif  
       endif
C ---------------------------------------------------------------
      Call Solution(NTENS, NDI, NSHR, T, STRESS, JAC
     1 , DDSDDE) 
	   if (isnan(stress(1))) then
	   term1=1.0d0
      endif
      
      END SUBROUTINE MCC
      
c------------------------------------------------------ 
c------------------------------------------------------           
      SUBROUTINE Initial(STRESS,T, DSTRAN, DEPS, NTENS,NDI, NSHR)      
      double precision STRESS(ntens), T(3,3)
     1 ,DSTRAN(ntens), DEPS(3,3) 
      Integer ntens, nshr, ndi   
      DEPS=0.0D0
      T=0.0D0
C     
      do i=1,ndi
      T(i,i)=stress(i)
      DEPS(i,i)=DSTRAN(i)
      enddo 
C     
      if (nshr.ge.1) then
      T(1,2)=stress(4)
      T(2,1)=stress(4)    
      DEPS(1,2)=0.5d0*DSTRAN(4)
      DEPS(2,1)=0.5d0*DSTRAN(4)           
      endif
      if (nshr.ge.2) then
      T(1,3)=stress(5)
      T(3,1)=stress(5)   
      DEPS(1,3)=0.5d0*DSTRAN(5)
      DEPS(3,1)=0.5d0*DSTRAN(5)          
      endif
      if (nshr.ge.3) then
      T(2,3)=stress(6)
      T(3,2)=stress(6)    
      DEPS(2,3)=0.5d0*DSTRAN(6)
      DEPS(3,2)=0.5d0*DSTRAN(6)         
      endif   
      return          
      END SUBROUTINE Initial  
      
c------------------------------------------------------ 
c------------------------------------------------------           
      SUBROUTINE Matrixtovector(vec1,mat1)      
      double precision vec1(6),  mat1(3,3)   
      vec1(1)=mat1(1,1)  
      vec1(2)=mat1(2,2) 
      vec1(3)=mat1(3,3) 
      vec1(4)=mat1(1,2)  
      vec1(5)=mat1(1,3) 
      vec1(6)=mat1(2,3)  
      return          
      END SUBROUTINE Matrixtovector
c------------------------------------------------------ 
c------------------------------------------------------ 
       SUBROUTINE Solution(NTENS, NDI, NSHR, T, STRESS, JAC, DDSDDE) 
       integer NTENS, NDI, NSHR, i, j, k, l
C     Subroutine for filling the stress and Jacobian matrix        
       double precision T(3,3), JAC(3,3,3,3), STRESS(NTENS),
     1 DDSDDE(NTENS,NTENS), JAC66(6,6)
      k=1
      l=1
c------------------------------------------------------
      do i=1,ndi
      stress(i)=T(i,i)
      enddo 
C     
      if (nshr.ge.1) then
      stress(ndi+1)=T(1,2)         
      endif
      if (nshr.ge.2) then
      stress(ndi+2)=T(1,3)         
      endif
      if (nshr.ge.3) then
      stress(ndi+3)=T(2,3)         
      endif   
      call tensortomatrix(jac,  jac66)   
        do i=1,ndi
        do j=1,ndi
          ddsdde(i,j)=jac66(i,j)
        enddo
      enddo 
      do i=ndi+1,ndi+nshr
        do j=1,ndi
          ddsdde(i,j)=jac66(3+k,j)
        enddo
        k=k+1
      enddo  
      do i=1,ndi
      l=1
        do j=ndi+1,ndi+nshr
          ddsdde(i,j)=jac66(i,3+l)
          l=l+1
        enddo
      enddo   
      k=1
      do i=ndi+1,ndi+nshr
        l=1
        do j=ndi+1,ndi+nshr
          ddsdde(i,j)=jac66(3+k,3+l)
          l=l+1
        enddo
        k=k+1
      enddo   
       Return       
      END SUBROUTINE Solution 
c------------------------------------------------------ 
c------------------------------------------------------       
      SUBROUTINE Vectortomatrix(vec1,mat1)      
      double precision vec1(6),  mat1(3,3)   
      mat1(1,1)=vec1(1)  
      mat1(2,2)=vec1(2)  
      mat1(3,3)=vec1(3)  
      mat1(1,2)=vec1(4)  
      mat1(2,1)=vec1(4)  
      mat1(1,3)=vec1(5)  
      mat1(3,1)=vec1(5)  
      mat1(2,3)=vec1(6)  
      mat1(3,2)=vec1(6)  
      return          
      END SUBROUTINE Vectortomatrix             
c------------------------------------------------------ 
c------------------------------------------------------ 
      subroutine tensortomatrix(a3333,  b66)   ! returns b(6,6)
      double precision a3333(3,3,3,3),b66(6,6)
      integer i,j,i9(6),j9(6)
      data i9/1,2,3,1,1,2/
     .     j9/1,2,3,2,3,3/
      do  i=1,6   !  switch to matrix notation
      do  j=1,6
      b66(i,j)=a3333(i9(i),j9(i),i9(j),j9(j))  
      enddo
      enddo   
      return
      end subroutine tensortomatrix
c------------------------------------------------------ 
c------------------------------------------------------ 
      subroutine pq(T,p,q)   ! p and q
      use tensor_tools
      double precision T(3,3),p,q, DevSTRESS(3,3),
     1 NormSt, Idev(3,3,3,3)
C      mean stress 
      p=-1.0d0/3.0d0*Tr(T)
C     deviator stress
      Call idev1(Idev)
      DevSTRESS=dev(T)
C     Norm of the deviator stress 
      NormSt=norm(DevSTRESS)
C     Trial q 
      q=sqrt(3.0d0/2.0d0)*NormSt    
      return
      end subroutine pq
      
c------------------------------------------------------ 
c------------------------------------------------------           
      SUBROUTINE Idev1(Idev) 
      use tensor_tools
      DOUBLE PRECISION Idev(3,3,3,3), 
     1  Ivol(3,3,3,3), Isym(3,3,3,3)
      INTEGER I,J,K,L
      Idev=0.0d0
      Ivol=0.0d0
      Isym=0.0D0
      Do i=1,3
      Do j=1,3
      Do k=1,3
      Do l=1,3
      Isym(i,j,k,l)=1.0d0/2.0d0*(delta(i,k)*delta(j,l)+
     1 delta(i,l)*delta(j,k))              
      Enddo
      Enddo
      Enddo
      Enddo
      Ivol=0.0D0
      Ivol=1.0d0/3.0d0*delta.out.delta    
      Idev=Isym-Ivol         
      END SUBROUTINE Idev1
c------------------------------------------------------ 
c------------------------------------------------------           
      SUBROUTINE ELMOD1(E, nu, ELMOD) 
      use tensor_tools
      DOUBLE PRECISION ELMOD(3,3,3,3), E, nu,
     1 Idev(3,3,3,3), Ivol(3,3,3,3), Isym(3,3,3,3)
      INTEGER I,J,K,L
      Idev=0.0d0
      Ivol=0.0d0
      ELMOD=0.0D0
      Isym=0.0D0
      Do i=1,3
      Do j=1,3
      Do k=1,3
      Do l=1,3
      Isym(i,j,k,l)=1.0d0/2.0d0*(delta(i,k)*delta(j,l)+
     1 delta(i,l)*delta(j,k))              
      Enddo
      Enddo
      Enddo
      Enddo
      Ivol=0.0D0
      Ivol=1.0d0/3.0d0*delta.out.delta    
      Idev=Isym-Ivol         
      ELMOD=E/(1.0D0-2.0D0*nu)*Ivol+E/(1.0d0+nu)*Idev
      END SUBROUTINE ELMOD1  