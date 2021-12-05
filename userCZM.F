*deck,userCZM     parallel       optimize  USERDISTRIB               gal
      subroutine userCZM (matId, elemId, kMatIntPt, ldstep,isubst,
     &                    keycut, ncomp,nProp, nstatev,
     &                    Time, dTime, Temp, dTemp,
     &                    coords, prop, Strain, dStrain, 
     &                    stress, dsdePl, sedEl, sedPl, statev,
     &                    var1, var2, var3, var4, var5)
c*************************************************************************
c     *** primary function ***
c
c           user cohesive zone model example
c
c           Commands
c             TB,CZM,mat,NTEMP,NPTS,user 
c                TBTEMP if mat. constants are temperature dependent
c                TBDATA define material constants
c*************************************************************************
c     input arguments
c     ===============
c      matId     (int,sc,in)              material #
c      elemId    (int,sc,in)              element #
c      kMatIntPt (int,sc,in)              material integration point #
c      ldstep    (int,sc,in)              load step number
c      isubst    (int,sc,in)              substep number
c      ncomp     (int,sc,in)              number of stress, strain components
c      nProp     (int,sc,in)              Number of material constants
c      nstatev   (int,sc,in)              Number of state variables
c
c      Temp      (dp ,sc,in)              temperature at beginning of time increment
c      dTemp     (dp ,sc,in)              temperature increment 
c      Time      (dp ,sc,in)              time at beginning of increment (t)
c      dTime     (dp ,sc,in)              time increment (dt)
c
c      prop     (dp,ar(nprop),i)          Material constants defined by TB command 
c      Strain   (dp,ar(ncomp),i)          Interface separation at beginning of time increment
c      dStrain  (dp,ar(ncomp),i)          Interface separation increment
c      coords   (dp,ar(3),i)              current coordinates
c
c     output arguments              
c     ======================             
c      stress   (dp,ar(nTesn),io)         Traction stress
c      sedEl    (dp,sc,io)                elastic work
c      sedPl    (dp,sc,io)                plastic work
c      keycut   (int,sc,io)               loading bisect/cut control
c                                         0 - no bisect/cut
c                                         1 - bisect/cut 
c                                         (factor will be determined by ANSYS solution control)
c      dsdePl   (dp,ar(ncomp,ncomp),io)   consistent tangent jacobian matrix
c
c     input output arguments              
c     ======================             
c      statev   (dp,ar(nstatev,io)        user defined solution state variables
c
c     misc.
c     ======================             
c      var1, var2, var3, var4, var5       currently not used
c
c     local variables
c     ======================             
c
c      debugflag (in,sc, l)                debugflag to print debug information
c*****************************************************************************************************************************************************************************
c                                                 USER'S INDICATIONS
c MATERIAL CONSTANTS (imput data):
c In TBDATA there are 4+(2*n)+(2*t1)+(2*t2) material constants: 
c  -the first 4 props are: prop(1)=t_a, prop(2)=n-1, prop(3)=t1-1, prop(4)=t2-1
c  -after there are separations (d) and tractions (sig, tau) for each direction ; i.e  t_a,n,t1,t2,d_1,....,d_n,sig_1,....,sig_n,d_1,...d_t1,tau_1,...tau_t1,...etc.
c
c  prop(2) = number points of normal traction-separation law (ModeI) (excluding first point in origin of axes)
c  prop(3) = number points of first tangential t-s law (ModeII) (excluding first point in origin of axes)
c  prop(4) = number points of second tangential t-s law (ModeIII) (excluding first point in origin of axes)
c
c  t_a        (int,sc) = "type analysis" determines the kind of interpolation of point's tractions-separation laws 
c             if t_a =1  Linear interpolation, otherwise Polynomial interpolation.The polynomial interpolation is calculated with Newton's method whith Divided Differences. 
c
c MATERIAL STATE VARIABLES (imput data):
c In TBDATA there are 4+(2*ncomp) material state variables:
c  -the first 4 state variable are internal Ansys state variable (not specified): kchstat=statev(1), kestat0=statev(2), kestat=statev(3), stiffi=statev(4).
c  -after there are user subroutine state variable:
c     if ncomp=2  state(5)=d_old(1), statev(6)=d_old(2), statev(7)=G_last_old(1), statev(8)=G_last_old(2).
c     if ncomp=3  state(5)=d_old(1), statev(6)=d_old(2), statev(7)=d_old(3), statev(8)=G_last_old(1), statev(9)=G_last_old(2), statev(10)=G_last_old(3).
c     Please note; for undamaged CZM material: d_old(i)=0, G_last_old(i)=G0(i) for i=,ncomp
c
c LOCAL USER SUBROUTINE VARIABLES
c  n          (int,sc) = prop(2)+1 number points of normal traction-separation law (ModeI) (including first point in origin of axes)
c  t1         (int,sc) = prop(3)+1 number points of first tangential t-s law (ModeII) (including first point in origin of axes)
c  t2         (int,sc) = prop(4)+1 number points of second tangential t-s law (ModeIII) (including first point in origin of axes)
c 
c  v          (dp,ar(i,j,k)) = array that contains points coordinates of cohesive laws, it has 3 dimensional:
c             where: i=3, it is number of rows: i=1 there are separations of t-s laws; i=2 there are tractions of t-s laws; i=3 there are areas below segments of the t-s laws.
c                    j=m_p+1, it is number of columns refer to max number of points of t-s laws. 
c                    k=ncomp, it is number of sheets: k=1 for normal t-s law; k=2 for first tangential t-s law; k=3 for second tangential t-s law.
c
c  T_inter    (dp,ar,(i,k)) = array that contains cohesive tractions, and local derivatives of cohesive laws when damage process accours.
c             where: i=2, it is number of rows: i=1 there are cohesive tractions, i=2 there are local derivatives of cohesive laws.  
c                    k=ncomp, it is number of columns: k=1 for normal t-s law; k=2 for first tangential t-s law; k=3 for second tangential t-s law. 
c
c  K0         (dp,ar(ncomp)) = array that contains initial penalty stiffness. 
c  d          (dp,ar(ncomp)) = array that contains damage parameters.
c  d_old      (dp,ar(ncomp)) = array that contains max historical damage parameters (state variables).
c  G0         (dp,ar(ncomp)) = array that contains total initial area below t-s laws i.e., critical strain energy release rate: GIc,GIIc,GIIIc.
c  G_last     (dp,ar(ncomp)) = array that contains remaining area under t-s laws curve after damage.
c  G_last_old (dp,ar(ncomp)) = array that contains min historical remaining area below t-s laws (state variables).
c   
c 
c                                                USER PARAMETER SETUP
c  m_p        (int,sc) ="max_points" is max number of interpolating points of t-s laws (excluding first point in origin of axes),thus m_p>=max(prop(2),prop(3),prop(4))
c  i_a        (int,sc) ="integral_accuracy"; greater is this parameter best will be approssimation of G0 and G_last, calculated with Newton-Cotes "composite trapezoidal rule",during polynomial interpolation. 
c  s_t        (logical)= "same t-s law" is a binary parameter used to distinguish the cases the cohesine laws of mode 2 and 3 are same or different: if s_t = true so are same, otherwise no.
c
c Subroutine's author: student Francesco Starace, Master's degree Mechanical Engineering for Design and Production, Federico II University of Naples, academic year 2020-2021 
c Author's contacts: fran.starace@protonmail.com         
c                   
c*****************************************************************************************************************************************************************************
c
#include "impcom.inc"
c
      interface
        function S(V,delta) 
        double precision, dimension(:,:,:), intent(in) :: V
        double precision, dimension(3), intent(in) :: delta
        double precision, dimension(2,3) :: S
        endfunction
      endinterface
      INTEGER          matId,elemId,kMatIntPt,ldstep,isubst,keycut,     
     &                 ncomp,nProp,nstatev,t_a,n,t1,t2,i,j,k,l,z,nco  
c
c --- COMMON
c          
      COMMON /INDICES/ t_a,n,t1,t2,nco 
c    
      DOUBLE PRECISION Time,dTime,Temp,dTemp,sedEl,sedPl,coords(3), 
     &                 prop(nProp), Strain(ncomp), dStrain(ncomp), 
     &                 stress(ncomp), dsdePl(ncomp,ncomp),
     &                 statev(nstatev), var1, var2, var3, var4, var5 
c
c --- PARAMETERS
c
      INTEGER          m_p,i_a
      LOGICAL          s_t
      PARAMETER       (m_p=5,       !"max_points"
     &                 i_a=50,      !"integral_accuracy"
     &                 s_t=.true.)  !"same t-s law" 
c               
      DOUBLE PRECISION, DIMENSION(ncomp) :: k0,d_old,G0,G_last,
     &                                      G_last_old,X
      DOUBLE PRECISION v(3,m_p+1,ncomp),T_inter(2,ncomp),d(3),delta(3)
c
      INTEGER          kchstat,debugflag,kestat,kestat0,kwr
      DOUBLE PRECISION stiffi,dsigCap,kdsigCap 
c
c --- debug includes
#include "locknm.inc"
      INTEGER          wrinqr, iott
      EXTERNAL         pplock, ppunlock, wrinqr
c
c*************************************************************************
c
!     INITIALIZATION

      debugflag = 1
      keycut   = 0
      kwr = 0
      sedPl    = 0.0
      stiffi   = 0.0
      sedEl    = 0.0
      sedPl    = 0.0
      keycut   = 0
      kchstat = nint(statev(1))
      kestat0 = nint(statev(2))
      kestat  = kestat0

      dsdePl=0.0
      nco=ncomp
      k0=0.0
      G0=0.0
      G_last=0.0
      X=0.0
      v=0.0
      T_inter=0.0
      d=0.0
      delta=0.0
!     STATE VARIABLES ASSIGNATION IN ARRAYS: d_old(k),G_last_old(k)
      do i=1,(2*ncomp)
        if(i.le.ncomp)then
          d_old(i)=statev(i+4)
        else
          G_last_old(i-ncomp)=statev(i+4)
        endif
      enddo
!     DATA ALLOCATION IN ARRAYS: delta(k),v(i,j,k),G0(k),K0(k)
      do k=1,ncomp
        if(k.eq.1)then
          delta(k)=0.5*(Strain(k)+abs(Strain(k)))
        else
          delta(k)=abs(Strain(k))
        endif
      enddo
      t_a=nint(prop(1)) 
      n=nint(prop(2))+1
      t1=nint(prop(3))+1
      if(ncomp.eq.3)then
        if(s_t.eq..true.)then
          t2=nint(prop(3))+1
        else  
          t2=nint(prop(4))+1
        endif
      else
        t2=0
      endif
      do k=1,ncomp
        z= 0.5*(k-1)*(k-2)*t2+(k-1)*(3-k)*t1+0.5*(2-k)*(3-k)*n   
        do i=1,2
          do j=2,n
            v(i,j,1)= prop(4+(i-1)*(n-1)+j-1)
          enddo
          do j=2,t1
            v(i,j,2)= prop(4+2*(n-1)+(i-1)*(t1-1)+j-1) 
          enddo
          if(ncomp.eq.3)then
            do j=2,t2 
              if(s_t.eq..true.)then !HP mode 2 ande mode 3 same czm law
                v(i,j,3)=v(i,j,2)
              else                  !HP mode 2 ande mode 3 different czm law           
                v(i,j,3)= prop(4+2*(n-1)+2*(t1-1)+(i-1)*(t2-1)+j-1)
              endif
            enddo
          else
          endif
        enddo 
        IF(t_a.eq.1)then 
          do j=2,z
          v(3,j,k)=0.5*(v(1,j,k)-v(1,j-1,k))*(v(2,j,k)+v(2,j-1,k))
          G0(k)=G0(k)+v(3,j,k)
          enddo
        ELSE
          do j=2,z 
            if(j.eq.2) then
              v(3,j,k)=0.5*v(2,j,k)*v(1,j,k)
            else
              do l=1,i_a-1    
              X(K)=v(1,j-1,k)+((v(1,j,k)-v(1,j-1,k))/i_a)*l
              T_inter=S(V,X)
              v(3,j,k)=v(3,j,k)+(2*T_inter(1,K))
              enddo
              v(3,j,k)=v(3,j,k)+v(2,j-1,k)+v(2,j,k)
              v(3,j,k)=v(3,j,k)*(v(1,j,k)-v(1,j-1,k))/(2*i_a)
            endif  
            G0(k)=G0(k)+v(3,j,k)
          enddo
        ENDIF
        K0(k)=v(2,2,k)/v(1,2,k) 
      enddo
!     CALCULATION OF ARRAYS VALUES: d(k),dsdePl(k,k),G_last(k)  
      T_inter=S(V,delta)
      do k=1,ncomp
        z= 0.5*(k-1)*(k-2)*t2+(k-1)*(3-k)*t1+0.5*(2-k)*(3-k)*n    
        if(delta(k).lt.v(1,2,k))then
          d(k)=0
          G_last(k)=G0(k)
          if(d(k).lt.d_old(k))then
            if((k.eq.1).and.(Strain(1).lt.0))then
              dsdePl(k,k)=K0(k)
            else
              dsdePl(k,k)=K0(k)*(1-d_old(k))
            endif     
          else
          dsdePl(k,k)=K0(k)
          endif
        elseif(delta(k).ge.v(1,z,k))then
          d(k)=1
          G_last(k)=0
          dsdePl(k,k)=0
        else
          d(k)=1-(T_inter(1,k))/(delta(k)*k0(k))
          IF(t_a.eq.1)then          !linear interpolation 
            do j=2,z-1
              if(delta(k).ge.v(1,j,k).and.delta(k).lt.v(1,j+1,k))then 
                G_last(k)=0.5*(T_inter(1,k)*delta(k)+(v(1,j+1,k)-
     &                     delta(k))*(T_inter(1,k)+v(2,j+1,k)))  
              endif
              if(j.le.z-2) G_last(k)=G_last(k)+v(3,j+2,k)  
            enddo
          ELSE                      !polynomial interpolation  
            do j=2,z-1
              if(delta(k).ge.v(1,j,k).and.delta(k).le.v(1,j+1,k))then 
                G_last(k)=0.5*T_inter(1,k)*delta(k)                 
                do l=1,i_a-1    
                  X(k)=delta(k)+((v(1,j+1,k)-delta(k))/i_a)*l
                  T_inter=S(V,X)
                  v(3,1,k)=v(3,1,k)+2*T_inter(1,K)
                enddo
                T_inter=S(V,delta)
                v(3,1,k)=v(3,1,k)+T_inter(1,k)+v(2,j+1,k)
                v(3,1,k)=v(3,1,k)*(v(1,j+1,k)-delta(k))/(2*i_a) 
                G_last(k)=G_last(k)+v(3,1,k)
              endif
              if(j.le.z-2) G_last(k)=G_last(k)+v(3,j+2,k)
            enddo
          ENDIF
          if(d(k).lt.d_old(k))then
            dsdePl(k,k)=K0(k)*(1-d_old(k))
          else 
            dsdePl(k,k)=T_inter(2,k)
          endif 
        endif
        d(k)=max(d(k),d_old(k))
        G_last(k)=min(G_last_old(k),G_last(k))
      enddo
!     CONDITION FOR COMPLETE FAILURE
      if((d(1).eq.1).or.(d(2).eq.1).or.(d(3).eq.1))then
        do k=1,ncomp
          d(k)=1
          G_last(k)=0
          dsdePl(k,k)=0
        enddo
      endif      
!     DETERMINATION OF STATE VARIABLES EACH SUBSTEP: d_old(k),G_last_old(k)
      do k=1,ncomp     
        d_old(k)=d(k)
        G_last_old(k)=G_last(k)
      enddo
      do i=1,(2*ncomp)
        if(i.le.ncomp)then
          statev(4+i)=d_old(i)
        else
          statev(4+i)=G_last_old(i-ncomp)
        endif
      enddo  
!     CALCULATION OF STRESS VALUES
      do k=1,ncomp
        if(k.eq.1)then
          Stress(k)=(1-d(k))*k0(k)*Strain(k)-d(k)*k0(k)*(0.5*
     &              (-Strain(k)+abs(Strain(k))))
        else
          Stress(k)=(1-d(k))*k0(k)*Strain(k)
        endif
      enddo     
c
  900 continue
      statev(2) = kchstat
      statev(3) = kestat 
      statev(4) = stiffi 
      sedEl = 0.d0
      sedPl = 0.d0     
c
  990 continue
  
      if (debugflag .gt. 0) then
         kwr = 1
      end if
c
      if(kwr .gt. 0) then
         call pplock(LOCKOT)
         iott = wrinqr(2)
         write(iott,1000) 'userCZM debug :',
     &   ' elemId=',elemId,
     &   ' kMatIntPt=',kMatIntPt,
     &   ' ldstep=',ldstep,
     &   ' isubst=',isubst,
     &   ' time=',Time,
     &   ' dtime=',dTime,
     &   ' Temp=',Temp,
     &   ' dTemp=',dTemp,
     &   ' coord=',coords
         write(iott,'(a,8e12.5/(11x,8e12.5))') 'prop      :',prop
         write(iott,'(a,6e12.5)') 'stress    :',stress
         write(iott,'(a,6e12.5)') 'Strain    :',Strain
         write(iott,'(a,6e12.5)') 'dsdePl(i,i):',
     &    (dsdePl(i,i),i=1,ncomp)
         write(iott,'(a,6e12.5/(10x,6e12.5))')
     &                 'statev    :',(statev(i),i=1,nstatev)
         call ppunlock(LOCKOT)
 1000    format(a/4(a,i5,1x),4(a,i1,1x)/5x,7(a,e12.5,1x),a,3e12.4)
      end if

      return
      end
!*****************************************************************************************************************************************************************************
! THIS FUNCTION S RETURNS AN ARRAY 3x2 THAT CONTAINS VALUES OF: INTERPOLATING FUNCTION AND DERIVATIVE, OF EACH COHESIVE LAWS IN THE 3 DIRECTIONS. 
      pure function S(V,delta)       
      implicit none
      integer i,j,k,l,z,t_a,n,t1,t2,nco
      COMMON /INDICES/ t_a,n,t1,t2,nco 
      double precision, dimension(:,:,:), intent(in) :: V
      double precision, dimension(3), intent(in) :: delta
      double precision, dimension(2,3) :: S
      double precision, dimension(max(n,t1,t2),max(n,t1,t2),nco) :: A ! Array that contains Divided Differences, used to calculate coefficient of polynomial.
      double precision, dimension((max(n,t1,t2)-1),
     &                           (max(n,t1,t2)-1),nco) :: A1,A2       ! Arrays used to calculate array C.
      double precision, dimension(max(n,t1,t2),nco) :: B,C            ! Arrays used to calculate powers of the terms of the polynomial function and its derivative. 
      S=0.0
      A=0.0
      A1=0.0
      A2=0.0
      B=0.0
      C=0.0
!************************************************LINEAR INPERPOLATION*********************************************************************************************************
      IF(t_a.eq.1)then 
!     CALCULATION OF TRACTIONS  
        do k=1,nco
          z= 0.5*(k-1)*(k-2)*t2+(k-1)*(3-k)*t1+0.5*(2-k)*(3-k)*n 
          do j=2,z-1
            if(delta(k).ge.v(1,j,k).and.delta(k).lt.v(1,j+1,k))then
              S(1,k)=(v(2,j,k)*(v(1,j+1,k)-delta(k))+v(2,j+1,k)*
     &               (delta(k)-v(1,j,k)))/(v(1,j+1,k)-v(1,j,k))
!     CALCULATION OF DERIVATIVES OF T-S LAWS (dsdePl consistent tangent jacobian matrix)
              S(2,k)=(v(2,j+1,k)-v(2,j,k))/(v(1,j+1,k)-v(1,j,k))
            endif 
          enddo  
        enddo
!************************************************POLYNOMIAL INTERPOLATION*****************************************************************************************************
      ELSE
!     CALCULATION OF TRACTIONS 
        do k=1,nco  
          z= 0.5*(k-1)*(k-2)*t2+(k-1)*(3-k)*t1+0.5*(2-k)*(3-k)*n 
          do i=1,z
            do j=i,z
              if(i.eq.1)then
                A(i,j,k)=v(2,j,k)
                B(i,k)=1
              else
                A(i,j,k)=(A(i-1,j-1,k)-A(i-1,j,k))/
     &                   (v(1,j-i+1,k)-v(1,j,k)) 
                B(i,k)=B(i-1,k)*(delta(k)-v(1,i-1,k))  
              endif
            enddo
            S(1,k)=S(1,k)+A(i,i,k)*B(i,k)
          enddo 
!     CALCULATION OF DERIVATIVES OF T-S LAWS (dsdePl consistent tangent jacobian matrix) 
          do i=1,z-1
            do j=1,z-1
              if(i.eq.j)then
                A1(i,j,k)=1
              else
                A1(i,j,k)=delta(k)-v(1,j,k)
              endif
            enddo          
            do j=1,i
              A2(i,j,k)=1  
              do l=1,i
                A2(i,j,k)=A2(i,j,k)*A1(j,l,k)
              enddo
            enddo
          enddo
          do i=1,z
            do j=1,i-1
              if(i.eq.1)then
                C(i,k)=0
              else  
                C(i,k)=C(i,k)+A2(i-1,j,k) 
              endif
            enddo
            S(2,k)=S(2,k)+A(i,i,k)*C(i,k)
          enddo 
        enddo
      ENDIF
      endfunction S
      
      