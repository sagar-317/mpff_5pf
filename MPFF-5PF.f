! MPFF - framework - Sagar (latest edit - 24th June 2025)
! Hybrid multi-phase-field fracture formulation based on 
! Cuntze failure mode concept for progressive damage 
! modelling of elasto-plastic FRP composites
! 
! Base Code:
! Wei Tan, Emilio Martinez-Paneda. 
! Phase field predictions of microscopic fracture and R-curve behaviour of fibre-reinforced composites
! (2021):108539 
! doi: https://doi.org/10.1016/j.compscitech.2020.108539












      module ktransf
      implicit none
C       real*8 UserVar(4,2,400000) !CPE4 or CPE8R 
C       real*8 UserVar(8,2,400000) !C3D8 or C3D20R
      real*8 UserVar(8,10,400000) !C3D8 or C3D20R - 2MPF
      !real*8 UserVar(4,2,400000) !C3D10
      integer nelem,kincK,kkiter,kflagE
      save
      end module

!***********************************************************************
      subroutine uexternaldb(lop,lrestart,time,dtime,kstep,kinc)
      include 'aba_param.inc' !implicit real(a-h o-z)
      dimension time(2)
      if (lop.eq.0) then !start of analysis
       call mutexinit(1)
      endif
      return
      end
      
      subroutine uel(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,
     1 props,nprops,coords,mcrd,nnode,u,du,v,a,jtype,time,dtime,
     2 kstep,kinc,jelem,params,ndload,jdltyp,adlmag,predef,npredf,
     3 lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,njpro,period)

      use ktransf
      include 'aba_param.inc' !implicit real(a-h o-z)

      dimension rhs(mlvarx,*),amatrx(ndofel,ndofel),props(*),svars(*),
     1 energy(*),coords(mcrd,nnode),u(ndofel),du(mlvarx,*),v(ndofel),
     2 a(ndofel),time(2),params(*),jdltyp(mdload,*),adlmag(mdload,*),
     3 ddlmag(mdload,*),predef(2,npredf,nnode),lflags(*),jprops(*)

C       parameter(ndim=2,ninpt=4,nsvint=2,ndof=1) !CPE4 or CPE8R
C       parameter(ndim=3,ninpt=8,nsvint=2,ndof=1) !C3D8 or C3D20R
C       parameter(ndim=3,ninpt=8,nsvint=4,ndof=2) ! 2MPF
      parameter(ndim=3,ninpt=8,nsvint=10,ndof=5) ! 3MPF
      !parameter(ndim=3,ninpt=4,nsvint=2,ndof=1) !C3D10
      
      dimension wght(ninpt),dN(nnode,1),dNdz(ndim,nnode),
     1 dNdx(ndim,nnode),statevLocal(nsvint)

      real*8 
     1 phinodal1(8), phinodal3(8), phinodal4(8), phinodal2(8), 
     2 amatrx_phi1(8, 8), amatrx_phi3(8, 8), amatrx_phi4(8, 8),  
     3 rhs_phi1(8,1), rhs_phi3(8,1), rhs_phi4(8,1), rhs_phi2(8,1), 
     4 rhs_phi5(8,1), phinodal5(8), amatrx_phi2(8, 8), amatrx_phi5(8, 8),
     5 M(3), ply_angle, costheta, sintheta, beta, Struct_Tensor(3,3),
     6 Iden(3,3), XX(8,1), YY1(3), YY2(3), YY3(3), YY4(3), YY5(3),
     7 phinodal1_mat(8,1), phinodal2_mat(8,1), phinodal3_mat(8,1),
     8 MdyadM(3,3), Struct_Tensor_Fibre(3,3),
     9 xr1(8,1), xr2(8,1), xr3(8,1), xr4(8,1), pi, lamina_ele, coh_ele

      integer
     1 rsame1, rsame2, rsame3, rsame4, rsame5, total_elements,
     2 elements_per_ply, LB(16), UB(16), total_plies,
     3 total_elements_overall, lyr_number, glob_num
        
!     initialising
      do k1=1,ndofel
       rhs(k1,1)=0.d0
      end do
      amatrx=0.d0
      rhs_phi1=0.d0; rhs_phi3=0.d0; rhs_phi4=0.d0; 
      rhs_phi2=0.d0; rhs_phi5=0.d0
      amatrx_phi1=0.d0; amatrx_phi3=0.d0; amatrx_phi4=0.d0; 
      amatrx_phi2=0.d0; amatrx_phi5=0.d0
      if (ninpt.eq.4.and.nnode.eq.10.and.ndim.eq.3) then
       wght=0.25d0/6.d0
      else
       wght=1.d0
      endif     


      ! assign ply angles according to element numbering:
      ! assign ply angles according to element numbering:

      pi=4.0d0*atan(1.0d0)
      elements_per_ply = 5157
      total_plies = 8  ! composite plies
      total_elements = elements_per_ply*total_plies  ! pf
      total_elements_overall = elements_per_ply*(2*total_plies-1) ! coh + pf
C       elements_per_ply = total_elements/total_plies


      do i=1,total_plies
        LB(i) = total_elements_overall + 
     +          1 + (i-1)*elements_per_ply
        UB(i) = total_elements_overall + 
     +          i*elements_per_ply
      end do


      if (jelem.ge.LB(1) .and.  jelem.le.UB(1)) then
        ply_angle=    0.d0 ; lyr_number = 1
      end if  
      if (jelem.ge.LB(2) .and.  jelem.le.UB(2)) then
        ply_angle= pi/4.d0 ; lyr_number = 2
      end if
      if (jelem.ge.LB(3) .and.  jelem.le.UB(3)) then
        ply_angle= pi/2.d0 ; lyr_number = 3
      end if 
      if (jelem.ge.LB(4) .and.  jelem.le.UB(4)) then
        ply_angle=-pi/4.d0 ; lyr_number = 4
      end if
      if (jelem.ge.LB(5) .and.  jelem.le.UB(5)) then
        ply_angle=    0.d0 ; lyr_number = 5
      end if 
      if (jelem.ge.LB(6) .and.  jelem.le.UB(6)) then
        ply_angle= pi/4.d0 ; lyr_number = 6
      end if
      if (jelem.ge.LB(7) .and.  jelem.le.UB(7)) then
        ply_angle= pi/2.d0 ; lyr_number = 7
      end if 
      if (jelem.ge.LB(8) .and.  jelem.le.UB(8)) then
        ply_angle=-pi/4.d0 ; lyr_number = 8
      end if


      beta=55.d0
C       ply_angle = props(9)
      costheta = cos(ply_angle)
      sintheta = sin(ply_angle)
      M = (/costheta,sintheta,0.d0/)

      call onem(Iden)
      call dyadic(M, M, 3, MdyadM)

      Struct_Tensor = Iden + beta*MdyadM
      Struct_Tensor_Fibre = Iden + beta*(Iden - MdyadM)
C       Struct_Tensor(3,3) = 0.d0
C       Struct_Tensor_Fibre(3,3) = 0.d0
C       print*, "Struct_Tensor = ", Struct_Tensor
C       print*, "MdyadM = ", MdyadM
C       print*, "****************************** time(2) = ", time(2)
C       print*, "jelem, nelem = ", jelem, nelem
      xelem = jelem - nelem
C       print*, "xelem = ", xelem  
C       print*, "nelem, jelem, ply_angle = ", nelem, jelem, ply_angle
C       print*, "coords = ", coords
C       if (time(2).gt.0.02d0) call XIT
      
      !     reading parameters
      xlc=props(1)
      Gc1=props(2) ! iff1
      Gc2=props(3) ! iff2
      Gc3=props(4) ! iff3
      Gc4=props(5) ! ff1
      Gc5=props(6) ! ff2
      zeta=props(7)*0.d0
      kflag=props(8)
C       xlc_ff1  = 0.418d0 
C       xlc_ff2  = 1.12d0
C       xlc_iff1 = 8.708d0
C       xlc_iff2 = 2.416d0
C       xlc_iff3 = 9.6558d0
      ! minimum value of xlc required = 1.8mm
C       xlc_ff1 = xlc
C       xlc_ff2 = xlc

!     AT1 or AT2 flag (if flag=1, then AT1, else AT2)
      kflagAT=2
      
!     viscous dissipation and iteration counter
      if (kflag.eq.1) then
       if (jelem.eq.1) then
        if (kinc.ne.kincK) then    
         kincK=kinc
         kkiter=1
        else
         kkiter=kkiter+1 
        endif
       endif
      endif 

      do kintk=1,ninpt
!     evaluate shape functions and derivatives
       call kshapefcn(kintk,ninpt,nnode,ndim,dN,dNdz)
       call kjacobian(jelem,ndim,nnode,coords,dNdz,djac,dNdx,mcrd)
       dvol=wght(kintk)*djac

!     recover and assign state variables       
       call kstatevar(kintk,nsvint,svars,statevLocal,1)
       Hn1=statevLocal(1)
       H1=statevLocal(2)
       Hn2=statevLocal(3)
       H2=statevLocal(4)
       Hn3=statevLocal(5)
       H3=statevLocal(6)
       Hn4=statevLocal(7)
       H4=statevLocal(8)
       Hn5=statevLocal(9)
       H5=statevLocal(10)

!     calculate phase field and incremental strains from nodal values
       phi1=0.d0 ;  phi3=0.d0; phi4=0.d0;  phi2=0.d0; phi5=0.d0
       dphi1=0.d0; dphi3=0.d0; dphi4=0.d0; dphi2=0.d0; dphi5=0.d0
       do inod=1,nnode
        phi1=phi1+dN(inod,1)*u(5*inod-4)
        dphi1=dphi1+dN(inod,1)*du(5*inod-4,1)
        phinodal1(inod) = u(5*inod-4)

        phi2=phi2+dN(inod,1)*u(5*inod-3)
        dphi2=dphi2+dN(inod,1)*du(5*inod-3,1)
        phinodal2(inod) = u(5*inod-3)

        phi3=phi3+dN(inod,1)*u(5*inod-2)
        dphi3=dphi3+dN(inod,1)*du(5*inod-2,1)
        phinodal3(inod) = u(5*inod-2)
        
        phi4=phi4+dN(inod,1)*u(5*inod-1)
        dphi4=dphi4+dN(inod,1)*du(5*inod-1,1)
        phinodal4(inod) = u(5*inod-1)        
        
        phi5=phi5+dN(inod,1)*u(5*inod)
        dphi5=dphi5+dN(inod,1)*du(5*inod,1)
        phinodal5(inod) = u(5*inod)
       end do


!      AT1 / AT2 variables
       if (kflagAT.eq.1) then ! AT1
        Hmin=3.d0*Gc/(16.d0*xlc)
        ATpar=1.d0
        c0 = 8.d0/3.d0
       else ! AT2
        Hmin=0.d0
        ATpar=0.d0
        c0= 2.d0
       endif

       rsame1=0; rsame2=0; rsame3=0; rsame4=0; rsame5=0
!     enforcing Karush-Kuhn-Tucker conditions
       H1=max(H1,Hmin)
       if (H1.le.Hn1) then
         H1=Hn1;  ! rsame1 = 1
       end if

       H2=max(H2,Hmin)
       if (H2.le.Hn2) then
         H2=Hn2;  ! rsame2 = 1
       end if

       H3=max(H3,Hmin)
       if (H3.le.Hn3) then
         H3=Hn3;  ! rsame3 = 1
       end if

       H4=max(H4,Hmin)
       if (H4.le.Hn4) then
         H4=Hn4;  ! rsame4 = 1
       end if

       H5=max(H5,Hmin)
       if (H5.le.Hn5) then
         H5=Hn5;  ! rsame5 = 1
       end if
       
!     collect information from UMAT and store state variables
       statevLocal(1)=H1  ! history
       statevLocal(3)=H2  ! history
       statevLocal(5)=H3  ! history
       statevLocal(7)=H4  ! history
       statevLocal(9)=H5  ! history 

       glob_num = jelem - nelem + (lyr_number-1)*elements_per_ply
       statevLocal(2)=UserVar(kintk,2,glob_num)! this is psip from umat
       statevLocal(4)=UserVar(kintk,4,glob_num)! this is psip from umat
       statevLocal(6)=UserVar(kintk,6,glob_num)
       statevLocal(8)=UserVar(kintk,8,glob_num)
       statevLocal(10)=UserVar(kintk,10,glob_num)
        
       call kstatevar(kintk,nsvint,svars,statevLocal,0)
       
C        if (kflag.eq.1) then
C         if (kkiter.le.3) zeta=0.d0
C        elseif (kflag.eq.2) then
C         if (dphi1.ge.0.d0) zeta=0.d0           
C        endif

C        if (kflag.eq.1) then
C         if (kkiter.le.3) zeta=0.d0
C        elseif (kflag.eq.2) then
C         if (dphi3.ge.0.d0) zeta=0.d0           
C        endif
       
       dalph1=2.d0*phi1; ddalph1=2.d0; dGdeg1=-2.d0+2.d0*phi1
       dalph2=2.d0*phi2; ddalph2=2.d0; dGdeg2=-2.d0+2.d0*phi2
       dalph3=2.d0*phi3; ddalph3=2.d0; dGdeg3=-2.d0+2.d0*phi3
       dalph4=2.d0*phi4; ddalph4=2.d0; dGdeg4=-2.d0+2.d0*phi4
       dalph5=2.d0*phi5; ddalph5=2.d0; dGdeg5=-2.d0+2.d0*phi5


       phinodal1_mat(1:8,1) = phinodal1(1:8) ! change to matrix form
       phinodal2_mat(1:8,1) = phinodal2(1:8)
       phinodal3_mat(1:8,1) = phinodal3(1:8)


       YY1 = MATMUL(Struct_Tensor, matmul(dNdx, phinodal1) ) ! for iff1
       YY2 = MATMUL(Struct_Tensor, matmul(dNdx, phinodal2) ) ! for iff2
       YY3 = MATMUL(Struct_Tensor, matmul(dNdx, phinodal3) ) ! for iff3
       YY4 = MATMUL(Struct_Tensor_Fibre, matmul(dNdx, phinodal4) ) ! for ff1
       YY5 = MATMUL(Struct_Tensor_Fibre, matmul(dNdx, phinodal5) ) ! for ff2

       
       rhs_phi1(1:8,1)=rhs_phi1(1:8,1)-dvol*
     1 (matmul(transpose(dNdx),YY1)*2.d0*Gc1*xlc/c0
     2 +dN(:,1)*(dGdeg1*H1+Gc1/xlc*dalph1/c0+zeta*dphi1/dtime))

       amatrx_phi1(1:8,1:8)=amatrx_phi1(1:8,1:8)+
     1 dvol*(matmul(transpose(dNdx),matmul(Struct_Tensor, dNdx))*2.d0*Gc1*xlc/c0+
     2 matmul(dN,transpose(dN))*(Gc1/xlc*ddalph1/c0+zeta/dtime+2.d0*H1))

       rhs_phi2(1:8,1)=rhs_phi2(1:8,1)-dvol*
     1 (matmul(transpose(dNdx),YY2)*2.d0*Gc2*xlc/c0
     2 +dN(:,1)*(dGdeg2*H2+Gc2/xlc*dalph2/c0+zeta*dphi2/dtime))

       amatrx_phi2(1:8,1:8)=amatrx_phi2(1:8,1:8)+
     1 dvol*(matmul(transpose(dNdx),matmul(Struct_Tensor, dNdx))*2.d0*Gc2*xlc/c0+
     2 matmul(dN,transpose(dN))*(Gc2/xlc*ddalph2/c0+zeta/dtime+2.d0*H2))

       rhs_phi3(1:8,1)=rhs_phi3(1:8,1)-dvol*
     1 (matmul(transpose(dNdx),YY3)*2.d0*Gc3*xlc/c0
     2 +dN(:,1)*(dGdeg3*H3+Gc3/xlc*dalph3/c0+zeta*dphi3/dtime))

       amatrx_phi3(1:8,1:8)=amatrx_phi3(1:8,1:8)+
     1 dvol*(matmul(transpose(dNdx),matmul(Struct_Tensor, dNdx))*2.d0*Gc3*xlc/c0+
     2 matmul(dN,transpose(dN))*(Gc3/xlc*ddalph3/c0+zeta/dtime+2.d0*H3))

       rhs_phi4(1:8,1)=rhs_phi4(1:8,1)-dvol*
     1 (matmul(transpose(dNdx),YY4)*2.d0*Gc4*xlc/c0
     2 +dN(:,1)*(dGdeg4*H4+Gc4/xlc*dalph4/c0+zeta*dphi4/dtime))

       amatrx_phi4(1:8,1:8)=amatrx_phi4(1:8,1:8)+
     1 dvol*(matmul(transpose(dNdx),matmul(Struct_Tensor_Fibre, dNdx))*2.d0*Gc4*xlc/c0+
     2 matmul(dN,transpose(dN))*(Gc4/xlc*ddalph4/c0+zeta/dtime+2.d0*H4))

       rhs_phi5(1:8,1)=rhs_phi5(1:8,1)-dvol*
     1 (matmul(transpose(dNdx),YY5)*2.d0*Gc5*xlc/c0
     2 +dN(:,1)*(dGdeg5*H5+Gc5/xlc*dalph5/c0+zeta*dphi5/dtime))

       amatrx_phi5(1:8,1:8)=amatrx_phi5(1:8,1:8)+
     1 dvol*(matmul(transpose(dNdx),matmul(Struct_Tensor_Fibre, dNdx))*2.d0*Gc5*xlc/c0+
     2 matmul(dN,transpose(dN))*(Gc5/xlc*ddalph5/c0+zeta/dtime+2.d0*H5))


      if (rsame1.eq.1) rhs_phi1=0.d0
      if (rsame2.eq.1) rhs_phi2=0.d0
      if (rsame3.eq.1) rhs_phi3=0.d0
      if (rsame4.eq.1) rhs_phi4=0.d0
      if (rsame5.eq.1) rhs_phi5=0.d0


      xx1=-dGdeg4*H4; xx2=-Gc4/xlc*dalph4/c0; xx3=-zeta*dphi4/dtime
      xr1(1:8,1)=-matmul(transpose(dNdx),YY4)*2.d0*Gc4*xlc/c0
      xr2(1:8,1)=-dN(:,1)*(dGdeg4*H4+Gc4/xlc*dalph4/c0+zeta*dphi4/dtime)

C       xr3=
C       xr4=

      if ((jelem-nelem).eq.12510 .and. kintk.eq.1) then
C         print*, " ********** uel ******** "
C         print*, "jelem, kinc, kintk = ", 
C      +           jelem, kinc, kintk
C         print*, "statevLocal(:) = ", statevLocal(:)
C C         print*, "sintheta, costheta = ", sintheta, costheta
C         print*, "phi1,phi2,phi3,phi4,phi5 = ", 
C      +           phi1,phi2,phi3,phi4,phi5
C         print*, "rhs_phi4 = ", rhs_phi4
C         print*, "YY4, xx1, xx2, xx3 = ", YY4, xx1, xx2, xx3
C         print*, "xr1, xr2 = ", xr1, xr2

      end if 

C       Original code isotropic pff:
C       amatrx(1:ndofel,1:ndofel)=amatrx(1:ndofel,1:ndofel)+
C      1 dvol*(matmul(transpose(dNdx),dNdx)*2.d0*Gc*xlc/c0+
C      2 matmul(dN,transpose(dN))*(Gc/xlc*ddalph/c0+zeta/dtime+2.d0*H))

C        rhs(1:ndofel,1)=rhs(1:ndofel,1)-dvol*
C      1 (matmul(transpose(dNdx),matmul(dNdx,u(1:ndofel)))*2.d0*Gc*xlc/c0
C      2 +dN(:,1)*(dGdeg*H+Gc/xlc*dalph/c0+zeta*dphi/dtime))

!     information transfer to UMAT
       
        UserVar(kintk,1,glob_num)=phi1
        UserVar(kintk,3,glob_num)=phi2
        UserVar(kintk,5,glob_num)=phi3
        UserVar(kintk,7,glob_num)=phi4
        UserVar(kintk,9,glob_num)=phi5
       
      end do       ! end loop on material integration points


      ! ASSEMBLY OF AMATRX AND RHS:
      do i=1,8
        rhs(5*i-4,1) = rhs_phi1(i,1)
        rhs(5*i-3,1) = rhs_phi2(i,1)
        rhs(5*i-2,1) = rhs_phi3(i,1)
        rhs(5*i-1,1) = rhs_phi4(i,1)
        rhs(5*i  ,1) = rhs_phi5(i,1)
      enddo


      do i=1,8
        do j=1,8
          amatrx(5*i-4,5*j-4) = amatrx_phi1(i,j)
          amatrx(5*i-3,5*j-3) = amatrx_phi2(i,j)
          amatrx(5*i-2,5*j-2) = amatrx_phi3(i,j)
          amatrx(5*i-1,5*j-1) = amatrx_phi4(i,j)
          amatrx(5*i  ,5*j  ) = amatrx_phi5(i,j) 
        enddo
      enddo     

      RETURN
      END

!***********************************************************************      
      subroutine kshapefcn(kintk,ninpt,nnode,ndim,dN,dNdz)
c
      include 'aba_param.inc'
c
      parameter (gaussCoord=0.577350269d0)
      parameter (gca=0.5854101966d0, gcb=0.1381966012d0)      
      dimension dN(nnode,1),dNdz(ndim,*),coord24(2,4),coord34(3,4),
     * coord38(3,8)
      
      data  coord24 /-1.d0, -1.d0,
     2                1.d0, -1.d0,
     3               -1.d0,  1.d0,
     4                1.d0,  1.d0/  
      
      data  coord34 /gcb, gcb, gcb,
     1               gca, gcb, gcb,
     2               gcb, gca, gcb,
     3               gcb, gcb, gca/       
      
      data  coord38 /-1.d0, -1.d0, -1.d0,
     2                1.d0, -1.d0, -1.d0,
     3               -1.d0,  1.d0, -1.d0,
     4                1.d0,  1.d0, -1.d0, 
     5               -1.d0, -1.d0,  1.d0,
     6                1.d0, -1.d0,  1.d0,
     7               -1.d0,  1.d0,  1.d0,
     8                1.d0,  1.d0,  1.d0/
      
      if (ninpt.eq.4.and.nnode.eq.4.and.ndim.eq.2) then ! CPE4
          
!     determine (g,h)
       g=coord24(1,kintk)*gaussCoord
       h=coord24(2,kintk)*gaussCoord

!     shape functions 
       dN(1,1)=(1.d0-g)*(1.d0-h)/4.d0
       dN(2,1)=(1.d0+g)*(1.d0-h)/4.d0
       dN(3,1)=(1.d0+g)*(1.d0+h)/4.d0
       dN(4,1)=(1.d0-g)*(1.d0+h)/4.d0

!     derivative d(Ni)/d(g)
       dNdz(1,1)=-(1.d0-h)/4.d0
       dNdz(1,2)=(1.d0-h)/4.d0
       dNdz(1,3)=(1.d0+h)/4.d0
       dNdz(1,4)=-(1.d0+h)/4.d0

!     derivative d(Ni)/d(h)
       dNdz(2,1)=-(1.d0-g)/4.d0
       dNdz(2,2)=-(1.d0+g)/4.d0
       dNdz(2,3)=(1.d0+g)/4.d0
       dNdz(2,4)=(1.d0-g)/4.d0
       
      elseif (ninpt.eq.4.and.nnode.eq.8.and.ndim.eq.2) then ! CPE8R 
          
!     determine (g,h,r)
       g=coord24(1,kintk)*gaussCoord
       h=coord24(2,kintk)*gaussCoord

!     shape functions 
       dN(1,1)=-0.25d0*(1.d0-g)*(1.d0-h)*(1.d0+g+h)
       dN(2,1)=0.25d0*(1.d0+g)*(1.d0-h)*(g-h-1.d0)
       dN(3,1)=0.25d0*(1.d0+g)*(1.d0+h)*(g+h-1.d0)
       dN(4,1)=0.25d0*(1.d0-g)*(1.d0+h)*(h-g-1.d0)
       dN(5,1)=0.5d0*(1.d0-g*g)*(1.d0-h)
       dN(6,1)=0.5d0*(1.d0+g)*(1.d0-h*h)
       dN(7,1)=0.5d0*(1.d0-g*g)*(1.d0+h)
       dN(8,1)=0.5d0*(1.d0-g)*(1.d0-h*h)        

!     derivative d(Ni)/d(g)
       dNdz(1,1)=0.25d0*(1.d0-h)*(2.d0*g+h)
       dNdz(1,2)=0.25d0*(1.d0-h)*(2.d0*g-h)
       dNdz(1,3)=0.25d0*(1.d0+h)*(2.d0*g+h)
       dNdz(1,4)=0.25d0*(1.d0+h)*(2.d0*g-h)
       dNdz(1,5)=-g*(1.d0-h)
       dNdz(1,6)=0.5d0*(1.d0-h*h)
       dNdz(1,7)=-g*(1.d0+h)
       dNdz(1,8)=-0.5d0*(1.d0-h*h)      

!     derivative d(Ni)/d(h)
       dNdz(2,1)=0.25d0*(1.d0-g)*(g+2.d0*h)
       dNdz(2,2)=0.25d0*(1.d0+g)*(2.d0*h-g)
       dNdz(2,3)=0.25d0*(1.d0+g)*(2.d0*h+g)
       dNdz(2,4)=0.25d0*(1.d0-g)*(2.d0*h-g)
       dNdz(2,5)=-0.5d0*(1.d0-g*g) 
       dNdz(2,6)=-(1.d0+g)*h 
       dNdz(2,7)=0.5d0*(1.d0-g*g)
       dNdz(2,8)=-(1.d0-g)*h           
          
      elseif (ninpt.eq.4.and.nnode.eq.10.and.ndim.eq.3) then ! C3D10
       f=coord34(1,kintk) ! Multiply here by for consistency
       g=coord34(2,kintk)
       h=coord34(3,kintk)
       
!     shape functions 
       dN(1,1)=(1.d0-f-g-h)*(1.d0-2.d0*f-2.d0*g-2.d0*h)
       dN(2,1)=f*(2.d0*f-1.d0)
       dN(3,1)=g*(2.d0*g-1.d0)
       dN(4,1)=h*(2.d0*h-1.d0)
       dN(5,1)=4.d0*f*(1.d0-f-g-h)
       dN(6,1)=4.d0*f*g
       dN(7,1)=4.d0*g*(1.d0-f-g-h)
       dN(8,1)=4.d0*h*(1.d0-f-g-h)
       dN(9,1)=4.d0*h*f
       dN(10,1)=4.d0*g*h

!     derivative d(Ni)/d(f)       
       dNdz(1,1)=-3.d0+4.d0*f+4.d0*g+4.d0*h
       dNdz(1,2)=4.d0*f-1.d0
       dNdz(1,3)=0.d0
       dNdz(1,4)=0.d0
       dNdz(1,5)=4.d0*(1.d0-f-g-h)-4.d0*f
       dNdz(1,6)=4.d0*g
       dNdz(1,7)=-4.d0*g
       dNdz(1,8)=-4.d0*h
       dNdz(1,9)=4.d0*h
       dNdz(1,10)=0.d0

!     derivative d(Ni)/d(g)
       dNdz(2,1)=-3.d0+4.d0*f+4.d0*g+4.d0*h
       dNdz(2,2)=0.d0
       dNdz(2,3)=4.d0*g-1.d0
       dNdz(2,4)=0.d0
       dNdz(2,5)=-4.d0*f
       dNdz(2,6)=4.d0*f
       dNdz(2,7)=4.d0*(1.d0-f-g-h)-4.d0*g
       dNdz(2,8)=-4.d0*h
       dNdz(2,9)=0.d0
       dNdz(2,10)=4.d0*h
      
!     derivative d(Ni)/d(h)
       dNdz(3,1)=-3.d0+4.d0*f+4.d0*g+4.d0*h
       dNdz(3,2)=0.d0
       dNdz(3,3)=0.d0
       dNdz(3,4)=4.d0*h-1.d0
       dNdz(3,5)=-4.d0*f
       dNdz(3,6)=0.d0
       dNdz(3,7)=-4.d0*g
       dNdz(3,8)=4.d0*(1.d0-f-g-h)-4.d0*h
       dNdz(3,9)=4.d0*f
       dNdz(3,10)=4.d0*g   
          
      elseif (ninpt.eq.8.and.nnode.eq.8.and.ndim.eq.3) then ! C3D8

!     determine (g,h,r)     
       f=coord38(1,kintk)*gaussCoord
       g=coord38(2,kintk)*gaussCoord
       h=coord38(3,kintk)*gaussCoord

!     shape functions
       dN(1,1)=0.125d0*(1.d0-f)*(1.d0-g)*(1.d0-h)
       dN(2,1)=0.125d0*(1.d0+f)*(1.d0-g)*(1.d0-h)
       dN(3,1)=0.125d0*(1.d0+f)*(1.d0+g)*(1.d0-h)
       dN(4,1)=0.125d0*(1.d0-f)*(1.d0+g)*(1.d0-h)
       dN(5,1)=0.125d0*(1.d0-f)*(1.d0-g)*(1.d0+h)
       dN(6,1)=0.125d0*(1.d0+f)*(1.d0-g)*(1.d0+h)
       dN(7,1)=0.125d0*(1.d0+f)*(1.d0+g)*(1.d0+h)
       dN(8,1)=0.125d0*(1.d0-f)*(1.d0+g)*(1.d0+h)     

!     derivative d(Ni)/d(f)
       dNdz(1,1)=-0.125d0*(1.d0-g)*(1.d0-h)
       dNdz(1,2)= 0.125d0*(1.d0-g)*(1.d0-h)
       dNdz(1,3)= 0.125d0*(1.d0+g)*(1.d0-h)
       dNdz(1,4)=-0.125d0*(1.d0+g)*(1.d0-h)
       dNdz(1,5)=-0.125d0*(1.d0-g)*(1.d0+h)
       dNdz(1,6)= 0.125d0*(1.d0-g)*(1.d0+h)
       dNdz(1,7)= 0.125d0*(1.d0+g)*(1.d0+h)
       dNdz(1,8)=-0.125d0*(1.d0+g)*(1.d0+h)

!     derivative d(Ni)/d(g)
       dNdz(2,1)=-0.125d0*(1.d0-f)*(1.d0-h)
       dNdz(2,2)=-0.125d0*(1.d0+f)*(1.d0-h)
       dNdz(2,3)= 0.125d0*(1.d0+f)*(1.d0-h)
       dNdz(2,4)= 0.125d0*(1.d0-f)*(1.d0-h)
       dNdz(2,5)=-0.125d0*(1.d0-f)*(1.d0+h)
       dNdz(2,6)=-0.125d0*(1.d0+f)*(1.d0+h)
       dNdz(2,7)= 0.125d0*(1.d0+f)*(1.d0+h)
       dNdz(2,8)= 0.125d0*(1.d0-f)*(1.d0+h)
      
!     derivative d(Ni)/d(h)
       dNdz(3,1)=-0.125d0*(1.d0-f)*(1.d0-g)
       dNdz(3,2)=-0.125d0*(1.d0+f)*(1.d0-g)
       dNdz(3,3)=-0.125d0*(1.d0+f)*(1.d0+g)
       dNdz(3,4)=-0.125d0*(1.d0-f)*(1.d0+g)
       dNdz(3,5)= 0.125d0*(1.d0-f)*(1.d0-g)
       dNdz(3,6)= 0.125d0*(1.d0+f)*(1.d0-g)
       dNdz(3,7)= 0.125d0*(1.d0+f)*(1.d0+g)
       dNdz(3,8)= 0.125d0*(1.d0-f)*(1.d0+g)
       
      elseif (ninpt.eq.8.and.nnode.eq.20.and.ndim.eq.3) then ! C3D20R       
       
!     determine (g,h,r)  
       f=coord38(1,kintk)*gaussCoord
       g=coord38(2,kintk)*gaussCoord
       h=coord38(3,kintk)*gaussCoord

!     shape functions       
       dN(9,1)= 0.25d0*(1.d0-f**2.d0)*(1.d0-g)*(1.d0-h)
       dN(10,1)=0.25d0*(1.d0+f)*(1.d0-g**2.d0)*(1.d0-h)
       dN(11,1)=0.25d0*(1.d0-f**2.d0)*(1.d0+g)*(1.d0-h)
       dN(12,1)=0.25d0*(1.d0-f)*(1.d0-g**2.d0)*(1.d0-h)
       dN(13,1)=0.25d0*(1.d0-f**2.d0)*(1.d0-g)*(1.d0+h)
       dN(14,1)=0.25d0*(1.d0+f)*(1.d0-g**2.d0)*(1.d0+h)
       dN(15,1)=0.25d0*(1.d0-f**2.d0)*(1.d0+g)*(1.d0+h)
       dN(16,1)=0.25d0*(1.d0-f)*(1.d0-g**2.d0)*(1.d0+h)
       dN(17,1)=0.25d0*(1.d0-f)*(1.d0-g)*(1.d0-h**2.d0)
       dN(18,1)=0.25d0*(1.d0+f)*(1.d0-g)*(1.d0-h**2.d0)
       dN(19,1)=0.25d0*(1.d0+f)*(1.d0+g)*(1.d0-h**2.d0)
       dN(20,1)=0.25d0*(1.d0-f)*(1.d0+g)*(1.d0-h**2.d0)
       dN(1,1)=0.125d0*(1.d0-f)*(1.d0-g)*(1.d0-h)
     & -(dN(9,1)+dN(12,1)+dN(17,1))/2.d0
       dN(2,1)=0.125d0*(1.d0+f)*(1.d0-g)*(1.d0-h)
     & -(dN(9,1)+dN(10,1)+dN(18,1))/2.d0
       dN(3,1)=0.125d0*(1.d0+f)*(1.d0+g)*(1.d0-h)
     & -(dN(10,1)+dN(11,1)+dN(19,1))/2.d0
       dN(4,1)=0.125d0*(1.d0-f)*(1.d0+g)*(1.d0-h)
     & -(dN(11,1)+dN(12,1)+dN(20,1))/2.d0
       dN(5,1)=0.125d0*(1.d0-f)*(1.d0-g)*(1.d0+h)
     & -(dN(13,1)+dN(16,1)+dN(17,1))/2.d0
       dN(6,1)=0.125d0*(1.d0+f)*(1.d0-g)*(1.d0+h)
     & -(dN(13,1)+dN(14,1)+dN(18,1))/2.d0
       dN(7,1)=0.125d0*(1.d0+f)*(1.d0+g)*(1.d0+h)
     & -(dN(14,1)+dN(15,1)+dN(19,1))/2.d0
       dN(8,1)=0.125d0*(1.d0-f)*(1.d0+g)*(1.d0+h)  
     & -(dN(15,1)+dN(16,1)+dN(20,1))/2.d0

!     derivative d(Ni)/d(f)
       dNdz(1,9)=-0.5d0*f*(1.d0-g)*(1.d0-h)
       dNdz(1,10)=0.25d0*(1.d0-g**2.d0)*(1.d0-h)
       dNdz(1,11)=-0.5d0*f*(1.d0+g)*(1.d0-h)
       dNdz(1,12)=-0.25d0*(1.d0-g**2.d0)*(1.d0-h)
       dNdz(1,13)=-0.5d0*f*(1.d0-g)*(1.d0+h)
       dNdz(1,14)=0.25d0*(1.d0-g**2.d0)*(1.d0+h)
       dNdz(1,15)=-0.5d0*f*(1.d0+g)*(1.d0+h)
       dNdz(1,16)=-0.25d0*(1.d0-g**2.d0)*(1.d0+h)
       dNdz(1,17)=-0.25d0*(1.d0-g)*(1.d0-h**2.d0)
       dNdz(1,18)=0.25d0*(1.d0-g)*(1.d0-h**2.d0)
       dNdz(1,19)=0.25d0*(1.d0+g)*(1.d0-h**2.d0)
       dNdz(1,20)=-0.25d0*(1.d0+g)*(1.d0-h**2.d0)
       dNdz(1,1)=-0.125d0*(1.d0-g)*(1.d0-h)
     & -(dNdz(1,9)+dNdz(1,12)+dNdz(1,17))/2.d0
       dNdz(1,2)=0.125d0*(1.d0-g)*(1.d0-h)
     & -(dNdz(1,9)+dNdz(1,10)+dNdz(1,18))/2.d0
       dNdz(1,3)=0.125d0*(1.d0+g)*(1.d0-h)
     & -(dNdz(1,10)+dNdz(1,11)+dNdz(1,19))/2.d0
       dNdz(1,4)=-0.125d0*(1.d0+g)*(1.d0-h)
     & -(dNdz(1,11)+dNdz(1,12)+dNdz(1,20))/2.d0
       dNdz(1,5)=-0.125d0*(1.d0-g)*(1.d0+h)
     & -(dNdz(1,13)+dNdz(1,16)+dNdz(1,17))/2.d0
       dNdz(1,6)=0.125d0*(1.d0-g)*(1.d0+h)
     & -(dNdz(1,13)+dNdz(1,14)+dNdz(1,18))/2.d0
       dNdz(1,7)=0.125d0*(1.d0+g)*(1.d0+h)
     & -(dNdz(1,14)+dNdz(1,15)+dNdz(1,19))/2.d0
       dNdz(1,8)=-0.125d0*(1.d0+g)*(1.d0+h)
     & -(dNdz(1,15)+dNdz(1,16)+dNdz(1,20))/2.d0

!     derivative d(Ni)/d(g)
       dNdz(2,9)=-0.25d0*(1.d0-f**2.d0)*(1.d0-h)
       dNdz(2,10)=-0.5d0*g*(1.d0+f)*(1.d0-h)
       dNdz(2,11)=0.25d0*(1.d0-f**2.d0)*(1.d0-h)
       dNdz(2,12)=-0.5d0*g*(1.d0-f)*(1.d0-h)
       dNdz(2,13)=-0.25d0*(1.d0-f**2.d0)*(1.d0+h)
       dNdz(2,14)=-0.5d0*g*(1.d0+f)*(1.d0+h)
       dNdz(2,15)= 0.25d0*(1.d0-f**2.d0)*(1.d0+h)
       dNdz(2,16)=-0.5d0*g*(1.d0-f)*(1.d0+h)
       dNdz(2,17)=-0.25d0*(1.d0-f)*(1.d0-h**2.d0)
       dNdz(2,18)=-0.25d0*(1.d0+f)*(1.d0-h**2.d0)
       dNdz(2,19)=0.25d0*(1.d0+f)*(1.d0-h**2.d0)
       dNdz(2,20)=0.25d0*(1.d0-f)*(1.d0-h**2.d0)
       dNdz(2,1)=-0.125d0*(1.d0-f)*(1.d0-h)
     & -(dNdz(2,9)+dNdz(2,12)+dNdz(2,17))/2.d0
       dNdz(2,2)=-0.125d0*(1.d0+f)*(1.d0-h)
     & -(dNdz(2,9)+dNdz(2,10)+dNdz(2,18))/2.d0
       dNdz(2,3)=0.125d0*(1.d0+f)*(1.d0-h)
     & -(dNdz(2,10)+dNdz(2,11)+dNdz(2,19))/2.d0
       dNdz(2,4)=0.125d0*(1.d0-f)*(1.d0-h)
     & -(dNdz(2,11)+dNdz(2,12)+dNdz(2,20))/2.d0
       dNdz(2,5)=-0.125d0*(1.d0-f)*(1.d0+h)
     & -(dNdz(2,13)+dNdz(2,16)+dNdz(2,17))/2.d0
       dNdz(2,6)=-0.125d0*(1.d0+f)*(1.d0+h)
     & -(dNdz(2,13)+dNdz(2,14)+dNdz(2,18))/2.d0
       dNdz(2,7)=0.125d0*(1.d0+f)*(1.d0+h)
     & -(dNdz(2,14)+dNdz(2,15)+dNdz(2,19))/2.d0
       dNdz(2,8)=0.125d0*(1.d0-f)*(1.d0+h)
     & -(dNdz(2,15)+dNdz(2,16)+dNdz(2,20))/2.d0
      
!     derivative d(Ni)/d(h)
       dNdz(3,9)=-0.25d0*(1.d0-f**2.d0)*(1.d0-g)
       dNdz(3,10)=-0.25d0*(1.d0+f)*(1.d0-g**2.d0)
       dNdz(3,11)=-0.25d0*(1.d0-f**2.d0)*(1.d0+g)
       dNdz(3,12)=-0.25d0*(1.d0-f)*(1.d0-g**2.d0)
       dNdz(3,13)=0.25d0*(1.d0-f**2.d0)*(1.d0-g)
       dNdz(3,14)=0.25d0*(1.d0+f)*(1.d0-g**2.d0)
       dNdz(3,15)=0.25d0*(1.d0-f**2.d0)*(1.d0+g)
       dNdz(3,16)=0.25d0*(1.d0-f)*(1.d0-g**2.d0)
       dNdz(3,17)=-0.5d0*h*(1.d0-f)*(1.d0-g)
       dNdz(3,18)=-0.5d0*h*(1.d0+f)*(1.d0-g)
       dNdz(3,19)=-0.5d0*h*(1.d0+f)*(1.d0+g)
       dNdz(3,20)=-0.5d0*h*(1.d0-f)*(1.d0+g)
       dNdz(3,1)=-0.125d0*(1.d0-f)*(1.d0-g)
     & -(dNdz(3,9)+dNdz(3,12)+dNdz(3,17))/2.d0
       dNdz(3,2)=-0.125d0*(1.d0+f)*(1.d0-g)
     & -(dNdz(3,9)+dNdz(3,10)+dNdz(3,18))/2.d0
       dNdz(3,3)=-0.125d0*(1.d0+f)*(1.d0+g)
     & -(dNdz(3,10)+dNdz(3,11)+dNdz(3,19))/2.d0
       dNdz(3,4)=-0.125d0*(1.d0-f)*(1.d0+g)
     & -(dNdz(3,11)+dNdz(3,12)+dNdz(3,20))/2.d0
       dNdz(3,5)=0.125d0*(1.d0-f)*(1.d0-g)
     & -(dNdz(3,13)+dNdz(3,16)+dNdz(3,17))/2.d0
       dNdz(3,6)=0.125d0*(1.d0+f)*(1.d0-g)
     & -(dNdz(3,13)+dNdz(3,14)+dNdz(3,18))/2.d0
       dNdz(3,7)=0.125d0*(1.d0+f)*(1.d0+g)
     & -(dNdz(3,14)+dNdz(3,15)+dNdz(3,19))/2.d0
       dNdz(3,8)=0.125d0*(1.d0-f)*(1.d0+g)
     & -(dNdz(3,15)+dNdz(3,16)+dNdz(3,20))/2.d0       
       
      else
       write (6,*) '***ERROR: The shape fuctions cannot be found'   
      endif    
 
      return
      end

!***********************************************************************
      subroutine kjacobian(jelem,ndim,nnode,coords,dNdz,djac,dNdx,mcrd)
!     Notation: djac - Jac determinant; xjaci - inverse of Jac matrix
!     dNdx - shape functions derivatives w.r.t. global coordinates      
      include 'aba_param.inc'

      dimension xjac(ndim,ndim),xjaci(ndim,ndim),coords(mcrd,nnode),
     1 dNdz(ndim,nnode),dNdx(ndim,nnode)      

      xjac=0.d0

      do inod=1,nnode
       do idim=1,ndim
        do jdim=1,ndim
         xjac(jdim,idim)=xjac(jdim,idim)+
     1        dNdz(jdim,inod)*coords(idim,inod)      
        end do
       end do 
      end do

      if (ndim.eq.3) then
          
       djac=xjac(1,1)*xjac(2,2)*xjac(3,3)+xjac(2,1)*xjac(3,2)*xjac(1,3)
     & +xjac(3,1)*xjac(2,3)*xjac(1,2)-xjac(3,1)*xjac(2,2)*xjac(1,3)
     & -xjac(2,1)*xjac(1,2)*xjac(3,3)-xjac(1,1)*xjac(2,3)*xjac(3,2)
       if (djac.gt.0.d0) then ! jacobian is positive - o.k.
        xjaci(1,1)=(xjac(2,2)*xjac(3,3)-xjac(2,3)*xjac(3,2))/djac
        xjaci(1,2)=(xjac(1,3)*xjac(3,2)-xjac(1,2)*xjac(3,3))/djac
        xjaci(1,3)=(xjac(1,2)*xjac(2,3)-xjac(1,3)*xjac(2,2))/djac
        xjaci(2,1)=(xjac(2,3)*xjac(3,1)-xjac(2,1)*xjac(3,3))/djac
        xjaci(2,2)=(xjac(1,1)*xjac(3,3)-xjac(1,3)*xjac(3,1))/djac
        xjaci(2,3)=(xjac(1,3)*xjac(2,1)-xjac(1,1)*xjac(2,3))/djac
        xjaci(3,1)=(xjac(2,1)*xjac(3,2)-xjac(2,2)*xjac(3,1))/djac
        xjaci(3,2)=(xjac(1,2)*xjac(3,1)-xjac(1,1)*xjac(3,2))/djac
        xjaci(3,3)=(xjac(1,1)*xjac(2,2)-xjac(1,2)*xjac(2,1))/djac
       else ! negative or zero jacobian
        write(7,*)'WARNING: element',jelem,'has neg. Jacobian'
       endif 
          
      else if (ndim.eq.2) then 
          
       djac=xjac(1,1)*xjac(2,2)-xjac(1,2)*xjac(2,1)
       if (djac.gt.0.d0) then ! jacobian is positive - o.k.
        xjaci(1,1)=xjac(2,2)/djac
        xjaci(2,2)=xjac(1,1)/djac
        xjaci(1,2)=-xjac(1,2)/djac
        xjaci(2,1)=-xjac(2,1)/djac
       else ! negative or zero jacobian
        write(7,*)'WARNING: element',jelem,'has neg. Jacobian'
       endif
       
      endif
      
      dNdx=matmul(xjaci,dNdz)

      return
      end

!***********************************************************************
      subroutine kstatevar(npt,nsvint,statev,statev_ip,icopy)
c
c     Transfer data to/from element-level state variable array from/to
c     material-point level state variable array.
c
      include 'aba_param.inc'

      dimension statev(*),statev_ip(*)

      isvinc=(npt-1)*nsvint     ! integration point increment

      if (icopy.eq.1) then ! Prepare arrays for entry into umat
       do i=1,nsvint
        statev_ip(i)=statev(i+isvinc)
       enddo
      else ! Update element state variables upon return from umat
       do i=1,nsvint
        statev(i+isvinc)=statev_ip(i)
       enddo
      end if

      return
      end

!***********************************************************************
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt,
     1 drplde,drpldt,stran,dstran,time,dtime,temp2,dtemp,predef,dpred,
     2 cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     3 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)

      use ktransf
      include 'aba_param.inc' !implicit real(a-h o-z)

      character*8 cmname

      PARAMETER (ZERO= 0.D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,SIX=6.0D0,
     + NEWTON=10,TOLER=1.D-6, NINPT = 8,HALF=0.5D0
     1, FOURTH=.25D0)

      PARAMETER (YSTRESS_MIN = 60.D0, MAXITER = 300, TOL = 1.D-3, 
     +       PHI_MAX = 0.99D0, VISCO = 0.05D0)
      dimension stress(ntens),statev(nstatv),ddsdde(ntens,ntens),
     1 ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens),
     2 time(2),predef(1),dpred(1),props(nprops),coords(3),drot(3,3),
     3 dfgrd0(3,3),dfgrd1(3,3),jstep(4)
      
      dimension Edev(ntens),eprin(3)

      REAL*8 EMOD,ENU,EBULK3,EBULK,EG,EG2,EG3,ELAM,TIMEZ  

      REAL*8 P_STR(6),P(6,6),Z(6,6), ESTRESS(6),M(6,6),M_INV(6,6),
     4 IDENT(6,6),DS_DDLAMBDA(6),F_DL_T1_1(6,6),F_DL_T1(6,6),PSTRAIN(6),
     5 F_DL_T2(6),DS_DDL_T1(6,6),DEQS_DDL_T2(6,6),NF(6),DEQS_DDL_T3(6),
     6 MAT(6,6),TV1(6),V2(6),C2(6),AIN(6,6),C(6,6),DDS_T3_1(6,6),T_INV(6,6),
     7 ESTRAIN(6),TSTRAIN(6),TSTRESS(6),CF(6,6),SE(6,6), CFDEG(6,6),
     8 LSTRESS(6),LCF(6,6),LDSTRAIN(6),LSTRAIN(6),T(6,6),LSTRESSOLD(6)

      REAL*8    
     1 E1,E2,E3,G12,G13,G23,V12,V13,V23,V21,V31,V32,IFF1,IFF2,IFF3,
     2 XT,XC,YT,YC,VF12,EF1,EQPSTRAIN, DLAMBDA,YSTRESS,F_DLAMBDA,
     3 S21,FFT,FFC,MFT,MFC,NMP,NW1,G23O,DYS_DDLAMBDA_1,TXI,TH,
     4 MFLC,DFT,DFC,DMT,DMC,DMLC,STRANT(6),CD(6,6),
     5 T21C,RVVA,ZERO,ONE,MAXIM,E2O,G12O,LC,THETA,
     6 CFULL(6,6),MFMV,DMG,COUNTER,XX,VAR,E1O,XTF,A,
     7 DEQS_DDLAMBDA_N,DEQS_DDLAMBDA_D,DEQS_DDLAMBDA,DF_DDLAMBDA,C1,
     8 DYS_DDLAMBDA,DLAMBDA_N,F_DL(1,1),F_DLAMDA,DEQS_DDL_N,PXI,FAILSURF,FS1

      REAL*8
     1 CF_INV(6,6), H_TILDE, NF_MAT(6,1), NF_MAT_TRANS(1,6),
     2 XI_TILDE, XI_TILDE_SQUARE(1,1), NF_NFTRANS(6,6), D_EP_INV(6,6) ,
     3 DPSTRAIN(6), ESTRESS_MAT(6,1), Hc1, Hc2, Hc3, Hc4, Hc5,
     4 hel1, hel2, hel3, hel4, hel5
      REAL*8 
     1 I_PLUS_DLAMBDA_DE_P(6,6), I_PLUS_DLAMBDA_DE_P_INV(6,6), 
     2 I_PLUS_DLAMBDA_DE_P_INV_SIG(6,1), ESTRESS_NPLUS1(6,1), ESTRESS_NPLUS1_TR(1,6),
     3 I_PLUS_DLAMBDA_DE_P_INV_SIG_TR(1,6), DSIG_BY_DDLAMBDA(6,1),
     4 DF_BY_DDLAMBDA_T1(1,1), DF_BY_DDLAMBDA_T1_SCALAR, P_SIG_NPLUS1(6,1),
     5 P_SIG_NPLUS1_TR(1,6), DP_BY_DDLAMBDA_T4(6,1), DP_BY_DDLAMBDA(1,1),
     6 DSIGY_BY_DDLAMBDA, PSIE_PLUS_F, PSIE_PLUS_M, PSIE_PLUS_S, INC_PCD_M, 
     7 INC_PCD_S, PHI_F, PHI_M, PHI_S, PCD_S, PCD_M, H_PCD_M, H_PCD_S,
     8 H_PLUS_F, H_PLUS_M, H_PLUS_S, YSTRESS_MIN, DEELAS(1:6), CF_NFMAT(6,1),
     9 X_B, X_A, NFMAT_CF_NFMAT, GD_F, GD_M, GD_S, G1, G2, PHI1, PHI2, PHI3

      REAL*8
     1 PHI_FT, PHI_FC, PHI_MTT, PHI_MTC, PHI_MLS, C_FF(6,6), C_IFF(6,6),
     2 GD_FT, GD_FC, GD_MTT, GD_MTC, GD_MLS, GD_FF, GD_IFF,
     3 ALPHA_FT, ALPHA_FC, ALPHA_MTT, ALPHA_MTC, ALPHA_MLS,
     4 PSIE_PLUS_FT, PSIE_PLUS_FC, PSIE_PLUS_MTT, PSIE_PLUS_MTC, PSIE_PLUS_MLS,
     5 HFTPLUS, HFCPLUS, G_M1, G_M2, G_M3,
     8 DMINV_BY_DPHI3(6,6), DMINV_BY_DPHI2(6,6), DMINV_BY_DPHI1(6,6),
     9 DMINV0_BY_DPHI3(6,6), DMINV0_BY_DPHI2(6,6), DMINV0_BY_DPHI1(6,6),
     8 w_ff1, w_ff2, w_iff1, w_iff2, w_iff3, psi_iff1, psi_iff2, psi_iff3,
     9 psi_ff1, psi_ff2, eelas(ntens), eplas(ntens)

      REAL*8
     1 DMINV0_BY_DPHIFC(6,6), DMINV0_BY_DPHIFT(6,6), DMINV_BY_DPHIFC(6,6), DMINV_BY_DPHIFT(6,6),
     2 COMPL(6,6), H1,H2,H3,H4,H5, H1PLUS, H2PLUS, H3PLUS,
     3 IFF1_COUNTER, IFF2_COUNTER, IFF3_COUNTER, H1_C, H2_C, H3_C,
     4 GC_S, GC_C, GC_T, MDOT, B1, B2, I_23_5, PHI_MAX, 
     5 INC_PCD_1, INC_PCD_2, INC_PCD_3, XSIG1(6), XSIG2(6), XSIG3(6),
     6 FFT_COUNTER, FFC_COUNTER, HC_FC, HC_FT, HFC(1,1), HFT(1,1),
     7 GC_FT, GC_FC, INC_PCD_FT, INC_PCD_FC, HPLUSFC, HPLUSFT,
     8 XSIGFT(6), XSIGFC(6), PCD1, PCD2, PCD3, PCD_FC, PCD_FT,
     9 H1V, H2V, H3V, HFTV, HFCV, VISCO, FE_EFF, gd1, gd2, gd3, gd4, gd5

      real*8
     1 IFF1n, IFF2n, IFF3n, FF1n, FF2n, iff1_hist, iff2_hist, iff3_hist,
     2 ff1_hist, ff2_hist, FE_EFFn, FE_EFF_hist, DELTA, xeelas11, 
     3 gd_ff_cap, gd_iff_cap, cff(6,6), ciff3(6,6), ciff(6,6),
     4 estress_ff(6), estress_iff(6), estress_iff3(6), gs_cap,
     5 ff1ind,ff2ind,iff1ind,iff2ind,iff3ind,fmax, PsifByPsif0, PsimByPsim0,
     3 XGDF, xpsi_m, xpsi_s
     4
     5 

      INTEGER MAXITER, ITER, VISCO_CHOICE, ISCOUPLED
     1 dcff1,dcff2,dciff1,dciff2,dciff3

      
!     find number of elements 
      if (npt.eq.1) then
       if ((time(1)-dtime).lt.0.d0) then
        if (kflagE.ne.13) then
         kincK=0 
         kkiter=1 
         nelem=noel
         UserVar=0.d0
         kflagE=13
        else
         CALL MutexLock(1)   
         if (noel.gt.nelem) nelem=noel
         CALL MutexUnlock(1)
        endif 
       endif
      endif


!     Initialization
      ddsdde=0.d0
C       E=props(1) ! Young's modulus
C       xnu=props(2) ! Poisson's ratio
      xk=1.d-07


!---------------------READ PROPERTIES--------------------------------------
      E1 = PROPS(1)           !YOUNG'S MODULUS IN DIRECTION 1 (L)
      E2 = PROPS(2)           !YOUNG'S MODULUS IN DIRECTION 2 (T) 
      E3=E2               
      G12 = PROPS(3)          !SHEAR MODULUS IN 12 PLANE
      G13=G12                 !SHEAR MODULUS IN 13 PLANE
      V12=PROPS(4)            !POISSON RATIO IN 12
      V23=PROPS(5)            !POISSON RATIO IN 23
      V13=V12                 !POISSON RATIO IN 13 
      EF1=PROPS(6)            !MODULUS OF FIBER PARALLEL TO FIBER
      VF12 = PROPS(7)         !POISSON RATIO OF FIBER
      XT = PROPS(8)           !TENSILE STRENGTH PARALLEL TO FIBER
      XC = PROPS(9)           !COMPRESSIVE STRENGTH PARALLEL TO FIBER
      YT = PROPS(10)          !TENSILE STRENGTH PERPENDICULAR TO FIBER
      YC = PROPS(11)          !COMPRESSIVE STRENGTH PERPENDICULAR TO FIBER
      S21 = PROPS(12)         !IN PLANE SHEAR STRENGTH
      MAT = PROPS(13)         !MATERIAL TYPE FOR INCLINATION PARAMETERS
      GC_T = PROPS(14)        !IFF1 GC
      GC_C = PROPS(15)        !IFF2 GC
      GC_S = PROPS(16)        !IFF3 GC
      GC_FT = PROPS(17)       !IFF2 GC
      GC_FC = PROPS(18)       !IFF3 GC
      ANGLE = PROPS(19)* ATAN(1.0) / 45.D0
      ETA = 0.0002D0
      G23 = E2/(2.D0*(1.D0+V23))     !SHEAR MODULUS IN 23 PLANE
      XTF=XT*EF1/E1           !EFFECTIVE TENSILE STRENGTH OF FIBER
      E11=E1; E22=E2; E33=E3
C
      V21=(E2/E1)*V12
      V31=(E3/E1)*V13
      V32=(E3/E2)*V23

      ! DEFINE CONSTANTS
      A = 1.5d0  ! as4 peek from u.din
      ALPHA = 0.21D0 ! AS4/PEEK changing to calibrate
      BETA = 280.D0  ! AS4/PEEK 
      DLAMBDA = 1.0D-20  ! INNITIALIZE MULTIPLIER AS VERY SMALL


!------------ INITIALIZATION OF MATRICES-------------------------------- 
      DO K1=1,6
          DO K2=1,6
              CF(K1,K2)=0.D0
              C_FF(K1,K2)=0.D0
              C_IFF(K1,K2)=0.D0
              P(K1,K2)=0.D0
              Z(K1,K2)=0.D0 
              IDENT(K1,K2) = 0.0D0     
          ENDDO
          IDENT(K1,K1) = 1.0D0
      ENDDO 
      M = IDENT
      M_INV = IDENT    
      DO I=1,NTENS  ! INITIALIZATION OF JACOBIAN MATRIX
          DO J=1,NTENS
              DDSDDE(I,J) = ZERO
          ENDDO
      ENDDO  


      ! UMAT CONVENTIONS TO CONSIDER - 
      ! UMAT TAKES ENGINEERING STRAINS WHILE VUMAT TAKES TENSORIAL
      ! UMAT HAS ORDER DIFFERENT FROM VUMAT, WHICH IS:
      ! UMAT - 11 22 33 12 13 23
      ! VUMAT- 11 22 33 12 23 13


      DELTA=1.D0/(1.D0-V12*V21-V23*V32-V13*V31-2.D0*V21*V32*V13)
      CF(1,1)=E11*DELTA*(ONE-V23*V32)
      CF(1,2)=E11*DELTA*(V21+V31*V23)
      CF(1,3)=E11*DELTA*(V31+V21*V32)
      CF(2,1)=E11*DELTA*(V21+V31*V23)
      CF(2,2)=E22*DELTA*(ONE-V13*V31)
      CF(2,3)=E22*DELTA*(V32+V12*V31)
      CF(3,1)=E11*DELTA*(V31+V21*V32)
      CF(3,2)=E22*DELTA*(V32+V12*V31)
      CF(3,3)=E33*DELTA*(ONE-V12*V21)
      CF(4,4)=G12
      CF(5,5)=G13
      CF(6,6)=G23
      DDSDDE(1:6,1:6) = CF(1:6,1:6)


      cff=0.d0; ciff=0.d0; ciff3=0.d0  ! all 6x6 matrices in umat notation
      cff  (1,1) = CF(1,1)
      ciff3(4,4) = CF(4,4)
      ciff3(5,5) = CF(5,5)
      ciff       = CF
      ciff(1,1)  = 0.d0
      ciff(4,4)  = 0.d0
      ciff(5,5)  = 0.d0


      ! DEFINE THE COMPLIANCE MATRIX, COMPL
      COMPL = 0.D0 ! INITIALIZE AS ZEROS
      COMPL(1,1) = 1.D0/E1
      COMPL(1,2) = -V12/E1
      COMPL(1,3) = -V13/E1
      COMPL(2,1) = COMPL(1,2)
      COMPL(2,2) = 1.D0/E2     
      COMPL(2,3) = -V23/E2      
      COMPL(3,1) = COMPL(1,3)    
      COMPL(3,2) = COMPL(2,3)    
      COMPL(3,3) = 1.D0/E3   
      COMPL(4,4) = 1.D0/G12   
      COMPL(5,5) = 1.D0/G13
      COMPL(6,6) = 1.D0/G23


      ! PROJECTION TESNOR
      P(2,2) = 3.0D0
      P(3,3) = 3.0D0
      P(2,3) = -3.0D0
      P(3,2) = -3.0D0
      P(4,4) = 6.0D0 * A
      P(5,5) = 6.0D0 * A
      P(6,6) = 12.0D0

      ! Z TENSOR FOR PLASTC RETURN MAPPING
      Z(2,2) = 2.0D0/3.0D0
      Z(3,3) = 2.0D0/3.0D0
      Z(2,3) = -1.0D0/6.0D0
      Z(3,2) = -1.0D0/6.0D0
      Z(4,4) = 1.0D0/(3.0D0 * A)
      Z(5,5) = 1.0D0/(3.0D0 * A)
      Z(6,6) = 1.0D0/6.0D0


!     Amor (kflagD=1) or Miehe [only 2D] (kflagD=2) decomposition flag      
      kflagD=0
      

      ! GET VARIABLES FROM STATEVARS
      phi1=statev(ntens+1)
      phi2=statev(ntens+3)
      phi3=statev(ntens+5)
      phi4=statev(ntens+7)
      phi5=statev(ntens+9) 
      dam_counter = statev(17)  ! damage counter. its initial value is 0
      ! previous value of loading functions


      dcff1=statev(65)
      dcff2=statev(66)
      dciff1=statev(67)
      dciff2=statev(68)
      dciff3=statev(69)
      IFF1n = statev(25)
      IFF2n = statev(26)
      IFF3n = statev(27)
      FF1n  = statev(28)
      FF2n  = statev(29)
      estress=statev(1:6)
      eelas = statev(30:35)
      eplas = statev(36:41)
      pcd1  = statev(42)
      pcd2  = statev(43)
      pcd3  = statev(44)
      pcd4  = statev(45)
      pcd5  = statev(46)
      psi_iff1 = statev(47)
      psi_iff2 = statev(48)
      psi_iff3 = statev(49)
      psi_ff1 = statev(50)
      psi_ff2 = statev(51)
      EQPSTRAIN=statev(53)
      FE_EFFn = statev(23) ! previous value of feff


      ! degradatrion functions defined below
      gd1 = (1.D0 - (phi1))**2
      gd2 = (1.D0 - (phi2))**2
      gd3 = (1.D0 - (phi3))**2
      gd4 = (1.D0 - (phi4))**2
      gd5 = (1.D0 - (phi5))**2  


      ! calculate weights of the modes based on previous estress
C       if (FE_EFFn.ge.1.D0) then
C         w_ff1 = FF1n/(FF1n+FF2n) 
C         w_ff2 = FF2n/(FF1n+FF2n)
C         w_iff1 = IFF1n/(IFF1n+IFF2n+IFF3n) 
C         w_iff2 = IFF2n/(IFF1n+IFF2n+IFF3n)
C         w_iff3 = IFF3n/(IFF1n+IFF2n+IFF3n)
C         M_INV(1,1) = w_ff1*(gd4+xk) + w_ff2*(gd5+xk)
C         M_INV(2,2) = w_iff1*(gd1+xk) + w_iff2*(gd2+xk) + w_iff3*(gd3+xk)
C         M_INV(3,3) = w_iff1*(gd1+xk) + w_iff2*(gd2+xk) + w_iff3*(gd3+xk)
C         M_INV(4,4) = (gd3+xk)
C         M_INV(5,5) = (gd3+xk)
C         M_INV(6,6) = w_iff1*(gd1+xk) + w_iff2*(gd2+xk) + w_iff3*(gd3+xk)
C         do i=1,6
C           M(i,i) = 1.d0/M_INV(i,i)
C         end do              
C       end if  
      ! if damage not activated, then M and M_inv will be identity matrices

      ! ISCOUPLED = 1 FOR COUPLING
      ! ISCOUPLED = 0 FOR DECOUPLING
C       ISCOUPLED = 1
      

      stran=stran+dstran   !stran represents the total strain 
                           !at the beginning of increment
C       stress=matmul(ddsdde,stran)
C       estress = matmul(ddsdde,stran)   ! effective stress - store as sdv
      
      !!!!!!!!!!!!!!!!!BEGINNING OF PLASTICITY MODELLING!!!!!!!!!!!!!!!!


      !Calculate predictor stress and elastic strain

      EELAS = EELAS + DSTRAN                ! TRIAL ELASTIC STRAIN 


      DO I=1,NTENS
        DO J=1,NTENS
          ESTRESS(I)=ESTRESS(I)+DDSDDE(I,J)*DSTRAN(J)
        END DO
      END DO 


      YSTRESS = MAX(YSTRESS_MIN, BETA * (EQPSTRAIN+1.D-10)**ALPHA+YSTRESS_MIN) ! CALCULATE YIELD STRESS   
      

      CALL TV_MAT_V(ESTRESS,P,ESTRESS,FS1)
      ! [sigma(6x1)]^T[P(6x6)][sigma(6x1)] = [1x1] = Yield criteria's 1st term
      FAILSURF = 0.5D0 * FS1 - YSTRESS**2 ! yield criteria
      DO I = 1,6   ! CONVERT EFFECTIVE STRESS VECTOR INTO A MATRIX
        ESTRESS_MAT(I,1) = ESTRESS(I)
      END DO
      DPSTRAIN = 0.D0 ! INITIALIZE AS 0.D0

      
      IF (FAILSURF .LT. ZERO) THEN
        DPSTRAIN = 0.D0
        DEQPSTRAIN=0.D0

      ELSE
        ITER1 = ONE; MARK = ONE
        NF_MAT = MATMUL(P,ESTRESS_MAT)
        
        DO I = 1,6
            NF_MAT_TRANS(1,I) = NF_MAT(I,1)
        END DO
        
        XI_TILDE_SQUARE = MATMUL(NF_MAT_TRANS, MATMUL(Z, NF_MAT))
        XI_TILDE = SQRT(XI_TILDE_SQUARE(1,1))
        ESTRESS_NPLUS1 = ESTRESS_MAT         ! UPDATED STRESS INSIDE RMA

        DO WHILE (MARK .EQ. ONE)
          
          ! CALCULATE EFFECTIVE SIGMA_N+1 FROM UPDATED DLAMBDA VALUE   
          I_PLUS_DLAMBDA_DE_P = IDENT + DLAMBDA*MATMUL(CF, P)
          
          CALL INVERSE(I_PLUS_DLAMBDA_DE_P, I_PLUS_DLAMBDA_DE_P_INV, 6)
          
          I_PLUS_DLAMBDA_DE_P_INV_SIG = MATMUL(I_PLUS_DLAMBDA_DE_P_INV, 
     +                                         ESTRESS_MAT)
          ! ESTRESS_MAT IS THE TRIAL VALUE, NOT THE UPDATED VALUE
          ESTRESS_NPLUS1 = I_PLUS_DLAMBDA_DE_P_INV_SIG
          
          ! CALCULATE YIELD SURFACE VALUE FOR THIS UPDATED EFF. STRESS
          NF_MAT = MATMUL(P,ESTRESS_NPLUS1)
          
          DO I = 1,6
            NF_MAT_TRANS(1,I) = NF_MAT(I,1)
            
            ESTRESS_NPLUS1_TR(1,I) = ESTRESS_NPLUS1(I,1)
          END DO
          
          XI_TILDE_SQUARE = MATMUL(NF_MAT_TRANS, MATMUL(Z, NF_MAT))
          
          XI_TILDE = SQRT(XI_TILDE_SQUARE(1,1))  

          DEQPSTRAIN = DLAMBDA * XI_TILDE
          
          EQPSTRAIN = EQPSTRAIN + DEQPSTRAIN  !update plastic strain 

          YSTRESS = MAX(YSTRESS_MIN, BETA * (EQPSTRAIN+1.D-10)**ALPHA+YSTRESS_MIN)    

          DO I=1,6 ! SAME AS ESTRESS_NPLUS1
            I_PLUS_DLAMBDA_DE_P_INV_SIG_TR(1,I) = 
     +      I_PLUS_DLAMBDA_DE_P_INV_SIG(I,1)
          END DO
          
          F_DL = MATMUL(I_PLUS_DLAMBDA_DE_P_INV_SIG_TR, MATMUL(P,
     +                  I_PLUS_DLAMBDA_DE_P_INV_SIG   ))            
          
          FAILSURF = 0.5D0 * F_DL(1,1) - YSTRESS**2.0D0

          ! CHECK IF YIELD SURFACE IS REACHED UNDER TOLERANCE
          
          IF (FAILSURF .LE. F_DL(1,1)*TOL)  THEN
C             PRINT*, '******* PASS', "DLAMBDA, FAILSURF = ", 
C      +      DLAMBDA, FAILSURF
            GOTO 200
          ENDIF
          ! YIELD SURFACE NOT REACHED.
          ! CALCULATE RESIDUAL DERIVATIVE TO UPDATE VALUE OF DLAMBDA AND TRY AGAIN

          DSIG_BY_DDLAMBDA = - MATMUL(I_PLUS_DLAMBDA_DE_P_INV, MATMUL(
     +    CF, MATMUL(P, MATMUL(
     +    I_PLUS_DLAMBDA_DE_P_INV, ESTRESS_MAT))) )
          
          ! ABOVE IS 6X1
          DF_BY_DDLAMBDA_T1 = MATMUL(ESTRESS_NPLUS1_TR,
     +                        MATMUL(P, DSIG_BY_DDLAMBDA)) ! 1X1
          
          DF_BY_DDLAMBDA_T1_SCALAR = DF_BY_DDLAMBDA_T1(1,1) ! SCALAR TERM 1 OF DDF_BY_DDLAMBDA
          
          P_SIG_NPLUS1 = MATMUL(P,ESTRESS_NPLUS1) !6X1
          
          DO I=1,6
            P_SIG_NPLUS1_TR(1,I) = P_SIG_NPLUS1(I,1) !1X6
          END DO
          
          DP_BY_DDLAMBDA_T4 = ESTRESS_NPLUS1 + DLAMBDA*DSIG_BY_DDLAMBDA ! 6X1
          
          DP_BY_DDLAMBDA = MATMUL(P_SIG_NPLUS1_TR,MATMUL(Z,MATMUL(
     +                             P, DP_BY_DDLAMBDA_T4 )) ) / XI_TILDE  ! 1X1
          
          DSIGY_BY_DDLAMBDA = ALPHA*BETA*(EQPSTRAIN+1.d-10)**(ALPHA-1.d0)
     +                        *DP_BY_DDLAMBDA(1,1) ! SCALAR
C           x_temp1 = (EQPSTRAIN+0.00000001d0)**(ALPHA-1)   ! to avoid numerical errors
C           x_temp2 = DP_BY_DDLAMBDA(1,1)
          RESDER = DF_BY_DDLAMBDA_T1_SCALAR - 2.D0*YSTRESS*DSIGY_BY_DDLAMBDA ! SCALAR
          DLAMBDA_N = DLAMBDA - FAILSURF/RESDER   ! LAMBDA UPDATED
          
          IF (DLAMBDA_N .NE. DLAMBDA_N) THEN      ! SEE IF NAN HAS APPEARED
           PRINT*,"DLAMBDA IS NOT DEFINED", DLAMBDA, DLAMBDA_N
           print*, "noel, npt, KINC = ", 
     +         noel, npt, KINC
           XTEMP1 = (EQPSTRAIN)**ALPHA
           PRINT*, "RESDER, BETA, EQPSTRAIN, ALPHA, XTEMP1",
     +              RESDER, BETA, EQPSTRAIN, ALPHA, XTEMP1      
             
           PRINT*, "ITER1 = ", ITER1
           PRINT*, "EQPSTRAIN, DEQPSTRAIN = ", EQPSTRAIN, DEQPSTRAIN 
           PRINT*,F_DL,YSTRESS**2
           PRINT*,"YSTRESS",YSTRESS
           PRINT*,"FAILSURF",FAILSURF

           PNEWDT = 0.5D0
           WRITE(*,*) 'DTIME IS VERY LARGE. PLASTIC RM DID NOT CONVERGE', DTIME
           WRITE(*,*) 'REQUESTING FINER TIME STEP: PNEWDT = ', PNEWDT
           goto 201
C            CALL XIT
          END IF

          ITER1  = ITER1 + ONE    ! UPDATE ITERATION NUMBER
          
          IF(ITER1 .GE. MAXITER) THEN ! CHECK IF MAX ITERATIONS HAVE BEEN PERFORMED
            PRINT*, 'TOO MANY ITERATIONS, ITER = ', ITER1, noel, npt, kinc
            CALL XIT
          END IF
          
          IF((DLAMBDA_N - DLAMBDA) .LE. DLAMBDA*5.0D-3) THEN
             DLAMBDA = DLAMBDA_N*1.20D0 
          ELSE 
             DLAMBDA = DLAMBDA_N 
          END IF
        END DO  
 201    pnewdt = 0.5d0       


 200    do I=1,6
          DPSTRAIN(I) = DLAMBDA*NF_MAT(I,1)! INC OF PLASTIC STRAIN TENSOR
        END DO
C         EQPSTRAIN  = EQPSTRAIN + DEQPSTRAIN

        
        ESTRESS(:) = ESTRESS(:) - MATMUL(CF,DPSTRAIN) ! ESTRESS WAS TRIAL PREVIOUSLY

        
        DO I = 1,6  ! CONVERT TO MATRIX FORMAT FOR FINDING HELMOLTZ FREE ENERGY
          ESTRESS_MAT(I,1) = ESTRESS(I)
        END DO

        
        DEELAS(1:6)= DSTRAN(:) - DPSTRAIN(:)   ! INCREMENT IN ELASTIC STRAIN
        EELAS(1:6) = EELAS(1:6) - DPSTRAIN(:)   ! ELASTIC STRAIN UPDATED
        EPLAS(1:6) = EPLAS(1:6) + DPSTRAIN(1:6)
        ! EELAS WAS ACTUALLY THE TRIAL ELASTIC STRAIN VALUE, NOT PREVIOUS
        

        CALL INVERSE(CF, CF_INV, 6)
        
        H_TILDE = ALPHA*BETA*(EQPSTRAIN+1.D-10)**(ALPHA-1.D0)
        
        NF_NFTRANS = MATMUL(NF_MAT, NF_MAT_TRANS)
        ! TANGENT WITH THE M MATRIX
C         D_EP_INV = CF_INV + DLAMBDA*P + 
C      +   1.D0/(2.D0*YSTRESS*H_TILDE*(XI_TILDE+1.D-6))*MATMUL(M, NF_NFTRANS) 
C      +   - DLAMBDA/(XI_TILDE+1.D-6)**2*MATMUL(M, MATMUL(NF_NFTRANS,MATMUL(Z,P))) ! EQ (A.5)[U.DIN]
        ! TANGENT WITHOUT THE M MATRIX
C         D_EP_INV = CF_INV + DLAMBDA*P + 
C      +   1.D0/(2.D0*YSTRESS*H_TILDE*(XI_TILDE+1.D-6))*NF_NFTRANS
C      +   - DLAMBDA/(XI_TILDE+1.D-6)**2*MATMUL(NF_NFTRANS,MATMUL(Z,P)) ! EQ (A.5)[U.DIN]
C          CALL INVERSE(D_EP_INV, DDSDDE, 6) ! ELASTO PLASTIC TANGENT MATRIX    
        ! AFTER DERIVING ELASTO-PLASTIC [NOT FROM U. DIN'S WORK]

        
        CF_NFMAT = MATMUL(CF, NF_MAT)   ! SIZE 6X1
        
        NFMAT_CF_NFMAT = 0.D0
        
        DO I=1,6
          NFMAT_CF_NFMAT = NFMAT_CF_NFMAT + CF_NFMAT(I,1)*NF_MAT(I,1)  
        ENDDO
C         X_B = 1.D0/XI_TILDE * NFMAT_CF_NFMAT
        
        X_A = 1.D0/(XI_TILDE * (1.D0/XI_TILDE *
     +         NFMAT_CF_NFMAT - 2.D0*YSTRESS*H_TILDE) )
        
        DDSDDE = CF - X_A * MATMUL(CF,MATMUL(NF_NFTRANS, CF))

C         SPD=SPD+DEQPSTRAIN*YSTRESS     ! update specific plastic dissipation    
      END IF

      !!!!!!!!!!!!!!!!!END OF PLASTICITY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      x_ind = 0.d0  ! damage indicator is 0 initially
      ! check if damage has been initiated:
C       CALL CUNTZE(
C      1     estress,S21,XT,XC,YT,YC,FF1,FF2,IFF1,IFF2,IFF3,FE_EFF,MDOT,
C      2     B1,B2,I_23_5)

      xeelas11 = eelas(1)
      
      call sample_sub(estress,xeelas11,E1,S21,XT,XC,YT,YC,
     1 FF1,FF2,IFF1,IFF2,IFF3,FE_EFF )

      
      ff1ind=0.d0; ff2ind=0.d0; iff1ind=0.d0; iff2ind=0.d0; iff3ind=0.d0
      
      IF (FE_EFF.gt.1.d0) then  
        fmax = max(FF1,FF2,IFF1,IFF2,IFF3)
        
        if ((ff1.ge.1.d0 .and. ff1.gt.ff1n) .or. 
     +       (abs(ff1-fmax).le.0.001d0 .and. ff1.gt.ff1n)) then
             dcff1=dcff1+1; ff1ind=1.d0
        end if     
        if ((ff2.ge.1.d0 .and. ff.gt.ffn) .or. 
     +       (abs(ff2-fmax).le.0.001d0 .and. ff2.gt.ff2n)) then
            dcff2=dcff2+1; ff2ind=1.d0
        end if    
        if ((iff1.ge.1.d0 .and. iff1.gt.iff1n) .or. 
     +       (abs(iff1-fmax).le.0.001d0 .and. iff1.gt.iff1n)) then
            dciff1=dciff1+1; iff1ind=1.d0
        end if    
        if ((iff2.ge.1.d0 .and. iff2.gt.iff2n) .or. 
     +       (abs(iff2-fmax).le.0.001d0 .and. iff2.gt.iff2n)) then
            dciff2=dciff2+1; iff2ind=1.d0
        end if    
        if ((iff3.ge.1.d0 .and. iff3.gt.iff3n) .or. 
     +       (abs(iff3-fmax).le.0.001d0 .and. iff3.gt.iff3n)) then
            dciff3=dciff3+1; iff3ind=1.d0
        end if   
C         if ((noel.eq.12510 .and. npt.eq.1).or.(noel.eq.3226 .and. npt.eq.1)) then
        if (noel.eq.12510 .and. npt.eq.1) then    
C             print*, "***********damage activated *****"
C             print*, "noel, kinc, npt = ", noel, kinc, npt
C             print*, "fmax, FE_EFF = ", 
C      +               fmax, FE_EFF
C             print*, "ff1, ff2, iff1, iff2, iff3 = ", 
C      +               ff1, ff2, iff1, iff2, iff3
C             print*, "ff1_hist, ff2_hist, iff1_hist, iff2_hist, iff3_hist = ", 
C      +               ff1_hist, ff2_hist, iff1_hist, iff2_hist, iff3_hist    
C             x1=abs(ff1-fmax); x2=abs(ff2-fmax); x3=abs(iff1-fmax)
C             x4=abs(iff1-fmax);x5=abs(iff1-fmax)
C             print*, "x1,x2,x3,x4,x5 = ", 
C      +               x1,x2,x3,x4,x5
C             print*, "dcff1,dcff2, dciff1, dciff2, dciff3 = ", 
C      +               dcff1,dcff2, dciff1, dciff2, dciff3
        end if 
        
      end if 
      
      iff1_hist=max(IFF1, IFF1n)
      iff2_hist=max(IFF2, IFF2n)
      iff3_hist=max(IFF3, IFF3n)
      ff1_hist =max(FF1,  FF1n)
      ff2_hist =max(FF2,  FF2n)
      FE_EFF_hist = max(FE_EFF, FE_EFFn)
      
            

C       IF (FE_EFF.gt.FE_EFFn .and. FE_EFF.gt.1.d0) then  ! added damcounter condition
C         x_ind = 1.d0 ! becomes 1 if activated.
C         dam_counter = dam_counter + 1.d0
C       end if      


        w_ff1 = ff1/(ff1+ff2) 
        w_ff2 = ff2/(ff1+ff2)
        w_iff1 = iff1/(iff1+iff2+iff3) 
        w_iff2 = iff2/(iff1+iff2+iff3)
        w_iff3 = iff3/(iff1+iff2+iff3)

        
        if ((ff1+ff2).eq.0.d0) then
          w_ff1=0.d0; w_ff2=0.d0
        end if
        
        if ((iff1+iff2+iff3).eq.0.d0) then
          w_iff1=0.d0; w_iff2=0.d0; w_iff3=0.d0
        end if 


        estress_ff     = matmul(cff,   eelas)
        estress_iff    = matmul(ciff,  eelas)
        estress_iff3   = matmul(ciff3, eelas)


        ! elastic driving energies
        psi_ff1  = w_ff1*0.5d0*dot_product(estress_ff,     eelas)
        psi_ff2  = w_ff2*0.5d0*dot_product(estress_ff,     eelas)
        psi_iff1 = w_iff1*0.5d0*dot_product(estress_iff,   eelas)
        psi_iff2 = w_iff2*0.5d0*dot_product(estress_iff,   eelas)
        psi_iff3 = w_iff3*( 0.5d0*dot_product(estress_iff, eelas)) 
     +                    + 0.5d0*dot_product(estress_iff3,eelas)

        ! DRIVING 2
C         psi_ff1  = gd5*0.5d0*dot_product(estress_ff,  eelas)
C         psi_ff2  = gd4*0.5d0*dot_product(estress_ff,  eelas)
C         psi_iff1 = gd4*gd5*gd2*gd3*0.5d0*dot_product(estress_iff, eelas)
C         psi_iff2 = gd4*gd5*gd1*gd3*0.5d0*dot_product(estress_iff, eelas)
C         psi_iff3 = gd4*gd5*gd1*gd2*( 0.5d0*dot_product(estress_iff, eelas)) 
C      +                    + gd4*gd5*0.5d0*dot_product(estress_iff3,eelas)
        if (estress(1).lt.0.d0) psi_ff1=0.d0
C         if (estress(2).lt.0.d0) psi_iff1=0.d0
        if (estress(1).gt.0.d0) psi_ff2=0.d0
C         if (estress(2).gt.0.d0) psi_iff2=0.d0
        ! incremental plastic driving density

        INC_PCD_4 = w_ff1*0.5d0*dot_product(estress_ff,  DPSTRAIN)
        INC_PCD_5 = w_ff2*0.5d0*dot_product(estress_ff,  DPSTRAIN)
        INC_PCD_1 = w_iff1*0.5d0*dot_product(estress_iff, DPSTRAIN)
        INC_PCD_2 = w_iff2*0.5d0*dot_product(estress_iff, DPSTRAIN)
        INC_PCD_3 = w_iff3*( 0.5d0*dot_product(estress_iff, DPSTRAIN)) 
     +                    + 0.5d0*dot_product(estress_iff3,DPSTRAIN)


        ! driving 2 plastic:
C         INC_PCD_4 = gd5*0.5d0*dot_product(estress_ff,  DPSTRAIN)
C         INC_PCD_5 = gd4*0.5d0*dot_product(estress_ff,  DPSTRAIN)
C         INC_PCD_1 = gd4*gd5*gd2*gd3*0.5d0*dot_product(estress_iff, DPSTRAIN)
C         INC_PCD_2 = gd4*gd5*gd1*gd3*0.5d0*dot_product(estress_iff, DPSTRAIN)
C         INC_PCD_3 = gd4*gd5*gd1*gd2*( 0.5d0*dot_product(estress_iff, DPSTRAIN)) 
C      +                    + gd4*gd5*0.5d0*dot_product(estress_iff3,DPSTRAIN)
        if (estress(1).lt.0.d0) INC_PCD_4=0.d0
C         if (estress(2).lt.0.d0) INC_PCD_1=0.d0
        if (estress(1).gt.0.d0) INC_PCD_5=0.d0
C         if (estress(2).gt.0.d0) INC_PCD_2=0.d0


C         INC_PCD_1 = w_iff1*(estress(2)*DPSTRAIN(2) +  
C      +  estress(3)*DPSTRAIN(3) + estress(6)*DPSTRAIN(6))
C         INC_PCD_2 = w_iff2*(estress(2)*DPSTRAIN(2) + 
C      +  estress(3)*DPSTRAIN(3) + estress(6)*DPSTRAIN(6))
C         INC_PCD_3 = w_iff3*(estress(2)*DPSTRAIN(2) +  
C      +  estress(3)*DPSTRAIN(3) +estress(6)*DPSTRAIN(6)) + 
C      +  (estress(4)*DPSTRAIN(4) + estress(5)*DPSTRAIN(5))   
C         INC_PCD_4 = w_ff1*(estress(1)*DPSTRAIN(1))
C         INC_PCD_5 = w_ff2*(estress(1)*DPSTRAIN(1))
        ! total plastic crack driving density
        

        pcd1 = pcd1 + INC_PCD_1
        pcd2 = pcd2 + INC_PCD_2
        pcd3 = pcd3 + INC_PCD_3
        pcd4 = pcd4 + INC_PCD_4
        pcd5 = pcd5 + INC_PCD_5
C       end if

      !
      psip1=0.d0; psip3=0.d0; psip4=0.d0; psip2=0.d0; psip5=0.d0


      psip1 = psi_iff1 + pcd1
      psip2 = psi_iff2 + pcd2
      psip3 = psi_iff3 + pcd3
      psip4 = psi_ff1  + pcd4
      psip5 = psi_ff2  + pcd5


      if (dcff1.le.1) then
        Hc4 = psip4
      else 
        Hc4 = statev(21)
      end if

      if (dcff2.le.1) then
        Hc5 = psip5
      else 
        Hc5 = statev(22)
      end if

      if (dciff1.le.1) then
        Hc1 = psip1
      else 
        Hc1 = statev(18)
      end if

      if (dciff2.le.1) then
        Hc2 = psip2
      else 
        Hc2 = statev(19)
      end if

      if (dciff3.le.1) then
        Hc3 = psip3
      else 
        Hc3 = statev(20)
      end if    


C       if (dam_counter.eq.1.d0) then ! fist time activated
C         !      Hc_{n+1} = (H_{n+1} + H_{n})/2
C         Hc1 = psip1
C         Hc2 = psip2
C         Hc3 = psip3
C         Hc4 = psip4
C         Hc5 = psip5
C C         Hc1 = (psip1 + statev(8))/2.d0
C C         Hc2 = (psip2 + statev(10))/2.d0
C C         Hc3 = (psip3 + statev(12))/2.d0
C C         Hc4 = (psip4 + statev(14))/2.d0
C C         Hc5 = (psip5 + statev(16))/2.d0
C       else if (dam_counter.lt.1.d0) then  ! not activated
C         !      Hc_{n+1} = H_{n+1}
C         Hc1 = psip1
C         Hc2 = psip2
C         Hc3 = psip3
C         Hc4 = psip4
C         Hc5 = psip5
C       else ! already activated (dam_count > 1)
C         !      Hc_{n+1} = Hc_{n}
C         Hc1 = statev(18) 
C         Hc2 = statev(19)
C         Hc3 = statev(20)
C         Hc4 = statev(21)
C         Hc5 = statev(22)  
C       end if
      

      H1 = max(0.d0, statev(8 ), psip1-Hc1)*iff1ind
      H2 = max(0.d0, statev(10), psip2-Hc2)*iff2ind
      H3 = max(0.d0, statev(12), psip3-Hc3)*iff3ind
      H4 = max(0.d0, statev(14), psip4-Hc4)*ff1ind
      H5 = max(0.d0, statev(16), psip5-Hc5)*ff2ind
      

      ! calculate weights from history values
      
      w_ff1_hist = ff1_hist/(ff1_hist+ff2_hist) 
      w_ff2_hist = ff2_hist/(ff1_hist+ff2_hist)
      w_iff1_hist = iff1_hist/(iff1_hist+iff2_hist+iff3_hist) 
      w_iff2_hist = iff2_hist/(iff1_hist+iff2_hist+iff3_hist)
      w_iff3_hist = iff3_hist/(iff1_hist+iff2_hist+iff3_hist)

      
      stress(1)   = (gd4+xk)*(gd5+xk)*CF(1,1)*EELAS(1) + 
     +              (gd3+xk)*(gd1+xk)*(gd2+xk)*
     +              (EELAS(2)*CF(1,2) + EELAS(3)*CF(1,3))   
      stress(2)   = (gd3+xk)*(gd1+xk)*(gd2+xk)*estress(2)
      stress(3)   = (gd3+xk)*(gd1+xk)*(gd2+xk)*estress(3)
      stress(4)   = (gd3+xk)*estress(4)
      stress(5)   = (gd3+xk)*estress(5)
      stress(6)   = (gd3+xk)*(gd1+xk)*(gd2+xk)*estress(6)
      

      XGDF=(gd4+xk)*(gd5+xk)
      

      xpsi_m = dot_product(estress_iff, eelas)*0.5D0
      xpsi_s = dot_product(estress_iff3,eelas)*0.5D0
      

      PsimByPsim0 = (gd3*gd1*gd2*xpsi_m + gd3*xpsi_s)/(xpsi_m + xpsi_s)

      
      ddsdde(1,1) = (gd4+xk)*(gd5+xk)*ddsdde(1,1)
      ddsdde(2,2) = (gd3+xk)*(gd1+xk)*(gd2+xk)*ddsdde(2,2)
      ddsdde(3,3) = (gd3+xk)*(gd1+xk)*(gd2+xk)*ddsdde(3,3)
      ddsdde(4,4) = (gd3+xk)*ddsdde(4,4)
      ddsdde(5,5) = (gd3+xk)*ddsdde(5,5)
      ddsdde(6,6) = (gd3+xk)*(gd1+xk)*(gd2+xk)*ddsdde(6,6)
      ddsdde(1,2) = (gd3+xk)*(gd1+xk)*(gd2+xk)*ddsdde(1,2)
      ddsdde(2,1) = (gd3+xk)*(gd1+xk)*(gd2+xk)*ddsdde(2,1)
      ddsdde(1,3) = (gd3+xk)*(gd1+xk)*(gd2+xk)*ddsdde(1,3)
      ddsdde(3,1) = (gd3+xk)*(gd1+xk)*(gd2+xk)*ddsdde(3,1)
      ddsdde(3,2) = (gd3+xk)*(gd1+xk)*(gd2+xk)*ddsdde(3,2)
      ddsdde(2,3) = (gd3+xk)*(gd1+xk)*(gd2+xk)*ddsdde(2,3) 
      

      ! REGULARIZING THE DRIVING ENERGY
C         H1V = (VISCO/(VISCO+DTIME)*STATEV(55)+
C      +       DTIME/(VISCO+DTIME)*H1)
C         H2V = (VISCO/(VISCO+DTIME)*STATEV(56)+
C      +       DTIME/(VISCO+DTIME)*H2)
C         H3V = (VISCO/(VISCO+DTIME)*STATEV(57)+
C      +       DTIME/(VISCO+DTIME)*H3)
C         VISCO_FF = VISCO ! INCREASE VISCOSITY FOR FIRBRE FRACTURE
C         H4V= (VISCO_FF/(VISCO_FF+DTIME)*STATEV(58)+
C      +       DTIME/(VISCO_FF+DTIME)*H4)
C         H5V= (VISCO_FF/(VISCO_FF+DTIME)*STATEV(59)+
C      +       DTIME/(VISCO_FF+DTIME)*H5)

        
        H1V=H1; H2V=H2; H3V=H3; H4V=H4; H5V=H5

C       if (dam_counter.ge.1.d0) then
C C         stress(1) = ((gd4+xk)*w_ff1 + (gd5+xk)*w_ff2)*estress(1) 
C C         stress(2) = ((gd1+xk)*w_iff1+(gd2+xk)*w_iff2+(gd3+xk)*w_iff3)
C C      +                *estress(2)
C C         stress(3) = ((gd1+xk)*w_iff1+(gd2+xk)*w_iff2+(gd3+xk)*w_iff3)
C C      +                *estress(3)
C C         stress(4) = (gd3+xk)*estress(4)
C C         stress(5) = (gd3+xk)*estress(5)
C C         stress(6) = ((gd1+xk)*w_iff1+(gd2+xk)*w_iff2+(gd3+xk)*w_iff3)
C C      +                *estress(6)

C         stress(1)   = (gd4+xk)*(gd5+xk)*estress(1)
C         stress(2)   = (gd3+xk)*(gd1+xk)*(gd2+xk)*estress(2)
C         stress(3)   = (gd3+xk)*(gd1+xk)*(gd2+xk)*estress(3)
C         stress(4)   = (gd3+xk)*estress(4)
C         stress(5)   = (gd3+xk)*estress(5)
C         stress(6)   = (gd3+xk)*(gd1+xk)*(gd2+xk)*estress(6)
C         ! REGULARIZING THE DRIVING ENERGY
C         H1V = (VISCO/(VISCO+DTIME)*STATEV(55)+
C      +       DTIME/(VISCO+DTIME)*H1)
C         H2V = (VISCO/(VISCO+DTIME)*STATEV(56)+
C      +       DTIME/(VISCO+DTIME)*H2)
C         H3V = (VISCO/(VISCO+DTIME)*STATEV(57)+
C      +       DTIME/(VISCO+DTIME)*H3)
C         VISCO_FF = VISCO ! INCREASE VISCOSITY FOR FIRBRE FRACTURE
C         H4V= (VISCO_FF/(VISCO_FF+DTIME)*STATEV(58)+
C      +       DTIME/(VISCO_FF+DTIME)*H4)
C         H5V= (VISCO_FF/(VISCO_FF+DTIME)*STATEV(59)+
C      +       DTIME/(VISCO_FF+DTIME)*H5)
C       else
C         stress = estress 
C         H1V=H1; H2V=H2; H3V=H3; H4V=H4; H5V=H5
C       end if  
C       if (dam_counter.ge.1.d0) then
C C         DDSDDE(1,1) = DDSDDE(1,1)*((gd4+xk)*w_ff1 + (gd5+xk)*w_ff2)
C C         DDSDDE(1,2) = DDSDDE(1,2)*((gd4+xk)*w_ff1 + (gd5+xk)*w_ff2)
C C         DDSDDE(1,3) = DDSDDE(1,3)*((gd4+xk)*w_ff1 + (gd5+xk)*w_ff2)

C C         DDSDDE(2,1) = DDSDDE(2,1)
C C      +               *((gd1+xk)*w_iff1+(gd2+xk)*w_iff2+(gd3+xk)*w_iff3)
C C         DDSDDE(2,2) = DDSDDE(2,2)
C C      +               *((gd1+xk)*w_iff1+(gd2+xk)*w_iff2+(gd3+xk)*w_iff3)
C C         DDSDDE(2,3) = DDSDDE(2,3)
C C      +               *((gd1+xk)*w_iff1+(gd2+xk)*w_iff2+(gd3+xk)*w_iff3)  

C C         DDSDDE(3,1) = DDSDDE(3,1)
C C      +               *((gd1+xk)*w_iff1+(gd2+xk)*w_iff2+(gd3+xk)*w_iff3)    
C C         DDSDDE(3,2) = DDSDDE(3,2)
C C      +               *((gd1+xk)*w_iff1+(gd2+xk)*w_iff2+(gd3+xk)*w_iff3)    
C C         DDSDDE(3,3) = DDSDDE(3,3)
C C      +               *((gd1+xk)*w_iff1+(gd2+xk)*w_iff2+(gd3+xk)*w_iff3)

C C         DDSDDE(4,4) = DDSDDE(4,4)*(gd3+xk)   
C C         DDSDDE(5,5) = DDSDDE(5,5)*(gd3+xk)   
C C         DDSDDE(6,6) = DDSDDE(6,6)
C C      +               *((gd1+xk)*w_iff1+(gd2+xk)*w_iff2+(gd3+xk)*w_iff3)
        
C         ddsdde(1,1) = (gd4+xk)*(gd5+xk)*ddsdde(1,1)
C         ddsdde(2,2) = (gd3+xk)*(gd1+xk)*(gd2+xk)*ddsdde(2,2)
C         ddsdde(3,3) = (gd3+xk)*(gd1+xk)*(gd2+xk)*ddsdde(3,3)
C         ddsdde(4,4) = (gd3+xk)*ddsdde(4,4)
C         ddsdde(5,5) = (gd3+xk)*ddsdde(5,5)
C         ddsdde(6,6) = (gd3+xk)*(gd1+xk)*(gd2+xk)*ddsdde(6,6)
C         ddsdde(1,2) = (gd3+xk)*(gd1+xk)*(gd2+xk)*ddsdde(1,2)
C         ddsdde(2,1) = (gd3+xk)*(gd1+xk)*(gd2+xk)*ddsdde(2,1)
C         ddsdde(1,3) = (gd3+xk)*(gd1+xk)*(gd2+xk)*ddsdde(1,3)
C         ddsdde(3,1) = (gd3+xk)*(gd1+xk)*(gd2+xk)*ddsdde(3,1)
C         ddsdde(3,2) = (gd3+xk)*(gd1+xk)*(gd2+xk)*ddsdde(3,2)
C         ddsdde(2,3) = (gd3+xk)*(gd1+xk)*(gd2+xk)*ddsdde(2,3) 
C       end if

C       ddsdde = CFDEG
C       DDSDDE = matmul(M_INV, DDSDDE)      


!     Collect information from UEL 
      phi1 = UserVar(npt,1,noel)
      phi2 = UserVar(npt,3,noel)
      phi3 = UserVar(npt,5,noel)
      phi4 = UserVar(npt,7,noel)
      phi5 = UserVar(npt,9,noel)


!     information transfer to UEL - SEND REGULARIZED DRIVING
      UserVar(npt,2,noel) = H1V
      UserVar(npt,4,noel) = H2V
      UserVar(npt,6,noel) = H3V
      UserVar(npt,8,noel) = H4V
      UserVar(npt,10,noel)= H5V
C       print*, "noel in umat = ", noel
C       print*, "coords = ", coords


!     output      
      statev(1:ntens) = estress
      statev(7)       = phi1
      statev(8)       = H1
      statev(9)       = phi2
      statev(10)      = H2
      statev(11)      = phi3
      statev(12)      = H3
      statev(13)      = phi4
      statev(14)      = H4
      statev(15)      = phi5
      statev(16)      = H5 ! 16th sdv here.
      statev(17)      = dam_counter
      statev(18)      = Hc1
      statev(19)      = Hc2
      statev(20)      = Hc3
      statev(21)      = Hc4
      statev(22)      = Hc5
      statev(23)      = max(FE_EFF, FE_EFFn)
      statev(24)      = x_ind
      statev(25)      = iff1_hist
      statev(26)      = iff2_hist
      statev(27)      = iff3_hist
      statev(28)      = ff1_hist
      statev(29)      = ff2_hist
      statev(30:35)   = EELAS
      statev(36:41)   = EPLAS
      statev(42)      = pcd1
      statev(43)      = pcd2
      statev(44)      = pcd3
      statev(45)      = pcd4
      statev(46)      = pcd5
      statev(47)      = psi_iff1
      statev(48)      = psi_iff2
      statev(49)      = psi_iff3
      statev(50)      = psi_ff1
      statev(51)      = psi_ff2
      statev(52)      = FAILSURF
      statev(53)      = EQPSTRAIN
      statev(54)      = YSTRESS
      STATEV(55)      = H1V
      STATEV(56)      = H2V
      STATEV(57)      = H3V
      STATEV(58)      = H4V
      STATEV(59)      = H5V
      STATEV(60)      = w_iff1
      STATEV(61)      = w_iff2
      STATEV(62)      = w_iff3
      STATEV(63)      = w_ff1
      STATEV(64)      = w_ff2
      statev(65)      = dcff1
      statev(66)      = dcff2
      statev(67)      = dciff1
      statev(68)      = dciff2
      statev(69)      = dciff3
      statev(70)      = XGDF
      statev(71)      = g1*g2*g3
      statev(72)      = g3
      statev(73)      = xpsi_m
      statev(74)      = xpsi_s


C       statev(71)      = (gd3*gd1*gd2*xpsi_m + gd3*xpsi_s)/(xpsi_m + xpsi_s)
C       statev(72)      = xpsi_m + xpsi_s

      
      if (noel.eq.12510 .and. npt.eq.1 .and. FE_EFF.gt.1.d0) then
C         print*, "*********** umat end *****"
C         print*, "noel, kinc, npt = ", noel, kinc, npt
C         print*, "H1,H2,H3,H4,H5 = ", H1,H2,H3,H4,H5
C         print*, "fmax = ", fmax
C         print*, "phi1,phi2,phi3,phi4,phi5 = ", 
C      +           phi1,phi2,phi3,phi4,phi5
C         print*, "ff1ind,ff2ind, iff1ind, iff2ind, iff3ind = ", 
C      +           ff1ind,ff2ind, iff1ind, iff2ind, iff3ind
C         print*, "dcff1,dcff2, dciff1, dciff2, dciff3 = ", 
C      +           dcff1,dcff2, dciff1, dciff2, dciff3
C         print*, "ff1,ff2, iff1, iff2, iff3 = ", 
C      +           ff1,ff2, iff1, iff2, iff3
C         print*, "FE_EFF = ", FE_EFF


      end if 
     
      return
      end


      subroutine sample_sub(S,xeelas11,E1,S21,XT,XC,YT,YC,
     1 FF1,FF2,IFF1,IFF2,IFF3,FE_EFF)

      real*8 xeelas11, S(6), E1, S21,XT,XC,YT,YC,
     1 FF1,FF2,IFF1,IFF2,IFF3,FE_EFF,B1,B2,MDOT,I_23_5

      PARAMETER (ZERO=0.D0, ONE=1.D0, B1=1.266D0, B2=0.1D0, MDOT=2.99D0,
     1           VF=0.61D0  )


      SIG11=S(1); SIG22=S(2); SIG33=S(3)
      SIG12=S(4); SIG13=S(5); SIG23=S(6)
      SIG11_STAR = xeelas11*E1


C       print*, "inside the subroutine - ", S,xeelas11,E1,S21,XT,XC,YT,YC
C       print*, SIG11,SIG22,SIG33,SIG12,SIG13,SIG23,SIG11_STAR
      

      IF(SIG11_STAR .GE. ZERO) THEN
          FF1 = SIG11_STAR/XT;      FF2 = 0.D0
C           FF1 = SIG11_STAR/XT;      FF2 = 0.D0
      ELSE
          FF2 = ABS(SIG11_STAR)/XC; FF1 = 0.D0
C           FF2 = ABS(SIG11_STAR)/XC; FF1 = 0.D0
      END IF
C       ff1 = SIG11_STAR/xt; ff2 = -1.d0 * SIG11_STAR/xc
      

      IFF1 = (SIG22 + SIG33 + SQRT((SIG22 - SIG33)**2 
     1  + 4. * SIG23**2)) / (2.0d0*YT)
      

      IFF2 = (B1*SQRT((SIG22 - SIG33)**2 + 4.D0*SIG23**2) 
     1     + (B1-1.D0)*(SIG22+SIG33))/(-YC)
C       print*, "FF1, FF2 = ", FF1, FF2
C       IF(SIG22 .GE. ZERO) THEN
C         IFF1 = (SIG22 + SIG33 + SQRT((SIG22 - SIG33)**2 
C      1  + 4. * SIG23**2)) / (2.0d0*YT)
C         IFF2 = 0.D0
C       ELSE 
C         IFF2 = ABS((B1 * SQRT((SIG22 - SIG33)**2 + 4.D0 * SIG23**2) 
C      1  + (B1 - 1.D0) * (SIG22 + SIG33)) /YC)
C         IFF1 = 0.D0
C       END IF 
        

        I_23_5 = 2.d0*SIG22*SIG12**2 + 2.d0*SIG33*SIG13**2 
     1  + 4.d0*SIG23*SIG13*SIG12
C         IFF3 = SQRT((SQRT(B2**2 * I_23_5**2 + 4. * S21**2 *
C      1  (SIG13**2 +SIG12**2)**2) + B2*I_23_5)/(2.*S21**3))
        X1 = B2**2 * I_23_5**2 + 4.d0*S21**2 *(SIG13**2 +SIG12**2)**2
        X2 = B2*I_23_5
        IF (ABS(X2).LT.1.D-8) X2 = 0.D0
        X3 = (SQRT(X1) + X2) / (2.D0 * S21**3)
        IF (ABS(X3).LT.1.D-8) X3 = 0.D0
        IF (X3.LT.0.D0) X3 = 0.D0
        IFF3 = SQRT(X3)
        if (IFF1.lt.0.d0) IFF1=0.d0
        if (IFF2.lt.0.d0) IFF2=0.d0
        if (IFF3.lt.0.d0) IFF3=0.d0
        if (FF1.lt.0.d0) FF1=0.d0
        if (FF2.lt.0.d0) FF2=0.d0

        
        FE_EFF = max(FF1, 0.d0)**MDOT + max(FF2, 0.d0)**MDOT +
     +   max(IFF1, 0.d0)**MDOT + max(IFF2, 0.d0)**MDOT +
     +   max(IFF3, 0.d0)**MDOT
C         FE_EFF = IFF1**MDOT + IFF2**MDOT + IFF3**MDOT
        IF(FE_EFF.NE.FE_EFF) THEN
          PRINT*, "INSIDE THE CUNTZE FUNCTION! THIS MEANS SOMETHING HAS GONE WRONG"
          PRINT*, "FE_EFF IS NOT DEFINED"
          print*, "FF1, FF2 = ", FF1, FF2
          PRINT*, "IFF1, IFF2, IFF3",  IFF1, IFF2, IFF3
          PRINT*, "MDOT, B1, B2, FE_EFF = ", MDOT, B1, B2, FE_EFF
          PRINT*, "S = ", S
          PRINT*, "I_23_5 = ", I_23_5
          PRINT*, "X1, X2, X3 = ", X1, X2, X3
        ENDIF
      return
      end
      

C******************************************************************************
C CUNTZE FAILRE CRITERIA******************************************************
C******************************************************************************
      
      SUBROUTINE CUNTZE(
     1      S,S21,XT,XC,YT,YC,FF1,FF2,IFF1,IFF2,IFF3,FE_EFF,MDOT
     2      B1,B2, I_23_5)
   
C       INCLUDE 'VABA_PARAM.INC'
      DOUBLE PRECISION
     1 S21,XT,XC,YT,YC,S(6),
     2 FF1,FF2,IFF1,IFF2,IFF3,FE_EFF,
     3 B1,B2, MDOT,I_23_5,TSTRAIN(6),STRESS(6)    
      PARAMETER (ZERO=0.D0, ONE=1.D0)
C       estress was defined in abaqus notation - 11, 22, 33, 12, 23, 13
C       tot_strain also
      SIG11=S(1)
      SIG22=S(2)
      SIG33=S(3)
      SIG12=S(4)
      SIG13=S(5)
      SIG23=S(6)      
C       TSN11=TSTRAIN(1)
C       TSN22=TSTRAIN(2)
C       TSN33=TSTRAIN(3)
C       TSN12=TSTRAIN(4)
C       TSN13=TSTRAIN(5)
C       TSN23=TSTRAIN(6)
      
      B1 = 1.09D0
      B2 = 0.13D0
      MDOT = 3.1D0
      ! INITIALIZE ALL AS ZERO
C       FF1=0.D0; FF2=0.D0
      
      IF(SIG11 .GE. ZERO) THEN
          FF1 = SIG11/XT;      FF2 = 0.D0
      ELSE
          FF2 = ABS(SIG11)/XC; FF1 = 0.D0
      END IF
C       Wrong expressions. Check paper Petersen 2016 
      IF(SIG22 .GE. ZERO) THEN
        IFF1 = (SIG22 + SIG33 + SQRT((SIG22 - SIG33)**2 
     1  + 4. * SIG23**2)) / (2.0d0*YT)
        IFF2 = 0.D0
      ELSE 
        IFF2 = ABS((B1 * SQRT((SIG22 - SIG33)**2 + 4.D0 * SIG23**2) 
     1  + (B1 - 1.D0) * (SIG22 + SIG33)) /YC)
        IFF1 = 0.D0
      END IF 
        I_23_5 = 2.*SIG22*SIG12**2 + 2.*SIG33*SIG13**2 
     1  + 4.*SIG23*SIG13*SIG12
C         IFF3 = SQRT((SQRT(B2**2 * I_23_5**2 + 4. * S21**2 *
C      1  (SIG13**2 +SIG12**2)**2) + B2*I_23_5)/(2.*S21**3))
        X1 = B2**2 * I_23_5**2 + 4. * S21**2 *(SIG13**2 +SIG12**2)**2
        X2 = B2*I_23_5
        IF (ABS(X2).LT.1.D-8) X2 = 0.D0
        X3 = (SQRT(X1) + X2) / (2.D0 * S21**3)
        IF (ABS(X3).LT.1.D-8) X3 = 0.D0
        IF (X3.LT.0.D0) X3 = 0.D0
        IFF3 = SQRT(X3)

   
        FE_EFF = FF1**MDOT+FF2**MDOT+IFF1**MDOT+IFF2**MDOT+IFF3**MDOT
C         FE_EFF = IFF1**MDOT + IFF2**MDOT + IFF3**MDOT
        IF(FE_EFF.NE.FE_EFF) THEN
          PRINT*, "INSIDE THE CUNTZE FUNCTION! THIS MEANS SOMETHING HAS GONE WRONG"
          PRINT*, "FE_EFF IS NOT DEFINED"
          PRINT*, "IFF1, IFF2, IFF3",  IFF1, IFF2, IFF3
C           PRINT*, "MDOT, B1, B2 = ", MDOT, B1, B2
          PRINT*, "S = ", S
          PRINT*, "I_23_5 = ", I_23_5
          PRINT*, "X1, X2, X3 = ", X1, X2, X3
        ENDIF

        

      RETURN
      END SUBROUTINE CUNTZE

C************************************************************
      SUBROUTINE TV_MAT_V(TV1,MAT,V2,C1)
        REAL(8), INTENT(OUT) :: C1
        REAL(8), INTENT(IN) :: TV1(6),V2(6),MAT(6,6)
        REAL(8) :: F(1,1),MV1(1,6),MV2(6,1)

        DO I = 1,6
          MV1(1,I) = TV1(I)
          MV2(I,1) = V2(I)
        ENDDO
      
        F = MATMUL(MATMUL(MV1,MAT),MV2)
        C1 = F(1,1)

        RETURN
      END SUBROUTINE TV_MAT_V

C**********************************************************       
      SUBROUTINE INVERSE(AA,C,N)
!============================================================
! INVERSE MATRIX
! METHOD: BASED ON DOOLITTLE LU FACTORIZATION FOR AX=B
! ALEX G. DECEMBER 2009
!-----------------------------------------------------------
! INPUT ...
! A(N,N) - ARRAY OF COEFFICIENTS FOR MATRIX A
! N      - DIMENSION
! OUTPUT ...
! C(N,N) - INVERSE MATRIX OF A
! COMMENTS ...
! THE ORIGINAL MATRIX A(N,N) WILL BE DESTROYED 
! DURING THE CALCULATION
!===========================================================
        IMPLICIT NONE 
        INTEGER N
        DOUBLE PRECISION A(N,N), C(N,N), AA(N,N)
        DOUBLE PRECISION L(N,N), U(N,N), B(N), D(N), X(N)
        DOUBLE PRECISION COEFF
        INTEGER I, J, K

! STEP 0: INITIALIZATION FOR MATRICES L AND U AND B
! FORTRAN 90/95 ALOOWS SUCH OPERATIONS ON MATRICES
        L=0.0
        U=0.0
        B=0.0
        A = AA
! STEP 1: FORWARD ELIMINATION
        DO K=1, N-1
            DO I=K+1,N
              COEFF=A(I,K)/A(K,K)
              L(I,K) = COEFF
                DO J=K+1,N
                   A(I,J) = A(I,J)-COEFF*A(K,J)
                END DO
             END DO
        END DO

! STEP 2: PREPARE L AND U MATRICES 
! L MATRIX IS A MATRIX OF THE ELIMINATION COEFFICIENT
! + THE DIAGONAL ELEMENTS ARE 1.0
       DO I=1,N
         L(I,I) = 1.0
       END DO
! U MATRIX IS THE UPPER TRIANGULAR PART OF A
        DO J=1,N
          DO I=1,J
              U(I,J) = A(I,J)
          END DO
        END DO

! STEP 3: COMPUTE COLUMNS OF THE INVERSE MATRIX C
          DO K=1,N
              B(K)=1.0
              D(1) = B(1)
! STEP 3A: SOLVE LD=B USING THE FORWARD SUBSTITUTION
                  DO I=2,N
                      D(I)=B(I)
                      DO J=1,I-1
                          D(I) = D(I) - L(I,J)*D(J)
                      END DO
                  END DO
! STEP 3B: SOLVE UX=D USING THE BACK SUBSTITUTION
                  X(N)=D(N)/U(N,N)
                  DO I = N-1,1,-1
                       X(I) = D(I)
                       DO J=N,I+1,-1
                          X(I)=X(I)-U(I,J)*X(J)
                       END DO
                       X(I) = X(I)/U(I,I)
                  END DO
! STEP 3C: FILL THE SOLUTIONS X(N) INTO COLUMN K OF C
                  DO I=1,N
                      C(I,K) = X(I)
                  END DO
                  B(K)=0.0
          END DO
       END SUBROUTINE INVERSE      
      
      
      SUBROUTINE TRANSFORMATION_MATRIX(T,T_INV,PROPS,NPROPS)
        
C       INCLUDE 'VABA_PARAM.INC'
      DOUBLE PRECISION T(6,6),T_INV(6,6),PROPS(NPROPS)
      DOUBLE PRECISION THETA
      
        THETA = PROPS(17) * ATAN(1.0_8) / 45.0
        
                
        ! COMPUTE TRANSFORMATION MATRIX
        ! T(1,1) = (COS(THETA))**2
        ! T(1,2) = (SIN(THETA))**2
        ! T(1,4) = SIN(THETA) * COS(THETA)
        ! T(2,1) = T(1,2)
        ! T(2,2) = T(1,1)
        ! T(2,4) = -SIN(THETA) * COS(THETA)
        ! T(3,3) = 1.0D0
        ! T(4,1) = -2. * SIN(THETA) * COS(THETA) 
        ! T(4,2) = 2. * SIN(THETA) * COS(THETA)
        ! T(4,4) = (COS(THETA))**2 - (SIN(THETA))**2
        ! T(5,5) = COS(THETA)
        ! T(5,6) = -SIN(THETA)
        ! T(6,5) = SIN(THETA)
        ! T(6,6) = COS(THETA)
        
        T(1,1) = 0.50D0
        T(1,2) = 0.50D0
        T(1,4) = 0.50D0
        T(2,1) = T(1,2)
        T(2,2) = T(1,1)
        T(2,4) = -0.50D0
        T(3,3) = 1.0D0
        T(4,1) = -1.0D0 
        T(4,2) = 1.0D0
        T(4,4) = 0.0d0
        T(5,5) = 1/sqrt(2.)
        T(5,6) = -1/sqrt(2.)
        T(6,5) = 1/sqrt(2.)
        T(6,6) = 1/sqrt(2.)


        
        ! INVERSE TRANSFORMATION_MATRIX
        ! T_INV(1,1) = (COS(THETA))**2
        ! T_INV(1,2) = (SIN(THETA))**2
        ! T_INV(1,4) = -SIN(THETA) * COS(THETA)
        ! T_INV(2,1) = T(1,2)
        ! T_INV(2,2) = T(1,1)
        ! T_INV(2,4) = SIN(THETA) * COS(THETA)
        ! T_INV(3,3) = 1.0D0
        ! T_INV(4,1) = 2. * SIN(THETA) * COS(THETA) 
        ! T_INV(4,2) = -2. * SIN(THETA) * COS(THETA)
        ! T_INV(4,4) = (COS(THETA))**2 - (SIN(THETA))**2
        ! T_INV(5,5) = COS(THETA)
        ! T_INV(5,6) = SIN(THETA)
        ! T_INV(6,5) = -SIN(THETA)
        ! T_INV(6,6) = COS(THETA)
        
        T_INV(1,1) = 0.50D0
        T_INV(1,2) = 0.50D0
        T_INV(1,4) = -0.50D0
        T_INV(2,1) = T_INV(1,2)
        T_INV(2,2) = T_INV(1,1)
        T_INV(2,4) = 0.50D0
        T_INV(3,3) = 1.0D0
        T_INV(4,1) = 1.0D0 
        T_INV(4,2) = -1.0D0
        T_INV(4,4) = 0.0d0
        T_INV(5,5) = 1/sqrt(2.)
        T_INV(5,6) = 1/sqrt(2.)
        T_INV(6,5) = -1/sqrt(2.)
        T_INV(6,6) = 1/sqrt(2.)
        
        
        ! PRINT*, 'T IN LOOP',T
        ! print*, "theta",THETA
        
      RETURN  
      END SUBROUTINE TRANSFORMATION_MATRIX

C**********************************************************************
      SUBROUTINE onem(A)

C     THIS SUBROUTINE STORES THE IDENTITY MATRIX IN THE 
C     3 BY 3 MATRIX [A]
C**********************************************************************

        REAL*8 A(3,3)
        DATA ZERO/0.D0/
        DATA ONE/1.D0/

      DO 1 I=1,3
        DO 1 J=1,3
          IF (I .EQ. J) THEN
              A(I,J) = 1.0
            ELSE
              A(I,J) = 0.0
            ENDIF
1       CONTINUE

      RETURN
      END

      subroutine dyadic(vector1,vector2, vlen, dyadicprod)
                  
                        integer  vlen, i, j
                        real*8 vector1(vlen),vector2(vlen)
                        real*8 dyadicprod(vlen,vlen)
                
                        do i = 1, vlen
                              do j = 1, vlen
                              dyadicprod(i,j) = vector1(i) * vector2(j)
                              end do
                        end do

                        return
      end subroutine dyadic