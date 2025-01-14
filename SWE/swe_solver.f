      program SWE
      !*******************************************************************************
!!     This code is solving the SWEs coupled with the suspendend load equation
!     After solving the flow field the bed elevation will be updated using Exner equ.
!     Finite volume is used, with Roe's solver to calculate the fluxes at the cell faces
!     MUSCL scheme is used to extend the solution to second order accuracy
!     This code is made to work on structure grid and the used schemes are Exiplicite schems
!     For more than one grain size for the sediment we can use D50 or the code must go through some modifications
!********************************************************************************************  
!     read the grid
      integer i,j,k,np,N,M,nt,t,kk,XX,QQ,no_of_sec,file_no,dum1,npt
      integer nprint,Allow_erosion,salin_flow,depo,d1d2
      Real  dtmax,dt_inmax,beta,c_initial_max,dtgmax,tt,dum2,Limiter,pi
      Allocatable x(:),y(:),z(:),xo(:,:),yo(:,:),zo(:,:),zo_old(:,:)
      Allocatable xco(:,:),yco(:,:),zco(:,:),xc(:,:),yc(:,:),zc(:,:)
      Allocatable ds(:,:,:),sphi(:,:,:),cphi(:,:,:),CellV(:,:)
      Allocatable U(:,:),V(:,:),H(:,:),C(:,:),zp(:,:),U_old(:,:)
      Allocatable V_old(:,:),C_old(:,:),H_old(:,:)
      Allocatable Vmag(:,:),Vmag_old(:,:),Zeta(:,:),zp_old(:,:)
      Allocatable Upre(:,:),Vpre(:,:),Hpre(:,:),Cpre(:,:),Vmag_pre(:,:)
      Allocatable xcgi(:,:),ycgi(:,:),xcgo(:,:),ycgo(:,:)      !coordinates of the gohst cells for the input and output BCs
      Allocatable xcgw(:,:),ycgw(:,:)                          !coordinates of the gohst cells for the wall BCs
      Allocatable Hgi(:,:),Ugi(:,:),Vgi(:,:),Cgi(:,:),H_inlet(:,:)          !input boundary condtions
      Allocatable Hgo(:,:),Ugo(:,:),Vgo(:,:),Cgo(:,:)          !output boundary condtions
      Allocatable Hgw(:,:),Ugw(:,:),Vgw(:,:),Cgw(:,:)          !Wall boundary condtions
      Allocatable dxi(:,:,:),deta(:,:,:)
      Allocatable dt_in(:,:),c_initail(:,:),dtg(:,:)
      Allocatable USG1_H(:,:),DSG1_H(:,:) ! for the flux limiter "the U/S gradiant & the D/S Gradiant for H IN ZAI DIR.
      Allocatable USG1_U(:,:),DSG1_U(:,:) ! for the flux limiter "the U/S gradiant & the D/S Gradiant for U IN ZAI DIR.
      Allocatable USG1_V(:,:),DSG1_V(:,:) ! for the flux limiter "the U/S gradiant & the D/S Gradiant for V IN ZAI DIR.
      Allocatable USG1_C(:,:),DSG1_C(:,:) ! for the flux limiter "the U/S gradiant & the D/S Gradiant for C IN ZAI DIR.
      Allocatable USG2_H(:,:),DSG2_H(:,:) ! for the flux limiter "the U/S gradiant & the D/S Gradiant for H IN ETA DIR.
      Allocatable USG2_U(:,:),DSG2_U(:,:) ! for the flux limiter "the U/S gradiant & the D/S Gradiant for U IN ETA DIR.
      Allocatable USG2_V(:,:),DSG2_V(:,:) ! for the flux limiter "the U/S gradiant & the D/S Gradiant for V IN ETA DIR.
      Allocatable USG2_C(:,:),DSG2_C(:,:) ! for the flux limiter "the U/S gradiant & the D/S Gradiant for C IN ETA DIR.
      Allocatable Lim_H(:,:,:),Lim_U(:,:,:),Lim_V(:,:,:),Lim_C(:,:,:)
      Allocatable hl_zai(:,:,:),hr_zai(:,:,:),hl_eta(:,:,:)
      Allocatable hr_eta(:,:,:),ul_zai(:,:,:),ur_zai(:,:,:)
      Allocatable ul_eta(:,:,:),ur_eta(:,:,:),vl_zai(:,:,:)
      Allocatable vr_zai(:,:,:),vl_eta(:,:,:),vr_eta(:,:,:)
      Allocatable cl_zai(:,:,:),cr_zai(:,:,:),cl_eta(:,:,:)
      Allocatable cr_eta(:,:,:)
      Allocatable H_hat(:,:,:),U_hat(:,:,:),V_hat(:,:,:),C_hat(:,:,:)           ! The Roe average values
      Allocatable a_hat(:,:,:),U_hat_perp(:,:,:)                                ! The Roe average values
      Allocatable Eign(:,:,:,:)           !Eign values for Roe solver
      Allocatable a_hat_L(:,:,:),a_hat_R(:,:,:)
      Allocatable U_hat_perp_L(:,:,:),U_hat_perp_R(:,:,:)
      Allocatable Eign1_L(:,:,:),Eign3_L(:,:,:),Eign1_R(:,:,:)                  ! Left and right Eign valus for Entropy fix
      Allocatable Eign3_R(:,:,:),dEign1(:,:,:),dEign3(:,:,:)
      Allocatable R_Eign1(:,:,:,:),R_Eign2(:,:,:,:),R_Eign3(:,:,:,:)            !Right Eign vectors
      Allocatable R_Eign4(:,:,:,:)
      Allocatable dh(:,:,:),du(:,:,:),dv(:,:,:),dc(:,:,:),du_perp(:,:,:)        ! Wave Strength Calculations
      Allocatable du_par(:,:,:),Ws(:,:,:,:)                                     ! Wave Strength
      Allocatable df(:,:,:,:)! ddd is the multiplication of the eign matrix by the wave strength one
      Allocatable ul_perp(:,:,:),ur_perp(:,:,:)
      Allocatable f1_R(:,:,:),f1_L(:,:,:),f2_L(:,:,:),f2_R(:,:,:)
      Allocatable f3_R(:,:,:),f3_L(:,:,:),f4_L(:,:,:),f4_R(:,:,:)                 !Flux calcualtion in Roe solver
      Allocatable dh_corec2(:,:,:),dh_corec3(:,:,:)
      Allocatable Flux1(:,:,:),Flux2(:,:,:),Flux3(:,:,:)
      Allocatable Flux4(:,:,:)
      Allocatable cb(:,:),cd(:,:),Ustar(:,:),Vstar(:,:),E(:,:)
      Allocatable cb_new(:,:),cd_new(:,:),Ustar_new(:,:)
      Allocatable Vstar_new(:,:),E_new(:,:)
      Allocatable q1(:,:),q2(:,:),q3(:,:),q4(:,:)                           ! Source terms
      Allocatable q1_new(:,:),q2_new(:,:),q3_new(:,:),q4_new(:,:)
      Allocatable dz_x(:,:),dz_y(:,:),dz_x_old(:,:),dz_y_old(:,:)
      Allocatable W1(:,:),W2(:,:),W3(:,:),W4(:,:)
      Allocatable W1_old(:,:),W2_old(:,:),W3_old(:,:),W4_old(:,:)
      Allocatable c_no(:,:),c_no_inlet(:,:)
      Allocatable dxcell(:,:),dycell(:,:)
      Allocatable xcgiw(:,:),ycgiw(:,:),Hgiw(:,:,:),cgiw(:,:,:)                 !new ghost cells at the domian corners
      Allocatable Ugiw(:,:,:),Vgiw(:,:,:)
      Allocatable Hc(:,:),Uc(:,:),Vc(:,:),Cc(:,:),Fr(:,:)                       !the new values at the original grid points
      Allocatable Hc_intial(:,:),Uc_intial(:,:),Vc_intial(:,:)
      Allocatable Cc_intial(:,:),Fr_intial(:,:)
      Allocatable theta11(:,:),theta12(:,:),theta21(:,:),theta22(:,:)
      Real Hinlet,Houtlet,Uinlet,Vinlet,Cinlet,Lim_H,Lim_U,Lim_V,Lim_C
      Real epso,D50,R,neu,ks,cb_method,constant_inflow_depth_option
      Real dxdxi,dxdeta,dydxi,dydeta,Uxi,Ueta                               ! Data reconstruction
      Real Roh,source_option,Output_option,output_file_formate
      Real Rep,Rf,Vfall,C_no_max,C_no_domain,C_no_inlet_max,widthratio
      Real intial_cs_option,L_WET,hieght_of_water1,hieght_of_water2
      real Concenteration
      Real BC1,BC2,BC3,BC4,Boundary_condition_option,addsource,order
      Real slope_factor,H_correction,lamda,tstp_adj,E_method,H_cri,s
      real*4 x,y,z,CellV,cphi,sphi
      Character ZZ
      CHARACTER*13 YY
      tt=0
      pi=3.14159265359
      print *, 'Enter the grid file name'
      read *,YY
      open (unit=1,file=YY)
      read (1,*) ZZ
      read (1,*) N, M
      np=N*M
      kk=4
      Allocate (x(np), STAT=keep)
      Allocate (y(np), STAT=keep)
      Allocate (z(np), STAT=keep)
      read (1,*) x
      read (1,*) z
      read (1,*) y
      close (1)
      i=1
      print *,'Enter the slope factor'
      read *, slope_factor
      allocate (xco(M,N), STAT=keep)
      allocate (yco(M,N), STAT=keep)
      allocate (zco(M,N), STAT=keep)
      do 1 j=1,M
      do 2 k=1,N
      xco(j,k)=x(i)
      yco(j,k)=y(i)
      zco(j,k)=z(i)*slope_factor
      i=i+1
   2  continue
   1  continue 
      allocate (xo(M,N), STAT=keep)
      allocate (yo(M,N), STAT=keep)
      allocate (zo(M,N),zo_old(M,N), STAT=keep)
      do 115 j=1,M
      do 116 k=1,N
      xo(j,k)=xco(j,N-(k-1))
      yo(j,k)=yco(j,N-(k-1))
      zo(j,k)=zco(j,N-(k-1))
      zo_old(j,k)=zo(j,k)
 116  continue
 115  continue
      deallocate(x,y,z,xco,yco,zco)
!     compute the cells center
      MM=M-1
      NN=N-1
      allocate (xc(MM,NN), STAT=keep)
      allocate (yc(MM,NN), STAT=keep)
      allocate (zc(MM,NN), STAT=keep)
      do 3 j=1,MM
      do 4 k=1,NN
      xc(j,k)=0.250*(xo(j,k)+xo(j+1,k)+xo(j,k+1)+xo(j+1,k+1))
      yc(j,k)=0.250*(yo(j,k)+yo(j+1,k)+yo(j,k+1)+yo(j+1,k+1))
      zc(j,k)=0.250*(zo(j,k)+zo(j+1,k)+zo(j,k+1)+zo(j+1,k+1))                       % the original bed elevation
   4  continue
   3  continue   
!     compute the length of the cell faces
      allocate (ds(MM,NN,4), STAT=keep)
      do 5 j=1,MM
      do 6 k=1,NN
      ds(j,k,1)=((xo(j+1,k)-xo(j,k))**2+(yo(j+1,k)-yo(j,k))**2)**0.5
      ds(j,k,2)=((xo(j+1,k+1)-xo(j+1,k))**2+ 
     &          (yo(j+1,k+1)-yo(j+1,k))**2)**0.5
      ds(j,k,3)=((xo(j,k+1)-xo(j+1,k+1))**2+ 
     &          (yo(j,k+1)-yo(j+1,k+1))**2)**0.5
      ds(j,k,4)=((xo(j,k)-xo(j,k+1))**2+(yo(j,k)-yo(j,k+1))**2)**0.5
   6  continue
   5  continue
 !    compute the cell volume
      allocate (CellV(MM,NN), STAT=keep)
      do 111 j=1,MM
      do 112 k=1,NN
      CellV(j,k)=0.50*((xo(j,k)-xo(j+1,k+1))*(yo(j+1,k)-yo(j,k+1))-
     &(xo(j+1,k)-xo(j,k+1))*(yo(j,k)-yo(j+1,k+1)))
 112  continue
 111  continue
 !    compute the angle of the cell faces with the the normal direction on the cell face
      allocate (cphi(MM,NN,4), STAT=keep)
      allocate (sphi(MM,NN,4), STAT=keep)
      do 7 j=1,MM
      do 8 k=1,NN
      cphi(j,k,1)=-(yo(j+1,k)-yo(j,k))/ds(j,k,1)                         cos phi hl face
      cphi(j,k,2)=(yo(j+1,k+1)-yo(j+1,k))/ds(j,k,2)                     cos phi vl face
      cphi(j,k,3)=-(yo(j+1,k+1)-yo(j,k+1))/ds(j,k,3)
      cphi(j,k,4)=(yo(j,k+1)-yo(j,k))/ds(j,k,4)
      sphi(j,k,1)=(xo(j+1,k)-xo(j,k))/ds(j,k,1)                         sin phi hl face
      sphi(j,k,2)=-(xo(j+1,k+1)-xo(j+1,k))/ds(j,k,2)                     sin phi vl face
      sphi(j,k,3)=(xo(j+1,k+1)-xo(j,k+1))/ds(j,k,3)
      sphi(j,k,4)=-(xo(j,k+1)-xo(j,k))/ds(j,k,4)
   8  continue
   7  continue
      file_no=417
      open (unit=file_no,file="cphi.txt",action="write",
     &      status="unknown")
      do 418 k=1,NN
      write (file_no,*) k
      do 419 kk=1,4
      write (file_no,*) 'face no',kk
      write (file_no,*) cphi(:,k,kk)
 419  continue
 418  continue   
      close (file_no)
      file_no=418
      open (unit=file_no,file="sphi.txt",action="write",
     &      status="unknown")
      do 420 k=1,NN
      write (file_no,*) k
      do 421 kk=1,4
      write (file_no,*) 'face no',kk
      write (file_no,*) sphi(:,k,kk)
 421  continue
 420  continue   
      close (file_no)    
******************************************************************************************************************************   
   !  Boundary condition
   !    initialization of gohst cells   "the BCs will be located at this points" "needs more work for complex geometery"
      allocate (xcgi(1,NN), STAT=keep)
      allocate (ycgi(1,NN), STAT=keep)
      allocate (xcgo(1,NN), STAT=keep)
      allocate (ycgo(1,NN), STAT=keep)
      do 200 k=1,NN                             !input and output gohst cells
      xcgi(1,k)=xc(1,k)-ds(1,k,1)
      ycgi(1,k)=yc(1,k)
      xcgo(1,k)=xc(MM,k)+ds(MM,k,1)
      ycgo(1,k)=yc(MM,k)
  200 continue 
      allocate (xcgw(MM,2), STAT=keep)
      allocate (ycgw(MM,2), STAT=keep)
      do 201 j=1,MM                                                             !Must be corrected
      ycgw(j,1)=yc(j,1)-(yc(j,2)-yc(j,1))
      xcgw(j,1)=xc(j,1)-(xc(j,2)-xc(j,1))
      ycgw(j,2)=yc(j,NN)-(yc(j,NN-1)-yc(j,NN))
      xcgw(j,2)=xc(j,NN)-(xc(j,NN-1)-xc(j,NN))
  201 continue
  !   BC s for supercritical " intially this model will solve a dam break problem" 30 April 2010 "This has been updated"
   !  inflow BC
******************************************************************************************************************
!     Reading the input file data
      print *,'read the input data from the input file  "(1)Y or (2)N"'
      read *,QQ
      if (QQ.EQ.1) then
      open (unit=2,file='Input.dat')
      read (2,*) ZZ
      read (2,*) d1d2
      read (2,*)ZZ
      read (2,*) roh
      read (2,*) ZZ
      read (2,*) nt
      read (2,*) ZZ,ZZ,ZZ,ZZ
      read (2,*) Hinlet,Uinlet,Vinlet,Cinlet
      read (2,*) ZZ
      read (2,*) constant_inflow_depth_option
      read (2,*) ZZ,ZZ,ZZ
      read (2,*) D50,R,lamda
      read (2,*) ZZ,ZZ
      read (2,*) neu,ks
      read (2,*) ZZ
      read (2,*) beta
      read (2,*) ZZ
      read (2,*) cb_method
      read (2,*) ZZ
      read (2,*)E_method
      read (2,*) ZZ
      read (2,*) source_option
      read (2,*) ZZ
      read (2,*) Output_option
      read (2,*) ZZ
      read (2,*) epso
      read (2,*) ZZ
      read (2,*) no_of_sec
      read (2,*) ZZ
      read (2,*) output_file_formate
      read (2,*) ZZ
      read (2,*) intial_cs_option
      read (2,*) ZZ
      read (2,*) Boundary_condition_option
      read (2,*) ZZ
      read (2,*) Houtlet
      read (2,*) ZZ
      read (2,*) nprint
      read (2,*) ZZ
      read (2,*) addsource
      read (2,*) ZZ
      read (2,*) H_correction
      read (2,*) ZZ
      read (2,*) order
      read (2,*) ZZ
      read (2,*) Allow_erosion
      read (2,*) ZZ
      read (2,*) salin_flow
      read (2,*) ZZ
      read (2,*) depo
      close (2)
      else                                      !needs to be like the input file it does not has all the options
      print *,'Enter the inflow conditions'
      print *,'H inlet (m)'
      read *,Hinlet
      print *,'U inlet (m/sec)'
      read *,Uinlet
      print *,'V inlet (m/sec)'
      read *,Vinlet
      print *,'C inlet (%)'
      read *,Cinlet
      print *,'The total time steps'
      read *,nt
      endif
      if (no_of_sec.GT.NN) then
      no_of_sec=NN
      endif
      Cinlet=Cinlet/100
      print *,'Enter the maximum courant number"Less than 1"'
      read *, tstp_adj

********************************************************************************************************************      
      file_no=1200
      open (unit=file_no,file="CellV.txt",action="write",
     &      status="unknown")
      write (file_no,*) MM,no_of_sec
      do 5000 k=1,no_of_sec
      write (file_no,*) k
      write (file_no,*) CellV(:,k)  
5000   continue         
      close (file_no)
      file_no=1201
      open (unit=file_no,file="ds.txt",action="write",
     &      status="unknown")
      write (file_no,*) MM,no_of_sec
      do 5002 k=1,no_of_sec
      write (file_no,*) k
      do 5001 kk=1,4
      write (file_no,*) 'Face number)',kk
      write (file_no,*) ds(:,k,kk) 
5001  continue       
5002  continue         
      close (file_no)
**********************************************************************************************************************
!     Intial conditions
!     " the in the code intial conditions are for the case of dam break problem unless you chosse read from input file" 
      Allocate (U(MM,NN),U_old(MM,NN),STAT=keep)
      Allocate (V(MM,NN),V_old(MM,NN), STAT=keep)
      Allocate (H(MM,NN),H_old(MM,NN), STAT=keep)
      Allocate (C(MM,NN),C_old(MM,NN), STAT=keep)
      Allocate (zp(MM,NN),Zeta(MM,NN),zp_old(MM,NN),STAT=keep)
      Allocate (Vmag(MM,NN),Vmag_old(MM,NN), STAT=Keep)
      Allocate (W1(MM,NN),W2(MM,NN),W3(MM,NN),W4(MM,NN),STAT=keep)
      Allocate (W1_old(MM,NN),W2_old(MM,NN),STAT=keep)
      Allocate (W3_old(MM,NN),W4_old(MM,NN),STAT=keep)
      zp=zc
      if (intial_cs_option.EQ.1) then
      open (unit=3,file='intial&BC_conditions.dat')
      read (3,*) ZZ
      read (3,*) L_WET
      read (3,*) ZZ
      read (3,*) hieght_of_water1
      read (3,*) zz
      read (3,*) hieght_of_water2
      read (3,*) ZZ
      read (3,*) Concenteration
      read (3,*) ZZ
      read (3,*) BC1
      read (3,*) ZZ
      read (3,*) BC2
      read (3,*) ZZ
      read (3,*) BC3
      read (3,*) ZZ
      read (3,*) BC4
      close (3)
      do 190 j=1,M
      dum2=xo(j,1)
      dum1=j                                                                    !no of wet points in X dirction
      if (dum2.GE.L_WET) then
      go to 1234
      endif
  190  continue
 1234 print *,dum1     
      do 191 j=1,MM
      do 192 k=1,NN
      if (j.LE.dum1) then
      H(j,k)=hieght_of_water1
      C(j,k)=Concenteration/100
      else
      H(j,k)=hieght_of_water2
      C(j,k)=Concenteration/100
      endif
      U(j,k)=0
      V(j,k)=0
      Vmag(j,k)=(U(j,k))**2+(V(j,k))**2
      W1(j,k)=H(j,k)
      W2(j,k)=H(j,k)*U(j,k)
      W3(j,k)=H(j,k)*V(j,k)
      W4(j,k)=H(j,k)*C(j,k)
      zeta(j,k)=zp(j,k)+H(j,k)
 192  continue
 191  continue          
      else
      open (unit=3,file='intial&BC_conditions.dat')
      read (3,*) ZZ
      read (3,*) L_WET
      read (3,*) ZZ
      read (3,*) hieght_of_water1
      read (3,*) zz
      read (3,*) hieght_of_water2
      read (3,*) ZZ
      read (3,*) Concenteration
      read (3,*) ZZ
      read (3,*) BC1
      read (3,*) ZZ
      read (3,*) BC2
      read (3,*) ZZ
      read (3,*) BC3
      read (3,*) ZZ
      read (3,*) BC4
      close (3)
      do 9 j=1,MM                                                           Dry bed as intial condition
      do 10 k=1,NN                    
      U(j,k)=0
      V(j,k)=0
      H(j,k)=0
      C(j,k)=0
      Vmag(j,k)=(U(j,k))**2+(V(j,k))**2
      W1(j,k)=H(j,k)
      W2(j,k)=H(j,k)*U(j,k)
      W3(j,k)=H(j,k)*V(j,k)
      W4(j,k)=H(j,k)*C(j,k)
      zeta(j,k)=zp(j,k)+H(j,k)
   10 continue
   9  continue
      endif
*************************************************************************************************************************     
!     Boundary conditions      
      allocate (Hgi(1,NN), STAT=keep)
      allocate (Ugi(1,NN), STAT=keep)
      allocate (Vgi(1,NN), STAT=keep)
      allocate (Cgi(1,NN), STAT=keep)
      allocate (Hgo(1,NN), STAT=keep)
      allocate (Ugo(1,NN), STAT=keep)
      allocate (Vgo(1,NN), STAT=keep)
      allocate (Cgo(1,NN), STAT=keep)
      allocate (Hgw(MM,2), STAT=keep)
      allocate (Ugw(MM,2), STAT=keep)
      allocate (Vgw(MM,2), STAT=keep)
      allocate (Cgw(MM,2), STAT=keep)
      Allocate (H_inlet(1,NN), STAT=keep)
      t=1
      if (H(MM,1).EQ.0.0) then
      H_cri=0.9*((Hinlet*Uinlet)**2.0/9.81)**0.333333333
      else
      H_cri=0.9*((H(MM,1)*U(MM,1))**2.0/9.81)**0.333333333
      endif
 1500 if (Boundary_condition_option.EQ.1) then
      ! BC for the wall no.1 "inlet"
      open (unit=4,file='inflow2.txt')
      read (4,*) ZZ
      do 2001 k=1,NN
      if (constant_inflow_depth_option .EQ.2) then
      read (4,*) H_inlet(1,k)
      endif
      H_inlet(1,k)=H_inlet(1,k)/100
      if (constant_inflow_depth_option .EQ.2) then
      Hgi(1,k)=Bound_h(H_inlet(1,k),H(1,k),H(MM,k),H_cri,Houtlet,BC1)
      else
      Hgi(1,k)=Bound_h(Hinlet,H(1,k),H(MM,k),H_cri,Houtlet,BC1)
      endif
      Cgi(1,k)=Bound_c(Cinlet,C(1,k),C(MM,k),BC1)
      Ugi(1,k)=Bound_U(Uinlet,U(1,k),V(1,k),cphi(1,k,4),
     &          sphi(1,k,4),BC1)
      Vgi(1,k)=Bound_V(Vinlet,U(1,k),V(1,k),cphi(1,k,4),
     &          sphi(1,k,4),BC1)
 2001 continue
      close (4)
      ! BC for the wall no.2 "outlet"
      do 2002 k=1,NN
      Hgo(1,k)=Bound_h(Hinlet,H(MM,k),H(MM,k),H_cri,Houtlet,BC2)
      Cgo(1,k)=Bound_c(Cinlet,C(MM,k),C(MM,k),BC2)
      Ugo(1,k)=Bound_U(Uinlet,U(MM,k),V(MM,k),cphi(MM,k,2),
     &          sphi(MM,k,2),BC2)
      Vgo(1,k)=Bound_V(Vinlet,U(MM,k),V(MM,k),cphi(MM,k,2),
     &          sphi(MM,k,2),BC2)
 2002 continue    
      ! BC for the wall no.3 "Lower wall"
      k=1
      do 2003 j=1,MM
      Hgw(j,1)=Bound_h(Hinlet,H(j,1),H(j,1),H_cri,Houtlet,BC3)
      Cgw(j,1)=Bound_c(Cinlet,C(1,k),C(j,1),BC3)
      Ugw(j,1)=Bound_U(Uinlet,U(j,1),V(j,1),cphi(j,1,1),
     &          sphi(j,1,1),BC3)
      Vgw(j,1)=Bound_V(Vinlet,U(j,1),V(j,1),cphi(j,1,1),
     &          sphi(j,1,1),BC3)
      ! BC for the wall no.4 "Top wall"
      Hgw(j,2)=Bound_h(Hinlet,H(j,NN),H(j,NN),H_cri,Houtlet,BC4)
      Cgw(j,2)=Bound_c(Cinlet,C(1,k),C(j,NN),BC3)
      Ugw(j,2)=Bound_U(Uinlet,U(j,NN),V(j,NN),cphi(j,NN,3),
     &          sphi(j,NN,3),BC4)
      Vgw(j,2)=Bound_V(Vinlet,U(j,NN),V(j,NN),cphi(j,NN,3),
     &          sphi(j,NN,3),BC4)
      
2003  continue
      if (t.GT.1) then
      goto 1501
      endif
      else
      ! BC for the wall no.1
      do 35 t=1,nt
      do 11 k=1,NN
      Hgi(1,k)=Hinlet
      Ugi(1,k)=Uinlet
      Vgi(1,k)=Vinlet
      Cgi(1,k)=Cinlet
   11 continue
   35 continue
      ! BC for the wall no.2
   !  outflow BC
      do 12 k=1,NN
      if (BC2.EQ.3) then
      Hgo(1,k)=H(MM,k)
      elseif(H(MM,k).LE.H_cri.AND.BC.EQ.4) then
      Hgo(1,k)=H(MM,k)
      else
      Hgo(1,k)=Houtlet
      endif
      Ugo(1,k)=U(MM,k)
      Vgo(1,k)=V(MM,k)
      Cgo(1,k)=C(MM,k)
   12 continue
      ! BC for the wall no.3
      
   !  Wall BCs          the net perpendicular velocity to the wall is zero
      do 14 j=1,MM                      
      Hgw(j,1)=H(j,1)                                       !For the lower wall
      Cgw(j,1)=C(j,1)                                       !For the lower wall
      Ugw(j,1)=U(j,1)*((sphi(j,1,1))**2-(cphi(j,1,1))**2)-  !For the lower wall
     &          2.0*V(j,1)*sphi(j,1,1)*cphi(j,1,1)
      Vgw(j,1)=-2.0*U(j,1)*sphi(j,1,1)*cphi(j,1,1)+         !For the lower wall
     &              +V(j,1)*((cphi(j,1,1))**2-(sphi(j,1,1))**2)
      ! BC for the wall no.4
      Hgw(j,2)=H(j,NN)                                       !For the Upper wall
      Cgw(j,2)=C(j,NN)                                       !For the Upper wall
      Ugw(j,2)=U(j,NN)*((sphi(j,NN,3))**2-(cphi(j,NN,3))**2)-  !For the Upper wall
     &          2.0*V(j,NN)*sphi(j,NN,3)*cphi(j,NN,3)
      Vgw(j,2)=-2.0*U(j,NN)*sphi(j,NN,3)*cphi(j,NN,3)+         !For the Upper wall
     &              +V(j,NN)*((cphi(j,NN,3))**2-(sphi(j,NN,3))**2)
   14 continue
      endif
!*************************************************************************************************************
!     time step calculation for the first time step
      Allocate (dtg(1,NN),c_initail(1,NN),dt_in(MM,NN),STAT=Keep)
      do 125 k=1,NN
      dtg=((Uinlet+sqrt(9.81*Hinlet))/(abs(xo(1,1)-xo(2,1)))+
     &(Vinlet+sqrt(9.81*Hinlet))/(abs(yo(1,1)-yo(1,2))))**(-1)
      c_initail(1,k)=(Uinlet+sqrt(9.81*Hinlet))/
     &(0.5*(abs(xo(1+1,k)-xo(1,k))+abs(xo(1,k+1)-xo(1+1,k+1))))
     &             +(Vinlet+sqrt(9.81*Hinlet))/
     &         (0.5*(abs(yo(1+1,k+1)-yo(1+1,k))+abs(yo(1,k)-yo(1,k+1))))
  125 continue   
      do 40 j=1,MM
      do 41 k=1,NN
      dt_in(j,k)=((U(j,k)+(9.81*H(j,k))**0.5)/
     &      (abs(xo(1,1)-xo(2,1)))+(V(j,k)+
     &    (9.81*H(j,k))**0.5)/(abs(yo(1,1)-yo(1,2))))**(-1)
  41  continue
  40  continue
      dt_inmax=MINVAL(dt_in)
      dtgmax=MINVAL(dtg)
      dtmax=min(dt_inmax,dtgmax)
  !    print *,'The time step for the first time step is=',dtmax,'sec'
!*************************************************************************************************************
   !  Predictor
   !  Data reconstruction
      allocate (dxi(MM,NN,2),deta(MM,NN,2),STAT=keep)
      file_no=505
      open (unit=file_no,file="area.txt",action="write",
     &      status="unknown")
      do 32 j=1,MM
      do 33 k=1,NN
      dxdxi=0.5*(-xo(j,k)+xo(j+1,k)+xo(j+1,k+1)-xo(j,k+1))
      dxdeta=0.5*(-xo(j,k)-xo(j+1,k)+xo(j+1,k+1)+xo(j,k+1))
      dydxi=0.5*(-yo(j,k)+yo(j+1,k)+yo(j+1,k+1)-yo(j,k+1))
      dydeta=0.5*(-yo(j,k)-yo(j+1,k)+yo(j+1,k+1)+yo(j,k+1))
      area=dxdxi*dydeta-dxeta*dydxi
      write (file_no,*) k
      write (file_no,*) area
      dxi(j,k,1)=dydeta/area
      dxi(j,k,2)=-dxdeta/area
      deta(j,k,1)=-dydxi/area
      deta(j,k,2)=dxdxi/area
  33  continue
  32  continue
      close (file_no)
  !**************************************
  ! output for dxi and deta
       file_no=501
      open (unit=file_no,file="dxi1.txt",action="write",
     &      status="unknown")
      do 301 k=1,NN
      write (file_no,*) k
      write (file_no,*) dxi(:,k,1)
 301  continue   
      close (file_no)   
      file_no=502
      open (unit=file_no,file="dxi2.txt",action="write",
     &      status="unknown")
      do 302 k=1,NN
      write (file_no,*) k
      write (file_no,*) dxi(:,k,2)
 302  continue   
      close (file_no) 
      file_no=503
      open (unit=file_no,file="deta1.txt",action="write",
     &      status="unknown")
      do 303 k=1,NN
      write (file_no,*) k
      write (file_no,*) deta(:,k,1)
 303  continue   
      close (file_no) 
      file_no=504
      open (unit=file_no,file="deta2.txt",action="write",
     &      status="unknown")
      do 304 k=1,NN
      write (file_no,*) k
      write (file_no,*) deta(:,k,2)
 304  continue   
      close (file_no)
       
*==========================================================================================================  
*     Allocate Memory for each parameter in the solution
      allocate (Upre(MM,NN),Vpre(MM,NN),Hpre(MM,NN), STAT=keep)
      allocate (Vmag_pre(MM,NN),Cpre(MM,NN), STAT=keep)
      Allocate (USG1_H(MM,NN),DSG1_H(MM,NN), STAT=Keep)
      Allocate (USG1_U(MM,NN),DSG1_U(MM,NN), STAT=Keep)
      Allocate (USG1_V(MM,NN),DSG1_V(MM,NN), STAT=Keep)
      Allocate (USG1_C(MM,NN),DSG1_C(MM,NN), STAT=Keep)
      Allocate (USG2_H(MM,NN),DSG2_H(MM,NN), STAT=Keep)
      Allocate (USG2_U(MM,NN),DSG2_U(MM,NN), STAT=Keep)
      Allocate (USG2_V(MM,NN),DSG2_V(MM,NN), STAT=Keep)
      Allocate (USG2_C(MM,NN),DSG2_C(MM,NN), STAT=Keep)
      Allocate (Lim_H(MM,NN,2),Lim_U(MM,NN,2),Lim_V(MM,NN,2), STAT=keep)
      Allocate (Lim_C(MM,NN,2), STAT=keep)
      Allocate (hl_zai(MM,NN,4),hr_zai(MM,NN,4),STAT=keep)
      Allocate (hl_eta(MM,NN,4),hr_eta(MM,NN,4),STAT=keep)
      Allocate (ul_zai(MM,NN,4),ur_zai(MM,NN,4),STAT=keep)
      Allocate (ul_eta(MM,NN,4),ur_eta(MM,NN,4),STAT=kepp)
      Allocate (vl_zai(MM,NN,4),vr_zai(MM,NN,4),STAT=kepp)
      Allocate (vl_eta(MM,NN,4),vr_eta(MM,NN,4),STAT=kepp)
      Allocate (cl_zai(MM,NN,4),cr_zai(MM,NN,4),STAT=kepp)
      Allocate (cl_eta(MM,NN,4),cr_eta(MM,NN,4),STAT=kepp)
      Allocate (H_hat(MM,NN,4),U_hat(MM,NN,4),V_hat(MM,NN,4),STAT=Keep)
      Allocate (C_hat(MM,NN,4),a_hat(MM,NN,4),STAT=Keep)
      Allocate (U_hat_perp(MM,NN,4),STAT=Keep)
      Allocate (Eign(MM,NN,4,4),STAT=keep)
      Allocate (Eign1_L(MM,NN,4),Eign3_L(MM,NN,4),STAT=keep)
      Allocate (Eign1_R(MM,NN,4),Eign3_R(MM,NN,4),STAT=keep)
      Allocate (dEign1(MM,NN,4),dEign3(MM,NN,4),STAT=keep)
      Allocate (U_hat_perp_L(MM,NN,4),U_hat_perp_R(MM,NN,4),STAT=keep)
      Allocate (a_hat_L(MM,NN,4),a_hat_R(MM,NN,4),STAT=Keep)
      Allocate (R_Eign1(MM,NN,4,4),R_Eign2(MM,NN,4,4),STAT=keep)
      Allocate (R_Eign3(MM,NN,4,4),R_Eign4(MM,NN,4,4),STAT=keep)
      Allocate (dh(MM,NN,4),du(MM,NN,4),dv(MM,NN,4),STAT=Keep)
      Allocate (dc(MM,NN,4),du_perp(MM,NN,4),du_par(MM,NN,4),STAT=Keep)
      Allocate (Ws(MM,NN,4,4),df(MM,NN,4,4),STAT=Keep)
      Allocate (ul_perp(MM,NN,4),ur_perp(MM,NN,4),STAT=Keep)
      Allocate (f1_L(MM,NN,4),f1_R(MM,NN,4),f2_L(MM,NN,4),STAT=Keep)
      Allocate (f2_R(MM,NN,4),f3_L(MM,NN,4),f3_R(MM,NN,4),STAT=Keep)
      Allocate (f4_L(MM,NN,4),f4_R(MM,NN,4),Flux1(MM,NN,4),STAT=Keep)
      Allocate (dh_corec2(MM,NN,4),dh_corec3(MM,NN,4),STAT=Keep)
      Allocate (Flux2(MM,NN,4),Flux3(MM,NN,4),STAT=Keep)
      Allocate (Flux4(MM,NN,4),STAT=Keep)
      Allocate (cd(MM,NN),cb(MM,NN),STAT=Keep)
      Allocate (cb_new(MM,NN),STAT=Keep)
      Allocate (Ustar(MM,NN),Vstar(MM,NN),E(MM,NN),STAT=Keep)
      Allocate (q1(MM,NN),q2(MM,NN),q3(MM,NN),STAT=Keep)                           ! Source terms
      Allocate (q4(MM,NN),dz_x(MM,NN),dz_y(MM,NN),STAT=Keep)
      Allocate (dz_x_old(MM,NN),dz_y_old(MM,NN), STAT=Keep)
      Allocate (Ustar_new(MM,NN),Vstar_new(MM,NN),STAT=Keep)
      Allocate (E_new(MM,NN),q1_new(MM,NN),STAT=Keep)
      Allocate (cd_new(MM,NN),q4_new(MM,NN),STAT=Keep)
      Allocate (q2_new(MM,NN),q3_new(MM,NN),STAT=Keep)
      Allocate (c_no(MM,NN),c_no_inlet(MM,NN),STAT=Keep)
      Allocate (dxcell(MM,NN),dycell(MM,NN),STAT=Keep)
      Allocate (xcgiw(2,2),ycgiw(2,2),Hgiw(2,2,nt),STAT=Keep)
      Allocate (Ugiw(2,2,nt),Vgiw(2,2,nt),cgiw(2,2,nt),STAT=Keep)
      Allocate (Hc(M,N),Uc(M,N),Vc(M,N),Cc(M,N),STAT=Keep)
      Allocate (Fr(M,N),STAT=Keep)
      Allocate (Hc_intial(M,N),Uc_intial(M,N),Vc_intial(M,N),STAT=Keep)
      Allocate (Cc_intial(M,N),Fr_intial(M,N),STAT=Keep)
**********************************************************************************************************************
**********************************************************************************************************************
      t=1
      goto 1501
*****************************************************************************************************************************      
 !   Computation of the primitive variables
1502  do 100 t=2,nt
      H_old=H
      U_old=U
      V_old=V
      C_old=C
      Vmag_old=Vmag                                         ! The beginning of the calculations
      tt=tt+dtmax
*     THE CALCULATION FOR THE U/S AND D/S GRADIENT in x direction
      do 47 j=1,MM
      do 48 k=1,NN  
      if (j.EQ.1) then
      USG1_H(j,k)=H_old(j,k)-Hgi(1,k)
      USG1_U(j,k)=U_old(j,k)-Ugi(1,k)
      USG1_V(j,k)=V_old(j,k)-Vgi(1,k)
      USG1_C(j,k)=C_old(j,k)-Cgi(1,k)
      DSG1_H(j,k)=H_old(j+1,k)-H_old(j,k)
      DSG1_U(j,k)=U_old(j+1,k)-U_old(j,k)
      DSG1_V(j,k)=V_old(j+1,k)-V_old(j,k)
      DSG1_C(j,k)=C_old(j+1,k)-C_old(j,k)
      else if (j.EQ.MM) then
      DSG1_H(j,k)=Hgo(1,k)-H_old(j,k)
      DSG1_U(j,k)=Ugo(1,k)-U_old(j,k)
      DSG1_V(j,k)=Vgo(1,k)-V_old(j,k)
      DSG1_C(j,k)=Cgo(1,k)-C_old(j,k)
      USG1_H(j,k)=H_old(j,k)-H_old(j-1,k)
      USG1_U(j,k)=U_old(j,k)-U_old(j-1,k)
      USG1_V(j,k)=V_old(j,k)-V_old(j-1,k)
      USG1_C(j,k)=C_old(j,k)-C_old(j-1,k)
      else
      USG1_H(j,k)=H_old(j,k)-H_old(j-1,k)
      USG1_U(j,k)=U_old(j,k)-U_old(j-1,k)
      USG1_V(j,k)=V_old(j,k)-V_old(j-1,k)
      USG1_C(j,k)=C_old(j,k)-C_old(j-1,k)
      DSG1_H(j,k)=H_old(j+1,k)-H_old(j,k)
      DSG1_U(j,k)=U_old(j+1,k)-U_old(j,k)
      DSG1_V(j,k)=V_old(j+1,k)-V_old(j,k)
      DSG1_C(j,k)=C_old(j+1,k)-C_old(j,k)
      endif
  48  continue
  47  continue
*     THE CALCULATION FOR THE U/S AND D/S GRADIENT in y direction
      do 49 k=1,NN
      do 50 j=1,MM  
      if (k.EQ.1) then
      USG2_H(j,k)=H_old(j,k)-Hgw(j,1)
      USG2_U(j,k)=U_old(j,k)-Ugw(j,1)
      USG2_V(j,k)=V_old(j,k)-Vgw(j,1)
      USG2_C(j,k)=C_old(j,k)-Cgw(j,1)
      DSG2_H(j,k)=H_old(j,k+1)-H_old(j,k)
      DSG2_U(j,k)=U_old(j,k+1)-U_old(j,k)
      DSG2_V(j,k)=V_old(j,k+1)-V_old(j,k)
      DSG2_C(j,k)=C_old(j,k+1)-C_old(j,k) 
      else if (k.EQ.NN) then
      DSG2_H(j,k)=Hgw(j,2)-H_old(j,k)
      DSG2_U(j,k)=Ugw(j,2)-U_old(j,k)
      DSG2_V(j,k)=Vgw(j,2)-V_old(j,k)
      DSG2_C(j,k)=Cgw(j,2)-C_old(j,k)
      USG2_H(j,k)=H_old(j,k)-H_old(j,k-1)
      USG2_U(j,k)=U_old(j,k)-U_old(j,k-1)
      USG2_V(j,k)=V_old(j,k)-V_old(j,k-1)
      USG2_C(j,k)=C_old(j,k)-C_old(j,k-1)
      else
      USG2_H(j,k)=H_old(j,k)-H_old(j,k-1)
      USG2_U(j,k)=U_old(j,k)-U_old(j,k-1)
      USG2_V(j,k)=V_old(j,k)-V_old(j,k-1)
      USG2_C(j,k)=C_old(j,k)-C_old(j,k-1)
      DSG2_H(j,k)=H_old(j,k+1)-H_old(j,k)
      DSG2_U(j,k)=U_old(j,k+1)-U_old(j,k)
      DSG2_V(j,k)=V_old(j,k+1)-V_old(j,k)
      DSG2_C(j,k)=C_old(j,k+1)-C_old(j,k) 
      endif
  50  continue
  49  continue   
*     CELL AVERAGE GRADIENT ACROSS THE CELL FACE USING THE FLUX LIMITER 
      do 51 j=1,MM
      do 52 k=1,NN
      Lim_H(j,k,1)=homar(order,beta,USG1_H(j,k),DSG1_H(j,k))
      Lim_U(j,k,1)=homar(order,beta,USG1_U(j,k),DSG1_U(j,k))
      Lim_V(j,k,1)=homar(order,beta,USG1_V(j,k),DSG1_V(j,k))
      Lim_C(j,k,1)=homar(order,beta,USG1_C(j,k),DSG1_C(j,k))
      Lim_H(j,k,2)=homar(order,beta,USG2_H(j,k),DSG2_H(j,k))
      Lim_U(j,k,2)=homar(order,beta,USG2_U(j,k),DSG2_U(j,k))
      Lim_V(j,k,2)=homar(order,beta,USG2_V(j,k),DSG2_V(j,k))
      Lim_C(j,k,2)=homar(order,beta,USG2_C(j,k),DSG2_C(j,k))
   52 continue
   51 continue   
*==================================================================================================================
      !=======================================================================================
*     PREDICTOR CALCULATION "Data Reconstrubtion"
*     This calculation are at time=n+1/2            see sec. 3.4 in shock caputering methods
      do 53 j=1,MM
      do 54 k=1,NN
      Uxi=U_old(j,k)*dxi(j,k,1)+V_old(j,k)*dxi(j,k,2)
      Ueta=U_old(j,k)*deta(j,k,1)+V_old(j,k)*deta(j,k,2)
      if (addsource.EQ.1) then
      q1(j,k)=0.0                                                   
      q2(j,k)=0.0
      q3(j,k)=0.0
      q4(j,k)=0.0
      endif
      Hpre(j,k)=H_old(j,k)-0.5*dtmax*(Uxi*Lim_H(j,k,1)+Ueta*Lim_H(j,k,2)
     &          +H_old(j,k)*(dxi(j,k,1)*Lim_U(j,k,1)+
     &          dxi(j,k,2)*Lim_V(j,k,1)+deta(j,k,1)*Lim_U(j,k,2)+
     &          deta(j,k,2)*Lim_V(j,k,2)))
      Upre(j,k)=U_old(j,k)-0.5*dtmax*(Uxi*Lim_U(j,k,1)+Ueta*Lim_U(j,k,2)
     &          +9.81*(dxi(j,k,1)*Lim_h(j,k,1)+
     &          deta(j,k,1)*Lim_H(j,k,2)))
      Vpre(j,k)=V_old(j,k)-0.5*dtmax*(Uxi*Lim_V(j,k,1)+Ueta*Lim_V(j,k,2)
     &          +9.81*(dxi(j,k,2)*Lim_H(j,k,1)+
     &          deta(j,k,2)*Lim_H(j,k,2)))
      Cpre(j,k)=C_old(j,k)-0.5*dtmax*(Uxi*Lim_C(j,k,1)+
     &          Ueta*Lim_C(j,k,2))
      Vmag_pre(j,k)=SQRT(Upre(j,k)**2+Vpre(j,k)**2)
      if(Hpre(j,k).LE.0.0) then
      Hpre(j,k)=0.0
      Upre(j,k)=0.0
      Vpre(j,k)=0.0
      Cpre(j,k)=0.0
      Vmag_pre(j,k)=0.0
      endif
  54  continue
  53  continue 
*********************************************************************************************************   
*     Recalculate the bounadry condition for the new time step "n+1/2"
      if (H(MM,1).EQ.0.0) then
      H_cri=0.9*((Hinlet*Uinlet)**2.0/9.81)**0.333333333
      else
      H_cri=0.9*((Hpre(MM,1)*Upre(MM,1))**2.0/9.81)**0.333333333
      endif
      if (Boundary_condition_option.EQ.1) then
      ! BC for the wall no.1
      do 2004 k=1,NN
      if (constant_inflow_depth_option .EQ.2) then
      Hgi(1,k)=Bound_h(H_inlet(1,k),H(1,k),H(MM,k),H_cri,Houtlet,BC1)
      else
      Hgi(1,k)=Bound_h(Hinlet,H(1,k),H(MM,k),H_cri,Houtlet,BC1)
      endif
      Cgi(1,k)=Bound_C(Cinlet,Cpre(1,k),Cpre(MM,k),BC1)
      Ugi(1,k)=Bound_U(Uinlet,Upre(1,k),Vpre(1,k),cphi(1,k,4),
     &          sphi(1,k,4),BC1)
      Vgi(1,k)=Bound_V(Vinlet,Upre(1,k),Vpre(1,k),cphi(1,k,4),
     &          sphi(1,k,4),BC1)
 2004 continue
      ! BC for the wall no.2
      do 2005 k=1,NN
      Hgo(1,k)=Bound_h(Hinlet,Hpre(MM,k),Hpre(MM,k),H_cri,Houtlet,BC2)
      Cgo(1,k)=Bound_c(Cinlet,Cpre(MM,k),Cpre(MM,k),BC2)
      Ugo(1,k)=Bound_U(Uinlet,Upre(MM,k),Vpre(MM,k),
     &          cphi(MM,k,2),sphi(MM,k,2),BC2)
      Vgo(1,k)=Bound_V(Vinlet,Upre(MM,k),Vpre(MM,k),
     &          cphi(MM,k,2),sphi(MM,k,2),BC2)
 2005 continue
      ! BC for the wall no.3
      k=1
      do 2006 j=1,MM
      Hgw(j,1)=Bound_h(Hinlet,Hpre(j,1),Hpre(j,1),H_cri,Houtlet,BC3)
      Cgw(j,1)=Bound_c(Cinlet,Cpre(1,k),Cpre(j,1),BC3)
      Ugw(j,1)=Bound_U(Uinlet,Upre(j,1),Vpre(j,1),cphi(j,1,1),
     &          sphi(j,1,1),BC3)
      Vgw(j,1)=Bound_V(Vinlet,Upre(j,1),Vpre(j,1),cphi(j,1,1),
     &          sphi(j,1,1),BC3)
      ! BC for the wall no.4
      Hgw(j,2)=Bound_h(Hinlet,Hpre(j,NN),Hpre(j,NN),H_cri,Houtlet,BC4)
      Cgw(j,2)=Bound_c(Cinlet,Cpre(1,k),Cpre(j,NN),BC4)
      Ugw(j,2)=Bound_U(Uinlet,Upre(j,NN),Vpre(j,NN),
     &          cphi(j,NN,3),sphi(j,NN,3),BC4)
      Vgw(j,2)=Bound_V(Vinlet,Upre(j,NN),Vpre(j,NN),
     &          cphi(j,NN,3),sphi(j,NN,3),BC4)
2006  continue
      else
*     Outflow BCs
      do 102 k=1,NN
      if (BC2.EQ.3) then
      Hgo(1,k)=H(MM,k)
      elseif(H(MM,k).LE.H_cri.AND.BC.EQ.4) then
      Hgo(1,k)=H(MM,k)
      else
      Hgo(1,k)=Houtlet
      endif
      Ugo(1,k)=Upre(MM,k)
      Vgo(1,k)=Vpre(MM,k)
      Cgo(1,k)=Cpre(MM,k)
  102 continue
*     Wall BCs in the cell local axies  "the net perpendicular velocity to the wall is zero"
      do 104 j=1,MM                      
      Hgw(j,1)=Hpre(j,1)                                       !For the lower wall
      Cgw(j,1)=Cpre(j,1)                                       !For the lower wall
      Ugw(j,1)=Upre(j,1)*((sphi(j,1,1))**2-(cphi(j,1,1))**2)-  !For the lower wall 
     &          2.0*Vpre(j,1)*sphi(j,1,1)*cphi(j,1,1)            
      Vgw(j,1)=-2.0*Upre(j,1)*sphi(j,1,1)*cphi(j,1,1)+         !For the lower wall
     &              +Vpre(j,1)*((cphi(j,1,1))**2-(sphi(j,1,1))**2)
      Hgw(j,2)=Hpre(j,NN)                                       !For the Upper wall
      Cgw(j,2)=Cpre(j,NN)                                       !For the Upper wall
      Ugw(j,2)=Upre(j,NN)*((sphi(j,NN,3))**2-(cphi(j,NN,3))**2)-  !For the Upper wall
     &          2.0*Vpre(j,NN)*sphi(j,NN,3)*cphi(j,NN,3)
      Vgw(j,2)=-2.0*Upre(j,NN)*sphi(j,NN,3)*cphi(j,NN,3)+         !For the Upper wall
     &            +Vpre(j,NN)*((cphi(j,NN,3))**2-(sphi(j,NN,3))**2)
 104  continue 
      endif
*=============================================================================      
*     Data Projection "using the predictor values to extrapolate the cell interfaces"         
*     calculation of the interfaces values for the four faces
*     THE RIGHT FACE "we will calculate the values for the permitive values on the left and right of the RIGHT FACE"
*     the left and right of the face is depending on the direction of the calculation
*     Faces 1 and 3 are the lower and upper faces
*     Faces 2 and 4 are the Right and Left faces
*     If the walls could be perpendicular to the zai dirction the the following part will need some modification for 
*     the gradiant across the the cell face at the wall
1509  if (output_option.EQ.2) then
      Print *,'Interface values for the Lower FACE'
      endif
      do 55 k=1,NN
      do 56 j=1,MM
      if (k.EQ.1) then                   
      hl_eta(j,k,1)=Hgw(j,1)+0.5*(-Lim_H(j,k,2))        !the gradients in eta dirctions are enforced to be zero
      ul_eta(j,k,1)=Ugw(j,1)+0.5*(Lim_U(j,k,2)*(cphi(j,k,1)*
     &   cphi(j,k,1)-sphi(j,k,1)*sphi(j,k,1))+2*Lim_V(j,k,2)*sphi(j,k,1)
     &              *cphi(j,k,1))
      vl_eta(j,k,1)=Vgw(j,1)+0.5*(Lim_V(j,k,2)*(sphi(j,k,1)*
     &   sphi(j,k,1)-cphi(j,k,1)*cphi(j,k,1))+2*Lim_U(j,k,2)*sphi(j,k,1)
     &              *cphi(j,k,1))
      cl_eta(j,k,1)=Cgw(j,1)+0.5*(-Lim_C(j,k,2))
      hr_eta(j,k,1)=Hpre(j,k)-0.5*Lim_H(j,k,2)            
      ur_eta(j,k,1)=Upre(j,k)-0.5*Lim_U(j,k,2)
      vr_eta(j,k,1)=Vpre(j,k)-0.5*Lim_V(j,k,2)
      cr_eta(j,k,1)=Cpre(j,k)-0.5*Lim_C(j,k,2)
      endif
      if (k.GT.1) then
      hl_eta(j,k,1)=Hpre(j,k-1)+0.5*Lim_H(j,k-1,2)
      ul_eta(j,k,1)=Upre(j,k-1)+0.5*Lim_U(j,k-1,2)
      vl_eta(j,k,1)=Vpre(j,k-1)+0.5*Lim_V(j,k-1,2)
      cl_eta(j,k,1)=Cpre(j,k-1)+0.5*Lim_C(j,k-1,2)
      hr_eta(j,k,1)=Hpre(j,k)-0.5*Lim_H(j,k,2)            
      ur_eta(j,k,1)=Upre(j,k)-0.5*Lim_U(j,k,2)
      vr_eta(j,k,1)=Vpre(j,k)-0.5*Lim_V(j,k,2)
      cr_eta(j,k,1)=Cpre(j,k)-0.5*Lim_C(j,k,2)
      endif
56    continue
55    continue 
      if (output_option.EQ.2) then
      Print *,'Interface values for the Upper FACE'
      endif
 !     if (t.EQ.220) then
 !     pause
 !     endif
       do 57 k=1,NN
      do 58 j=1,MM
      if (k.EQ.NN) then                    
      hr_eta(j,k,3)=Hgw(j,2)-0.5*(-Lim_H(j,k,2))
      ur_eta(j,k,3)=Ugw(j,2)-0.5*(Lim_U(j,k,2)*(cphi(j,k,3)*
     &   cphi(j,k,3)-sphi(j,k,3)*sphi(j,k,3))+2*Lim_V(j,k,2)*sphi(j,k,3)
     &              *cphi(j,k,3))
      vr_eta(j,k,3)=Vgw(j,2)-0.5*(Lim_V(j,k,2)*(sphi(j,k,3)*
     &   sphi(j,k,3)-cphi(j,k,3)*cphi(j,k,3))+2*Lim_U(j,k,2)*sphi(j,k,3)
     &              *cphi(j,k,3))
      cr_eta(j,k,3)=cgw(j,2)-0.5*(-Lim_C(j,k,2))
      hl_eta(j,k,3)=Hpre(j,k)+0.5*Lim_H(j,k,2)
      ul_eta(j,k,3)=Upre(j,k)+0.5*Lim_U(j,k,2)
      vl_eta(j,k,3)=Vpre(j,k)+0.5*Lim_V(j,k,2)
      cl_eta(j,k,3)=Cpre(j,k)+0.5*Lim_C(j,k,2)
      else
      hl_eta(j,k,3)=Hpre(j,k)+0.5*Lim_H(j,k,2)
      ul_eta(j,k,3)=Upre(j,k)+0.5*Lim_U(j,k,2)
      vl_eta(j,k,3)=Vpre(j,k)+0.5*Lim_V(j,k,2)
      cl_eta(j,k,3)=Cpre(j,k)+0.5*Lim_C(j,k,2)
      hr_eta(j,k,3)=Hpre(j,k+1)-0.5*Lim_H(j,k+1,2)            
      ur_eta(j,k,3)=Upre(j,k+1)-0.5*Lim_U(j,k+1,2)
      vr_eta(j,k,3)=Vpre(j,k+1)-0.5*Lim_V(j,k+1,2)
      cr_eta(j,k,3)=Cpre(j,k+1)-0.5*Lim_C(j,k+1,2)
      endif
58    continue
57    continue           
      if (output_option.EQ.2) then
      Print *,'Interface values for the Right FACE'
      endif
      do 59 j=1,MM
      do 60 k=1,NN
      if (j.EQ.MM.AND.BC2.EQ.4) then                         ! For the right state the gradiant is zero see Bradford paper BCs Section
      hl_zai(j,k,2)=Hpre(j,k)+0.5*Lim_H(j,k,1)
      ul_zai(j,k,2)=Upre(j,k)+0.5*Lim_U(j,k,1)
      vl_zai(j,k,2)=Vpre(j,k)+0.5*Lim_V(j,k,1)
      cl_zai(j,k,2)=Cpre(j,k)+0.5*Lim_C(j,k,1)
      hr_zai(j,k,2)=Hgo(1,k)+0.0
      if (Hgo(1,k).GT.h_cri) then
      ur_zai(j,k,2)=ul_zai(j,k,2)+(9.81**0.5)*((hl_zai(j,k,2))**0.5-
     &                                      (hr_zai(j,k,2))**0.5)                       ! this should be hr_zai
      else
      ur_zai(j,k,2)=ul_zai(j,k,2)
      endif
      vr_zai(j,k,2)=vl_zai(j,k,2)
      cr_zai(j,k,2)=cl_zai(j,k,2)
      elseif (j.EQ.MM.AND.BC2.NE.4) then
      hl_zai(j,k,2)=Hpre(j,k)+0.5*Lim_H(j,k,1)
      ul_zai(j,k,2)=Upre(j,k)+0.5*Lim_U(j,k,1)
      vl_zai(j,k,2)=Vpre(j,k)+0.5*Lim_V(j,k,1)
      cl_zai(j,k,2)=Cpre(j,k)+0.5*Lim_C(j,k,1)
      hr_zai(j,k,2)=Hgo(1,k)
      ur_zai(j,k,2)=Ugo(1,k)
      vr_zai(j,k,2)=Vgo(1,k)
      cr_zai(j,k,2)=Cgo(1,k)
      else
      hl_zai(j,k,2)=Hpre(j,k)+0.5*Lim_H(j,k,1)
      ul_zai(j,k,2)=Upre(j,k)+0.5*Lim_U(j,k,1)
      vl_zai(j,k,2)=Vpre(j,k)+0.5*Lim_V(j,k,1)
      cl_zai(j,k,2)=Cpre(j,k)+0.5*Lim_C(j,k,1)
      hr_zai(j,k,2)=Hpre(j+1,k)-0.5*Lim_H(j+1,k,1)            
      ur_zai(j,k,2)=Upre(j+1,k)-0.5*Lim_U(j+1,k,1)
      vr_zai(j,k,2)=Vpre(j+1,k)-0.5*Lim_V(j+1,k,1)
      cr_zai(j,k,2)=Cpre(j+1,k)-0.5*Lim_C(j+1,k,1)
      endif
60    continue
59    continue
      if (output_option.EQ.2) then
      Print *,'Interface values for the Left FACE'
      endif
      do 61 j=1,MM
      do 62 k=1,NN
      if (j.EQ.1) then
      hl_zai(j,k,4)=Hgi(1,k)
      ul_zai(j,k,4)=Ugi(1,k)
      vl_zai(j,k,4)=Vgi(1,k)
      cl_zai(j,k,4)=Cgi(1,k)
      hr_zai(j,k,4)=Hpre(j,k)-0.5*Lim_H(j,k,1)            
      ur_zai(j,k,4)=Upre(j,k)-0.5*Lim_U(j,k,1)
      vr_zai(j,k,4)=Vpre(j,k)-0.5*Lim_V(j,k,1)
      cr_zai(j,k,4)=Cpre(j,k)-0.5*Lim_C(j,k,1)
      else
      hl_zai(j,k,4)=Hpre(j-1,k)+0.5*Lim_H(j-1,k,1)
      ul_zai(j,k,4)=Upre(j-1,k)+0.5*Lim_U(j-1,k,1)
      vl_zai(j,k,4)=Vpre(j-1,k)+0.5*Lim_V(j-1,k,1)
      cl_zai(j,k,4)=Cpre(j-1,k)+0.5*Lim_C(j-1,k,1)
      hr_zai(j,k,4)=Hpre(j,k)-0.5*Lim_H(j,k,1)            
      ur_zai(j,k,4)=Upre(j,k)-0.5*Lim_U(j,k,1)
      vr_zai(j,k,4)=Vpre(j,k)-0.5*Lim_V(j,k,1)
      cr_zai(j,k,4)=Cpre(j,k)-0.5*Lim_C(j,k,1)
      endif
62    continue
61    continue
*******************************************************************************************************************************
*     ROE SOLVER    "Data Evolving" "Approximate Riemann problem solution"  -Roe-Pike method 1984-
*     Computing the Roe average values 
      do 66 j=1,MM
      do 67 k=1,NN
      ! Lower Face
      H_hat(j,k,1)=(hl_eta(j,k,1)*hr_eta(j,k,1))**(0.5)
      U_hat(j,k,1)=(((hl_eta(j,k,1))**(0.5))*ul_eta(j,k,1)+
     &              ((hr_eta(j,k,1))**(0.5))*ur_eta(j,k,1))/
     &              ((hl_eta(j,k,1))**(0.5)+(hr_eta(j,k,1))**(0.5))
      V_hat(j,k,1)=(((hl_eta(j,k,1))**(0.5))*vl_eta(j,k,1)+
     &              ((hr_eta(j,k,1))**(0.5))*vr_eta(j,k,1))/
     &              ((hl_eta(j,k,1))**(0.5)+(hr_eta(j,k,1))**(0.5))
      C_hat(j,k,1)=(((hl_eta(j,k,1))**(0.5))*cl_eta(j,k,1)+
     &              ((hr_eta(j,k,1))**(0.5))*cr_eta(j,k,1))/
     &              ((hl_eta(j,k,1))**(0.5)+(hr_eta(j,k,1))**(0.5))
      a_hat(j,k,1)=(0.5*9.81*(hl_eta(j,k,1)+hr_eta(j,k,1)))**0.5
      ! Right Face
      H_hat(j,k,2)=(hl_zai(j,k,2)*hr_zai(j,k,2))**(0.5)
      U_hat(j,k,2)=(((hl_zai(j,k,2))**(0.5))*ul_zai(j,k,2)+
     &              ((hr_zai(j,k,2))**(0.5))*ur_zai(j,k,2))/
     &              ((hl_zai(j,k,2))**(0.5)+(hr_zai(j,k,2))**(0.5))
      V_hat(j,k,2)=(((hl_zai(j,k,2))**(0.5))*vl_zai(j,k,2)+
     &              ((hr_zai(j,k,2))**(0.5))*vr_zai(j,k,2))/
     &              ((hl_zai(j,k,2))**(0.5)+(hr_zai(j,k,2))**(0.5))
      C_hat(j,k,2)=(((hl_zai(j,k,2))**(0.5))*cl_zai(j,k,2)+
     &              ((hr_zai(j,k,2))**(0.5))*cr_zai(j,k,2))/
     &              ((hl_zai(j,k,2))**(0.5)+(hr_zai(j,k,2))**(0.5))
      a_hat(j,k,2)=(0.5*9.81*(hl_zai(j,k,2)+hr_zai(j,k,2)))**0.5
      ! Top Face
      H_hat(j,k,3)=(hl_eta(j,k,3)*hr_eta(j,k,3))**(0.5)
      U_hat(j,k,3)=(((hl_eta(j,k,3))**(0.5))*ul_eta(j,k,3)+
     &              ((hr_eta(j,k,3))**(0.5))*ur_eta(j,k,3))/
     &              ((hl_eta(j,k,3))**(0.5)+(hr_eta(j,k,3))**(0.5))
      V_hat(j,k,3)=(((hl_eta(j,k,3))**(0.5))*vl_eta(j,k,3)+
     &              ((hr_eta(j,k,3))**(0.5))*vr_eta(j,k,3))/
     &              ((hl_eta(j,k,3))**(0.5)+(hr_eta(j,k,3))**(0.5))
      C_hat(j,k,3)=(((hl_eta(j,k,3))**(0.5))*cl_eta(j,k,3)+
     &              ((hr_eta(j,k,3))**(0.5))*cr_eta(j,k,3))/
     &              ((hl_eta(j,k,3))**(0.5)+(hr_eta(j,k,3))**(0.5))
      a_hat(j,k,3)=(0.5*9.81*(hl_eta(j,k,3)+hr_eta(j,k,3)))**0.5
      ! Left Face
      H_hat(j,k,4)=(hl_zai(j,k,4)*hr_zai(j,k,4))**(0.5)
      U_hat(j,k,4)=(((hl_zai(j,k,4))**(0.5))*ul_zai(j,k,4)+
     &              ((hr_zai(j,k,4))**(0.5))*ur_zai(j,k,4))/
     &              ((hl_zai(j,k,4))**(0.5)+(hr_zai(j,k,4))**(0.5))
      V_hat(j,k,4)=(((hl_zai(j,k,4))**(0.5))*vl_zai(j,k,4)+
     &              ((hr_zai(j,k,4))**(0.5))*vr_zai(j,k,4))/
     &              ((hl_zai(j,k,4))**(0.5)+(hr_zai(j,k,4))**(0.5))
      C_hat(j,k,4)=(((hl_zai(j,k,4))**(0.5))*cl_zai(j,k,4)+
     &              ((hr_zai(j,k,4))**(0.5))*cr_zai(j,k,4))/
     &              ((hl_zai(j,k,4))**(0.5)+(hr_zai(j,k,4))**(0.5))
      a_hat(j,k,4)=(0.5*9.81*(hl_zai(j,k,4)+hr_zai(j,k,4)))**0.5
  67  continue
  66  continue
*     Computing the Eign values      
      do 68 j=1,MM
      do 69 k=1,NN
      U_hat_perp(j,k,1)=U_hat(j,k,1)*cphi(j,k,1)+
     &                  V_hat(j,k,1)*sphi(j,k,1)
      U_hat_perp(j,k,2)=U_hat(j,k,2)*cphi(j,k,2)+
     &                  V_hat(j,k,2)*sphi(j,k,2)
      U_hat_perp(j,k,3)=U_hat(j,k,3)*cphi(j,k,3)+
     &                  V_hat(j,k,3)*sphi(j,k,3)
      U_hat_perp(j,k,4)=U_hat(j,k,4)*cphi(j,k,4)+
     &                  V_hat(j,k,4)*sphi(j,k,4)
      do 120 kk=1,4
      Eign(j,k,kk,1)=abs(U_hat_perp(j,k,kk)-a_hat(j,k,kk))
      Eign(j,k,kk,3)=abs(U_hat_perp(j,k,kk)+a_hat(j,k,kk))
      Eign(j,k,kk,2)=abs(U_hat_perp(j,k,kk))
      Eign(j,k,kk,4)=abs(U_hat_perp(j,k,kk))
 120  continue     
 69   continue
 68   continue
*     Entropy Fix for the Eign values
      do 70 j=1,MM
      do 71 k=1,NN
      U_hat_perp_L(j,k,1)=ul_eta(j,k,1)*cphi(j,k,1)+
     &                    vl_eta(j,k,1)*sphi(j,k,1)
      a_hat_L(j,k,1)=(hl_eta(j,k,1)*9.81)**0.5
      U_hat_perp_L(j,k,2)=ul_zai(j,k,2)*cphi(j,k,2)+
     &                    vl_zai(j,k,2)*sphi(j,k,2)
      a_hat_L(j,k,2)=(hl_zai(j,k,2)*9.81)**0.5
      U_hat_perp_L(j,k,3)=ul_eta(j,k,3)*cphi(j,k,3)+
     &                    vl_eta(j,k,3)*sphi(j,k,3)
      a_hat_L(j,k,3)=(hl_eta(j,k,3)*9.81)**0.5
      U_hat_perp_L(j,k,4)=ul_zai(j,k,4)*cphi(j,k,4)+
     &                    vl_zai(j,k,4)*sphi(j,k,4)
      a_hat_L(j,k,4)=(hl_zai(j,k,4)*9.81)**0.5
      U_hat_perp_R(j,k,1)=ur_eta(j,k,1)*cphi(j,k,1)+
     &                    vr_eta(j,k,1)*sphi(j,k,1)
      a_hat_R(j,k,1)=(hr_eta(j,k,1)*9.81)**0.5
      U_hat_perp_R(j,k,2)=ur_zai(j,k,2)*cphi(j,k,2)+
     &                    vr_zai(j,k,2)*sphi(j,k,2)
      a_hat_R(j,k,2)=(hr_zai(j,k,2)*9.81)**0.5
      U_hat_perp_R(j,k,3)=ur_eta(j,k,3)*cphi(j,k,3)+
     &                    vr_eta(j,k,3)*sphi(j,k,3)
      a_hat_R(j,k,3)=(hr_eta(j,k,3)*9.81)**0.5
      U_hat_perp_R(j,k,4)=ur_zai(j,k,4)*cphi(j,k,4)+
     &                    vr_zai(j,k,4)*sphi(j,k,4)
      a_hat_R(j,k,4)=(hr_zai(j,k,4)*9.81)**0.5
      do 72 kk=1,4
      Eign1_L(j,k,kk)=abs(U_hat_perp_L(j,k,kk)-a_hat_L(j,k,kk))
      Eign3_L(j,k,kk)=abs(U_hat_perp_L(j,k,kk)+a_hat_L(j,k,kk))
      Eign1_R(j,k,kk)=abs(U_hat_perp_R(j,k,kk)-a_hat_R(j,k,kk))
      Eign3_R(j,k,kk)=abs(U_hat_perp_R(j,k,kk)+a_hat_R(j,k,kk))
      dEign1(j,k,kk)=4.0*(Eign1_R(j,k,kk)-Eign1_L(j,k,kk))
      dEign3(j,k,kk)=4.0*(Eign3_R(j,k,kk)-Eign3_L(j,k,kk))
      if ((Eign(j,k,kk,1).GT.(-0.5*dEign1(j,k,kk))).AND.
     &   (Eign(j,k,kk,1).LT.(0.5*dEign1(j,k,kk)))) then
      Eign(j,k,kk,1)=((Eign(j,k,kk,1))**(0.5))/dEign1(j,k,kk)+
     &                dEign1(j,k,kk)*0.25
      endif
      if ((Eign(j,k,kk,3).GT.(-0.5*dEign3(j,k,kk))).AND.
     &   (Eign(j,k,kk,3).LT.(0.5*dEign3(j,k,kk)))) then
      Eign(j,k,kk,3)=((Eign(j,k,kk,3))**(0.5))/dEign3(j,k,kk)+
     &                dEign3(j,k,kk)*0.25
      endif
  72  continue  
  71  continue
  70  continue 
*     Computing the right Eign vectors
      do 73 j=1,MM
      do 74 k=1,NN
      do 75 kk=1,4
*     1st Right Eign Vector
      R_Eign1(j,k,kk,1)=1
      R_Eign1(j,k,kk,2)=U_hat(j,k,kk)-a_hat(j,k,kk)*cphi(j,k,kk)
      R_Eign1(j,k,kk,3)=V_hat(j,k,kk)-a_hat(j,k,kk)*sphi(j,k,kk)
      R_Eign1(j,k,kk,4)=C_hat(j,k,kk)
*     2nd Right Eign Vector
      R_Eign2(j,k,kk,1)=0
      R_Eign2(j,k,kk,2)=-sphi(j,k,kk)
      R_Eign2(j,k,kk,3)=cphi(j,k,kk)
      R_Eign2(j,k,kk,4)=0
*     3ed Right Eign Vector
      R_Eign3(j,k,kk,1)=1
      R_Eign3(j,k,kk,2)=U_hat(j,k,kk)+a_hat(j,k,kk)*cphi(j,k,kk)
      R_Eign3(j,k,kk,3)=V_hat(j,k,kk)+a_hat(j,k,kk)*sphi(j,k,kk)
      R_Eign3(j,k,kk,4)=C_hat(j,k,kk)
*     4th Right Eign Vector
      R_Eign4(j,k,kk,1)=1
      R_Eign4(j,k,kk,2)=U_hat(j,k,kk)
      R_Eign4(j,k,kk,3)=V_hat(j,k,kk)
      R_Eign4(j,k,kk,4)=-C_hat(j,k,kk)
  75  continue
  74  continue
  73  continue
*     Computing the wave strength
      do 76 j=1,MM
      do 77 k=1,NN
*     Lower Face
      dh(j,k,1)=hr_eta(j,k,1)-hl_eta(j,k,1)
      du(j,k,1)=ur_eta(j,k,1)-ul_eta(j,k,1)
      dv(j,k,1)=vr_eta(j,k,1)-vl_eta(j,k,1)
      du_perp(j,k,1)=du(j,k,1)*cphi(j,k,1)+dv(j,k,1)*sphi(j,k,1)
      du_par(j,k,1)=-du(j,k,1)*sphi(j,k,1)+dv(j,k,1)*cphi(j,k,1)
      dc(j,k,1)=cr_eta(j,k,1)-cl_eta(j,k,1)
*     Right Face
      dh(j,k,2)=hr_zai(j,k,2)-hl_zai(j,k,2)
      du(j,k,2)=ur_zai(j,k,2)-ul_zai(j,k,2)
      dv(j,k,2)=vr_zai(j,k,2)-vl_zai(j,k,2)
      du_perp(j,k,2)=du(j,k,2)*cphi(j,k,2)+dv(j,k,2)*sphi(j,k,2)
      du_par(j,k,2)=-du(j,k,2)*sphi(j,k,2)+dv(j,k,2)*cphi(j,k,2)
      dc(j,k,2)=cr_zai(j,k,2)-cl_zai(j,k,2)
*     Top Face
      dh(j,k,3)=hr_eta(j,k,3)-hl_eta(j,k,3)
      du(j,k,3)=ur_eta(j,k,3)-ul_eta(j,k,3)
      dv(j,k,3)=vr_eta(j,k,3)-vl_eta(j,k,3)
      du_perp(j,k,3)=du(j,k,3)*cphi(j,k,3)+dv(j,k,3)*sphi(j,k,3)
      du_par(j,k,3)=-du(j,k,3)*sphi(j,k,3)+dv(j,k,3)*cphi(j,k,3)
      dc(j,k,3)=cr_eta(j,k,3)-cl_eta(j,k,3)
*     Left Face
      dh(j,k,4)=hr_zai(j,k,4)-hl_zai(j,k,4)
      du(j,k,4)=ur_zai(j,k,4)-ul_zai(j,k,4)
      dv(j,k,4)=vr_zai(j,k,4)-vl_zai(j,k,4)
      du_perp(j,k,4)=du(j,k,4)*cphi(j,k,4)+dv(j,k,4)*sphi(j,k,4)
      du_par(j,k,4)=-du(j,k,4)*sphi(j,k,4)+dv(j,k,4)*cphi(j,k,4)
      dc(j,k,4)=cr_zai(j,k,4)-cl_zai(j,k,4)
  77  continue
  76  continue
*     Wave Strengths                                        NEEDS MORE ATTENTION"Done"
      do 78 j=1,MM
      do 79 k=1,NN
      do 80 kk=1,4
      Ws(j,k,kk,1)=0.5*(dh(j,k,kk)-H_hat(j,k,kk)*du_perp(j,k,kk)/
     &             a_hat(j,k,kk)+H_hat(j,k,kk)*dc(j,k,kk)/
     &              (2*C_hat(j,k,kk)))
      Ws(j,k,kk,2)=H_hat(j,k,kk)*du_par(j,k,kk)
      Ws(j,k,kk,3)=0.5*(dh(j,k,kk)+H_hat(j,k,kk)*du_perp(j,k,kk)/
     &             a_hat(j,k,kk)+H_hat(j,k,kk)*dc(j,k,kk)/
     &             (2*C_hat(j,k,kk)))
      Ws(j,k,kk,4)=H_hat(j,k,kk)*(-dc(j,k,kk))/(2*C_hat(j,k,kk))
  80  continue
  79  continue
  78  continue
*     Computing the left and right interfacial fluxes  
      ! Compution of the flux change increment for all the faces
      XX=1
      do 81 j=1,MM
      do 82 k=1,NN
      do 83 kk=1,4
      df(j,k,kk,1)=R_Eign1(j,k,kk,1)*Eign(j,k,kk,1)*Ws(j,k,kk,1)+
     &             R_Eign2(j,k,kk,1)*Eign(j,k,kk,2)*Ws(j,k,kk,2)+
     &             R_Eign3(j,k,kk,1)*Eign(j,k,kk,3)*Ws(j,k,kk,3)+
     &             R_Eign4(j,k,kk,1)*Eign(j,k,kk,4)*Ws(j,k,kk,4)
      df(j,k,kk,2)=R_Eign1(j,k,kk,2)*Eign(j,k,kk,1)*Ws(j,k,kk,1)+
     &             R_Eign2(j,k,kk,2)*Eign(j,k,kk,2)*Ws(j,k,kk,2)+
     &             R_Eign3(j,k,kk,2)*Eign(j,k,kk,3)*Ws(j,k,kk,3)+
     &             R_Eign4(j,k,kk,2)*Eign(j,k,kk,4)*Ws(j,k,kk,4)
      df(j,k,kk,3)=R_Eign1(j,k,kk,3)*Eign(j,k,kk,1)*Ws(j,k,kk,1)+
     &             R_Eign2(j,k,kk,3)*Eign(j,k,kk,2)*Ws(j,k,kk,2)+
     &             R_Eign3(j,k,kk,3)*Eign(j,k,kk,3)*Ws(j,k,kk,3)+
     &             R_Eign4(j,k,kk,3)*Eign(j,k,kk,4)*Ws(j,k,kk,4)
      df(j,k,kk,4)=R_Eign1(j,k,kk,4)*Eign(j,k,kk,1)*Ws(j,k,kk,1)+
     &             R_Eign2(j,k,kk,4)*Eign(j,k,kk,2)*Ws(j,k,kk,2)+
     &             R_Eign3(j,k,kk,4)*Eign(j,k,kk,3)*Ws(j,k,kk,3)+
     &             R_Eign4(j,k,kk,4)*Eign(j,k,kk,4)*Ws(j,k,kk,4)
   83 continue
   82 continue
   81 continue
*============================The dry wet tracking needs to be added======================================
* the fluxes on left and right of each face will be calculated then the new flux at time step n+1 will be determined
      do 85 j=1,MM
      do 86 k=1,NN
      ul_perp(j,k,1)=ul_eta(j,k,1)*cphi(j,k,1)+vl_eta(j,k,1)*sphi(j,k,1)
      ul_perp(j,k,2)=ul_zai(j,k,2)*cphi(j,k,2)+vl_zai(j,k,2)*sphi(j,k,2)
      ul_perp(j,k,3)=ul_eta(j,k,3)*cphi(j,k,3)+vl_eta(j,k,3)*sphi(j,k,3)
      ul_perp(j,k,4)=ul_zai(j,k,4)*cphi(j,k,4)+vl_zai(j,k,4)*sphi(j,k,4)
      ur_perp(j,k,1)=ur_eta(j,k,1)*cphi(j,k,1)+vr_eta(j,k,1)*sphi(j,k,1)
      ur_perp(j,k,2)=ur_zai(j,k,2)*cphi(j,k,2)+vr_zai(j,k,2)*sphi(j,k,2)
      ur_perp(j,k,3)=ur_eta(j,k,3)*cphi(j,k,3)+vr_eta(j,k,3)*sphi(j,k,3)
      ur_perp(j,k,4)=ur_zai(j,k,4)*cphi(j,k,4)+vr_zai(j,k,4)*sphi(j,k,4)
*     The lower face
      F1_L(j,k,1)=f1(hl_eta(j,k,1),ul_perp(j,k,1))
      F1_R(j,k,1)=f1(hr_eta(j,k,1),ur_perp(j,k,1))
      F2_L(j,k,1)=f2(hl_eta(j,k,1),ul_eta(j,k,1),ul_perp(j,k,1),
     &              cphi(j,k,1))
      F2_R(j,k,1)=f2(hr_eta(j,k,1),ur_eta(j,k,1),ur_perp(j,k,1),
     &              cphi(j,k,1))
      F3_L(j,k,1)=f3(hl_eta(j,k,1),vl_eta(j,k,1),ul_perp(j,k,1),
     &              sphi(j,k,1))
      F3_R(j,k,1)=f3(hr_eta(j,k,1),vr_eta(j,k,1),ur_perp(j,k,1),
     &              sphi(j,k,1))
      F4_L(j,k,1)=f4(hl_eta(j,k,1),cl_eta(j,k,1),ul_perp(j,k,1))
      F4_R(j,k,1)=f4(hr_eta(j,k,1),cr_eta(j,k,1),ur_perp(j,k,1))
*     The Right face
      F1_L(j,k,2)=f1(hl_zai(j,k,2),ul_perp(j,k,2))
      F1_R(j,k,2)=f1(hr_zai(j,k,2),ur_perp(j,k,2))
      F2_L(j,k,2)=f2(hl_zai(j,k,2),ul_zai(j,k,2),ul_perp(j,k,2),
     &              cphi(j,k,2))
      F2_R(j,k,2)=f2(hr_zai(j,k,2),ur_zai(j,k,2),ur_perp(j,k,2),
     &              cphi(j,k,2))
      F3_L(j,k,2)=f3(hl_zai(j,k,2),vl_zai(j,k,2),ul_perp(j,k,2),
     &              sphi(j,k,2))
      F3_R(j,k,2)=f3(hr_zai(j,k,2),vr_zai(j,k,2),ur_perp(j,k,2),
     &              sphi(j,k,2))
      F4_L(j,k,2)=f4(hl_zai(j,k,2),cl_zai(j,k,2),ul_perp(j,k,2))
      F4_R(j,k,2)=f4(hr_zai(j,k,2),cr_zai(j,k,2),ur_perp(j,k,2))
*     The Top face
      F1_L(j,k,3)=f1(hl_eta(j,k,3),ul_perp(j,k,3))
      F1_R(j,k,3)=f1(hr_eta(j,k,3),ur_perp(j,k,3))
      F2_L(j,k,3)=f2(hl_eta(j,k,3),ul_eta(j,k,3),ul_perp(j,k,3),
     &              cphi(j,k,3))
      F2_R(j,k,3)=f2(hr_eta(j,k,3),ur_eta(j,k,3),ur_perp(j,k,3),
     &              cphi(j,k,3))
      F3_L(j,k,3)=f3(hl_eta(j,k,3),vl_eta(j,k,3),ul_perp(j,k,3),
     &              sphi(j,k,3))
      F3_R(j,k,3)=f3(hr_eta(j,k,3),vr_eta(j,k,3),ur_perp(j,k,3),
     &              sphi(j,k,3))
      F4_L(j,k,3)=f4(hl_eta(j,k,3),cl_eta(j,k,3),ul_perp(j,k,3))
      F4_R(j,k,3)=f4(hr_eta(j,k,3),cr_eta(j,k,3),ur_perp(j,k,3))
*     The Left face      
      F1_L(j,k,4)=f1(hl_zai(j,k,4),ul_perp(j,k,4))
      F1_R(j,k,4)=f1(hr_zai(j,k,4),ur_perp(j,k,4))
      F2_L(j,k,4)=f2(hl_zai(j,k,4),ul_zai(j,k,4),ul_perp(j,k,4),
     &              cphi(j,k,4))
      F2_R(j,k,4)=f2(hr_zai(j,k,4),ur_zai(j,k,4),ur_perp(j,k,4),
     &              cphi(j,k,4))
      F3_L(j,k,4)=f3(hl_zai(j,k,4),vl_zai(j,k,4),ul_perp(j,k,4),
     &              sphi(j,k,4))
      F3_R(j,k,4)=f3(hr_zai(j,k,4),vr_zai(j,k,4),ur_perp(j,k,4),
     &              sphi(j,k,4))
      F4_L(j,k,4)=f4(hl_zai(j,k,4),cl_zai(j,k,4),ul_perp(j,k,4))
      F4_R(j,k,4)=f4(hr_zai(j,k,4),cr_zai(j,k,4),ur_perp(j,k,4))
*=================================================================================      
*     Pressure term correction " For more detailes see Sanders Paper 2002" THIS PART NEEDS TO BE CORRECTED TO WORK WITH TC
      dh_corec2(j,k,1)=0.5*9.81*cphi(j,k,1)*
     &                  ((Hc(j+1,k)-Hc(j,k))**2)/12
      dh_corec3(j,k,1)=0.5*9.81*sphi(j,k,1)*
     &                  ((Hc(j+1,k)-Hc(j,k))**2)/12
      dh_corec2(j,k,2)=0.5*9.81*cphi(j,k,2)*
     &                  ((Hc(j+1,k+1)-Hc(j+1,k))**2)/12
      dh_corec3(j,k,2)=0.5*9.81*sphi(j,k,2)*
     &                  ((Hc(j+1,k+1)-Hc(j+1,k))**2)/12
      dh_corec2(j,k,3)=0.5*9.81*cphi(j,k,3)*
     &                  ((Hc(j+1,k+1)-Hc(j,k+1))**2)/12
      dh_corec3(j,k,3)=0.5*9.81*sphi(j,k,3)*
     &                  ((Hc(j+1,k+1)-Hc(j,k+1))**2)/12
      dh_corec2(j,k,4)=0.5*9.81*cphi(j,k,4)*
     &                  ((Hc(j,k+1)-Hc(j,k))**2)/12
      dh_corec3(j,k,4)=0.5*9.81*sphi(j,k,4)*
     &                  ((Hc(j,k+1)-Hc(j,k))**2)/12
*===================================================================================    
*     Computing the Fluxes at cell faces
      do 87 kk=1,4
      Flux1(j,k,kk)=0.5*(F1_R(j,k,kk)+F1_L(j,k,kk)-df(j,k,kk,1))
      Flux2(j,k,kk)=0.5*(F2_R(j,k,kk)+F2_L(j,k,kk)-df(j,k,kk,2)
     &                      +dh_corec2(j,k,kk)*H_correction)
      Flux3(j,k,kk)=0.5*(F3_R(j,k,kk)+F3_L(j,k,kk)-df(j,k,kk,3)
     &                      +dh_corec2(j,k,kk)*H_correction)
      Flux4(j,k,kk)=0.5*(F4_R(j,k,kk)+F4_L(j,k,kk)-df(j,k,kk,4))
          ! Dry Bed Compution
      
      if (hl_eta(j,k,1).LE.epso .AND.hr_eta(j,k,1).LE.epso)then
      Flux1(j,k,1)=0.0
      Flux2(j,k,1)=0.0
      Flux3(j,k,1)=0.0
      Flux4(j,k,1)=0.0
      endif
      if (hl_eta(j,k,3).LE.epso .AND.hr_eta(j,k,3).LE.epso)then
      Flux1(j,k,3)=0.0
      Flux2(j,k,3)=0.0
      Flux3(j,k,3)=0.0
      Flux4(j,k,3)=0.0
      endif
      if (hl_zai(j,k,2).LE.epso .AND.hr_zai(j,k,2).LE.epso)then
      Flux1(j,k,2)=0.0
      Flux2(j,k,2)=0.0
      Flux3(j,k,2)=0.0
      Flux4(j,k,2)=0.0
      endif
      if (hl_zai(j,k,4).LE.epso .AND.hr_zai(j,k,4).LE.epso)then
      Flux1(j,k,4)=0.0
      Flux2(j,k,4)=0.0
      Flux3(j,k,4)=0.0
      Flux4(j,k,4)=0.0
      endif
  87  continue    
  86  continue
  85  continue   
      
*=============================================================================================================================
*     Source term Calculations "we want to make external functions for all the source term parameters"
      ! I am not sure if we need to calculate the source term based on the previouse time step"t=n" or based on 
      ! the predicted solution "t=n+1/2"
*     Fall velocity "Dietrich 1982"
      vfall=v_fall(R,D50,neu)
      Rep=sqrt(R*9.81*D50)*D50/neu
      do 88 j=1,MM
      do 89 k=1,NN
      if (source_option.EQ.1) then                          !source term is calculated based on values at t=n
*     Near bed concentration        "Van Rijn 1984" &"Parker"
      if (cb_method.EQ.1.0) then
      cb(j,k)=c_bed1(C_old(j,k))
      else
      cb(j,k)=c_bed2(D50,R,neu,Vmag_old(j,k),D50)
      endif
*     Drag coeeficient                  "Haaland 1983"
      cd(j,k)=drag(Vmag_old(j,k),H_old(j,k),neu,ks)
*     Shear velocity
      Ustar(j,k)=cd(j,k)*Vmag_old(j,k)*U_old(j,k)
      Vstar(j,k)=cd(j,k)*Vmag_old(j,k)*V_old(j,k)
*     Bed slop term discretization
      dz_x_old(j,k)=((zo(j+1,k)-zo(j,k+1))*(yo(j+1,k+1)-yo(j,k))-
     &              ((zo(j+1,k+1)-zo(j,k))*(yo(j+1,k)-yo(j,k+1))))/
     &              (2*CellV(j,k))
      dz_y_old(j,k)=((zo(j+1,k+1)-zo(j,k))*(xo(j+1,k)-xo(j,k+1))-
     &              ((zo(j+1,k)-zo(j,k+1))*(xo(j+1,k+1)-xo(j,k))))/
     &              (2*CellV(j,k))
*     Dimensionless Entrainment rate
      if (E_method.EQ.1) then
      s=SQRT(dz_x_old(j,k)**2+dz_y_old(j,k)**2)
      E(j,k)=Es(Ustar(j,k),Vstar(j,k),vfall,s,Rep)
      else
      E(j,k)=E_van(Ustar(j,k),Vstar(j,k),D50,neu,H_old(j,k),R)
      endif
*     Source term calculation
      q1(j,k)=0.0             ! mas conservation equation
      q2(j,k)=-9.81*H_old(j,k)*dz_x_old(j,k)*CellV(j,k)-Ustar(j,k)
      q3(j,k)=-9.81*H_old(j,k)*dz_y_old(j,k)*CellV(j,k)-Vstar(j,k)
      q4(j,k)=vfall*(E(j,k)-cb(j,k))
      else
*===============================================================================
*===============================================================================
      !source term is calculated based on values at t=n+1/2
*      Near bed concentration        "Van Rijn 1984" &"Parker"
      if (cb_method.EQ.1.0) then
      cb_new(j,k)=c_bed1(Cpre(j,k))
      else
      cb_new(j,k)=c_bed2(D50,R,neu,Vmag_pre(j,k),D50)
      endif
*     Drag coeeficient                  "Haaland 1983"
      cd_new(j,k)=drag(Vmag_pre(j,k),Hpre(j,k),neu,ks)
*     Shear velocity
      Ustar_new(j,k)=cd_new(j,k)*Vmag_pre(j,k)*Upre(j,k)
      Vstar_new(j,k)=cd_new(j,k)*Vmag_pre(j,k)*Vpre(j,k)
*     Bed slop term discretization                  "if the bed is allowed to erode be sure that the zo is updated"
      dz_x(j,k)=((zo(j+1,k)-zo(j,k+1))*(yo(j+1,k+1)-yo(j,k))-
     &              ((zo(j+1,k+1)-zo(j,k))*(yo(j+1,k)-yo(j,k+1))))/
     &              (2*CellV(j,k))
      dz_y(j,k)=((zo(j+1,k+1)-zo(j,k))*(xo(j+1,k)-xo(j,k+1))-
     &              ((zo(j+1,k)-zo(j,k+1))*(xo(j+1,k+1)-xo(j,k))))/
     &              (2*CellV(j,k))
*     Dimensionless Entrainment rate
      if (E_method.EQ.1) then
      s=SQRT(dz_x_old(j,k)**2+dz_y_old(j,k)**2)
      E_new(j,k)=Es(Ustar_new(j,k),Vstar_new(j,k),vfall,s,Rep)
      else
      E_new(j,k)=E_van(Ustar_new(j,k),Vstar_new(j,k),D50,neu,Hpre(j,k)
     &          ,R)
      endif
*     Source term calculation
      q1_new(j,k)=0.0                                ! mas conservation equation
      q2_new(j,k)=-9.81*Hpre(j,k)*dz_x(j,k)*CellV(j,k)-Ustar_new(j,k)
      q3_new(j,k)=-9.81*Hpre(j,k)*dz_y(j,k)*CellV(j,k)-Vstar_new(j,k)
      q4_new(j,k)=vfall*(E_new(j,k)-cb_new(j,k))
      endif
      if (output_option.EQ.2) then
      print *,Flux1(1,1,1),Flux1(1,1,2),Flux1(1,1,3),
     &          Flux1(1,1,4)
      endif
      if (Hpre(j,k).LE.espo) then
      q1(j,k)=0.0                                                           ! mas conservation equation
      q2(j,k)=0.0
      q3(j,k)=0.0
      q4(j,k)=0.0
      q1_new(j,k)=0.0                                                           ! mas conservation equation
      q2_new(j,k)=0.0
      q3_new(j,k)=0.0
      q4_new(j,k)=0.0
      endif
*=============================================================================================================================
*     Computing the new Values for the consrvative parameters [H HUU HUV HUC]     
      W1_old(j,k)=W1(j,k)
      W2_old(j,k)=W2(j,k)
      W3_old(j,k)=W3(j,k)
      W4_old(j,k)=W4(j,k)
      if (source_option.EQ.1) then
      W1(j,k)=(dtmax/CellV(j,k))*(Flux1(j,k,1)*ds(j,k,1)-
     &                          Flux1(j,k,2)*ds(j,k,2)-
     &                          Flux1(j,k,3)*ds(j,k,3)+
     &                          Flux1(j,k,4)*ds(j,k,4)+
     &                          q1(j,k)*CellV(j,k))+W1_old(j,k)
      W2(j,k)=(dtmax/CellV(j,k))*(Flux2(j,k,1)*ds(j,k,1)-
     &                          Flux2(j,k,2)*ds(j,k,2)-
     &                          Flux2(j,k,3)*ds(j,k,3)+
     &                          Flux2(j,k,4)*ds(j,k,4)+
     &                          q2(j,k)*CellV(j,k))+W2_old(j,k)
      W3(j,k)=(dtmax/CellV(j,k))*(Flux3(j,k,1)*ds(j,k,1)-
     &                          Flux3(j,k,2)*ds(j,k,2)-
     &                          Flux3(j,k,3)*ds(j,k,3)+
     &                          Flux3(j,k,4)*ds(j,k,4)+
     &                          q3(j,k)*CellV(j,k))+W3_old(j,k)
      W4(j,k)=(dtmax/CellV(j,k))*(Flux4(j,k,1)*ds(j,k,1)-
     &                          Flux4(j,k,2)*ds(j,k,2)-
     &                          Flux4(j,k,3)*ds(j,k,3)+
     &                          Flux4(j,k,4)*ds(j,k,4)+
     &                          q4(j,k)*CellV(j,k))+W4_old(j,k)
      else
      W1(j,k)=(dtmax/CellV(j,k))*(Flux1(j,k,1)*ds(j,k,1)-
     &                          Flux1(j,k,2)*ds(j,k,2)-
     &                          Flux1(j,k,3)*ds(j,k,3)+
     &                          Flux1(j,k,4)*ds(j,k,4)+
     &                          q1_new(j,k)*CellV(j,k))+W1_old(j,k)
      W2(j,k)=(dtmax/CellV(j,k))*(Flux2(j,k,1)*ds(j,k,1)-
     &                          Flux2(j,k,2)*ds(j,k,2)-
     &                          Flux2(j,k,3)*ds(j,k,3)+
     &                          Flux2(j,k,4)*ds(j,k,4)+
     &                          q2_new(j,k)*CellV(j,k))+W2_old(j,k)
      W3(j,k)=(dtmax/CellV(j,k))*(Flux3(j,k,1)*ds(j,k,1)-
     &                          Flux3(j,k,2)*ds(j,k,2)-
     &                          Flux3(j,k,3)*ds(j,k,3)+
     &                          Flux3(j,k,4)*ds(j,k,4)+
     &                          q3_new(j,k)*CellV(j,k))+W3_old(j,k)
      W4(j,k)=(dtmax/CellV(j,k))*(Flux4(j,k,1)*ds(j,k,1)-
     &                          Flux4(j,k,2)*ds(j,k,2)-
     &                          Flux4(j,k,3)*ds(j,k,3)+
     &                          Flux4(j,k,4)*ds(j,k,4)+
     &                          q4_new(j,k)*CellV(j,k))+W4_old(j,k)
      endif
*=============================================================================================================================
*     Computing the new Values for the consrvative parameters [H U V C]
      if (W1(j,k).LE.epso) then
      H(j,k)=0.0
      U(j,k)=0.0
      V(j,k)=0.0
      C(j,k)=0.0
      else
      H(j,k)=W1(j,k)
      U(j,k)=W2(j,k)/H(j,k)
      V(j,k)=W3(j,k)/H(j,k)
      C(j,k)=W4(j,k)/H(j,k)
      endif
      if (d1d2.eq.1) then
      V=0
      endif
      Vmag(j,k)=sqrt((U(j,k))**2+(V(j,k))**2)
      dxcell(j,k)=0.5*(abs(xo(j+1,k)-xo(j,k))+
     &          abs(xo(j,k+1)-xo(j+1,k+1)))
      dycell(j,k)=0.5*(abs(yo(j+1,k+1)-yo(j+1,k))+
     &              abs(yo(j,k)-yo(j,k+1)))
      c_no(j,k)=(abs(U(j,k))+sqrt(9.81*H(j,k)))/
     &          dxcell(j,k)+(abs(V(j,k))+
     &          sqrt(9.81*H(j,k)))/dycell(j,k)

  
*********************************************************************************************************************************
*     Exner's Equation solution
 !     if (j.EQ.32.AND.k.EQ.3.AND.t.EQ.78) then
 !     pause
   !   endif
      cb_n=cb_new(j,k)
      Es_n=E_new(j,k)
      ! calculating the new cb and Es
      if (cb_method.EQ.1.0) then
      cb_n_1=c_bed1(C(j,k))
      else
      cb_n_1=c_bed2(D50,R,neu,Vmag(j,k),D50)
      endif
      cd_n_1=drag(Vmag(j,k),H(j,k),neu,ks)
      Ustar_n_1=(cd_n_1)*Vmag(j,k)*U(j,k)
      Vstar_n_1=(cd_n_1)*Vmag(j,k)*V(j,k)
      Es_n_1=E_van(Ustar_n_1,Vstar_n_1,D50,neu,H(j,k),R)
*     Bed load

*     New bed elevation
      if (t.EQ.2) then
      zp(j,k)=zc(j,k)
      endif
      zp_old(j,k)=zp(j,k)
      zp(j,k)=zp_old(j,k)+vfall*dtmax/(2*(1-lamda))*                     ! New bed elevation
     &      ((cb_n-Es_n)+((cb_n_1)-(Es_n_1)))
      if (zp(j,k).LT.zc(j,k).AND.Allow_erosion.EQ.1) then
      zp(j,k)=zc(j,k)
      endif
      if (depo.EQ.2) then
      zp(j,k)=0.0
      endif
      Zeta(j,k)=H(j,k)+zp(j,k)                                    ! New water elevation
  89  continue
  88  continue
      do 95 j=1,M                 !********Updating the bed elevation*********
      do 96 k=1,N
      if (j.EQ.1.AND.k.EQ.1) then
      zp_g1=zp(j,k)+(zp(j,k)-zp(j+1,k+1))/
     &  (sqrt((xc(j,k)-xc(j+1,k+1))**2+(yc(j,k)-yc(j+1,k+1))**2))*
     &  2*(sqrt((xc(j,k)-xo(j,k))**2+(yc(j,k)-yo(j,k))**2))
      zp_g2=zp(j,k)+(zp(j,k)-zp(j,k+1))/
     &  (sqrt((xc(j,k)-xc(j,k+1))**2+(yc(j,k)-yc(j,k+1))**2))*
     &  2*(yc(j,k)-yo(j,k))
      zp_g3=zp(j,k)+(zp(j,k)-zp(j+1,k))/
     &  (sqrt((xc(j,k)-xc(j+1,k))**2+(yc(j,k)-yc(j+1,k))**2))
     &  *2*(xc(j,k)-xo(j,k))
      zo(j,k)=0.25*(zp(j,k)+zp_g1+zp_g2+zp_g3)
      elseif (j.EQ.1.AND.k.EQ.N) then
      zp_g1=zp(j,k-1)+(zp(j,k-1)-zp(j+1,k-2))/
     &  (sqrt((xc(j,k-1)-xc(j+1,k-2))**2+(yc(j,k-1)-yc(j+1,k-2))**2))*
     &  2*(sqrt((xc(j,k-1)-xo(j,k))**2+(yc(j,k-1)-yo(j,k))**2))
      zp_g2=zp(j,k-1)+(zp(j,k-1)-zp(j,k-2))/
     &  (sqrt((xc(j,k-1)-xc(j,k-2))**2+(yc(j,k-1)-yc(j,k-2))**2))*
     &  2*(yc(j,k-1)-yo(j,k))
      zp_g3=zp(j,k-1)+(zp(j,k-1)-zp(j+1,k-1))/
     &  (sqrt((xc(j,k-1)-xc(j+1,k-1))**2+(yc(j,k-1)-yc(j+1,k-1))**2))
     &  *2*(xc(j,k-1)-xo(j,k))
      zo(j,k)=0.25*(zp(j,k-1)+zp_g1+zp_g2+zp_g3)
      elseif (j.EQ.1.And.K.GT.1.AND.k.LT.N) then
      zp_g1=zp(j,k-1)+(zp(j,k-1)-zp(j+1,k-1))/
     &  (sqrt((xc(j,k-1)-xc(j+1,k-1))**2+(yc(j,k-1)-yc(j+1,k-1))**2))
     &  *2*(xc(j,k-1)-xo(j,k))
      zp_g1=zp(j,k)+(zp(j,k)-zp(j+1,k))/
     &  (sqrt((xc(j,k)-xc(j+1,k))**2+(yc(j,k)-yc(j+1,k))**2))
     &  *2*(xc(j,k)-xo(j,k))
      zo(j,k)=0.25*(zp_g1+zp_g2+zp(j,k-1)+zp(j,k))
      elseif (k.EQ.1.And.j.GT.1.AND.j.LT.M) then
      zp_g1=zp(j-1,k)+(zp(j-1,k)-zp(j-1,k+1))/
     &  (sqrt((xc(j-1,k)-xc(j-1,k+1))**2+(yc(j-1,k)-yc(j-1,k+1))**2))
     &  *2*(yc(j-1,k)-yo(j,k))
      zp_g2=zp(j,k)+(zp(j,k)-zp(j,k+1))/
     &  (sqrt((xc(j,k)-xc(j,k+1))**2+(yc(j,k)-yc(j,k+1))**2))
     &  *2*(yc(j,k)-yo(j,k))
      zo(j,k)=0.25*(zp_g1+zp_g2+zp(j-1,k)+zp(j,k))
      elseif (k.EQ.N.And.j.GT.1.AND.j.LT.M) then
      zp_g1=zp(j-1,k-1)+(zp(j-1,k-1)-zp(j-1,k-2))/
     & (sqrt((xc(j-1,k-1)-xc(j-1,k-2))**2+(yc(j-1,k-1)-yc(j-1,k-2))**2))
     &  *2*(yc(j-1,k-1)-yo(j,k))
      zp_g2=zp(j,k-1)+(zp(j,k-1)-zp(j,k-2))/
     &  (sqrt((xc(j,k-1)-xc(j,k-2))**2+(yc(j,k-1)-yc(j,k-2))**2))
     &  *2*(yc(j,k-1)-yo(j,k))
      zo(j,k)=0.25*(zp_g1+zp_g2+zp(j-1,k-1)+zp(j,k-1))
      elseif (j.EQ.M.AND.k.EQ.1) then
      zp_g1=zp(j-1,k)+(zp(j-1,k)-zp(j-2,k+1))/
     &  (sqrt((xc(j-1,k)-xc(j-2,k+1))**2+(yc(j-1,k)-yc(j-2,k+1))**2))*
     &  2*(sqrt((xc(j-1,k)-xo(j,k))**2+(yc(j-1,k)-yo(j,k))**2))
      zp_g2=zp(j-1,k)+(zp(j-1,k)-zp(j-1,k+1))/
     &  (sqrt((xc(j-1,k)-xc(j-1,k+1))**2+(yc(j-1,k)-yc(j-1,k+1))**2))*
     &  2*(yc(j-1,k)-yo(j,k))
      zp_g3=zp(j-1,k)+(zp(j-1,k)-zp(j-2,k))/
     &  (sqrt((xc(j-1,k)-xc(j-2,k))**2+(yc(j-1,k)-yc(j-2,k))**2))
     &  *2*(xc(j-1,k)-xo(j,k))
      zo(j,k)=0.25*(zp(j-1,k)+zp_g1+zp_g2+zp_g3)
      elseif (j.EQ.M.AND.k.EQ.N) then
      zp_g1=zp(j-1,k-1)+(zp(j-1,k-1)-zp(j-2,k-2))/
     &(sqrt((xc(j-1,k-1)-xc(j-2,k-2))**2+(yc(j-1,k-1)-yc(j-2,k-2))**2))*
     &  2*(sqrt((xc(j-1,k-1)-xo(j,k))**2+(yc(j-1,k-1)-yo(j,k))**2))
      zp_g2=zp(j-1,k-1)+(zp(j-1,k-1)-zp(j-1,k-2))/
     &(sqrt((xc(j-1,k-1)-xc(j-1,k-2))**2+(yc(j-1,k-1)-yc(j-1,k-2))**2))*
     &  2*(yc(j-1,k-1)-yo(j,k))
      zp_g3=zp(j-1,k-1)+(zp(j-1,k-1)-zp(j-2,k-1))/
     & (sqrt((xc(j-1,k-1)-xc(j-1,k-2))**2+(yc(j-1,k-1)-yc(j-1,k-2))**2))
     &  *2*(xc(j-1,k-1)-xo(j,k))
      zo(j,k)=0.25*(zp(j-1,k-1)+zp_g1+zp_g2+zp_g3)
      elseif (j.EQ.M.AND.k.GT.1.AND.k.LT.N) then
      zp_g1=zp(j-1,k-1)+(zp(j-1,k-1)-zp(j-2,k-1))/
     & (sqrt((xc(j-1,k-1)-xc(j-2,k-1))**2+(yc(j-1,k-1)-yc(j-2,k-1))**2))
     &  *2*(xc(j-1,k-1)-xo(j,k))
      zp_g1=zp(j-1,k)+(zp(j-1,k)-zp(j-2,k))/
     &  (sqrt((xc(j-1,k)-xc(j-2,k))**2+(yc(j-1,k)-yc(j-2,k))**2))
     &  *2*(xc(j-1,k)-xo(j,k))
      zo(j,k)=0.25*(zp_g1+zp_g2+zp(j-1,k-1)+zp(j-1,k))
      else
      zo(j,k)=(zp(j-1,k-1)*
     &         (sqrt((xc(j-1,k-1)-xo(j,k))**2+(yc(j-1,k-1)-yo(j,k))**2))
     &  +zp(j-1,k)*(sqrt((xc(j-1,k)-xo(j,k))**2+(yc(j-1,k)-yo(j,k))**2))
     &  +zp(j,k-1)*(sqrt((xc(j,k-1)-xo(j,k))**2+(yc(j,k-1)-yo(j,k))**2))
     &  +zp(j,k)*(sqrt((xc(j,k)-xo(j,k))**2+(yc(j,k)-yo(j,k))**2)))/
     &  (sqrt((xc(j-1,k-1)-xo(j,k))**2+(yc(j-1,k-1)-yo(j,k))**2)+
     &  sqrt((xc(j-1,k)-xo(j,k))**2+(yc(j-1,k)-yo(j,k))**2)+
     &  sqrt((xc(j,k-1)-xo(j,k))**2+(yc(j,k-1)-yo(j,k))**2)+
     &  sqrt((xc(j,k)-xo(j,k))**2+(yc(j,k)-yo(j,k))**2))
      endif
 96   continue
 95   continue

******************************************************************************************************************************** 
******************************************************************************************************************************
*     call new boundary conditions and intial condition for the new time step
      if (H(MM,1).EQ.0.0) then
      H_cri=0.9*((Hinlet*Uinlet)**2.0/9.81)**0.333333333
      else
      H_cri=0.9*((H(MM,1)*U(MM,1))**2.0/9.81)**0.333333333
      endif
      if (Boundary_condition_option.EQ.1) then
      goto 1500
      else
      do 91 k=1,NN
      if (BC2.EQ.3) then
      Hgo(1,k)=H(MM,k)
      elseif(H(MM,k).LE.H_cri.AND.BC.EQ.4) then
      Hgo(1,k)=H(MM,k)
      else
      Hgo(1,k)=Houtlet
      endif
      Ugo(1,k)=U(MM,k)
      Vgo(1,k)=V(MM,k)
      Cgo(1,k)=C(MM,k)
   91 continue
   !  Wall BCs          the net perpendicular velocity to the wall is zero
      do 92 j=1,MM                      
      Hgw(j,1)=H(j,1)                                       !For the lower wall
      Cgw(j,1)=C(j,1)                                       !For the lower wall
      Ugw(j,1)=U(j,1)*((sphi(j,1,1))**2-(cphi(j,1,1))**2)-  !For the lower wall
     &          2.0*V(j,1)*sphi(j,1,1)*cphi(j,1,1)
      Vgw(j,1)=-2.0*U(j,1)*sphi(j,1,1)*cphi(j,1,1)+         !For the lower wall
     &              +V(j,1)*((cphi(j,1,1))**2-(sphi(j,1,1))**2)
      Hgw(j,2)=H(j,NN)                                       !For the Upper wall
      Cgw(j,2)=C(j,NN)                                       !For the Upper wall
      Ugw(j,2)=U(j,NN)*((sphi(j,NN,3))**2-(cphi(j,NN,3))**2)-  !For the Upper wall
     &          2.0*V(j,NN)*sphi(j,NN,3)*cphi(j,NN,3)
      Vgw(j,2)=-2.0*U(j,NN)*sphi(j,NN,3)*cphi(j,NN,3)+         !For the Upper wall
     &              +V(j,NN)*((cphi(j,NN,3))**2-(sphi(j,NN,3))**2)
   92 continue
      endif
******************************************************************************************************************************
*     Update the values at the cells corners
      ! Adding some values in the new four ghost cells at the corners
*    ===============================================================
      ! ***AT THE ENTRANCE***(location,wall number,time)
1501  Hgiw(1,1,t)=Hgi(1,1)
      Hgiw(1,2,t)=Hgi(1,NN)
      Cgiw(1,1,t)=Cgi(1,1)
      Cgiw(1,2,t)=Cgi(1,NN)
      Ugiw(1,1,t)=Ugi(1,1)*((sphi(1,1,1))**2-(cphi(1,1,1))**2)-  !For the lower wall
     &          2.0*Vgi(1,1)*sphi(1,1,1)*cphi(1,1,1)
      Vgiw(1,1,t)=-2.0*Ugi(1,1)*sphi(1,1,1)*cphi(1,1,1)+         !For the lower wall
     &              +Vgi(1,1)*((cphi(1,1,1))**2-(sphi(1,1,1))**2)
      Ugiw(1,2,t)=Ugi(1,NN)*((sphi(1,NN,3))**2-(cphi(1,NN,3))**2)-  !For the upper wall
     &          2.0*Vgi(1,NN)*sphi(1,NN,3)*cphi(1,NN,3)
      Vgiw(1,2,t)=-2.0*Ugi(1,NN)*sphi(1,NN,3)*cphi(1,NN,3)+         !For the upper wall
     &              +Vgi(1,NN)*((cphi(1,NN,3))**2-(sphi(1,NN,3))**2)
      ! ***AT THE EXIT***
      Hgiw(2,1,t)=Hgo(1,1)
      Hgiw(2,2,t)=Hgo(1,NN)
      Cgiw(2,1,t)=Cgo(1,1)
      Cgiw(2,2,t)=Cgo(1,NN)
      Ugiw(2,1,t)=Ugo(1,1)*((sphi(MM,1,1))**2-(cphi(MM,1,1))**2)-  !For the lower wall
     &          2.0*Vgo(1,1)*sphi(MM,1,1)*cphi(MM,1,1)
      Vgiw(2,1,t)=-2.0*Ugo(1,1)*sphi(MM,1,1)*cphi(MM,1,1)+         !For the lower wall
     &              +Vgo(1,1)*((cphi(MM,1,1))**2-(sphi(MM,1,1))**2)
      Ugiw(2,2,t)=Ugo(1,NN)*((sphi(MM,NN,3))**2-(cphi(MM,NN,3))**2)-  !For the upper wall
     &          2.0*Vgo(1,NN)*sphi(MM,NN,3)*cphi(MM,NN,3)
      Vgiw(2,2,t)=-2.0*Ugo(1,NN)*sphi(MM,NN,3)*cphi(MM,NN,3)+         !For the upper wall
     &              +Vgo(1,NN)*((cphi(MM,NN,3))**2-(sphi(MM,NN,3))**2)
*     Coordinats of the four ghost corners(xcgiw,ycgiw)
      xcgiw(1,1)=xcgi(1,1)
      ycgiw(1,1)=ycgw(1,1)
      xcgiw(1,2)=xcgi(1,NN)
      ycgiw(1,2)=ycgw(1,2)
      xcgiw(2,1)=xcgo(1,1)
      ycgiw(2,1)=ycgw(MM,1)
      xcgiw(2,2)=xcgo(1,1)
      ycgiw(2,2)=ycgw(MM,2)
*     Interpolation for H
      do 500 j=1,M
      do 501 k=1,N
      if ((k.EQ.1).AND.(j.EQ.1)) then
      Hc(j,k)=(Hgiw(1,1,t)*sqrt((xo(j,k)-xcgiw(1,1))**2+
     &          (yo(j,k)-ycgiw(1,1))**2)+Hgi(1,k)
     &         *sqrt((xo(j,k)-xcgi(1,k))**2+(yo(j,k)-ycgi(1,k))**2)+
     &          H(j,k)*sqrt((xo(j,k)-xc(j,k))**2+(yo(j,k)-yc(j,k))**2)
     & +Hgw(j,1)*sqrt((xo(j,k)-xcgw(j,1))**2+(yo(j,k)-ycgw(j,1))**2))/
     &      (sqrt((xo(j,k)-xcgiw(1,1))**2+(yo(j,k)-ycgiw(1,1))**2)+
     &      sqrt((xo(j,k)-xcgi(1,k))**2+(yo(j,k)-ycgi(1,k))**2)+
     &      sqrt((xo(j,k)-xc(j,k))**2+(yo(j,k)-yc(j,k))**2)+
     &      sqrt((xo(j,k)-xcgw(j,1))**2+(yo(j,k)-ycgw(j,1))**2))
     
      Uc(j,k)=(Ugiw(1,1,t)*sqrt((xo(j,k)-xcgiw(1,1))**2+
     &          (yo(j,k)-ycgiw(1,1))**2)+Ugi(1,k)
     &         *sqrt((xo(j,k)-xcgi(1,k))**2+(yo(j,k)-ycgi(1,k))**2)+
     &          U(j,k)*sqrt((xo(j,k)-xc(j,k))**2+(yo(j,k)-yc(j,k))**2)
     & +Ugw(j,1)*sqrt((xo(j,k)-xcgw(j,1))**2+(yo(j,k)-ycgw(j,1))**2))/
     &      (sqrt((xo(j,k)-xcgiw(1,1))**2+(yo(j,k)-ycgiw(1,1))**2)+
     &      sqrt((xo(j,k)-xcgi(1,k))**2+(yo(j,k)-ycgi(1,k))**2)+
     &      sqrt((xo(j,k)-xc(j,k))**2+(yo(j,k)-yc(j,k))**2)+
     &      sqrt((xo(j,k)-xcgw(j,1))**2+(yo(j,k)-ycgw(j,1))**2))
     
      Vc(j,k)=(Vgiw(1,1,t)*sqrt((xo(j,k)-xcgiw(1,1))**2+
     &          (yo(j,k)-ycgiw(1,1))**2)+Vgi(1,k)
     &         *sqrt((xo(j,k)-xcgi(1,k))**2+(yo(j,k)-ycgi(1,k))**2)+
     &          V(j,k)*sqrt((xo(j,k)-xc(j,k))**2+(yo(j,k)-yc(j,k))**2)
     & +Vgw(j,1)*sqrt((xo(j,k)-xcgw(j,1))**2+(yo(j,k)-ycgw(j,1))**2))/
     &      (sqrt((xo(j,k)-xcgiw(1,1))**2+(yo(j,k)-ycgiw(1,1))**2)+
     &      sqrt((xo(j,k)-xcgi(1,k))**2+(yo(j,k)-ycgi(1,k))**2)+
     &      sqrt((xo(j,k)-xc(j,k))**2+(yo(j,k)-yc(j,k))**2)+
     &      sqrt((xo(j,k)-xcgw(j,1))**2+(yo(j,k)-ycgw(j,1))**2))
     
      Cc(j,k)=C(j,k)
     
      elseif ((k.EQ.N).AND.(j.EQ.1)) then
      Hc(j,k)=(Hgiw(1,2,t)*sqrt((xo(j,k)-xcgiw(1,2))**2+
     &          (yo(j,k)-ycgiw(1,2))**2)+Hgi(1,k-1)
     &         *sqrt((xo(j,k)-xcgi(1,k-1))**2+(yo(j,k)-ycgi(1,k-1))**2)+
     &  H(j,k-1)*sqrt((xo(j,k)-xc(j,k-1))**2+(yo(j,k)-yc(j,k-1))**2)
     & +Hgw(j,2)*sqrt((xo(j,k)-xcgw(j,2))**2+(yo(j,k)-ycgw(j,2))**2))/
     & (sqrt((xo(j,k)-xcgiw(1,2))**2+(yo(j,k)-ycgiw(1,2))**2)+
     & sqrt((xo(j,k)-xcgi(1,k-1))**2+(yo(j,k)-ycgi(1,k-1))**2)+
     & sqrt((xo(j,k)-xc(j,k-1))**2+(yo(j,k)-yc(j,k-1))**2)+
     & sqrt((xo(j,k)-xcgw(j,2))**2+(yo(j,k)-ycgw(j,2))**2))
     
      Uc(j,k)=(Ugiw(1,2,t)*sqrt((xo(j,k)-xcgiw(1,2))**2+
     &          (yo(j,k)-ycgiw(1,2))**2)+Ugi(1,k-1)
     &         *sqrt((xo(j,k)-xcgi(1,k-1))**2+(yo(j,k)-ycgi(1,k-1))**2)+
     &  U(j,k-1)*sqrt((xo(j,k)-xc(j,k-1))**2+(yo(j,k)-yc(j,k-1))**2)
     & +Ugw(j,2)*sqrt((xo(j,k)-xcgw(j,2))**2+(yo(j,k)-ycgw(j,2))**2))/
     & (sqrt((xo(j,k)-xcgiw(1,2))**2+(yo(j,k)-ycgiw(1,2))**2)+
     & sqrt((xo(j,k)-xcgi(1,k-1))**2+(yo(j,k)-ycgi(1,k-1))**2)+
     & sqrt((xo(j,k)-xc(j,k-1))**2+(yo(j,k)-yc(j,k-1))**2)+
     & sqrt((xo(j,k)-xcgw(j,2))**2+(yo(j,k)-ycgw(j,2))**2))
     
      Vc(j,k)=(Vgiw(1,2,t)*sqrt((xo(j,k)-xcgiw(1,2))**2+
     &          (yo(j,k)-ycgiw(1,2))**2)+Vgi(1,k-1)
     &         *sqrt((xo(j,k)-xcgi(1,k-1))**2+(yo(j,k)-ycgi(1,k-1))**2)+
     &  V(j,k-1)*sqrt((xo(j,k)-xc(j,k-1))**2+(yo(j,k)-yc(j,k-1))**2)
     & +Vgw(j,2)*sqrt((xo(j,k)-xcgw(j,2))**2+(yo(j,k)-ycgw(j,2))**2))/
     & (sqrt((xo(j,k)-xcgiw(1,2))**2+(yo(j,k)-ycgiw(1,2))**2)+
     & sqrt((xo(j,k)-xcgi(1,k-1))**2+(yo(j,k)-ycgi(1,k-1))**2)+
     & sqrt((xo(j,k)-xc(j,k-1))**2+(yo(j,k)-yc(j,k-1))**2)+
     & sqrt((xo(j,k)-xcgw(j,2))**2+(yo(j,k)-ycgw(j,2))**2))
     
      Cc(j,k)=C(j,k-1)
     
      elseif ((k.EQ.1).AND.(j.EQ.M)) then
      Hc(j,k)=(Hgiw(2,1,t)*sqrt((xo(j,k)-xcgiw(2,1))**2+
     &          (yo(j,k)-ycgiw(2,1))**2)+Hgo(1,k)
     &         *sqrt((xo(j,k)-xcgo(1,k))**2+(yo(j,k)-ycgo(1,k))**2)+
     & H(j-1,k)*sqrt((xo(j,k)-xc(j-1,k))**2+
     &         (yo(j,k)-yc(j-1,k))**2)+Hgw(j-1,1)*
     &         sqrt((xo(j,k)-xcgw(j-1,1))**2+(yo(j,k)-ycgw(j-1,1))**2))/
     & (sqrt((xo(j,k)-xcgiw(2,1))**2+(yo(j,k)-ycgiw(2,1))**2)+
     & sqrt((xo(j,k)-xcgo(1,k))**2+(yo(j,k)-ycgo(1,k))**2)+
     & sqrt((xo(j,k)-xc(j-1,k))**2+(yo(j,k)-yc(j-1,k))**2)+
     & sqrt((xo(j,k)-xcgw(j-1,1))**2+(yo(j,k)-ycgw(j-1,1))**2))
     
      Uc(j,k)=(Ugiw(2,1,t)*sqrt((xo(j,k)-xcgiw(2,1))**2+
     &          (yo(j,k)-ycgiw(2,1))**2)+Ugo(1,k)
     &         *sqrt((xo(j,k)-xcgo(1,k))**2+(yo(j,k)-ycgo(1,k))**2)+
     & U(j-1,k)*sqrt((xo(j,k)-xc(j-1,k))**2+
     &         (yo(j,k)-yc(j-1,k))**2)+Ugw(j-1,1)*
     &         sqrt((xo(j,k)-xcgw(j-1,1))**2+(yo(j,k)-ycgw(j-1,1))**2))/
     & (sqrt((xo(j,k)-xcgiw(2,1))**2+(yo(j,k)-ycgiw(2,1))**2)+
     & sqrt((xo(j,k)-xcgo(1,k))**2+(yo(j,k)-ycgo(1,k))**2)+
     & sqrt((xo(j,k)-xc(j-1,k))**2+(yo(j,k)-yc(j-1,k))**2)+
     & sqrt((xo(j,k)-xcgw(j-1,1))**2+(yo(j,k)-ycgw(j-1,1))**2))
     
      Vc(j,k)=(Vgiw(2,1,t)*sqrt((xo(j,k)-xcgiw(2,1))**2+
     &          (yo(j,k)-ycgiw(2,1))**2)+Vgo(1,k)
     &         *sqrt((xo(j,k)-xcgo(1,k))**2+(yo(j,k)-ycgo(1,k))**2)+
     & V(j-1,k)*sqrt((xo(j,k)-xc(j-1,k))**2+
     &         (yo(j,k)-yc(j-1,k))**2)+Vgw(j-1,1)*
     &         sqrt((xo(j,k)-xcgw(j-1,1))**2+(yo(j,k)-ycgw(j-1,1))**2))/
     & (sqrt((xo(j,k)-xcgiw(2,1))**2+(yo(j,k)-ycgiw(2,1))**2)+
     & sqrt((xo(j,k)-xcgo(1,k))**2+(yo(j,k)-ycgo(1,k))**2)+
     & sqrt((xo(j,k)-xc(j-1,k))**2+(yo(j,k)-yc(j-1,k))**2)+
     & sqrt((xo(j,k)-xcgw(j-1,1))**2+(yo(j,k)-ycgw(j-1,1))**2))
     
      Cc(j,k)=(Cgiw(2,1,t)*sqrt((xo(j,k)-xcgiw(2,1))**2+
     &          (yo(j,k)-ycgiw(2,1))**2)+Cgo(1,k)
     &         *sqrt((xo(j,k)-xcgo(1,k))**2+(yo(j,k)-ycgo(1,k))**2)+
     & C(j-1,k)*sqrt((xo(j,k)-xc(j-1,k))**2+
     &         (yo(j,k)-yc(j-1,k))**2)+Cgw(j-1,1)*
     &         sqrt((xo(j,k)-xcgw(j-1,1))**2+(yo(j,k)-ycgw(j-1,1))**2))/
     & (sqrt((xo(j,k)-xcgiw(2,1))**2+(yo(j,k)-ycgiw(2,1))**2)+
     & sqrt((xo(j,k)-xcgo(1,k))**2+(yo(j,k)-ycgo(1,k))**2)+
     & sqrt((xo(j,k)-xc(j-1,k))**2+(yo(j,k)-yc(j-1,k))**2)+
     & sqrt((xo(j,k)-xcgw(j-1,1))**2+(yo(j,k)-ycgw(j-1,1))**2))
     
      elseif ((k.EQ.N).AND.(j.EQ.M)) then
      Hc(j,k)=(Hgiw(2,2,t)*sqrt((xo(j,k)-xcgiw(2,2))**2+
     &          (yo(j,k)-ycgiw(2,2))**2)+Hgo(1,k-1)
     &         *sqrt((xo(j,k)-xcgo(1,k-1))**2+(yo(j,k)-ycgo(1,k-1))**2)+
     &          H(j-1,k-1)*sqrt((xo(j,k)-xc(j-1,k-1))**2+
     &          (yo(j,k)-yc(j-1,k-1))**2)+Hgw(j-1,2)*
     &         sqrt((xo(j,k)-xcgw(j-1,2))**2+(yo(j,k)-ycgw(j-1,2))**2))/
     &  (sqrt((xo(j,k)-xcgiw(2,2))**2+(yo(j,k)-ycgiw(2,2))**2)+
     &  sqrt((xo(j,k)-xcgo(1,k-1))**2+(yo(j,k)-ycgo(1,k-1))**2)+
     &  sqrt((xo(j,k)-xc(j-1,k-1))**2+(yo(j,k)-yc(j-1,k-1))**2)+
     &  sqrt((xo(j,k)-xcgw(j-1,2))**2+(yo(j,k)-ycgw(j-1,2))**2))
     
      Uc(j,k)=(Ugiw(2,2,t)*sqrt((xo(j,k)-xcgiw(2,2))**2+
     &          (yo(j,k)-ycgiw(2,2))**2)+Ugo(1,k-1)
     &         *sqrt((xo(j,k)-xcgo(1,k-1))**2+(yo(j,k)-ycgo(1,k-1))**2)+
     &          U(j-1,k-1)*sqrt((xo(j,k)-xc(j-1,k-1))**2+
     &          (yo(j,k)-yc(j-1,k-1))**2)+Ugw(j-1,2)*
     &         sqrt((xo(j,k)-xcgw(j-1,2))**2+(yo(j,k)-ycgw(j-1,2))**2))/
     &  (sqrt((xo(j,k)-xcgiw(2,2))**2+(yo(j,k)-ycgiw(2,2))**2)+
     &  sqrt((xo(j,k)-xcgo(1,k-1))**2+(yo(j,k)-ycgo(1,k-1))**2)+
     &  sqrt((xo(j,k)-xc(j-1,k-1))**2+(yo(j,k)-yc(j-1,k-1))**2)+
     &  sqrt((xo(j,k)-xcgw(j-1,2))**2+(yo(j,k)-ycgw(j-1,2))**2))
     
      Vc(j,k)=(Hgiw(2,2,t)*sqrt((xo(j,k)-xcgiw(2,2))**2+
     &          (yo(j,k)-ycgiw(2,2))**2)+Vgo(1,k-1)
     &         *sqrt((xo(j,k)-xcgo(1,k-1))**2+(yo(j,k)-ycgo(1,k-1))**2)+
     &          V(j-1,k-1)*sqrt((xo(j,k)-xc(j-1,k-1))**2+
     &          (yo(j,k)-yc(j-1,k-1))**2)+Vgw(j-1,2)*
     &         sqrt((xo(j,k)-xcgw(j-1,2))**2+(yo(j,k)-ycgw(j-1,2))**2))/
     &  (sqrt((xo(j,k)-xcgiw(2,2))**2+(yo(j,k)-ycgiw(2,2))**2)+
     &  sqrt((xo(j,k)-xcgo(1,k-1))**2+(yo(j,k)-ycgo(1,k-1))**2)+
     &  sqrt((xo(j,k)-xc(j-1,k-1))**2+(yo(j,k)-yc(j-1,k-1))**2)+
     &  sqrt((xo(j,k)-xcgw(j-1,2))**2+(yo(j,k)-ycgw(j-1,2))**2))
     
      Cc(j,k)=(Cgiw(2,2,t)*sqrt((xo(j,k)-xcgiw(2,2))**2+
     &          (yo(j,k)-ycgiw(2,2))**2)+Cgo(1,k-1)
     &         *sqrt((xo(j,k)-xcgo(1,k-1))**2+(yo(j,k)-ycgo(1,k-1))**2)+
     &          C(j-1,k-1)*sqrt((xo(j,k)-xc(j-1,k-1))**2+
     &          (yo(j,k)-yc(j-1,k-1))**2)+Cgw(j-1,2)*
     &         sqrt((xo(j,k)-xcgw(j-1,2))**2+(yo(j,k)-ycgw(j-1,2))**2))/
     &  (sqrt((xo(j,k)-xcgiw(2,2))**2+(yo(j,k)-ycgiw(2,2))**2)+
     &  sqrt((xo(j,k)-xcgo(1,k-1))**2+(yo(j,k)-ycgo(1,k-1))**2)+
     &  sqrt((xo(j,k)-xc(j-1,k-1))**2+(yo(j,k)-yc(j-1,k-1))**2)+
     &  sqrt((xo(j,k)-xcgw(j-1,2))**2+(yo(j,k)-ycgw(j-1,2))**2))
     
      elseif ((j.EQ.1).AND.(k.LT.N).AND.(k.GT.1)) then                          !First column
      Hc(j,k)=(Hgi(1,k-1)*sqrt((xo(j,k)-xcgi(1,k-1))**2+
     &          (yo(j,k)-ycgi(1,k-1))**2)+Hgi(1,k)
     &         *sqrt((xo(j,k)-xcgi(1,k))**2+(yo(j,k)-ycgi(1,k))**2)+
     &          H(j,k)*sqrt((xo(j,k)-xc(j,k))**2+(yo(j,k)-yc(j,k))**2)
     & +H(j,k-1)*sqrt((xo(j,k)-xc(j,k-1))**2+(yo(j,k)-yc(j,k-1))**2))/
     & (sqrt((xo(j,k)-xcgi(1,k-1))**2+(yo(j,k)-ycgi(1,k-1))**2)+
     &  sqrt((xo(j,k)-xcgi(1,k))**2+(yo(j,k)-ycgi(1,k))**2)+
     &  sqrt((xo(j,k)-xc(j,k))**2+(yo(j,k)-yc(j,k))**2)+
     &  sqrt((xo(j,k)-xc(j,k-1))**2+(yo(j,k)-yc(j,k-1))**2))
     
      Uc(j,k)=(Ugi(1,k-1)*sqrt((xo(j,k)-xcgi(1,k-1))**2+
     &          (yo(j,k)-ycgi(1,k-1))**2)+Ugi(1,k)
     &         *sqrt((xo(j,k)-xcgi(1,k))**2+(yo(j,k)-ycgi(1,k))**2)+
     &          U(j,k)*sqrt((xo(j,k)-xc(j,k))**2+(yo(j,k)-yc(j,k))**2)
     & +U(j,k-1)*sqrt((xo(j,k)-xc(j,k-1))**2+(yo(j,k)-yc(j,k-1))**2))/
     & (sqrt((xo(j,k)-xcgi(1,k-1))**2+(yo(j,k)-ycgi(1,k-1))**2)+
     &  sqrt((xo(j,k)-xcgi(1,k))**2+(yo(j,k)-ycgi(1,k))**2)+
     &  sqrt((xo(j,k)-xc(j,k))**2+(yo(j,k)-yc(j,k))**2)+
     &  sqrt((xo(j,k)-xc(j,k-1))**2+(yo(j,k)-yc(j,k-1))**2))
     
      Vc(j,k)=(Vgi(1,k-1)*sqrt((xo(j,k)-xcgi(1,k-1))**2+
     &          (yo(j,k)-ycgi(1,k-1))**2)+Vgi(1,k)
     &         *sqrt((xo(j,k)-xcgi(1,k))**2+(yo(j,k)-ycgi(1,k))**2)+
     &          V(j,k)*sqrt((xo(j,k)-xc(j,k))**2+(yo(j,k)-yc(j,k))**2)
     & +V(j,k-1)*sqrt((xo(j,k)-xc(j,k-1))**2+(yo(j,k)-yc(j,k-1))**2))/
     & (sqrt((xo(j,k)-xcgi(1,k-1))**2+(yo(j,k)-ycgi(1,k-1))**2)+
     &  sqrt((xo(j,k)-xcgi(1,k))**2+(yo(j,k)-ycgi(1,k))**2)+
     &  sqrt((xo(j,k)-xc(j,k))**2+(yo(j,k)-yc(j,k))**2)+
     &  sqrt((xo(j,k)-xc(j,k-1))**2+(yo(j,k)-yc(j,k-1))**2))
     
      Cc(j,k)=C(j,k)
     
      elseif ((j.EQ.M).AND.(k.LT.N).AND.(k.GT.1)) then                          !Last column
      Hc(j,k)=(Hgo(1,k-1)*sqrt((xo(j,k)-xcgo(1,k-1))**2+
     &  (yo(j,k)-ycgo(1,k-1))**2)+Hgo(1,k)*sqrt((xo(j,k)-xcgo(1,k))**2
     & +(yo(j,k)-ycgo(1,k))**2)+H(j-1,k)*sqrt((xo(j,k)-xc(j-1,k))**2+
     &(yo(j,k)-yc(j-1,k))**2)+H(j-1,k-1)*sqrt((xo(j,k)-xc(j-1,k-1))**2
     &  +(yo(j,k)-yc(j-1,k-1))**2))/(sqrt((xo(j,k)-xcgo(1,k-1))**2+
     &  (yo(j,k)-ycgo(1,k-1))**2)+sqrt((xo(j,k)-xcgo(1,k))**2
     & +(yo(j,k)-ycgo(1,k))**2)+sqrt((xo(j,k)-xc(j-1,k))**2+
     &(yo(j,k)-yc(j-1,k))**2)+sqrt((xo(j,k)-xc(j-1,k-1))**2
     &  +(yo(j,k)-yc(j-1,k-1))**2))
     
      Uc(j,k)=(Ugo(1,k-1)*sqrt((xo(j,k)-xcgo(1,k-1))**2+
     &  (yo(j,k)-ycgo(1,k-1))**2)+Ugo(1,k)*sqrt((xo(j,k)-xcgo(1,k))**2
     & +(yo(j,k)-ycgo(1,k))**2)+U(j-1,k)*sqrt((xo(j,k)-xc(j-1,k))**2+
     &(yo(j,k)-yc(j-1,k))**2)+U(j-1,k-1)*sqrt((xo(j,k)-xc(j-1,k-1))**2
     &  +(yo(j,k)-yc(j-1,k-1))**2))/(sqrt((xo(j,k)-xcgo(1,k-1))**2+
     &  (yo(j,k)-ycgo(1,k-1))**2)+sqrt((xo(j,k)-xcgo(1,k))**2
     & +(yo(j,k)-ycgo(1,k))**2)+sqrt((xo(j,k)-xc(j-1,k))**2+
     &(yo(j,k)-yc(j-1,k))**2)+sqrt((xo(j,k)-xc(j-1,k-1))**2
     &  +(yo(j,k)-yc(j-1,k-1))**2))
     
      Vc(j,k)=(Vgo(1,k-1)*sqrt((xo(j,k)-xcgo(1,k-1))**2+
     &  (yo(j,k)-ycgo(1,k-1))**2)+Vgo(1,k)*sqrt((xo(j,k)-xcgo(1,k))**2
     & +(yo(j,k)-ycgo(1,k))**2)+V(j-1,k)*sqrt((xo(j,k)-xc(j-1,k))**2+
     &(yo(j,k)-yc(j-1,k))**2)+V(j-1,k-1)*sqrt((xo(j,k)-xc(j-1,k-1))**2
     &  +(yo(j,k)-yc(j-1,k-1))**2))/(sqrt((xo(j,k)-xcgo(1,k-1))**2+
     &  (yo(j,k)-ycgo(1,k-1))**2)+sqrt((xo(j,k)-xcgo(1,k))**2
     & +(yo(j,k)-ycgo(1,k))**2)+sqrt((xo(j,k)-xc(j-1,k))**2+
     &(yo(j,k)-yc(j-1,k))**2)+sqrt((xo(j,k)-xc(j-1,k-1))**2
     &  +(yo(j,k)-yc(j-1,k-1))**2))
     
      Cc(j,k)=(Cgo(1,k-1)*sqrt((xo(j,k)-xcgo(1,k-1))**2+
     &  (yo(j,k)-ycgo(1,k-1))**2)+Cgo(1,k)*sqrt((xo(j,k)-xcgo(1,k))**2
     & +(yo(j,k)-ycgo(1,k))**2)+C(j-1,k)*sqrt((xo(j,k)-xc(j-1,k))**2+
     &(yo(j,k)-yc(j-1,k))**2)+C(j-1,k-1)*sqrt((xo(j,k)-xc(j-1,k-1))**2
     &  +(yo(j,k)-yc(j-1,k-1))**2))/(sqrt((xo(j,k)-xcgo(1,k-1))**2+
     &  (yo(j,k)-ycgo(1,k-1))**2)+sqrt((xo(j,k)-xcgo(1,k))**2
     & +(yo(j,k)-ycgo(1,k))**2)+sqrt((xo(j,k)-xc(j-1,k))**2+
     &(yo(j,k)-yc(j-1,k))**2)+sqrt((xo(j,k)-xc(j-1,k-1))**2
     &  +(yo(j,k)-yc(j-1,k-1))**2))
     
      elseif ((k.EQ.1).AND.(j.LT.M).AND.(j.GT.1)) then                          !right wall
      Hc(j,k)=(Hgw(j,1)*sqrt((xo(j,k)-xcgw(j,1))**2+
     &(yo(j,k)-ycgw(j,1))**2)+Hgw(j-1,1)*sqrt((xo(j,k)-xcgw(j-1,1))**2
     &+(yo(j,k)-ycgw(j-1,1))**2)+H(j-1,k)*sqrt((xo(j,k)-xc(j-1,k))**2+
     &(yo(j,k)-yc(j-1,k))**2)+H(j,k)*sqrt((xo(j,k)-xc(j,k))**2
     & +(yo(j,k)-yc(j,k))**2))/(sqrt((xo(j,k)-xcgw(j,1))**2+
     & (yo(j,k)-ycgw(j,1))**2)+sqrt((xo(j,k)-xcgw(j-1,1))**2
     &+(yo(j,k)-ycgw(j-1,1))**2)+sqrt((xo(j,k)-xc(j-1,k))**2+
     &(yo(j,k)-yc(j-1,k))**2)+sqrt((xo(j,k)-xc(j,k))**2
     &+(yo(j,k)-yc(j,k))**2))
     
      Uc(j,k)=(Ugw(j,1)*sqrt((xo(j,k)-xcgw(j,1))**2+
     &(yo(j,k)-ycgw(j,1))**2)+Ugw(j-1,1)*sqrt((xo(j,k)-xcgw(j-1,1))**2
     &+(yo(j,k)-ycgw(j-1,1))**2)+U(j-1,k)*sqrt((xo(j,k)-xc(j-1,k))**2+
     &(yo(j,k)-yc(j-1,k))**2)+U(j,k)*sqrt((xo(j,k)-xc(j,k))**2
     &+(yo(j,k)-yc(j,k))**2))/(sqrt((xo(j,k)-xcgw(j,1))**2+
     & (yo(j,k)-ycgw(j,1))**2)+sqrt((xo(j,k)-xcgw(j-1,1))**2
     &+(yo(j,k)-ycgw(j-1,1))**2)+sqrt((xo(j,k)-xc(j-1,k))**2+
     &(yo(j,k)-yc(j-1,k))**2)+sqrt((xo(j,k)-xc(j,k))**2
     &+(yo(j,k)-yc(j,k))**2))
     
      Vc(j,k)=(Vgw(j,1)*sqrt((xo(j,k)-xcgw(j,1))**2+
     &(yo(j,k)-ycgw(j,1))**2)+Vgw(j-1,1)*sqrt((xo(j,k)-xcgw(j-1,1))**2
     &+(yo(j,k)-ycgw(j-1,1))**2)+V(j-1,k)*sqrt((xo(j,k)-xc(j-1,k))**2+
     &(yo(j,k)-yc(j-1,k))**2)+V(j,k)*sqrt((xo(j,k)-xc(j,k))**2
     &+(yo(j,k)-yc(j,k))**2))/(sqrt((xo(j,k)-xcgw(j,1))**2+
     & (yo(j,k)-ycgw(j,1))**2)+sqrt((xo(j,k)-xcgw(j-1,1))**2
     &+(yo(j,k)-ycgw(j-1,1))**2)+sqrt((xo(j,k)-xc(j-1,k))**2+
     &(yo(j,k)-yc(j-1,k))**2)+sqrt((xo(j,k)-xc(j,k))**2
     &+(yo(j,k)-yc(j,k))**2))
     
      Cc(j,k)=(Cgw(j,1)*sqrt((xo(j,k)-xcgw(j,1))**2+
     &(yo(j,k)-ycgw(j,1))**2)+Cgw(j-1,1)*sqrt((xo(j,k)-xcgw(j-1,1))**2
     &+(yo(j,k)-ycgw(j-1,1))**2)+C(j-1,k)*sqrt((xo(j,k)-xc(j-1,k))**2+
     &(yo(j,k)-yc(j-1,k))**2)+C(j,k)*sqrt((xo(j,k)-xc(j,k))**2
     &+(yo(j,k)-yc(j,k))**2))/(sqrt((xo(j,k)-xcgw(j,1))**2+
     & (yo(j,k)-ycgw(j,1))**2)+sqrt((xo(j,k)-xcgw(j-1,1))**2
     &+(yo(j,k)-ycgw(j-1,1))**2)+sqrt((xo(j,k)-xc(j-1,k))**2+
     &(yo(j,k)-yc(j-1,k))**2)+sqrt((xo(j,k)-xc(j,k))**2
     &+(yo(j,k)-yc(j,k))**2))
     
      elseif ((k.EQ.N).AND.(j.LT.M).AND.(j.GT.1)) then                          !left wall
      Hc(j,k)=(Hgw(j,2)*sqrt((xo(j,k)-xcgw(j,2))**2+
     &(yo(j,k)-ycgw(j,2))**2)+Hgw(j-1,2)*sqrt((xo(j,k)-xcgw(j-1,2))**2
     & +(yo(j,k)-ycgw(1,2))**2)+H(j,k-1)*sqrt((xo(j,k)-xc(j,k-1))**2+
     &(yo(j,k)-yc(j,k-1))**2)+H(j-1,k-1)*sqrt((xo(j,k)-xc(j-1,k-1))**2
     &+(yo(j,k)-yc(j-1,k-1))**2))/(sqrt((xo(j,k)-xcgw(j,2))**2+
     & (yo(j,k)-ycgw(j,2))**2)+sqrt((xo(j,k)-xcgw(j-1,2))**2
     & +(yo(j,k)-ycgw(1,2))**2)+sqrt((xo(j,k)-xc(j,k-1))**2+
     &(yo(j,k)-yc(j,k-1))**2)+sqrt((xo(j,k)-xc(j-1,k-1))**2
     &+(yo(j,k)-yc(j-1,k-1))**2))
     
      Uc(j,k)=(Ugw(j,2)*sqrt((xo(j,k)-xcgw(j,2))**2+
     &(yo(j,k)-ycgw(j,2))**2)+Ugw(j-1,2)*sqrt((xo(j,k)-xcgw(j-1,2))**2
     & +(yo(j,k)-ycgw(1,2))**2)+U(j,k-1)*sqrt((xo(j,k)-xc(j,k-1))**2+
     &(yo(j,k)-yc(j,k-1))**2)+U(j-1,k-1)*sqrt((xo(j,k)-xc(j-1,k-1))**2
     &+(yo(j,k)-yc(j-1,k-1))**2))/(sqrt((xo(j,k)-xcgw(j,2))**2+
     & (yo(j,k)-ycgw(j,2))**2)+sqrt((xo(j,k)-xcgw(j-1,2))**2
     & +(yo(j,k)-ycgw(1,2))**2)+sqrt((xo(j,k)-xc(j,k-1))**2+
     &(yo(j,k)-yc(j,k-1))**2)+sqrt((xo(j,k)-xc(j-1,k-1))**2
     &+(yo(j,k)-yc(j-1,k-1))**2))
     
      Vc(j,k)=(Vgw(j,2)*sqrt((xo(j,k)-xcgw(j,2))**2+
     &(yo(j,k)-ycgw(j,2))**2)+Vgw(j-1,2)*sqrt((xo(j,k)-xcgw(j-1,2))**2
     & +(yo(j,k)-ycgw(1,2))**2)+V(j,k-1)*sqrt((xo(j,k)-xc(j,k-1))**2+
     &(yo(j,k)-yc(j,k-1))**2)+V(j-1,k-1)*sqrt((xo(j,k)-xc(j-1,k-1))**2
     &+(yo(j,k)-yc(j-1,k-1))**2))/(sqrt((xo(j,k)-xcgw(j,2))**2+
     & (yo(j,k)-ycgw(j,2))**2)+sqrt((xo(j,k)-xcgw(j-1,2))**2
     & +(yo(j,k)-ycgw(1,2))**2)+sqrt((xo(j,k)-xc(j,k-1))**2+
     &(yo(j,k)-yc(j,k-1))**2)+sqrt((xo(j,k)-xc(j-1,k-1))**2
     &+(yo(j,k)-yc(j-1,k-1))**2))
     
      Cc(j,k)=(Cgw(j,2)*sqrt((xo(j,k)-xcgw(j,2))**2+
     &(yo(j,k)-ycgw(j,2))**2)+Cgw(j-1,2)*sqrt((xo(j,k)-xcgw(j-1,2))**2
     & +(yo(j,k)-ycgw(1,2))**2)+C(j,k-1)*sqrt((xo(j,k)-xc(j,k-1))**2+
     &(yo(j,k)-yc(j,k-1))**2)+C(j-1,k-1)*sqrt((xo(j,k)-xc(j-1,k-1))**2
     &+(yo(j,k)-yc(j-1,k-1))**2))/(sqrt((xo(j,k)-xcgw(j,2))**2+
     & (yo(j,k)-ycgw(j,2))**2)+sqrt((xo(j,k)-xcgw(j-1,2))**2
     & +(yo(j,k)-ycgw(1,2))**2)+sqrt((xo(j,k)-xc(j,k-1))**2+
     &(yo(j,k)-yc(j,k-1))**2)+sqrt((xo(j,k)-xc(j-1,k-1))**2
     &+(yo(j,k)-yc(j-1,k-1))**2))
     
      else
      Hc(j,k)=(H(j-1,k-1)*sqrt((xo(j,k)-xc(j-1,k-1))**2+
     &          (yo(j,k)-yc(j-1,k-1))**2)+H(j-1,k)
     &         *sqrt((xo(j,k)-xc(j-1,k))**2+(yo(j,k)-yc(j-1,k))**2)+
     & H(j,k-1)*sqrt((xo(j,k)-xc(j,k-1))**2+(yo(j,k)-yc(j,k-1))**2)+
     &  H(j,k)*sqrt((xo(j,k)-xc(j,k))**2+(yo(j,k)-yc(j,k))**2))/
     & (sqrt((xo(j,k)-xc(j-1,k-1))**2+(yo(j,k)-yc(j-1,k-1))**2)+
     &  sqrt((xo(j,k)-xc(j-1,k))**2+(yo(j,k)-yc(j-1,k))**2)+
     &  sqrt((xo(j,k)-xc(j,k-1))**2+(yo(j,k)-yc(j,k-1))**2)+
     &  sqrt((xo(j,k)-xc(j,k))**2+(yo(j,k)-yc(j,k))**2))
     
      Uc(j,k)=(U(j-1,k-1)*sqrt((xo(j,k)-xc(j-1,k-1))**2+
     &          (yo(j,k)-yc(j-1,k-1))**2)+U(j-1,k)
     &         *sqrt((xo(j,k)-xc(j-1,k))**2+(yo(j,k)-yc(j-1,k))**2)+
     & U(j,k-1)*sqrt((xo(j,k)-xc(j,k-1))**2+(yo(j,k)-yc(j,k-1))**2)+
     &  U(j,k)*sqrt((xo(j,k)-xc(j,k))**2+(yo(j,k)-yc(j,k))**2))/
     & (sqrt((xo(j,k)-xc(j-1,k-1))**2+(yo(j,k)-yc(j-1,k-1))**2)+
     &  sqrt((xo(j,k)-xc(j-1,k))**2+(yo(j,k)-yc(j-1,k))**2)+
     &  sqrt((xo(j,k)-xc(j,k-1))**2+(yo(j,k)-yc(j,k-1))**2)+
     &  sqrt((xo(j,k)-xc(j,k))**2+(yo(j,k)-yc(j,k))**2))
     
      Vc(j,k)=(V(j-1,k-1)*sqrt((xo(j,k)-xc(j-1,k-1))**2+
     &          (yo(j,k)-yc(j-1,k-1))**2)+V(j-1,k)
     &         *sqrt((xo(j,k)-xc(j-1,k))**2+(yo(j,k)-yc(j-1,k))**2)+
     & V(j,k-1)*sqrt((xo(j,k)-xc(j,k-1))**2+(yo(j,k)-yc(j,k-1))**2)+
     &  V(j,k)*sqrt((xo(j,k)-xc(j,k))**2+(yo(j,k)-yc(j,k))**2))/
     & (sqrt((xo(j,k)-xc(j-1,k-1))**2+(yo(j,k)-yc(j-1,k-1))**2)+
     &  sqrt((xo(j,k)-xc(j-1,k))**2+(yo(j,k)-yc(j-1,k))**2)+
     &  sqrt((xo(j,k)-xc(j,k-1))**2+(yo(j,k)-yc(j,k-1))**2)+
     &  sqrt((xo(j,k)-xc(j,k))**2+(yo(j,k)-yc(j,k))**2))
     
      Cc(j,k)=(C(j-1,k-1)*sqrt((xo(j,k)-xc(j-1,k-1))**2+
     &          (yo(j,k)-yc(j-1,k-1))**2)+C(j-1,k)
     &         *sqrt((xo(j,k)-xc(j-1,k))**2+(yo(j,k)-yc(j-1,k))**2)+
     & c(j,k-1)*sqrt((xo(j,k)-xc(j,k-1))**2+(yo(j,k)-yc(j,k-1))**2)+
     &  C(j,k)*sqrt((xo(j,k)-xc(j,k))**2+(yo(j,k)-yc(j,k))**2))/
     & (sqrt((xo(j,k)-xc(j-1,k-1))**2+(yo(j,k)-yc(j-1,k-1))**2)+
     &  sqrt((xo(j,k)-xc(j-1,k))**2+(yo(j,k)-yc(j-1,k))**2)+
     &  sqrt((xo(j,k)-xc(j,k-1))**2+(yo(j,k)-yc(j,k-1))**2)+
     &  sqrt((xo(j,k)-xc(j,k))**2+(yo(j,k)-yc(j,k))**2))
      endif
      Fr(j,k)=(sqrt((Uc(j,k))**2.0+(Vc(j,k))**2.0))/
     &          sqrt(9.81*Hc(j,k))
  501 continue
  500 continue
      if (t.EQ.1) then
      Hc_intial=Hc
      Uc_intial=Uc
      Vc_intial=Vc
      Cc_intial=Cc
      goto 1502
      endif
******************************************************************************************************************************
*     Checking the time step and CLF number
      do 93 j=1,1
      do 90 k=1,NN
      dxcell=0.5*(abs(xo(1+1,k)-xo(1,k))+abs(xo(1,k+1)-xo(1+1,k+1)))
      dycell=0.5*(abs(yo(1+1,k+1)-yo(1+1,k))+abs(yo(1,k)-yo(1,k+1)))
      c_no_inlet(1,k)=(abs(Ugi(1,k))+
     &              sqrt(9.81*Hgi(1,k)))/dxcell(j,k)+
     &  (abs(Vgi(1,k))+sqrt(9.81*Hgi(1,k)))/dycell(j,k)
 90   continue
 93   continue
      C_no_domain=maxval(c_no)
      C_no_inlet_max=maxval(c_no_inlet)
      C_no_max=max(C_no_domain,C_no_inlet_max)
      if (maxval(c_initail).EQ.0) then
      dtmax=((C_no_max)**(-1.0))*tstp_adj
      else
      dtmax=((C_no_max)**(-1.0))*tstp_adj
      endif
      
******************************************************************************************************************************
******************************************************************************************************************************
*     Output files
      print *,'Time step',t
      print *, 'dt=', dtmax
      print *, 'time now is', tt
      print *,'============================='
      print *,' '
*     output files for the debuging
      if (output_file_formate.EQ.2) then
      goto 15004
      else
      goto 1504
      endif
**********************************************************************************************************************      
!=======================================================================
15004 file_no=13
      if (t.EQ.2) then
      open (unit=file_no,file='USG1_H.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) USG1_H(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'USG1_H.txt',MM,NN,nt-1,1,no_of_sec)
      endif
      !=======================================================================
      file_no=14
      if (t.EQ.2) then
      open (unit=file_no,file='USG2_H.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) USG2_H(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'USG2_H.txt',MM,NN,nt-1,1,no_of_sec)
      endif
      !=======================================================================
      file_no=15
      if (t.EQ.2) then
      open (unit=file_no,file='USG1_U.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) USG1_U(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'USG1_U.txt',MM,NN,nt-1,1,no_of_sec)
      endif
      !=======================================================================
      file_no=16
      if (t.EQ.2) then
      open (unit=file_no,file='USG2_U.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) USG2_U(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'USG2_U.txt',MM,NN,nt-1,1,no_of_sec)
      endif
      !=======================================================================
      file_no=17
      if (t.EQ.2) then
      open (unit=file_no,file='USG1_V.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) USG1_V(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'USG1_V.txt',MM,NN,nt-1,1,no_of_sec)
      endif
      !=======================================================================
      file_no=18
      if (t.EQ.2) then
      open (unit=file_no,file='USG2_V.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) USG2_V(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'USG2_V.txt',MM,NN,nt-1,1,no_of_sec)
      endif
      !=======================================================================
      file_no=19
      if (t.EQ.2) then
      open (unit=file_no,file='USG1_C.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) USG1_C(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'USG1_C.txt',MM,NN,nt-1,1,no_of_sec)
      endif
      !=======================================================================
      file_no=20
      if (t.EQ.2) then
      open (unit=file_no,file='USG2_C.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) USG2_C(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'USG2_C.txt',MM,NN,nt-1,1,no_of_sec)
      endif
      !=======================================================================
      file_no=21
      if (t.EQ.2) then
      open (unit=file_no,file='Lim_H.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) Lim_H(:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'Lim_H.txt',MM,NN,nt-1,2,no_of_sec)
      endif
      !=======================================================================
      file_no=22
      if (t.EQ.2) then
      open (unit=file_no,file='Lim_U.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) Lim_U(:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'Lim_U.txt',MM,NN,nt-1,2,no_of_sec)
      endif
      !=======================================================================
      file_no=23
      if (t.EQ.2) then
      open (unit=file_no,file='Lim_V.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) Lim_V(:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'Lim_V.txt',MM,NN,nt-1,2,no_of_sec)
      endif
      !=======================================================================
      file_no=24
      if (t.EQ.2) then
      open (unit=file_no,file='Lim_C.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) Lim_C(:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'Lim_C.txt',MM,NN,nt-1,2,no_of_sec)
      endif
      !=======================================================================
      file_no=25
      if (t.EQ.2) then
      open (unit=file_no,file='Hpre.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) Hpre(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'Hpre.txt',MM,NN,nt-1,1,no_of_sec)
      endif
      !=======================================================================
      file_no=26
      if (t.EQ.2) then
      open (unit=file_no,file='Upre.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) Upre(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'Upre.txt',MM,NN,nt-1,1,no_of_sec)
      endif
      !===================================================================
      file_no=27
      if (t.EQ.2) then
      open (unit=file_no,file='Vpre.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) Vpre(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'Vpre.txt',MM,NN,nt-1,1,no_of_sec)
      endif
      !=======================================================================
      file_no=28
      if (t.EQ.2) then
      open (unit=file_no,file='Cpre.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) Cpre(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'Cpre.txt',MM,NN,nt-1,1,no_of_sec)
      endif
      !=======================================================================
      file_no=30
      if (t.EQ.2) then
      open (unit=file_no,file='Hgi.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) Hgi(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'Hgi.txt',1,NN,nt-1,1,no_of_sec)
      endif
!=======================================================================
      file_no=31
      if (t.EQ.2) then
      open (unit=file_no,file='Ugi.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) Ugi(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'Ugi.txt',1,NN,nt-1,1,no_of_sec)
      endif
!=======================================================================
      file_no=32
      if (t.EQ.2) then
      open (unit=file_no,file='Vgi.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) Vgi(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'Vgi.txt',1,NN,nt-1,1,no_of_sec)
      endif
!=======================================================================
      file_no=33
      if (t.EQ.2) then
      open (unit=file_no,file='Cgi.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) Cgi(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'Cgi.txt',1,NN,nt-1,1,no_of_sec)
      endif
!=======================================================================
      file_no=34
      if (t.EQ.2) then
      open (unit=file_no,file='Hgo.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) Hgo(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'Hgo.txt',1,NN,nt-1,1,no_of_sec)
      endif
!=======================================================================
      file_no=35
      if (t.EQ.2) then
      open (unit=file_no,file='Ugo.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) Ugo(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'Ugo.txt',1,NN,nt-1,1,no_of_sec)
      endif
!=======================================================================
      file_no=36
      if (t.EQ.2) then
      open (unit=file_no,file='Vgo.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) Vgo(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'Vgo.txt',1,NN,nt-1,1,no_of_sec)
      endif
!=======================================================================
      file_no=37
      if (t.EQ.2) then
      open (unit=file_no,file='Cgo.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) Cgo(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'Cgo.txt',1,NN,nt-1,1,no_of_sec)
      endif
!=======================================================================
      file_no=38
      if (t.EQ.2) then
      open (unit=file_no,file='Hgw.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) Hgw(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'Hgw.txt',MM,2,nt-1,1,no_of_sec)
      endif
!=======================================================================
      file_no=39
      if (t.EQ.2) then
      open (unit=file_no,file='Ugw.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) Ugw(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'Ugw.txt',MM,2,nt-1,1,no_of_sec)
      endif
!=======================================================================
      file_no=40
      if (t.EQ.2) then
      open (unit=file_no,file='Vgw.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) Vgw(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'Vgw.txt',MM,2,nt-1,1,no_of_sec)
      endif
!=======================================================================
      file_no=41
      if (t.EQ.2) then
      open (unit=file_no,file='Cgw.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) Cgw(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'Cgw.txt',MM,2,nt-1,1,no_of_sec)
      endif
!=======================================================================
      file_no=42
      if (t.EQ.2) then
      open (unit=file_no,file='hl_eta.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) hl_eta(:,:,1),hl_eta(:,:,3)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'hl_eta.txt',MM,NN,nt-1,2,no_of_sec)
      endif
!=======================================================================
      file_no=43
      if (t.EQ.2) then
      open (unit=file_no,file='hr_eta.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) hr_eta(:,:,1),hr_eta(:,:,3)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'hr_eta.txt',MM,NN,nt-1,2,no_of_sec)
      endif
!=======================================================================
      file_no=44
      if (t.EQ.2) then
      open (unit=file_no,file='ul_eta.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) ul_eta(:,:,1),ul_eta(:,:,3)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'ul_eta.txt',MM,NN,nt-1,2,no_of_sec)
      endif
!=======================================================================
      file_no=45
      if (t.EQ.2) then
      open (unit=file_no,file='ur_eta.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) ur_eta(:,:,1),ur_eta(:,:,3)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'ur_eta.txt',MM,NN,nt-1,2,no_of_sec)
      endif
!=======================================================================
      file_no=46
      if (t.EQ.2) then
      open (unit=file_no,file='vl_eta.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) vl_eta(:,:,1),vl_eta(:,:,3)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'vl_eta.txt',MM,NN,nt-1,2,no_of_sec)
      endif
!=======================================================================
      file_no=47
      if (t.EQ.2) then
      open (unit=file_no,file='vr_eta.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) vr_eta(:,:,1),vr_eta(:,:,3)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'vr_eta.txt',MM,NN,nt-1,2,no_of_sec)
      endif
!=======================================================================
      file_no=48
      if (t.EQ.2) then
      open (unit=file_no,file='cl_eta.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) cl_eta(:,:,1),cl_eta(:,:,3)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'cl_eta.txt',MM,NN,nt-1,2,no_of_sec)
      endif
!=======================================================================
      file_no=49
      if (t.EQ.2) then
      open (unit=file_no,file='cr_eta.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) cr_eta(:,:,1),cr_eta(:,:,3)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'cr_eta.txt',MM,NN,nt-1,2,no_of_sec)
      endif
!=======================================================================
      file_no=50
      if (t.EQ.2) then
      open (unit=file_no,file='hl_zai.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) hl_zai(:,:,2),hl_zai(:,:,4)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'hl_zai.txt',MM,NN,nt-1,2,no_of_sec)
      endif
!=======================================================================
      file_no=51
      if (t.EQ.2) then
      open (unit=file_no,file='hr_zai.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) hr_zai(:,:,2),hr_zai(:,:,4)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'hr_zai.txt',MM,NN,nt-1,2,no_of_sec)
      endif
!=======================================================================
      file_no=52
      if (t.EQ.2) then
      open (unit=file_no,file='ul_zai.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) ul_zai(:,:,2),ul_zai(:,:,4)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'ul_zai.txt',MM,NN,nt-1,2,no_of_sec)
      endif
!=======================================================================
      file_no=53
      if (t.EQ.2) then
      open (unit=file_no,file='ur_zai.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) ur_zai(:,:,2),ur_zai(:,:,4)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'ur_zai.txt',MM,NN,nt-1,2,no_of_sec)
      endif
!=======================================================================
      file_no=54
      if (t.EQ.2) then
      open (unit=file_no,file='vl_zai.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) vl_zai(:,:,2),vl_zai(:,:,4)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'vl_zai.txt',MM,NN,nt-1,2,no_of_sec)
      endif
!=======================================================================
      file_no=55
      if (t.EQ.2) then
      open (unit=file_no,file='vr_zai.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) vr_zai(:,:,2),vr_zai(:,:,4)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'vr_zai.txt',MM,NN,nt-1,2,no_of_sec)
      endif
!=======================================================================
      file_no=56
      if (t.EQ.2) then
      open (unit=file_no,file='cl_zai.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) cl_zai(:,:,2),cl_zai(:,:,4)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'cl_zai.txt',MM,NN,nt-1,2,no_of_sec)
      endif
!=======================================================================
      file_no=57
      if (t.EQ.2) then
      open (unit=file_no,file='cr_zai.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) cr_zai(:,:,2),cr_zai(:,:,4)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'cr_zai.txt',MM,NN,nt-1,2,no_of_sec)
      endif
!=======================================================================
      file_no=58
      if (t.EQ.2) then
      open (unit=file_no,file='H_hat.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) H_hat(:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'H_hat.txt',MM,NN,nt-1,4,no_of_sec)
      endif
!=======================================================================
      file_no=59
      if (t.EQ.2) then
      open (unit=file_no,file='U_hat.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) U_hat(:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'U_hat.txt',MM,NN,nt-1,4,no_of_sec)
      endif
!=======================================================================
      file_no=60
      if (t.EQ.2) then
      open (unit=file_no,file='V_hat.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) V_hat(:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'V_hat.txt',MM,NN,nt-1,4,no_of_sec)
      endif
!=======================================================================
      file_no=61
      if (t.EQ.2) then
      open (unit=file_no,file='C_hat.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) C_hat(:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'C_hat.txt',MM,NN,nt-1,4,no_of_sec)
      endif
!=======================================================================
      file_no=62
      if (t.EQ.2) then
      open (unit=file_no,file='a_hat.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) a_hat(:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'a_hat.txt',MM,NN,nt-1,4,no_of_sec)
      endif
!=======================================================================
      file_no=63
      if (t.EQ.2) then
      open (unit=file_no,file='Eign.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) Eign(:,:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output4D(file_no,'Eign.txt',MM,NN,nt-1,4,4,no_of_sec)
      endif
!=======================================================================
      file_no=64
      if (t.EQ.2) then
      open (unit=file_no,file='R_Eign.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) Eign(:,:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output4D(file_no,'R_Eign.txt',MM,NN,nt-1,4,4,no_of_sec)
      endif
!=======================================================================
      file_no=65
      if (t.EQ.2) then
      open (unit=file_no,file='dh__.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) dh(:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'dh__.txt',MM,NN,nt-1,4,no_of_sec)
      endif
!=======================================================================
      file_no=66
      if (t.EQ.2) then
      open (unit=file_no,file='du__.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) du(:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'du__.txt',MM,NN,nt-1,4,no_of_sec)
      endif
!=======================================================================
      file_no=67
      if (t.EQ.2) then
      open (unit=file_no,file='dv__.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) dv(:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'dv__.txt',MM,NN,nt-1,4,no_of_sec)
      endif
!=======================================================================
      file_no=68
      if (t.EQ.2) then
      open (unit=file_no,file='dc__.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) dc(:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'dc__.txt',MM,NN,nt-1,4,no_of_sec)
      endif
!=======================================================================
      file_no=69
      if (t.EQ.2) then
      open (unit=file_no,file='du_perp.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) du_perp(:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'du_perp.txt',MM,NN,nt-1,4,no_of_sec)
      endif
!=======================================================================
      file_no=70
      if (t.EQ.2) then
      open (unit=file_no,file='du_par.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) du_par(:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'du_par.txt',MM,NN,nt-1,4,no_of_sec)
      endif
!=======================================================================
      file_no=71
      if (t.EQ.2) then
      open (unit=file_no,file='Ws__.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) Ws(:,:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output4D(file_no,'ws__.txt',MM,NN,nt-1,4,4,no_of_sec)
      endif
!=======================================================================
      file_no=72
      if (t.EQ.2) then
      open (unit=file_no,file='df__.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) df(:,:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output4D(file_no,'df__.txt',MM,NN,nt-1,4,4,no_of_sec)
      endif
!=======================================================================
      file_no=73
      if (t.EQ.2) then
      open (unit=file_no,file='F1_L.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) F1_L(:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'F1_L.txt',MM,NN,nt-1,4,no_of_sec)
      endif
!=======================================================================
      file_no=74
      if (t.EQ.2) then
      open (unit=file_no,file='F2_L.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) F2_L(:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'F2_L.txt',MM,NN,nt-1,4,no_of_sec)
      endif
!=======================================================================
      file_no=75
      if (t.EQ.2) then
      open (unit=file_no,file='F3_L.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) F3_L(:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'F3_L.txt',MM,NN,nt-1,4,no_of_sec)
      endif
!=======================================================================
      file_no=76
      if (t.EQ.2) then
      open (unit=file_no,file='F4_L.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) F4_L(:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'F4_L.txt',MM,NN,nt-1,4,no_of_sec)
      endif
!=======================================================================
      file_no=77
      if (t.EQ.2) then
      open (unit=file_no,file='F1_R.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) F1_R(:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'F1_R.txt',MM,NN,nt-1,4,no_of_sec)
      endif
!=======================================================================
      file_no=78
      if (t.EQ.2) then
      open (unit=file_no,file='F2_R.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) F2_R(:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'F2_R.txt',MM,NN,nt-1,4,no_of_sec)
      endif
!=======================================================================
      file_no=79
      if (t.EQ.2) then
      open (unit=file_no,file='F3_R.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) F3_R(:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'F3_R.txt',MM,NN,nt-1,4,no_of_sec)
      endif
!=======================================================================
      file_no=80
      if (t.EQ.2) then
      open (unit=file_no,file='F4_R.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) F4_R(:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'F4_R.txt',MM,NN,nt-1,4,no_of_sec)
      endif
!=======================================================================
      file_no=81
      if (t.EQ.2) then
      open (unit=file_no,file='Flux1.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) Flux1(:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'Flux1.txt',MM,NN,nt-1,4,no_of_sec)
      endif
!=======================================================================
      file_no=82
      if (t.EQ.2) then
      open (unit=file_no,file='Flux2.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) Flux2(:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'Flux2.txt',MM,NN,nt-1,4,no_of_sec)
      endif
!=======================================================================
      file_no=83
      if (t.EQ.2) then
      open (unit=file_no,file='Flux3.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) Flux3(:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'Flux3.txt',MM,NN,nt-1,4,no_of_sec)
      endif
!=======================================================================
      file_no=84
      if (t.EQ.2) then
      open (unit=file_no,file='Flux4.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) Flux4(:,:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'Flux4.txt',MM,NN,nt-1,4,no_of_sec)
      endif
!=======================================================================
      file_no=85
      if (t.EQ.2) then
      open (unit=file_no,file='cb_new.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) cb_new(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'cb_new.txt',MM,NN,nt-1,1,no_of_sec)
      endif
!=======================================================================
      file_no=86
      if (t.EQ.2) then
      open (unit=file_no,file='cd_new.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) cd_new(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'cd_new.txt',MM,NN,nt-1,1,no_of_sec)
      endif
!=======================================================================
      file_no=87
      if (t.EQ.2) then
      open (unit=file_no,file='Ustar_new.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) Ustar_new(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'Ustar_new.txt',MM,NN,nt-1,1,no_of_sec)
      endif
!=======================================================================
      file_no=88
      if (t.EQ.2) then
      open (unit=file_no,file='Vstar_new.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) Vstar_new(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'Vstar_new.txt',MM,NN,nt-1,1,no_of_sec)
      endif
!=======================================================================
      file_no=89
      if (t.EQ.2) then
      open (unit=file_no,file='E_new.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) E_new(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'E_new.txt',MM,NN,nt-1,1,no_of_sec)
      endif
!=======================================================================
      file_no=90
      if (t.EQ.2) then
      open (unit=file_no,file='dz_x.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) dz_x(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'dz_x.txt',MM,NN,nt-1,1,no_of_sec)
      endif
!=======================================================================
      file_no=91
      if (t.EQ.2) then
      open (unit=file_no,file='dz_y.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) dz_y(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'dz_y.txt',MM,NN,nt-1,1,no_of_sec)
      endif
!=======================================================================
      file_no=92
      if (t.EQ.2) then
      open (unit=file_no,file='q1_new.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) q1_new(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'q1_new.txt',MM,NN,nt-1,1,no_of_sec)
      endif
!=======================================================================
      file_no=93
      if (t.EQ.2) then
      open (unit=file_no,file='q2_new.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) q2_new(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'q2_new.txt',MM,NN,nt-1,1,no_of_sec)
      endif
!=======================================================================
      file_no=94
      if (t.EQ.2) then
      open (unit=file_no,file='q3_new.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) q3_new(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'q3_new.txt',MM,NN,nt-1,1,no_of_sec)
      endif
!=======================================================================
      file_no=95
      if (t.EQ.2) then
      open (unit=file_no,file='q4_new.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) q4_new(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'q4_new.txt',MM,NN,nt-1,1,no_of_sec)
      endif
!=======================================================================
      file_no=96
      if (t.EQ.2) then
      open (unit=file_no,file='W1__.txt',action="write",
     &      status="unknown")
      write (file_no,*) W1_old(:,:),W1(:,:)
      else
      write (file_no,*) W1(:,:)
      endif
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'W1__.txt',MM,NN,nt,1,no_of_sec)
      endif
!=======================================================================
      file_no=97
      if (t.EQ.2) then
      open (unit=file_no,file='W2__.txt',action="write",
     &      status="unknown")
      write (file_no,*) W2_old(:,:),W2(:,:)
      else
      write (file_no,*) W2(:,:)
      endif
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'W2__.txt',MM,NN,nt,1,no_of_sec)
      endif
!=======================================================================
      file_no=98
      if (t.EQ.2) then
      open (unit=file_no,file='W3__.txt',action="write",
     &      status="unknown")
      write (file_no,*) W3_old(:,:),W3(:,:)
      else
      write (file_no,*) W3(:,:)
      endif
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'W3__.txt',MM,NN,nt,1,no_of_sec)
      endif
!=======================================================================
      file_no=99
      if (t.EQ.2) then
      open (unit=file_no,file='W4__.txt',action="write",
     &      status="unknown")
      write (file_no,*) W4_old(:,:),W4(:,:)
      else
      write (file_no,*) W4(:,:)
      endif
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'W4__.txt',MM,NN,nt,1,no_of_sec)
      endif
!=======================================================================
      file_no=100
      if (t.EQ.2) then
      open (unit=file_no,file='H__.txt',action="write",
     &      status="unknown")
      write (file_no,*) H_old(:,:),H(:,:)
      else
      write (file_no,*) H(:,:)
      endif
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'H__.txt',MM,NN,nt,1,no_of_sec)
      endif
      !=======================================================================
      file_no=101
      if (t.EQ.2) then
      open (unit=file_no,file='U__.txt',action="write",
     &      status="unknown")
      write (file_no,*) U_old(:,:),U(:,:)
      else
      write (file_no,*) U(:,:)
      endif
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'U__.txt',MM,NN,nt,1,no_of_sec)
      endif
      !=======================================================================
      file_no=102
      if (t.EQ.2) then
      open (unit=file_no,file='V__.txt',action="write",
     &      status="unknown")
      write (file_no,*) V_old(:,:),V(:,:)
      else
      write (file_no,*) V(:,:)
      endif
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'V__.txt',MM,NN,nt,1,no_of_sec)
      endif
      !=======================================================================
      file_no=103
      if (t.EQ.2) then
      open (unit=file_no,file='C__.txt',action="write",
     &      status="unknown")
      write (file_no,*) C_old(:,:),C(:,:)
      else
      write (file_no,*) C(:,:)
      endif
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'C__.txt',MM,NN,nt,1,no_of_sec)
      endif
!=======================================================================
      file_no=104
      if (t.EQ.2) then
      open (unit=file_no,file='zp__.txt',action="write",
     &      status="unknown")
      write (file_no,*) zp_old(:,:),zp(:,:)
      else
      write (file_no,*) zp(:,:)
      endif
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'zp__.txt',MM,NN,nt,1,no_of_sec)
      endif
!=======================================================================
      file_no=106
      if (t.EQ.2) then
      open (unit=file_no,file='dz_x_old.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) dz_x_old(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'dz_x_old.txt',MM,NN,nt-1,1,no_of_sec)
      endif
!=======================================================================
      file_no=107
      if (t.EQ.2) then
      open (unit=file_no,file='dz_y_old.txt',action="write",
     &      status="unknown")
      endif
      write (file_no,*) dz_y_old(:,:)
      if (t.EQ.nt) then
      close (file_no)
      call output3D(file_no,'dz_y_old.txt',MM,NN,nt-1,1,no_of_sec)
      endif
!=======================================================================
******************************************************************************************************
!=======================================option 2 in the output file "prefered for final results"
1504  if (npt.GE.nprint) then
      npt=1
      else
      npt=npt+1
      endif
      if ((t.EQ.2).OR.(npt.EQ.nprint)) then
      file_no=258
      if (t.EQ.2) then
      open (unit=file_no,file='time.txt',action="write",
     &      status="unknown")
      Write (file_no,*) 0
      endif
      Write (file_no,*) tt
!      endif
      if (t.EQ.nt) then
      close (file_no)
      endif
*********************************************************************************************      
   !   if ((t.EQ.2).OR.(npt.EQ.nprint)) then
      file_no=108
      if (t.EQ.2) then
      npt=2
      endif
      if (t.EQ.2) then
      open (unit=file_no,file="Discharge per width.txt",action="write",
     &      status="unknown")
      write (file_no,*) M,no_of_sec+1
      write (file_no,*) nt/nprint+2
      endif
      if(t.EQ.2) then
      do 475 k=1,no_of_sec
      write (file_no,*) 1
      write (file_no,*) H_old(:,k)*U_old(:,k)
 475  continue
      do 476 k=1,no_of_sec
      write (file_no,*) t
      write (file_no,*) H(:,k)*U(:,k)  
476   continue
      else
      do 477 k=1,no_of_sec
      write (file_no,*) t
      write (file_no,*) H(:,k)*U(:,k) 
477   continue
      endif
      if (t.EQ.nt) then           
      close (file_no)
      endif
      !=================================================== General outputfile
      file_no=1000
      if (t.EQ.2) then
      open (unit=file_no,file="Hc.txt",action="write",
     &      status="unknown")
      write (file_no,*) M,no_of_sec+1
      write (file_no,*) nt/nprint+2
      endif
      if(t.EQ.2) then
      do 410 k=1,no_of_sec+1
      write (file_no,*) Hc_intial(:,k)
 410  continue
      do 415 k=1,no_of_sec+1
      write (file_no,*) Hc(:,k)  
415   continue
      else
      do 416 k=1,no_of_sec+1
      write (file_no,*) Hc(:,k)  
416   continue
      endif
      if (t.EQ.nt) then           
      close (file_no)
      endif
!=================================================== General outputfile
      file_no=1001
      if (t.EQ.2) then
      open (unit=file_no,file="Uc.txt",action="write",
     &      status="unknown")
      write (file_no,*) M,no_of_sec+1
      write (file_no,*) nt/nprint+2
      endif
      if(t.EQ.2) then
      do 430 k=1,no_of_sec+1
      write (file_no,*) Uc_intial(:,k)
 430  continue
      do 431 k=1,no_of_sec+1
      write (file_no,*) Uc(:,k)  
431   continue
      else
      do 432 k=1,no_of_sec+1
      write (file_no,*) Uc(:,k)  
432   continue
      endif
      if (t.EQ.nt) then           
      close (file_no)
      endif
  !=================================================== General outputfile
      file_no=1002
      if (t.EQ.2) then
      open (unit=file_no,file="Vc.txt",action="write",
     &      status="unknown")
      write (file_no,*) M,no_of_sec+1
      write (file_no,*) nt/nprint+2
      endif
      if(t.EQ.2) then
      do 433 k=1,no_of_sec+1
      write (file_no,*) Vc_intial(:,k)
 433  continue
      do 434 k=1,no_of_sec+1
      write (file_no,*) Vc(:,k)  
434   continue
      else
      do 435 k=1,no_of_sec+1
      write (file_no,*) Vc(:,k)  
435   continue
      endif
      if (t.EQ.nt) then           
      close (file_no)
      endif
      !=================================================== General outputfile
      file_no=1003
      if (t.EQ.2) then
      open (unit=file_no,file="Cc.txt",action="write",
     &      status="unknown")
      write (file_no,*) M,no_of_sec+1
      write (file_no,*) nt/nprint+2
      endif
      if(t.EQ.2) then
      do 436 k=1,no_of_sec+1
      write (file_no,*) Cc_intial(:,k)
 436  continue
      do 437 k=1,no_of_sec+1
      write (file_no,*) Cc(:,k)  
437   continue
      else
      do 438 k=1,no_of_sec+1
      write (file_no,*) Cc(:,k)  
438   continue
      endif
      if (t.EQ.nt) then           
      close (file_no)
      endif
!=================================================== General outputfile
      file_no=1004
      do 442 j=1,MM
      do 443 k=1,NN
      Fr(j,k)=U(j,k)/SQRT(H(j,k)*9.81)
 443  continue     
 442  continue
      if(t.EQ.2) then
      open (unit=file_no,file="Fr.txt",action="write",
     &      status="unknown")
      write (file_no,*) MM,no_of_sec
      write (file_no,*) nt/nprint+2-1
      endif
      do 441 k=1,no_of_sec
      write (file_no,*) t
      write (file_no,*) Fr(:,k)
441   continue
      if (t.EQ.nt) then           
      close (file_no)
      endif
!      file_no=1004
!      if (t.EQ.2) then
!      open (unit=file_no,file="zp.txt",action="write",
!     &      status="unknown")
!      write (file_no,*) MM,no_of_sec
!      write (file_no,*) nt/nprint+2
!      endif
!      if(t.EQ.2) then
!      do 439 k=1,no_of_sec
!      write (file_no,*) zp_old(:,k)
! 439  continue
!      do 440 k=1,no_of_sec+1
!      write (file_no,*) zp(:,k)  
!440   continue
!      else
!      do 441 k=1,no_of_sec+1
!      write (file_no,*) zp(:,k)  
!441   continue
!      endif
!      if (t.EQ.nt) then           
!      close (file_no)
!      endif
!=================================================== General outputfile
      if (t.EQ.2) then
      open (unit=500,file="Xo.txt",action="write",
     &      status="unknown")
      write (500,*) M,no_of_sec+1
      do 412 k=1,no_of_sec+1
      write (500,*) Xo(:,k)  
412   continue         
      close (500)
      endif
!=================================================== General outputfile
      if (t.EQ.2) then
      open (unit=501,file="Yo.txt",action="write",
     &      status="unknown")
      write (501,*) M,no_of_sec+1
      do 414 k=1,no_of_sec+1
      write (501,*) Yo(:,k)  
414   continue         
      close (501)
      endif
!====================================================================      
      if (t.EQ.2) then
      open (unit=502,file="Zo.txt",action="write",
     &      status="unknown")
      write (502,*) M,no_of_sec+1
      write (502,*) nt/nprint+2
      do 413 k=1,no_of_sec+1
      write (502,*) Zo_old(:,k)  
413   continue
      do 425 k=1,no_of_sec+1
      write (502,*) Zo(:,k)  
425   continue
      else
      do 426 k=1,no_of_sec+1
      write (502,*) Zo(:,k)  
426   continue  
      endif
      if (t.EQ.nt) then        
      close (502)
      endif
!=======================================================================
*******************************************************************************************************************************
      endif
  100 continue
      pause
      end
********************************************************************************************************************************
!==================================================================================================================================
*     Output files "Formating the output files for 3D vectors"
      subroutine output3D(fileno,filename,MM,NN,nt,no_of_faces,no_secs)
      Real ZZ
      Integer fileno,MM,NN,nt,no_of_faces,size,t,NNN,n
      Character*(*) filename
      Allocatable ZZ(:),z(:,:,:,:)
      print *, 'Formatting=================>  ',filename
      size=MM*NN*nt*no_of_faces
      Allocate (ZZ(size),z(MM,NN,no_of_faces,nt), STAT=Keep)
      open (unit=fileno,file=filename)
      read (fileno,*) ZZ
      close (fileno)
      i=1
      do 4 t=1,nt
      do 7 n=1,no_of_faces
      do 5 k=1,NN
      do 6 j=1,MM
      z(j,k,n,t)=ZZ(i)
      i=i+1
  6   continue
  5   continue
  7   continue
  4   continue
      if (no_secs.GT.NN) then
      NNN=NN
      else
      NNN=no_secs
      endif   
      open (unit=fileno,file=filename,action="write",
     &      status="unknown")

      do 1 t=1,nt
      write (fileno,*) 'Time step no=',t
      do 8 n=1,no_of_faces
      write (fileno,*) 'kk face number or dir. 1 for zai&2 for eta',n
      do 2 k=1,NNN
      write (fileno,*) 'k=',k
      write (fileno,*) z(:,k,n,t)   
 2    continue
 8    continue
 1    continue
      close (fileno)
      return
      end
************************************************************************************************************************ 
*********************************************************************************************************************************
*     Output files "Formating the output files for 4D vectors"
      subroutine output4D(fileno,filename,MM,NN,nt,no_of_faces,No_coeff,
     &                    no_secs)
      Real ZZ
      Integer fileno,MM,NN,nt,no_of_faces,size,t,No_coeff,n,nnn,NNNN
      Character*(*) filename
      Allocatable ZZ(:),z(:,:,:,:,:)
      size=MM*NN*nt*no_of_faces*No_coeff
      Allocate (ZZ(size),z(MM,NN,No_coeff,no_of_faces,nt), STAT=Keep)
      print *, 'Formatting=================>  ',filename
      open (unit=fileno,file=filename)
      read (fileno,*) ZZ
      close (fileno)
      i=1
      do 4 t=1,nt
      do 7 nnn=1,No_coeff
      do 9 n=1,no_of_faces
      do 5 k=1,NN
      do 6 j=1,MM
      z(j,k,n,nnn,t)=ZZ(i)
      i=i+1
  6   continue
  5   continue
  9   continue
  7   continue
  4   continue 
      if (no_secs.GT.NN) then
      NNNN=NN
      else
      NNNN=no_secs
      endif   
      open (unit=fileno,file=filename,action="write",
     &      status="unknown")

      do 1 t=1,nt
      write (fileno,*) 'Time step no=',t
      do 8 nnn=1,no_of_faces
      write (fileno,*) 'coeff.no',nnn
      do 10 n=1,No_coeff
      write (fileno,*) 'face no.',n
      do 2 k=1,NNNN
      write (fileno,*) 'k=',k
      write (fileno,*) z(:,k,n,nnn,t)   
 2    continue
 10   continue
 8    continue
 1    continue
      close (fileno)
      return
      end
************************************************************************************************************************ 
*********************************************************************************************************************************
!   FLUX LIMITER FUNCTION             
      Function  homar(order,beta,USG,DSG)
      Real beta,USG,DSG,qq1,qq2,ab,order
      homar=0.00
      if (order.EQ.2) then
      ab=USG*DSG
      qq1=max(abs(USG),abs(DSG))
      qq2=beta*min(abs(USG),abs(DSG))
      if (ab.GT.0.0.AND.USG.LT.0.0) then
      homar=-1*min(qq1,qq2)
      elseif (ab.GT.0.0.AND.USG.GT.0.0) then
      homar=min(qq1,qq2)
      endif
      endif
      return
      end      
*********************************************************************************************************************************      
!     LEFT AND RIGHT FLUX CALCULATOR 
      Function f1(h,u)
      real  h,u
      f1=h*u
      return
      end
*==============================================================================================================================
      Function f2(h,u,uper,cosin)
      real  h,u,uper,cosin
      f2=h*u*uper+0.5*9.81*cosin*h**2
      return
      end
*============================================================================================================================  
      Function f3(h,v,uper,sin)
      real  h,v,uper,sin
      f3=h*v*uper+0.5*9.81*sin*h**2
      return
      end
*============================================================================================================================    
      Function f4(h,c,uper)
      real  h,u,uper
      f4=h*c*uper
      return
      end
*=========================================================================================================================== 
      Function  v_fall(R,D,neu)                                         !Fall velocity "Dietrich 1982"
      real R,D,neu,Rep,Rf,b1,b2,b3,b4
      Rep=sqrt(R*9.81*D)*D/neu
      b1=2.891394
      b2=0.95296
      b3=0.056835
      b4=0.002892
      b5=0.000245
      Rf=exp(-b1+b2*log(Rep)-b3*(log(Rep))**2-b4*(log(Rep))**3+
     &      b5*(log(Rep))**4)
      v_fall=Rf*sqrt(R*9.81*D)
      return
      end
*=============================================================================================================================
!     Near bed concentration "Parker"
      Function c_bed1(C)
      real C,cb
      c_bed1=2.0*C
      return
      end
*============================================================================================================================
      Function c_bed2(D,R,neu,u,d90)        !needs more work
      real C,cb,D,R,neu,u,d90
      c_bed2=2.0*C
      return
      end
***************************************************************************************************************************** 
*     Drag coeeficient                  "Haaland 1983"
      Function drag(u,h,neu,ks)
      real   u,h,neu,Re,CD,ks
      Re=u*h/neu
      if (Re.GT.1000) then
      if (ks/h.LT.0.2) then
      drag=0.204/(log(1.725/Re+(ks/(14.8*h))**1.11))**2.0
      else
      drag=1.56*.01*(ks/h)**(0.33333333333333)
      endif
      else
      drag=0.204/(log(1.725/1000+(ks/(14.8*h))**1.11))**2.0
      endif
      return
      end
*******************************************************************************************************************************
*     Nondimensional etrainmentrate "Wright &Parker 2004" 
      Function Es(Ustar,Vstar,vf,s,Rep)                                         !which slope should used
      real Ustar,Vstar,vf,s,E,B,Zu,Vmagstar
      B=7.8e-7
      Vmagstar=(Ustar**2+Vstar**2)**0.5
      Zu=(Vmagstar*Rep**0.6)*(s**0.08)/vf
      Es=(B*Zu**5)/(1+B*(Zu**5)/0.3)
      return
      end
**********************************************************************************************************************************
*     Nondimensional etrainmentrate "Van Rijn 1984"
      Function E_van(Ustar,Vstar,D,neu,h,R)
      real Ustar,Vstar,vf,D,neu,E,ks,h,tao_stac,tao_sta,R,Vmagstar,Rep
      Vmagstar=(Ustar**2+Vstar**2)**0.5
      tao_sta=Vmagstar/(R*9.81*D)
      Rep=sqrt(R*9.81*D)*D/neu
      tao_stac=0.22*Rep**(-0.6)+0.06*10**(-7.7*Rep**(-0.6))
      if (3*D.GT.0.01*h) then
      b=3*D
      else
      b=0.01*h
      endif
      if (tao_sta.LE.tao_stac) then
      E_van=0.0
      else
      E_van=0.015*(D/b)*(tao_sta/tao_stac-1)**(1.5)*Rep**(-0.2)
      endif
      return
      end
********************************************************************************************************************************
*     Courant number function
      Function courant(u,v,h,dx,dy)
      real u,v,h,dx,dy,c,R
      cournat=((abs(u)+sqrt(9.81*h))/dx+
     &          (abs(v)+sqrt(9.81*h)))/dy
      return
      end
********************************************************************************************************************************
!     Time step Function
      Function time_step(ui,vi,hi,u,v,h,dx,dy,MM,NN)
      Allocatable dt_domain(:,:)
      real ui,vi,hi,u,v,h,dx,dy,dt_in,dt_domain
      integer MM,NN
      Allocate (dt_domain(MM,NN),STAT=Keep)
      dt_in=(((abs(ui)+sqrt(9.81*hi))/dx+
     &          (abs(vi)+sqrt(9.81*hi)))/dy)**(-1)
      do 1011 j=1,MM
      do 1012 k=1,NN
      dt_domain(j,k)=(((abs(u)+sqrt(9.81*h))/dx+
     &          (abs(v)+sqrt(9.81*h)))/dy)**(-1)
1012  continue
1011  continue
      dt=MINVAL(dt_domain(:,:))
      time_step=min(dt_in,dt)
      return
      end
*******************************************************************************************************************************
*     Boundary condition Function for water depth
      Function Bound_H(hi,hdom1,hdom2,h_c,ho,bc)
      Real hi,hdom1,hdom2,h_c,ho,bc
      if (bc.EQ.1) then         !open boundary condition inflow supercritical
      Bound_h=hi
      elseif (bc.EQ.2) then         !open boundary condition inflow subcritical
      Bound_h=hdom1
      elseif (bc.EQ.3) then         ! open BC supercritical outflow
      Bound_h=hdom2
      elseif (bc.EQ.4.And.hdom.LE.h_c) then         ! Open BC subcritical outflow
      Bound_h=hdom2
      elseif (bc.EQ.4.And.h_in.GT.h_c) then
      Bound_h=ho
      elseif (bc.EQ.5) then         !soild wall boundary condition "free slip" HL Wall
      Bound_h=hdom2
      elseif (bc.EQ.6) then          !soild wall boundary condition "free slip" VL Wall
      Bound_h=hdom1
      endif
      return
      end
*******************************************************************************************************************************
*     Boundary condition Function for concentration
      Function Bound_C(ci,cdom1,cdom2,bc)
      Real ci,cdom1,cdom2,bc
      if (bc.EQ.1) then         !open boundary condition inflow supercritical
      Bound_c=ci
      elseif (bc.EQ.2) then         !open boundary condition inflow subcritical
      Bound_c=ci
      elseif (bc.EQ.3) then         ! open BC supercritical outflow
      Bound_c=cdom2
      elseif (bc.EQ.4) then
      Bound_c=cdom2
      elseif (bc.EQ.5) then         !soild wall boundary condition "free slip" HL Wall
      Bound_c=cdom2
      elseif (bc.EQ.6) then          !soild wall boundary condition "free slip" VL Wall
      Bound_c=cdom1
      endif
      return
      end
*******************************************************************************************************************************
*     Boundary condition Function for velocity"U"
      Function Bound_U(ui,udom,vdom,cphi,sphi,bc)
      Real ui,udom,vdom,cphi,sphi,bc
      if (bc.EQ.1) then         !open boundary condition inflow supercritical
      Bound_U=ui
      elseif (bc.EQ.2) then         !open boundary condition inflow subcritical
      Bound_U=ui
      elseif (bc.EQ.3) then         ! open BC supercritical outflow
      Bound_U=udom
      elseif (bc.EQ.4) then
      Bound_U=udom
      elseif (bc.EQ.5) then         !soild wall boundary condition "free slip" HL Wall
      Bound_U=udom*(sphi**2-cphi**2)-2*vdom*sphi*cphi
      elseif (bc.EQ.6) then          !soild wall boundary condition "free slip" VL Wall
      Bound_U=-2.0*udom*sphi*cphi+vdom*(cphi**2-sphi**2)
      endif
      return
      end
******************************************************************************************************************************
*     Boundary condition Function for velocity"V"
      Function Bound_V(vi,udom,vdom,cphi,sphi,bc)
      Real vi,udom,vdom,cphi,sphi,bc
      if (bc.EQ.1) then         !open boundary condition inflow supercritical
      Bound_V=vi
      elseif (bc.EQ.2) then         !open boundary condition inflow subcritical
      Bound_V=vi
      elseif (bc.EQ.3) then         ! open BC supercritical outflow
      Bound_V=vdom
      elseif (bc.EQ.4) then
      Bound_V=vdom
      elseif (bc.EQ.5) then         !soild wall boundary condition "free slip" HL Wall
      Bound_V=-2.0*udom*sphi*cphi+vdom*(cphi**2-sphi**2)
      elseif (bc.EQ.6) then          !soild wall boundary condition "free slip" VL Wall
      Bound_V=udom*(sphi**2-cphi**2)-2*vdom*sphi*cphi
      endif
      return
      end
******************************************************************************************************************************