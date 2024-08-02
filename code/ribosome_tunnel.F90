!******************************************************************************
!
!  Module: ribosome_tunnel_mod
!
!  Description: This module will build a rigid model of ribosome exit
!  tunnel, in which a protein can cotranslational folding.
! 
!  Authors: Peng Tao (current)
!           Ercheng Wang (revised)
!           Changjun Chen (original)
!
!  Organization: Institute of Biophysics, School of Physics, Huazhong
!  University of Science and Technology, Wuhan 430074, China
!  
!  Note: 1. The initial configuration of protein should be extended and
!           be put along the x axis. Besides, the N terminal of protein 
!           should locate in the positive direction of the x axis.
!        2. Since this module is not a part of pmemd, before compliering
!           this module you should add the module information before module
!           irunmd_mod'into 'Makefile' located in this directory.
!        3. This tunnel model only applied to protein, it does not 
!           interact with other molecules like ions or waters!
!  
!  Tunnel parameters:
!                                                    |     |         
!                                                    cone- #
!                                                   height #
!          tunnelminpos                tunnelmidpos  |     #------
!                     |                            \ |   ## cone-
!                     |                             \|### maxradius
!    -----------------################################------------ 
!    cylinder- -------#
!    radius           ################################   tunnelmaxpos
!                     |                               ###   /
!                     |                                  ##/
!                     |                                    #
!                     |                                    #
!                     |                                    #
!                     |-----------tunnellength-------------|
!   
!******************************************************************************

module ribosome_tunnel_mod

use mdin_ctrl_dat_mod

  implicit none
  integer,private              :: nres,natom,nresmol,natommol
  integer,private              :: relax_res_cnt
  integer,private              :: rpanel_rejected,lpanel_rejected,&
                                  cylinder_rejected,cone_rejected
  real*8,private               :: tunnelminpos,tunnelmidpos,tunnelmaxpos
  !real*8,private               :: ctdz
  !real*8,private               :: xminpos
  
  public  :: parse_tunnel_params ,build_tunnel ,apply_tunnel  
  !private :: centroid_of_protein ,normalize

contains

subroutine parse_tunnel_params(nres,natom,res_atms)

  implicit none
  integer,intent(in)           :: nres,natom
  integer,intent(in)           :: res_atms(nres)
  integer                      :: i
  !real*8,intent(in)            :: crd(3,natom)
  
  nresmol=0
  do i=1,size(res_atms)
    if (res_atms(i+1)-res_atms(i) > 3) then
      nresmol=nresmol+1
    else
      exit
    end if
  end do
  
  natommol=res_atms(nresmol+1)-1

  !xminpos=minval(crd(1,1:natommol))
   
  !call centroid_of_protein(crd(1:3,1:natommol))

  write(*,*) ""
  write(*,*) "**************************************************************"
  write(*,*) "*                 Cotranslational folding                    *"
  write(*,*) "**************************************************************"
  write(*,*) ""
  write(*,*) "Simulated system information:"
  write(*,'(a,i10)') "  residues(system)    :",nres
  write(*,'(a,i10)') "  residues(protein)   :",nresmol
  write(*,'(a,i10)') "  atoms(system)       :",natom
  write(*,'(a,i10)') "  atoms(protein)      :",natommol
  write(*,*) ""
  write(*,*) "Tunnel model parameters:"
  write(*,'(a,f10.2)') "  tunnel length       :",tunnellength
  write(*,'(a,f10.2)') "  cylinder radius     :",cylinderradius
  ! If all residues have been released, the cylinder minimum position will set 
  ! to be xminpos.
  write(*,'(a,f10.2)') "  cylinder lowpos     :",xminpos
  write(*,'(a,f10.2)') "  cone max radius     :",conemaxradius
  write(*,'(a,f10.2)') "  cone height         :",coneheight
  write(*,'(a,i10)') "  release speed       :",staysteps ! per residue
  write(*,*) ""
  write(*,*) "Output parameters:"
  write(*,'(a,i10)') "  print interval      :",ntpr
  write(*,'(a,2(f10.2))') "  protein centroid    :",ctdy,ctdz
  write(*,*) ""

end subroutine

subroutine build_tunnel(crd,res_atms)

  implicit none
  real*8,intent(in)         :: crd(3,natom)
  integer,intent(in)        :: res_atms(nres)
  real*8                    :: resmaxpos(nresmol),resminpos(nresmol)
  integer                   :: i,j
  
  cylinder_rejected=0;rpanel_rejected=0;lpanel_rejected=0;cone_rejected=0

  resmaxpos(1:nresmol)=-1000.d0;resminpos(1:nresmol)=1000.d0
  do i=1,nresmol
    do j=res_atms(i),res_atms(i+1)-1
       if (crd(1,j) < resminpos(i)) resminpos(i)=crd(1,j)
       if (crd(1,j) > resmaxpos(i)) resmaxpos(i)=crd(1,j)     
    end do
  end do
  
  write(*,*) "Debug information:"
  write(*,'(5x,a,2(f8.2),a,2(f8.2))') "resminpos:",resminpos(1),resminpos(2),&
           "   ......",resminpos(nresmol-1),resminpos(nresmol)
  write(*,'(5x,a,2(f8.2),a,2(f8.2))') "resmaxpos:",resmaxpos(1),resmaxpos(2),&
           "   ......",resmaxpos(nresmol-1),resmaxpos(nresmol)
  write(*,*) ""
  write(*,'(5x,a,f12.2,a)') "Current time is:",t,"ps."
  write(*,*) ""
  relax_res_cnt=int( int((t+1)/1000) / (dt*staysteps/1000) ) + 1
  if (relax_res_cnt <= nresmol) then
    write(*,'(5x,i5,a)') relax_res_cnt,"th residue is relaxing into the tunnel..."
    tunnelminpos=resminpos(relax_res_cnt)-1
  else
    write(*,'(5x,a)') "  All residues have relaxed into the tunnel."
    write(*,'(5x,a)') "  In this situation, the location of the tunnel  will"
    write(*,'(5x,a)') "  not change over time."
    tunnelminpos=xminpos-1
    relax_res_cnt=nresmol
  end if
  write(*,*) ""
  
  tunnelmaxpos=tunnelminpos+tunnellength
  tunnelmidpos=tunnelmaxpos-coneheight
  
  write(*,*) ""
  write(*,'(10x,a,3(f8.2))') "tunnel position:",tunnelminpos,&
                 tunnelmidpos,tunnelmaxpos
  write(*,*) "**************************************************************"

end subroutine

subroutine apply_tunnel(nstep,crd,vel,res_atms)

  implicit none
  integer,intent(in)        :: nstep
  integer,intent(in)        :: res_atms(nres)
  real*8, intent(in)        :: crd(3,natom)
  real*8, intent(inout)     :: vel(3,natom)
  integer                   :: i,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13
  integer                   :: counter1,counter2,counter3,counter4
  real*8                    :: temp,temp2,temp3
  real*8                    :: tpcrd(3)

  counter1=0;counter2=0;counter3=0;counter4=0
  c1=0;c2=0;c3=0;c4=0;c5=0;c6=0;c7=0
  c8=0;c9=0;c10=0;c11=0;c12=0;c13=0
  
  ! The following loop will confine the protein inside the tunnel model.  
  do i=1,res_atms(relax_res_cnt+1)-1
    temp = sqrt((crd(2,i)-ctdy)**2 + (crd(3,i)-ctdz)**2)
    !left
    if ( crd(1,i) < tunnelminpos ) then
      if ( vel(1,i) > 0 ) cycle
      counter1=counter1+1
      vel(1,i)=-vel(1,i)
    !cylinder
    else if ( crd(1,i) >= tunnelminpos .and. crd(1,i) < tunnelmidpos ) then
      if ( temp > cylinderradius ) then
        tpcrd(1) = 0.0d0                                                             
        tpcrd(2)=crd(2,i)-ctdy                                                    
        tpcrd(3)=crd(3,i)-ctdz                                                    
        call normalize(tpcrd)                                                    
        if (dot_product(vel(1:3,i),tpcrd) < 0) cycle
        vel(1:3,i) = vel(1:3,i) - 2 * (dot_product(vel(1:3,i),tpcrd) * tpcrd(1:3))                                             
        counter2 = counter2 + 1                                                         
      end if                                                                             
    else if ( crd(1,i)>=tunnelmidpos .and. crd(1,i)<tunnelmaxpos ) then
      !cone
      temp2=( crd(1,i)-tunnelmidpos)*conemaxradius/coneheight
      if ( temp > cylinderradius+temp2 .and. temp <= cylinderradius &
          +conemaxradius) then
        tpcrd(1) = -temp*conemaxradius/coneheight
        tpcrd(2) = crd(2,i)-ctdy
        tpcrd(3) = crd(3,i)-ctdz
        call normalize(tpcrd)
        if (dot_product(vel(1:3,i),tpcrd) < 0) cycle
        vel(1:3,i) = vel(1:3,i) - 2 * (dot_product(vel(1:3,i),tpcrd) * tpcrd(1:3))
        counter3 = counter3 + 1
      !right panel
      else if ( temp > cylinderradius+conemaxradius ) then
        if ( vel(1,i)>0 ) cycle
        vel(1,i)=-vel(1,i)
        counter4=counter4+1
      end if
    end if
  end do
  cylinder_rejected=cylinder_rejected+counter2
  rpanel_rejected=rpanel_rejected+counter4
  lpanel_rejected=lpanel_rejected+counter1
  cone_rejected=cone_rejected+counter3

  ! Every ntw step the distribution of atoms of molecule will be printed on the
  ! screen.
  if (mod(nstep,ntpr)==0) then
    write(*,'(a,3(f10.2))') "tunnel position: ",tunnelminpos,&
          tunnelmidpos,tunnelmaxpos
    do i=1,natommol
      temp = sqrt((crd(2,i)-ctdy)**2 + (crd(3,i)-ctdz)**2)
      if ( crd(1,i) < tunnelminpos )  then
        if ( temp > cylinderradius + conemaxradius)  c1=c1+1
        if ( temp <= cylinderradius + conemaxradius .and. temp &
            > cylinderradius)  c2=c2+1
        if ( temp <= cylinderradius ) c3=c3+1
      else if ( crd(1,i) >= tunnelminpos .and. crd(1,i) < tunnelmidpos)  then
        if ( temp > cylinderradius + conemaxradius)  c4=c4+1
        if ( temp <= cylinderradius + conemaxradius .and. temp &
            > cylinderradius) c5=c5+1
        if ( temp <= cylinderradius ) c6=c6+1
      else if ( crd(1,i) >= tunnelmidpos .and. crd(1,i) < tunnelmaxpos)  then
        temp2=( crd(1,i)-tunnelmidpos)*conemaxradius/coneheight
        if ( temp > cylinderradius + conemaxradius) c7=c7+1
        if ( temp <= cylinderradius + conemaxradius .and. temp &
            > cylinderradius + temp2 ) c8=c8+1
        if ( temp > cylinderradius .and. temp <= cylinderradius + temp2 &
            ) c9=c9+1
        if ( temp <= cylinderradius ) c10=c10+1
      else
        if ( temp > cylinderradius + conemaxradius) c11=c11+1
        if ( temp <= cylinderradius + conemaxradius .and. temp &
            > cylinderradius) c12=c12+1
        if ( temp <= cylinderradius ) c13=c13+1
      end if
    end do

    write(*,'(a,4(i11))') "rejected counter:",lpanel_rejected,&
          cylinder_rejected,cone_rejected,rpanel_rejected

    ! '#' represent the wall of tunnel
    ! '-', '+' and '|' make nonsense
    write(*,*)   "+----------+--------------------------+---------#----------+"
    write(*,*)   "|          |                          |         #          |"
    write(*,100) "|  ",c1,"  |          ",c4,"          | ",c7,"  #   ",c11,"|"
    write(*,*)   "+----------+--------------------------+---------#----------+"
    write(*,200) "|          |                          |",c8," ##|          |"
    write(*,300) "|  ",c2,"  |          ",c5,"          |   ###   |   ",c12,"|"
    write(*,200) "|          |                          |## ",c9,"|          |"
    write(*,*)   "+----------############################---------+----------+"
    write(*,*)   "|          #                          |         |          |"
    write(*,400) "|  ",c3,"  #          ",c6,"          |  ",c10,"|   ",c13,"|"
    write(*,*)   "+----------#--------------------------+---------+----------+"

    if ( (c11+c12+c13) == natommol ) then 
      itunnel=0
      nscm=1000
      write(*,*) "************************************************************"
      write(*,*) "*                       CAUTION                            *"
      write(*,*) "*     Since the whole protein has got out of the exit      *"
      write(*,*) "*   tunnel, this tunnel will be turned off in the rest     *"
      write(*,*) "*   MD simulation.                                         *"
      write(*,*) "************************************************************"
    end if
  end if

100 format(1x,a,i6,a,i6,a,i6,a,i7,a)
200 format(1x,a,i6,a)
300 format(1x,a,i6,a,i6,a,i7,a)
400 format(1x,a,i6,a,i6,a,i7,a,i7,a)

end subroutine
 
!subroutine centroid_of_protein(atm_crd)  
!
!  implicit none
!  integer               :: i,j
!  real*8,intent(in)     :: atm_crd(3,natommol)
!
!  ctd(1:3)=(/ 0.0d0,0.0d0,0.0d0 /)
!  do i=1,natommol
!    do j=1,3
!      ctd(j)=ctd(j)+atm_crd(j,i)
!    end do
!  end do
!  do i=1,3
!    ctd(i)=ctd(i)/natommol
!  enddo
!
!end subroutine

subroutine normalize(vector)

  implicit none
  real*8,intent(inout):: vector(3)
  real*8 :: temp

  temp=sqrt(dot_product(vector,vector))
  vector(1:3)=vector(1:3)/temp

end subroutine

end module
