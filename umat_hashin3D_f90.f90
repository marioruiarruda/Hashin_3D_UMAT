!     ------------------------------------------------------------------
!     |     DATE           PROGRAMMER            DESCRIPTION OF CHANGE |  
!     |     ====           ==========            ===================== |
!     |  09/05/2024      MÁRIO RUI ARRUDA         SUBROUTINES ADDED    |
!     |  LISBON/TORINO                                                 |
!     ------------------------------------------------------------------

!     #COPYRIGHT 2022 BY MÁRIO RUI TIAGO ARRUDA
!     ALL RIGHTS RESERVED. NO PART OF THIS SUBROUTINE MAY BE REPRODUCED OR 
!     USED IN ANY MANNER WITHOUT WRITTEN PERMISSION OF THE COPYRIGHT OWNER.

!     If using this UMAT in future papers please cite https://doi.org/10.3390/app13021155
!     and give credit to the original authors of this UMAT.


!     --------UMAT 2D AND 3D ELEMENTS WITH HASHIN LINEAR DAMAGE---------

!     ------------------------------------------------------------------
!     ----------ABAQUS INPUT VARIABLES IN UMAT SUBROUTINE---------------
!     ------------------------------------------------------------------

subroutine umat(stress, statev, ddsdde, sse, spd, scd, rpl, &
           ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, &
           dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatv, props, & 
           nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, noel, &
           npt, layer, kspt, kstep, kinc)

IMPLICIT NONE

! Used this to debbug the UMAT, by commenting "include 'aba_param.inc'"
  
! include 'aba_param.inc'
! WARNING the aba_param.inc file declares Implicit real*8(a-h,o-z) !!!!!!!!!!!!!!!!!!!!!!!
! This means that, by default, any variables with first letter between a-h or o-z are real.
! The rest are integers.
! Note that this also means that if you type a variable name incorrectly, the compiler won't catch your typo. 
 

!     ------------------------------------------------------------------
!     -------------------ABAQUS DIMENSION VARIABLES---------------------
!     ------------------------------------------------------------------
 
!     -----------------ABAQUS CHARACTER VARIABLES----------------------
character(len=80) :: cmname 

!     -------------------ABAQUS INTEGER VARIABLES-----------------------  
integer :: ntens,ndi,nshr,nstatv,nprops,noel,npt,kspt,kstep,kinc,nprecd,layer

!     -------------------ABAQUS REAL VARIABLES--------------------------
real(kind=8) :: celent,sse,spd,scd,rpl,drpldt,dtime,temp,dtemp,pnewdt 

!     -------------------ABAQUS ARRAY VARIABLES-------------------------
real(kind=8) :: stress(ntens),statev(nstatv),ddsdde(ntens,ntens),&
                ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens),&
                predef(1),dpred(1),props(nprops),coords(3),drot(3,3),&
                dfgrd0(3,3),dfgrd1(3,3),time(2)

					 
!     ------------------------------------------------------------------    
!     -----------------DECLARATION OF VARIABLES-------------------------
!     ------------------------------------------------------------------

integer :: i,j,k,l,m,n,kiter,ktotal

real(kind=8) :: E1,E2,E3,G12,G13,G23,v12,v21,v13,v31,v23,v32,Xt,Xc,Yt,Yc,Zt,Zc,Sl,St,Si,&
								Gft,Gfc,Gmt,Gmc,Git,Gic,etaft,etafc,etamt,etamc,etait,etaic,dmax,&
								dfto,dfco,dmto,dmco,dito,dico,Fafto,Fafco,Famto,Famco,Faito,Faico,&
								dvft,dvfc,dvmt,dvmc,dvit,dvic,dvfto,dvfco,dvmto,dvmco,dvito,dvico,dvso,&
	              seqft0,seqfc0,seqmt0,seqmc0,seqit0,seqic0,alpha_f,alpha_i,&
								df,dvf,dm,dvm,di,dvi,Faft,Fafc,Famt,Famc,Fait,Faic,Dcoef,catLc,&
								ueqft,ueqft0,ueqftu,ueqfc,ueqfc0,ueqfcu,ueqmt,ueqmt0,ueqmtu,ueqmc,ueqmc0,ueqmcu,&
								ueqit,ueqit0,ueqitu,ueqic,ueqic0,ueqicu,dft,dfc,dmt,dmc,dit,dic,&
								ds12,ds13,ds23,ds12o,ds13o,ds23o,dvs12,dvs13,dvs23,dvs12o,dvs13o,dvs23o,&
                presft,presfc,presmt,presmc,presit,presic,presshear,Flagshear12,Flagshear13,Flagshear23,&
                slip12eq0,slip13eq0,slip23eq0,shear12eq0,shear13eq0,shear23eq0

!     ----------------------DUMMY VARIABLES-----------------------------
real(kind=8) :: const1,const2,const3,const4

!     -------------INITIATIONS OF ARRAYS-------------------------------- 
real(kind=8) :: eij(ntens),deij(ntens),eijo(ntens),sij(ntens),sijo(ntens), & 
                sije(ntens),sijeo(ntens),ddsddee(ntens,ntens),eij_e(ntens), &
                nsij(ntens),nsijo(ntens),nddsdde(ntens,ntens) !!! sij,sijo,ddsdde VISCOUS nsij,nsijo,nddsdde NON-VISCOUS

!     -------------DOUBLE PRECISION VALUES OF 0,1,2,3-------------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0, FOUR=4.d0,&
                           HALF=0.5d0,TOL0=1.d-8			
						
!     --------------ELASTIC AND MECHANICAL PROPERTIES-------------------

E1=props(1)     ! LONGITUDINAL ELASTIC MODULUS
E2=props(2)     ! TRANSVERSAL ELASTIC MODULUS
E3=props(3)     ! OUT-PLANE ELASTIC MODULUS

G12=props(4)    ! SHEAR ELASTIC MODULUS
G13=props(5)    ! SHEAR ELASTIC MODULUS
G23=props(6)    ! SHEAR ELASTIC MODULUS

v12=props(7)    ! POISSON COEFICIENT
v21=v12*E2/E1   ! SYMMETRIC POISSON COEFICIENT
v13=props(8)    ! POISSON COEFICIENT
v31=v13*E3/E1   ! SYMMETRIC POISSON COEFICIENT
v23=props(9)    ! POISSON COEFICIENT
v32=v23*E3/E2   ! SYMMETRIC POISSON COEFICIENT

Xt=props(10)    ! LONGITUDINAL TENSILE RESISTENCE
Xc=props(11)    ! LONGITUDINAL COMPRESSIVE RESISTENCE
Yt=props(12)    ! TRANSVERSAL TENSILE RESISTENCE
Yc=props(13)    ! TRANSVERSAL COMPRESSIVE RESISTENCE
Zt=props(14)    ! OUT-PLANE TENSILE RESISTENCE
Zc=props(15)    ! OUT-PLANE COMPRESSIVE RESISTENCE

Sl=props(16)    ! LONGITUDINAL SHEAR RESISTENCE
St=props(17)    ! TRANSVERSAL SHEAR RESISTENCE
Si=props(18)    ! OUT-PLANE SHEAR RESISTENCE

Gft=props(19)   ! FRACTURE ENERGY FOR FIBRE TENSION
Gfc=props(20)   ! FRACTURE ENERGY FOR FIBRE COMPRESSION
Gmt=props(21)   ! FRACTURE ENERGY FOR MATRIX TENSION
Gmc=props(22)   ! FRACTURE ENERGY FOR MATRIX COMPRESSION
Git=props(23)   ! FRACTURE ENERGY FOR INTERLAMINAR TENSION
Gic=props(24)   ! FRACTURE ENERGY FOR INTERLAMINAR COMPRESSION

etaft=props(25) ! VISCOUS REGULARIZATION COEFICIENT FOR FIBRE TENSION
etafc=etaft     ! VISCOUS REGULARIZATION COEFICIENT FOR FIBRE COMPRESSION
etamt=etaft     ! VISCOUS REGULARIZATION COEFICIENT FOR MATRIX TENSION
etamc=etaft     ! VISCOUS REGULARIZATION COEFICIENT FOR MATRIX COMPRESSION
etait=etaft     ! VISCOUS REGULARIZATION COEFICIENT FOR INTERLAMINAR TENSION
etaic=etaft     ! VISCOUS REGULARIZATION COEFICIENT FOR INTERLAMINAR COMPRESSION

dmax=props(26)    ! MAXIMUM ALLOWED DAMAGE FOR CONVERGENCY 
alpha_f=props(27) ! INTERACTION COEFICIENT FOR FIBRE TO SIMULATE 2D BEHAVIOUR
alpha_i=props(28) ! INTERACTION COEFICIENT FOR INTERLAMINAR TO SIMULATE 2D BEHAVIOUR

presft=props(29)    ! FIBRE TENSILE RESIDUAL STRESS
presfc=props(30)    ! FIBRE COMPRESSIVE RESIDUAL STRESS
presmt=props(31)    ! MATRIX TENSILE RESIDUAL STRESS
presmc=props(32)    ! MATRIX COMPRESSIVE RESIDUAL STRESS
presit=props(33)    ! INTERLAMINAR TENSILE RESIDUAL STRESS
presic=props(34)    ! INTERLAMINAR COMPRESSIVE RESIDUAL STRESS

!     ---------------STATE FIELD VARIABLES FOR ABAQUS-------------------
if (kinc == 1) then 
	do i=1, nstatv
		statev(i)=ZERO  ! VARIABLES INITIATION (if not=0, depending on compiler LINUX/WINDOWS)
  end do
end if

dfto=statev(1)      ! FIBRE TENSION DAMAGE FROM PREVIOUS INCREMENT
dfco=statev(2)      ! FIBRE COMPRESSION DAMAGE FROM PREVIOUS INCREMENT
dmto=statev(3)      ! MATRIX TENSION DAMAGE FROM PREVIOUS INCREMENT
dmco=statev(4)      ! MATRIX COMPRESSION DAMAGE FROM PREVIOUS INCREMENT
dito=statev(5)      ! INTERLAMINAR TENSION DAMAGE FROM PREVIOUS INCREMENT
dico=statev(6)      ! INTERLAMINAR COMPRESSION DAMAGE FROM PREVIOUS INCREMENT
ds12o=statev(7)     ! SHEAR 12 DAMAGE FROM PREVIOUS INCREMENT
ds13o=statev(8)     ! SHEAR 13 DAMAGE FROM PREVIOUS INCREMENT
ds23o=statev(9)     ! SHEAR 23 DAMAGE FROM PREVIOUS INCREMENT

Fafto=statev(10)    ! FIBRE TENSION CRITERIA FROM PREVIOUS INCREMENT
Fafco=statev(11)    ! FIBRE COMPRESSION CRITERIA FROM PREVIOUS INCREMENT
Famto=statev(12)    ! MATRIX TENSION CRITERIA FROM PREVIOUS INCREMENT
Famco=statev(13)    ! MATRIX COMPRESSION CRITERIA FROM PREVIOUS INCREMENT
Faito=statev(14)    ! INTERLAMINAR TENSION CRITERIA FROM PREVIOUS INCREMENT
Faico=statev(15)    ! INTERLAMINAR COMPRESSION CRITERIA FROM PREVIOUS INCREMENT

dvfto=statev(16)    ! VISCOUS FIBRE TENSION DAMAGE FROM PREVIOUS INCREMENT
dvfco=statev(17)    ! VISCOUS FIBRE COMRESSION DAMAGE FROM PREVIOUS INCREMENT
dvmto=statev(18)    ! VISCOUS MATRIX TENSION DAMAGE FROM PREVIOUS INCREMENT
dvmco=statev(19)    ! VISCOUS MATRIX COMPRESSION DAMAGE FROM PREVIOUS INCREMENT
dvito=statev(20)    ! VISCOUS INTERLAMINAR TENSION DAMAGE FROM PREVIOUS INCREMENT
dvico=statev(21)    ! VISCOUS INTERLAMINAR COMPRESSION DAMAGE FROM PREVIOUS INCREMENT
dvs12o=statev(22)   ! VISCOUS SHEAR 12 DAMAGE FROM PREVIOUS INCREMENT
dvs13o=statev(23)   ! VISCOUS SHEAR 13 DAMAGE FROM PREVIOUS INCREMENT
dvs23o=statev(24)   ! VISCOUS SHEAR 23 DAMAGE FROM PREVIOUS INCREMENT

seqft0=statev(25)   ! EQUIVALENTE FIBRE TENSION AT THE ONSET OF DAMAGE FROM PREVIOUS INCREMENT
seqfc0=statev(26)   ! EQUIVALENTE FIBRE COMPRESSION AT THE ONSET OF DAMAGE FROM PREVIOUS INCREMENT
seqmt0=statev(27)   ! EQUIVALENTE MATRIX TENSION AT THE ONSET OF DAMAGE FROM PREVIOUS INCREMENT
seqmc0=statev(28)   ! EQUIVALENTE MATRIX COMPRESSION AT THE ONSET OF DAMAGE FROM PREVIOUS INCREMENT
seqit0=statev(29)   ! EQUIVALENTE INTERLAMINAR TENSION AT THE ONSET OF DAMAGE FROM PREVIOUS INCREMENT
seqic0=statev(30)   ! EQUIVALENTE INTERLAMINAR COMPRESSION AT THE ONSET OF DAMAGE FROM PREVIOUS INCREMENT

ueqft0=statev(31)   ! EQUIVALENTE FIBRE TENSION DISPLACEMENT AT THE ONSET OF DAMAGE FROM PREVIOUS INCREMENT
ueqfc0=statev(32)   ! EQUIVALENTE FIBRE COMPRESSION DISPLACEMENT AT THE ONSET OF DAMAGE FROM PREVIOUS INCREMENT
ueqmt0=statev(33)   ! EQUIVALENTE MATRIX TENSION DISPLACEMENT AT THE ONSET OF DAMAGE FROM PREVIOUS INCREMENT
ueqmc0=statev(34)   ! EQUIVALENTE MATRIX COMPRESSON DISPLACEMENT AT THE ONSET OF DAMAGE FROM PREVIOUS INCREMENT
ueqit0=statev(35)   ! EQUIVALENTE INTERLAMINAR TENSION DISPLACEMENT AT THE ONSET OF DAMAGE FROM PREVIOUS INCREMENT
ueqic0=statev(36)   ! EQUIVALENTE INTERLAMINAR COMPRESSON DISPLACEMENT AT THE ONSET OF DAMAGE FROM PREVIOUS INCREMENT

if (ndi==3) then !!!! FOR 3D ANALYSIS
sijeo(1)=statev(37) ! EFFECTIVE STRESS 11 FROM PREVIOUS INCREMENT
sijeo(2)=statev(38) ! EFFECTIVE STRESS 22 FROM PREVIOUS INCREMENT
sijeo(3)=statev(39) ! EFFECTIVE STRESS 33 FROM PREVIOUS INCREMENT
sijeo(4)=statev(40) ! EFFECTIVE STRESS 12 FROM PREVIOUS INCREMENT
sijeo(5)=statev(41) ! EFFECTIVE STRESS 13 FROM PREVIOUS INCREMENT
sijeo(6)=statev(42) ! EFFECTIVE STRESS 23 FROM PREVIOUS INCREMENT
else             !!!! FOR 2D ANALYSIS
sijeo(1)=statev(37) ! EFFECTIVE STRESS 11 FROM PREVIOUS INCREMENT
sijeo(2)=statev(38) ! EFFECTIVE STRESS 22 FROM PREVIOUS INCREMENT
sijeo(3)=statev(39) ! EFFECTIVE STRESS 12 FROM PREVIOUS INCREMENT
end if


!     ------------------------------------------------------------------
!     ---------FIELD PREDICTOR FOR VARIABLES AND STIFFNESS -------------
!     ------------------------------------------------------------------

!     --------------FINAL STRAIN INCREMENT------------------------------

do i=1, ntens
  eij(i)=stran(i)+dstran(i)       ! TOTAL STRAIN IN THE BEGGIND OF THE INCREMENT
  eijo(i)=stran(i)                ! STRAIN FROM PREVIOUS INCREMENT
	deij(i)=dstran(i)               ! STRAIN INCREMENT
  sijo(i)=stress(i)               ! STRESS FROM PREVIOUS INCREMENT
  !const4=ZERO                   ! FULL EXPLICIT INTEGRATION
  !const4=ONE                    ! FULL IMPLICIT INTEGRATION
  const4=dtime/(dtime+etaft)      ! FROM IMPLICIT TO EXPLICIT WHEN dtime<<eta
  eij_e(i)=stran(i)+const4*dstran(i) ! FOR EULER INTEGRATION (ZERO,HALF,ONE)
end do

!     ----STABILITY CONDITION DURING EXPLICIT MATERIAL INTEGRATION------
! do i=1,ntens
	! if (abs(dstran(i)) > (Xt/E1)) then
		! pnewdt=0.5
	! end if
! end do

!     -------VERIFICATION OF FIBRE AND MATRIX TENSION/COMPRESSION-------

call calc_damage_verification(ndi,ntens,sijeo,df,dm,di,dvf,dvm,dvi,dfto,dfco,dmto,dmco,dito,dico,dvfto,dvfco,dvmto,dvmco,dvito,dvico)

Faft=Fafto ! INITIATION FIBRE TENSION FAILURE
Fafc=Fafco ! INITIATION FIBRE COMPRESSION FAILURE
Famt=Famto ! INITIATION MATRIX TENSION FAILURE
Famc=Famco ! INITIATION MATRIX COMPRESSION FAILURE
Fait=Faito ! INITIATION INTERLAMINAR TENSION FAILURE
Faic=Faico ! INITIATION INTERLAMINAR COMPRESSION FAILURE

!     -----SECANT DAMAGED MATRIX AT THE BEGINNING OF THE INCREMENT-------
call calc_stiffness(ntens,ndi,ddsdde,E1,E2,E3,G12,G13,G23,v12,v21,v23,v32,v13,v31,df,dm,di,ds12o,ds23o,ds13o)

!     -----SECANT EFFECTIVE MATRIX AT THE BEGINNING OF THE INCREMENT-------
call calc_stiffness_effective(ntens,ndi,ddsddee,E1,E2,E3,G12,G13,G23,v12,v21,v23,v32,v13,v31,df,dm,di)

!    ---------CALCULATION OF THE STRESS FROM THE PREVIOUS INCREMENT--------
call calc_stress(ndi,ntens,nsijo,eijo,ddsdde)

!    ---------CALCULATION OF THE STRESS AT THE BEGINNING OF THE INCREMENT-------
call calc_stress(ndi,ntens,sij,eij_e,ddsdde)

!    ---CALCULATION OF THE EFFECTIVE STRESS AT THE BEGINNING OF THE INCREMENT---
call calc_stress(ndi,ntens,sije,eij_e,ddsddee)


!     ------------------------------------------------------------------
!     --------------FAILURE CRITERIA VERIFICATION-----------------------
!     ------------------------------------------------------------------

!     -------------CARACTERISTIC LENGTH FROM ABAQUS---------------------
catLc=celent

!     --------VERIFICATION OF FIBRE TENSION/COMPRESSION FAILURE---------
call calc_fibre_failure(nstatv,ntens,ndi,catLc,Faft,Fafc,eijo,sijo,sije,sijeo,statev,ueqft0,ueqfc0,seqft0,seqfc0,Xt,Xc,Sl,St,Si,alpha_f)

!     --------VERIFICATION OF MATRIX TENSION/COMPRESSION FAILURE---------
call calc_matrix_failure(nstatv,ntens,ndi,catLc,Famt,Famc,eijo,sijo,sije,sijeo,statev,ueqmt0,ueqmc0,seqmt0,seqmc0,Yt,Yc,Sl,St,Si,ZERO)

!     -----VERIFICATION OF INTERLAMINAR TENSION/COMPRESSION FAILURE------
if (ndi==3) then
call calc_interlaminar_failure(nstatv,ntens,ndi,catLc,Fait,Faic,eijo,sijo,sije,sijeo,statev,ueqit0,ueqic0,seqit0,seqic0,Zt,Zc,Sl,St,Si,alpha_i)
end if


!     ------------------------------------------------------------------
!     -----------------HASHIN DAMAGE EVOLUTION--------------------------
!     ------------------------------------------------------------------

dft=dfto;dfc=dfco;dmt=dmto;dmc=dmco;dit=dito;dic=dico
ds12=ds12o;ds13=ds13o;ds23=ds23o;dvs12=dvs12o;dvs13=dvs13o;dvs23=dvs23o
dvft=dvfto;dvfc=dvfco;dvmt=dvmto;dvmc=dvmco;dvit=dvito;dvic=dvico

!     ------------------FIBER DAMAGE EVOLUTION--------------------------
call calc_fibre_damage(statev,eij_e,nstatv,ntens,ndi,catLc,Faft,Fafc,alpha_f,Gft,Gfc,ueqft0,ueqfc0, & 
                       seqft0,seqfc0,dft,dfc,dvft,dvfc,dmax,dfto,dfco,dvfto,dvfco,etaft,etafc,dtime,ueqft, &
                       ueqftu,ueqfc,ueqfcu,presft,presfc)

!     ------------------MATRIX DAMAGE EVOLUTION-------------------------
call calc_matrix_damage(statev,eij_e,nstatv,ntens,ndi,catLc,Famt,Famc,ZERO,Gmt,Gmc,ueqmt0,ueqmc0, & 
                        seqmt0,seqmc0,dmt,dmc,dvmt,dvmc,dmax,dmto,dmco,dvmto,dvmco,etamt,etamc,dtime,ueqmt, &
                        ueqmtu,ueqmc,ueqmcu,presmt,presmc)
                        
!     ------------------INTERLAMINAR DAMAGE EVOLUTION-------------------
if (ndi==3) then
call calc_interlaminar_damage(statev,eij_e,nstatv,ntens,ndi,catLc,Fait,Faic,alpha_i,Git,Gic,ueqit0,ueqic0, & 
                              seqit0,seqic0,dit,dic,dvit,dvic,dmax,dito,dico,dvito,dvico,etait,etaic,dtime,ueqit, &
                              ueqitu,ueqic,ueqicu,presit,presic)
end if                              
                              
!     --------------------SHEAR DAMAGE EVOLUTION------------------------
ds12=ONE-(ONE-dft)*(ONE-dfc)*(ONE-dmt)*(ONE-dmc)
ds13=ONE-(ONE-dft)*(ONE-dfc)*(ONE-dit)*(ONE-dic)
ds23=ONE-(ONE-dit)*(ONE-dic)*(ONE-dmt)*(ONE-dmc)

dvs12=ONE-(ONE-dvft)*(ONE-dvfc)*(ONE-dvmt)*(ONE-dvmc)
dvs13=ONE-(ONE-dvft)*(ONE-dvfc)*(ONE-dvit)*(ONE-dvic)
dvs23=ONE-(ONE-dvit)*(ONE-dvic)*(ONE-dvmt)*(ONE-dvmc)

statev(7)=ds12
statev(8)=ds13
statev(9)=ds23
statev(22)=dvs12
statev(23)=dvs13
statev(24)=dvs23


!     ------------------------------------------------------------------
!     -----------------FINAL DAMAGE BEHAVIOUR---------------------------
!     ------------------------------------------------------------------

call calc_damage_verification(ndi,ntens,sije,df,dm,di,dvf,dvm,dvi,dft,dfc,dmt,dmc,dit,dic,dvft,dvfc,dvmt,dvmc,dvit,dvic)


!     ------------------------------------------------------------------
!     -------FIELD CORRECTOR FOR STRESS AND SECANT STIFFNESS -----------
!     ------------------------------------------------------------------

!     -----VISCOUS SECANT DAMAGED MATRIX FOR THE ITERATION PROCESS------
call calc_stiffness(ntens,ndi,ddsdde,E1,E2,E3,G12,G13,G23,v12,v21,v23,v32,v13,v31,dvf,dvm,dvi,dvs12,dvs23,dvs13)

!     ---------SECANT DAMAGED MATRIX FOR THE ITERATION PROCESS----------
call calc_stiffness(ntens,ndi,nddsdde,E1,E2,E3,G12,G13,G23,v12,v21,v23,v32,v13,v31,df,dm,di,ds12,ds23,ds13)

!     ----------SECANT EFFECTIVE MATRIX IN THE ACTUAL LOAD STEP---------
call calc_stiffness_effective(ntens,ndi,ddsddee,E1,E2,E3,G12,G13,G23,v12,v21,v23,v32,v13,v31,df,dm,di)

!    -------------------------STRESS CORRECTOR--------------------------
call calc_stress(ndi,ntens,nsij,eij,nddsdde)

!    ---------------------VISCOUS STRESS CORRECTOR----------------------
call calc_stress(ndi,ntens,stress,eij,ddsdde)

!    ----------------VISCOUS EFFECTIVE STRESS CORRECTOR-----------------
call calc_stress(ndi,ntens,sije,eij,ddsddee)

!     ---------------------UPDATED STRAIN ENERGY------------------------
call calc_energy(ndi,ntens,sse,stress,sijo,dstran,scd,nsij,nsijo)

if (ndi==3) then 
statev(37)=sije(1);statev(38)=sije(2);statev(39)=sije(3);statev(40)=sije(4);statev(41)=sije(5);statev(42)=sije(6)
else
statev(37)=sije(1);statev(38)=sije(2);statev(39)=sije(3)
end if

return
end




!     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!     ||||||||||||||||||||||||USER SUBROUTINES||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


!     __________________________________________________________________
!     ________________SUBROUTINE DAMAGE VERIFICATION____________________
!     __________________________________________________________________
subroutine calc_damage_verification(ndi,ntens,sij,df,dm,di,dvf,dvm,dvi,dft,dfc,dmt,dmc,dit,dic,dvft,dvfc,dvmt,dvmc,dvit,dvic)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer, intent(in):: ntens,ndi

!     ------------------------REAL VARIABLES---------------------------- 
real(kind=8), intent(in)  :: sij(ntens),dft,dfc,dmt,dmc,dit,dic,dvft,dvfc,dvmt,dvmc,dvit,dvic
real(kind=8), intent(inout)  :: df,dm,di,dvf,dvm,dvi

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

!     -----------------------FIBRE DAMAGE-------------------------------
if (sij(1) >= ZERO) then
	df=dft;dvf=dvft
else
	df=dfc;dvf=dvfc
endif
		
!     -----------------------MATRIX DAMAGE------------------------------  
if (ndi==3) then !!!!!!!!!!!! FOR 3D ANALYSIS 
  if ((sij(2)+sij(3)) >= ZERO) then
    dm=dmt;dvm=dvmt
  else
    dm=dmc;dvm=dvmc
  endif
!     --------------------INTERLAMINAR DAMAGE---------------------------
  if (sij(3) >= ZERO) then
    di=dit;dvi=dvit
  else
    di=dic;dvi=dvic
  endif
else             !!!!!!!!!!!! FOR 2D ANALYSIS
!     -----------------------MATRIX DAMAGE------------------------------ 
  if (sij(2) >= ZERO) then
    dm=dmt;dvm=dvmt
  else
    dm=dmc;dvm=dvmc
  endif
end if

return
end subroutine


!     __________________________________________________________________
!     ________________SUBROUTINE STIFFNESS ASSEMBLER____________________
!     __________________________________________________________________

subroutine calc_stiffness(ntens,ndi,Cij,E1,E2,E3,G12,G13,G23,v12,v21,v23,v32,v13,v31,df,dm,di,ds12,ds23,ds13)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer, intent(in):: ntens,ndi

!     ------------------------REAL VARIABLES---------------------------- 
real(kind=8) :: Dvar
real(kind=8), intent(in)  :: E1,E2,E3,G12,G13,G23,v12,v21,v23,v32,v13,v31,df,dm,di,ds12,ds23,ds13
real(kind=8), intent(inout) :: Cij(ntens,ntens)

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

if (ndi==3) then ! FOR 3D SOLID ANALYSIS
  Dvar=ONE/(ONE-(ONE-df)*(ONE-dm)*v12*v21-(ONE-dm)*(ONE-di)*v23*v32-(ONE-di)*(ONE-df)*v31*v13-TWO*(ONE-df)*(ONE-dm)*(ONE-di)*v12*v23*v31)
  Cij(1,1)=(ONE-df)*E1*(ONE-(ONE-dm)*(ONE-di)*v23*v32)*Dvar
  Cij(1,2)=(ONE-df)*(ONE-dm)*E1*(v21+(ONE-di)*v31*v23)*Dvar
  Cij(1,3)=(ONE-df)*(ONE-di)*E1*(v31+(ONE-dm)*v21*v32)*Dvar
  Cij(2,1)=(ONE-df)*(ONE-dm)*E2*(v12+(ONE-di)*v13*v32)*Dvar
  Cij(2,2)=(ONE-dm)*E2*(ONE-(ONE-df)*(ONE-di)*v31*v13)*Dvar
  Cij(2,3)=(ONE-dm)*(ONE-di)*E2*(v32+(ONE-df)*v31*v12)*Dvar
  Cij(3,1)=(ONE-df)*(ONE-di)*E3*(v13+(ONE-dm)*v12*v23)*Dvar
  Cij(3,2)=(ONE-dm)*(ONE-di)*E3*(v23+(ONE-df)*v13*v21)*Dvar
  Cij(3,3)=(ONE-di)*E3*(ONE-(ONE-df)*(ONE-dm)*v12*v21)*Dvar
  Cij(4,4)=(ONE-ds12)*G12
  Cij(5,5)=(ONE-ds13)*G13
  Cij(6,6)=(ONE-ds23)*G23
  Cij(4,1)=ZERO;Cij(4,2)=ZERO;Cij(4,3)=ZERO;Cij(4,5)=ZERO;Cij(4,6)=ZERO
  Cij(5,1)=ZERO;Cij(5,2)=ZERO;Cij(5,3)=ZERO;Cij(5,4)=ZERO;Cij(5,6)=ZERO
  Cij(6,1)=ZERO;Cij(6,2)=ZERO;Cij(6,3)=ZERO;Cij(6,4)=ZERO;Cij(6,5)=ZERO
else            ! FOR 2D SOLID ANALYSIS
  Dvar=ONE-(ONE-df)*(ONE-dm)*v12*v21 
  Cij(1,1)=(ONE-df)*E1*Dvar
  Cij(2,2)=(ONE-dm)*E2*Dvar
  Cij(3,3)=(ONE-ds12)*G12 
  Cij(1,2)=(ONE-df)*(ONE-dm)*v21*E1*Dvar
  Cij(2,1)=(ONE-df)*(ONE-dm)*v12*E2*Dvar
  Cij(1,3)=ZERO;Cij(3,1)=ZERO
  Cij(2,3)=ZERO;Cij(3,2)=ZERO
end if

return
end subroutine


!     __________________________________________________________________
!     ___________SUBROUTINE EFFECTIVE STIFFNESS ASSEMBLER_______________
!     __________________________________________________________________

subroutine calc_stiffness_effective(ntens,ndi,Cij,E1,E2,E3,G12,G13,G23,v12,v21,v23,v32,v13,v31,df,dm,di)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer, intent(in):: ntens,ndi

!     ------------------------REAL VARIABLES---------------------------- 
real (kind=8) :: Dvar
real(kind=8), intent(in)  :: E1,E2,E3,G12,G13,G23,v12,v21,v23,v32,v13,v31,df,dm,di
real(kind=8), intent(inout) :: Cij(ntens,ntens)

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

if (ndi==3) then ! FOR 3D SOLID ANALYSIS
  Dvar=ONE/(ONE-(ONE-df)*(ONE-dm)*v12*v21-(ONE-dm)*(ONE-di)*v23*v32-(ONE-di)*(ONE-df)*v31*v13-TWO*(ONE-df)*(ONE-dm)*(ONE-di)*v12*v23*v31)
  Cij(1,1)=E1*(ONE-(ONE-dm)*(ONE-di)*v23*v32)*Dvar
  Cij(1,2)=(ONE-dm)*E1*(v21+(ONE-di)*v31*v23)*Dvar
  Cij(1,3)=(ONE-di)*E1*(v31+(ONE-dm)*v21*v32)*Dvar
  Cij(2,1)=(ONE-df)*E2*(v12+(ONE-di)*v13*v32)*Dvar
  Cij(2,2)=E2*(ONE-(ONE-df)*(ONE-di)*v31*v13)*Dvar
  Cij(2,3)=(ONE-di)*E2*(v32+(ONE-df)*v31*v12)*Dvar
  Cij(3,1)=(ONE-df)*E3*(v13+(ONE-dm)*v12*v23)*Dvar
  Cij(3,2)=(ONE-dm)*E3*(v23+(ONE-df)*v13*v21)*Dvar
  Cij(3,3)=E3*(ONE-(ONE-df)*(ONE-dm)*v12*v21)*Dvar
  Cij(4,4)=G12
  Cij(5,5)=G13
  Cij(6,6)=G23
  Cij(4,1)=ZERO;Cij(4,2)=ZERO;Cij(4,3)=ZERO;Cij(4,5)=ZERO;Cij(4,6)=ZERO
  Cij(5,1)=ZERO;Cij(5,2)=ZERO;Cij(5,3)=ZERO;Cij(5,4)=ZERO;Cij(5,6)=ZERO
  Cij(6,1)=ZERO;Cij(6,2)=ZERO;Cij(6,3)=ZERO;Cij(6,4)=ZERO;Cij(6,5)=ZERO
else            ! FOR 2D SOLID ANALYSIS
  Dvar=ONE-(ONE-df)*(ONE-dm)*v12*v21 
  Cij(1,1)=E1*Dvar
  Cij(2,2)=E2*Dvar
  Cij(3,3)=G12
  Cij(1,2)=(ONE-dm)*v21*E1*Dvar
  Cij(2,1)=(ONE-df)*v12*E2*Dvar
  Cij(1,3)=ZERO;Cij(3,1)=ZERO
  Cij(2,3)=ZERO;Cij(3,2)=ZERO
end if

return
end subroutine


!     __________________________________________________________________
!     ___________SUBROUTINE FOR MACAULY BRACKET OPERATION_______________
!     __________________________________________________________________

subroutine calc_macaulay_bracket(ndi,mat,mat_b,signal)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer :: i
integer, intent(in):: ndi

!     ------------------------REAL VARIABLES---------------------------- 
real(kind=8), intent(in)  :: mat(ndi),signal
real(kind=8), intent(inout) :: mat_b(ndi)

!     -------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2---------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

do i=1,ndi
	if (signal>=ZERO) then
		mat_b(i)=abs(max(ZERO,mat(i)))
	else
		mat_b(i)=abs(min(ZERO,mat(i)))
	end if
end do

return
end subroutine


!     __________________________________________________________________
!     ___________SUBROUTINE FOR MACAULY BRACKET OPERATION_______________
!     __________________________________________________________________

subroutine calc_macaulay_bracket_scalar(x,x_b,signal)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer :: i

!     ------------------------REAL VARIABLES---------------------------- 
real(kind=8), intent(in)  :: x,signal
real(kind=8), intent(inout) :: x_b

!     -------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2---------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

if (signal>=ZERO) then
	x_b=abs(max(ZERO,x))
else
	x_b=abs(min(ZERO,x))
end if

return
end subroutine


!     __________________________________________________________________
!     _________SUBROUTINE FOR STRAIN STIFFNESS MULTIPLICATION___________
!     __________________________________________________________________

subroutine calc_stress(ndi,ntens,mat_si,mat_ej,mat_Cij)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer :: i,j
integer, intent(in):: ndi,ntens

!     ------------------------REAL VARIABLES---------------------------- 
real(kind=8), intent(in)  :: mat_Cij(ntens,ntens), mat_ej(ntens)
real(kind=8), intent(inout) :: mat_si(ntens)

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

do i=1,ndi
  mat_si(i)=ZERO
  do j=1,ndi ! NORMAL STRAIN
	  mat_si(i)=mat_si(i)+mat_Cij(i,j)*mat_ej(j)
  end do
end do
do i=ndi+1,ntens ! SHEAR STRAIN
  mat_si(i)=mat_Cij(i,i)*mat_ej(i) 
end do

return
end subroutine


!     __________________________________________________________________
!     _______________SUBROUTINE FOR DEFORMATION ENERGY__________________
!     __________________________________________________________________

subroutine calc_energy(ndi,ntens,energy,sij,sij_old,deij,nenergy,nsij,nsij_old)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer :: i,j
integer, intent(in):: ndi,ntens

!     ------------------------REAL VARIABLES---------------------------- 
real(kind=8), intent(in)  :: sij(ntens), sij_old(ntens), deij(ntens),&
                             nsij(ntens), nsij_old(ntens)
real(kind=8), intent(inout) :: energy,nenergy

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0

do i=1,ntens
  energy=energy+HALF*(sij_old(i)+sij(i))*deij(i)
  nenergy=nenergy+HALF*(sij_old(i)+sij(i)-nsij_old(i)-nsij(i))*deij(i)
end do

return
end subroutine


!     __________________________________________________________________
!     _______________SUBROUTINE FIBRE FAILURE CRITERIA__________________
!     __________________________________________________________________

subroutine calc_fibre_failure(nstatv,ntens,ndi,catLc,Fat,Fac,eijo,sijo,sije,sijeo,statev,ueqt0,ueqc0,seqt0,seqc0,Rt,Rc,Tl,Ts,Ti,a)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer :: i,j
integer, intent(in):: ntens,ndi,nstatv

!     ------------------------REAL VARIABLES---------------------------- 
real (kind=8) :: s11e,s12e,s13e,s11eo,s12eo,s13eo,s11o,s12o,s13o,e11o,e12o,e13o,&
                 s11obp,e11obp,s11obn,e11obn,slip12eq0,slip13eq0,shear12eq0,shear13eq0

real(kind=8), intent(in)  :: Rt,Rc,Tl,Ts,Ti,a,catLc
real(kind=8), intent(inout) :: statev(nstatv),sijo(ntens),eijo(ntens),sijeo(ntens),sije(ntens), &
                               ueqt0,ueqc0,seqt0,seqc0,Fat,Fac

!     ----------------------DUMMY VARIABLES-----------------------------
real(kind=8) :: const1,const2,const3,const4,const5,const6

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0, mONE=-1.d0

if (ndi == 3) then
  s11e=sije(1);s12e=sije(4);s13e=sije(5)
  s11eo=sijeo(1);s12eo=sijeo(4);s13eo=sijeo(5)
  s11o=sijo(1);s12o=sijo(4);s13o=sijo(5)
  e11o=eijo(1);e12o=eijo(4);e13o=eijo(5)
else
  s11e=sije(1);s12e=sije(3);s13e=ZERO
  s11eo=sijeo(1);s12eo=sijeo(3);s13eo=ZERO
  s11o=sijo(1);s12o=sijo(3);s13o=ZERO
  e11o=eijo(1);e12o=eijo(3);e13o=ZERO
end if

call calc_macaulay_bracket_scalar(s11o,s11obp,ONE)
call calc_macaulay_bracket_scalar(e11o,e11obp,ONE)
call calc_macaulay_bracket_scalar(s11o,s11obn,mONE)
call calc_macaulay_bracket_scalar(e11o,e11obn,mONE)

if (s11e >= ZERO .and. Fat < ONE) then
	Fat=(s11e/Rt)**2+a*(s12e**2+s13e**2)/Tl**2
	if (Fat >= ONE) then
	  Fat=(s11eo/Rt)**2+a*(s12eo**2+s13eo**2)/Tl**2
		const1=ONE/sqrt(Fat)
    
		ueqt0=const1*sqrt(e11obp**2+a*(e12o**2+e13o**2))*catLc
		seqt0=const1**2*(s11obp*e11obp+a*(s12o*e12o+s13o*e13o))/(ueqt0/catLc)
    
		statev(31)=ueqt0
		statev(25)=seqt0
    
		Fat=1.01d0
	end if
	statev(10)=Fat
else if (s11e < ZERO .and. Fac < ONE) then
	Fac=(s11e/Rc)**2
	if (Fac>=ONE) then
	  Fac=(s11eo/Rc)**2
		const1=ONE/sqrt(Fac)
		ueqc0=const1*e11obn*catLc
		seqc0=const1*s11obn
    
		statev(32)=ueqc0
		statev(26)=seqc0
    
		Fac=1.01d0
	end if
	statev(11)=Fac
end if

return
end subroutine


!     __________________________________________________________________
!     _______________SUBROUTINE MATRIX FAILURE CRITERIA__________________
!     __________________________________________________________________

subroutine calc_matrix_failure(nstatv,ntens,ndi,catLc,Fat,Fac,eijo,sijo,sije,sijeo,statev,ueqt0,ueqc0,seqt0,seqc0,Rt,Rc,Tl,Ts,Ti,a)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer :: i,j
integer, intent(in):: ntens,ndi,nstatv

!     ------------------------REAL VARIABLES---------------------------- 
real (kind=8) :: s22e,s33e,s12e,s13e,s23e,s22eo,s33eo,s12eo,s13eo,s23eo, &
                 s22o,s33o,s12o,s13o,s23o,e22o,e33o,e12o,e13o,e23o,&
                 s22obp,s33obp,e22obp,e33obp,s22obn,s33obn,e22obn,e33obn,&
                 slip12eq0,slip13eq0,slip23eq0,shear12eq0,shear13eq0,shear23eq0

real(kind=8), intent(in)  :: Rt,Rc,Tl,Ts,Ti,a,catLc
real(kind=8), intent(inout) :: statev(nstatv),sijo(ntens),eijo(ntens),sijeo(ntens),sije(ntens), &
                               ueqt0,ueqc0,seqt0,seqc0,Fat,Fac

!     ----------------------DUMMY VARIABLES-----------------------------
real(kind=8) :: const1,const2,const3,const4,const5,const6

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0, mONE=-1.d0

if (ndi == 3) then
  s22e=sije(2);s33e=sije(3);s12e=sije(4);s13e=sije(5);s23e=sije(6)
  s22eo=sijeo(2);s33eo=sijeo(3);s12eo=sijeo(4);s13eo=sijeo(5);s23eo=sijeo(6)
  s22o=sijo(2);s33o=sijo(3);s12o=sijo(4);s13o=sijo(5);s23o=sijo(6)
  e22o=eijo(2);e33o=eijo(3);e12o=eijo(4);e13o=eijo(5);e23o=eijo(6)
else
  s22e=sije(2);s33e=ZERO;s12e=sije(3);s13e=ZERO;s23e=ZERO
  s22eo=sijeo(2);s33eo=ZERO;s12eo=sijeo(3);s13eo=ZERO;s23eo=ZERO
  s22o=sijo(2);s33o=ZERO;s12o=sijo(3);s13o=ZERO;s23o=ZERO
  e22o=eijo(2);e33o=ZERO;e12o=eijo(3);e13o=ZERO;e23o=ZERO
end if

call calc_macaulay_bracket_scalar(s22o,s22obp,ONE)
call calc_macaulay_bracket_scalar(s33o,s33obp,ONE)
call calc_macaulay_bracket_scalar(e22o,e22obp,ONE)
call calc_macaulay_bracket_scalar(e33o,e33obp,ONE)
call calc_macaulay_bracket_scalar(s22o,s22obn,mONE)
call calc_macaulay_bracket_scalar(s33o,s33obn,mONE)
call calc_macaulay_bracket_scalar(e22o,e22obn,mONE)
call calc_macaulay_bracket_scalar(e33o,e33obn,mONE)

if ((s22e+s33e) >= ZERO .and. Fat < ONE) then
	Fat=((s22e+s33e)/Rt)**2+(s23e**2-s22e*s33e)/Ti**2+(s12e**2+s13e**2)/Tl**2
	if (Fat >= ONE) then
	  Fat=((s22eo+s33eo)/Rt)**2+(s23eo**2-s22eo*s33eo)/Ti**2+(s12eo**2+s13eo**2)/Tl**2
		const1=ONE/sqrt(Fat)
    
		ueqt0=const1*sqrt(e22obp**2+e33obp**2+e12o**2+e13o**2+e23o**2)*catLc
		seqt0=const1**2*(s22obp*e22obp+s33obp*e33obp+s12o*e12o+s13o*e13o+s23o*e23o)/(ueqt0/catLc)
    
		statev(33)=ueqt0
		statev(27)=seqt0
    
		Fat=1.01d0
	end if
	statev(12)=Fat
else if ((s22e+s33e) < ZERO .and. Fac < ONE) then
	Fac=((Rc*HALF/Ts)**2-ONE)*(s22e+s33e)/Rc+((s22e+s33e)*HALF/Ts)**2+(s23e**2-s22e*s33e)/Ti**2+(s12e**2+s13e**2)/Tl**2
	if (Fac>=ONE) then
		const2=((s22eo+s33eo)*HALF/Ts)**2+(s23eo**2-s22eo*s33eo)/Ti**2+(s12eo**2+s13eo**2)/Tl**2
		const3=((Rc*HALF/Ts)**2-ONE)*(s22eo+s33eo)/Rc
		const1=(-const3+sqrt(const3**2+4*const2))*HALF/const2
    
		ueqc0=const1*sqrt(e22obn**2+e33obn**2+e12o**2+e13o**2+e23o**2)*catLc
		seqc0=const1**2*(s22obn*e22obn+s33obn*e33obn+s12o*e12o+s13o*e13o+s23o*e23o)/(ueqc0/catLc)
    
		statev(34)=ueqc0
		statev(28)=seqc0
    
		Fac=1.01d0
	end if
	statev(13)=Fac
end if

return
end subroutine


!     __________________________________________________________________
!     ___________SUBROUTINE INTERLAMINAR FAILURE CRITERIA_______________
!     __________________________________________________________________

subroutine calc_interlaminar_failure(nstatv,ntens,ndi,catLc,Fat,Fac,eijo,sijo,sije,sijeo,statev,ueqt0,ueqc0,seqt0,seqc0,Rt,Rc,Tl,Ts,Ti,a)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer :: i,j
integer, intent(in):: ntens,ndi,nstatv

!     ------------------------REAL VARIABLES---------------------------- 
real (kind=8) :: s33e,s13e,s23e,s33eo,s13eo,s23eo, &
                 s33o,s13o,s23o,e33o,e13o,e23o,&
                 s33obp,e33obp,s33obn,e33obn,slip13eq0,slip23eq0,shear13eq0,shear23eq0

real(kind=8), intent(in)  :: Rt,Rc,Tl,Ts,Ti,a,catLc
real(kind=8), intent(inout) :: statev(nstatv),sijo(ntens),eijo(ntens),sijeo(ntens),sije(ntens), &
                               ueqt0,ueqc0,seqt0,seqc0,Fat,Fac

!     ----------------------DUMMY VARIABLES-----------------------------
real(kind=8) :: const1,const2,const3,const4,const5,const6

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0, mONE=-1.d0

if (ndi == 3) then
  s33e=sije(3);s13e=sije(5);s23e=sije(6) 
  s33eo=sijeo(3);s13eo=sijeo(5);s23eo=sijeo(6)
  s33o=sijo(3);s13o=sijo(5);s23o=sijo(6)
  e33o=eijo(3);e13o=eijo(5);e23o=eijo(6)
else
  s33e=ZERO;s13e=ZERO;s23e=ZERO 
  s33eo=ZERO;s13eo=ZERO;s23eo=ZERO
  s33o=ZERO;s13o=ZERO;s23o=ZERO
  e33o=ZERO;e13o=ZERO;e23o=ZERO
end if

call calc_macaulay_bracket_scalar(s33o,s33obp,ONE)
call calc_macaulay_bracket_scalar(e33o,e33obp,ONE)
call calc_macaulay_bracket_scalar(s33o,s33obn,mONE)
call calc_macaulay_bracket_scalar(e33o,e33obn,mONE)

if (s33e >= ZERO .and. Fat < ONE) then
	Fat=(s33e/Rt)**2+a*((s13e/Tl)**2+(s23e/Ti)**2)
	if (Fat >= ONE) then
	  Fat=(s33eo/Rt)**2+a*((s13eo/Tl)**2+(s23eo/Ti)**2)
		const1=ONE/sqrt(Fat)
    
		ueqt0=const1*sqrt(e33obp**2+a*(e13o**2+e23o**2))*catLc
		seqt0=const1**2*(s33obp*e33obp+a*(s13o*e13o+s23o*e23o))/(ueqt0/catLc)
    
		statev(35)=ueqt0
		statev(29)=seqt0
    
		Fat=1.01d0
	end if
	statev(14)=Fat
else if (s33e < ZERO .and. Fac < ONE) then
	Fac=(s33e/Rc)**2
	if (Fac>=ONE) then
	  Fac=(s33eo/Rc)**2
		const1=ONE/sqrt(Fac)
    
		ueqc0=const1*e33obn*catLc
		seqc0=const1*s33obn
    
		statev(36)=ueqc0
		statev(30)=seqc0
    
		Fac=1.01d0
	end if
	statev(15)=Fac
end if

return
end subroutine


!     __________________________________________________________________
!     _______________SUBROUTINE FIBRE DAMAGE EVOLUTION__________________
!     __________________________________________________________________

subroutine calc_fibre_damage(statev,eij,nstatv,ntens,ndi,catLc,Fat,Fac,a,Gt,Gc,ueqt0,ueqc0, & 
                             seqt0,seqc0,dt,dc,dvt,dvc,dmax,dto,dco,dvto,dvco,etat,etac,dtime,ueqt, &
                             ueqtu,ueqc,ueqcu,prest,presc)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer :: i,j
integer, intent(in):: ntens,ndi,nstatv

!     ------------------------REAL VARIABLES---------------------------- 
real (kind=8) :: e11,e12,e13,e11bp,e11bn,ueqtp,ueqcp

real(kind=8), intent(in)  :: Fat,Fac,a,catLc,Gt,Gc,ueqt0,ueqc0,seqt0,seqc0,dmax,&
                             dto,dco,dvto,dvco,etat,etac,dtime,&
                             eij(ntens),prest,presc
real(kind=8), intent(inout) :: statev(nstatv),dt,dc,dvt,dvc,ueqt,ueqtu,ueqc,ueqcu

!     ----------------------DUMMY VARIABLES-----------------------------
real(kind=8) :: const1,const2,const3,const4,const5,const6

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0, mONE=-1.d0

if (ndi == 3) then
  e11=eij(1);e12=eij(4);e13=eij(5)
else
  e11=eij(1);e12=eij(3);e13=ZERO
end if

call calc_macaulay_bracket_scalar(e11,e11bp,ONE)
call calc_macaulay_bracket_scalar(e11,e11bn,mONE)

!     -------------FIBER TENSION DAMAGE EVOLUTION-----------------------
if (Fat >= ONE) then
	ueqt=sqrt(e11bp**2+a*(e12**2+e13**2))*catLc
	ueqtu=TWO*Gt/seqt0
  ueqtp=ueqtu-prest*(ueqtu-ueqt0)
	if (ueqt > ueqt0 .and. ueqt < ueqtp) then
		dt=ueqtu*(ueqt-ueqt0)/(ueqt*(ueqtu-ueqt0))
		dvt=etat/(etat+dtime)*dvto+dtime/(etat+dtime)*dt
  else if (ueqt >= ueqtp) then
    dt=ONE-prest*ueqt0/ueqt
    dvt=etat/(etat+dtime)*dvto+dtime/(etat+dtime)*dt
  end if
  
  dt=min(dt,dmax)
  dvt=min(dvt,dmax)
  dt=max(dt,dto)
  dvt=max(dvt,dvto)
  
  statev(1)=dt
  statev(16)=dvt
  
end if 

!     -------------FIBER COMPRESSION DAMAGE EVOLUTION-------------------
if (Fac >= ONE) then
	ueqc=e11bn*catLc
	ueqcu=TWO*Gc/seqc0
  ueqcp=ueqcu-presc*(ueqcu-ueqc0)
	if (ueqc > ueqc0 .and. ueqc < ueqcp) then
		dc=ueqcu*(ueqc-ueqc0)/(ueqc*(ueqcu-ueqc0))
		dvc=etac/(etac+dtime)*dvco+dtime/(etac+dtime)*dc
  else if (ueqc >= ueqcp) then
    dc=ONE-presc*ueqc0/ueqc
    dvc=etac/(etac+dtime)*dvco+dtime/(etac+dtime)*dc
  end if  
  
  dc=min(dc,dmax)
  dvc=min(dvc,dmax)
  dc=max(dc,dco)
  dvc=max(dvc,dvco)

  statev(2)=dc
	statev(17)=dvc
end if

return
end subroutine


!     __________________________________________________________________
!     _______________SUBROUTINE MATRIX DAMAGE EVOLUTION_________________
!     __________________________________________________________________

subroutine calc_matrix_damage(statev,eij,nstatv,ntens,ndi,catLc,Fat,Fac,a,Gt,Gc,ueqt0,ueqc0, & 
                              seqt0,seqc0,dt,dc,dvt,dvc,dmax,dto,dco,dvto,dvco,etat,etac,dtime,ueqt, &
                              ueqtu,ueqc,ueqcu,prest,presc)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer :: i,j
integer, intent(in):: ntens,ndi,nstatv

!     ------------------------REAL VARIABLES---------------------------- 
real (kind=8) :: e22,e33,e12,e13,e23,e22bp,e33bp,e22bn,e33bn,ueqtp,ueqcp

real(kind=8), intent(in)  :: Fat,Fac,a,catLc,Gt,Gc,ueqt0,ueqc0,seqt0,seqc0,dmax,&
                             dto,dco,dvto,dvco,etat,etac,dtime,&
                             eij(ntens),prest,presc
real(kind=8), intent(inout) :: statev(nstatv),dt,dc,dvt,dvc,ueqt,ueqtu,ueqc,ueqcu

!     ----------------------DUMMY VARIABLES-----------------------------
real(kind=8) :: const1,const2,const3,const4,const5,const6

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0, mONE=-1.d0

if (ndi == 3) then
  e22=eij(2);e33=eij(3);e12=eij(4);e13=eij(5);e23=eij(6)
else
  e22=eij(2);e33=ZERO;e12=eij(3);e13=ZERO;e23=ZERO
end if

call calc_macaulay_bracket_scalar(e22,e22bp,ONE)
call calc_macaulay_bracket_scalar(e33,e33bp,ONE)
call calc_macaulay_bracket_scalar(e22,e22bn,mONE)
call calc_macaulay_bracket_scalar(e33,e33bn,mONE)

!     -------------MATRIX TENSION DAMAGE EVOLUTION----------------------
if (Fat >= ONE) then
	ueqt=sqrt(e22bp**2+e33bp**2+e12**2+e13**2+e23**2)*catLc
	ueqtu=TWO*Gt/seqt0
  ueqtp=ueqtu-prest*(ueqtu-ueqt0)
	if (ueqt > ueqt0 .and. ueqt < ueqtp) then
		dt=ueqtu*(ueqt-ueqt0)/(ueqt*(ueqtu-ueqt0))
		dvt=etat/(etat+dtime)*dvto+dtime/(etat+dtime)*dt
  else if (ueqt >= ueqtp) then
    dt=ONE-prest*ueqt0/ueqt
    dvt=etat/(etat+dtime)*dvto+dtime/(etat+dtime)*dt
  end if
  
  dt=min(dt,dmax)
  dvt=min(dvt,dmax)
  dt=max(dt,dto)
  dvt=max(dvt,dvto)

  statev(3)=dt
	statev(18)=dvt
end if

!     -----------MATRIX COMPRESSION DAMAGE EVOLUTION--------------------
if (Fac >= ONE) then
	ueqc=sqrt(e22bn**2+e33bn**2+e12**2+e13**2+e23**2)*catLc
	ueqcu=TWO*Gc/seqc0
  ueqcp=ueqcu-presc*(ueqcu-ueqc0)
	if (ueqc > ueqc0 .and. ueqc < ueqcp) then
		dc=ueqcu*(ueqc-ueqc0)/(ueqc*(ueqcu-ueqc0))
		dvc=etac/(etac+dtime)*dvco+dtime/(etac+dtime)*dc
  else if (ueqc >= ueqcp) then
    dc=ONE-presc*ueqc0/ueqc
    dvc=etac/(etac+dtime)*dvco+dtime/(etac+dtime)*dc
  end if  
  
  dc=min(dc,dmax)
  dvc=min(dvc,dmax)
  dc=max(dc,dco)
  dvc=max(dvc,dvco)

  statev(4)=dc
	statev(19)=dvc
end if

return
end subroutine


!     __________________________________________________________________
!     ____________SUBROUTINE INTERLAMINAR DAMAGE EVOLUTION______________
!     __________________________________________________________________

subroutine calc_interlaminar_damage(statev,eij,nstatv,ntens,ndi,catLc,Fat,Fac,a,Gt,Gc,ueqt0,ueqc0, & 
                                    seqt0,seqc0,dt,dc,dvt,dvc,dmax,dto,dco,dvto,dvco,etat,etac,dtime,ueqt, &
                                    ueqtu,ueqc,ueqcu,prest,presc)

IMPLICIT NONE

!     -----------------------INTEGER VARIABLES-------------------------- 
integer :: i,j
integer, intent(in):: ntens,ndi,nstatv

!     ------------------------REAL VARIABLES---------------------------- 
real (kind=8) :: e33,e13,e23,e33bp,e33bn,ueqtp,ueqcp

real(kind=8), intent(in)  :: Fat,Fac,a,catLc,Gt,Gc,ueqt0,ueqc0,seqt0,seqc0,dmax,&
                             dto,dco,dvto,dvco,etat,etac,dtime,&
                             eij(ntens),prest,presc
real(kind=8), intent(inout) :: statev(nstatv),dt,dc,dvt,dvc,ueqt,ueqtu,ueqc,ueqcu

!     ----------------------DUMMY VARIABLES-----------------------------
real(kind=8) :: const1,const2,const3,const4,const5,const6

!     ------------DOUBLE PRECISION VALUES OF 0,1,2,3,1/2----------------
real(kind=8), parameter :: ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,&
                           FOUR=4.d0, HALF=0.5d0, mONE=-1.d0

if (ndi == 3) then
  e33=eij(3);e13=eij(5);e23=eij(6)
else
  e33=ZERO;e13=ZERO;e23=ZERO
end if

call calc_macaulay_bracket_scalar(e33,e33bp,ONE)
call calc_macaulay_bracket_scalar(e33,e33bn,mONE)

!     -----------INTERLAMINAR TENSION DAMAGE EVOLUTION------------------
if (Fat >= ONE) then
	ueqt=sqrt(e33bp**2+a*(e13**2+e23**2))*catLc
	ueqtu=TWO*Gt/seqt0
  ueqtp=ueqtu-prest*(ueqtu-ueqt0)
	if (ueqt > ueqt0 .and. ueqt < ueqtp) then
		dt=ueqtu*(ueqt-ueqt0)/(ueqt*(ueqtu-ueqt0))
		dvt=etat/(etat+dtime)*dvto+dtime/(etat+dtime)*dt
  else if (ueqt >= ueqtp) then
    dt=ONE-prest*ueqt0/ueqt
    dvt=etat/(etat+dtime)*dvto+dtime/(etat+dtime)*dt
  end if
  
  dt=min(dt,dmax)
  dvt=min(dvt,dmax)
  dt=max(dt,dto)
  dvt=max(dvt,dvto)

  statev(5)=dt
	statev(20)=dvt
end if

!     ---------INTERLAMINAR COMPRESSION DAMAGE EVOLUTION----------------
if (Fac >= ONE) then
	ueqc=e33bn*catLc
	ueqcu=TWO*Gc/seqc0
  ueqcp=ueqcu-presc*(ueqcu-ueqc0)
	if (ueqc > ueqc0 .and. ueqc < ueqcp) then
		dc=ueqcu*(ueqc-ueqc0)/(ueqc*(ueqcu-ueqc0))
		dvc=etac/(etac+dtime)*dvco+dtime/(etac+dtime)*dc
  else if (ueqc >= ueqcp) then
    dc=ONE-presc*ueqc0/ueqc
    dvc=etac/(etac+dtime)*dvco+dtime/(etac+dtime)*dc
  end if  
  
  dc=min(dc,dmax)
  dvc=min(dvc,dmax)
  dc=max(dc,dco)
  dvc=max(dvc,dvco)
  
  statev(6)=dc
	statev(21)=dvc

end if

return
end subroutine
