Module pec_module

  Implicit None
  
  Private
  Public::PEC_init,potential,R_e,Morse_energy,Mor_E,D_e
    
  Real*8 D_0(5),R_e(5),ome_e(5),D_e(5),a_0(5)
  Real*8 au_eV,bohr_2_ang,au_cm,mass,shift(5)
  Real*8, allocatable :: Mor_E(:)

  contains
!--------------------------------------------------------------------------
Subroutine PEC_init

! Data are collected from J. Phys. Chem. A 2003, 107, 8521-8529.
!  ^2 Pi_u(O2-) data are collected from J. Chem. Phys. 105, 1807 (1996)
! Array, entries 1-5: 1. ^2 Pi_g (O2-), 2. ^2 Pi_u(O2-) 
! 3. ^3 Sig_g^- (O2) 4. a ^1 Delta_g (O2) 5. b ^1 Sig_g^+ (O2)
  
  D_0=[4.104d0,0.777d0,5.11665d0,4.13935d0,3.48985d0] 
  !eV, D_0(4:5) are shifted using T_0, @Table2
  
  R_e=[1.348d0,1.69d0,1.2075d0, 1.2155d0,1.2268d0]	!Angstrom
  
  ome_e=[1108.d0,592.d0,1580.1d0, 1509.8d0,1432.7d0]	!cm-1
    
  au_eV=27.211386245988d0
  bohr_2_ang=5.291772109038d-1
  au_cm=2.194746313632d5
  
  mass=14582.6d0
  
  D_0=D_0/au_eV
  R_e=R_e/bohr_2_ang
  ome_e=ome_e/au_cm
  
  D_e=ome_e/2.d0+D_0	!Add the zpe
  shift(2)=-D_e(2)+D_e(1)	!3.39 eV is the T_e from Ref.
  a_0=dsqrt(0.5d0*mass/D_e)*ome_e
  
  shift(1)=-0.409/au_eV
  shift(2)= shift(2) + shift(1)
  shift(3)=0.d0
  shift(4)=0.9773d0/au_eV-0.5d0*(ome_e(4)-ome_e(3))
  shift(5)=1.6268/au_eV-0.5d0*(ome_e(5)-ome_e(3))  

End Subroutine PEC_init
!--------------------------------------------------------------------------
Pure function potential(R)

  Integer i
  Real*8 potential(5)  
  Real*8, intent(in) :: R
     
    Do i=1,5
    
      potential(i)=D_e(i)*(1.d0-dexp(-a_0(i)*(R-R_e(i))))**2
      potential(i)=potential(i)+shift(i)
      
    End do
      
End function potential
!------------------------------------------------------------------------------
Subroutine Morse_energy(numE,potindx)
    Integer, intent(in) :: numE,potindx 	!Number of spectrum
    Integer :: n
    Real*8  :: v
    
    If (potindx>5) stop 'Maximum potindx is 5'
    Allocate( Mor_E(numE))
     
    v=a_0(potindx)*dsqrt(2.d0*D_e(potindx)/mass)
    Do n=0,numE-1
        
       Mor_E(n+1)= v*( float(n) +0.5d0) - &
                 & ( v*(float(n) +0.5d0))**2/(4.d0*D_e(potindx))    
       
    End do
    
End Subroutine Morse_energy
!------------------------------------------------------------------------------

End Module pec_module
