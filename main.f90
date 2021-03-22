! To access the ith basis function, use D^idx u_i(x) = bspl(i,idx,x)
! To inspect basis function, call printbasisfunction(ideriv,ith)
! To inspect potential energy curve, call print_pec
! To inspect matrix elements, & 
! & call print_matrix_elements('H ',numB,H) (char*2 string)
Program Bspline_SE

  Use bspline_module
  Use basis_module
  Use pec_module
  Use schroe_module
  
  Implicit none
  
  Integer                 ::  i,j,iR,num_wf
  Real*8                  ::  Rmin,Rmax,valu,m1,m2,amu,au_eV
  character*5             ::  gtype,ctype,stype
  character*100            ::  filename  
  
! Read input file
  Open(99,file='bse.inp')
  Read(99,*)
  Read(99,*)
  Read(99,*) order
  Read(99,*) 
  Read(99,*)
  Read(99,*) gtype,Rmin,Rmax,Rsteps
  Read(99,*)
  Read(99,*)
  Read(99,*) m1,m2,potindx
  Read(99,*)
  Read(99,*) ctype,stype
  Read(99,*)
  Read(99,*) num_wf
  Allocate(wfnum(num_wf))
  Read(99,*) wfnum   
  Close(99)
  
  Write(*,*)'Order of B-spline =',order
  Write(*,*)'Type of R-grid = ',gtype
  Write(*,*)'mass of 1st and 2nd atoms (in amu)',m1,m2
  Write(*,*)'Index of potential to be used is',potindx 
  Write(*,*)'Type of calculation is ',stype,' and ',ctype
  Write(*,*)'The following wave functions will be printed:',wfnum
   
  amu=1822.888486209d0
  au_eV=27.211386245988d0
  mu=m1*m2/(m1+m2)*amu
  
  Write(*,*)'Reduced mass of the system is (in a.u.)',mu
   
! Create grid points  
  call RgridMaker(Rmin,Rmax,gtype)
  
! Create knot sequence which ensure zero at the left and right boundary  
  call RknotMaker
  
  nB=Rsteps+order-2
  Write(*,*)'Total number of B-Spline functions is',nB
  
! Initiate parameters for potential energy curves  
  call PEC_init
   
  If(trim(stype)=='band') then
   
    call Calc_banded_Hamiltonian(ctype)  
    call banded_eigensolver
    
  Else
    
    call Calc_full_Hamiltonian(ctype)
    call full_eigensolver
    
  End if
  
  write(filename,'("erg_order_",i0,"_Rmax_",i0,"_Rsteps_",i0,".dat")') order,ceiling(Rmax),Rsteps
  
  Open(79,file=trim(filename)) 
    
  Write(*,*) 'Eigenvalues, Morse_energy, and diff are printed at '//trim(filename)//'(in eV)'
  
  call Morse_energy(numB,potindx)

  Do i=1,numB
  
    if (Erg(i) > D_e(potindx)) exit
    Write(79,*) i,Erg(i)*au_eV,Mor_E(i)*au_eV,dabs(Mor_E(i)-Erg(i))*au_eV
    
  End do

  Close(79)
  
!  Call Norm_wavefunction
  
  Call plot_wavefunction
  
  Deallocate(R, R_knots)
  Deallocate(H,S,Erg,Mor_E)
  If(allocated(c)) Deallocate(c)
   
End program Bspline_SE
