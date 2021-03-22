Module schroe_module

  Use bspline_module
  Use basis_module
  Use pec_module
  
  Implicit none
  
  Integer                :: potindx,numB
  Real*8                 :: mu
  Real*8, allocatable    :: KE(:,:),PE(:,:),H(:,:),S(:,:)	! Matrices for SE
  Real*8, allocatable    :: Erg(:),c(:,:) !Eigenvalue and Eigenvector
  Integer, allocatable   :: wfnum(:)    
  
  Contains
  
!-----------------------------------------------------------------------  
Subroutine Calc_full_Hamiltonian(ctype)

  character*5,intent(in) :: ctype
  Integer                :: i,j
    
  If(ctype /= 'bound') stop 'This subroutine is only for bound state'
  
  Write(*,*)'Constructing the full Hamiltonian and the overlap matrix...'  
  
  numB=nB-2	!Number of basis is nB-2, since the 1st and last are excluded
  
  Allocate(KE(numB,numB),PE(numB,numB),H(numB,numB),S(numB,numB))
  
  Do i=2,nB-1
    Do j=i,nB-1
      call Calc_overlap(KE(i-1,j-1),i,j,1,1)
      call Calc_overlap(PE(i-1,j-1),i,j,0,0,pot)
      call Calc_overlap(S(i-1,j-1),i,j,0,0)
    End do    
  End do
  
  KE=0.5d0/mu*KE  
  
  H = KE + PE
! The above calculates the upper trianglar part only.  
! Symmetrize the matrix  

  H = H + transpose(H)
  S = S + transpose(S)
  
! The diagonal elements were doubled by above operations.  
  Do i=1,numB
  
    H(i,i)=H(i,i)/2.d0
    S(i,i)=S(i,i)/2.d0
  
  End do
    
  Write(*,*)'Done.'     
  
End Subroutine Calc_full_Hamiltonian  
!-----------------------------------------------------------------------
Subroutine full_eigensolver

  Integer             :: Info,lwork,i,j
  Real*8, allocatable :: work(:)
  
  Write(*,*)'Solving the generalized eigenvalue problem...'
  
  lwork=3*numB
  Allocate(Erg(numB),work(lwork))

  call dsygv(1,'V','U',numB,H,numB,S,numB,Erg,work,lwork,info)
  Deallocate(work)
  
  Write(*,*)'Done.'
  
End Subroutine full_eigensolver
!-----------------------------------------------------------------------    
Subroutine Calc_banded_Hamiltonian(ctype)

  character*5,intent(in) :: ctype
  Integer                :: i,j
    
  If(ctype /= 'bound') stop 'This subroutine is only for bound state.'
  
  Write(*,*)'Constructing the banded Hamiltonian and the overlap matrix...'  
  
  numB=nB-2	!Number of basis is nB-2, since the 1st and last are excluded
  
  Allocate(KE(order,numB),PE(order,numB),H(order,numB),S(order,numB))
  
  Do i=2,nB-1
    Do j=i,nB-1    
      if ( i>=max(1,j-order)+1 ) then
        call Calc_overlap(KE(order+i-j,j-1),i,j,1,1)
        call Calc_overlap(PE(order+i-j,j-1),i,j,0,0,pot)
        call Calc_overlap(S(order+i-j,j-1),i,j,0,0)        
      end if  
    End do    
  End do
  
  KE=0.5d0/mu*KE  
  
  H = KE + PE
  
  Deallocate(KE,PE)
    
  Write(*,*)'Done.'     
  
End Subroutine Calc_banded_Hamiltonian  
!-----------------------------------------------------------------------
Subroutine banded_eigensolver

  Integer             :: Info,lwork
  Real*8, allocatable :: work(:)
  
  Write(*,*)'Solving the generalized eigenvalue problem for banded matrix...'
  
  lwork=3*numB
  Allocate(Erg(numB),c(numB,numB),work(lwork))

  call dsbgv('V','U',numB,order-1,order-1,H,order,S,order,Erg,c,numB,work,info)

  Deallocate(work)
  
  Write(*,*)'Done.'
  
End Subroutine banded_eigensolver
!----------------------------------------------------------------------- 
! Interface for picking out the potential energy curve needed  
Pure Function pot(x)
  
    Real*8,intent(in) :: x
    Real*8            :: pot,V(5)
    
    V(:)=potential(x)
    pot=V(potindx)
  
End function
!-----------------------------------------------------------------------
Function wavefn(wnum,x)
  Integer              :: inbv,iflag
  Integer,intent(in)   :: wnum
  Real*8, intent(in)   :: x
  Real*8               :: wavefn
  Real*8, allocatable  :: bcoef(:),work(:)
  
  Allocate( bcoef(nB),work(3*order) )
  
  bcoef=0.d0
  bcoef(2:nB-1)=c(1:numB,wnum)
  inbv=1
  
  call dbvalu(R_knots,bcoef,nB,order,0,x,inbv,work,iflag,wavefn)
  
  if (iflag/=0) then
  
    write(*,*) 'iflag=',iflag
    stop 'iflag/=0 in dbvalu!'
    
  end if  
  
  Deallocate(work,bcoef)
  
End Function wavefn
!-----------------------------------------------------------------------
! Calculate the overlap matrix element of \int D^idx[B_i(R)] D^jdx[B_j(R)] (pot(R)) dR
! Potential function is optional. If not present, this subroutine calculate without pot(R).  
Subroutine Calc_wf_overlap(norm,i,j)
  Integer,intent(in)        :: i,j  !The i/jth wave function
  Real*8,intent(out)        :: norm  !Value of the integral
  Integer                   :: iflag
  Real*8, parameter         :: tol = 1.d2*epsilon(1.d0)	!Don't change. Errors occur when tol is too small. 
  Real*8                    :: x1,x2
  Real*8, allocatable       :: bcoef(:),work(:)
  
  Allocate( bcoef(nB),work(3*order) )
 
! Initialize for jth wavefn
  bcoef=0.d0;bcoef(2:nB-1)=c(1:numB,j)
       
  x1=R_knots(order);x2=R_knots(order+Rsteps)	
! dbfqad will determine a smaller interval for integration
 
  call dbfqad(fun,R_knots,bcoef,nB,order,0,x1,x2,tol,norm,iflag,work)
  
  if (iflag/=0) then
  
    write(*,*) 'iflag=',iflag
    stop 'iflag/=0 in dbfqad!'
    
  end if
  
  Deallocate(work,bcoef)
  
    contains
  
  Function fun(x)
    
    Real*8, intent(in) :: x
    Real*8             :: fun
      
      fun=wavefn(i,x)
      
  End function  

End Subroutine Calc_wf_overlap
!-----------------------------------------------------------------------
!Normalize the wave function (in terms of c(1:numB,wfnum))
Subroutine Norm_wavefunction

  Integer :: i
  Real*8  :: norm
  Open(99,file='norm.dat')
  
  Do i=1,numB
    
    call Calc_wf_overlap(norm,i,i)
    write(99,*) i,norm
    c(1:numB,i)=c(1:numB,i)/dsqrt(norm)
    
  End do
  
  Close(99)
  
End Subroutine Norm_wavefunction
!-----------------------------------------------------------------------
Subroutine plot_wavefunction  

  Integer                         :: iR,nR,i
  Real*8                          :: dR,Rtemp
  character*16                    :: filename
  
  dR=1.d-2
  nR=ceiling(R(Rsteps)-R(1))/dR
  
      Write(*,*) 'Printing the wave functions #',wfnum    
      
  Do i=1,size(wfnum)
    
    Write(filename,'("wf_",i0,".dat")') wfnum(i)
    Open(71,file=trim(filename))
    
    Do iR=1, nR
      
        Rtemp= R(1) + float(iR-1) * dR
        write(71,*) Rtemp,wavefn(wfnum(i),Rtemp)
    
    End do

    Close(71) 
  
  End do
  
End Subroutine plot_wavefunction
!-----------------------------------------------------------------------
Subroutine print_pec
  Integer :: iR

  Write(*,*) 'Printing potential energy curve'  
    
  Open(85,file='pec.dat')
  
  Do iR=1,Rsteps
  
    write(85,*) R(iR),pot(R(iR))
    
  End do
  
  Close(85)
  
End Subroutine print_pec
!-----------------------------------------------------------------------
Subroutine print_matrix_elements(ele,dim1,dim2,A)

  character*2,intent(in)          :: ele
  Real*8, dimension(:),intent(in) :: A(:,:)
  character*16                    :: filename
  Integer,intent(in)              :: dim1,dim2
  Integer                         :: i
    
  Write(*,*) 'Printing matrix element for ',trim(ele)   
  
  filename=trim(ele)//".dat"
  Open(85,file=filename)
  
  Do i=1,dim1
  
    Write(85,'(100e20.12)') A(i,1:dim2)
  
  End do
  
  Close(85)
  
End Subroutine print_matrix_elements
!-----------------------------------------------------------------------   

End Module schroe_module
