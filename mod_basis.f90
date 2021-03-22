Module basis_module
  Use bspline_module
  Implicit none
  
  Integer             :: order,Rsteps,nB
  Real*8, allocatable :: R(:),R_knots(:)
       
  contains
!-----------------------------------------------------------------------
Subroutine RgridMaker(Rmin,Rmax,gtype)
  
  Real*8, intent(in)      :: Rmin, Rmax
  Character*5, intent(in) :: gtype
  Integer                 :: iR
  Real*8                  :: dR  
  
  Allocate(R(Rsteps))
  
  Write(*,*) 'Making a ',gtype,' grid for R, with Rmin,Rmax,Rsteps=' 
  Write(*,'(3x,f5.3,3x,f10.3,3x,I0)')   Rmin,Rmax,Rsteps
  
  Select case(gtype)
  
    Case('const')
    
      dR=(Rmax-Rmin)/float(Rsteps-1)
    
      Do iR=1, Rsteps
         R(iR) = Rmin + float(iR-1) * dR    
      End do
  
    Case default
      Write(*,*)'This grid type is not available.'  
        
  End select
  
End Subroutine RgridMaker
!-----------------------------------------------------------------------
Subroutine RknotMaker
  Integer :: nt
  
  nt=Rsteps + 2*(Order-1)
  Allocate(R_knots(nt))
  R_knots(1:Order-1) = R(1)
  R_knots(nt-Order+1:nt) = R(Rsteps)
  R_knots(Order:Order+Rsteps) = R(1:Rsteps)  
  
End Subroutine RknotMaker
!-----------------------------------------------------------------------
Subroutine printbasisfunction(ideriv,ith)
  Integer, intent(in) ::  ideriv,ith
  Integer             ::  iR,nR
  Real*8              ::  dR,Rtemp
  character*16        ::  filename   
  
  nR=1000
  dR=( R(Rsteps) - R(1) ) / nR  
  
  Write(*,*) 'Printing the ',ideriv,' derivative of ',ith,' Bspline functions'     
    
  Write(filename,'("D",i0,"_spline_",i0,".dat")')ideriv,ith
  Open(71,file=trim(filename))
    
  Do iR=1, nR
      
      Rtemp= R(1) + float(iR-1) * dR
      write(71,*) Rtemp,bspl(ith,ideriv,Rtemp)
    
  End do

  Close(71) 

End Subroutine printbasisfunction
!-----------------------------------------------------------------------
Function bspl(ith,ideriv,x)
  Integer              :: inbv,iflag
  Integer,intent(in)   :: ith,ideriv
  Real*8, intent(in)   :: x
  Real*8               :: bspl  
  Real*8, allocatable  :: bcoef(:),work(:)
  
  Allocate( bcoef(nB),work(3*order) )
  
  bcoef=0.d0
  bcoef(ith)=1.d0
  inbv=1
  
  call dbvalu(R_knots,bcoef,nB,order,ideriv,x,inbv,work,iflag,bspl)
  
  if (iflag/=0) then
  
    write(*,*) 'iflag=',iflag
    stop 'iflag/=0 in dbvalu!'
    
  end if  
  
  Deallocate(work,bcoef)
  
End Function bspl
!-----------------------------------------------------------------------
! Calculate the overlap matrix element of \int D^idx[B_i(R)] D^jdx[B_j(R)] (pot(R)) dR
! Potential function is optional. If not present, this subroutine calculate without pot(R).  
Subroutine Calc_overlap(valu,i,j,idx,jdx,pot)
  Integer,intent(in)        :: i,j,idx,jdx  !The i/jth i/jdx th deriv B-spline
  Real*8,intent(out)        :: valu  !Value of the integral
  Real*8,optional, external :: pot   !Optional, potential function pot(x)
  Integer                   :: iflag
  Real*8, parameter         :: tol = 1.d2*epsilon(1.d0)	!Don't change. Errors occur when tol is too small. 
  Real*8                    :: x1,x2
  Real*8, allocatable       :: bcoef(:),work(:)
  
  Allocate( bcoef(nB),work(3*order) )
 
! Initialize for jth Bspline 
  bcoef=0.d0;bcoef(j)=1.d0     
  x1=R_knots(order);x2=R_knots(order+Rsteps)	
! dbfqad will determine a smaller interval for integration 
  call dbfqad(fun,R_knots,bcoef,nB,order,jdx,x1,x2,tol,valu,iflag,work)
  
  if (iflag/=0) then
  
    write(*,*) 'iflag=',iflag
    stop 'iflag/=0 in dbfqad!'
    
  end if
  
  Deallocate(work,bcoef)
  
    contains
  
  Function fun(x)
    
    Real*8, intent(in) :: x
    Real*8             :: fun
      
    If(present(pot)) then
      
      fun=bspl(i,idx,x)*pot(x)
        
    Else
      
      fun=bspl(i,idx,x)
        
    End if
      
  End function  

End Subroutine Calc_overlap
!-----------------------------------------------------------------------

End module basis_module
