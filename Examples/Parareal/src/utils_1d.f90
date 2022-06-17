!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>  Some helpful routines which are problem/DIM dependent  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module pf_mod_zutils
  use pfasst
  use pf_mod_solutions
  use pf_mod_zndarray
  use pf_mod_fftpackage
  use fnpy
  implicit none
contains  

  !> Routine to return the exact solution
  subroutine exact(fft,t, y_exact)  
    type(pf_fft_t), pointer, intent(in) :: fft
    real(pfdp), intent(in)  :: t
    type(pf_zndarray_t), intent(inout) :: y_exact
    
    complex(pfdp), pointer :: yex(:)
    complex(pfdp), pointer :: yreal(:)  !  Real space exact solution

    allocate(yreal(y_exact%arr_shape(1)))
    call y_exact%get_array(yex)    
    call exact_realspace(t,yreal)
    call fft%fft(yreal,yex)

    deallocate(yreal)
    
  end subroutine exact
  
  !> Routine to return the exact solution
  subroutine exact_realspace(t, yex)
    use probin, only: eq_type,lam1,lam2,nu, a, t00, kfreq,dom_size,beta,ic_type
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(out) :: yex(:)

    yex = cmplx(0.0,0.0,pfdp)

    select case (eq_type)
    case (0) 
       yex=exp(zi*(lam1+lam2)*t)
    case (1) ! Nonlinear Schrodinger Equation
        select case (ic_type)
            case (1)
                call ic_nls_ppw_smooth(yex, dom_size(1))
            case (2)
                call ic_nls_ppw_oscillatory(yex, dom_size(1))
            case (3)  
                call ic_nls_ppw_fullspectrum(yex, dom_size(1))
        end select
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',eq_type)
    end select
    
  end subroutine exact_realspace
  
  ! perturbed plane wave initial condition for NLS 

  subroutine ic_nls_ppw_smooth(u0, Lx)
  
    complex(pfdp),  intent(inout) :: u0(:)
    real(pfdp),     intent(in) :: Lx  
  
    real(pfdp), dimension(1) :: amps
    integer,    dimension(1) :: modes

    modes(1) = 1
    amps(1)  = 0.01_pfdp

    call ic_nls_ppw(u0, Lx, modes, amps)
  
  end subroutine ic_nls_ppw_smooth

  subroutine ic_nls_ppw_oscillatory(u0, Lx)
  
    complex(pfdp),  intent(inout) :: u0(:)
    real(pfdp),     intent(in) :: Lx

    real(pfdp) :: amps(2)  = [0.01_pfdp, 0.01_pfdp]
    integer    :: modes(2) = [1, 45]

    call ic_nls_ppw(u0, Lx, modes, amps)
  
  end subroutine ic_nls_ppw_oscillatory

  subroutine ic_nls_ppw_fullspectrum(u0, Lx)
  
    complex(pfdp),  intent(inout) :: u0(:)
    real(pfdp),     intent(in) :: Lx  
  
    integer :: i
    real(pfdp) :: amps(45)  = 0.01_pfdp
    integer    :: modes(45) = [(i, i=1, 45, 1)]

    call ic_nls_ppw(u0, Lx, modes, amps)
  
  end subroutine ic_nls_ppw_fullspectrum
  
  subroutine ic_nls_ppw(u0, Lx, modes, amps)
    
    complex(pfdp),  intent(inout) :: u0(:)
    real(pfdp),     intent(in) :: Lx
    integer,        intent(in) :: modes(:)
    real(pfdp),     intent(in) :: amps(:)
    
    integer    :: Nx, Nm, i, j
    real(pfdp) :: x, amp

    Nx  = SIZE(u0)
    Nm  = SIZE(modes)
    
    do i = 1, Nx
       x = Lx * REAL(i-1,pfdp) / REAL(nx,pfdp) - Lx/2.0_pfdp
       u0(i) = 1.0_pfdp
       do j = 1, Nm
        u0(i) = u0(i) + amps(j) * cos(two_pi * x  * modes(j) / Lx)
       enddo
    end do

  end subroutine ic_nls_ppw

  !> Routine to return set the linear and nonlinear operators
  subroutine set_ops(opL, opR, opNL, ddx,lap, fft)
    use probin, only: eq_type, dealias, rho
    complex(pfdp), intent(inout) :: opL(:)
    complex(pfdp), intent(inout) :: opNL(:)
    complex(pfdp), intent(inout) :: opR(:)
    complex(pfdp), intent(in) :: ddx(:),lap(:)
    type(pf_fft_t), intent(in),  pointer :: fft
    
    select case (eq_type)
        case (1)  ! Nonlinear Schrodinger Equation
            opL  = zi * lap
            if ( rho .ne. 0.0_pfdp ) then
                opR  = -1.0_pfdp * abs(opL) / tan(two_pi/4.0_pfdp - rho)    
                opL  = opL + opR
            endif
            opNL = 0.0_pfdp
            
            if(dealias) then
                call fft%dealias(opL, 2)
                call fft%dealias(opR, 2)
                call fft%dealias(opNL, 2)                
            endif
        case DEFAULT
            call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',eq_type)
    end select

  end subroutine set_ops

  !> Routine to compute the nonlinear operators
  subroutine f_NL(yvec,fvec,opR,opNL,tmp,fft)
    use probin, only: eq_type,gamma,beta, dealias, rho
    complex(pfdp), intent(in) :: yvec(:)
    complex(pfdp), intent(inout) :: fvec(:)
    complex(pfdp), intent(in) :: opNL(:)
    complex(pfdp), intent(in) :: opR(:)
    complex(pfdp), intent(inout) :: tmp(:)
    type(pf_fft_t), intent(in),    pointer :: fft
 
    fvec=yvec
    tmp=yvec
    select case (eq_type)
    case (1) ! Nonlinear Schrodinger Equation
        if (dealias) call fft%dealias(tmp,2)
        call fft%ifft(tmp,tmp)
        fvec=conjg(tmp)*tmp*tmp
        call fft%fft(fvec,fvec)
        fvec=2.0_pfdp*zi*fvec
        if ( rho .ne. 0.0_pfdp ) then
             fvec = fvec - OpR * yvec
        endif
        if (dealias) call fft%dealias(fvec,2)
    case DEFAULT
       call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',eq_type)
    end select

  end subroutine f_NL
  
end module pf_mod_zutils

!  Module to hold the spectral operators
module pf_mod_fftops
  use pf_mod_dtype
  use pf_mod_fftpackage  
  use pf_mod_zutils
  implicit none

  type, public  :: pf_fft_ops_t
     complex(pfdp), allocatable :: lap(:) ! Laplacian operators
     complex(pfdp), allocatable :: ddx(:)  ! first derivative operator
     complex(pfdp), allocatable :: opL(:)  ! implcit operator
     complex(pfdp), allocatable :: opNL(:) ! explicit operator
     complex(pfdp), allocatable :: opR(:)  ! Repartitioning operator
   contains
        procedure :: init  =>  fftops_init
        procedure :: destroy  =>  fftops_destroy
  end type pf_fft_ops_t

  contains

    subroutine fftops_init(this,fft,grid_size)
      use probin, only: d0,d1,r0,r1,rho
      class(pf_fft_ops_t), intent(inout)    :: this
      type(pf_fft_t), pointer, intent(in) :: fft
      integer, intent(in) :: grid_size(1)

      integer :: istat, nx
      nx = grid_size(1)

      allocate(this%lap(nx),STAT=istat)
      if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
      allocate(this%ddx(nx),STAT=istat)
      if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
      allocate(this%opL(nx),STAT=istat)
      if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
      allocate(this%opNL(nx),STAT=istat)
      if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
      
      call fft%make_deriv(this%ddx) !  First derivative
      call fft%make_lap(this%lap)   !  Second derivative
      allocate(this%opR(nx),STAT=istat)

      ! initialize  operators
      call set_ops(this%opL, this%opR, this%opNL, this%ddx, this%lap, fft)

      deallocate(this%lap)
      deallocate(this%ddx)
    end subroutine fftops_init

    subroutine fftops_destroy(this)
      class(pf_fft_ops_t), intent(inout)    :: this

      deallocate(this%opL)
      deallocate(this%opNL)
      deallocate(this%opR)
      
    end subroutine fftops_destroy
    
  !> Routine to return out put the solution to numpy (dimension dependent)
  subroutine numpy_dump(fft,t, y_out,fname, do_complex_in)
    type(pf_fft_t), pointer, intent(in) :: fft
    real(pfdp), intent(in)  :: t
    type(pf_zndarray_t), intent(inout) :: y_out
    character(len=*),  intent(in   ) :: fname
    logical,           intent(in   ), optional :: do_complex_in

    complex(pfdp), pointer :: y(:)
    complex(pfdp), pointer :: yreal(:)  !  Real space exact solution
    logical :: do_complex
    
    do_complex=.false.
    if (present(do_complex_in)) do_complex=do_complex_in

    call y_out%get_array(y)  !  Grab the solution from encapsulationi
    
    if (do_complex) then
       call save_complex_double( fname, shape(y), y)
    else  ! output real solution
       allocate(yreal(y_out%arr_shape(1)))
       call fft%ifft(y,yreal)  !  compute the solution in real space
       call save_complex_double( fname, shape(yreal), yreal)
       deallocate(yreal)
    end if
    
  end subroutine numpy_dump
  
    
end module pf_mod_fftops
