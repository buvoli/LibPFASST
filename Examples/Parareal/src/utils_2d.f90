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
  subroutine exact(fft, fft1d, t, y_exact)
    use probin, only: dealias
    type(pf_fft_t), pointer, intent(in) :: fft
    type(pf_fft_t), pointer, intent(in) :: fft1d(:)
    real(pfdp), intent(in)  :: t
    type(pf_zndarray_t), intent(inout) :: y_exact
    
    complex(pfdp), pointer :: yex(:,:)
    complex(pfdp), pointer :: yreal(:,:)  !  Real space exact solution

    allocate(yreal(y_exact%arr_shape(1),y_exact%arr_shape(2)))

    call y_exact%get_array(yex)    
    call exact_realspace(t,yreal)

    call fft%fft(yreal,yex)
    if(dealias) then
        call fft%dealias(yex, 2)
    endif

    deallocate(yreal)
    
  end subroutine exact
  
  !> Routine to return the exact solution
  subroutine exact_realspace(t, yex)
    use probin, only: eq_type, dom_size
    real(pfdp), intent(in)  :: t
    complex(pfdp), intent(out) :: yex(:,:)

    yex=cmplx(0.0,0.0,pfdp)

    select case (eq_type)
        case (1) ! KP equation
            call ic_kp(yex, dom_size)
        case DEFAULT
            call pf_stop(__FILE__,__LINE__,'Bad case  for eq_type ',eq_type)
    end select
    
  end subroutine exact_realspace
  
  subroutine ic_kp(u0,dom_size)
    complex(pfdp), intent(inout) :: u0(:,:)
    real(pfdp), intent(in) :: dom_size(2)
    
    integer    :: Nx, Ny, i, j
    real(pfdp) :: x, y, Lx, Ly, x0, m, delta, ky
    
    Nx = SIZE(u0, 1)
    Ny = SIZE(u0, 2)

    Lx = dom_size(1)
    Ly = dom_size(2)
    
    ! parameters
    x0 = -Lx / 4.0_pfdp 
    m = 1.0_pfdp
    delta = 0.2_pfdp; 
    ky = two_pi * m / Ly;
    
    do j = 1, ny
        y = Ly*REAL(j-1,pfdp)/REAL(ny,pfdp)
        do i = 1, nx
            x = Lx*REAL(i-1,pfdp)/REAL(nx,pfdp) - Lx/2.0_pfdp
            u0(i,j) = 2 * ((1.0_pfdp / cosh((x - x0) + delta * cos(ky * y))) ** 2)
        end do
    end do

  end subroutine ic_kp
  
  !> Routine to return set the linear and nonlinear operators
  subroutine set_ops(opL, opR, opNL, ddx, ddy, lap, fft, fft1d)
    use probin, only: eq_type, dealias, rho
    complex(pfdp), intent(inout) :: opL(:,:)
    complex(pfdp), intent(inout) :: opR(:, :)
    complex(pfdp), intent(inout) :: opNL(:,:)
    complex(pfdp), intent(in) :: ddx(:,:),ddy(:,:),lap(:,:)
    type(pf_fft_t), intent(in), pointer :: fft
    type(pf_fft_t), pointer, intent(in) :: fft1d(:)

    select case (eq_type)
        
        case (1) ! KP Equation
            
            opL = 3 * (ddy ** 2) / ddx
            opL(1, :) = 0; ! zero our first x mode
            opL = opL - 1.0_pfdp * (ddx ** 3)
            
            if ( rho .ne. 0.0_pfdp ) then
                opR  = -1.0_pfdp * abs(opL) / tan(two_pi/4.0_pfdp - rho)    
                opL  = opL + opR
            endif

            opNL = - 3 * ddx
            
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
  subroutine f_NL(yvec, fvec, opR, opNL, tmp, fft, fft1d)
    use probin, only: eq_type, dealias, rho
    complex(pfdp), intent(in) :: yvec(:,:)
    complex(pfdp), intent(inout) :: fvec(:,:)
    complex(pfdp), intent(in) :: opR(:,:)
    complex(pfdp), intent(in) :: opNL(:,:)
    complex(pfdp), intent(inout) :: tmp(:,:)
    type(pf_fft_t), intent(in),    pointer :: fft
    type(pf_fft_t), pointer, intent(in) :: fft1d(:)
     
    fvec = yvec
    tmp  = yvec

    select case (eq_type)
        
        case (1) ! KP Equation

            if (dealias) call fft%dealias(tmp,2)
            call fft%ifft(tmp, tmp)
            tmp = tmp * tmp ! u^2
            call fft%fft(tmp, fvec)
            fvec = opNL * fvec ! -3 d/dx u^2
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
     complex(pfdp), allocatable :: lap(:,:) ! Laplacian operators
     complex(pfdp), allocatable :: ddx(:,:) ! first derivative operator
     complex(pfdp), allocatable :: ddy(:,:) ! first derivative operator
     complex(pfdp), allocatable :: opL(:,:) ! implcit operator
     complex(pfdp), allocatable :: opR(:,:) ! Repartitioning operator
     complex(pfdp), allocatable :: opNL(:,:) ! explicit operator
   contains
        procedure :: init  =>  fftops_init
        procedure :: destroy  =>  fftops_destroy
  end type pf_fft_ops_t

  contains

    subroutine fftops_init(this,fft,fft1d,grid_size)
      use probin, only: rho
      class(pf_fft_ops_t), intent(inout)  :: this
      type(pf_fft_t), pointer, intent(in) :: fft
      type(pf_fft_t), pointer, intent(in) :: fft1d(:)
      integer, intent(in) :: grid_size(2)
      
      integer :: istat, nx, ny
      nx = grid_size(1)
      ny = grid_size(2)

      allocate(this%lap(nx,ny),STAT=istat)
      if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
      allocate(this%ddx(nx,ny),STAT=istat)
      if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
      allocate(this%ddy(nx,ny),STAT=istat)
      if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
      allocate(this%opL(nx,ny),STAT=istat)
      if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
      allocate(this%opNL(nx,ny),STAT=istat)
      if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
      allocate(this%opR(nx,ny),STAT=istat)
      if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
      
      call fft%make_deriv(this%ddx,1) !  First derivative
      call fft%make_deriv(this%ddy,2) !  First derivative
      call fft%make_lap(this%lap)  !  Second derivative

      ! initialize  operators
      call set_ops(this%opL,this%opR,this%opNL,this%ddx,this%ddy,this%lap,fft,fft1d)
      
      deallocate(this%lap)
      deallocate(this%ddx)
      deallocate(this%ddy)
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

    complex(pfdp), pointer :: y(:,:)
    complex(pfdp), pointer :: yreal(:,:)  !  Real space exact solution
    logical :: do_complex
    
    do_complex=.false.
    if (present(do_complex_in)) do_complex=do_complex_in

    call y_out%get_array(y)  !  Grab the solution from encapsulationi
    
    if (do_complex) then
       call save_complex_double( fname, shape(y), y)
    else  ! output real solution
       allocate(yreal(y_out%arr_shape(1),y_out%arr_shape(2)))
       call fft%ifft(y,yreal)  !  compute the solution in real space
       call save_complex_double( fname, shape(yreal), yreal)
       deallocate(yreal)
    end if
    
  end subroutine numpy_dump
  
    
end module pf_mod_fftops
