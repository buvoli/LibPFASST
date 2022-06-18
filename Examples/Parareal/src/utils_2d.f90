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
    use probin, only: eq_type, dealias
    type(pf_fft_t), pointer, intent(in) :: fft
    type(pf_fft_t), pointer, intent(in) :: fft1d(:)
    real(pfdp), intent(in)  :: t
    type(pf_zndarray_t), intent(inout) :: y_exact
    
    complex(pfdp), pointer :: yex(:,:)
    complex(pfdp), pointer :: yreal(:,:)  !  Real space exact solution
    integer :: i, Nv

    allocate(yreal(y_exact%arr_shape(1),y_exact%arr_shape(2)))

    call y_exact%get_array(yex)    
    call exact_realspace(t,yreal)

    select case (eq_type)
        case (1) ! KP equation (transform x,y into Fourier space)
            call fft%fft(yreal,yex)
            if(dealias) then
                call fft%dealias(yex, 2)
            endif
        case (2) ! VP equation (transform only x into Fourier space)
            Nv = SIZE(yex,2)
            do i = 1, Nv
                call fft1d(1)%fft(yreal(:,i), yex(:,i))
            enddo
    end select

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
        case (2) ! VP equation
            call ic_vp(yex, dom_size)
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
    
  subroutine ic_vp(u0, dom_size)
    complex(pfdp), intent(inout) :: u0(:,:)
    real(pfdp), intent(in) :: dom_size(2)
    
    integer    :: Nx, Nv, i, j
    real(pfdp) :: x, v, Lx, Lv
    
    Nx = SIZE(u0, 1)
    Nv = SIZE(u0, 2)

    Lx = dom_size(1)
    Lv = dom_size(2)
    
    do j = 1, Nv
        v = Lv*REAL(j-1,pfdp)/REAL(Nv,pfdp)  - Lv / 2.0_pfdp
        do i = 1, nx
            x = Lx*REAL(i-1,pfdp)/REAL(Nx,pfdp)
            u0(i,j) = ( 0.9_pfdp / sqrt(two_pi) * exp(-v**2.0_pfdp  / 2.0_pfdp) + 0.2_pfdp / sqrt(two_pi) * exp(-2.0_pfdp * (v - 4.5_pfdp) ** 2.0_pfdp)) * ( 1.0_pfdp + 0.04_pfdp * cos(0.3_pfdp * x))
        end do
    end do

  end subroutine ic_vp

  !> Routine to return set the linear and nonlinear operators
  subroutine set_ops(opL, opR, opNL1, opNL2, ddx, ddy, lap, fft, fft1d)
    use probin, only: eq_type, dealias, rho, dom_size
    complex(pfdp), intent(inout) :: opL(:,:)
    complex(pfdp), intent(inout) :: opR(:, :)
    complex(pfdp), intent(inout) :: opNL1(:,:)
    complex(pfdp), intent(inout) :: opNL2(:,:)
    complex(pfdp), intent(in) :: ddx(:,:),ddy(:,:),lap(:,:)
    type(pf_fft_t), intent(in), pointer :: fft
    type(pf_fft_t), pointer, intent(in) :: fft1d(:)

    integer :: Nv, j
    real(pfdp) :: Lv, v

    select case (eq_type)
        
        case (1) ! KP Equation
            
            opL = 3 * (ddy ** 2) / ddx
            opL(1, :) = 0; ! zero out first x mode
            opL = opL - 1.0_pfdp * (ddx ** 3)
            
            if ( rho .ne. 0.0_pfdp ) then
                opR  = -1.0_pfdp * abs(opL) / tan(two_pi/4.0_pfdp - rho)    
                opL  = opL + opR
            endif

            opNL1 = - 3 * ddx
            
            if(dealias) then
                call fft%dealias(opL, 2)
                call fft%dealias(opR, 2)
                call fft%dealias(opNL1, 2)                
            endif

        case (2) ! VP Equation

            Nv = SIZE(opL, 2)
            Lv = dom_size(2)
    
            opL = ddx
            do j = 1, Nv
                v = Lv*REAL(j-1,pfdp)/REAL(Nv,pfdp)  - Lv / 2.0_pfdp
                opL(:,j) = -1.0_pfdp * v * opL(:,j)               
            end do
        
            if ( rho .ne. 0.0_pfdp ) then
                opR  = -1.0_pfdp * abs(opL) / tan(two_pi/4.0_pfdp - rho)    
                opL  = opL + opR
            endif

            opNl1 = 1 / ddx ! IDX
            opNl1(1, :) = 0 ! zero out first x mode

            opNl2 = ddy ! DV

            if(dealias) then
                do j = 1, Nv
                    call fft1d(1)%dealias(opL(:, j), 2)
                    call fft1d(1)%dealias(opR(:, j), 2)
                    call fft1d(1)%dealias(opNL1(:, j), 2)
                    call fft1d(1)%dealias(opNL2(:, j), 2)
                end do              
            endif

        case DEFAULT
        
            call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',eq_type)

    end select
  
  end subroutine set_ops

  !> Routine to compute the nonlinear operators
  subroutine f_NL(yvec, fvec, opR, opNL1, opNL2, tmp, fft, fft1d)
    use probin, only: eq_type, dealias, rho, dom_size
    complex(pfdp), intent(in) :: yvec(:,:)
    complex(pfdp), intent(inout) :: fvec(:,:)
    complex(pfdp), intent(in) :: opR(:,:)
    complex(pfdp), intent(in) :: opNL1(:,:)
    complex(pfdp), intent(in) :: opNL2(:,:)
    complex(pfdp), intent(inout) :: tmp(:,:)
    type(pf_fft_t), intent(in),    pointer :: fft
    type(pf_fft_t), pointer, intent(in) :: fft1d(:)

    integer :: Nx, Nv, i
    real(pfdp) :: Lv, dv

    select case (eq_type)
        
        case (1) ! KP Equation

            fvec = yvec
            tmp  = yvec
        
            if (dealias) call fft%dealias(tmp,2)
            call fft%ifft(tmp, tmp)
            tmp = tmp * tmp ! u^2
            call fft%fft(tmp, fvec)
            fvec = opNL1 * fvec ! -3 d/dx u^2
            if ( rho .ne. 0.0_pfdp ) then
                fvec = fvec - OpR * yvec
            endif
            if (dealias) call fft%dealias(fvec,2)

        case (2) ! VP Equation

            Nx = SIZE(opNL1, 1)
            Nv = SIZE(opNL1, 2)
            Lv = dom_size(2)
            dv = Lv / Nv

            ! --> E Term  (stored in tmp)     
            tmp = opNL1 ! IDX
            tmp(1, :) = (dv * sum(yvec(1, :)) - Nv) * tmp(1, :)
            do i = 2, Nx
                tmp(i, :) = (dv * sum(yvec(i, :))) * tmp(i, :)
            enddo
            do i = 1, Nv
                call fft1d(1)%ifft(tmp(:, i), tmp(:, i))
            enddo

            !--> F_v term (stored in fvec)
            fvec = yvec
            do i = 1, Nx
                call fft1d(2)%fft(fvec(i, :), fvec(i, :)) ! transform in v
            enddo
            fvec = opNL2 * fvec ! DV * \hat{y}
            call fft%ifft(fvec, fvec)

            ! --> E * F_v
            fvec = -1.0_pfdp * fvec * tmp

            if (dealias) then
                call fft%fft(fvec, fvec)
                call fft%dealias(fvec, 2)
                do i = 1, Nx ! inverse v transform
                    call fft1d(2)%ifft(fvec(i, :), fvec(i, :))
                enddo
            else
                do i = 1, Nv ! forwards x
                    call fft1d(1)%fft(fvec(:,i), fvec(:,i))
                enddo
            endif

            ! Repartitioning
            if ( rho .ne. 0.0_pfdp ) then
                fvec = fvec - OpR * yvec
            endif

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
     complex(pfdp), allocatable :: opNL1(:,:) ! explicit operator
     complex(pfdp), allocatable :: opNL2(:,:) ! explicit operator
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
      allocate(this%opNL1(nx,ny),STAT=istat)
      if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
      allocate(this%opNL2(nx,ny),STAT=istat)
      if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
      allocate(this%opR(nx,ny),STAT=istat)
      if (istat .ne. 0)  call pf_stop(__FILE__,__LINE__,'Allocate failed ',istat)
      
      call fft%make_deriv(this%ddx,1) !  First derivative
      call fft%make_deriv(this%ddy,2) !  First derivative
      call fft%make_lap(this%lap)  !  Second derivative

      ! initialize  operators
      call set_ops(this%opL,this%opR,this%opNL1,this%opNL2, this%ddx,this%ddy,this%lap,fft,fft1d)
      
      deallocate(this%lap)
      deallocate(this%ddx)
      deallocate(this%ddy)
    end subroutine fftops_init

    subroutine fftops_destroy(this)
      class(pf_fft_ops_t), intent(inout)    :: this

      deallocate(this%opL)
      deallocate(this%opNL1)
      deallocate(this%opNL2)
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
    
    do_complex=.true.
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
