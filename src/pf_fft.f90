!!  Module of FFT based routines using fftpack
!
! This file is part of LIBPFASST.
!
!>  Module for using fftpack
module pf_mod_fft_abs
  use pf_mod_dtype
  use pf_mod_utils
  implicit none
  
  real(pfdp), parameter :: two_pi = 6.2831853071795862_pfdp
  !>  Variables and storage for FFT
  type, abstract  :: pf_fft_abs_t
     integer ::   nx,ny,nz                         !! grid sizes
     integer ::   dim                              !! spatial dimension
     real(pfdp) :: Lx, Ly, Lz                      !! domain size
     complex(pfdp), pointer :: wk_1d(:)            ! work space
     complex(pfdp), pointer :: wk_2d(:,:)          ! work space
     complex(pfdp), pointer :: wk_3d(:,:,:)        ! work space                    
   contains
     procedure(pf_fft_s_p),deferred :: fft_setup
     procedure(pf_fft_p),deferred :: fft_destroy
     procedure(pf_fft_p),deferred :: fftf   !  Forward FFT
     procedure(pf_fft_p),deferred :: fftb   !  Inverse (backward)  FFT
     !  FFT
     procedure, private  :: fft_1d, fft_2d, fft_3d, zfft_1d, zfft_2d, zfft_3d  
     generic :: fft => fft_1d, fft_2d, fft_3d, zfft_1d, zfft_2d, zfft_3d  
     !  Inverse FFT
     procedure, private  :: ifft_1d, ifft_2d, ifft_3d,izfft_1d, izfft_2d, izfft_3d  
     generic :: ifft => ifft_1d, ifft_2d, ifft_3d,izfft_1d, izfft_2d, izfft_3d    

     !  Convolution in spectral space     
     procedure, private  :: conv_1d, conv_2d, conv_3d  
     generic :: conv => conv_1d, conv_2d, conv_3d  
     !  Complex convolution in real space     
     procedure, private :: zconv_1d, zconv_2d, zconv_3d  
     generic :: zconv => zconv_1d, zconv_2d, zconv_3d  
     !  Convenience function to grab pointer to workspace
     procedure , private :: get_wk_ptr_1d, get_wk_ptr_2d, get_wk_ptr_3d  
     generic   :: get_wk_ptr =>get_wk_ptr_1d,get_wk_ptr_2d,get_wk_ptr_3d
     !  Construct spectral Laplacian
     procedure , private :: make_lap_1d,make_lap_2d, make_lap_3d
     generic   :: make_lap =>make_lap_1d,make_lap_2d,make_lap_3d
     !  Construct spectral derivative
     procedure , private :: make_deriv_1d,make_deriv_2d, make_deriv_3d
     generic   :: make_deriv =>make_deriv_1d,make_deriv_2d,make_deriv_3d
     !  Restrict in spectral space
     procedure , private :: restrict_1d,restrict_2d, restrict_3d,zrestrict_1d,zrestrict_2d, zrestrict_3d
     generic   :: restrict =>restrict_1d,restrict_2d,restrict_3d,zrestrict_1d,zrestrict_2d,zrestrict_3d

     
  end type pf_fft_abs_t

  interface
     subroutine pf_fft_s_p(this, grid_shape, dim, grid_size)
       import pf_fft_abs_t,pfdp
       class(pf_fft_abs_t), intent(inout) :: this
       integer,              intent(in   ) :: dim
       integer,              intent(in   ) :: grid_shape(dim)
       real(pfdp), optional, intent(in   ) :: grid_size(dim)
     end subroutine pf_fft_s_p
     subroutine pf_fft_p(this)
       import pf_fft_abs_t,pfdp
       class(pf_fft_abs_t), intent(inout) :: this
     end subroutine pf_fft_p
  end interface
  
contains
  ! START private convolution procedures
  
  !++++++++++ Forward FFT real to complex  ++++++++++++++++
  subroutine fft_1d(this, g,ghat)
    class(pf_fft_abs_t), intent(inout) :: this
    real(pfdp), intent(in) :: g(:)
    complex(pfdp), intent(inout) :: ghat(:)

    this%wk_1d=g
    call this%fftf()
    ghat=this%wk_1d
  end subroutine fft_1d

  subroutine fft_2d(this, g,ghat)
    class(pf_fft_abs_t), intent(inout) :: this
    real(pfdp), intent(in) :: g(:,:)
    complex(pfdp), intent(inout) :: ghat(:,:)

    this%wk_2d=g
    call this%fftf()
    ghat=this%wk_2d
  end subroutine fft_2d

    subroutine fft_3d(this, g,ghat)
    class(pf_fft_abs_t), intent(inout) :: this
    real(pfdp), intent(in) :: g(:,:,:)
    complex(pfdp), intent(inout) :: ghat(:,:,:)

    this%wk_3d=g
    call this%fftf()
    ghat=this%wk_3d
  end subroutine fft_3d

  !++++++++++ Backward FFT complex to real   ++++++++++++++++
  subroutine ifft_1d(this,ghat,g)
    class(pf_fft_abs_t), intent(inout) :: this
    real(pfdp), intent(inout) :: g(:)
    complex(pfdp), intent(in) :: ghat(:)

    this%wk_1d=ghat
    call this%fftb()
    g=real(this%wk_1d)
  end subroutine ifft_1d

  subroutine ifft_2d(this, ghat,g)
    class(pf_fft_abs_t), intent(inout) :: this
    real(pfdp), intent(inout) :: g(:,:)
    complex(pfdp), intent(in) :: ghat(:,:)

    this%wk_2d=ghat
    call this%fftb()
    g=real(this%wk_2d)
  end subroutine ifft_2d

    subroutine ifft_3d(this, ghat,g)
    class(pf_fft_abs_t), intent(inout) :: this
    real(pfdp), intent(inout) :: g(:,:,:)
    complex(pfdp), intent(in) :: ghat(:,:,:)

    this%wk_3d=ghat
    call this%fftb()
    g=real(this%wk_3d)
  end subroutine ifft_3d

  !++++++++++ Forward FFT complex to complex   ++++++++++++++++
  subroutine zfft_1d(this, g,ghat)
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp), intent(in) :: g(:)
    complex(pfdp), intent(inout) :: ghat(:)

    this%wk_1d=g
    call this%fftf()
    ghat=this%wk_1d
  end subroutine zfft_1d

    subroutine zfft_2d(this, g,ghat)
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp), intent(in) :: g(:,:)
    complex(pfdp), intent(inout) :: ghat(:,:)

    this%wk_2d=g
    call this%fftf()
    ghat=this%wk_2d
  end subroutine zfft_2d

  subroutine zfft_3d(this, g,ghat)
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp), intent(in) :: g(:,:,:)
    complex(pfdp), intent(inout) :: ghat(:,:,:)

    this%wk_3d=g
    call this%fftf()
    ghat=this%wk_3d
  end subroutine zfft_3d

  !++++++++++ Backward FFT complex to complex   ++++++++++++++++
  subroutine izfft_1d(this,ghat,g)
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp), intent(inout) :: g(:)
    complex(pfdp), intent(in) :: ghat(:)

    
    this%wk_1d=ghat
    call this%fftb()
    g=this%wk_1d
  end subroutine izfft_1d

  ! Take forward FFT
  subroutine izfft_2d(this, ghat,g)
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp), intent(inout) :: g(:,:)
    complex(pfdp), intent(in) :: ghat(:,:)

    this%wk_2d=ghat
    call this%fftb()
    g=this%wk_2d
  end subroutine izfft_2d

  subroutine izfft_3d(this, ghat,g)
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp), intent(inout) :: g(:,:,:)
    complex(pfdp), intent(in) :: ghat(:,:,:)

    this%wk_3d=ghat
    call this%fftb()
    g=this%wk_3d
  end subroutine izfft_3d

  ! Convolve g with spectral op and return in c
  subroutine conv_1d(this, g,op,c)
    class(pf_fft_abs_t), intent(inout) :: this
    real(pfdp), intent(inout) :: g(:)
    complex(pfdp), intent(in) :: op(:)
    real(pfdp), intent(inout) :: c(:)

    this%wk_1d=g
    call this%fftf()
    this%wk_1d = this%wk_1d * op
    call this%fftb()
    c=real(this%wk_1d,pfdp)
  end subroutine conv_1d

  ! Convolve g with spectral op and return in c
  subroutine conv_2d(this, g,op,c)
    class(pf_fft_abs_t), intent(inout) :: this
    real(pfdp), intent(in) :: g(:,:)
    complex(pfdp), intent(in) :: op(:,:)
    real(pfdp), intent(inout) :: c(:,:)
    this%wk_2d=g
    ! Compute Convolution
    call this%fftf()
    this%wk_2d = this%wk_2d * op
    call this%fftb()
    c=real(this%wk_2d,pfdp)        
  end subroutine conv_2d
  
  subroutine conv_3d(this, g,op,c)
        class(pf_fft_abs_t), intent(inout) :: this
    real(pfdp), intent(in) :: g(:,:,:)
    complex(pfdp), intent(in) :: op(:,:,:)
    real(pfdp), intent(inout) :: c(:,:,:)
    this%wk_3d=g
    ! Compute Convolution
    call this%fftf()
    this%wk_3d = this%wk_3d * op
    call this%fftb()
    c=real(this%wk_3d,pfdp)            
  end subroutine conv_3d

 subroutine get_wk_ptr_1d(this,wk) 
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp), pointer,intent(inout) :: wk(:)              ! work space
    wk=>this%wk_1d
  end subroutine get_wk_ptr_1d
  subroutine get_wk_ptr_2d(this,wk) 
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp), pointer,intent(inout) :: wk(:,:)              ! work space
    wk=>this%wk_2d
  end subroutine get_wk_ptr_2d
  subroutine get_wk_ptr_3d(this,wk) 
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp), pointer,intent(inout) :: wk(:,:,:)              ! work space
    wk=>this%wk_3d
  end subroutine get_wk_ptr_3d
    
  subroutine zconv_1d(this, g)
    ! Variable Types
        class(pf_fft_abs_t), intent(inout) :: this
        real(pfdp), intent(in) :: g(:)
        ! Compute Convolution
        call this%fftb()
        this%wk_1d = this%wk_1d * g
        call this%fftf()
    end subroutine zconv_1d

    subroutine zconv_2d(this, g)
        ! Variable Types
        class(pf_fft_abs_t), intent(inout) :: this
        real(pfdp), intent(in) :: g(:,:)
        ! Compute Convolution
        call this%fftb()
        this%wk_2d = this%wk_2d * g
        call this%fftf()
    end subroutine zconv_2d

    subroutine zconv_3d(this, g)
        ! Variable Types
        class(pf_fft_abs_t), intent(inout) :: this
        real(pfdp), intent(in) :: g(:,:,:)
        ! Compute Convolution
        call this%fftb()
        this%wk_3d = this%wk_3d * g
        call this%fftf()
    end subroutine zconv_3d
    
    subroutine make_lap_1d(this, lap)
      class(pf_fft_abs_t), intent(inout) :: this
      complex(pfdp), intent(inout) :: lap(:)
      
      integer     :: i,nx
      real(pfdp)  :: kx, Lx
      
      nx=this%nx
      Lx=this%Lx
      do i = 1, nx
         if (i <= nx/2+1) then
            kx = two_pi / Lx * dble(i-1)
         else
            kx = two_pi / Lx * dble(-nx + i - 1)
         end if
         lap(i) = -kx**2
      end do
      
    end subroutine make_lap_1d
    subroutine make_deriv_1d(this, ddx)
      class(pf_fft_abs_t), intent(inout) :: this
      complex(pfdp), intent(inout) :: ddx(:)
      
      integer     :: i,nx
      real(pfdp)  :: kx, Lx
      
      nx=this%nx
      Lx=this%Lx
      do i = 1, nx
         if (i <= nx/2+1) then
            kx = two_pi / Lx * dble(i-1)
         else
            kx = two_pi / Lx * dble(-nx + i - 1)
         end if
         
         ddx(i) = (0.0_pfdp, 1.0_pfdp) * kx
      end do
    end subroutine make_deriv_1d
     
    subroutine make_lap_2d(this, lap)
      class(pf_fft_abs_t), intent(inout) :: this
      complex(pfdp), intent(inout) :: lap(:,:)
      
      integer     :: i,j,nx,ny
      real(pfdp)  :: kx,ky,Lx,Ly
      
      nx=this%nx
      ny=this%ny
      Lx=this%Lx
      Ly=this%Ly
      
      do j = 1, ny
         if (j <= ny/2+1) then
            ky = two_pi / Ly * dble(j-1)
         else
            ky = two_pi / Ly * dble(-ny + j - 1)
         end if
         do i = 1, nx
            if (i <= nx/2+1) then
               kx = two_pi / Lx * dble(i-1)
            else
               kx = two_pi / Lx * dble(-nx + i - 1)
            end if
            
            lap(i,j) = -(kx**2+ky**2)
         end do
      end do
    end subroutine make_lap_2d

    subroutine make_deriv_2d(this, deriv,dir)
      class(pf_fft_abs_t), intent(inout) :: this
      complex(pfdp), intent(inout) :: deriv(:,:)
      integer, intent(in) :: dir
      
      integer     :: i,j,nx,ny
      real(pfdp)  :: kx,ky,Lx,Ly
      
      nx=this%nx
      ny=this%ny
      Lx=this%Lx
      Ly=this%Ly
      
      do j = 1, ny
         if (j <= ny/2+1) then
            ky = two_pi / Ly * dble(j-1)
         else
            ky = two_pi / Ly * dble(-ny + j - 1)
         end if
         do i = 1, nx
            if (i <= nx/2+1) then
               kx = two_pi / Lx * dble(i-1)
            else
               kx = two_pi / Lx * dble(-nx + i - 1)
            end if
            
            if (dir .eq. 1) then
               deriv(i,j) = (0.0_pfdp,1.0_pfdp)*kx
            else
               deriv(i,j) = (0.0_pfdp,1.0_pfdp)*ky
            endif
         end do
      end do
    end subroutine make_deriv_2d


    subroutine make_lap_3d(this, lap)
      class(pf_fft_abs_t), intent(inout) :: this
      complex(pfdp), intent(inout) :: lap(:,:,:)
      
      integer     :: i,j,k,nx,ny,nz
      real(pfdp)  :: kx,ky,kz,Lx,Ly,Lz
      
      nx=this%nx
      ny=this%ny
      nz=this%nz
      Lx=this%Lx
      Ly=this%Ly
      Lz=this%Lz
      do k = 1,nz
         if (k <= nz/2+1) then
            kz = two_pi / Lz * dble(k-1)
         else
            kz = two_pi / Ly * dble(-nz + k - 1)
         end if
         do j = 1, ny
            if (j <= ny/2+1) then
               ky = two_pi / Ly * dble(j-1)
            else
               ky = two_pi / Ly * dble(-ny + j - 1)
            end if
            do i = 1, nx
               if (i <= nx/2+1) then
                  kx = two_pi / Lx * dble(i-1)
               else
                  kx = two_pi / Lx * dble(-nx + i - 1)
               end if
               lap(i,j,k) = -(kx**2+ky**2+kz**2)
            end do
         end do
      end do
      
    end subroutine make_lap_3d

    subroutine make_deriv_3d(this, deriv,dir)
      class(pf_fft_abs_t), intent(inout) :: this
      complex(pfdp), intent(inout) :: deriv(:,:,:)
      integer, intent(in) :: dir
      
      integer     :: i,j,k,nx,ny,nz
      real(pfdp)  :: kx,ky,kz,Lx,Ly,Lz
      
      nx=this%nx
      ny=this%ny
      nz=this%nz       
      Lx=this%Lx
      Ly=this%Ly
      Lz=this%Lz
      
      select case (dir)
      case (1)  
         do i = 1, nx
            if (i <= nx/2+1) then
               kx = two_pi / Lx * dble(i-1)
            else
               kx = two_pi / Lx * dble(-nx + i - 1)
            end if
            deriv(i,:,:) = (0.0_pfdp,1.0_pfdp)*kx
         end do
      case (2)
         do j = 1, ny
            if (j <= ny/2+1) then
               ky = two_pi / Ly * dble(j-1)
            else
               ky = two_pi / Ly * dble(-ny + j - 1)
            end if
            deriv(:,j,:) = (0.0_pfdp,1.0_pfdp)*ky
         end do
      case (3)
         do k = 1, nz
            if (k <= nz/2+1) then
               kz = two_pi / Lz * dble(k-1)
            else
               kz = two_pi / Ly * dble(-nz + k - 1)
            end if
            deriv(:,:,k) = (0.0_pfdp,1.0_pfdp)*kz
         end do
      case DEFAULT
         call pf_stop(__FILE__,__LINE__,'Bad case in SELECT',dir)
      end select

    end subroutine make_deriv_3d
  !  Interp routines that take a coarse vector and produce a fine    
  
  

  subroutine restrict_1d(this, yvec_f, yvec_c)
    class(pf_fft_abs_t), intent(inout) :: this
    real(pfdp),         pointer :: yvec_f(:), yvec_c(:)
    integer :: nx_f, nx_c,irat
    nx_f = size(yvec_f)
    nx_c = size(yvec_c)
    irat  = nx_f/nx_c

    yvec_c = yvec_f(::irat)
  end subroutine restrict_1d
  subroutine restrict_2d(this, yvec_f, yvec_c)
    class(pf_fft_abs_t), intent(inout) :: this
    real(pfdp),         pointer :: yvec_f(:,:), yvec_c(:,:)
    integer :: nx_f(2), nx_c(2), irat, jrat

    nx_f = shape(yvec_f)
    nx_c = shape(yvec_c)

    irat  = nx_f(1)/nx_c(1)
    jrat  = nx_f(2)/nx_c(2)
    
    yvec_c = yvec_f(::irat,::jrat)           
  end subroutine restrict_2d
  subroutine restrict_3d(this, yvec_f, yvec_c)
    class(pf_fft_abs_t), intent(inout) :: this
    real(pfdp),         pointer :: yvec_f(:,:,:), yvec_c(:,:,:)
    integer :: nx_f(3), nx_c(3)
    integer :: irat, jrat,krat

    nx_f = shape(yvec_f)
    nx_c = shape(yvec_c)

    irat  = nx_f(1)/nx_c(1)
    jrat  = nx_f(2)/nx_c(2)
    krat  = nx_f(3)/nx_c(3)
    
    yvec_c = yvec_f(::irat,::jrat,::krat)           
  end subroutine restrict_3d
  
  subroutine zrestrict_1d(this, yhat_f, yhat_c)
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp),  pointer :: yhat_f(:), yhat_c(:)
    integer :: nx_f, nx_c

    nx_f = size(yhat_f)
    nx_c = size(yhat_c)

    yhat_c(1:nx_c/2) = yhat_f(1:nx_c/2)
    yhat_c(nx_c/2+2:nx_c) = yhat_f(nx_f-nx_c/2+2:nx_f)
    
  end subroutine zrestrict_1d
  subroutine zrestrict_2d(this, yhat_f, yhat_c)
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp),         pointer :: yhat_f(:,:), yhat_c(:,:)
    integer :: nx_f(2), nx_c(2)
  end subroutine zrestrict_2d
  subroutine zrestrict_3d(this, yhat_f, yhat_c)
    class(pf_fft_abs_t), intent(inout) :: this
    complex(pfdp),         pointer :: yhat_f(:,:,:), yhat_c(:,:,:)
    integer :: nx_f(3), nx_c(3)
  end subroutine zrestrict_3d

  end module pf_mod_fft_abs
  
   



