!
! Copyright (C) 2012, 2013 Matthew Emmett and Michael Minion.
!
! This file is part of LIBPFASST.
!
! LIBPFASST is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! LIBPFASST is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with LIBPFASST.  If not, see <http://www.gnu.org/licenses/>.
!

module pf_mod_imex
  use pf_mod_dtype
  use pf_mod_utils
  implicit none

  type, extends(pf_sweeper_t), abstract :: pf_imex_t
     real(pfdp), allocatable :: SdiffE(:,:)
     real(pfdp), allocatable :: SdiffI(:,:)
   contains
     procedure :: sweep        => imex_sweep
     procedure :: initialize   => imex_initialize
     procedure :: evaluate     => imex_evaluate
     procedure :: integrate    => imex_integrate
     procedure :: residual     => imex_residual
     procedure :: evaluate_all => imex_evaluate_all
     procedure :: destroy      => imex_destroy
     procedure :: imex_destroy
     procedure :: f3eval => f3eval_dummy
     procedure :: f3comp => f3comp_dummy
  end type pf_imex_t

  interface

  end interface

contains

  ! Perform on SDC sweep on level lev and set qend appropriately.
  subroutine imex_sweep(this, pf, lev, t0, dt)
    use pf_mod_timer

    class(pf_imex_t),  intent(inout) :: this
    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in   ) :: dt, t0
    class(pf_level_t), intent(inout) :: lev

    integer     :: m, n
    real(pfdp)  :: t
    real(pfdp)  :: dtsdc(1:lev%nnodes-1)
    class(pf_encap_t), allocatable :: rhs
    call start_timer(pf, TLEVEL+lev%level-1)

    ! compute integrals and add fas correction
    do m = 1, lev%nnodes-1
       call lev%S(m)%setval(0.0_pfdp)
       do n = 1, lev%nnodes
          call lev%S(m)%axpy(dt*this%SdiffE(m,n), lev%F(n,1))
          call lev%S(m)%axpy(dt*this%SdiffI(m,n), lev%F(n,2))
       end do
       if (allocated(lev%tau)) then
          call lev%S(m)%axpy(1.0_pfdp, lev%tau(m))
       end if
    end do

    ! do the time-stepping
    call lev%Q(1)%unpack(lev%q0)

    !    call this%f1eval(lev%Q(1), t0, lev%level, lev%F(1,1))
    call this%f1eval(lev%Q(1), t0, lev%level, lev%F(1,1))    
    call this%f2eval(lev%Q(1), t0, lev%level, lev%F(1,2))

    call lev%ulevel%factory%create_single(rhs, lev%level, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)

    t = t0
    dtsdc = dt * (lev%nodes(2:lev%nnodes) - lev%nodes(1:lev%nnodes-1))
    do m = 1, lev%nnodes-1
       t = t + dtsdc(m)

       call rhs%copy(lev%Q(m))
       call rhs%axpy(dtsdc(m), lev%F(m,1))
       call rhs%axpy(1.0_pfdp, lev%S(m))

       call this%f2comp(lev%Q(m+1), t, dtsdc(m), rhs, lev%level, lev%F(m+1,2))
       call this%f1eval(lev%Q(m+1), t, lev%level, lev%F(m+1,1))
    end do

    call lev%qend%copy(lev%Q(lev%nnodes))

    ! done
    call lev%ulevel%factory%destroy_single(rhs, lev%level, SDC_KIND_SOL_FEVAL, lev%nvars, lev%shape)
   
    call end_timer(pf, TLEVEL+lev%level-1)
  end subroutine imex_sweep

  ! Evaluate function values
  subroutine imex_evaluate(this, lev, t, m)
    class(pf_imex_t),  intent(inout) :: this
    real(pfdp),        intent(in   ) :: t
    integer,           intent(in   ) :: m
    class(pf_level_t), intent(inout) :: lev

    call this%f1eval(lev%Q(m), t, lev%level, lev%F(m,1))
    call this%f2eval(lev%Q(m), t, lev%level, lev%F(m,2))
  end subroutine imex_evaluate

  ! Initialize matrices
  subroutine imex_initialize(this, lev)
    class(pf_imex_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev

    real(pfdp) :: dsdc(lev%nnodes-1)
    integer    :: m, nnodes

    this%npieces = 2

    nnodes = lev%nnodes
    allocate(this%SdiffE(nnodes-1,nnodes))  !  S-FE
    allocate(this%SdiffI(nnodes-1,nnodes))  !  S-BE

    this%SdiffE = lev%s0mat
    this%SdiffI = lev%s0mat

    dsdc = lev%nodes(2:nnodes) - lev%nodes(1:nnodes-1)
    do m = 1, nnodes-1
       this%SdiffE(m,m)   = this%SdiffE(m,m)   - dsdc(m)
       this%SdiffI(m,m+1) = this%SdiffI(m,m+1) - dsdc(m)
    end do
  end subroutine imex_initialize

  ! Destroy the matrices
  subroutine imex_destroy(this, lev)
    class(pf_imex_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    
    deallocate(this%SdiffE)
    deallocate(this%SdiffI)
  end subroutine imex_destroy

  ! Compute SDC integral
  subroutine imex_integrate(this, lev, qSDC, fSDC, dt, fintSDC)
    class(pf_imex_t),  intent(inout) :: this
    class(pf_level_t), intent(in   ) :: lev
    class(pf_encap_t), intent(in   ) :: qSDC(:), fSDC(:, :)
    real(pfdp),        intent(in   ) :: dt
    class(pf_encap_t), intent(inout) :: fintSDC(:)

    integer :: n, m, p

    do n = 1, lev%nnodes-1
       call fintSDC(n)%setval(0.0_pfdp)
       do m = 1, lev%nnodes
          do p = 1, this%npieces
             call fintSDC(n)%axpy(dt*lev%s0mat(n,m), fSDC(m,p))
          end do
       end do
    end do
  end subroutine imex_integrate

  subroutine imex_residual(this, lev, dt)
    class(pf_imex_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    real(pfdp),        intent(in   ) :: dt
    call pf_generic_residual(this, lev, dt)
  end subroutine imex_residual

  subroutine imex_evaluate_all(this, lev, t)
    class(pf_imex_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    real(pfdp),        intent(in   ) :: t(:)
    call pf_generic_evaluate_all(this, lev, t)
  end subroutine imex_evaluate_all

  subroutine f3eval_dummy(this, y,t, level, f)
    !  Dummy routine
    class(pf_imex_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: y
    class(pf_encap_t),   intent(inout) :: f
    real(pfdp),          intent(in   ) :: t
    integer,             intent(in   ) :: level
  end subroutine f3eval_dummy
  subroutine f3comp_dummy(this, y,t,dt,rhs, level, f)
    !  Dummy routine
    class(pf_imex_t), intent(inout) :: this
    class(pf_encap_t),   intent(inout) :: y
    real(pfdp),          intent(in   ) :: t
    real(pfdp),          intent(in   ) :: dt
    class(pf_encap_t),   intent(in   ) :: rhs
    integer,             intent(in   ) :: level
    class(pf_encap_t),   intent(inout) :: f
  end subroutine f3comp_dummy

end module pf_mod_imex
