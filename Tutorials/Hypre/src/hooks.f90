!
! This file is part of LIBPFASST.
!
!>  User defined routines that can be called from inside libpfasst using hooks
module hooks
  use encap
  use pf_my_sweeper
  use probin
  implicit none

  interface
     function HypreMaxErr(x, t, init_cond) result(max_err) bind(c, name="HypreMaxErr")
        use iso_c_binding
        type(c_ptr), value :: x
        real(c_double), value :: t
        real(c_double), value :: init_cond
        real(c_double) :: max_err
     end function
  end interface
contains

  !>  Output the error and residual in the solution
  subroutine echo_error(pf, level_index)
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    real(pfdp) :: yexact
    real(pfdp) :: maxerr
    real(pfdp) :: residual
    class(hypre_vector_encap), pointer :: y_end
    integer :: nproc, rank, error

    !> Get the solution at the end of this step
    y_end => cast_as_hypre_vector(pf%levels(level_index)%qend)

    !>  compute error
    maxerr = HypreMaxErr(y_end%c_hypre_vector_ptr, Tfin, init_cond)
    residual=pf%levels(level_index)%residual
   
    call mpi_comm_rank(pf%comm%comm, rank, error)
    call mpi_comm_size(pf%comm%comm, nproc, error)
    

    if (rank == nproc-1 .and. level_index == pf%nlevels .and. pf%state%iter > 0) then
       print '("error: step: ",i3.3," iter: ",i4.3," level: ",i2.2," error: ",es14.7," res: ",es18.10e4)', &
            pf%state%step+1, pf%state%iter,level_index, maxerr,residual
       call flush(6)
    end if
  end subroutine echo_error


end module hooks
