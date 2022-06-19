!
! This file is part of LIBPFASST.
!
!>  User defined routines that can be called from inside libpfasst using hooks
module hooks
  use pfasst
  use pf_mod_zndarray
  use pf_mod_zutils
  implicit none
contains

  !>  Output the error and residual in the solution
  subroutine dump_sol(pf, level_index)
    use probin, only: grid_size
    use pf_mod_fftops
    use pf_my_stepper, only: my_stepper_t, as_my_stepper
    type(pf_pfasst_t), intent(inout) :: pf
    character(len = 10) :: yname   !!  output file name for delta_q0
    real(pfdp) :: t    
    integer, intent(in) :: level_index

    class(my_stepper_t), pointer :: stepper
    class(pf_zndarray_t), pointer :: yend
    stepper => as_my_stepper(pf%levels(level_index)%ulevel%stepper)    
    yend => cast_as_zndarray(pf%levels(level_index)%qend)
    !>  compute the exact solution
    t=pf%state%t0+pf%state%dt

    !> save solution at end
    write (yname, "(A1,I0.5,A4)") 'y',pf%state%step+1,'.npy'
    call numpy_dump(stepper%fft_tool,t,yend, (trim(pf%results%datpath) // '/' // trim(yname)))
    print *,'saving ',yname
  end subroutine dump_sol

  !>  Output the Pararael solution for the last processor
  subroutine saveLastProcIteration(pf, level_index)
    use probin, only: save_last_proc_iters
    use pf_mod_fftops
    use pf_my_stepper, only: my_stepper_t, as_my_stepper
    type(pf_pfasst_t), intent(inout) :: pf
    character(len = 24) :: filename
    real(pfdp) :: t    
    integer, intent(in) :: level_index

    integer :: parareal_iteration, proc_rank, num_procs, coarse_step
    class(my_stepper_t), pointer :: stepper
    class(pf_zndarray_t), pointer :: yend

    parareal_iteration = pf%state%iter
    coarse_step = pf%state%step+1
    proc_rank = pf%rank
    num_procs = pf%comm%nproc

    if( (save_last_proc_iters) .and. ( proc_rank .eq. num_procs - 1 ) ) then
        
        stepper => as_my_stepper(pf%levels(level_index)%ulevel%stepper)    
        yend => cast_as_zndarray(pf%levels(level_index)%qend)
    
        write (filename, "(A4,I0.5,A6,I0.5,A4)") 'sol-', coarse_step , '-iter-', parareal_iteration, '.npy'
        call numpy_dump(stepper%fft_tool, t, yend, (trim(pf%results%datpath) // '/' // trim(filename)))

    endif
    
  end subroutine saveLastProcIteration
  
  !>  Output the error and residual in the solution
  subroutine echo_error(pf, level_index)
    use probin, only: grid_size
    use pf_my_stepper, only: my_stepper_t, as_my_stepper
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    class(my_stepper_t), pointer :: stepper
    real(pfdp) :: maxerr,t
    type(pf_zndarray_t) :: y_ex      !<  the initial condition
    
    call zndarray_build(y_ex, pf%levels(level_index)%lev_shape)
    stepper => as_my_stepper(pf%levels(level_index)%ulevel%stepper)    

    
    !>  compute the exact solution
    t=pf%state%t0+pf%state%dt
    call exact(stepper%fft_tool, stepper%fft1d_tools, t, y_ex)

    !>  compute error
    call y_ex%axpy(-1.0d0,pf%levels(level_index)%qend)

    !>  compute error
    maxerr = y_ex%norm()
    print '("error: time: ", f8.4," step: ",i8.1," rank: ",i3.3," iter: ",i4.3," level: ",i2.2," error: ",es14.7," resid: ",es14.7," deltaq0: ",es14.7)', &
         t,pf%state%step+1, pf%rank, pf%state%iter,level_index,maxerr,pf%levels(level_index)%residual,pf%levels(level_index)%max_delta_q0

    call flush(6)

    call zndarray_destroy(y_ex)
    
  end subroutine echo_error

   subroutine set_error(pf, level_index)
    use probin, only: grid_size
    use pf_my_stepper, only: my_stepper_t, as_my_stepper
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    class(my_stepper_t), pointer :: stepper
    real(pfdp) :: maxerr,t
    type(pf_zndarray_t) :: y_ex      !<  the initial condition
    
    call zndarray_build(y_ex, pf%levels(level_index)%lev_shape)
    stepper => as_my_stepper(pf%levels(level_index)%ulevel%stepper)    

    !>  compute the exact solution
    t=pf%state%t0+pf%state%dt
    call exact(stepper%fft_tool, stepper%fft1d_tools, t, y_ex)

    !>  compute error
    call y_ex%axpy(-1.0d0,pf%levels(level_index)%qend)

    !>  compute error
    maxerr = y_ex%norm()

    call zndarray_destroy(y_ex)

    call pf_set_error(pf,level_index,maxerr)
  end subroutine set_error
end module hooks
