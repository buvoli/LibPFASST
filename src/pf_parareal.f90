!!  Routines that run the parareal algorithm
!
! This file is part of LIBPFASST.
!

!> Module of routines to run parareal
module pf_mod_parareal
  use pf_mod_interpolate
  use pf_mod_restrict
  use pf_mod_utils
  use pf_mod_timer
  use pf_mod_dtype
  use pf_mod_hooks
  use pf_mod_pfasst
  use pf_mod_comm
  use pf_mod_results
  
  implicit none
  
contains
  !>  Do the parareal algorithm
  subroutine pf_parareal_run(pf, q0, dt, tend, nsteps, qend, pred_group_size)
    type(pf_pfasst_t), intent(inout), target   :: pf   !!  The complete PFASST structure
    class(pf_encap_t), intent(in   )           :: q0   !!  The initial condition
    real(pfdp),        intent(in   )           :: dt   !!  The time step for each processor
    real(pfdp),        intent(in   )           :: tend !!  The final time of run
    integer,           intent(in   ), optional :: nsteps  !!  The number of time steps
    class(pf_encap_t), intent(inout), optional :: qend    !!  The computed solution at tend
    integer,           intent(in   ), optional :: pred_group_size !! processor group size for Parareal predictor (see pf_parareal_predictor)
    
    !  Local variables
    integer :: nproc  !!  Total number of processors
    integer :: nsteps_loc  !!  local number of time steps
    real(pfdp) :: tend_loc !!  The final time of run
    integer :: ierr, pred_group_size_

    if(present(pred_group_size)) then
        pred_group_size_ = pred_group_size
    else
        pred_group_size_ = 1
    end if

    ! make a local copy of nproc
    nproc = pf%comm%nproc

    !>  Set the number of time steps to do
    !!  The user can either pass in the number of time steps or
    !!  pass in the time step size and length of run
    if (present(nsteps)) then
      nsteps_loc = nsteps
      tend_loc=real(nsteps_loc*dt,pfdp)
    else
      nsteps_loc = ceiling(tend/dt)
      !  Do  sanity check on steps
      if (abs(real(nsteps_loc,pfdp)-tend/dt) > dt/1d-7) then
        print *,'dt=',dt
        print *,'nsteps=',nsteps_loc
        print *,'tend=',tend
       call pf_stop(__FILE__,__LINE__,'Invalid nsteps ,nsteps=',nsteps)
      end if
    end if
    pf%state%nsteps = nsteps_loc


    !  do sanity checks on Nproc
    if (mod(nsteps,nproc) > 0)  call pf_stop(__FILE__,__LINE__,'nsteps must be multiple of nproc ,nsteps=',nsteps)

    ! do sanity checks on group_size
    if(mod(nproc, pred_group_size_) .ne. 0) call pf_stop(__FILE__,__LINE__,'pred_group_size must be a divisor of nproc ,nproc=',nproc)

    !>  Allocate stuff for holding results 
    call initialize_results(pf)
    
    !>  Try to sync everyone
    call mpi_barrier(pf%comm%comm, ierr)

    !> Start timer
    call pf_start_timer(pf, T_TOTAL)
    if (present(qend)) then
       call pf_parareal_block_run(pf, q0, dt, nsteps_loc, pred_group_size_, qend=qend)
    else
       call pf_parareal_block_run(pf, q0, dt,  nsteps_loc, pred_group_size_)
    end if

    !> End timer    
    call pf_stop_timer(pf, T_TOTAL)

    !>  Output stats
    call pf_dump_stats(pf)
    call mpi_barrier(pf%comm%comm, ierr)
  end subroutine pf_parareal_run

  !>  parareal controller for block mode
  subroutine pf_parareal_block_run(pf, q0, dt, nsteps, pred_group_size, qend, flags)
    use pf_mod_mpi, only: MPI_REQUEST_NULL
    type(pf_pfasst_t), intent(inout), target   :: pf
    class(pf_encap_t), intent(in   )           :: q0
    real(pfdp),        intent(in   )           :: dt
    integer,           intent(in   )           :: nsteps
    integer,           intent(in   )           :: pred_group_size
    class(pf_encap_t), intent(inout), optional :: qend
    integer,           intent(in   ), optional :: flags(:)

    class(pf_level_t), pointer :: lev  !!  pointer to the one level we are operating on
    integer                   :: j, k, ierr
    integer                   :: nblocks !!  The number of blocks of steps to do
    integer                   :: nproc   !!  The number of processors being used
    integer                   :: level_index_c !!  Coarsest level in V (Lambda)-cycle
    integer                   :: level_max_depth !!  Finest level in V-cycle

    pf%state%dt      = dt
    pf%state%proc    = pf%rank+1
    pf%state%step    = pf%rank
    pf%state%t0      = pf%state%step * dt

    ! set finest level to visit in the following run
    pf%state%finest_level = pf%nlevels

    !  pointer to finest  level to start
    lev => pf%levels(pf%state%finest_level)

    !  Stick the initial condition into q0 (will happen on all processors)
    call lev%q0%copy(q0, flags=0)


    nproc = pf%comm%nproc
    nblocks = nsteps/nproc

    !  Decide what the coarsest level in the V-cycle is
    level_index_c=1
    if (.not. pf%Vcycle)     level_index_c=pf%state%finest_level

    do k = 1, nblocks   !  Loop over blocks of time steps
       call pf_start_timer(pf, T_BLOCK)
       call call_hooks(pf, -1, PF_PRE_BLOCK)
       ! print *,'Starting  step=',pf%state%step,'  block k=',k
       ! Each block will consist of
       !  1.  predictor
       !  2.  Vcycle until max iterations, or tolerances met
       !  3.  Move solution to next block

       !  Reset some flags
       pf%state%iter    = 0
       pf%state%status  = PF_STATUS_PREDICTOR
       pf%state%pstatus = PF_STATUS_PREDICTOR
       pf%comm%statreq  = MPI_REQUEST_NULL
       pf%comm%sendreq  = MPI_REQUEST_NULL
       pf%state%pfblock = k
       pf%state%sweep = 1   !  Needed for compatibility of residual storage       


       if (k > 1) then
          !>  When starting a new block, broadcast new initial conditions to all procs
          if (pf%debug) print *,'DEBUG-rank=',pf%rank, ' at barrier at k=',k
          call mpi_barrier(pf%comm%comm, ierr)
          if (pf%debug) print *,'DEBUG-rank=',pf%rank, ' past barrier at k=',k
          if (nproc > 1)  then
             call lev%qend%pack(lev%send)    !!  Pack away your last solution
             call pf_broadcast(pf, lev%send, lev%mpibuflen, pf%comm%nproc-1)
             call lev%q0%unpack(lev%send)    !!  Everyone resets their q0
          else
             call lev%q0%copy(lev%qend, flags=0)    !!  Just stick qend in q0
          end if

          !>  Update the step and t_n variables for new block
          pf%state%step = pf%state%step + pf%comm%nproc
          pf%state%t0   = pf%state%step * dt
       end if

       !> Call the predictor to get an initial guess on all levels and all processors
       call pf_parareal_predictor(pf, dt, pred_group_size, flags)
       ! After the predictor, the residual and delta_q0 are just zero
       if (pf%save_delta_q0) call pf_set_delta_q0(pf,1,0.0_pfdp)       
       call pf_set_resid(pf,pf%nlevels,0.0_pfdp)       
       call call_hooks(pf, -1, PF_POST_ITERATION)       !  This is the zero iteration
       
       if (pf%nlevels > 1) then
          !>  Start the parareal iterations
           do j = 1, pf%niters
              call call_hooks(pf, -1, PF_PRE_ITERATION)
             call pf_start_timer(pf, T_ITERATION)
             
             pf%state%iter = j
             
             !  Do a v_cycle
             call pf_parareal_v_cycle(pf, k, dt, 1,2)
             
             !  Check for convergence
             !call pf_check_convergence_block(pf, pf%state%finest_level, send_tag=1111*k+j)
             
              call pf_stop_timer(pf, T_ITERATION)
              call call_hooks(pf, -1, PF_POST_ITERATION)
             
             !  If we are converged, exit block
             if (pf%state%status == PF_STATUS_CONVERGED)  then
                 call call_hooks(pf, -1, PF_POST_CONVERGENCE)
                 call pf_set_iter(pf,j)                 
                 exit
              end if
          end do  !  Loop over j, the iterations in this block
!          if(pf%rank .ne. nproc-1) then
!             print *, pf%rank, nproc-1
!             call pf_start_timer(pf, T_AUX)
!             call pf%comm%wait(pf, 2, ierr)       

!             call pf_stop_timer(pf, T_AUX)
!          endif
       call pf_stop_timer(pf, T_BLOCK)
       call call_hooks(pf, -1, PF_POST_BLOCK)
    end if
    
    
    end do !  Loop over the blocks
    call call_hooks(pf, -1, PF_POST_ALL)

    !  Grab the last solution for return (if wanted)
    if (present(qend)) then
       call qend%copy(lev%qend, flags=0)
    end if
  end subroutine pf_parareal_block_run 

  !>  The parareal predictor computes the initial solution values values at
  !>  each processor.  It sets each of the following values:
  !>
  !>    pf%levels(1)%q0   - coarse solution at beginning of coarse timestep
  !>    pf%levels(1)%qend - coarse solution at end of coarse timestep
  !>    f%levels(2)%q0    - fine solution at beginning of coarse timestep
  !>    pf%levels(2)%qend - fine solution at end of coarse timestep
  !>
  !> This function divides all processors into groups of size group_size and
  !> then uses a single processor to compute the initial values for the entire
  !> group. The value group_size must be a divisor of num_procs. For example,
  !> if groupsize=3 and num_procs=9 then:
  !>
  !>   - Proc 2 computes initial conditions for Procs 0, 1, and 2
  !>   - Proc 5 computes initial conditions for Procs 3, 4, and 5
  !>   - Proc 8 computes initial conditions for Procs 6, 7, and 8
  !>
  !> For problems where the computational cost of a coarse step is 
  !> significantly larger than the computational cost of communication,
  !> the groupsize should be set to the number of mpi tasks per node. This
  !> ensures that (1) any communication remains local to the node and (2) only
  !> one processor per node is used to compute the initial conditions (this
  !> minimizes slowdown caused by running multiple cores on the same node)

  subroutine pf_parareal_predictor(pf, dt, group_size, flags)
    type(pf_pfasst_t), intent(inout), target :: pf     !! PFASST main data structure
    real(pfdp),        intent(in)            :: dt     !! time step
    integer,           intent(in)            :: group_size !! size of processor groups           
    integer,           intent(in), optional  :: flags(:)  !!  User defined flags

    class(pf_level_t), pointer :: c_lev
    class(pf_level_t), pointer :: f_lev     
    integer                   :: k,n             !!  Loop indices
    integer                   :: nsteps_c        !!  Number of RK  steps on coarse level
    integer                   :: level_index     !!  Local variable for looping over levels
    real(pfdp)                :: t_n             !!  Initial time of this processor
    real(pfdp)                :: t_k             !!  Initial time at time step k
    real(pfdp)                :: dt_all          !!  Length of time to integrate 

    integer :: i, factor, leader_rank
    
    call call_hooks(pf, 1, PF_PRE_PREDICTOR)
    call pf_start_timer(pf, T_PREDICTOR)

    !  This is for one two levels only or one if only RK is done
    c_lev => pf%levels(1)
    f_lev => pf%levels(pf%state%finest_level)

    if (pf%debug) print*, 'DEBUG --', pf%rank, 'beginning parareal predictor'
    
    !! Step 1. Getting the initial condition on the coarsest level
    if (pf%state%finest_level > 1) then
       if (pf%q0_style < 2) then  !  Copy coarse
          call c_lev%q0%copy(f_lev%q0)
       end if
    end if

    nsteps_c= c_lev%ulevel%stepper%nsteps
    t_n=pf%state%t0
    
    if(MOD(pf%rank + 1, group_size) .eq. 0) then ! -- group leader --------------------

        ! compute solution for self and all other processors in group
        factor = pf%rank - MOD(pf%rank, group_size)
        
        ! serially advance solution to previous leader (no communication)
        if(factor .ne. 0) then
            call c_lev%ulevel%stepper%do_n_steps(pf, 1, t_n, c_lev%q0, c_lev%qend, dt * real(factor, pfdp), nsteps_c * factor)
        else
            call c_lev%qend%copy(c_lev%q0) ! simply copy info since num_steps = 0
        end if

        ! compute q0 for all previous processors in group
        do i = 1, group_size - 1
            call pf_sendTo(pf, c_lev, 1, pf%rank - group_size + i, 1000, .false.) ! send qend value as q0 for proc with rank - group_size + i
            call c_lev%ulevel%stepper%do_n_steps(pf, 1, t_n + real(i, pfdp) * dt, c_lev%qend, c_lev%qend, dt, nsteps_c) ! single coarse step
        end do
        if (group_size > 1) then
            call pf_sendTo(pf, c_lev, 1, pf%rank - 1, 1001, .false., wait_after=.TRUE.) ! send qend value as qend of processor rank - group_size + i; NOTE: wait_after true since this may be last non-blocking call (avoid warning object was not returned to mpool ucp_requests )
        end if
        
        ! store fine q0 value 
        call f_lev%q0%copy(c_lev%qend)        
        
        ! advance one step, and store fine qend value 
        call c_lev%ulevel%stepper%do_n_steps(pf, 1, t_n + dt * real(group_size, pfdp), c_lev%qend, c_lev%qend, dt, nsteps_c)
        call f_lev%qend%copy(c_lev%qend)
    
    else ! -- follower --------------------------------------------------------
        
        ! receive q0 from node leader
        leader_rank = pf%rank - MOD(pf%rank, group_size) + group_size - 1               
        call pf_recvFrom(pf, c_lev, 0, leader_rank, 1000, .true.) ! recieve q0 from leader
        call f_lev%q0%copy(c_lev%q0)

        ! send q0 to previous neighbor (this provides qend for previous neighbor) - do not send if you are the first processor in the group
        if(MOD(pf%rank, group_size) .ne. 0) then 
            call pf_sendTo(pf, c_lev, 0, pf%rank - 1, 1001, .true.)
        end if

        ! recieve qend from subsequent neighbor
        call pf_recvFrom(pf, c_lev, 1, pf%rank + 1, 1001, .true.)
        call f_lev%qend%copy(c_lev%qend)

    endif

    ! Save the coarse level value to be used in parareal iteration
    call c_lev%Q(1)%copy(f_lev%qend, flags=0)     
   
    call pf_stop_timer(pf, T_PREDICTOR)

    call call_hooks(pf, -1, PF_POST_PREDICTOR)

    pf%state%status = PF_STATUS_ITERATING
    pf%state%pstatus = PF_STATUS_ITERATING
    if (pf%debug) print*,  'DEBUG --', pf%rank, 'ending predictor'

  end subroutine pf_parareal_predictor

  !> Execute a parareal V-cycle (iteration)
  !!  It is assumed that we have two levels and two nodes here
  !!  When this is called the previous coarse integrator result should be stored in Q(1)
  !!  and the parareal iteration in qend (both on coarse level).  If this is called
  !!  directly after the predictor, these will be the same thing
  subroutine pf_parareal_v_cycle(pf, iteration, dt,level_index_c,level_index_f, flags)


    type(pf_pfasst_t), intent(inout), target :: pf
    real(pfdp),        intent(in)    :: dt
    integer,           intent(in)    :: iteration
    integer,           intent(in)    :: level_index_c  !! Coarsest level of V-cycle (not supported)
    integer,           intent(in)    :: level_index_f  !! Finest level of V-cycle (not supported)
    integer, optional, intent(in)    :: flags

    type(pf_level_t), pointer :: f_lev, c_lev
    integer :: level_index, j,nsteps_f,nsteps_c

    if (pf%nlevels <2) return     !  This is for two levels only

    c_lev => pf%levels(1)
    f_lev => pf%levels(2)
    nsteps_c= c_lev%ulevel%stepper%nsteps
    nsteps_f= f_lev%ulevel%stepper%nsteps  

    !  Save the old value of q0 and qend so that we can compute difference

     if (pf%save_delta_q0) call c_lev%delta_q0%copy(f_lev%q0, flags=0) !  Prime the delta_q0 stored in c_lev%delta_q0

     call f_lev%delta_q0%copy(f_lev%qend, flags=0) !  Holding delta_qend in f_lev%delta_q0

    !  Step on fine and store in  fine qend 
    level_index=2
    call f_lev%ulevel%stepper%do_n_steps(pf, level_index,pf%state%t0, f_lev%q0,f_lev%qend, dt, nsteps_f)

    !  Subtract the old coarse to get parareal correction  in f_lev%qend
    call f_lev%qend%axpy(-1.0_pfdp,c_lev%Q(1))
    
    ! Get a new initial condition on fine (will be put in q0)
    call pf_recv(pf, f_lev, 10000+iteration, .true.)

    !  Step on coarse and save in Q(1) for next iteration
    level_index=1    
    call c_lev%ulevel%stepper%do_n_steps(pf, level_index,pf%state%t0, f_lev%q0,c_lev%Q(1), dt, nsteps_c)
!    call c_lev%qend%copy(c_lev%Q(1), flags=0) !  save in Q(1) for next iteration
    !  Finish the parareal update (store in fine qend) F_old-G_old+G_new
    call f_lev%qend%axpy(1.0_pfdp,c_lev%Q(1))        

    !  Send new solution  forward  (nonblocking)
!    call pf_send(pf, f_lev, 10000+iteration, .false.)
    call pf_send(pf, f_lev, 10000+iteration, .true.)


    !  Complete the delta_q0 on coarse with new initial condition
    call c_lev%delta_q0%axpy(-1.0_pfdp,f_lev%q0, flags=0) !  Complete delta_q0

    !  Complete the jump at the end
     if (pf%save_delta_q0) call f_lev%delta_q0%axpy(-1.0_pfdp,f_lev%qend)

    !  Save residual
    if (pf%save_residuals) then    
        f_lev%residual=f_lev%delta_q0%norm(flags=0)     ! max jump in qend
        call pf_set_resid(pf,1,f_lev%residual)
        call pf_set_resid(pf,2,f_lev%residual)
     end if
    !  Save jumps
     if (pf%save_delta_q0) then
        c_lev%max_delta_q0=c_lev%delta_q0%norm(flags=0) ! max jump in q0
        call pf_set_delta_q0(pf,1,c_lev%max_delta_q0)
        call pf_set_delta_q0(pf,2,c_lev%max_delta_q0)
     end if

  end subroutine pf_parareal_v_cycle
  

  !> Subroutine to check if the current processor has converged and
  !> to update the next processor on the status
  !> Note that if the previous processor hasn't converged yet
  !> (pstatus), the current processor can't be converged yet either
  subroutine pf_check_convergence_block(pf, level_index, send_tag)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: level_index
    integer,           intent(in)    :: send_tag  !! identifier for status send and receive

    logical           :: residual_converged, converged


    ! Shortcut for fixed iteration mode
    if (pf%abs_res_tol == 0.0 .and. pf%rel_res_tol == 0.0) then
       pf%state%pstatus = PF_STATUS_ITERATING
       pf%state%status  = PF_STATUS_ITERATING
       return
    end if

    call call_hooks(pf, 1, PF_PRE_CONVERGENCE)

    !> Check to see if tolerances are met
    call pf_check_residual(pf, level_index, residual_converged)

    !>  Until I hear the previous processor is done, recieve it's status
    if (pf%state%pstatus /= PF_STATUS_CONVERGED) call pf_recv_status(pf, send_tag)

    !>  Check to see if I am converged
    converged = .false.
    if (residual_converged) then
       if (pf%rank == 0) then
          converged = .true.
       else  !  I am not the first processor, so I need to check the previous one
          if (pf%state%pstatus == PF_STATUS_CONVERGED) converged = .true.
       end if
    end if ! (residual_converged)

    !>  For parareal, Proc N is converged after iteration N
    if (pf%rank .lt. pf%state%iter) then
       converged = .true.
    end if
    
    !> Assign status and send it forward
    if (converged) then
       if (pf%state%status == PF_STATUS_ITERATING) then
          !  If I am converged for the first time
          !  then flip my flag and send the last status update
          pf%state%status = PF_STATUS_CONVERGED
          call pf_send_status(pf, send_tag)
       end if
    else
       !  I am not converged, send the news
       pf%state%status = PF_STATUS_ITERATING
       call pf_send_status(pf, send_tag)
    end if

  end subroutine pf_check_convergence_block

  !> Subroutine to test residuals to determine if the current processor has converged.
  subroutine pf_check_residual(pf, level_index, residual_converged)
    type(pf_pfasst_t), intent(inout) :: pf
    integer,           intent(in)    :: level_index
    logical,           intent(out)   :: residual_converged  !! Return true if residual is below tolerances

    residual_converged = .false.

    ! Check to see if absolute tolerance is met
    if   (pf%levels(level_index)%residual     < pf%abs_res_tol)  then
       if (pf%debug) print*, 'DEBUG --',pf%rank, 'residual tol met',pf%levels(level_index)%residual
       residual_converged = .true.
    end if

  end subroutine pf_check_residual

end module pf_mod_parareal
