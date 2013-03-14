!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

program fpfasst
  use pfasst
  use feval
  use hooks
  use transfer
  use encap_array1d
  use pf_mod_mpi, only: MPI_COMM_WORLD

  implicit none

  type(pf_pfasst_t)  :: pf
  type(pf_comm_t)    :: comm
  type(pf_sweeper_t), target :: sweeper
  type(pf_encap_t),   target :: encap
  integer            :: ierror, nlevs, nvars(3), nnodes(3), l
  double precision   :: dt

  type(array1d), target :: q0

  !!!! initialize mpi
  call mpi_init(ierror)
  if (ierror .ne. 0) &
       stop "ERROR: Can't initialize MPI."


  !!!! initialize pfasst
  nvars  = [ 32,64,128 ]
  nnodes = [ 3, 5, 9 ]
  dt     = 0.01_pfdp
  nlevs  = 3

  call array1d_encap_create(encap)
  call pf_mpi_create(comm, MPI_COMM_WORLD)
  call pf_imex_create(sweeper, eval_f1, eval_f2, comp_f2)
  call pf_pfasst_create(pf, comm, nlevs)

  pf%niters  = 12
  pf%qtype   = 1

  pf%echo_timings = .false.
  if (nlevs > 1) then
     pf%levels(1)%nsweeps = 2
  end if

  do l = 1, nlevs
     pf%levels(l)%nvars  = nvars(l)
     pf%levels(l)%nnodes = nnodes(l)

     call feval_create_workspace(pf%levels(l)%ctx, nvars(l))

     pf%levels(l)%interpolate => interpolate
     pf%levels(l)%restrict    => restrict
     pf%levels(l)%encap       => encap
     pf%levels(l)%sweeper     => sweeper
  end do

  call pf_mpi_setup(comm, pf)
  call pf_pfasst_setup(pf)

  call add_hook(pf, nlevs, PF_POST_ITERATION, echo_error)

  ! do l = 1, nlevs
  !    call add_hook(pf, l, PF_POST_PREDICTOR, echo_error)
  ! end do


  !!!! initialize advection/diffusion
  allocate(q0%array(nvars(nlevs)))
  call initial(q0)

  !!!! run
  call pf_pfasst_run(pf, c_loc(q0), dt, 0.0_pfdp, 2*comm%nproc)


  !!!! done
  deallocate(q0%array)

  do l = 1, nlevs
     call feval_destroy_workspace(pf%levels(l)%ctx)
  end do

  call pf_imex_destroy(sweeper)
  call pf_pfasst_destroy(pf)
  call pf_mpi_destroy(comm)
  call mpi_finalize(ierror)
  call fftw_cleanup()

end program fpfasst
