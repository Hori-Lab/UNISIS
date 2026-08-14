program unisis

   use, intrinsic :: iso_fortran_env, only : iostat_end, compiler_version, compiler_options, output_unit
   use const, only : CHAR_FILE_PATH, init_const, PREC, FILENAME_DIGIT_REPLICA
   use const_phys, only : BOLTZ_KCAL_MOL
   use const_idx, only : ENE, JOBT, REPT, RSTBLK
   use var_potential, only : flg_ele, flg_twz, flg_bias_rg, flg_restraint
   use var_state, only : restarted, xyz, tempK, kT, job, opt_anneal, anneal_tempK, reset_step
   use var_top, only : nmp, flg_freeze
   use var_io, only : flg_out_bp, flg_out_bpall, flg_out_bpe, hdl_in_rst, &
                      hdl_out, hdl_bp, hdl_bpall, hdl_bpe
   use var_parallel, only : init_parallel, end_parallel
   use var_replica, only : nrep_proc, flg_replica
   use pbc, only : flg_pbc, init_pbc
   use unisis_version, only : UNISIS_RELEASE_LABEL
   !use mt19937_64, only : init_genrand64

   implicit none

   character(len=CHAR_FILE_PATH) :: cfile_inp, cfile_rst
   character(len=CHAR_FILE_PATH) :: carg

   integer :: nargs
   integer :: istat
   integer :: irep
   integer :: rst_status
   logical :: stat

   call init_const()

   nargs = command_argument_count()
   if (nargs == 1) then
      call get_command_argument(1, carg)
      if (trim(carg) == '--version') then
         call print_version()
         stop
      endif
   endif

   call init_parallel()

   call print_program_info()

   if (nargs < 1 .or. 3 < nargs) then
      call print_usage()
      stop
   endif

   reset_step = .False.
   if (nargs == 3) then
      !! Check the third argument
      call get_command_argument(3, carg)
      if (trim(carg) == '--reset-step') then
         reset_step = .True.
      else
         call print_usage()
         stop
      endif
   endif

   !! Read input file
   call get_command_argument(1, cfile_inp)
   call read_input(cfile_inp)


   !! Open restart file if given
   restarted = .False.

   if (nargs >= 2) then
      call get_command_argument(2, cfile_rst)

      open(hdl_in_rst, file=cfile_rst, status='old', action='read', iostat=istat, form='unformatted', access='stream')

      if (istat > 0) then
         print '(2a)', 'Error: cannot open the restart file ', trim(cfile_rst)
         error stop
      endif

      restarted = .True.
   end if

   !! Initialise replicas
   call init_replica()

   !! Initialize the random number generator
   !call init_genrand64(rng_seed)
   call init_mts()

   !! Load force field
   call read_force_field()

   if (flg_pbc) then
      call init_pbc(restarted)
      ! RSTBLK%PBC and RSTBLK%PBCRESIZE will be read in the subroutine
   endif

   !! Read annealing file
   if (opt_anneal > 0) then
      call read_anneal(stat)
      if (.not. stat) then
         print *, 'Error in reading annealing-schedule file'
         stop (2)
      endif

      tempK = anneal_tempK(1)
   endif

   !! Prepare output files
   call init_out_files()


   !! Construct the sequences
   call init_sequence()

   !! Read coordinates (xyz)
   if (job /= JOBT%DCD) then

      allocate(xyz(3, nmp, nrep_proc))

      if (restarted) then
         call read_rst(RSTBLK%XYZ, rst_status)
      else
         call init_structure()
      endif
   endif

   print '(a)', 'Temperature'
   print '(a,f7.3)', '# T(K): ', tempK
   print '(a,f7.5)', '# kT(kcal/mol): ', kT
   print *
   flush(output_unit)

   !! Construct possible combinations of basepairs
   call init_bp()

   !! Allocation and initialisation of electrostatics
   if (flg_ele) call init_ele()

   if (flg_freeze) call init_freeze()

   if (flg_twz) call init_tweezers()

   if (flg_restraint) call read_restraint()

   if (flg_bias_rg) call init_bias_rg()

   !! Construct pair lists of local potentials
   call list_local()

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if (job == JOBT%MD) then
      ! Will use neighbor list
      continue

   else
      ! No neighbor list
      call list_bp()
      call list_exv()
      if (flg_ele) call list_ele()

   endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Main jobs
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (job == JOBT%DCD) then
      print '(a)', 'Starting DCD'
      call job_dcd()

   else if (job == JOBT%CHECK_FORCE) then

      call job_check_force()

   else if (job == JOBT%MD) then

      if (flg_replica) call write_rep_table()

      print '(a)', 'Starting MD'
      call job_md()

   endif

   do irep = 1, nrep_proc
      if (flg_out_bp) close(hdl_bp(irep))
      if (flg_out_bpe) close(hdl_bpe(irep))
      if (flg_out_bpall) close(hdl_bpall(irep))
   enddo

   deallocate(xyz)

   do irep = 1, nrep_proc
      close(hdl_out(irep))
   enddo
   deallocate(hdl_out)

   flush(output_unit)

   call end_parallel()

   stop

contains

   subroutine print_usage()

#ifdef PAR_MPI
      print *, 'Usage: unisis_mpi input.toml [restart file (.rst)] [--reset-step]'
#else
      print *, 'Usage: unisis input.toml [restart file (.rst)] [--reset-step]'
#endif

   end subroutine print_usage

   subroutine print_version()

#ifdef PAR_MPI
      print '(3a)', 'UNISIS ', UNISIS_RELEASE_LABEL, ' (MPI)'
#else
      print '(3a)', 'UNISIS ', UNISIS_RELEASE_LABEL, ' (serial)'
#endif

   end subroutine print_version

   subroutine print_program_info()

      use var_parallel, only : nprocs, nthreads, ncores
#ifdef PAR_MPI
      use var_parallel, only : myrank, MPI_VERSION, MPI_SUBVERSION
#endif

      character(len=40) :: githash, git
      character(len=8) :: date
      character(len=10) :: time
      character(len=5) :: zone
      character(len=500) :: com
      integer :: length

      call date_and_time(date, time, zone)
      call get_command(com, length)
      git = githash()

      print '(a)', '############ Program information ############'
      print '(a)', 'UNISIS model simulation program'
      print '(a)', 'Authors: N. Hori, H.T. Vu, J.A. Robins'
      print '(a)', 'Source: https://github.com/Hori-Lab/UNISIS'
      print '(2a)', 'UNISIS version: ', UNISIS_RELEASE_LABEL
#ifdef PAR_MPI
      print '(a)', 'Build variant: MPI'
#else
      print '(a)', 'Build variant: serial'
#endif
      if (git(1:1) == '?') then
         continue
      else
         print '(2a)', 'Git commit: ', git
      endif
      print '(2a)', 'Compiler version: ', compiler_version()
      print '(2a)', 'Compiler options: ', compiler_options()
      print '(2a)', 'Command: ', com(1:length)
      print '(15a)', 'Executed at ', date(1:4), '-', date(5:6), '-', date(7:8), &  ! Date
                     'T', time(1:2), ':', time(3:4), ':', time(5:6), zone(1:3), ':', zone(4:5)  ! ISO 8601 Time
      print '(a,i6)', 'Number of cores: ', ncores
      print '(a,i6)', 'Number of OpenMP threads: ', nthreads
      print '(a,i6)', 'Number of MPI processes: ', nprocs
#ifdef PAR_MPI
      print '(a,i2,a,i2)', 'MPI version: ', MPI_VERSION, '.', MPI_SUBVERSION
      print '(a,i6)', 'MPI rank: ', myrank
#endif
      print '(a)', '#############################################'
      print *

   end subroutine print_program_info

end program unisis
