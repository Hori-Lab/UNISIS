subroutine job_md()

   use, intrinsic :: iso_fortran_env, Only : output_unit, INT64
   use const, only : PREC
   use const_phys, only : KCAL2JOUL, N_AVO, PI, BOLTZ_KCAL_MOL, JOUL2KCAL_MOL
   use const_idx, only : ENE, SEQT, RSTBLK, BPT, REPT, INTGRT
   use progress, only : progress_init, progress_update, wall_time_sec
   use pbc, only : flg_pbc, pbc_box, flg_pbc_resize, pbc_resize_step, pbc_resize
   use var_top, only : nmp, seq, mass, lmp_mp, ichain_mp, flg_freeze, is_frozen
   use var_state, only : restarted, reset_step, flg_bp_energy, ionic_strength, integrator, &
                         viscosity_Pas, xyz,  energies, dt, velos, accels, tempK, &
                         nstep, nstep_save, nstep_save_rst, &
                         nstep_check_stop, stop_wall_time_sec, fix_com_origin, &
                         nl_margin, Ekinetic, &
                         opt_anneal, nanneal, anneal_tempK, anneal_step, &
                         istep, ianneal, istep_anneal_next, &
                         nstep_bp_MC, flg_bp_MC, bp_status_MC, bp_status, rg
   use var_io, only : flg_progress, step_progress, hdl_dcd, cfile_dcd, hdl_rep, cfile_prefix
   use var_potential, only : stage_sigma, wca_sigma, bp_paras, bp_cutoff_energy, bp_cutoff_dist, &
                             ele_cutoff, flg_stage, flg_ele, &
                             flg_twz, ntwz_DCF, twz_DCF_forces, twz_DCF_direction, &
                             flg_bias_rg, flg_timed_bias_rg, ntimed_bias_rg, timed_bias_rg_k, timed_bias_rg_0, &
                             istep_timed_bias_rg_next, timed_bias_rg_step, itimed_bias_rg, &
                             bias_rg_k, bias_rg_0, flg_restraint
   use var_replica, only : nrep_all, nrep_proc, flg_replica, rep2val, irep2grep, rep2lab, &
                           nstep_rep_exchange, nstep_rep_save, nrep_all, flg_repvar, flg_exchange
   use var_parallel
   use dcd, only : file_dcd, DCD_OPEN_MODE
#ifdef DUMPFORCE
   use var_state, only : flg_step_dump_force
#endif

   implicit none

   integer :: i, irep, imp
   integer :: grep, rep_label
   integer :: rst_status
   real(PREC) :: tK, k, ionstr
   real(PREC) :: dxyz(3)
   real(PREC) :: xyz_move(3, nmp, nrep_proc)
   type(file_dcd) :: fdcd(nrep_proc)
   real(PREC) :: fric(nmp), radius
   real(PREC) :: v, c1, c2, md_coef(3, nmp), md_coef_rep(1, nmp, nrep_proc)
   real(PREC) :: rnd_bm(3, nmp)
   real(PREC) :: accels_pre(3)
   real(PREC) :: d2, d2max, d2max_2nd
   real(PREC), allocatable :: forces(:, :)
   real(PREC) :: replica_energies(2, nrep_all)
   logical :: flg_stop, flg_stop_file
   logical :: flg_step_save, flg_step_rep_save, flg_step_rep_exchange
   character(len=:), allocatable :: cfile_stop
#ifdef PAR_MPI
   integer :: istat
   real(PREC) :: replica_energies_l(2, nrep_all)
#endif

   ! Function
   real(PREC) :: rnd_boxmuller

   interface
   subroutine energy_sumup(irep, tempK_in, energies)
      use const, only : PREC
      use const_idx, only : ENE
      implicit none
      integer,    intent(in)  :: irep
      real(PREC), intent(in)  :: tempK_in
      real(PREC), intent(out) :: energies(0:ENE%MAX)
   endsubroutine energy_sumup
   subroutine set_ele(irep, tempk, ionstr, out_lb, out_Zp)
      use const, only : PREC
      integer, intent(in) :: irep
      real(PREC), intent(in) :: tempk
      real(PREC), intent(in) :: ionstr
      real(PREC), intent(out), optional :: out_lb
      real(PREC), intent(out), optional :: out_Zp
   endsubroutine set_ele
   endinterface

   cfile_stop = trim(cfile_prefix(:index(cfile_prefix, '/', .True.))) // 'STOP_UNISIS'

   allocate(mass(nmp))
   allocate(accels(3, nmp, nrep_proc))
   allocate(velos(3, nmp, nrep_proc))
   allocate(forces(3, nmp))
   allocate(energies(0:ENE%MAX, nrep_proc))
   allocate(Ekinetic(nrep_proc))
   !allocate(replica_energies(2, nrep_all))
   allocate(rg(nrep_proc))

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Set up constants for the dynamics
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if (integrator == INTGRT%LD_GJF2GJ) then
      print '(a)', 'Langevin dynamics parameters'
      print '(a)', 'Integrator: GFJ-2GJ'

      radius = 10.0_PREC
      v = viscosity_Pas * sqrt(1.0e3_PREC / KCAL2JOUL) * N_AVO * 1.0e-20_PREC
      fric(:) = 6.0e0_PREC * PI * v * radius

   else if (integrator == INTGRT%BD_EM) then
      print '(a)', 'Brownian dynamics parameters'
      print '(a)', 'Integrator: Ermak-McCammon'

      radius = 10.0_PREC
      v = (JOUL2KCAL_MOL * 1.0e-18_PREC) * viscosity_Pas
      fric(:) = 6.0e0_PREC * PI * v * radius
   endif

   print '(a,f10.3)', 'Stokes radius = ', radius
   print '(a,f10.3)', 'Viscosity = ', v
   print '(a,f10.3)', 'Friction coefficient = ', fric(1)

   ! Coefficients that do not depend on replicas
   do imp = 1, nmp
      ! Mass hard coded
      select case (seq(lmp_mp(imp), ichain_mp(imp)))
         case (SEQT%A)
            mass(imp) = 328.212_PREC
         case (SEQT%G)
            mass(imp) = 344.212_PREC
         case (SEQT%C)
            mass(imp) = 304.182_PREC
         case (SEQT%U)
            mass(imp) = 305.164_PREC
         case (SEQT%D)
            mass(imp) = 320.443_PREC
         case default
            write(*,*) 'Error: Unknown seq type, ', seq(lmp_mp(imp), ichain_mp(imp)), ' imp=',imp
            call sis_abort()
      endselect
   enddo

   print '(a)', 'mass(A) = 328.212'
   print '(a)', 'mass(G) = 344.212'
   print '(a)', 'mass(C) = 304.182'
   print '(a)', 'mass(U) = 305.164'
   print *
   flush(output_unit)

   call set_md_coef()


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! set up the initial state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (restarted) then
      call read_rst(RSTBLK%VELO, rst_status)
      call read_rst(RSTBLK%ACCEL, rst_status)

   else

      if (flg_stage) call check_int_stage()

      ! Set up initial velocities by Maxwell–Boltzmann distribution
      do irep = 1, nrep_proc
         if (flg_repvar(REPT%TEMP)) then
            tK = rep2val(irep2grep(irep), REPT%TEMP)
         else
            tK = tempK
         endif
         do imp = 1, nmp
            c1 = sqrt(tK * BOLTZ_KCAL_MOL / mass(imp))
            do i = 1, 3
               velos(i, imp, irep) = c1 * rnd_boxmuller(irep)
            enddo
         enddo

         ! Set up initial accerelations
         do imp = 1, nmp
            do i = 1, 3
               accels(i, imp, irep) = md_coef_rep(1, imp, irep) * rnd_boxmuller(irep)
            enddo
         end do
      enddo
   endif


   print '(a)', 'Potential and neighbor list cutoffs'
   print '(a,f10.3)', 'bp_cutoff_energy = ', bp_cutoff_energy
   print '(a,f10.3)', 'bp_cutoff_ddist(GC) = ', bp_paras(BPT%GC)%cutoff_ddist
   print '(a,f10.3)', 'bp_cutoff_ddist(AU) = ', bp_paras(BPT%AU)%cutoff_ddist
   print '(a,f10.3)', 'bp_cutoff_ddist(GU) = ', bp_paras(BPT%GU)%cutoff_ddist
   print '(a,f10.3)', 'bp_cutoff_dist(for neighbor list) = ', bp_cutoff_dist
   print '(a,f10.3)', 'wca_sigma = ', wca_sigma
   print '(a,f10.3)', 'stage_sigma = ', stage_sigma
   if (flg_ele) then
      do irep = 1, nrep_proc
         print '(a,i4,a,f10.3)', 'ele_cutoff(irep=',irep,') = ', ele_cutoff(irep)
      enddo
   endif
   print '(a,f10.3)', 'nl_margin = ', nl_margin
   print *

   do irep = 1, nrep_proc
      call neighbor_list(irep)
   enddo
   xyz_move(:,:,:) = 0.0e0_PREC

   if (restarted .and. .not. reset_step) then
      call read_rst(RSTBLK%STEP, rst_status)
   else
      istep = 0_INT64
   endif

   ! Setting up Simulated Annealing
   ianneal = 1
   if (opt_anneal > 0) then

      if (istep >= anneal_step(nanneal)) then
         ianneal = nanneal
         opt_anneal = 0

      else if (restarted) then
         call read_rst(RSTBLK%ANNEAL, rst_status)

      endif

      tempK = anneal_tempK(ianneal)

      call set_md_coef()
      !call set_bp_map()

      if (ianneal < nanneal) then
         istep_anneal_next = anneal_step(ianneal+1)
      endif
   endif

   ! Setting up Timed bias-Rg
   if (flg_timed_bias_rg) then

      if (flg_repvar(REPT%RG)) then
         write(*,*) 'Error: timed-bias-rg does not work with nrep_rg > 1.'
         call sis_abort()
      endif

      call read_timed_bias_rg()

      itimed_bias_rg = 1

      if (istep >= timed_bias_rg_step(ntimed_bias_rg)) then
         itimed_bias_rg = ntimed_bias_rg
         flg_timed_bias_rg = .False.

      else if (restarted) then

         do while (timed_bias_rg_step(itimed_bias_rg + 1) <= istep)

            itimed_bias_rg = itimed_bias_rg + 1

            if (itimed_bias_rg == ntimed_bias_rg) then
               flg_timed_bias_rg = .False.
               exit
            endif
         enddo
      endif

      bias_rg_k = timed_bias_rg_k(itimed_bias_rg)
      bias_rg_0(:) = timed_bias_rg_0(itimed_bias_rg)

      if (itimed_bias_rg < ntimed_bias_rg) then
         istep_timed_bias_rg_next = timed_bias_rg_step(itimed_bias_rg + 1)
      endif
   endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Initial energies
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (nstep_bp_MC > 0) then
      flg_bp_MC = .True.
   else
      flg_bp_MC = .False.
   endif

   do irep = 1, nrep_proc

      if (flg_replica) then
         grep = irep2grep(irep)
         rep_label = rep2lab(grep)
      endif

      if (flg_repvar(REPT%TEMP)) then
         tK = rep2val(grep, REPT%TEMP)
      else
         tK = tempK
      endif

      flg_bp_energy = .False.
      call energy_sumup(irep, tK, energies(0:ENE%MAX, irep))
      call energy_kinetic(irep, Ekinetic(irep))

      ! Open DCD file and write the header
      print '(a,i5,2a)', '# Opening dcd file for irep = ', irep, ', ', trim(cfile_dcd(irep))
      fdcd(irep) = file_dcd(hdl_dcd(irep), cfile_dcd(irep), DCD_OPEN_MODE%WRITE)

      ! set PBC box
      if (flg_pbc) then
         fdcd(irep)%flg_unitcell = .True.
         fdcd(irep)%box(:) = pbc_box(:)
      endif

      call fdcd(irep)%write_header(nmp)

      call write_out_header(irep)

      ! Output initial states

      ! STDOUT
      print '(a)', '##### Initial state'
      print '(a, i12)', 'Step ', istep
      if (flg_replica) then
         print '(a, i4)', 'Replica ID (in-process) ', irep
         print '(a, i4)', 'Replica ID (global) ', grep
         print '(a, i4)', 'Replica label ', rep_label
      endif
      print '(a, f6.2)',  'Temperature  ', tK
      print '(a, g13.6)', 'E_kinetic ', Ekinetic(irep)
      print '(a, g13.6)', 'E_total   ', energies(ENE%TOTAL, irep)
      print '(a, g13.6)', 'E_bond    ', energies(ENE%BOND, irep)
      print '(a, g13.6)', 'E_angl    ', energies(ENE%ANGL, irep)
      print '(a, g13.6)', 'E_dih     ', energies(ENE%DIHE, irep)
      print '(a, g13.6)', 'E_bp      ', energies(ENE%BP, irep)
      if (flg_ele  ) print '(a, g13.6)', 'E_ele     ', energies(ENE%ELE, irep)
      if (flg_stage) print '(a, g13.6)', 'E_stage   ', energies(ENE%STAGE, irep)
      if (flg_twz  ) print '(a, g13.6)', 'E_tweezers', energies(ENE%TWZ, irep)
      if (flg_bias_rg) print '(a, g13.6)', 'E_rg      ', energies(ENE%RG, irep)
      if (flg_bias_rg) print '(a, g13.6)', 'Rg        ', rg(irep)
      if (flg_restraint) print '(a, g13.6)', 'E_rest    ', energies(ENE%REST, irep)
      print *

      if (.not. restarted) then
         ! At istep = 0 (not restarted), write both .out and DCD
         call write_out(irep, istep, rep_label, tK)

         call write_bp(irep, tK)

         call fdcd(irep)%write_onestep(nmp, xyz(:,:,irep), fix_com_origin)
      endif
   enddo
   print *
   flush(output_unit)

   if (nstep_bp_MC > 0) then
      bp_status_MC(:,:) = bp_status(:,:)
   else
      bp_status_MC(:,:) = .True.
   endif

   if (myrank == 0 .and. flg_progress) then
      call progress_init(istep)
   endif

   flg_stop = .False.
   flg_stop_file = .False.

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Main loop for time integration
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do while (istep < nstep)

      istep = istep + 1

      flg_bp_MC = .False.
      if (nstep_bp_MC > 0) then
         if (mod(istep, nstep_bp_MC) == 0) flg_bp_MC = .True.
      endif

#ifdef DUMPFORCE
      flg_step_dump_force = .False.
      if (mod(istep-1, nstep_save) == 0) flg_step_dump_force = .True.
#endif

      do irep = 1, nrep_proc

         d2max = 0.0e0_PREC
         d2max_2nd = 0.0e0_PREC
         do imp = 1, nmp
            d2 = dot_product(xyz_move(1:3, imp, irep), xyz_move(1:3, imp, irep))
            if (d2 > d2max) then
               d2max_2nd = d2max
               d2max = d2
            else if (d2 > d2max_2nd) then
               d2max_2nd = d2
            endif
         enddo

         if (sqrt(d2max) + sqrt(d2max_2nd) > nl_margin) then
            call neighbor_list(irep)
            xyz_move(:, :, irep) = 0.0e0_PREC

            if (nstep_bp_MC > 0) then
               flg_bp_MC = .True.  ! Enforce MC after updating the neighbor list
            endif
         endif

         call force(irep, forces(:,:))

         if (flg_bp_MC) then
            bp_status_MC(:,irep) = bp_status(:,irep)
         endif

         if (maxval(forces(:,:)) > 100.0) then
            print *, 'Warning: maxval(forces) > 100.0. (irep, grep, istep, force) =', &
                     irep, irep2grep(irep), istep, maxval(forces)
            flush(output_unit)
         endif

         do imp= 1, nmp
            do i = 1, 3
               rnd_bm(i, imp) = rnd_boxmuller(irep)
            enddo
         enddo

         if (integrator == INTGRT%LD_GJF2GJ) then
            !$omp parallel do private(accels_pre, dxyz)
            do imp = 1, nmp
               !if(fix_mp(imp)) cycle

               !! a = (1 - gamma h / 2m) / (1 + gamma h / 2m)
               !! b = 1 / (1 + gamma h / 2m)

               ! beta(t + h) with the associated coefficient
               accels_pre(1:3) = md_coef_rep(1, imp, irep) * rnd_bm(1:3, imp)
               ! md_coef_rep(1) = sqrt(b) / 2m * sqrt(2 gamma kT h)

               ! v(t + 1/2h) update the half-step velocity
               velos(1:3, imp, irep) =  md_coef(1, imp) * velos(1:3, imp, irep)  &
                                + md_coef(2, imp) * forces(1:3, imp) &
                                + (accels(1:3, imp, irep) + accels_pre(1:3))
               ! md_coef(1) = a
               ! md_coef(2) = sqrt(b) h / m

               ! beta(t) <= beta(t+h) (incluing the coefficient) save for the next iteration
               accels(1:3, imp, irep) = accels_pre(1:3)

               dxyz(1:3) =  md_coef(3, imp) * velos(1:3, imp, irep)
               ! md_coef(3) = sqrt(b) h
               !           (= 0 if frozen)

               xyz(1:3, imp, irep) = xyz(1:3, imp, irep) + dxyz(1:3)
               xyz_move(1:3, imp, irep) = xyz_move(1:3, imp, irep) + dxyz(1:3)
            end do
            !$omp end parallel do

         else if (integrator == INTGRT%BD_EM) then
            !$omp parallel do private(dxyz)
            do imp = 1, nmp

               dxyz(1:3) = md_coef(1, imp) * forces(1:3, imp) &
                          + md_coef_rep(1, imp, irep) * rnd_bm(1:3, imp)

               xyz(1:3, imp, irep) = xyz(1:3, imp, irep) + dxyz(1:3)
               xyz_move(1:3, imp, irep) = xyz_move(1:3, imp, irep) + dxyz(1:3)
            enddo
            !$omp end parallel do

         endif ! integrator

      enddo ! irep

      flg_step_save = .False.
      flg_step_rep_save = .False.
      flg_step_rep_exchange = .False.
      if (mod(istep, nstep_save) == 0) flg_step_save = .True.
      if (flg_replica .and. flg_exchange) then
         if (mod(istep, nstep_rep_exchange) == 0) flg_step_rep_exchange = .True.
         if (mod(istep, nstep_rep_save) == 0) flg_step_rep_save = .True.
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Energy calculation
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (flg_step_save .or. flg_step_rep_exchange) then
         replica_energies(:, :) = 0.0_PREC
         if (flg_step_rep_exchange .and. nstep_bp_MC > 0) then
            flg_bp_MC = .True.
            flg_bp_energy = .False.
         endif
         call energy_replica(energies, replica_energies, flg_step_rep_exchange, flg_step_save)
      endif


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Save energies and coordinates to out, dcd, bp files
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (flg_step_save) then
         do irep = 1, nrep_proc
            if (flg_replica) then
               grep = irep2grep(irep)
               rep_label = rep2lab(grep)
            endif

            if (flg_repvar(REPT%TEMP)) then
               tK = rep2val(grep, REPT%TEMP)
            else
               tK = tempK
            endif

            !call energy_sumup(irep, energies(1:ENE%MAX, irep))
            call energy_kinetic(irep, Ekinetic(irep))

            call write_out(irep, istep, rep_label, tK)

            call write_bp(irep, tK)

            call fdcd(irep)%write_onestep(nmp, xyz(:,:,irep), fix_com_origin)
         enddo
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Replica exchange
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (flg_replica) then

         ! Write rep file
         if (myrank == 0 .and. flg_step_rep_save) then
            write(hdl_rep, '(i12,1x)', ADVANCE='no') istep

            do irep = 1, nrep_all - 1
               write(hdl_rep, '(i5,1x)', ADVANCE='no') rep2lab(irep)
            enddo

            write(hdl_rep, '(i5)') rep2lab(nrep_all)
         endif

         ! Exchange
         if (flg_step_rep_exchange) then

#ifdef PAR_MPI
            replica_energies_l(:,:) = replica_energies(:,:)
            call MPI_ALLREDUCE(replica_energies_l, replica_energies, 2*nrep_all, PREC_MPI, &
                               MPI_SUM, MPI_COMM_WORLD, istat)
#endif

            call replica_exchange(velos, replica_energies, tempk)

            if (flg_repvar(REPT%TEMP)) call set_md_coef()

            do irep = 1, nrep_proc
               grep = irep2grep(irep)

               if (flg_ele) then
                  if (flg_repvar(REPT%TEMP)) then
                     tK = rep2val(grep, REPT%TEMP)
                  else
                     tK = tempK
                  endif

                  if (flg_repvar(REPT%ION)) then
                     ionstr = rep2val(grep, REPT%ION)
                  else
                     ionstr = ionic_strength
                  endif

                  call set_ele(irep, tK, ionstr)
               endif

               if (flg_repvar(REPT%TWZDCF)) then
                  k = rep2val(grep, REPT%TWZDCF)
                  twz_DCF_forces(:, 1:ntwz_DCF, irep) = k * twz_DCF_direction(:, 1:ntwz_DCF)
               endif

               if (flg_repvar(REPT%RG)) then
                  bias_rg_0(irep) = rep2val(grep, REPT%RG)
               endif
            enddo

         endif
      endif


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Update the box size
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (flg_pbc_resize) then
         if (mod(istep, pbc_resize_step) == 0) then

            call pbc_resize()

            do irep = 1, nrep_proc
               fdcd(irep)%box(:) = pbc_box(:)

               call neighbor_list(irep)
               xyz_move(:,:,irep) = 0.0e0_PREC
            enddo
            print '(a,i12,a,3(1x,f8.3))', 'PBC resized: step = ',istep, ', box size =', pbc_box(1:3)

            ! If the flag is set to False, this is the last update.
            if (.not. flg_pbc_resize) then
               print '(a)', 'The target box size has been reached.'
            endif
         endif
      endif


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Annealing
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (opt_anneal > 0 .and. istep + 1 == istep_anneal_next) then
         ianneal = ianneal + 1
         tempK = anneal_tempK(ianneal)

         call set_md_coef()
         !call set_bp_map()

         if (ianneal < nanneal) then
            istep_anneal_next = anneal_step(ianneal+1)
         else
            opt_anneal = 0
         endif
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Timed bias-Rg
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (flg_timed_bias_rg .and. istep + 1 == istep_timed_bias_rg_next) then
         itimed_bias_rg = itimed_bias_rg + 1
         bias_rg_k = timed_bias_rg_k(itimed_bias_rg)
         bias_rg_0(:) = timed_bias_rg_0(itimed_bias_rg)

         if (itimed_bias_rg < ntimed_bias_rg) then
            istep_timed_bias_rg_next = timed_bias_rg_step(itimed_bias_rg + 1)
         else
            flg_timed_bias_rg = .False.
         endif
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Write progress
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (myrank == 0 .and. flg_progress) then
         if (mod(istep, step_progress) == 0) then
            call progress_update(istep, nstep)
         endif
      endif


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Check wall time
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (mod(istep, nstep_check_stop) == 0) then
         if (myrank == 0) then
            inquire(file=trim(cfile_stop), exist=flg_stop_file)
            if (flg_stop_file) flg_stop = .True.

            if (stop_wall_time_sec > 0) then
               if (wall_time_sec() > stop_wall_time_sec) then
                  flg_stop = .True.
               endif
            endif
         endif
#ifdef PAR_MPI
         call MPI_BCAST(flg_stop, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
         if (flg_stop) then
            call MPI_BCAST(flg_stop_file, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
         endif
#endif
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Write restart file (.rst)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (mod(istep, nstep_save_rst) == 0 .or. flg_stop) then
         call write_rst()
      endif

      if (flg_stop) then
         do irep = 1, nrep_proc
            call fdcd(irep)%close()
         enddo

         if (flg_stop_file) then
            print '(a,i13,a)', 'STOP_UNISIS file was found at step ', istep, '. Stop the job.'
         else
            print '(a,i13,a)', 'Wall-clock time limit reached at step ', istep, '. Stop the job.'
         endif
         flush(output_unit)

         exit
      endif

   enddo  ! <--- Main loop for time integration

#ifdef DUMPFORCE
   if (mod(nstep, nstep_save) == 0) then
      flg_step_dump_force = .True.
      do irep = 1, nrep_proc
         call force(irep, forces(:,:))
      enddo
   endif
#endif

   do irep = 1, nrep_proc
      call fdcd(irep)%close()
   enddo

contains

   subroutine set_md_coef()

      logical, save :: flg_first = .True.

      ! Coefficients that do not depend on replicas
      if (flg_first) then

         if (integrator == INTGRT%LD_GJF2GJ) then
            do imp = 1, nmp
               !! sqrt(b) = sqrt(1 / (1 + gamma h / 2m))
               c1 = 0.5 * dt * fric(imp) / mass(imp)
               c2 = sqrt(1.0e0_PREC / (1.0_PREC + c1))

               ! md_coef(1) = a
               md_coef(1, imp) = (1.0_PREC - c1) / (1.0_PREC + c1)
               ! md_coef(2) = sqrt(b) h / m
               md_coef(2, imp) = c2 * dt / mass(imp)
               ! md_coef(3) = sqrt(b) h
               md_coef(3, imp) = c2 * dt

               if (flg_freeze) then
                  if (is_frozen(imp)) md_coef(3, imp) = 0.0_PREC
               endif
            enddo

         else if (integrator == INTGRT%BD_EM) then
            do imp = 1, nmp
               md_coef(1, imp) = dt / fric(imp)

               if (flg_freeze) then
                  if (is_frozen(imp)) md_coef(1, imp) = 0.0_PREC
               endif
            enddo
         endif

         flg_first = .False.
      endif

      ! Coefficients depending on replicas
      do irep = 1, nrep_proc

         if (flg_repvar(REPT%TEMP)) then
            tK = rep2val(irep2grep(irep), REPT%TEMP)
         else
            tK = tempK
         endif

         if (integrator == INTGRT%LD_GJF2GJ) then
            do imp = 1, nmp
               c1 = 0.5 * dt * fric(imp) / mass(imp)
               c2 = sqrt(1.0e0_PREC / (1.0_PREC + c1))
               ! md_coef_rep(1) = sqrt(b) / 2m * sqrt(2 gamma kT h)
               md_coef_rep(1, imp, irep) = 0.5_PREC * c2 / mass(imp) * sqrt(2.0_PREC * fric(imp) * BOLTZ_KCAL_MOL * tK * dt)
            enddo

         else if (integrator == INTGRT%BD_EM) then
            do imp = 1, nmp
               md_coef_rep(1, imp, irep) = sqrt(2.0_PREC * BOLTZ_KCAL_MOL * tK * dt / fric(imp))

               if (flg_freeze) then
                  if (is_frozen(imp)) md_coef_rep(1, imp, irep) = 0.0_PREC
               endif
            enddo
         endif

      enddo

   endsubroutine set_md_coef

endsubroutine job_md
