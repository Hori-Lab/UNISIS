subroutine read_force_field()

   use,intrinsic :: ISO_FORTRAN_ENV, only: OUTPUT_UNIT
   use tomlf

   use const_idx, only : NNT, nnt2char, NNENDT
   use const_phys, only : INVALID_VALUE, INVALID_JUDGE
   use var_potential
   use var_io, only : iopen_hdl, cfile_ff
   use var_parallel

   implicit none

   integer :: i
   integer :: istat
   integer :: hdl
   logical :: flg_angl
#ifdef PAR_MPI
   integer :: sz
#endif

   character(len=:), allocatable :: cline, csource

   !======= TOML
   type(toml_table), allocatable :: table
   type(toml_table), pointer :: group, node, subnode
   !type(toml_array), pointer :: array
   type(toml_context) :: context
   type(toml_error), allocatable :: tm_err
   !======= 

   call set_invalid()

   if (myrank == 0) then

      print '(2a)', 'Reading force-field file: ', trim(cfile_ff)
      flush(OUTPUT_UNIT)

      iopen_hdl = iopen_hdl + 1
      hdl = iopen_hdl

      open(hdl, file=cfile_ff, status='old', action='read', iostat=istat)

      if (istat /= 0) then
         print '(2a)', 'Error: failed to open the force-field file. ', trim(cfile_ff)
         flush(6)
         error stop
      endif

      call toml_load(table, hdl, context=context, error=tm_err)

      close(hdl)
      iopen_hdl = iopen_hdl - 1

      if (allocated(tm_err)) then
         print '(a)', tm_err%message
         call sis_abort()
      endif

      call get_value(table, "title", cline)
      print '(2a)', '# title: ', trim(cline)
      flush(OUTPUT_UNIT)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      !! Potential
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      call get_value(table, "potential", group)
      if (.not. associated(group)) then
         print '(a)', 'Error: [potential] required in FF file'
         call sis_abort()
      endif

      call get_value(group, "bond", node)
      if (associated(node)) then
         call get_value(node, "k", bond_k)
         call get_value(node, "r0", bond_r0)

      else
         print '(a)', 'Error: [bond] parameters required in FF file'
         call sis_abort()
      endif

      flg_angl = .True.
      call get_value(group, "angle", node, requested=.False.)
      if (associated(node)) then
         call get_value(node, "k", angl_k)
         call get_value(node, "a0", angl_t0)
      else
         flg_angl = .False.
      endif

      flg_angl_ReB = .True.
      call get_value(group, "angle_ReB", node, requested=.False.)
      if (associated(node)) then
         call get_value(node, "k", angl_k)
         call get_value(node, "a0", angl_t0)
      else
         flg_angl_ReB = .False.
      endif

      if (flg_angl .and. flg_angl_ReB) then
         print '(a)', 'Error: [angle] and [angle_ReB] cannot be specified together in FF file.'
         call sis_abort()

      else if (.not. flg_angl .and. .not. flg_angl_ReB) then
         print '(a)', 'Error: Either [angle] or [angle_ReB] parameters are required in FF file'
         call sis_abort()
      endif

      flg_dih_cos = .True.
      call get_value(group, "dihedral", node, requested=.False.)

      if (associated(node)) then
         call get_value(node, "k", dih_k)
         call get_value(node, "phi0", dih_p0)

      else
         flg_dih_cos = .False.
      endif

      flg_dih_exp = .True.
      call get_value(group, "dihedral_exp", node, requested=.False.)

      if (associated(node)) then
         call get_value(node, "k", dih_k)
         call get_value(node, "w", dih_w)
         call get_value(node, "phi0", dih_p0)

      else
         flg_dih_exp = .False.
      endif

      if (flg_angl .and. flg_angl_ReB) then
         print '(a)', 'Error: [dihedral] and [dihedral_exp] cannot be specified together in FF file.'
         call sis_abort()

      else if (.not. flg_angl .and. .not. flg_angl_ReB) then
         print '(a)', 'Neither [potential.dihedral] or [potential.dihedral_exp] was found, thus no dihedral potential set.'
         call sis_abort()
      endif

      call get_value(group, "basepair", node)
      bp_cutoff_energy = 0.001_PREC  ! Default /kcal/mol
      coef_dG = 1.0_PREC
      dH0 = 0.0_PREC
      dS0 = 0.0_PREC
      flg_NNend = .False.
      dHend0 = 0.0_PREC
      dSend0 = 0.0_PREC
      if (associated(node)) then
         !call get_value(node, "min_loop", bp_min_loop)
         !call get_value(node, "cutoff", bp_cutoff_dist)
         call get_value(node, "cutoff_energy", bp_cutoff_energy)

         call get_value(node, "coef_dG", coef_dG)
         call get_value(node, "dH0", dH0)
         call get_value(node, "dS0", dS0)

         call get_value(node, "NNend", flg_NNend)
         call get_value(node, "dHend0", dHend0)
         call get_value(node, "dSend0", dSend0)

         ! older format
         !call get_value(node, "bond_k", bp_bond_k)
         !call get_value(node, "bond_r", bp_bond_r)
         !call get_value(node, "angl_k", bp_angl_k)
         !call get_value(node, "angl_theta1", bp_angl_theta1)
         !call get_value(node, "angl_theta2", bp_angl_theta2)
         !call get_value(node, "dihd_k", bp_dihd_k)
         !call get_value(node, "dihd_phi1", bp_dihd_phi1)
         !call get_value(node, "dihd_phi2", bp_dihd_phi2)

         !bp_seqdep = 0
         !! = 0 (Default): No sequence dependence. Only U0_GC, U0_AU, U0_GU are required.
         !! = 1: Sequence dependent parameters. All possible combinations of trinucleotide-dimer are required.
         !call get_value(node, "seqdep", bp_seqdep)

         !if (bp_seqdep == 0) then
         !   call get_value(node, "U0_GC", bp_U0_GC)
         !   call get_value(node, "U0_AU", bp_U0_AU)
         !   call get_value(node, "U0_GU", bp_U0_GU)
         !endif

         call get_value(node, "GC", subnode)
         if (associated(subnode)) then
            call get_value(subnode, "U0", bp_paras(BPT%GC)%U0)
            call get_value(subnode, "bond_k", bp_paras(BPT%GC)%bond_k)
            call get_value(subnode, "bond_r", bp_paras(BPT%GC)%bond_r)
            call get_value(subnode, "angl_k1", bp_paras(BPT%GC)%angl_k1)
            call get_value(subnode, "angl_k2", bp_paras(BPT%GC)%angl_k2)
            call get_value(subnode, "angl_k3", bp_paras(BPT%GC)%angl_k3)
            call get_value(subnode, "angl_k4", bp_paras(BPT%GC)%angl_k4)
            call get_value(subnode, "angl_theta1", bp_paras(BPT%GC)%angl_theta1)
            call get_value(subnode, "angl_theta2", bp_paras(BPT%GC)%angl_theta2)
            call get_value(subnode, "angl_theta3", bp_paras(BPT%GC)%angl_theta3)
            call get_value(subnode, "angl_theta4", bp_paras(BPT%GC)%angl_theta4)
            call get_value(subnode, "dihd_k1", bp_paras(BPT%GC)%dihd_k1)
            call get_value(subnode, "dihd_k2", bp_paras(BPT%GC)%dihd_k2)
            call get_value(subnode, "dihd_phi1", bp_paras(BPT%GC)%dihd_phi1)
            call get_value(subnode, "dihd_phi2", bp_paras(BPT%GC)%dihd_phi2)

            bp_paras(BPT%CG) = bp_paras(BPT%GC)
            bp_paras(BPT%CG)%angl_k1 = bp_paras(BPT%GC)%angl_k2
            bp_paras(BPT%CG)%angl_k2 = bp_paras(BPT%GC)%angl_k1
            bp_paras(BPT%CG)%angl_k3 = bp_paras(BPT%GC)%angl_k4
            bp_paras(BPT%CG)%angl_k4 = bp_paras(BPT%GC)%angl_k3
         endif

         call get_value(node, "AU", subnode)
         if (associated(subnode)) then
            call get_value(subnode, "U0", bp_paras(BPT%AU)%U0)
            call get_value(subnode, "bond_k", bp_paras(BPT%AU)%bond_k)
            call get_value(subnode, "bond_r", bp_paras(BPT%AU)%bond_r)
            call get_value(subnode, "angl_k1", bp_paras(BPT%AU)%angl_k1)
            call get_value(subnode, "angl_k2", bp_paras(BPT%AU)%angl_k2)
            call get_value(subnode, "angl_k3", bp_paras(BPT%AU)%angl_k3)
            call get_value(subnode, "angl_k4", bp_paras(BPT%AU)%angl_k4)
            call get_value(subnode, "angl_theta1", bp_paras(BPT%AU)%angl_theta1)
            call get_value(subnode, "angl_theta2", bp_paras(BPT%AU)%angl_theta2)
            call get_value(subnode, "angl_theta3", bp_paras(BPT%AU)%angl_theta3)
            call get_value(subnode, "angl_theta4", bp_paras(BPT%AU)%angl_theta4)
            call get_value(subnode, "dihd_k1", bp_paras(BPT%AU)%dihd_k1)
            call get_value(subnode, "dihd_k2", bp_paras(BPT%AU)%dihd_k2)
            call get_value(subnode, "dihd_phi1", bp_paras(BPT%AU)%dihd_phi1)
            call get_value(subnode, "dihd_phi2", bp_paras(BPT%AU)%dihd_phi2)

            bp_paras(BPT%UA) = bp_paras(BPT%AU)
            bp_paras(BPT%UA)%angl_k1 = bp_paras(BPT%AU)%angl_k2
            bp_paras(BPT%UA)%angl_k2 = bp_paras(BPT%AU)%angl_k1
            bp_paras(BPT%UA)%angl_k3 = bp_paras(BPT%AU)%angl_k4
            bp_paras(BPT%UA)%angl_k4 = bp_paras(BPT%AU)%angl_k3
         endif

         call get_value(node, "GU", subnode)
         if (associated(subnode)) then
            call get_value(subnode, "U0", bp_paras(BPT%GU)%U0)
            call get_value(subnode, "bond_k", bp_paras(BPT%GU)%bond_k)
            call get_value(subnode, "bond_r", bp_paras(BPT%GU)%bond_r)
            call get_value(subnode, "angl_k1", bp_paras(BPT%GU)%angl_k1)
            call get_value(subnode, "angl_k2", bp_paras(BPT%GU)%angl_k2)
            call get_value(subnode, "angl_k3", bp_paras(BPT%GU)%angl_k3)
            call get_value(subnode, "angl_k4", bp_paras(BPT%GU)%angl_k4)
            call get_value(subnode, "angl_theta1", bp_paras(BPT%GU)%angl_theta1)
            call get_value(subnode, "angl_theta2", bp_paras(BPT%GU)%angl_theta2)
            call get_value(subnode, "angl_theta3", bp_paras(BPT%GU)%angl_theta3)
            call get_value(subnode, "angl_theta4", bp_paras(BPT%GU)%angl_theta4)
            call get_value(subnode, "dihd_k1", bp_paras(BPT%GU)%dihd_k1)
            call get_value(subnode, "dihd_k2", bp_paras(BPT%GU)%dihd_k2)
            call get_value(subnode, "dihd_phi1", bp_paras(BPT%GU)%dihd_phi1)
            call get_value(subnode, "dihd_phi2", bp_paras(BPT%GU)%dihd_phi2)

            bp_paras(BPT%UG) = bp_paras(BPT%GU)
            bp_paras(BPT%UG)%angl_k1 = bp_paras(BPT%GU)%angl_k2
            bp_paras(BPT%UG)%angl_k2 = bp_paras(BPT%GU)%angl_k1
            bp_paras(BPT%UG)%angl_k3 = bp_paras(BPT%GU)%angl_k4
            bp_paras(BPT%UG)%angl_k4 = bp_paras(BPT%GU)%angl_k3
         endif

      else
         print '(a)', 'Error: [basepair] parameters required in FF file'
         call sis_abort()
      endif

      call get_value(group, "wca", node)
      if (associated(node)) then
         call get_value(node, "sigma", wca_sigma)
         call get_value(node, "epsilon", wca_eps)

      else
         print '(a)', 'Error: [wca] parameters required in FF file'
         call sis_abort()
      endif


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      !! NN
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      if (bp_model == 4 .or. bp_model == 5) then
         call get_value(table, "NN", group)
         if (.not. associated(group)) then
            print '(a)', 'Error: [NN] required in FF file'
            call sis_abort()
         endif

         if (bp_model == 4) then

            call get_value(group, "dG", node)
            if (.not. associated(group)) then
               print '(a)', 'Error: [NN.dG] required in FF file'
               call sis_abort()
            endif

            allocate(NN_dG(NNT%MAX))
            NN_dG(:) = INVALID_VALUE

            do i = 1, NNT%MAX
               call get_value(node, nnt2char(i), NN_dG(i))
            enddo

         else if (bp_model == 5) then

            call get_value(group, "dH", node)
            if (.not. associated(group)) then
               print '(a)', 'Error: [NN.dH] required in FF file'
               call sis_abort()
            endif

            allocate(NN_dH(NNT%MAX))
            NN_dH(:) = INVALID_VALUE

            do i = 1, NNT%MAX
               call get_value(node, nnt2char(i), NN_dH(i))
            enddo


            call get_value(group, "dS", node)
            if (.not. associated(group)) then
               print '(a)', 'Error: [NN.dS] required in FF file'
               call sis_abort()
            endif

            allocate(NN_dS(NNT%MAX))
            NN_dS(:) = INVALID_VALUE

            do i = 1, NNT%MAX
               call get_value(node, nnt2char(i), NN_dS(i))
            enddo
         endif ! bp_model == 5
      endif ! NN

      if (bp_model == 5 .and. flg_NNend) then
         call get_value(table, "NNend", group)
         if (.not. associated(group)) then
            print '(a)', 'Error: [NNend] required in FF file'
            call sis_abort()
         endif

         call get_value(group, "dH", node)
         if (.not. associated(group)) then
            print '(a)', 'Error: [NNend.dH] required in FF file'
            call sis_abort()
         endif

         allocate(NNend_dH(1:6))
         NNend_dH(:) = INVALID_VALUE

         call get_value(node, 'AUonAU', NNend_dH(NNENDT%AUonAU))
         call get_value(node, 'AUonCG', NNend_dH(NNENDT%AUonCG))
         call get_value(node, 'AUonGU', NNend_dH(NNENDT%AUonGU))
         call get_value(node, 'GUonCG', NNend_dH(NNENDT%GUonCG))
         call get_value(node, 'GUonAU', NNend_dH(NNENDT%GUonAU))
         call get_value(node, 'GUonGU', NNend_dH(NNENDT%GUonGU))

         call get_value(group, "dS", node)
         if (.not. associated(group)) then
            print '(a)', 'Error: [NNend.dS] required in FF file'
            call sis_abort()
         endif

         allocate(NNend_dS(1:6))
         NNend_dS(:) = INVALID_VALUE

         call get_value(node, 'AUonAU', NNend_dS(NNENDT%AUonAU))
         call get_value(node, 'AUonCG', NNend_dS(NNENDT%AUonCG))
         call get_value(node, 'AUonGU', NNend_dS(NNENDT%AUonGU))
         call get_value(node, 'GUonCG', NNend_dS(NNENDT%GUonCG))
         call get_value(node, 'GUonAU', NNend_dS(NNENDT%GUonAU))
         call get_value(node, 'GUonGU', NNend_dS(NNENDT%GUonGU))

      endif ! NNend

   endif ! myrank == 0


#ifdef PAR_MPI

   if (myrank == 0) then
      print '(2a)', 'Sending force-field data via MPI.'
      flush(OUTPUT_UNIT)
   else
      print '(2a)', 'Receiving force-field data via MPI.'
      flush(OUTPUT_UNIT)
   endif

   !! Potentials

   call MPI_BCAST(bond_k, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(bond_r0, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)

   call MPI_BCAST(flg_angl, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(flg_angl_ReB, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(angl_k, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(angl_t0, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)

   call MPI_BCAST(flg_dih_cos, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(flg_dih_exp, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(dih_k, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(dih_w, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(dih_p0, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)

   call MPI_BCAST(bp_cutoff_energy, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(coef_dG, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(flg_NNend, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(dH0, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(dS0, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(dHend0, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(dSend0, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)

   sz = sizeof(bp_paras(1))
   call MPI_BCAST(bp_paras, sz*BPT%MAX, MPI_BYTE, 0, MPI_COMM_WORLD, istat)

   call MPI_BCAST(wca_sigma, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(wca_eps, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)


   !! NN
   if (bp_model == 4) then
      
      if (myrank /= 0) then
         allocate(NN_dG(NNT%MAX))
         NN_dG(:) = INVALID_VALUE
      endif

      call MPI_BCAST(NN_dG, NNT%MAX, PREC_MPI, 0, MPI_COMM_WORLD, istat)

   else if (bp_model == 5) then

      if (myrank /= 0) then
         allocate(NN_dH(NNT%MAX))
         NN_dH(:) = INVALID_VALUE

         allocate(NN_dS(NNT%MAX))
         NN_dS(:) = INVALID_VALUE
      endif

      call MPI_BCAST(NN_dH, NNT%MAX, PREC_MPI, 0, MPI_COMM_WORLD, istat)
      call MPI_BCAST(NN_dS, NNT%MAX, PREC_MPI, 0, MPI_COMM_WORLD, istat)

      !! NNend
      if (flg_NNend) then
         if (myrank /= 0) then
            allocate(NNend_dH(6))
            NNend_dH(:) = INVALID_VALUE

            allocate(NNend_dS(6))
            NNend_dS(:) = INVALID_VALUE
         endif

         call MPI_BCAST(NNend_dH, 6, PREC_MPI, 0, MPI_COMM_WORLD, istat)
         call MPI_BCAST(NNend_dS, 6, PREC_MPI, 0, MPI_COMM_WORLD, istat)
      endif

   endif

#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   !! Check
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

   if (myrank == 0) then
      csource = 'in the FF file.'
   else
      csource = 'received via MPI.'
   endif

   !! Potentials

   if (bond_k  > INVALID_JUDGE) then
      print '(a)', "INVALID bond_k " // trim(csource)
      call sis_abort()
   else
      print '(a,g15.8)', "# bond_k: ", bond_k
   endif

   if (bond_r0 > INVALID_JUDGE) then
      print '(a)', "INVALID bond_r0 " // trim(csource)
      call sis_abort()
   else
      print '(a,g15.8)', "# bond_r0: ", bond_r0
   endif

   if (angl_k  > INVALID_JUDGE) then
      print '(a)', "INVALID angl_k " // trim(csource)
      call sis_abort()
   else
      print '(a,g15.8)', "# angl_k: ", angl_k
   endif

   if (angl_t0 > INVALID_JUDGE) then
      print '(a)', "INVALID angl_t0 " // trim(csource)
      call sis_abort()
   else
      print '(a,g15.8)', "# angl_t0: ", angl_t0 
   endif

   if ((flg_dih_cos .or. flg_dih_exp) .and. dih_k  > INVALID_JUDGE) then
      print '(a)', "INVALID dih_k " // trim(csource)
      call sis_abort()
   else
      print '(a,g15.8)', "# dih_k: ", dih_k
   endif

   if (flg_dih_exp .and. dih_w  > INVALID_JUDGE) then
      print '(a)', "INVALID dih_w " // trim(csource)
      call sis_abort()
   else if (flg_dih_exp) then
      print '(a,g15.8)', "# dih_w: ", dih_w
   endif

   if ((flg_dih_cos .or. flg_dih_exp) .and. dih_p0 > INVALID_JUDGE) then
      print '(a)', "INVALID dih_p0 " // trim(csource)
      call sis_abort()
   else
      print '(a,g15.8)', "# dih_p0: ", dih_p0 
   endif

   print '(a,g15.8)', '# bp_cutoff_energy: ', bp_cutoff_energy
   print '(a,g15.8)', '# coef_dG: ', coef_dG
   print '(a,g15.8)', '# dH0: ', dH0
   print '(a,g15.8)', '# dS0: ', dS0
   if (flg_NNend) then
      print '(a)', '# NNend: True'
      print '(a,g15.8)', '# dHend0: ', dHend0
      print '(a,g15.8)', '# dSend0: ', dSend0
   endif

   call check_bp_paras(BPT%GC)
   call check_bp_paras(BPT%AU)
   call check_bp_paras(BPT%GU)
   call check_bp_paras(BPT%CG)
   call check_bp_paras(BPT%UA)
   call check_bp_paras(BPT%UG)

   if (wca_sigma > INVALID_JUDGE) then
      print '(a)', "INVALID wca_sigma " // trim(csource)
      call sis_abort()
   else
      print '(a,g15.8)', "# wca_sigma: ", wca_sigma
   endif

   if (wca_eps > INVALID_JUDGE) then
      print '(a)', "INVALID wca_eps " // trim(csource)
      call sis_abort()
   else
      print '(a,g15.8)', "# wca_eps: ", wca_eps
   endif

   !! NN

   if (bp_model == 4 .or. bp_model == 5) then
      print '(a)', '# NN parameters:'
   endif

   if (bp_model == 4) then

      print '(a)', '# dG:'
      do i = 1, NNT%MAX
         if (NN_dG(i) > INVALID_JUDGE) then
            print '(a)', 'Error: invalid NN.dG value for ' // nnt2char(i) // ' ' // trim(csource)
            call sis_abort()
         else
            print '(a3, a5, 2x, f7.3)', '#  ', nnt2char(i), NN_dG(i)
         endif
      enddo

   else if (bp_model == 5) then

      print '(a)', '# dH:'
      do i = 1, NNT%MAX
         if (NN_dH(i) > INVALID_JUDGE) then
            print '(a)', 'Error: invalid NN.dG value for ' // nnt2char(i) // ' ' // trim(csource)
            call sis_abort()
         else
            print '(a3, a5, 2x, f7.3)', '#  ', nnt2char(i), NN_dH(i)
         endif
      enddo

      print '(a)', '# dS:'
      do i = 1, NNT%MAX
         if (NN_dS(i) > INVALID_JUDGE) then
            print '(a)', 'Error: invalid NN.dS value for ' // nnt2char(i) // ' ' // trim(csource)
            call sis_abort()
         else
            print '(a3, a5, 2x, f7.3)', '#  ', nnt2char(i), NN_dS(i)
         endif
      enddo

      if (flg_NNend) then

         print '(a)', '# dHend:'
         do i = 1, 6
            if (NNend_dH(i) > INVALID_JUDGE) then
               print '(a)', 'Error: invalid NNend.dH value ' // trim(csource)
               call sis_abort()
            endif
         enddo
         print '(a8, 2x, f7.3)', '# AUonAU', NNend_dH(NNENDT%AUonAU)
         print '(a8, 2x, f7.3)', '# AUonCG', NNend_dH(NNENDT%AUonCG)
         print '(a8, 2x, f7.3)', '# AUonGU', NNend_dH(NNENDT%AUonGU)
         print '(a8, 2x, f7.3)', '# GUonCG', NNend_dH(NNENDT%GUonCG)
         print '(a8, 2x, f7.3)', '# GUonAU', NNend_dH(NNENDT%GUonAU)
         print '(a8, 2x, f7.3)', '# GUonGU', NNend_dH(NNENDT%GUonGU)

         print '(a)', '# dSend:'
         do i = 1, 6
            if (NNend_dS(i) > INVALID_JUDGE) then
               print '(a)', 'Error: invalid NNend.dS value ' // trim(csource)
               call sis_abort()
            endif
         enddo
         print '(a8, 2x, f7.3)', '# AUonAU', NNend_dS(NNENDT%AUonAU)
         print '(a8, 2x, f7.3)', '# AUonCG', NNend_dS(NNENDT%AUonCG)
         print '(a8, 2x, f7.3)', '# AUonGU', NNend_dS(NNENDT%AUonGU)
         print '(a8, 2x, f7.3)', '# GUonCG', NNend_dS(NNENDT%GUonCG)
         print '(a8, 2x, f7.3)', '# GUonAU', NNend_dS(NNENDT%GUonAU)
         print '(a8, 2x, f7.3)', '# GUonGU', NNend_dS(NNENDT%GUonGU)
      endif

   endif

   if (myrank == 0) then
      print '(a)', 'Done: reading force-field file'
   else
      print '(a)', 'Done: receiving force-field data'
   endif
   print *
   flush(OUTPUT_UNIT)

contains

   subroutine set_invalid()
      integer :: bptype

      bond_k  = INVALID_VALUE
      bond_r0 = INVALID_VALUE

      angl_k  = INVALID_VALUE
      angl_t0 = INVALID_VALUE

      dih_k  = INVALID_VALUE
      dih_w  = INVALID_VALUE
      dih_p0 = INVALID_VALUE

      !bp_min_loop = -1
      !bp_cutoff_dist = INVALID_VALUE
      !bp_cutoff_energy = INVALID_VALUE
      do bptype = 1, BPT%MAX
         bp_paras(bptype)%U0 = INVALID_VALUE
         bp_paras(bptype)%bond_k = INVALID_VALUE
         bp_paras(bptype)%bond_r = INVALID_VALUE
         bp_paras(bptype)%angl_k1 = INVALID_VALUE
         bp_paras(bptype)%angl_k2 = INVALID_VALUE
         bp_paras(bptype)%angl_k3 = INVALID_VALUE
         bp_paras(bptype)%angl_k4 = INVALID_VALUE
         bp_paras(bptype)%angl_theta1 = INVALID_VALUE
         bp_paras(bptype)%angl_theta2 = INVALID_VALUE
         bp_paras(bptype)%angl_theta3 = INVALID_VALUE
         bp_paras(bptype)%angl_theta4 = INVALID_VALUE
         bp_paras(bptype)%dihd_k1 = INVALID_VALUE
         bp_paras(bptype)%dihd_k2 = INVALID_VALUE
         bp_paras(bptype)%dihd_phi1 = INVALID_VALUE
         bp_paras(bptype)%dihd_phi2 = INVALID_VALUE
      enddo

      wca_sigma = INVALID_VALUE
      wca_eps = INVALID_VALUE
   endsubroutine set_invalid

   subroutine check_bp_paras(bptype)

      use const_idx, only : BPTYPE_CHAR
      integer, intent(in) :: bptype

      if (bp_paras(bptype)%U0 > INVALID_JUDGE) then
         print '(3a)', "INVALID value for bp ", BPTYPE_CHAR(bptype), " U0 " // trim(csource)
         call sis_abort()
      else
         print '(3a,g15.8)', "# ", BPTYPE_CHAR(bptype), " U0: ", bp_paras(bptype)%U0
      endif

      if (bp_paras(bptype)%bond_k > INVALID_JUDGE) then
         print '(3a)', "INVALID value for bp ", BPTYPE_CHAR(bptype), " bond_k " // trim(csource)
         call sis_abort()
      else
         print '(3a,g15.8)', "# ", BPTYPE_CHAR(bptype), " bond_k: ", bp_paras(bptype)%bond_k
      endif

      if (bp_paras(bptype)%bond_r > INVALID_JUDGE) then
         print '(3a)', "INVALID value for bp ", BPTYPE_CHAR(bptype), " bond_r " // trim(csource)
         call sis_abort()
      else
         print '(3a,g15.8)', "# ", BPTYPE_CHAR(bptype), " bond_r: ", bp_paras(bptype)%bond_r
      endif

      if (bp_paras(bptype)%angl_k1 > INVALID_JUDGE) then
         print '(3a)', "INVALID value for bp ", BPTYPE_CHAR(bptype), " angl_k1 " // trim(csource)
         call sis_abort()
      else
         print '(3a,g15.8)', "# ", BPTYPE_CHAR(bptype), " angl_k1: ", bp_paras(bptype)%angl_k1
      endif

      if (bp_paras(bptype)%angl_k2 > INVALID_JUDGE) then
         print '(3a)', "INVALID value for bp ", BPTYPE_CHAR(bptype), " angl_k2 " // trim(csource)
         call sis_abort()
      else
         print '(3a,g15.8)', "# ", BPTYPE_CHAR(bptype), " angl_k2: ", bp_paras(bptype)%angl_k2
      endif

      if (bp_paras(bptype)%angl_k3 > INVALID_JUDGE) then
         print '(2a)', "INVALID value in the force field file. bp angl_k3 for ", BPTYPE_CHAR(bptype)
         call sis_abort()
      else
         print '(3a,g15.8)', "# ", BPTYPE_CHAR(bptype), " angl_k3: ", bp_paras(bptype)%angl_k3
      endif

      if (bp_paras(bptype)%angl_k4 > INVALID_JUDGE) then
         print '(2a)', "INVALID value in the force field file. bp angl_k4 for ", BPTYPE_CHAR(bptype)
         call sis_abort()
      else
         print '(3a,g15.8)', "# ", BPTYPE_CHAR(bptype), " angl_k4: ", bp_paras(bptype)%angl_k4
      endif

      if (bp_paras(bptype)%angl_theta1 > INVALID_JUDGE) then
         print '(3a)', "INVALID value for bp ", BPTYPE_CHAR(bptype), " angl_theta1 " // trim(csource)
         call sis_abort()
      else
         print '(3a,g15.8)', "# ", BPTYPE_CHAR(bptype), " angl_theta1: ", bp_paras(bptype)%angl_theta1
      endif

      if (bp_paras(bptype)%angl_theta2 > INVALID_JUDGE) then
         print '(3a)', "INVALID value for bp ", BPTYPE_CHAR(bptype), " angl_theta2 " // trim(csource)
         call sis_abort()
      else
         print '(3a,g15.8)', "# ", BPTYPE_CHAR(bptype), " angl_theta2: ", bp_paras(bptype)%angl_theta2
      endif

      if (bp_paras(bptype)%angl_theta3 > INVALID_JUDGE) then
         print '(2a)', "INVALID value in the force field file. bp angl_theta3 for ", BPTYPE_CHAR(bptype)
         call sis_abort()
      else
         print '(3a,g15.8)', "# ", BPTYPE_CHAR(bptype), " angl_theta3: ", bp_paras(bptype)%angl_theta3
      endif

      if (bp_paras(bptype)%angl_theta4 > INVALID_JUDGE) then
         print '(2a)', "INVALID value in the force field file. bp angl_theta4 for ", BPTYPE_CHAR(bptype)
         call sis_abort()
      else
         print '(3a,g15.8)', "# ", BPTYPE_CHAR(bptype), " angl_theta4: ", bp_paras(bptype)%angl_theta4
      endif

      if (bp_paras(bptype)%dihd_k1 > INVALID_JUDGE) then
         print '(3a)', "INVALID value for bp ", BPTYPE_CHAR(bptype), " angl_dihd_k1 " // trim(csource)
         call sis_abort()
      else
         print '(3a,g15.8)', "# ", BPTYPE_CHAR(bptype), " dihd_k1: ", bp_paras(bptype)%dihd_k1
      endif

      if (bp_paras(bptype)%dihd_k2 > INVALID_JUDGE) then
         print '(3a)', "INVALID value for bp ", BPTYPE_CHAR(bptype), " angl_dihd_k2 " // trim(csource)
         call sis_abort()
      else
         print '(3a,g15.8)', "# ", BPTYPE_CHAR(bptype), " dihd_k2: ", bp_paras(bptype)%dihd_k2
      endif

      if (bp_paras(bptype)%dihd_phi1 > INVALID_JUDGE) then
         print '(3a)', "INVALID value for bp ", BPTYPE_CHAR(bptype), " angl_dihd_phi1 " // trim(csource)
         call sis_abort()
      else
         print '(3a,g15.8)', "# ", BPTYPE_CHAR(bptype), " dihd_phi1: ", bp_paras(bptype)%dihd_phi1
      endif

      if (bp_paras(bptype)%dihd_phi2 > INVALID_JUDGE) then
         print '(3a)', "INVALID value for bp ", BPTYPE_CHAR(bptype), " angl_dihd_phi2 " // trim(csource)
         call sis_abort()
      else
         print '(3a,g15.8)', "# ", BPTYPE_CHAR(bptype), " dihd_phi2: ", bp_paras(bptype)%dihd_phi2
      endif

   endsubroutine check_bp_paras

end subroutine read_force_field
