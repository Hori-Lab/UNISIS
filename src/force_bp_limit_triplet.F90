subroutine force_bp_limit_triplet(irep, forces)

   use mt_stream, only : genrand_double1, genrand_double3
   use const, only : PREC
   use const_idx, only : REPT
   use const_phys, only : BOLTZ_KCAL_MOL
   use pbc, only : pbc_vec_d
   use var_top, only : nmp
   use var_state, only : xyz, bp_status, ene_bp, for_bp, kT, flg_bp_energy, nt_bp_excess, mts, tempK
   use var_potential, only : max_bp_per_nt, nbp, bp_cutoff_energy, bp_mp, bp_paras, bp_coef, &
                             basepair_parameters, flg_bias_ss, bias_ss_force
   use var_replica, only : flg_repvar, rep2val, irep2grep
#ifdef DUMPFORCE
   use const_idx, only : ENE
   use var_io, only : hdl_force
   use var_state, only: flg_step_dump_force
#endif

   implicit none

   integer, intent(in) :: irep
   real(PREC), intent(inout) :: forces(3, nmp)

   integer :: i, ibp, jbp, ibp_delete, nt_delete
   integer :: imp1, imp2, imp3, imp4, imp5, imp6
   integer :: i_save, i_swap
   integer :: nbp_seq
   integer :: bp_seq(nbp(irep))
   integer :: nnt_bp_excess
   integer :: ntlist_excess(nmp)
   ! Per-nucleotide reverse index for O(degree) BP lookup
   integer, parameter :: MAX_BP_DEGREE = 32
   integer :: nt_bp_count(nmp)
   integer :: nt_bp_idx(MAX_BP_DEGREE, nmp)
   integer :: nt_excess_pos(nmp)
   integer :: imp_del, jmp_del, pos, k
   type(basepair_parameters) :: bpp
   real(PREC) :: tK, dG
   real(PREC) :: u, pre, beta, e, ratio, rnd
   real(PREC) :: d, cosine, dih
   real(PREC) :: f_i(3), f_j(3), f_k(3), f_l(3)
   real(PREC) :: v12(3), v13(3), v42(3), v15(3), v62(3)
   real(PREC) :: a12
   real(PREC) :: d1212, d1313, d4242, d1213, d1242, d1215, d1515, d6262, d1262
   real(PREC) :: d1213over1212, d1242over1212, d1215over1212, d1262over1212
   real(PREC) :: m(3), n(3)
#ifdef DUMPFORCE
   real(PREC) :: force_save(3, 1:nmp)

   if (flg_step_dump_force) then
      !$omp master
      force_save(:,:) = 0.0_PREC
      !$omp end master
      !$omp barrier
   endif
#endif

   !real(PREC) :: for_bp(3, 6, nbp)
   ! In this subroutine, for_bp will be used in more than one loop.
   ! Therefore for_bp cannot be private and it has to be stored in var_state module.

   !#######################################
   ! imp-1 (3) --- imp (1) --- imp+1 (5)
   !                ||
   ! jmp+1 (6) --- jmp (2) --- jmp-1 (4)
   !#######################################

   if (flg_repvar(REPT%TEMP)) then
      tK = rep2val(irep2grep(irep), REPT%TEMP)
      beta = 1.0_PREC / (BOLTZ_KCAL_MOL * tK)
   else
      tK = tempK
      beta = 1.0_PREC / kT
   endif

   !$omp master
   bp_status(1:nbp(irep), irep) = .False.
   ene_bp(1:nbp(irep), irep) = 0.0_PREC
   for_bp(1:3,1:6,1:nbp(irep)) = 0.0_PREC
   nt_bp_excess(1:nmp) = -max_bp_per_nt
   !$omp end master

   ! Wait until the master initializes the arrays
   !$omp barrier

   !$omp do private(bpp, imp1, imp2, imp3, imp4, imp5, imp6, d, u, e, pre, cosine, dih, &
   !$omp&           f_i, f_j, f_k, f_l, v12, v13, v42, v15, v62, a12, m, n, &
   !$omp&           d1212, d1313, d4242, d1213, d1242, d1215, d1515, d6262, d1262, &
   !$omp&           d1213over1212, d1242over1212, d1215over1212, d1262over1212, &
   !$omp&           dG)
   do ibp = 1, nbp(irep)

      imp1 = bp_mp(1, ibp, irep)
      imp2 = bp_mp(2, ibp, irep)
      imp3 = bp_mp(3, ibp, irep)
      imp4 = bp_mp(4, ibp, irep)
      imp5 = bp_mp(5, ibp, irep)
      imp6 = bp_mp(6, ibp, irep)
      bpp = bp_paras(bp_mp(7, ibp, irep))

      ! dG = dH - T * dS
      dG = bp_coef(1, ibp, irep) - tK * bp_coef(2, ibp, irep)
      if (dG >= 0.0_PREC) cycle

      v12(:) = pbc_vec_d(xyz(:, imp1, irep), xyz(:, imp2, irep))
      d1212 = dot_product(v12,v12)
      a12 = sqrt(d1212)
      d = a12 - bpp%bond_r

      if (flg_bias_ss) then
         if (abs(d) > bpp%cutoff_ddist) then
            v12(:) = bias_ss_force / a12 * v12(:)
            !$omp critical
            forces(:, imp1) = forces(:, imp1) - v12(:)
            forces(:, imp2) = forces(:, imp2) + v12(:)
            !$omp end critical
            cycle
         endif
      else
         if (abs(d) > bpp%cutoff_ddist) cycle
      endif

      !===== Distance =====
      !d = a12 - bpp%bond_r

      u = bpp%bond_k * d**2
      f_i(:) = (2.0e0_PREC * bpp%bond_k * d / a12) * v12(:)
      for_bp(:, 1, ibp) = + f_i(:)
      for_bp(:, 2, ibp) = - f_i(:)

      v13(:) = pbc_vec_d(xyz(:, imp1, irep), xyz(:, imp3, irep))
      v15(:) = pbc_vec_d(xyz(:, imp1, irep), xyz(:, imp5, irep))
      v42(:) = pbc_vec_d(xyz(:, imp4, irep), xyz(:, imp2, irep))
      v62(:) = pbc_vec_d(xyz(:, imp6, irep), xyz(:, imp2, irep))
 
      d1313 = dot_product(v13, v13)
      d4242 = dot_product(v42, v42)
      d1213 = dot_product(v12, v13)
      d1242 = dot_product(v12, v42)
      d1215 = dot_product(v12, v15)
      d1515 = dot_product(v15, v15)
      d1262 = dot_product(v12, v62)
      d6262 = dot_product(v62, v62)
      d1213over1212 = d1213 / d1212
      d1242over1212 = d1242 / d1212
      d1215over1212 = d1215 / d1212
      d1262over1212 = d1262 / d1212

      !===== Angle of 3-1=2 (imp-1 -- imp -- jmp) =====
      cosine = d1213 / (sqrt(d1313) * a12)
      d = acos(cosine) - bpp%angl_theta2
      pre = 2.0e0_PREC * bpp%angl_k2 * d / sqrt(d1313*d1212 - d1213**2)
      u = u + bpp%angl_k2 * d**2

      f_i(:) = pre * (v12(:) - (d1213 / d1313 * v13(:)))
      f_k(:) = pre * (v13(:) - (d1213over1212 * v12(:)))

      for_bp(:, 1, ibp) = for_bp(:, 1, ibp) - f_i(:) - f_k(:)
      for_bp(:, 2, ibp) = for_bp(:, 2, ibp) + f_k(:)
      for_bp(:, 3, ibp) = f_i(:)

      !===== Angle of 1=2-4 (imp -- jmp -- jmp-1) =====
      cosine = d1242 / (a12 * sqrt(d4242))
      d = acos(cosine) - bpp%angl_theta1
      u = u + bpp%angl_k1 * d**2
      pre = 2.0e0_PREC * bpp%angl_k1 * d / sqrt(d1212*d4242 - d1242**2)

      f_i(:) = - pre * (v42(:) - (d1242over1212 * v12(:)))
      f_k(:) = - pre * (v12(:) - (d1242 / d4242 * v42(:)))

      for_bp(:, 1, ibp) = for_bp(:, 1, ibp) + f_i(:)
      for_bp(:, 2, ibp) = for_bp(:, 2, ibp) - f_i(:) - f_k(:)
      for_bp(:, 4, ibp) = f_k(:)

      !===== Angle of 5-1=2 (imp+1 -- imp -- jmp) =====
      cosine = d1215 / (sqrt(d1515) * a12)
      d = acos(cosine) - bpp%angl_theta4
      u = u + bpp%angl_k4 * d**2
      pre = 2.0e0_PREC * bpp%angl_k4 * d / sqrt(d1515*d1212 - d1215**2)

      f_i(:) = pre * (v12(:) - (d1215 / d1515 * v15(:)))
      f_k(:) = pre * (v15(:) - (d1215over1212 * v12(:)))

      for_bp(:, 1, ibp) = for_bp(:, 1, ibp) - f_i(:) - f_k(:)
      for_bp(:, 2, ibp) = for_bp(:, 2, ibp) + f_k(:)
      for_bp(:, 5, ibp) = f_i(:)

      !===== Angle of 1=2-6 (imp -- jmp -- jmp+1) =====
      cosine = d1262 / (a12 * sqrt(d6262))
      d = acos(cosine) - bpp%angl_theta3
      u = u + bpp%angl_k3 * d**2
      pre = 2.0e0_PREC * bpp%angl_k3 * d / sqrt(d1212*d6262 - d1262**2)

      f_i(:) = - pre * (v62(:) - (d1262over1212 * v12(:)))
      f_k(:) = - pre * (v12(:) - (d1262 / d6262 * v62(:)))

      for_bp(:, 1, ibp) = for_bp(:, 1, ibp) + f_i(:)
      for_bp(:, 2, ibp) = for_bp(:, 2, ibp) - f_i(:) - f_k(:)
      for_bp(:, 6, ibp) = f_k(:)
 
      !===== Dihedral angle among 4-2=1=3 (jmp-1 -- jmp -- imp -- imp-1) =====
      m(1) = v42(2)*v12(3) - v42(3)*v12(2)
      m(2) = v42(3)*v12(1) - v42(1)*v12(3)
      m(3) = v42(1)*v12(2) - v42(2)*v12(1)
      n(1) = v12(2)*v13(3) - v12(3)*v13(2)
      n(2) = v12(3)*v13(1) - v12(1)*v13(3)
      n(3) = v12(1)*v13(2) - v12(2)*v13(1)

      dih = atan2(dot_product(v42,n)*a12, dot_product(m,n))
      d = dih + bpp%dihd_phi1
      u = u + bpp%dihd_k1 * (1.0 + cos(d))

      pre = -bpp%dihd_k1 * sin(d) * a12
      f_i(:) = + pre / dot_product(m, m) * m(:)
      f_l(:) = - pre / dot_product(n, n) * n(:)

      for_bp(:, 4, ibp) = for_bp(:, 4, ibp) + f_i(:)
      for_bp(:, 2, ibp) = for_bp(:, 2, ibp) + (-1.0e0_PREC + d1242over1212) * f_i(:) &
                                        - (              d1213over1212) * f_l(:)
      for_bp(:, 1, ibp) = for_bp(:, 1, ibp) + (-1.0e0_PREC + d1213over1212) * f_l(:) &
                                        - (              d1242over1212) * f_i(:)
      for_bp(:, 3, ibp) = for_bp(:, 3, ibp) + f_l(:)
 
      !===== Dihedral angle among 6-2=1=5 (jmp+1 -- jmp -- imp -- imp+1) =====
      m(1) = v62(2)*v12(3) - v62(3)*v12(2)
      m(2) = v62(3)*v12(1) - v62(1)*v12(3)
      m(3) = v62(1)*v12(2) - v62(2)*v12(1)
      n(1) = v12(2)*v15(3) - v12(3)*v15(2)
      n(2) = v12(3)*v15(1) - v12(1)*v15(3)
      n(3) = v12(1)*v15(2) - v12(2)*v15(1)

      dih = atan2(dot_product(v62,n)*a12, dot_product(m,n))
      d = dih + bpp%dihd_phi2
      u = u + bpp%dihd_k2 * (1.0 + cos(d))

      pre = -bpp%dihd_k2 * sin(d) * a12
      f_i(:) = + pre / dot_product(m, m) * m(:)
      f_l(:) = - pre / dot_product(n, n) * n(:)

      for_bp(:, 6, ibp) = for_bp(:, 6, ibp) + f_i(:)
      for_bp(:, 2, ibp) = for_bp(:, 2, ibp) + (-1.0e0_PREC + d1262over1212) * f_i(:) &
                                        - (              d1215over1212) * f_l(:)
      for_bp(:, 1, ibp) = for_bp(:, 1, ibp) + (-1.0e0_PREC + d1215over1212) * f_l(:) &
                                        - (              d1262over1212) * f_i(:)
      for_bp(:, 5, ibp) = for_bp(:, 5, ibp) + f_l(:)

      !===== Total =====
      !e = bpp%U0 * bp_map_dG(imp1, imp2, irep) * exp(-u)
      e = bpp%U0 * dG * exp(-u)

      if (e <= bp_cutoff_energy) then
         bp_status(ibp, irep) = .True.
         ene_bp(ibp, irep) = e

         for_bp(:, :, ibp) = e * for_bp(:, :, ibp)

         !$omp atomic
         nt_bp_excess(imp1) = nt_bp_excess(imp1) + 1
         !$omp atomic
         nt_bp_excess(imp2) = nt_bp_excess(imp2) + 1
      endif
   enddo
   !$omp end do

   !$omp master

   ! Build per-nucleotide reverse index of active BPs
   nt_bp_count(:) = 0
   do ibp = 1, nbp(irep)
      if (bp_status(ibp, irep)) then
         imp_del = bp_mp(1, ibp, irep)
         jmp_del = bp_mp(2, ibp, irep)
         nt_bp_count(imp_del) = nt_bp_count(imp_del) + 1
         if (nt_bp_count(imp_del) > MAX_BP_DEGREE) then
            write(*, '(a,i0,a,i0)') 'WARNING: nt_bp_count exceeds MAX_BP_DEGREE for nt ', imp_del, &
               ', count=', nt_bp_count(imp_del)
         endif
         nt_bp_idx(nt_bp_count(imp_del), imp_del) = ibp
         nt_bp_count(jmp_del) = nt_bp_count(jmp_del) + 1
         if (nt_bp_count(jmp_del) > MAX_BP_DEGREE) then
            write(*, '(a,i0,a,i0)') 'WARNING: nt_bp_count exceeds MAX_BP_DEGREE for nt ', jmp_del, &
               ', count=', nt_bp_count(jmp_del)
         endif
         nt_bp_idx(nt_bp_count(jmp_del), jmp_del) = ibp
      endif
   enddo

   ! Build initial excess list with position tracking
   nnt_bp_excess = 0
   nt_excess_pos(:) = 0
   do imp1 = 1, nmp
      if (nt_bp_excess(imp1) > 0) then
         nnt_bp_excess = nnt_bp_excess + 1
         ntlist_excess(nnt_bp_excess) = imp1
         nt_excess_pos(imp1) = nnt_bp_excess
      endif
   enddo

   do     ! loop while (nnt_bp_excess > 0)

      if (nnt_bp_excess == 0) exit

      ! Randomly choose one nucleotide (nt) that has more than one base pair
      rnd = genrand_double3(mts(irep))    ! (0,1)-real-interval

      nt_delete = ntlist_excess( ceiling(rnd * nnt_bp_excess) )

      ! Look up BPs involving nt_delete via reverse index
      nbp_seq = nt_bp_count(nt_delete)
      bp_seq(1:nbp_seq) = nt_bp_idx(1:nbp_seq, nt_delete)

      ! Shuffle
      do i = 1, nbp_seq
         rnd = genrand_double3(mts(irep))   ! (0,1)-real-interval
         i_swap = ceiling(rnd*nbp_seq)
         i_save = bp_seq(i)
         bp_seq(i) = bp_seq(i_swap)
         bp_seq(i_swap) = i_save
      enddo

      ! Randomly choose one "ibp" that will be deleted, depending on the energies
      ibp_delete = bp_seq(1)
      do i = 2, nbp_seq
         jbp = bp_seq(i)

         ratio = exp( (ene_bp(jbp, irep) - ene_bp(ibp_delete, irep)) * beta )
         rnd = genrand_double1(mts(irep))  ! [0,1]-real-interval

         if (rnd < ratio) then
            ibp_delete = jbp
         endif
      enddo

      ! Delete
      bp_status(ibp_delete, irep) = .False.

      ! Update reverse index: remove ibp_delete from both nucleotides' lists
      imp_del = bp_mp(1, ibp_delete, irep)
      jmp_del = bp_mp(2, ibp_delete, irep)

      do k = 1, nt_bp_count(imp_del)
         if (nt_bp_idx(k, imp_del) == ibp_delete) then
            nt_bp_idx(k, imp_del) = nt_bp_idx(nt_bp_count(imp_del), imp_del)
            nt_bp_count(imp_del) = nt_bp_count(imp_del) - 1
            exit
         endif
      enddo
      do k = 1, nt_bp_count(jmp_del)
         if (nt_bp_idx(k, jmp_del) == ibp_delete) then
            nt_bp_idx(k, jmp_del) = nt_bp_idx(nt_bp_count(jmp_del), jmp_del)
            nt_bp_count(jmp_del) = nt_bp_count(jmp_del) - 1
            exit
         endif
      enddo

      ! Update nt_bp_excess
      nt_bp_excess(imp_del) = nt_bp_excess(imp_del) - 1
      nt_bp_excess(jmp_del) = nt_bp_excess(jmp_del) - 1

      ! Incrementally update ntlist_excess (shift-remove to preserve sorted order)
      if (nt_bp_excess(imp_del) == 0) then
         pos = nt_excess_pos(imp_del)
         do k = pos, nnt_bp_excess - 1
            ntlist_excess(k) = ntlist_excess(k + 1)
            nt_excess_pos(ntlist_excess(k)) = k
         enddo
         nt_excess_pos(imp_del) = 0
         nnt_bp_excess = nnt_bp_excess - 1
      endif
      if (nt_bp_excess(jmp_del) == 0) then
         pos = nt_excess_pos(jmp_del)
         if (pos > 0) then
            do k = pos, nnt_bp_excess - 1
               ntlist_excess(k) = ntlist_excess(k + 1)
               nt_excess_pos(ntlist_excess(k)) = k
            enddo
            nt_excess_pos(jmp_del) = 0
            nnt_bp_excess = nnt_bp_excess - 1
         endif
      endif

   enddo

   !$omp end master

   ! Wait until the master finishes deletions
   !$omp barrier

   !$omp do private(imp1, imp2, imp3, imp4, imp5, imp6)
   do ibp = 1, nbp(irep)

      if (.not. bp_status(ibp, irep)) cycle

      imp1 = bp_mp(1, ibp, irep)
      imp2 = bp_mp(2, ibp, irep)
      imp3 = bp_mp(3, ibp, irep)
      imp4 = bp_mp(4, ibp, irep)
      imp5 = bp_mp(5, ibp, irep)
      imp6 = bp_mp(6, ibp, irep)

      forces(:, imp1) = forces(:, imp1) + for_bp(:, 1, ibp)
      forces(:, imp2) = forces(:, imp2) + for_bp(:, 2, ibp)
      forces(:, imp3) = forces(:, imp3) + for_bp(:, 3, ibp)
      forces(:, imp4) = forces(:, imp4) + for_bp(:, 4, ibp)
      forces(:, imp5) = forces(:, imp5) + for_bp(:, 5, ibp)
      forces(:, imp6) = forces(:, imp6) + for_bp(:, 6, ibp)
#ifdef DUMPFORCE
      if (flg_step_dump_force) then
         do i = 1, 3
            !$omp atomic update
            force_save(i, imp1) = force_save(i, imp1) + for_bp(i, 1,ibp)
            !$omp atomic update
            force_save(i, imp2) = force_save(i, imp2) + for_bp(i, 2,ibp)
            !$omp atomic update
            force_save(i, imp3) = force_save(i, imp3) + for_bp(i, 3,ibp)
            !$omp atomic update
            force_save(i, imp4) = force_save(i, imp4) + for_bp(i, 4,ibp)
            !$omp atomic update
            force_save(i, imp5) = force_save(i, imp5) + for_bp(i, 5,ibp)
            !$omp atomic update
            force_save(i, imp6) = force_save(i, imp6) + for_bp(i, 6,ibp)
         enddo
      endif
#endif
   end do
   !$omp end do nowait

   !$omp master
   flg_bp_energy = .True.
   !$omp end master


#ifdef DUMPFORCE
   if (flg_step_dump_force) then
      !$omp barrier
      !$omp master
      do imp1 = 1, nmp
         write(hdl_force(ENE%BP), '(3(1x,e10.4))', advance='no') force_save(1:3, imp1)
      enddo
      write(hdl_force(ENE%BP), '(a)') ''
      !$omp end master
   endif
#endif

end subroutine force_bp_limit_triplet
