subroutine neighbor_list(irep)

   use const, only : PREC
   use const_idx, only : SEQT, MOLT
   use pbc, only : flg_pbc, pbc_vec_d, pbc_wrap
   use var_top, only : nmp_chain, imp_chain, nchains, nmp, has_charge, seq, moltypes
   use var_state, only : xyz, bp_status, bp_status_MC, nstep_bp_MC, &
                         ene_bp, for_bp, nt_bp_excess, nl_margin, lambdaD
   use var_potential, only : wca_sigma, nwca, nwca_max, wca_mp, &
                             bp_cutoff_dist, bp_mp, bp_coef, nbp, nbp_max, bp_map, &
                             ele_cutoff_type, ele_cutoff, ele_cutoff_inp, &
                             nele, nele_max, ele_mp, flg_ele, &
                             ele_exclude_covalent_bond_pairs, ele_exclude_covalent_angle_pairs, &
                             bp3_dH, bp3_dS, bp3_map, &
                             flg_bias_ss
   use var_replica, only : nrep_proc

   implicit none

   integer, intent(in) :: irep

   integer :: ichain, jchain, i, imp, j, jmp, j_start, iseq, jseq
   integer :: iwca, ibp, iele
   logical :: flg_add
   real(PREC) :: d2, v(3)
   real(PREC) :: wca_nl_cut2, bp_nl_cut2, ele_nl_cut2

   if (flg_pbc) then
      call pbc_wrap(irep)
   end if

   if (allocated(wca_mp)) then
      wca_mp(:,:,irep) = 0
   else
      allocate(nwca(nrep_proc))
      nwca_max = 5 * nmp
      allocate(wca_mp(2, nwca_max, nrep_proc))
      nwca(:) = 0
      wca_mp(:,:,:) = 0
   endif

   if (allocated(bp_mp)) then
      bp_mp(:,:,irep) = 0
      bp_coef(:,:,irep) = 0.0_PREC
   else
      allocate(nbp(nrep_proc))
      nbp_max = nmp / 2
      allocate(bp_mp(7, nbp_max, nrep_proc))
      allocate(bp_coef(2, nbp_max, nrep_proc))
      allocate(bp_status(nbp_max, nrep_proc))
      allocate(ene_bp(nbp_max, nrep_proc))
      allocate(for_bp(3, 6, nbp_max))
      allocate(nt_bp_excess(nmp))
      allocate(bp_status_MC(nbp_max, nrep_proc))
      nbp(:) = 0
      bp_mp(:,:,:) = 0
      bp_coef(:,:,:) = 0.0_PREC
      bp_status(:,:) = .False.
      ene_bp(:,:) = 0.0_PREC
      for_bp(:,:,:) = 0.0_PREC
      nt_bp_excess(:) = 0
      bp_status_MC(:,:) = .False.
   endif

   ele_nl_cut2 = 0.0_PREC
   if (flg_ele) then
      if (allocated(ele_mp)) then
         ele_mp(:,:,irep) = 0
         !ele_coef(:) = 0.0_PREC
      else
         allocate(nele(nrep_proc))
         nele_max = 5 * nmp
         allocate(ele_mp(2, nele_max, nrep_proc))
         nele(:) = 0
         ele_mp(:,:,:) = 0
         !allocate(ele_coef(nele_max))
         !ele_coef(:) = 0.0_PREC
      endif

      if (ele_cutoff_type == 2) then
         ele_cutoff(irep) = ele_cutoff_inp * lambdaD(irep)
      endif
      ! Cutoff
      ele_nl_cut2 = (ele_cutoff(irep) + nl_margin) ** 2
   endif

   ! Cutoff
   wca_nl_cut2 = (wca_sigma + nl_margin) ** 2
   if (flg_bias_ss) then
      bp_nl_cut2 = 9999999999.9
   else
      bp_nl_cut2 = (bp_cutoff_dist + nl_margin) ** 2
   endif

   iwca = 0
   ibp = 0
   iele = 0

   do ichain = 1, nchains 

      do jchain = ichain, nchains 

         do i = 1, nmp_chain(ichain)

            imp = imp_chain(i, ichain)
            iseq = seq(i, ichain)
   
            if (ichain == jchain) then
               j_start = i + 1
            else
               j_start = 1
            endif

            do j = j_start, nmp_chain(jchain)

               jmp = imp_chain(j, jchain)
               jseq = seq(j, jchain)

               v(:) = pbc_vec_d(xyz(:,imp,irep), xyz(:,jmp,irep))
               d2 = dot_product(v,v)

               ! WCA
               if (iseq /= SEQT%D .and. jseq /= SEQT%D) then
                  if (d2 <= wca_nl_cut2) then

                     flg_add = .True.
                     if (ichain == jchain) then
                        if (i + 2 >= j) then
                           flg_add = .False.

                        else if (moltypes(ichain) == MOLT%CIRCRNA) then
                           if (j + 2 >= i + nmp_chain(ichain)) then
                              flg_add = .False.
                           endif
                        endif
                     endif

                     if (flg_add) then
                        iwca = iwca + 1
                        if (iwca > nwca_max) then
                           call reallocate_wca_mp()
                        endif

                        wca_mp(1, iwca, irep) = imp
                        wca_mp(2, iwca, irep) = jmp
                     endif
                  endif
               endif

               ! Electrostatic
               if (flg_ele .and. d2 <= ele_nl_cut2) then

                  if (has_charge(imp) .and. has_charge(jmp)) then

                     flg_add = .True.
                     if (ichain == jchain) then
                        if (j == i + 1 .and. ele_exclude_covalent_bond_pairs) flg_add = .False.
                        if (j == i + 2 .and. ele_exclude_covalent_angle_pairs) flg_add = .False.
                        if (moltypes(ichain) == MOLT%CIRCRNA) then
                           if (i == 1) then
                              if (j == nmp_chain(ichain) .and. ele_exclude_covalent_bond_pairs) flg_add = .False.
                              if (j == nmp_chain(ichain) - 1 .and. ele_exclude_covalent_angle_pairs) flg_add = .False.
                            else if (i == 2) then
                              if (j == nmp_chain(ichain) .and. ele_exclude_covalent_angle_pairs) flg_add = .False.
                           endif
                        endif
                     endif

                     if (flg_add) then
                        iele = iele + 1
                        if (iele > nele_max) then
                           call reallocate_ele_mp()
                        endif

                        ele_mp(1, iele, irep) = imp
                        ele_mp(2, iele, irep) = jmp
                        !ele_coef(iele) = charge(i) * charge(j)
                     endif
                  endif
               endif

               ! Basepairs
               if (bp3_map(imp, jmp) > 0 .and. d2 <= bp_nl_cut2) then

                  ibp = ibp + 1
                  if (ibp > nbp_max) then
                     !write(*,*) 'Error: ibp > nbp_max. ibp =', ibp, 'nbp_max = ', nbp_max
                     call reallocate_bp_mp()
                  endif

                  bp_mp(1, ibp, irep) = imp
                  bp_mp(2, ibp, irep) = jmp

                  if (i == 1) then
                     bp_mp(3, ibp, irep) = imp + nmp_chain(ichain) - 1
                  else
                     bp_mp(3, ibp, irep) = imp - 1
                  endif

                  if (j == 1) then
                     bp_mp(4, ibp, irep) = jmp + nmp_chain(jchain) - 1
                  else
                     bp_mp(4, ibp, irep) = jmp - 1
                  endif

                  if (i == nmp_chain(ichain)) then
                     bp_mp(5, ibp, irep) = imp + 1 - nmp_chain(ichain)
                  else
                     bp_mp(5, ibp, irep) = imp + 1
                  endif

                  if (j == nmp_chain(jchain)) then
                     bp_mp(6, ibp, irep) = jmp + 1 - nmp_chain(jchain)
                  else
                     bp_mp(6, ibp, irep) = jmp + 1
                  endif

                  bp_mp(7, ibp, irep) = bp_map(imp, jmp)

                  bp_coef(1, ibp, irep) = bp3_dH(bp3_map(imp, jmp))  ! dH
                  bp_coef(2, ibp, irep) = bp3_dS(bp3_map(imp, jmp))  ! dS (0.001 already multiplied so that the unit is kcal/mol/K)

               endif

            enddo
         enddo

      enddo
   enddo

   nwca(irep) = iwca
   nbp(irep) = ibp
   if (flg_ele) nele(irep) = iele

contains
   
   subroutine reallocate_wca_mp()
      
      integer :: old_max
      integer :: tmp(2, nwca_max, nrep_proc)

      old_max = nwca_max   
      tmp(1:2, 1:old_max, 1:nrep_proc) = wca_mp(1:2, 1:old_max, 1:nrep_proc)

      deallocate(wca_mp)

      nwca_max = int(nwca_max * 1.2)

      allocate(wca_mp(2, nwca_max, nrep_proc))

      wca_mp(:, :, :) = 0
      wca_mp(1:2, 1:old_max, 1:nrep_proc) = tmp(1:2, 1:old_max, 1:nrep_proc)

   endsubroutine reallocate_wca_mp

   subroutine reallocate_ele_mp()

      integer :: old_max
      integer :: tmp(2, nele_max, nrep_proc)
      !real(PREC) :: tmp2(nele_max)

      old_max = nele_max   
      tmp(1:2, 1:old_max, 1:nrep_proc) = ele_mp(1:2, 1:old_max, 1:nrep_proc)
      !tmp2(1:old_max) = ele_coef(1:old_max)

      deallocate(ele_mp)
      !deallocate(ele_coef)

      nele_max = int(nele_max * 1.2)

      allocate(ele_mp(2, nele_max, nrep_proc))
      !allocate(ele_coef(nele_max))

      ele_mp(:, :, :) = 0
      ele_mp(1:2, 1:old_max, 1:nrep_proc) = tmp(1:2, 1:old_max, 1:nrep_proc)

      !ele_coef(:) = 0.0_PREC
      !ele_coef(1:old_max) = tmp2(1:old_max)

   endsubroutine reallocate_ele_mp

   subroutine reallocate_bp_mp()
      
      integer :: old_max
      integer :: tmp(7, nbp_max, nrep_proc)
      real(PREC) :: tmp2(2, nbp_max, nrep_proc)

      old_max = nbp_max   
      tmp(1:7, 1:old_max, 1:nrep_proc) = bp_mp(1:7, 1:old_max, 1:nrep_proc)
      tmp2(1:2, 1:old_max, 1:nrep_proc) = bp_coef(1:2, 1:old_max, 1:nrep_proc)

      deallocate(bp_mp)
      deallocate(bp_coef)
      deallocate(bp_status)
      deallocate(ene_bp)
      deallocate(for_bp)

      nbp_max = int(nbp_max * 1.2)

      allocate(bp_mp(7, nbp_max, nrep_proc))
      allocate(bp_coef(2, nbp_max, nrep_proc))
      allocate(bp_status(nbp_max, nrep_proc))
      allocate(ene_bp(nbp_max, nrep_proc))
      allocate(for_bp(3, 6, nbp_max))

      bp_mp(:, :, :) = 0
      bp_mp(1:7, 1:old_max, 1:nrep_proc) = tmp(1:7, 1:old_max, 1:nrep_proc)

      bp_coef(:, :, :) = 0.0_PREC
      bp_coef(1:2, 1:old_max, 1:nrep_proc) = tmp2(1:2, 1:old_max, 1:nrep_proc)

      bp_status(:,:) = .False.
      ene_bp(:,:) = 0.0_PREC
      for_bp(:,:,:) = 0.0_PREC

      if (nstep_bp_MC > 0) then
         deallocate(bp_status_MC)
         allocate(bp_status_MC(nbp_max, nrep_proc))
         bp_status_MC(:,:) = .False.
      endif

   endsubroutine reallocate_bp_mp

endsubroutine neighbor_list
