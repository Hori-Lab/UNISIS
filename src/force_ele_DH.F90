subroutine force_ele_DH(irep, forces)

   use const, only : PREC
   use pbc, only : pbc_vec_d
   use var_top, only : nmp
   use var_state, only : xyz, lambdaD
   use var_potential, only : ele_coef, ele_cutoff, nele, ele_mp
#ifdef DUMPFORCE
   use const_idx, only : ENE
   use var_io, only : hdl_force
   use var_state, only: flg_step_dump_force
#endif

   implicit none

   integer, intent(in) :: irep
   real(PREC), intent(inout) :: forces(3, nmp)

   integer :: iele, imp1, imp2
   real(PREC) :: dist1, dist2, rdist1
   real(PREC) :: dvdw_dr, rcdist, cutoff2
   real(PREC) :: v21(3), for(3)
   !character(ARRAY_MSG_ERROR) :: error_message
#ifdef DUMPFORCE
   integer :: i
   real(PREC) :: force_save(3, 1:nmp)

   if (flg_step_dump_force) then
      !$omp master
      force_save(:,:) = 0.0_PREC
      !$omp end master
      !$omp barrier
   endif
#endif


   ! --------------------------------------------------------------------

   cutoff2 = ele_cutoff(irep) ** 2
   rcdist = 1.0_PREC / lambdaD(irep)

   !$omp do private(imp1,imp2,v21,dist2,dist1,rdist1,dvdw_dr,for)
   do iele = 1, nele(irep)
      imp1 = ele_mp(1, iele, irep)
      imp2 = ele_mp(2, iele, irep)

      v21(:) = pbc_vec_d(xyz(:, imp2, irep), xyz(:, imp1, irep))

      dist2 = dot_product(v21,v21)

      if(dist2 > cutoff2) cycle

      ! -----------------------------------------------------------------
      dist1 = sqrt(dist2)
      rdist1 = 1.0_PREC / dist1

      dvdw_dr = ele_coef(irep) * rdist1 * rdist1 &
               * (rdist1 + rcdist) * exp(-dist1 * rcdist)

     
      !if(dvdw_dr > DE_MAX) then
      !   write(error_message,*) 'force_ele_DH > DE_MAX', istep, grep, imp1, imp2, dist1, dvdw_dr, DE_MAX
      !   write(6, '(a)') trim(error_message)
      !   dvdw_dr = DE_MAX
      !end if

      for(1:3) = dvdw_dr * v21(1:3)
      forces(1:3, imp1) = forces(1:3, imp1) - for(1:3)
      forces(1:3, imp2) = forces(1:3, imp2) + for(1:3)

#ifdef DUMPFORCE
      if (flg_step_dump_force) then
         do i = 1, 3
            !$omp atomic update
            force_save(i, imp1) = force_save(i, imp1) - for(i)
            !$omp atomic update
            force_save(i, imp2) = force_save(i, imp2) + for(i)
         enddo
      endif
#endif
   end do
   !$omp end do nowait


#ifdef DUMPFORCE
   if (flg_step_dump_force) then
      !$omp barrier
      !$omp master
      do imp1 = 1, nmp
         write(hdl_force(ENE%ELE), '(3(1x,e10.4))', advance='no') force_save(1:3, imp1)
      enddo
      write(hdl_force(ENE%ELE), '(a)') ''
      !$omp end master
   endif
#endif

end subroutine force_ele_DH
