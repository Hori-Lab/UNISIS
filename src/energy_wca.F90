subroutine energy_wca(irep, Ewca)

   use const
   use pbc, only : pbc_vec_d
   use var_state, only : xyz
   use var_potential, only : wca_sigma, wca_eps, wca_mp, nwca

   implicit none

   integer, intent(in) :: irep
   real(PREC), intent(inout) :: Ewca

   integer :: iwca, imp1, imp2
   real(PREC) :: d
   real(PREC) :: e_wca

   e_wca = 0.0e0_PREC

   !$omp parallel do private(imp1,imp2,d) reduction(+:e_wca)
   do iwca = 1, nwca(irep)

      imp1 = wca_mp(1, iwca, irep)
      imp2 = wca_mp(2, iwca, irep)

      d = norm2( pbc_vec_d(xyz(:,imp1,irep), xyz(:, imp2,irep)) )

      if (d >= wca_sigma) cycle

      d = wca_sigma / d
      e_wca = e_wca + wca_eps * (d**12 - 2 * d**6 + 1.0_PREC)

   enddo
   !$omp end parallel do

   Ewca = e_wca

end subroutine energy_wca
