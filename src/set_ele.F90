subroutine set_ele(irep, tempk, ionstr, out_lb, out_Zp)

   use const
   use const_phys, only : PI, EPS0, BOLTZ_J, N_AVO, ELE, JOUL2KCAL_MOL
   use var_state, only : lambdaD, diele, length_per_charge, &
                         temp_independent, diele_dTcoef, tempK_ref
   use var_potential, only : ele_coef

   implicit none

   integer, intent(in) :: irep
   real(PREC), intent(in) :: tempk
   real(PREC), intent(in) :: ionstr
   real(PREC), intent(out), optional :: out_lb
   real(PREC), intent(out), optional :: out_Zp

   real(PREC) :: Tc, lb, Zp
   real(PREC), parameter ::  MM_A=87.740e0_PREC, MM_B=-0.40008e0_PREC
   real(PREC), parameter ::  MM_C=9.398e-4_PREC, MM_D=-1.410e-6_PREC

#ifdef _HTN_CONSISTENT
   real(PREC) :: temp_kT, lb_kT, kappaD, rho
#endif

   ! -----------------------------------------------------------------------
   ! Dielectric constant
   ! Temperature dependent (Malmberg and Maryott, 1956)

#ifdef _HTN_CONSISTENT
   !To be consistent with Hung's code
   temp_kT = tempk * BOLTZ_J * JOUL2KCAL_MOL
   diele(irep) = 296.0736276 - 619.2813716 * temp_kT + 531.2826741 * temp_kT**2 - 180.0369914 * temp_kT**3

   ! Bjerrum length
   lb_kT = 332.0637 / diele(irep)
   lb = lb_kT / temp_kT
   Zp = length_per_charge / lb

   if (temp_independent /= 0) then
      print '(a)', 'Error: temp_independnet /= 0 not implemented in _HTN_CONSISTENT (set_ele.F90)'
      flush(6)
      error stop
   endif

#else
   ! Default
   if (temp_independent == 0) then
      Tc = tempk - 273.15_PREC
      diele(irep) =  MM_A + MM_B*Tc + MM_C*Tc*Tc + MM_D*Tc*Tc*Tc

   else
      Tc = tempK_ref - 273.15_PREC
      diele(irep) =  MM_A + MM_B*Tc + MM_C*Tc*Tc + MM_D*Tc*Tc*Tc
      diele_dTcoef(irep) = 1.0_PREC + tempK_ref / diele(irep) &
                          * (MM_B + 2.0_PREC*MM_C*Tc + 3.0_PREC*MM_D*Tc*Tc)
   endif

   ! Bjerrum length
   lb = ELE * ELE / (4.0e0_PREC * PI * EPS0 * diele(irep) * BOLTZ_J * tempk) * 1.0e10_PREC

   ! Calculate phosphate charge taking into account monovalent salt condensation
   !xi = lb / length_per_charge
   !theta = 1.0 - 1.0/xi
   Zp = length_per_charge / lb
#endif

   !charge(:,irep) = 0.0_PREC

   ! Reflect Zp
   !do imp = 1, nmp
   !   if (.not. has_charge(imp)) cycle
   !   charge(imp, irep) = -Zp
   !enddo

   ! ----------------------------------------------------------------------
   ! coef: j_kcal * eq**2 / (4.0e0_PREC * PI * e0 * ek * rij)
   !   =  332.063713019 / ek
#ifdef _HTN_CONSISTENT
   ele_coef(irep) = lb_kT
#else
   ele_coef(irep) = JOUL2KCAL_MOL * 1.0e10_PREC * ELE**2 / (4.0e0_PREC * PI * EPS0 * diele(irep))
#endif

   ! Currently, all charges are phosphate with the same value of charge (-Zp).
   ! Therefore it can be included in coef.
   ele_coef(irep) = ele_coef(irep) * Zp**2

   ! Kd: sqrt(e0 * ek * RT / 2 * NA**2 * eq**2 * I)
   !   = sqrt(1.0e-3 * e0 * ek * kb / 2 * NA * eq**2) * sqrt(T(K) / I(M))

#ifdef _HTN_CONSISTENT
   rho = 2 * ionstr * 6.022e-4_PREC
   kappaD = sqrt(4 * 3.14159 * lb_kT * rho / temp_kT)
   lambdaD(irep) = 1.0_PREC / kappaD
#else
   lambdaD(irep) = 1.0e10_PREC * sqrt( (1.0e-3_PREC * EPS0 * diele(irep) * BOLTZ_J) &
                                  / (2.0_PREC * N_AVO * ELE**2)  )     &
                                * sqrt(tempk / ionstr)
#endif

   if (present(out_lb)) out_lb = lb
   if (present(out_Zp)) out_Zp = Zp

end subroutine set_ele
