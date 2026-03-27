module const_phys

   use const

   real(PREC), parameter :: PI = 4.D0 * atan(1.0_PREC)
   real(PREC), parameter :: EPS0 = 8.8541878128e-12_PREC    ! Electric constant [F/m]
   real(PREC), parameter :: ELE  = 1.602176634e-19_PREC     ! Elementary charge [C]
   real(PREC), parameter :: BOLTZ_J = 1.380649e-23_PREC     ! Boltzmann constant [J/K]
   real(PREC), parameter :: N_AVO   = 6.02214076e23_PREC    ! Avogadro constant [/mol]
   real(PREC), parameter :: KCAL2JOUL = 4184.0_PREC         ! (kcal -> J)  [J/kcal]

   real(PREC), parameter :: JOUL2KCAL = 1.0_PREC / KCAL2JOUL   !< (J -> kcal)  [kcal/J]
   real(PREC), parameter :: JOUL2KCAL_MOL  = JOUL2KCAL * N_AVO  !< (J -> kcal/mol)
   real(PREC), parameter :: BOLTZ_KCAL_MOL = BOLTZ_J * JOUL2KCAL_MOL   !< Boltzmann constant [kcal/mol/K]
   real(PREC), parameter :: SMALL_VALUE = epsilon(PI)
   real(PREC), parameter :: INVALID_JUDGE = 1.0e30_PREC
   real(PREC), parameter :: INVALID_VALUE = 1.0e31_PREC
   real(PREC), parameter :: ZERO_JUDGE    = 1.0e-6_PREC

   integer, parameter :: INVALID_INT_VALUE = huge(M_INT)
   integer, parameter :: INVALID_INT_JUDGE = INVALID_INT_VALUE - 1

contains

   real(PREC) function J2kcalmol(j)
      real(PREC), intent(in) :: j
      J2kcalmol = j * JOUL2KCAL_MOL
   endfunction J2kcalmol

endmodule const_phys
