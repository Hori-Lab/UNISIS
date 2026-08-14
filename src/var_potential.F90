module var_potential

   use, intrinsic :: ISO_FORTRAN_ENV, only: INT64
   use const, only : PREC
   use const_idx, only : BPT

   implicit none
  
   integer, save :: nbond
   integer, allocatable, save :: bond_mp(:,:)
   real(PREC), allocatable, save :: bond_para(:,:)
   real(PREC), save :: bond_k  ! = 15.0_PREC
   real(PREC), save :: bond_r0  ! = 5.9_PREC

   logical, save :: flg_angl_ReB
   integer, save :: nangl
   integer, allocatable, save :: angl_mp(:,:)
   real(PREC), allocatable, save :: angl_para(:,:)
   real(PREC), save :: angl_k  ! = 10.0_PREC
   real(PREC), save :: angl_t0  ! = 2.618_PREC
   
   integer, save :: ndih
   logical, save :: flg_dih_cos
   logical, save :: flg_dih_exp
   integer, allocatable, save :: dih_mp(:,:)
   real(PREC), allocatable, save :: dih_para(:,:)
   real(PREC), save :: dih_k
   real(PREC), save :: dih_w
   real(PREC), save :: dih_p0

   ! Nearest neighbour parameters
   real(PREC), allocatable :: NN_dG(:)  ! (1:NNT%MAX) if bp_model == 4
   real(PREC), allocatable :: NN_dH(:)  ! (1:NNT%MAX) if bp_model == 5
   real(PREC), allocatable :: NN_dS(:)  ! (1:NNT%MAX) if bp_model == 5
   logical, save :: flg_NNend
   real(PREC), allocatable :: NNend_dH(:)  ! (1:6) if flg_NNend
   real(PREC), allocatable :: NNend_dS(:)  ! (1:6) if flg_NNend

   ! Basepair
   integer,    save :: max_bp_per_nt
   integer,    save :: bp_min_loop  ! = 4 (in the original CAG work), = 3 (for mRNA)
   integer,    save :: bp_model

   !integer, allocatable :: bp_map_0(:,:)
   integer, allocatable :: bp_map(:,:)
   integer, allocatable :: bp3_map(:,:)
   !real(PREC), allocatable :: bp_map_dG(:,:,:)   ! (nmp, nmp, nrep_proc)
   real(PREC) :: bp_cutoff_energy  ! BP cutoff energy value in kcal/mol, e.g. -0.001.
                                   ! In the MC procedure, a base pair is
                                   ! regarded as formed only when its energy is below this value.
                                   ! If positive, no cutoff applies.
   real(PREC) :: bp_cutoff_dist
   integer,    save :: bp_seqdep
      ! = 0 (Default): No sequence dependence. Only U0_GC, U0_AU, U0_GU are required.
      ! = 1: Sequence dependent parameters. All possible combinations of trinucleotide-dimer are required.

   real(PREC), save :: coef_dG
   real(PREC), save :: dH0
   real(PREC), save :: dS0
   real(PREC), save :: dHend0
   real(PREC), save :: dSend0

   real(PREC) :: bp3_dH(1:468)
   real(PREC) :: bp3_dS(1:468)

   type basepair_parameters
      real(PREC) :: cutoff_ddist
      real(PREC) :: U0
      real(PREC) :: bond_k
      real(PREC) :: bond_r
      real(PREC) :: angl_k1
      real(PREC) :: angl_k2
      real(PREC) :: angl_k3
      real(PREC) :: angl_k4
      real(PREC) :: angl_theta1
      real(PREC) :: angl_theta2
      real(PREC) :: angl_theta3
      real(PREC) :: angl_theta4
      real(PREC) :: dihd_k1
      real(PREC) :: dihd_k2
      real(PREC) :: dihd_phi1
      real(PREC) :: dihd_phi2
   end type basepair_parameters

   type(basepair_parameters), save :: bp_paras(BPT%MAX)

   ! Basepair list
   integer, allocatable :: nbp(:)
   integer :: nbp_max  ! This defines the size of neighbor list
   integer, allocatable :: bp_mp(:,:,:)  ! (3, nbp, nrep_proc) 1: imp1, 2: imp2, 3: bp type
   real(PREC), allocatable :: bp_coef(:,:,: )  ! (2, nbp, nrep_proc) 1: dH, 2: dS

   ! WCA parameters
   real(PREC), save :: wca_sigma  ! = 10.0
   real(PREC), save :: wca_eps    ! = 2.0

   ! WCA list
   integer, allocatable :: nwca(:)
   integer, save :: nwca_max  ! This defines the size of neighbor list
   integer, allocatable, save :: wca_mp(:,:,:)  ! (1:2, nwca, nrep_proc) or (1:2, nwca_max, nrep_proc),
                                                ! 1:imp1, 2:imp2

   ! Electrostatic parameters
   logical, save :: flg_ele
   logical, save :: ele_exclude_covalent_bond_pairs
   logical, save :: ele_exclude_covalent_angle_pairs
   integer, save :: ele_cutoff_type
   real(PREC), save :: ele_cutoff_inp
   real(PREC), allocatable :: ele_cutoff(:)  ! (nrep_proc)
   real(PREC), allocatable :: ele_coef(:)  ! (nrep_proc)

   ! Electrostatic list
   integer, allocatable :: nele(:)
   integer, save :: nele_max
   integer, allocatable, save :: ele_mp(:, :, :)  ! (2, nele_max, nrep_proc)

   ! Stage potential parameters
   logical, save :: flg_stage
   real(PREC), save :: stage_sigma ! = wca_sigma = 10.0 Angstrom for now
   real(PREC), save :: stage_eps ! = 1.2 kcal/mol

   ! Tweezers
   logical, save :: flg_twz
   integer, save :: ntwz_DCF  !!! Dual Constant Force
   integer, allocatable :: twz_DCF_pairs(:,:)  ! (2, ntwz_DCF), pair IDs imp1 and imp2
   real(PREC), allocatable :: twz_DCF_direction(:,:)  ! (3, ntwz_DCF), force vectors given in input
                                                      ! In f-REMD, this will be normalised.
                                                      ! If non f-REMD, this will be copied to twz_DCF_forces.
   real(PREC), allocatable :: twz_DCF_forces(:,:,:)  ! (3, ntwz_DCF, nrep_proc), force vectors (magniture x direction)
   integer, save :: ntwz_FR  !!! Force Ramp
   integer, allocatable :: twz_FR_pairs(:,:)   ! (2, ntwz_FR)
   real(PREC), allocatable :: twz_FR_k(:,:)  ! (2, ntwz_FR)
   real(PREC), allocatable :: twz_FR_speed(:,:)  ! (2, ntwz_FR)
   real(PREC), allocatable :: twz_FR_init(:,:,:)  ! (3, 2, ntwz_FR)
   real(PREC), allocatable :: twz_FR_velo(:,:,:)  ! (3, 2, ntwz_FR)


   ! Bias_SS
   logical, save :: flg_bias_ss
   real(PREC), save :: bias_ss_force

   ! Bias_Rg
   logical, save :: flg_bias_rg
   integer, save :: bias_rg_pott  ! potential type (POTT%HARMONIC or POTT%FLATBOTTOMED)
   real(PREC), save :: bias_rg_k
   real(PREC), save :: bias_rg_0_inp        ! Value from the input file
   real(PREC), allocatable :: bias_rg_0(:)  ! (1:nrep_proc)
                                            ! bias_rg_0 = bias_rg_0_inp if not Rg-REMD

   ! Time-dependent Bias-Rg
   logical, save :: flg_timed_bias_rg
   integer, save :: ntimed_bias_rg
   integer(INT64), allocatable, save :: timed_bias_rg_step(:)
   real(PREC), allocatable, save :: timed_bias_rg_k(:)
   real(PREC), allocatable, save :: timed_bias_rg_0(:)
   integer :: itimed_bias_rg
   integer(INT64) :: istep_timed_bias_rg_next

   ! Restraint
   logical, save :: flg_restraint
   logical, save :: flg_rest_sigb
   integer, save :: nrest_sigb
   integer, allocatable :: rest_sigb_id(:,:)        ! (1:3, nrest_sig) 1=target ID, 2=reference ID, 3=type
                                                    ! type = 1, Sigmoid; 2, Sigmoid-to-bead
   real(PREC), allocatable :: rest_sigb_rcut(:)     ! (nrest_sig) Cutoff distance, r_cut
   real(PREC), allocatable :: rest_sigb_para(:, :)  ! (1:3, nrest_sig) 1=epsilon, 2=r_bound, 3=delta

end module var_potential
