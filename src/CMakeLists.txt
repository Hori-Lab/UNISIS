option(BP_HALT_IEEE_EXCEPTIONS "Use ieee halting mode functions in write_bp" OFF)
option(ALLOW_ARGUMENT_MISMATCH "Use -fallow-argument-mismatch option" OFF)

set(UNISIS_COMMON_SOURCES
    const.F90 const_idx.F90 const_phys.F90
    rnd_boxmuller.F90
    dcd.F90
    pbc.F90 progress.F90
    var_parallel.F90 var_replica.F90
    var_top.F90 var_state.F90 var_potential.F90 var_io.F90
    read_force_field.F90 read_input.F90 read_pdb.F90 read_xyz.F90 read_fasta.F90
    read_anneal.F90 read_ss.F90 read_timed_bias_rg.F90 read_restraint.F90
    write_out.F90
    write_rep_table.F90
    write_rst.F90 read_rst.F90
    write_bp.F90
    write_bpcoef.F90
    init_mts.F90 init_out_files.F90 init_replica.F90 init_ele.F90 init_freeze.F90
    init_tweezers.F90 init_bias_rg.F90
    init_sequence.F90 init_structure.F90 gen_random_coil.F90
    list_local.F90 list_bp.F90
    list_exv.F90 list_ele.F90
    init_bp.F90 set_ele.F90
    replica_exchange.F90
    energy_replica.F90 energy_sumup.F90
    energy_kinetic.F90
    energy_bond.F90 energy_angl.F90 energy_angl_ReB.F90
    energy_bp.F90 energy_bp_limit.F90 energy_bp_triplet.F90 energy_bp_limit_triplet.F90
    energy_dih_cos.F90 energy_dih_exp.F90
    energy_wca.F90 energy_ele_DH.F90 energy_stage.F90 energy_tweezers.F90
    energy_rg.F90 energy_restraint.F90
    force.F90
    force_bond.F90 force_angl.F90 force_angl_ReB.F90
    force_bp.F90 force_bp_limit.F90 force_bp_triplet.F90 force_bp_limit_triplet.F90
    force_dih_cos.F90 force_dih_exp.F90
    force_wca.F90 force_ele_DH.F90 force_stage.F90 force_tweezers.F90
    force_rg.F90 force_restraint.F90
    neighbor_list.F90
    job_md.F90 job_dcd.F90 job_check_force.F90
    sis_abort.F90
    warning_stage.F90
    main.F90
    ${CMAKE_BINARY_DIR}/generated/unisis_version.F90
)

function(configure_unisis_target target_name output_name)
    set_target_properties(${target_name} PROPERTIES
        OUTPUT_NAME ${output_name}
        Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules/${target_name}"
    )

    target_include_directories(${target_name} PRIVATE
        "${CMAKE_BINARY_DIR}/generated"
        "${CMAKE_BINARY_DIR}/modules/${target_name}"
    )

    if (BP_HALT_IEEE_EXCEPTIONS)
        target_compile_definitions(${target_name} PRIVATE BP_HALT_IEEE_EXCEPTIONS)
    endif()

    if (ALLOW_ARGUMENT_MISMATCH)
        target_compile_options(${target_name} PRIVATE -fallow-argument-mismatch)
    endif()

    if (DUMPFORCE)
        target_compile_definitions(${target_name} PRIVATE DUMPFORCE)
    endif()

    if (OpenMP_Fortran_FOUND)
        target_link_libraries(${target_name} PRIVATE OpenMP::OpenMP_Fortran)
    endif()

    target_link_libraries(${target_name} PRIVATE
        toml-f-lib
        githash
        mt_stream-lib
    )
endfunction()

if(UNISIS_BUILD_SERIAL)
    add_executable(unisis_serial ${UNISIS_COMMON_SOURCES})
    configure_unisis_target(unisis_serial unisis)
    install(TARGETS unisis_serial RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

if(UNISIS_BUILD_MPI)
    add_executable(unisis_mpi ${UNISIS_COMMON_SOURCES})
    configure_unisis_target(unisis_mpi unisis_mpi)
    target_compile_definitions(unisis_mpi PRIVATE PAR_MPI)
    target_link_libraries(unisis_mpi PRIVATE MPI::MPI_Fortran)
    install(TARGETS unisis_mpi RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

add_executable(show_rst
    const.F90 const_idx.F90
    show_rst.F90
)
add_dependencies(show_rst mt_stream-lib)
target_link_libraries(show_rst mt_stream-lib)
set_target_properties(show_rst PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules/show_rst"
)
