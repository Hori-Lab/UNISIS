cmake_minimum_required(VERSION 3.15)

project(
    UNISIS
    VERSION 1.0.0
    LANGUAGES Fortran
)

include(GNUInstallDirs)

set(UNISIS_RELEASE_LABEL "2026.08" CACHE STRING "User-facing UNISIS release label")

option(UNISIS_BUILD_SERIAL "Build the serial UNISIS executable" ON)
option(UNISIS_BUILD_MPI "Build the MPI UNISIS executable" OFF)
option(BUILD_MPI "Deprecated: use -DUNISIS_BUILD_MPI=ON instead" OFF)

if(BUILD_MPI)
    message(WARNING "BUILD_MPI is deprecated; use UNISIS_BUILD_MPI. Building the MPI variant only.")
    set(UNISIS_BUILD_SERIAL OFF CACHE BOOL "Build the serial UNISIS executable" FORCE)
    set(UNISIS_BUILD_MPI ON CACHE BOOL "Build the MPI UNISIS executable" FORCE)
endif()

if(NOT UNISIS_BUILD_SERIAL AND NOT UNISIS_BUILD_MPI)
    message(FATAL_ERROR "At least one of UNISIS_BUILD_SERIAL or UNISIS_BUILD_MPI must be enabled.")
endif()

#ifort
if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    #set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -traceback -check all -warn all -fp-model strict -standf18")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -traceback -check all -warn all -fp-model strict")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native")
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O3 -g -march=native")
endif()

#gfortran
if (CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    #set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -Wall -g -fbacktrace -ffpe-trap=zero,overflow,denormal,underflow -ffree-line-length-none -fbounds-check -fmax-errors=5")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -Wall -g -fbacktrace -ffpe-trap=zero,overflow,denormal -ffree-line-length-none -fbounds-check -fmax-errors=5")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -ffree-line-length-none")
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O3 -g -march=native -ffree-line-length-none")
endif()

#MPI
if (UNISIS_BUILD_MPI)
    find_package(MPI REQUIRED COMPONENTS Fortran)
    if(NOT MPI_Fortran_FOUND)
        message(FATAL_ERROR "Could not find Fortran MPI.  Please set MPI_Fortran_COMPILER to point to the mpifort wrapper.")
    endif()
endif (UNISIS_BUILD_MPI)

# Option to dumping force data
option(DUMPFORCE "Use DUMPFORCE option to output force data" OFF)

find_package(OpenMP)


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
foreach(config DEBUG RELEASE RELWITHDEBINFO MINSIZEREL)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${config} "${CMAKE_BINARY_DIR}/bin")
endforeach()

file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/generated")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/unisis_version.F90.in"
    "${CMAKE_BINARY_DIR}/generated/unisis_version.F90"
    @ONLY
)


# To output Git commit hash. Either CheckGitSetup or NoGitSetup will be called.
include(./CheckGit.cmake)

find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")

    # Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()

    # Output Git hash
    CheckGitSetup()

else()

    # Set githash = "?" and output version number (in main.F90)
    NoGitSetup()

endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/subprojects/toml-f/CMakeLists.txt")
    message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()

add_subdirectory(subprojects)
add_subdirectory(src)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_BINARY_DIR}/cmake_uninstall.cmake"
    @ONLY
)

add_custom_target(
    uninstall
    COMMAND "${CMAKE_COMMAND}" -P "${CMAKE_BINARY_DIR}/cmake_uninstall.cmake"
)
