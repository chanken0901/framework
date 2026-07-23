cmake_minimum_required(VERSION 3.24)

# Windowsの長いソースパスでは、CMakeにオブジェクト名の短縮を早めに行わせる。
set(CMAKE_OBJECT_PATH_MAX 240)

project(
  NSE
  VERSION 1.0.0
  DESCRIPTION "MPI/OpenMP Navier-Stokes solver"
  LANGUAGES Fortran
)

include(GNUInstallDirs)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(NSEMPI)

option(NSE_ENABLE_WARNINGS "Enable additional Fortran compiler warnings" ON)
set(
  NSE_MPI_PROVIDER "AUTO"
  CACHE STRING
  "MPI provider: AUTO, MSMPI, or SYSTEM"
)
set_property(CACHE NSE_MPI_PROVIDER PROPERTY STRINGS AUTO MSMPI SYSTEM)

nse_configure_mpi(nse_mpi)
find_package(OpenMP REQUIRED COMPONENTS Fortran)

set(NSE_MODULE_DIR "${CMAKE_CURRENT_BINARY_DIR}/mod")

set(
  NSE_CORE_SOURCES
  src/common/mod_precision.f90
  src/common/mod_constants.f90
  src/common/mod_common_config.f90
  src/common/mod_model_config.f90
  src/grid/mod_grid_fvm.f90
  src/field/mod_field_nse.f90
  src/io/mod_input_reader.f90
  src/mpi/module_mpi_TW_NEC.f90
  src/io/mod_slf_output.f90
)

add_library(nse_core STATIC ${NSE_CORE_SOURCES})
add_library(NSE::core ALIAS nse_core)

set_target_properties(
  nse_core
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}"
    POSITION_INDEPENDENT_CODE ON
)

target_include_directories(
  nse_core
  PUBLIC
    "$<BUILD_INTERFACE:${NSE_MODULE_DIR}>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>"
)

target_link_libraries(
  nse_core
  PUBLIC
    nse_mpi
    OpenMP::OpenMP_Fortran
)

add_executable(nse_solver src/main/main_nse.f90)
set_target_properties(
  nse_solver
  PROPERTIES
    OUTPUT_NAME solver
    Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
)
target_link_libraries(nse_solver PRIVATE NSE::core)

function(nse_apply_fortran_options target_name)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(
      ${target_name}
      PRIVATE
        -fallow-invalid-boz
        -fallow-argument-mismatch
        "$<$<BOOL:${NSE_ENABLE_WARNINGS}>:-Wall>"
        "$<$<BOOL:${NSE_ENABLE_WARNINGS}>:-Wextra>"
        "$<$<CONFIG:Debug>:-O0>"
        "$<$<CONFIG:Debug>:-g>"
        "$<$<CONFIG:Debug>:-fcheck=all>"
        "$<$<CONFIG:Debug>:-fbacktrace>"
        "$<$<CONFIG:Release>:-O2>"
        "$<$<CONFIG:RelWithDebInfo>:-O2>"
        "$<$<CONFIG:RelWithDebInfo>:-g>"
    )
  endif()
endfunction()

nse_apply_fortran_options(nse_core)
nse_apply_fortran_options(nse_solver)

install(
  TARGETS nse_solver
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)

message(STATUS "NSE MPI provider: ${NSE_MPI_PROVIDER_RESOLVED}")
message(STATUS "NSE Fortran modules: ${NSE_MODULE_DIR}")
message(STATUS "NSE executable: ${CMAKE_CURRENT_BINARY_DIR}/bin/solver${CMAKE_EXECUTABLE_SUFFIX}")
