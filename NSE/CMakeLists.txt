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
option(NSE_ENABLE_OPENMP "Enable OpenMP support" ON)
set(
  NSE_MPI_PROVIDER "AUTO"
  CACHE STRING
  "MPI provider: AUTO, MSMPI, or SYSTEM"
)
set_property(CACHE NSE_MPI_PROVIDER PROPERTY STRINGS AUTO MSMPI SYSTEM)

nse_configure_mpi(nse_mpi)
if(NSE_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

set(NSE_MODULE_DIR "${CMAKE_CURRENT_BINARY_DIR}/mod")

set(
  NSE_DEFAULT_CORE_SOURCES
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
set(
  NSE_CORE_SOURCES "${NSE_DEFAULT_CORE_SOURCES}"
  CACHE STRING
  "Semicolon-separated NSE library source list"
)
set(
  NSE_MAIN_SOURCE "src/main/main_nse.f90"
  CACHE FILEPATH
  "Fortran source containing the NSE main program"
)
set(NSE_OUTPUT_NAME "solver" CACHE STRING "NSE executable file name without suffix")

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(NSE_DEFAULT_COMMON_FLAGS -fallow-invalid-boz -fallow-argument-mismatch)
  set(NSE_DEFAULT_WARNING_FLAGS -Wall -Wextra)
  set(NSE_DEFAULT_DEBUG_FLAGS -O0 -g -fcheck=all -fbacktrace)
  set(NSE_DEFAULT_RELEASE_FLAGS -O2)
  set(NSE_DEFAULT_RELWITHDEBINFO_FLAGS -O2 -g)
else()
  set(NSE_DEFAULT_COMMON_FLAGS "")
  set(NSE_DEFAULT_WARNING_FLAGS "")
  set(NSE_DEFAULT_DEBUG_FLAGS "")
  set(NSE_DEFAULT_RELEASE_FLAGS "")
  set(NSE_DEFAULT_RELWITHDEBINFO_FLAGS "")
endif()

set(NSE_COMMON_FLAGS "${NSE_DEFAULT_COMMON_FLAGS}" CACHE STRING "Common Fortran flags")
set(NSE_WARNING_FLAGS "${NSE_DEFAULT_WARNING_FLAGS}" CACHE STRING "Fortran warning flags")
set(NSE_DEBUG_FLAGS "${NSE_DEFAULT_DEBUG_FLAGS}" CACHE STRING "Debug Fortran flags")
set(NSE_RELEASE_FLAGS "${NSE_DEFAULT_RELEASE_FLAGS}" CACHE STRING "Release Fortran flags")
set(
  NSE_RELWITHDEBINFO_FLAGS "${NSE_DEFAULT_RELWITHDEBINFO_FLAGS}"
  CACHE STRING
  "RelWithDebInfo Fortran flags"
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
)
if(NSE_ENABLE_OPENMP)
  target_link_libraries(nse_core PUBLIC OpenMP::OpenMP_Fortran)
endif()

add_executable(nse_solver "${NSE_MAIN_SOURCE}")
set_target_properties(
  nse_solver
  PROPERTIES
    OUTPUT_NAME "${NSE_OUTPUT_NAME}"
    Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
)
target_link_libraries(nse_solver PRIVATE NSE::core)

function(nse_apply_fortran_options target_name)
  foreach(_flag IN LISTS NSE_COMMON_FLAGS)
    target_compile_options(${target_name} PRIVATE "${_flag}")
  endforeach()
  if(NSE_ENABLE_WARNINGS)
    foreach(_flag IN LISTS NSE_WARNING_FLAGS)
      target_compile_options(${target_name} PRIVATE "${_flag}")
    endforeach()
  endif()
  foreach(_flag IN LISTS NSE_DEBUG_FLAGS)
    target_compile_options(${target_name} PRIVATE "$<$<CONFIG:Debug>:${_flag}>")
  endforeach()
  foreach(_flag IN LISTS NSE_RELEASE_FLAGS)
    target_compile_options(${target_name} PRIVATE "$<$<CONFIG:Release>:${_flag}>")
  endforeach()
  foreach(_flag IN LISTS NSE_RELWITHDEBINFO_FLAGS)
    target_compile_options(${target_name} PRIVATE "$<$<CONFIG:RelWithDebInfo>:${_flag}>")
  endforeach()
endfunction()

nse_apply_fortran_options(nse_core)
nse_apply_fortran_options(nse_solver)

install(
  TARGETS nse_solver
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)

message(STATUS "NSE MPI provider: ${NSE_MPI_PROVIDER_RESOLVED}")
message(STATUS "NSE Fortran modules: ${NSE_MODULE_DIR}")
message(STATUS "NSE executable: ${CMAKE_CURRENT_BINARY_DIR}/bin/${NSE_OUTPUT_NAME}${CMAKE_EXECUTABLE_SUFFIX}")
