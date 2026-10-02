cmake_minimum_required(VERSION 3.20)
project(ReactingFlow LANGUAGES Fortran)
include(CTest)
set(RF_LEGACY_NSE_SOURCE_DIR "" CACHE PATH "Optional NSE source root for migration comparison tests only")
set(RF_DVODE_SOURCE_DIR "" CACHE PATH "Offline DVODE source checkout (fixed revision documented in REACTORS.md)")
if(NOT RF_DVODE_SOURCE_DIR)
  include(FetchContent)
  FetchContent_Declare(rf_dvode
    GIT_REPOSITORY https://github.com/jacobwilliams/dvode.git
    GIT_TAG 59d2eedd1410a4bd6c196c24d08f95a7c20565c0)
  FetchContent_MakeAvailable(rf_dvode)
  set(RF_DVODE_SOURCE_DIR "${rf_dvode_SOURCE_DIR}")
endif()
configure_file("${RF_DVODE_SOURCE_DIR}/LICENSE.md" "${CMAKE_CURRENT_BINARY_DIR}/DVODE_LICENSE.md" COPYONLY)
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")
file(MAKE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}")
add_library(rf_dvode STATIC
  ${RF_DVODE_SOURCE_DIR}/src/dvode_kinds_module.F90
  ${RF_DVODE_SOURCE_DIR}/src/dvode_blas_module.F90
  ${RF_DVODE_SOURCE_DIR}/src/dvode_linpack_module.F90
  ${RF_DVODE_SOURCE_DIR}/src/dvode_module.F90)
target_compile_definitions(rf_dvode PRIVATE REAL64)
add_library(reactingflow STATIC src/fortran/mod_rf_mechanism.f90
  src/fortran/mod_rf_thermo.f90 src/fortran/mod_rf_kinetics.f90 src/fortran/mod_rf_reactor.f90
  src/fortran/mod_rf_units.f90 src/fortran/mod_rf_transport.f90 src/fortran/mod_rf_flow1d.f90
  src/fortran/mod_rf_shock.f90 src/fortran/mod_rf_znd.f90
  src/fortran/mod_rf_equilibrium.f90 src/fortran/mod_rf_cj.f90 src/fortran/mod_rf_profile.f90
  src/fortran/mod_rf_wave.f90 src/fortran/mod_rf_finite_volume.f90 src/fortran/mod_rf_grid2d.f90
  src/fortran/mod_rf_flow2d.f90)
target_link_libraries(reactingflow PUBLIC rf_dvode)
target_include_directories(reactingflow PUBLIC "$<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>")
add_executable(rf_reactor src/fortran/main_rf_reactor.f90)
target_link_libraries(rf_reactor PRIVATE reactingflow)
add_executable(rf_flow1d src/fortran/main_rf_flow1d.f90)
target_link_libraries(rf_flow1d PRIVATE reactingflow)
add_executable(rf_mesh2d src/fortran/main_rf_mesh2d.f90)
target_link_libraries(rf_mesh2d PRIVATE reactingflow)
add_executable(rf_flow2d src/fortran/main_rf_flow2d.f90)
target_link_libraries(rf_flow2d PRIVATE reactingflow)
add_executable(rf_flow_balance src/fortran/main_rf_flow_balance.f90)
target_link_libraries(rf_flow_balance PRIVATE reactingflow)
add_executable(rf_normal_shock src/fortran/main_rf_shock.f90)
target_link_libraries(rf_normal_shock PRIVATE reactingflow)
add_executable(rf_znd src/fortran/main_rf_znd.f90)
target_link_libraries(rf_znd PRIVATE reactingflow)
add_executable(rf_cj src/fortran/main_rf_cj.f90)
target_link_libraries(rf_cj PRIVATE reactingflow)
add_executable(rf_probe tests/fortran/rf_probe.f90)
target_link_libraries(rf_probe PRIVATE reactingflow)
add_executable(rf_equilibrium_probe tests/fortran/rf_equilibrium_probe.f90)
target_link_libraries(rf_equilibrium_probe PRIVATE reactingflow)
add_executable(rf_transport_probe tests/fortran/rf_transport_probe.f90)
target_link_libraries(rf_transport_probe PRIVATE reactingflow)
if(BUILD_TESTING)
  add_executable(rf_flow2d_unit tests/fortran/rf_flow2d_unit.f90)
  target_link_libraries(rf_flow2d_unit PRIVATE reactingflow)
  add_test(NAME rf_flow2d_unit COMMAND rf_flow2d_unit)
  add_executable(rf_grid2d_unit tests/fortran/rf_grid2d_unit.f90)
  target_link_libraries(rf_grid2d_unit PRIVATE reactingflow)
  add_test(NAME rf_grid2d_unit COMMAND rf_grid2d_unit)
  foreach(invalid_grid reversed closed concave)
    add_test(NAME rf_grid2d_reject_${invalid_grid} COMMAND rf_grid2d_unit ${invalid_grid})
    set_tests_properties(rf_grid2d_reject_${invalid_grid} PROPERTIES WILL_FAIL TRUE)
  endforeach()
  add_executable(rf_finite_volume_unit tests/fortran/rf_finite_volume_unit.f90)
  target_link_libraries(rf_finite_volume_unit PRIVATE reactingflow)
  add_test(NAME rf_finite_volume_unit COMMAND rf_finite_volume_unit)
  add_executable(rf_wave_unit tests/fortran/rf_wave_unit.f90)
  target_link_libraries(rf_wave_unit PRIVATE reactingflow)
  add_test(NAME rf_wave_unit COMMAND rf_wave_unit)
  add_executable(rf_unit tests/fortran/rf_unit.f90)
  target_link_libraries(rf_unit PRIVATE reactingflow)
  add_test(NAME rf_unit COMMAND rf_unit)
  # Keep scratch I/O inside this build tree, independent of system temp ACLs.
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Testing/tmp")
  set_tests_properties(rf_unit PROPERTIES ENVIRONMENT
    "TMPDIR=${CMAKE_CURRENT_BINARY_DIR}/Testing/tmp;TMP=${CMAKE_CURRENT_BINARY_DIR}/Testing/tmp;TEMP=${CMAKE_CURRENT_BINARY_DIR}/Testing/tmp")
  add_executable(rf_flow_unit tests/fortran/rf_flow_unit.f90)
  target_link_libraries(rf_flow_unit PRIVATE reactingflow)
  add_test(NAME rf_flow_unit COMMAND rf_flow_unit)
  add_executable(rf_transport_unit tests/fortran/rf_transport_unit.f90)
  target_link_libraries(rf_transport_unit PRIVATE reactingflow)
  add_test(NAME rf_transport_unit COMMAND rf_transport_unit)
  add_executable(rf_muscl_unit tests/fortran/rf_muscl_unit.f90)
  target_link_libraries(rf_muscl_unit PRIVATE reactingflow)
  add_test(NAME rf_muscl_unit COMMAND rf_muscl_unit)
  if(RF_LEGACY_NSE_SOURCE_DIR)
    set(legacy_mc "${RF_LEGACY_NSE_SOURCE_DIR}/src/extensions/multicomponent")
    add_library(rf_legacy_reference STATIC
      ${RF_LEGACY_NSE_SOURCE_DIR}/src/common/mod_precision.f90
      ${legacy_mc}/mod_mc_config.f90
      ${legacy_mc}/mod_mc_state_layout.f90
      ${legacy_mc}/mod_mc_euler_config.f90
      ${legacy_mc}/mod_mc_geometry.f90
      ${legacy_mc}/providers/mod_mc_thermodynamics_thermally_perfect.f90
      ${legacy_mc}/providers/mod_mc_chemistry_one_step_arrhenius.f90
      ${legacy_mc}/mod_mc_chemistry_integrator.f90
      ${legacy_mc}/mod_mc_euler_field.f90
      ${legacy_mc}/mod_mc_boundary.f90
      ${legacy_mc}/providers/mod_mc_mapped_transport_none.f90
      ${legacy_mc}/mod_mc_mapped_flux.f90
      ${legacy_mc}/mod_mc_euler_flux.f90)
    target_include_directories(rf_legacy_reference PUBLIC "${CMAKE_Fortran_MODULE_DIRECTORY}")
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(rf_legacy_reference PRIVATE -ffree-line-length-none $<$<CONFIG:Debug>:-fcheck=all>)
    endif()
    add_executable(rf_legacy_unit tests/fortran/rf_legacy_unit.f90)
    target_link_libraries(rf_legacy_unit PRIVATE reactingflow rf_legacy_reference)
    add_test(NAME rf_legacy_unit COMMAND rf_legacy_unit "${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/legacy_one_step.in")
  endif()
endif()
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_options(reactingflow PRIVATE -Wall -Wextra -ffree-line-length-none
    $<$<CONFIG:Debug>:-fcheck=all> $<$<CONFIG:Debug>:-fbacktrace>)
endif()
