cmake_minimum_required(VERSION 3.20)
project(ReactingFlow LANGUAGES Fortran)
include(CTest)
set(RF_DVODE_SOURCE_DIR "" CACHE PATH "Offline DVODE source checkout (fixed revision documented in FORTRAN.md)")
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
  src/fortran/mod_rf_units.f90)
target_link_libraries(reactingflow PUBLIC rf_dvode)
target_include_directories(reactingflow PUBLIC "$<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>")
add_executable(rf_reactor src/fortran/main_rf_reactor.f90)
target_link_libraries(rf_reactor PRIVATE reactingflow)
add_executable(rf_probe tests/fortran/rf_probe.f90)
target_link_libraries(rf_probe PRIVATE reactingflow)
if(BUILD_TESTING)
  add_executable(rf_unit tests/fortran/rf_unit.f90)
  target_link_libraries(rf_unit PRIVATE reactingflow)
  add_test(NAME rf_unit COMMAND rf_unit)
endif()
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_options(reactingflow PRIVATE -Wall -Wextra -ffree-line-length-none
    $<$<CONFIG:Debug>:-fcheck=all> $<$<CONFIG:Debug>:-fbacktrace>)
endif()
