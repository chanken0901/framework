cmake_minimum_required(VERSION 3.24)

# Windowsの長いソースパスでは、CMakeにオブジェクト名の短縮を早めに行わせる。
set(CMAKE_OBJECT_PATH_MAX 240)

project(
  NSE
  VERSION 1.0.0
  DESCRIPTION "MPI/OpenMP and CUDA Navier-Stokes solver"
  LANGUAGES Fortran
)

include(GNUInstallDirs)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(NSEMPI)

option(NSE_ENABLE_WARNINGS "Enable additional Fortran compiler warnings" ON)
option(NSE_ENABLE_OPENMP "Enable OpenMP support" ON)
option(NSE_USE_MPI "Build the MPI domain-decomposition backend" ON)
option(
  NSE_ENABLE_CUFFTMP_CONTRACT_TESTS
  "Validate the future MPI+CUDA cuFFTMp forcing interface contract"
  ON
)
set(
  NSE_GPU_BACKEND "none"
  CACHE STRING
  "GPU backend: none or cuda"
)
set_property(CACHE NSE_GPU_BACKEND PROPERTY STRINGS none cuda)
set(
  NSE_CUDA_ARCHITECTURES ""
  CACHE STRING
  "Optional CUDA architectures, for example 86 for an RTX 3050"
)
set(
  NSE_MPI_PROVIDER "AUTO"
  CACHE STRING
  "MPI provider: AUTO, MSMPI, or SYSTEM"
)
set_property(CACHE NSE_MPI_PROVIDER PROPERTY STRINGS AUTO MSMPI SYSTEM)

set(
  NSE_CONVECTIVE_BACKEND "keep_family"
  CACHE STRING
  "Convective-flux backend compiled into the NSE solver"
)
set_property(CACHE NSE_CONVECTIVE_BACKEND PROPERTY STRINGS keep_family)
set(
  NSE_VISCOUS_SCHEME "central6"
  CACHE STRING
  "Viscous-term backend compiled into the NSE solver"
)
set_property(CACHE NSE_VISCOUS_SCHEME PROPERTY STRINGS none central6)
set(
  NSE_BOUNDARY_SCHEME "periodic"
  CACHE STRING
  "Boundary-condition backend compiled into the NSE solver"
)
set_property(CACHE NSE_BOUNDARY_SCHEME PROPERTY STRINGS periodic)
set(
  NSE_INIT_FFT_BACKEND "none"
  CACHE STRING
  "Initial-condition FFT backend: none or 2decomp_fftw"
)
set_property(CACHE NSE_INIT_FFT_BACKEND PROPERTY STRINGS none 2decomp_fftw)
set(
  NSE_FORCING_FFT_BACKEND "none"
  CACHE STRING
  "Forcing Helmholtz FFT backend: none, 2decomp_fftw, or cufft"
)
set_property(
  CACHE NSE_FORCING_FFT_BACKEND
  PROPERTY STRINGS none 2decomp_fftw cufft
)
set(
  NSE_2DECOMP_ROOT "$ENV{DECOMP2D_ROOT}"
  CACHE PATH
  "2DECOMP&FFT installation prefix containing lib/decomp2d"
)

string(TOLOWER "${NSE_CONVECTIVE_BACKEND}" NSE_CONVECTIVE_BACKEND)
string(TOLOWER "${NSE_VISCOUS_SCHEME}" NSE_VISCOUS_SCHEME)
string(TOLOWER "${NSE_BOUNDARY_SCHEME}" NSE_BOUNDARY_SCHEME)
string(TOLOWER "${NSE_GPU_BACKEND}" NSE_GPU_BACKEND)
string(TOLOWER "${NSE_INIT_FFT_BACKEND}" NSE_INIT_FFT_BACKEND)
string(TOLOWER "${NSE_FORCING_FFT_BACKEND}" NSE_FORCING_FFT_BACKEND)

if(NOT NSE_GPU_BACKEND STREQUAL "none" AND
   NOT NSE_GPU_BACKEND STREQUAL "cuda")
  message(FATAL_ERROR
    "Unsupported NSE_GPU_BACKEND='${NSE_GPU_BACKEND}'. "
    "Available backends: none, cuda"
  )
endif()
if(NSE_GPU_BACKEND STREQUAL "cuda" AND NSE_USE_MPI)
  message(FATAL_ERROR
    "The NSE CUDA backend is single-GPU only. Configure with NSE_USE_MPI=OFF."
  )
endif()
if(NOT NSE_INIT_FFT_BACKEND STREQUAL "none" AND
   NOT NSE_INIT_FFT_BACKEND STREQUAL "2decomp_fftw")
  message(FATAL_ERROR
    "Unsupported NSE_INIT_FFT_BACKEND='${NSE_INIT_FFT_BACKEND}'. "
    "Available backends: none, 2decomp_fftw"
  )
endif()
if(NSE_INIT_FFT_BACKEND STREQUAL "2decomp_fftw" AND NOT NSE_USE_MPI)
  message(FATAL_ERROR
    "The 2decomp_fftw initialization backend requires NSE_USE_MPI=ON."
  )
endif()
if(NSE_INIT_FFT_BACKEND STREQUAL "2decomp_fftw" AND
   NOT NSE_GPU_BACKEND STREQUAL "none")
  message(FATAL_ERROR
    "The 2decomp_fftw initialization backend is currently CPU/MPI only."
  )
endif()
if(NOT NSE_FORCING_FFT_BACKEND STREQUAL "none" AND
   NOT NSE_FORCING_FFT_BACKEND STREQUAL "2decomp_fftw" AND
   NOT NSE_FORCING_FFT_BACKEND STREQUAL "cufft")
  message(FATAL_ERROR
    "Unsupported NSE_FORCING_FFT_BACKEND='${NSE_FORCING_FFT_BACKEND}'. "
    "Available runtime backends: none, 2decomp_fftw, cufft. "
    "cuFFTMp is represented by a test-only interface contract."
  )
endif()
if(NSE_FORCING_FFT_BACKEND STREQUAL "2decomp_fftw" AND
   (NOT NSE_USE_MPI OR NOT NSE_GPU_BACKEND STREQUAL "none"))
  message(FATAL_ERROR
    "The 2decomp_fftw forcing backend requires CPU execution with MPI."
  )
endif()
if(NSE_FORCING_FFT_BACKEND STREQUAL "cufft" AND
   NOT NSE_GPU_BACKEND STREQUAL "cuda")
  message(FATAL_ERROR
    "The cufft forcing backend requires NSE_GPU_BACKEND=cuda."
  )
endif()

if(NSE_CONVECTIVE_BACKEND STREQUAL "keep_family")
  set(NSE_CONVECTIVE_SOURCE src/numerics/convective/mod_convective_keep.f90)
else()
  message(FATAL_ERROR
    "Unsupported NSE_CONVECTIVE_BACKEND='${NSE_CONVECTIVE_BACKEND}'. "
    "Available backends: keep_family"
  )
endif()

if(NSE_VISCOUS_SCHEME STREQUAL "none")
  set(NSE_VISCOUS_SOURCE src/numerics/viscous/mod_viscous_none.f90)
elseif(NSE_VISCOUS_SCHEME STREQUAL "central6")
  set(NSE_VISCOUS_SOURCE src/numerics/viscous/mod_viscous_central6.f90)
else()
  message(FATAL_ERROR
    "Unsupported NSE_VISCOUS_SCHEME='${NSE_VISCOUS_SCHEME}'. "
    "Available backends: none, central6"
  )
endif()

if(NSE_BOUNDARY_SCHEME STREQUAL "periodic")
  set(NSE_BOUNDARY_SOURCE src/boundary/mod_boundary_periodic.f90)
else()
  message(FATAL_ERROR
    "Unsupported NSE_BOUNDARY_SCHEME='${NSE_BOUNDARY_SCHEME}'. "
    "Available backends: periodic"
  )
endif()

if(NSE_USE_MPI)
  nse_configure_mpi(nse_mpi)
else()
  add_library(nse_mpi INTERFACE)
  set(NSE_MPI_PROVIDER_RESOLVED "DISABLED")
endif()
if(NSE_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

set(NSE_INIT_FFT_SOURCE src/init/mod_init_hit_spectral_stub.f90)
if(NSE_INIT_FFT_BACKEND STREQUAL "2decomp_fftw" OR
   NSE_FORCING_FFT_BACKEND STREQUAL "2decomp_fftw")
  if(NSE_2DECOMP_ROOT)
    find_package(
      decomp2d CONFIG REQUIRED
      PATHS "${NSE_2DECOMP_ROOT}"
      NO_DEFAULT_PATH
    )
  else()
    find_package(decomp2d CONFIG REQUIRED)
  endif()
  if(NSE_INIT_FFT_BACKEND STREQUAL "2decomp_fftw")
    set(NSE_INIT_FFT_SOURCE src/init/mod_init_hit_spectral_2decomp.f90)
  endif()
endif()

set(NSE_FORCING_SOURCE src/forcing/mod_nse_forcing_stub.f90)
if(NSE_FORCING_FFT_BACKEND STREQUAL "2decomp_fftw")
  set(NSE_FORCING_SOURCE src/forcing/mod_nse_forcing_2decomp.f90)
endif()

set(NSE_MODULE_DIR "${CMAKE_CURRENT_BINARY_DIR}/mod")

set(NSE_MPI_SOURCE src/mpi/module_mpi_TW_NEC.f90)
if(NOT NSE_USE_MPI)
  set(NSE_MPI_SOURCE src/mpi/module_mpi_stub.f90)
endif()

set(NSE_OPENMP_RUNTIME_SOURCE
  src/parallel/mod_openmp_runtime_enabled.f90)
if(NOT NSE_ENABLE_OPENMP)
  set(NSE_OPENMP_RUNTIME_SOURCE
    src/parallel/mod_openmp_runtime_stub.f90)
endif()

set(NSE_GPU_FORTRAN_SOURCE "")
set(NSE_DEFAULT_MAIN_SOURCE src/main/main_nse.f90)
set(NSE_DEFAULT_OUTPUT_NAME solver)
if(NSE_GPU_BACKEND STREQUAL "cuda")
  set(NSE_GPU_FORTRAN_SOURCE src/gpu/mod_nse_gpu_cuda.f90)
  set(NSE_DEFAULT_MAIN_SOURCE src/main/main_nse_cuda.f90)
  set(NSE_DEFAULT_OUTPUT_NAME nse_cuda)

  # nvcc uses MSVC as its host compiler on Windows while the final executable
  # is linked by GNU Fortran. Build the CUDA C ABI as a DLL between them.
  set(NSE_FORTRAN_ARCHIVER "${CMAKE_AR}")
  if(WIN32)
    find_program(NSE_MSVC_LIBRARIAN NAMES lib.exe REQUIRED)
    get_filename_component(NSE_MSVC_BIN_DIR
      "${NSE_MSVC_LIBRARIAN}" DIRECTORY)
    find_program(NSE_MSVC_LINKER NAMES link.exe
      HINTS "${NSE_MSVC_BIN_DIR}" NO_DEFAULT_PATH REQUIRED)
    find_program(NSE_MSVC_RESOURCE_COMPILER NAMES rc.exe REQUIRED)
    set(CMAKE_AR "${NSE_MSVC_LIBRARIAN}")
  endif()
  set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)
  enable_language(CUDA)
  set(CMAKE_CUDA_RUNTIME_LIBRARY None)
  if(WIN32)
    set(CMAKE_LINKER "${NSE_MSVC_LINKER}")
    set(CMAKE_CUDA_HOST_LINK_LAUNCHER "${NSE_MSVC_LINKER}")
    set(CMAKE_RC_COMPILER "${NSE_MSVC_RESOURCE_COMPILER}")
    set(CMAKE_SHARED_LIBRARY_CREATE_CUDA_FLAGS "")
    set(CMAKE_CUDA_CREATE_STATIC_LIBRARY
      "\"${NSE_MSVC_LIBRARIAN}\" /nologo <LINK_FLAGS> /out:<TARGET> <OBJECTS>")
    set(CMAKE_AR "${NSE_FORTRAN_ARCHIVER}")
  endif()
  find_package(CUDAToolkit REQUIRED)
endif()

set(
  NSE_DEFAULT_CORE_SOURCES
  src/common/mod_precision.f90
  src/common/mod_constants.f90
  src/common/mod_common_config.f90
  src/common/mod_model_config.f90
  src/forcing/mod_nse_forcing_common.f90
  src/grid/mod_grid_fvm.f90
  src/field/mod_field_nse.f90
  src/io/mod_input_reader.f90
  ${NSE_MPI_SOURCE}
  ${NSE_OPENMP_RUNTIME_SOURCE}
  src/init/mod_init_taylor_green.f90
  ${NSE_INIT_FFT_SOURCE}
  src/init/mod_nse_initial_conditions.f90
  src/io/mod_slf_output.f90
)
if(NSE_GPU_BACKEND STREQUAL "cuda")
  list(APPEND NSE_DEFAULT_CORE_SOURCES ${NSE_GPU_FORTRAN_SOURCE})
else()
  list(APPEND NSE_DEFAULT_CORE_SOURCES
    ${NSE_BOUNDARY_SOURCE}
    ${NSE_CONVECTIVE_SOURCE}
    ${NSE_VISCOUS_SOURCE}
    ${NSE_FORCING_SOURCE}
    src/numerics/mod_nse_spatial_operator.f90
    src/time/mod_nse_time_integration.f90
  )
endif()
set(
  NSE_CORE_SOURCES "${NSE_DEFAULT_CORE_SOURCES}"
  CACHE STRING
  "Semicolon-separated NSE library source list"
)
set(
  NSE_MAIN_SOURCE "${NSE_DEFAULT_MAIN_SOURCE}"
  CACHE FILEPATH
  "Fortran source containing the NSE main program"
)
set(
  NSE_OUTPUT_NAME "${NSE_DEFAULT_OUTPUT_NAME}"
  CACHE STRING
  "NSE executable file name without suffix"
)

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
if(NSE_INIT_FFT_BACKEND STREQUAL "2decomp_fftw" OR
   NSE_FORCING_FFT_BACKEND STREQUAL "2decomp_fftw")
  target_link_libraries(nse_core PUBLIC decomp2d)
endif()

if(NSE_GPU_BACKEND STREQUAL "cuda")
  add_library(nse_cuda_bridge SHARED src/gpu/nse_cuda_bridge.cu)
  set_target_properties(
    nse_cuda_bridge
    PROPERTIES
      CUDA_STANDARD 17
      CUDA_STANDARD_REQUIRED ON
      CUDA_RUNTIME_LIBRARY None
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
  )
  if(NSE_CUDA_ARCHITECTURES)
    set_property(
      TARGET nse_cuda_bridge
      PROPERTY CUDA_ARCHITECTURES "${NSE_CUDA_ARCHITECTURES}"
    )
  endif()
  if(WIN32)
    target_compile_options(
      nse_cuda_bridge
      PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/utf-8>
        $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/Zc:preprocessor>
    )
  endif()
  target_link_libraries(nse_cuda_bridge PRIVATE CUDA::cudart)
  if(NSE_FORCING_FFT_BACKEND STREQUAL "cufft")
    target_compile_definitions(nse_cuda_bridge PRIVATE NSE_FORCING_CUFFT=1)
    target_link_libraries(nse_cuda_bridge PRIVATE CUDA::cufft)
  endif()
  target_link_libraries(nse_core PUBLIC nse_cuda_bridge)
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

include(CTest)
if(BUILD_TESTING AND NSE_ENABLE_CUFFTMP_CONTRACT_TESTS)
  add_test(
    NAME nse_cufftmp_forcing_contract
    COMMAND "${CMAKE_COMMAND}"
      "-DCONTRACT_FILE=${CMAKE_CURRENT_SOURCE_DIR}/src/forcing/cufftmp_forcing_contract.yaml"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/validate_cufftmp_forcing_contract.cmake"
  )
  set_tests_properties(
    nse_cufftmp_forcing_contract
    PROPERTIES PASS_REGULAR_EXPRESSION "cuFFTMp forcing contract test passed"
  )
endif()
if(BUILD_TESTING)
  add_executable(
    nse_test_forcing_coefficients
    tests/test_forcing_coefficients.f90
  )
  set_target_properties(
    nse_test_forcing_coefficients
    PROPERTIES
      Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}/forcing_coefficient_test"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
  )
  target_link_libraries(nse_test_forcing_coefficients PRIVATE NSE::core)
  nse_apply_fortran_options(nse_test_forcing_coefficients)
  add_test(
    NAME nse_forcing_coefficients
    COMMAND $<TARGET_FILE:nse_test_forcing_coefficients>
  )
  set_tests_properties(
    nse_forcing_coefficients
    PROPERTIES
      PASS_REGULAR_EXPRESSION
        "Petersen-Livescu forcing coefficient test passed"
  )
endif()
if(BUILD_TESTING AND NSE_GPU_BACKEND STREQUAL "none")
  add_executable(nse_test_keep6_accuracy tests/test_keep6_accuracy.f90)
  set_target_properties(
    nse_test_keep6_accuracy
    PROPERTIES
      Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}/keep6_test"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
  )
  target_link_libraries(nse_test_keep6_accuracy PRIVATE NSE::core)
  nse_apply_fortran_options(nse_test_keep6_accuracy)
  add_test(NAME nse_keep6_accuracy
    COMMAND $<TARGET_FILE:nse_test_keep6_accuracy>)
  set_tests_properties(
    nse_keep6_accuracy
    PROPERTIES
      PASS_REGULAR_EXPRESSION "KEEP selectable-order accuracy test passed"
  )
endif()
if(BUILD_TESTING AND NSE_GPU_BACKEND STREQUAL "none" AND
   NSE_VISCOUS_SCHEME STREQUAL "central6")
  add_executable(nse_test_viscous_central6 tests/test_viscous_central6.f90)
  set_target_properties(
    nse_test_viscous_central6
    PROPERTIES
      Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}/viscous_test"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
  )
  target_link_libraries(nse_test_viscous_central6 PRIVATE NSE::core)
  nse_apply_fortran_options(nse_test_viscous_central6)
  add_test(NAME nse_viscous_central6
    COMMAND $<TARGET_FILE:nse_test_viscous_central6>)
  set_tests_properties(
    nse_viscous_central6
    PROPERTIES PASS_REGULAR_EXPRESSION "central6 viscous test passed"
  )
endif()
if(BUILD_TESTING AND NSE_USE_MPI AND NSE_GPU_BACKEND STREQUAL "none")
  add_executable(nse_test_boundary_periodic tests/test_boundary_periodic.f90)
  set_target_properties(
    nse_test_boundary_periodic
    PROPERTIES
      Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}/boundary_test"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
  )
  target_link_libraries(nse_test_boundary_periodic PRIVATE NSE::core)
  nse_apply_fortran_options(nse_test_boundary_periodic)
  add_test(
    NAME nse_boundary_periodic
    COMMAND mpiexec -n 4 $<TARGET_FILE:nse_test_boundary_periodic>
  )
  set_tests_properties(
    nse_boundary_periodic
    PROPERTIES
      PASS_REGULAR_EXPRESSION
        "periodic boundary face/edge/corner test passed"
  )
  add_test(
    NAME nse_mpi_normal_completion
    COMMAND mpiexec -n 4 $<TARGET_FILE:nse_solver>
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_tgv_viscous_small.dat"
  )
  set_tests_properties(
    nse_mpi_normal_completion
    PROPERTIES
      PASS_REGULAR_EXPRESSION
        "NSE calculation completed successfully:"
      FAIL_REGULAR_EXPRESSION
        "floating-point exceptions are signalling"
  )
endif()
if(BUILD_TESTING AND NSE_USE_MPI AND
   NSE_FORCING_FFT_BACKEND STREQUAL "2decomp_fftw")
  add_test(
    NAME nse_forcing_2decomp_smoke
    COMMAND mpiexec -n 4 $<TARGET_FILE:nse_solver>
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_tgv_forcing_small.dat"
  )
  set_tests_properties(
    nse_forcing_2decomp_smoke
    PROPERTIES
      PASS_REGULAR_EXPRESSION
        "NSE calculation completed successfully:"
      FAIL_REGULAR_EXPRESSION
        "floating-point exceptions are signalling"
  )
  add_test(
    NAME nse_forcing_hit_2decomp_smoke
    COMMAND mpiexec -n 4 $<TARGET_FILE:nse_solver>
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_hit_forcing_small.dat"
  )
  set_tests_properties(
    nse_forcing_hit_2decomp_smoke
    PROPERTIES
      PASS_REGULAR_EXPRESSION
        "NSE calculation completed successfully:"
      FAIL_REGULAR_EXPRESSION
        "floating-point exceptions are signalling"
  )
endif()
if(BUILD_TESTING AND NSE_GPU_BACKEND STREQUAL "cuda")
  add_executable(nse_test_cuda_keep tests/test_cuda_keep.f90)
  set_target_properties(
    nse_test_cuda_keep
    PROPERTIES
      Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
  )
  target_link_libraries(nse_test_cuda_keep PRIVATE NSE::core)
  nse_apply_fortran_options(nse_test_cuda_keep)
  add_test(NAME nse_cuda_keep_uniform
    COMMAND $<TARGET_FILE:nse_test_cuda_keep>)
  set_tests_properties(
    nse_cuda_keep_uniform
    PROPERTIES PASS_REGULAR_EXPRESSION "CUDA KEEP test passed"
  )
  add_executable(
    nse_test_cuda_keep_compare
    tests/test_cuda_keep_compare.f90
    src/numerics/convective/mod_convective_keep.f90
  )
  set_target_properties(
    nse_test_cuda_keep_compare
    PROPERTIES
      Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}/cuda_compare"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
  )
  target_link_libraries(nse_test_cuda_keep_compare PRIVATE NSE::core)
  nse_apply_fortran_options(nse_test_cuda_keep_compare)
  add_test(NAME nse_cuda_keep_compare
    COMMAND $<TARGET_FILE:nse_test_cuda_keep_compare> 6)
  set_tests_properties(
    nse_cuda_keep_compare
    PROPERTIES PASS_REGULAR_EXPRESSION "CUDA KEEP comparison passed"
  )
  add_test(NAME nse_cuda_keep2_compare
    COMMAND $<TARGET_FILE:nse_test_cuda_keep_compare> 2)
  set_tests_properties(
    nse_cuda_keep2_compare
    PROPERTIES PASS_REGULAR_EXPRESSION "CUDA KEEP comparison passed"
  )
  add_executable(
    nse_test_cuda_viscous_compare
    tests/test_cuda_viscous_compare.f90
    src/numerics/convective/mod_convective_keep.f90
    src/numerics/viscous/mod_viscous_central6.f90
  )
  set_target_properties(
    nse_test_cuda_viscous_compare
    PROPERTIES
      Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}/cuda_viscous_compare"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
  )
  target_link_libraries(nse_test_cuda_viscous_compare PRIVATE NSE::core)
  nse_apply_fortran_options(nse_test_cuda_viscous_compare)
  add_test(NAME nse_cuda_viscous_compare
    COMMAND $<TARGET_FILE:nse_test_cuda_viscous_compare>)
  set_tests_properties(
    nse_cuda_viscous_compare
    PROPERTIES PASS_REGULAR_EXPRESSION "CUDA central6 comparison passed"
  )
  add_test(
    NAME nse_cuda_tgv_smoke
    COMMAND $<TARGET_FILE:nse_solver>
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_tgv_cuda_small.dat"
  )
  set_tests_properties(
    nse_cuda_tgv_smoke
    PROPERTIES PASS_REGULAR_EXPRESSION "# step time dt"
  )
  if(NSE_FORCING_FFT_BACKEND STREQUAL "cufft")
    add_test(
      NAME nse_cuda_forcing_smoke
      COMMAND $<TARGET_FILE:nse_solver>
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_tgv_cuda_forcing_small.dat"
    )
    set_tests_properties(
      nse_cuda_forcing_smoke
      PROPERTIES
        PASS_REGULAR_EXPRESSION
          "NSE calculation completed successfully:"
        FAIL_REGULAR_EXPRESSION
          "CUDA backend failure"
    )
  endif()
endif()

message(STATUS "NSE MPI provider: ${NSE_MPI_PROVIDER_RESOLVED}")
message(STATUS "NSE MPI enabled: ${NSE_USE_MPI}")
message(STATUS "NSE GPU backend: ${NSE_GPU_BACKEND}")
message(STATUS "NSE convective backend: ${NSE_CONVECTIVE_BACKEND}")
message(STATUS "NSE viscous scheme: ${NSE_VISCOUS_SCHEME}")
message(STATUS "NSE boundary scheme: ${NSE_BOUNDARY_SCHEME}")
message(STATUS "NSE init FFT backend: ${NSE_INIT_FFT_BACKEND}")
message(STATUS "NSE forcing FFT backend: ${NSE_FORCING_FFT_BACKEND}")
if(NSE_INIT_FFT_BACKEND STREQUAL "2decomp_fftw")
  message(STATUS "NSE 2DECOMP root: ${NSE_2DECOMP_ROOT}")
endif()
message(STATUS "NSE Fortran modules: ${NSE_MODULE_DIR}")
message(STATUS "NSE executable: ${CMAKE_CURRENT_BINARY_DIR}/bin/${NSE_OUTPUT_NAME}${CMAKE_EXECUTABLE_SUFFIX}")
