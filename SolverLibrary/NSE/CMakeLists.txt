cmake_minimum_required(VERSION 3.24)

# Windowsの長いソースパスでは、CMakeにオブジェクト名の短縮を早めに行わせる。
set(CMAKE_OBJECT_PATH_MAX 240)

project(
  NSE
  VERSION 1.0.0
  DESCRIPTION "MPI/OpenMP, single-GPU CUDA, and MPI+CUDA Navier-Stokes solver"
  LANGUAGES Fortran
)

include(GNUInstallDirs)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(NSEMPI)

option(NSE_ENABLE_WARNINGS "Enable additional Fortran compiler warnings" ON)
option(NSE_ENABLE_OPENMP "Enable OpenMP support" ON)
option(NSE_USE_MPI "Build the MPI domain-decomposition backend" ON)
option(NSE_ENABLE_CUDA_AWARE_MPI "Enable single-component CUDA device-buffer MPI halo exchange" OFF)
set(
  NSE_TEST_SUITE "legacy"
  CACHE STRING
  "CTest suite: legacy, multicomponent_foundation, multicomponent_passive_scalar, multicomponent_euler, multicomponent_thermally_perfect, multicomponent_viscous, multicomponent_reactor, multicomponent_reactive, multicomponent_reactive_boundaries, or none"
)
set_property(
  CACHE NSE_TEST_SUITE
  PROPERTY STRINGS legacy multicomponent_foundation
    multicomponent_passive_scalar multicomponent_euler
    multicomponent_thermally_perfect multicomponent_viscous
    multicomponent_reactor multicomponent_reactive
    multicomponent_reactive_boundaries none
)
option(
  NSE_ENABLE_MULTI_GPU_TESTS
  "Run MPI+CUDA tests that require one visible GPU per MPI rank"
  OFF
)
set(
  NSE_MULTI_GPU_TEST_PROCESSES 4
  CACHE STRING
  "MPI process/GPU count for the optional MPI+CUDA smoke test"
)
option(
  NSE_ENABLE_CUFFTMP_CONTRACT_TESTS
  "Validate the MPI+CUDA cuFFTMp forcing interface contract"
  ON
)
set(
  NSE_GPU_BACKEND "none"
  CACHE STRING
  "GPU backend: none or cuda"
)
set_property(CACHE NSE_GPU_BACKEND PROPERTY STRINGS none cuda)
option(NSE_MC_CUDA "Build the separate resident multicomponent CUDA bridge" OFF)
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
  NSE_CONVECTIVE_BACKEND "runtime"
  CACHE STRING
  "Convective-flux backend compiled into the NSE solver"
)
set_property(CACHE NSE_CONVECTIVE_BACKEND PROPERTY STRINGS runtime keep_family)
set(
  NSE_VISCOUS_SCHEME "central6"
  CACHE STRING
  "Viscous-term backend compiled into the NSE solver"
)
set_property(CACHE NSE_VISCOUS_SCHEME PROPERTY STRINGS none central6)
set(
  NSE_BOUNDARY_SCHEME "runtime"
  CACHE STRING
  "Boundary-condition backend compiled into the NSE solver"
)
set_property(CACHE NSE_BOUNDARY_SCHEME PROPERTY STRINGS runtime periodic)
set(
  NSE_INIT_FFT_BACKEND "none"
  CACHE STRING
  "Initial-condition FFT backend: none, 2decomp_fftw, cufft, or cufftmp"
)
set_property(
  CACHE NSE_INIT_FFT_BACKEND
  PROPERTY STRINGS none 2decomp_fftw cufft cufftmp
)
set(
  NSE_FORCING_FFT_BACKEND "none"
  CACHE STRING
  "Forcing Helmholtz FFT backend: none, 2decomp_fftw, cufft, or cufftmp"
)
set_property(
  CACHE NSE_FORCING_FFT_BACKEND
  PROPERTY STRINGS none 2decomp_fftw cufft cufftmp
)
set(
  NSE_2DECOMP_ROOT "$ENV{DECOMP2D_ROOT}"
  CACHE PATH
  "2DECOMP&FFT installation prefix containing lib/decomp2d"
)
set(
  NSE_CUFFTMP_ROOT "$ENV{CUFFTMP_ROOT}"
  CACHE PATH
  "cuFFTMp prefix containing include/cufftmp and lib64/libcufftMp"
)
set(
  NSE_NVSHMEM_ROOT "$ENV{NVSHMEM_ROOT}"
  CACHE PATH
  "NVSHMEM prefix matching the selected cuFFTMp installation"
)

string(TOLOWER "${NSE_CONVECTIVE_BACKEND}" NSE_CONVECTIVE_BACKEND)
string(TOLOWER "${NSE_VISCOUS_SCHEME}" NSE_VISCOUS_SCHEME)
string(TOLOWER "${NSE_BOUNDARY_SCHEME}" NSE_BOUNDARY_SCHEME)
string(TOLOWER "${NSE_GPU_BACKEND}" NSE_GPU_BACKEND)
string(TOLOWER "${NSE_INIT_FFT_BACKEND}" NSE_INIT_FFT_BACKEND)
string(TOLOWER "${NSE_FORCING_FFT_BACKEND}" NSE_FORCING_FFT_BACKEND)
string(TOLOWER "${NSE_TEST_SUITE}" NSE_TEST_SUITE)

if(NOT NSE_TEST_SUITE STREQUAL "legacy" AND
   NOT NSE_TEST_SUITE STREQUAL "multicomponent_foundation" AND
   NOT NSE_TEST_SUITE STREQUAL "multicomponent_passive_scalar" AND
   NOT NSE_TEST_SUITE STREQUAL "multicomponent_euler" AND
   NOT NSE_TEST_SUITE STREQUAL "multicomponent_thermally_perfect" AND
   NOT NSE_TEST_SUITE STREQUAL "multicomponent_viscous" AND
   NOT NSE_TEST_SUITE STREQUAL "multicomponent_reactor" AND
   NOT NSE_TEST_SUITE STREQUAL "multicomponent_reactive" AND
   NOT NSE_TEST_SUITE STREQUAL "multicomponent_reactive_boundaries" AND
   NOT NSE_TEST_SUITE STREQUAL "none")
  message(FATAL_ERROR
    "Unsupported NSE_TEST_SUITE='${NSE_TEST_SUITE}'. "
    "Available suites: legacy, multicomponent_foundation, "
    "multicomponent_passive_scalar, multicomponent_euler, "
    "multicomponent_thermally_perfect, multicomponent_viscous, "
    "multicomponent_reactor, multicomponent_reactive, "
    "multicomponent_reactive_boundaries, none"
  )
endif()

if(NOT NSE_GPU_BACKEND STREQUAL "none" AND
   NOT NSE_GPU_BACKEND STREQUAL "cuda")
  message(FATAL_ERROR
    "Unsupported NSE_GPU_BACKEND='${NSE_GPU_BACKEND}'. "
    "Available backends: none, cuda"
  )
endif()
if(NOT NSE_INIT_FFT_BACKEND STREQUAL "none" AND
   NOT NSE_INIT_FFT_BACKEND STREQUAL "2decomp_fftw" AND
   NOT NSE_INIT_FFT_BACKEND STREQUAL "cufft" AND
   NOT NSE_INIT_FFT_BACKEND STREQUAL "cufftmp")
  message(FATAL_ERROR
    "Unsupported NSE_INIT_FFT_BACKEND='${NSE_INIT_FFT_BACKEND}'. "
    "Available backends: none, 2decomp_fftw, cufft, cufftmp"
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
if(NSE_INIT_FFT_BACKEND STREQUAL "cufft" AND
   NOT NSE_GPU_BACKEND STREQUAL "cuda")
  message(FATAL_ERROR
    "The cufft initialization backend requires NSE_GPU_BACKEND=cuda."
  )
endif()
if(NSE_INIT_FFT_BACKEND STREQUAL "cufftmp" AND
   (NOT NSE_USE_MPI OR NOT NSE_GPU_BACKEND STREQUAL "cuda"))
  message(FATAL_ERROR
    "The cufftmp initialization backend requires NSE_USE_MPI=ON and "
    "NSE_GPU_BACKEND=cuda."
  )
endif()
if(NSE_INIT_FFT_BACKEND STREQUAL "cufft" AND NSE_USE_MPI)
  message(FATAL_ERROR
    "The cufft initialization backend is single-GPU only. "
    "Use NSE_INIT_FFT_BACKEND=none for MPI+CUDA."
  )
endif()
if(NOT NSE_FORCING_FFT_BACKEND STREQUAL "none" AND
   NOT NSE_FORCING_FFT_BACKEND STREQUAL "2decomp_fftw" AND
   NOT NSE_FORCING_FFT_BACKEND STREQUAL "cufft" AND
   NOT NSE_FORCING_FFT_BACKEND STREQUAL "cufftmp")
  message(FATAL_ERROR
    "Unsupported NSE_FORCING_FFT_BACKEND='${NSE_FORCING_FFT_BACKEND}'. "
    "Available runtime backends: none, 2decomp_fftw, cufft, cufftmp."
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
if(NSE_FORCING_FFT_BACKEND STREQUAL "cufftmp" AND
   (NOT NSE_USE_MPI OR NOT NSE_GPU_BACKEND STREQUAL "cuda"))
  message(FATAL_ERROR
    "The cufftmp forcing backend requires NSE_USE_MPI=ON and "
    "NSE_GPU_BACKEND=cuda."
  )
endif()
if((NSE_INIT_FFT_BACKEND STREQUAL "cufftmp" OR
    NSE_FORCING_FFT_BACKEND STREQUAL "cufftmp") AND WIN32)
  message(FATAL_ERROR
    "cuFFTMp is supported by this project on Linux HPC systems only."
  )
endif()
if(NSE_FORCING_FFT_BACKEND STREQUAL "cufft" AND NSE_USE_MPI)
  message(FATAL_ERROR
    "The cufft forcing backend is single-GPU only. "
    "Use NSE_FORCING_FFT_BACKEND=none for MPI+CUDA."
  )
endif()

if(NSE_CONVECTIVE_BACKEND STREQUAL "runtime" OR
   NSE_CONVECTIVE_BACKEND STREQUAL "keep_family")
  set(
    NSE_CONVECTIVE_SOURCES
    src/numerics/reconstruction/mod_reconstruction_weno5z.f90
    src/numerics/riemann/mod_riemann_roe.f90
    src/numerics/convective/mod_convective_keep.f90
    src/numerics/convective/mod_convective_weno5z_roe.f90
    src/numerics/convective/mod_convective_leaf_registry.f90
    src/numerics/convective/mod_convective_hybrid.f90
    src/numerics/convective/mod_convective_dispatch.f90
  )
else()
  message(FATAL_ERROR
    "Unsupported NSE_CONVECTIVE_BACKEND='${NSE_CONVECTIVE_BACKEND}'. "
    "Available backend: runtime"
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

if(NSE_BOUNDARY_SCHEME STREQUAL "runtime")
  set(NSE_BOUNDARY_SOURCE src/boundary/mod_boundary_runtime.f90)
elseif(NSE_BOUNDARY_SCHEME STREQUAL "periodic")
  set(NSE_BOUNDARY_SOURCE src/boundary/mod_boundary_periodic.f90)
else()
  message(FATAL_ERROR
    "Unsupported NSE_BOUNDARY_SCHEME='${NSE_BOUNDARY_SCHEME}'. "
    "Available backends: runtime, periodic"
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
set(NSE_INIT_FFT_SUPPORT_SOURCE "")
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
    set_source_files_properties(
      src/init/mod_init_hit_spectral_2decomp.f90
      PROPERTIES Fortran_PREPROCESS ON
    )
  endif()
endif()
if(NSE_INIT_FFT_BACKEND STREQUAL "cufft")
  set(NSE_INIT_FFT_SOURCE src/init/mod_init_hit_spectral_cufft.f90)
endif()
if(NSE_INIT_FFT_BACKEND STREQUAL "cufftmp")
  set(NSE_INIT_FFT_SUPPORT_SOURCE src/gpu/mod_nse_cufftmp_fft.f90)
  set(NSE_INIT_FFT_SOURCE src/init/mod_init_hit_spectral_2decomp.f90)
  set_source_files_properties(
    src/init/mod_init_hit_spectral_2decomp.f90
    PROPERTIES Fortran_PREPROCESS ON
  )
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
if(NSE_ENABLE_CUDA_AWARE_MPI AND NOT NSE_GPU_BACKEND STREQUAL "cuda")
  message(FATAL_ERROR "NSE_ENABLE_CUDA_AWARE_MPI requires NSE_GPU_BACKEND=cuda")
endif()
if(NSE_GPU_BACKEND STREQUAL "cuda")
  set(NSE_GPU_FORTRAN_SOURCE src/gpu/mod_nse_gpu_cuda.f90)
  if(NSE_USE_MPI)
    list(APPEND NSE_GPU_FORTRAN_SOURCE src/gpu/mod_nse_gpu_mpi.f90)
    set(NSE_DEFAULT_MAIN_SOURCE src/main/main_nse_mpi_cuda.f90)
    set(NSE_DEFAULT_OUTPUT_NAME nse_mpi_cuda)
  else()
    set(NSE_DEFAULT_MAIN_SOURCE src/main/main_nse_cuda.f90)
    set(NSE_DEFAULT_OUTPUT_NAME nse_cuda)
  endif()

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

set(NSE_CUFFTMP_ENABLED OFF)
if(NSE_INIT_FFT_BACKEND STREQUAL "cufftmp" OR
   NSE_FORCING_FFT_BACKEND STREQUAL "cufftmp")
  set(NSE_CUFFTMP_ENABLED ON)
  find_path(
    NSE_CUFFTMP_INCLUDE_DIR
    NAMES cufftMp.h
    HINTS
      "${NSE_CUFFTMP_ROOT}"
      "$ENV{CUFFT_INC}"
    PATH_SUFFIXES include include/cufftmp math_libs/include/cufftmp
    REQUIRED
  )
  find_library(
    NSE_CUFFTMP_LIBRARY
    NAMES cufftMp
    HINTS
      "${NSE_CUFFTMP_ROOT}"
      "$ENV{CUFFT_LIB}"
    PATH_SUFFIXES lib lib64 math_libs/lib64
    REQUIRED
  )
  find_path(
    NSE_NVSHMEM_INCLUDE_DIR
    NAMES nvshmem.h
    HINTS
      "${NSE_NVSHMEM_ROOT}"
      "$ENV{NVSHMEM_INC}"
    PATH_SUFFIXES include comm_libs/nvshmem/include
    REQUIRED
  )
  find_library(
    NSE_NVSHMEM_HOST_LIBRARY
    NAMES nvshmem_host
    HINTS
      "${NSE_NVSHMEM_ROOT}"
      "$ENV{NVSHMEM_LIB}"
    PATH_SUFFIXES lib lib64 comm_libs/nvshmem/lib
    REQUIRED
  )
  find_path(
    NSE_CUFFTMP_MPI_INCLUDE_DIR
    NAMES mpi.h
    HINTS "$ENV{MPI_HOME}" ${MPI_Fortran_INCLUDE_DIRS}
    PATH_SUFFIXES include
    REQUIRED
  )
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
  src/init/mod_init_imported_turbulence.f90
  src/init/mod_init_shock_turbulence.f90
  src/init/mod_init_shock_tube_turbulence.f90
  src/init/mod_hit_isotropy_math.f90
  ${NSE_INIT_FFT_SUPPORT_SOURCE}
  ${NSE_INIT_FFT_SOURCE}
  src/init/mod_nse_initial_conditions.f90
  src/io/mod_slf_output.f90
  ${NSE_BOUNDARY_SOURCE}
)
if(NSE_GPU_BACKEND STREQUAL "cuda")
  list(APPEND NSE_DEFAULT_CORE_SOURCES ${NSE_GPU_FORTRAN_SOURCE})
else()
  list(APPEND NSE_DEFAULT_CORE_SOURCES
    ${NSE_CONVECTIVE_SOURCES}
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

# Keep optional FH sources available even in previously configured build trees
# whose NSE_CORE_SOURCES cache predates this extension.
if("${NSE_CORE_SOURCES}" MATCHES "mod_nse_spatial_operator[.]f90|mod_nse_gpu_cuda[.]f90")
  if(NOT "${NSE_CORE_SOURCES}" MATCHES "mod_fh_random[.]f90")
    list(APPEND NSE_CORE_SOURCES src/extensions/fluctuating/mod_fh_random.f90)
  endif()
  if(NOT "${NSE_CORE_SOURCES}" MATCHES "mod_nse_fluctuating[.]f90")
    list(APPEND NSE_CORE_SOURCES src/extensions/fluctuating/mod_nse_fluctuating.f90)
  endif()
endif()
list(REMOVE_DUPLICATES NSE_CORE_SOURCES)
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
if(NSE_INIT_FFT_BACKEND STREQUAL "cufftmp")
  target_compile_definitions(nse_core PRIVATE NSE_INIT_CUFFTMP=1)
endif()

if(NSE_GPU_BACKEND STREQUAL "cuda")
  if(NSE_MC_CUDA)
    add_library(nse_cuda_bridge SHARED src/extensions/multicomponent/gpu/mc_cuda.cu)
    target_compile_options(nse_cuda_bridge PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--fmad=false>)
  else()
    add_library(nse_cuda_bridge SHARED src/gpu/nse_cuda_bridge.cu)
  endif()
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
  if(NSE_ENABLE_CUDA_AWARE_MPI)
    if(NOT NSE_USE_MPI OR NSE_MPI_PROVIDER_RESOLVED STREQUAL "MSMPI" OR
       NSE_TEST_SUITE MATCHES "^multicomponent")
      message(FATAL_ERROR "CUDA-aware halo exchange requires single-component NSE and SYSTEM MPI (not MS-MPI)")
    endif()
    enable_language(C)
    find_package(MPI REQUIRED COMPONENTS C)
    target_compile_definitions(nse_cuda_bridge PRIVATE NSE_CUDA_AWARE_MPI=1 OMPI_SKIP_MPICXX=1 MPICH_SKIP_MPICXX=1)
    target_link_libraries(nse_cuda_bridge PRIVATE MPI::MPI_C)
  endif()
  if(NSE_FORCING_FFT_BACKEND STREQUAL "cufft")
    target_compile_definitions(nse_cuda_bridge PRIVATE NSE_FORCING_CUFFT=1)
  endif()
  if(NSE_FORCING_FFT_BACKEND STREQUAL "cufftmp")
    target_compile_definitions(nse_cuda_bridge PRIVATE NSE_FORCING_CUFFTMP=1)
  endif()
  if(NSE_INIT_FFT_BACKEND STREQUAL "cufft")
    target_compile_definitions(nse_cuda_bridge PRIVATE NSE_INIT_CUFFT=1)
  endif()
  if(NSE_INIT_FFT_BACKEND STREQUAL "cufftmp")
    target_compile_definitions(nse_cuda_bridge PRIVATE NSE_INIT_CUFFTMP=1)
  endif()
  if(NSE_FORCING_FFT_BACKEND STREQUAL "cufft" OR
     NSE_INIT_FFT_BACKEND STREQUAL "cufft" OR
     NSE_FORCING_FFT_BACKEND STREQUAL "cufftmp" OR
     NSE_INIT_FFT_BACKEND STREQUAL "cufftmp")
    target_link_libraries(nse_cuda_bridge PRIVATE CUDA::cufft)
  endif()
  if(NSE_CUFFTMP_ENABLED)
    target_include_directories(
      nse_cuda_bridge
      PRIVATE
        "${NSE_CUFFTMP_INCLUDE_DIR}"
        "${NSE_NVSHMEM_INCLUDE_DIR}"
        "${NSE_CUFFTMP_MPI_INCLUDE_DIR}"
    )
    target_link_libraries(
      nse_cuda_bridge
      PRIVATE
        nse_mpi
        "${NSE_CUFFTMP_LIBRARY}"
        "${NSE_NVSHMEM_HOST_LIBRARY}"
    )
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
if(BUILD_TESTING AND
   (NSE_TEST_SUITE STREQUAL "multicomponent_foundation" OR
    NSE_TEST_SUITE STREQUAL "multicomponent_passive_scalar" OR
    NSE_TEST_SUITE STREQUAL "multicomponent_euler" OR
    NSE_TEST_SUITE STREQUAL "multicomponent_thermally_perfect" OR
    NSE_TEST_SUITE STREQUAL "multicomponent_viscous" OR
    NSE_TEST_SUITE STREQUAL "multicomponent_reactor" OR
    NSE_TEST_SUITE STREQUAL "multicomponent_reactive" OR
    NSE_TEST_SUITE STREQUAL "multicomponent_reactive_boundaries"))
  add_executable(
    nse_test_multicomponent_foundation
    tests/test_multicomponent_foundation.f90
  )
  set_target_properties(
    nse_test_multicomponent_foundation
    PROPERTIES
      Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}/multicomponent_foundation_test"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
  )
  target_link_libraries(nse_test_multicomponent_foundation PRIVATE NSE::core)
  nse_apply_fortran_options(nse_test_multicomponent_foundation)
  add_test(
    NAME nse_multicomponent_foundation
    COMMAND $<TARGET_FILE:nse_test_multicomponent_foundation>
  )
  set_tests_properties(
    nse_multicomponent_foundation
    PROPERTIES
      PASS_REGULAR_EXPRESSION
        "Multicomponent foundation contract tests passed"
  )
  if(NSE_TEST_SUITE STREQUAL "multicomponent_foundation")
    add_test(
      NAME nse_multicomponent_foundation_input
      COMMAND $<TARGET_FILE:nse_solver>
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_multicomponent_foundation.dat"
    )
    set_tests_properties(
      nse_multicomponent_foundation_input
      PROPERTIES
        PASS_REGULAR_EXPRESSION
          "Multicomponent NSE stage-0 foundation initialized successfully"
    )
  endif()
endif()

if(BUILD_TESTING AND NSE_TEST_SUITE STREQUAL "multicomponent_reactor")
  add_executable(
    nse_test_multicomponent_reactor
    tests/test_multicomponent_reactor.f90
  )
  set_target_properties(
    nse_test_multicomponent_reactor
    PROPERTIES
      Fortran_MODULE_DIRECTORY
        "${NSE_MODULE_DIR}/multicomponent_reactor_test"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
  )
  target_link_libraries(nse_test_multicomponent_reactor PRIVATE NSE::core)
  nse_apply_fortran_options(nse_test_multicomponent_reactor)
  add_test(
    NAME nse_multicomponent_reactor
    COMMAND $<TARGET_FILE:nse_test_multicomponent_reactor>
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_multicomponent_reactor.dat"
  )
  set_tests_properties(
    nse_multicomponent_reactor
    PROPERTIES
      PASS_REGULAR_EXPRESSION
        "Homogeneous finite-rate chemistry tests passed"
      FAIL_REGULAR_EXPRESSION
        "floating-point exceptions are signalling;ERROR STOP"
  )
  add_test(
    NAME nse_multicomponent_reactor_input
    COMMAND $<TARGET_FILE:nse_solver>
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_multicomponent_reactor.dat"
  )
  set_tests_properties(
    nse_multicomponent_reactor_input
    PROPERTIES
      PASS_REGULAR_EXPRESSION
        "Homogeneous reactor calculation completed successfully"
      FAIL_REGULAR_EXPRESSION
        "floating-point exceptions are signalling;ERROR STOP"
  )
endif()

if(BUILD_TESTING AND
   (NSE_TEST_SUITE STREQUAL "multicomponent_reactive" OR
    NSE_TEST_SUITE STREQUAL "multicomponent_reactive_boundaries"))
  add_executable(
    nse_test_multicomponent_reactive
    tests/test_multicomponent_reactive.f90
  )
  set_target_properties(
    nse_test_multicomponent_reactive
    PROPERTIES
      Fortran_MODULE_DIRECTORY
        "${NSE_MODULE_DIR}/multicomponent_reactive_test"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
  )
  target_link_libraries(nse_test_multicomponent_reactive PRIVATE NSE::core)
  nse_apply_fortran_options(nse_test_multicomponent_reactive)
  add_test(
    NAME nse_multicomponent_reactive
    COMMAND $<TARGET_FILE:nse_test_multicomponent_reactive>
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_multicomponent_reactive.dat"
  )
  set_tests_properties(
    nse_multicomponent_reactive
    PROPERTIES
      PASS_REGULAR_EXPRESSION
        "Reactive multicomponent Strang-splitting tests passed"
      FAIL_REGULAR_EXPRESSION
        "floating-point exceptions are signalling;ERROR STOP"
  )
  add_test(
    NAME nse_multicomponent_reactive_input
    COMMAND $<TARGET_FILE:nse_solver>
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_multicomponent_reactive.dat"
  )
  set_tests_properties(
    nse_multicomponent_reactive_input
    PROPERTIES
      PASS_REGULAR_EXPRESSION
        "Reactive multicomponent Navier-Stokes calculation completed successfully"
      FAIL_REGULAR_EXPRESSION
        "floating-point exceptions are signalling;ERROR STOP"
  )
endif()

if(BUILD_TESTING AND
   NSE_TEST_SUITE STREQUAL "multicomponent_reactive_boundaries")
  add_executable(nse_test_multicomponent_geometry tests/test_multicomponent_geometry.f90)
  set_target_properties(nse_test_multicomponent_geometry PROPERTIES
    Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}/mc_geometry_test"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
  target_link_libraries(nse_test_multicomponent_geometry PRIVATE NSE::core)
  nse_apply_fortran_options(nse_test_multicomponent_geometry)
  add_test(NAME nse_mc_geometry COMMAND $<TARGET_FILE:nse_test_multicomponent_geometry>
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_multicomponent_reactive_boundaries.dat")
  set_tests_properties(nse_mc_geometry PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")
  if(NSE_USE_MPI)
    add_executable(nse_test_multicomponent_mpi_pencil tests/test_multicomponent_mpi_pencil.f90)
    set_target_properties(nse_test_multicomponent_mpi_pencil PROPERTIES
      Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}/mc_mpi_test"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
    target_link_libraries(nse_test_multicomponent_mpi_pencil PRIVATE NSE::core)
    nse_apply_fortran_options(nse_test_multicomponent_mpi_pencil)
    foreach(ranks IN ITEMS 1 2 3 4)
      add_test(NAME nse_mc_pencil_${ranks}
        COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${ranks}
          ${MPIEXEC_PREFLAGS} $<TARGET_FILE:nse_test_multicomponent_mpi_pencil>
          ${MPIEXEC_POSTFLAGS}
          "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_multicomponent_reactive_boundaries.dat")
      set_tests_properties(nse_mc_pencil_${ranks} PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=2;OMP_DYNAMIC=FALSE"
        TIMEOUT 60 PROCESSORS ${ranks}
        PASS_REGULAR_EXPRESSION "MPI pencil transverse/corner and Strang regression passed"
        FAIL_REGULAR_EXPRESSION "ERROR STOP;Mismatch")
    endforeach()
  endif()
  add_executable(nse_test_multicomponent_rhs_workspace
    tests/test_multicomponent_rhs_workspace.f90)
  set_target_properties(nse_test_multicomponent_rhs_workspace PROPERTIES
    Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}/multicomponent_rhs_workspace_test"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
  target_link_libraries(nse_test_multicomponent_rhs_workspace PRIVATE NSE::core)
  nse_apply_fortran_options(nse_test_multicomponent_rhs_workspace)
  add_test(NAME nse_multicomponent_rhs_workspace
    COMMAND $<TARGET_FILE:nse_test_multicomponent_rhs_workspace>
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_multicomponent_reactive_boundaries.dat")
  set_tests_properties(nse_multicomponent_rhs_workspace PROPERTIES
    PASS_REGULAR_EXPRESSION "Multicomponent RHS workspace regression tests passed"
    FAIL_REGULAR_EXPRESSION "floating-point exceptions are signalling;ERROR STOP")
  add_executable(
    nse_test_multicomponent_reactive_boundaries
    tests/test_multicomponent_reactive_boundaries.f90
  )
  set_target_properties(
    nse_test_multicomponent_reactive_boundaries
    PROPERTIES
      Fortran_MODULE_DIRECTORY
        "${NSE_MODULE_DIR}/multicomponent_reactive_boundary_test"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
  )
  target_link_libraries(
    nse_test_multicomponent_reactive_boundaries PRIVATE NSE::core
  )
  nse_apply_fortran_options(nse_test_multicomponent_reactive_boundaries)
  add_test(
    NAME nse_multicomponent_reactive_boundaries
    COMMAND $<TARGET_FILE:nse_test_multicomponent_reactive_boundaries>
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_multicomponent_reactive_boundaries.dat"
  )
  set_tests_properties(
    nse_multicomponent_reactive_boundaries
    PROPERTIES
      PASS_REGULAR_EXPRESSION
        "Reactive multicomponent boundary and shock-tube tests passed"
      FAIL_REGULAR_EXPRESSION
        "floating-point exceptions are signalling;ERROR STOP"
  )
  add_test(
    NAME nse_multicomponent_reactive_boundaries_input
    COMMAND $<TARGET_FILE:nse_solver>
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_multicomponent_reactive_boundaries.dat"
  )
  set_tests_properties(
    nse_multicomponent_reactive_boundaries_input
    PROPERTIES
      PASS_REGULAR_EXPRESSION
        "Reactive multicomponent Navier-Stokes calculation completed successfully"
      FAIL_REGULAR_EXPRESSION
        "floating-point exceptions are signalling;ERROR STOP"
  )
endif()

if(BUILD_TESTING AND NSE_TEST_SUITE STREQUAL "multicomponent_viscous")
  add_executable(
    nse_test_multicomponent_viscous
    tests/test_multicomponent_viscous.f90
  )
  set_target_properties(
    nse_test_multicomponent_viscous
    PROPERTIES
      Fortran_MODULE_DIRECTORY
        "${NSE_MODULE_DIR}/multicomponent_viscous_test"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
  )
  target_link_libraries(nse_test_multicomponent_viscous PRIVATE NSE::core)
  nse_apply_fortran_options(nse_test_multicomponent_viscous)
  add_test(
    NAME nse_multicomponent_viscous
    COMMAND $<TARGET_FILE:nse_test_multicomponent_viscous>
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_multicomponent_viscous.dat"
  )
  set_tests_properties(
    nse_multicomponent_viscous
    PROPERTIES
      PASS_REGULAR_EXPRESSION
        "Multicomponent viscous transport tests passed"
      FAIL_REGULAR_EXPRESSION
        "floating-point exceptions are signalling;ERROR STOP"
  )
  add_test(
    NAME nse_multicomponent_viscous_input
    COMMAND $<TARGET_FILE:nse_solver>
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_multicomponent_viscous.dat"
  )
  set_tests_properties(
    nse_multicomponent_viscous_input
    PROPERTIES
      PASS_REGULAR_EXPRESSION
        "Multicomponent Navier-Stokes calculation completed successfully"
      FAIL_REGULAR_EXPRESSION
        "floating-point exceptions are signalling;ERROR STOP"
  )
endif()

if(BUILD_TESTING AND NSE_TEST_SUITE STREQUAL "multicomponent_euler")
  add_executable(
    nse_test_multicomponent_euler
    tests/test_multicomponent_euler.f90
  )
  set_target_properties(
    nse_test_multicomponent_euler
    PROPERTIES
      Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}/multicomponent_euler_test"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
  )
  target_link_libraries(nse_test_multicomponent_euler PRIVATE NSE::core)
  nse_apply_fortran_options(nse_test_multicomponent_euler)
  add_test(
    NAME nse_multicomponent_euler
    COMMAND $<TARGET_FILE:nse_test_multicomponent_euler>
  )
  set_tests_properties(
    nse_multicomponent_euler
    PROPERTIES
      PASS_REGULAR_EXPRESSION "Multicomponent Euler tests passed"
  )
  add_test(
    NAME nse_multicomponent_euler_input
    COMMAND $<TARGET_FILE:nse_solver>
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_multicomponent_euler.dat"
  )
  set_tests_properties(
    nse_multicomponent_euler_input
    PROPERTIES
      PASS_REGULAR_EXPRESSION
        "Multicomponent Euler calculation completed successfully"
      FAIL_REGULAR_EXPRESSION
        "floating-point exceptions are signalling;ERROR STOP"
  )
endif()

if(BUILD_TESTING AND
   NSE_TEST_SUITE STREQUAL "multicomponent_thermally_perfect")
  add_executable(
    nse_test_multicomponent_thermally_perfect
    tests/test_multicomponent_thermally_perfect.f90
  )
  set_target_properties(
    nse_test_multicomponent_thermally_perfect
    PROPERTIES
      Fortran_MODULE_DIRECTORY
        "${NSE_MODULE_DIR}/multicomponent_thermally_perfect_test"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
  )
  target_link_libraries(
    nse_test_multicomponent_thermally_perfect PRIVATE NSE::core
  )
  nse_apply_fortran_options(nse_test_multicomponent_thermally_perfect)
  add_test(
    NAME nse_multicomponent_thermally_perfect
    COMMAND $<TARGET_FILE:nse_test_multicomponent_thermally_perfect>
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_multicomponent_thermally_perfect.dat"
  )
  set_tests_properties(
    nse_multicomponent_thermally_perfect
    PROPERTIES
      PASS_REGULAR_EXPRESSION
        "Thermally-perfect multicomponent tests passed"
      FAIL_REGULAR_EXPRESSION
        "floating-point exceptions are signalling;ERROR STOP"
  )
  add_test(
    NAME nse_multicomponent_thermally_perfect_input
    COMMAND $<TARGET_FILE:nse_solver>
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_multicomponent_thermally_perfect.dat"
  )
  set_tests_properties(
    nse_multicomponent_thermally_perfect_input
    PROPERTIES
      PASS_REGULAR_EXPRESSION
        "Multicomponent Euler calculation completed successfully"
      FAIL_REGULAR_EXPRESSION
        "floating-point exceptions are signalling;ERROR STOP"
  )
endif()

if(BUILD_TESTING AND NSE_TEST_SUITE STREQUAL "multicomponent_passive_scalar")
  add_executable(
    nse_test_multicomponent_passive_scalar
    tests/test_multicomponent_passive_scalar.f90
  )
  set_target_properties(
    nse_test_multicomponent_passive_scalar
    PROPERTIES
      Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}/multicomponent_passive_scalar_test"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
  )
  target_link_libraries(nse_test_multicomponent_passive_scalar PRIVATE NSE::core)
  nse_apply_fortran_options(nse_test_multicomponent_passive_scalar)
  add_test(
    NAME nse_multicomponent_passive_scalar
    COMMAND $<TARGET_FILE:nse_test_multicomponent_passive_scalar>
  )
  set_tests_properties(
    nse_multicomponent_passive_scalar
    PROPERTIES
      PASS_REGULAR_EXPRESSION
        "Multicomponent passive-scalar tests passed"
  )
  add_test(
    NAME nse_multicomponent_passive_scalar_input
    COMMAND $<TARGET_FILE:nse_solver>
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_multicomponent_passive_scalar.dat"
  )
  set_tests_properties(
    nse_multicomponent_passive_scalar_input
    PROPERTIES
      PASS_REGULAR_EXPRESSION
        "Passive-scalar calculation completed successfully"
      FAIL_REGULAR_EXPRESSION
        "floating-point exceptions are signalling;ERROR STOP"
  )
endif()

if(BUILD_TESTING AND NSE_MC_CUDA)
  add_executable(nse_test_mc_cuda_species tests/test_multicomponent_cuda_species.f90)
  set_target_properties(nse_test_mc_cuda_species PROPERTIES
    Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}/mc_cuda_species_test"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
  target_link_libraries(nse_test_mc_cuda_species PRIVATE NSE::core)
  nse_apply_fortran_options(nse_test_mc_cuda_species)
  add_test(NAME nse_mc_cuda_species COMMAND $<TARGET_FILE:nse_test_mc_cuda_species>)
  set_tests_properties(nse_mc_cuda_species PROPERTIES
    PASS_REGULAR_EXPRESSION "CUDA variable species comparison passed"
    FAIL_REGULAR_EXPRESSION "ERROR STOP")
  add_executable(nse_test_mc_cuda tests/test_multicomponent_cuda.f90)
  set_target_properties(nse_test_mc_cuda PROPERTIES
    Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}/mc_cuda_test"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
  target_link_libraries(nse_test_mc_cuda PRIVATE NSE::core)
  nse_apply_fortran_options(nse_test_mc_cuda)
  add_test(NAME nse_mc_cuda_cpu_comparison COMMAND $<TARGET_FILE:nse_test_mc_cuda>
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_multicomponent_reactive_boundaries.dat")
  foreach(failure IN ITEMS bad_device bad_state bad_temperature substep_limit)
    add_test(NAME nse_mc_cuda_reject_${failure} COMMAND $<TARGET_FILE:nse_test_mc_cuda>
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_multicomponent_reactive_boundaries.dat" ${failure})
    set_tests_properties(nse_mc_cuda_reject_${failure} PROPERTIES WILL_FAIL TRUE TIMEOUT 30)
  endforeach()
  set_tests_properties(nse_mc_cuda_cpu_comparison PROPERTIES
    PASS_REGULAR_EXPRESSION "Resident CUDA CPU comparison passed"
    FAIL_REGULAR_EXPRESSION "ERROR STOP;MC CUDA .*error")
  if(NSE_USE_MPI)
    add_executable(nse_test_mc_cuda_mpi tests/test_multicomponent_cuda_mpi.f90)
    set_target_properties(nse_test_mc_cuda_mpi PROPERTIES
      Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}/mc_cuda_mpi_test"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
    target_link_libraries(nse_test_mc_cuda_mpi PRIVATE NSE::core)
    nse_apply_fortran_options(nse_test_mc_cuda_mpi)
    foreach(np IN ITEMS 1 2 3 4)
      add_test(NAME nse_mc_cuda_pencil_${np} COMMAND "${MPIEXEC_EXECUTABLE}"
        ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:nse_test_mc_cuda_mpi> ${MPIEXEC_POSTFLAGS}
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_multicomponent_reactive_boundaries.dat")
      set_tests_properties(nse_mc_cuda_pencil_${np} PROPERTIES
        TIMEOUT 120 RUN_SERIAL TRUE
        PASS_REGULAR_EXPRESSION "Resident CUDA MPI pencil comparison passed")
    endforeach()
  endif()
endif()

if(NSE_TEST_SUITE STREQUAL "legacy")
if(BUILD_TESTING AND NSE_USE_MPI AND
   (NSE_GPU_BACKEND STREQUAL "none" OR NSE_ENABLE_MULTI_GPU_TESTS))
  add_test(NAME nse_mpi_positivity_retry
    COMMAND "${MPIEXEC_EXECUTABLE}" ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:nse_solver> ${MPIEXEC_POSTFLAGS}
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_positivity_retry.dat")
  set_tests_properties(nse_mpi_positivity_retry PROPERTIES
    TIMEOUT 60)
endif()
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
  add_executable(nse_test_fluctuating tests/test_fluctuating.f90)
  target_link_libraries(nse_test_fluctuating PRIVATE NSE::core)
  set_target_properties(nse_test_fluctuating PROPERTIES
    Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}/fh_test"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
  nse_apply_fortran_options(nse_test_fluctuating)
  add_test(NAME nse_fluctuating COMMAND $<TARGET_FILE:nse_test_fluctuating>)
  if(NSE_VISCOUS_SCHEME STREQUAL "central6")
    file(READ "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_fluctuating.dat" fh_smoke_input)
    foreach(scheme IN ITEMS keep2 keep6 weno5z_roe hybrid)
      string(REPLACE "convective_scheme='keep6'" "convective_scheme='${scheme}'" fh_input "${fh_smoke_input}")
      string(REPLACE "use_fixed_dt=.true." "use_fixed_dt=.false." fh_input "${fh_input}")
      string(REPLACE "t_max=0.0003" "t_max=0.3" fh_input "${fh_input}")
      file(GENERATE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/fh_${scheme}_variable.dat" CONTENT "${fh_input}")
      if(NSE_USE_MPI)
        add_test(NAME nse_fh_${scheme}_variable COMMAND "${MPIEXEC_EXECUTABLE}"
          ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:nse_solver> ${MPIEXEC_POSTFLAGS}
          "${CMAKE_CURRENT_BINARY_DIR}/fh_${scheme}_variable.dat")
      else()
        add_test(NAME nse_fh_${scheme}_variable COMMAND $<TARGET_FILE:nse_solver>
          "${CMAKE_CURRENT_BINARY_DIR}/fh_${scheme}_variable.dat")
      endif()
      set_tests_properties(nse_fh_${scheme}_variable PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=2;NSE_CUDA_DEVICE_POLICY=fixed"
        PASS_REGULAR_EXPRESSION "completed successfully" TIMEOUT 60)
    endforeach()
  endif()
  set_tests_properties(nse_fluctuating PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2" TIMEOUT 120)
  foreach(mode IN ITEMS adaptive weno)
    add_test(NAME nse_fh_accept_${mode} COMMAND $<TARGET_FILE:nse_test_fluctuating> ${mode})
  endforeach()
  foreach(failure IN ITEMS dt transport boundary inviscid negative_beta seed)
    add_test(NAME nse_fh_reject_${failure} COMMAND $<TARGET_FILE:nse_test_fluctuating> ${failure})
    set_tests_properties(nse_fh_reject_${failure} PROPERTIES WILL_FAIL TRUE TIMEOUT 30)
  endforeach()
  if(NSE_GPU_BACKEND STREQUAL "none" AND NSE_VISCOUS_SCHEME STREQUAL "central6")
    if(NOT NSE_USE_MPI)
    add_test(NAME nse_fluctuating_smoke COMMAND $<TARGET_FILE:nse_solver>
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_fluctuating.dat")
    set_tests_properties(nse_fluctuating_smoke PROPERTIES
      ENVIRONMENT "OMP_NUM_THREADS=2" PASS_REGULAR_EXPRESSION "completed successfully" TIMEOUT 60)
    endif()
    if(NSE_USE_MPI)
      add_test(NAME nse_fluctuating_mpi_smoke COMMAND "${MPIEXEC_EXECUTABLE}"
        ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:nse_solver> ${MPIEXEC_POSTFLAGS}
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_fluctuating.dat")
      set_tests_properties(nse_fluctuating_mpi_smoke PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=2" PASS_REGULAR_EXPRESSION "completed successfully" TIMEOUT 60)
    endif()
  endif()
  if(NSE_GPU_BACKEND STREQUAL "cuda" AND NOT NSE_USE_MPI)
    add_test(NAME nse_fluctuating_cuda_smoke COMMAND $<TARGET_FILE:nse_solver>
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_fluctuating.dat")
    set_tests_properties(nse_fluctuating_cuda_smoke PROPERTIES
      PASS_REGULAR_EXPRESSION "completed successfully" TIMEOUT 30)
  endif()
  add_executable(
    nse_test_imported_turbulence
    tests/test_imported_turbulence.f90
  )
  set_target_properties(
    nse_test_imported_turbulence
    PROPERTIES
      Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}/imported_turbulence_test"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
  )
  target_link_libraries(nse_test_imported_turbulence PRIVATE NSE::core)
  nse_apply_fortran_options(nse_test_imported_turbulence)
  add_test(
    NAME nse_imported_turbulence
    COMMAND $<TARGET_FILE:nse_test_imported_turbulence>
  )
  set_tests_properties(
    nse_imported_turbulence
    PROPERTIES
      PASS_REGULAR_EXPRESSION "Imported turbulence initialization tests passed"
  )

  add_executable(
    nse_test_hit_parameter_resolution
    tests/test_hit_parameter_resolution.f90
  )
  set_target_properties(
    nse_test_hit_parameter_resolution
    PROPERTIES
      Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}/hit_parameter_test"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
  )
  target_link_libraries(nse_test_hit_parameter_resolution PRIVATE NSE::core)
  nse_apply_fortran_options(nse_test_hit_parameter_resolution)
  add_test(
    NAME nse_hit_parameter_resolution
    COMMAND $<TARGET_FILE:nse_test_hit_parameter_resolution>
  )
  set_tests_properties(
    nse_hit_parameter_resolution
    PROPERTIES PASS_REGULAR_EXPRESSION "HIT parameter resolution test passed"
  )

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
  if(NOT NSE_USE_MPI)
    add_executable(nse_test_positivity tests/test_nse_positivity.f90)
    set_target_properties(nse_test_positivity PROPERTIES
      Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
    target_link_libraries(nse_test_positivity PRIVATE NSE::core)
    nse_apply_fortran_options(nse_test_positivity)
    add_test(NAME nse_positivity COMMAND $<TARGET_FILE:nse_test_positivity>)
    add_test(NAME nse_positivity_fixed_reject COMMAND $<TARGET_FILE:nse_test_positivity> fixed-reject)
    set_tests_properties(nse_positivity_fixed_reject PROPERTIES WILL_FAIL TRUE)
  endif()
  add_executable(nse_test_hit_isotropy_math tests/test_hit_isotropy_math.f90)
  set_target_properties(
    nse_test_hit_isotropy_math
    PROPERTIES
      Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}/hit_isotropy_test"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
  )
  target_link_libraries(nse_test_hit_isotropy_math PRIVATE NSE::core)
  nse_apply_fortran_options(nse_test_hit_isotropy_math)
  add_test(NAME nse_hit_isotropy_math
    COMMAND $<TARGET_FILE:nse_test_hit_isotropy_math>)
  set_tests_properties(
    nse_hit_isotropy_math
    PROPERTIES PASS_REGULAR_EXPRESSION "HIT isotropy tensor test passed"
  )

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

  add_executable(nse_test_weno5z_roe tests/test_weno5z_roe.f90)
  set_target_properties(
    nse_test_weno5z_roe
    PROPERTIES
      Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}/weno5z_roe_test"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
  )
  target_link_libraries(nse_test_weno5z_roe PRIVATE NSE::core)
  nse_apply_fortran_options(nse_test_weno5z_roe)
  add_test(NAME nse_weno5z_roe
    COMMAND $<TARGET_FILE:nse_test_weno5z_roe>)
  set_tests_properties(
    nse_weno5z_roe
    PROPERTIES
      PASS_REGULAR_EXPRESSION "WENO5-Z/Roe tests passed"
  )

  add_executable(nse_test_hybrid_flux tests/test_hybrid_flux.f90)
  set_target_properties(
    nse_test_hybrid_flux
    PROPERTIES
      Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}/hybrid_flux_test"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
  )
  target_link_libraries(nse_test_hybrid_flux PRIVATE NSE::core)
  nse_apply_fortran_options(nse_test_hybrid_flux)
  add_test(NAME nse_hybrid_flux
    COMMAND $<TARGET_FILE:nse_test_hybrid_flux>)
  set_tests_properties(
    nse_hybrid_flux
    PROPERTIES PASS_REGULAR_EXPRESSION "Hybrid KEEP/WENO flux tests passed"
  )

  if(NSE_USE_MPI)
    add_test(
      NAME nse_weno5z_roe_smoke
      COMMAND mpiexec -n 4 $<TARGET_FILE:nse_solver>
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_tgv_weno5z_roe_small.dat"
    )
  else()
    add_test(
      NAME nse_weno5z_roe_smoke
      COMMAND $<TARGET_FILE:nse_solver>
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_tgv_weno5z_roe_small.dat"
    )
  endif()
  set_tests_properties(
    nse_weno5z_roe_smoke
    PROPERTIES
      PASS_REGULAR_EXPRESSION "NSE calculation completed successfully:"
      FAIL_REGULAR_EXPRESSION "floating-point exceptions are signalling"
  )
  if(NSE_USE_MPI)
    add_test(
      NAME nse_hybrid_smoke
      COMMAND mpiexec -n 4 $<TARGET_FILE:nse_solver>
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_tgv_hybrid_small.dat"
    )
  else()
    add_test(
      NAME nse_hybrid_smoke
      COMMAND $<TARGET_FILE:nse_solver>
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_tgv_hybrid_small.dat"
    )
  endif()
  set_tests_properties(
    nse_hybrid_smoke
    PROPERTIES
      PASS_REGULAR_EXPRESSION "NSE calculation completed successfully:"
      FAIL_REGULAR_EXPRESSION "floating-point exceptions are signalling"
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
  if(NSE_BOUNDARY_SCHEME STREQUAL "runtime")
    add_executable(
      nse_test_boundary_non_reflecting
      tests/test_boundary_non_reflecting.f90
    )
    set_target_properties(
      nse_test_boundary_non_reflecting
      PROPERTIES
        Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}/boundary_non_reflecting_test"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
    )
    target_link_libraries(nse_test_boundary_non_reflecting PRIVATE NSE::core)
    nse_apply_fortran_options(nse_test_boundary_non_reflecting)
    add_test(
      NAME nse_boundary_non_reflecting
      COMMAND mpiexec -n 4 $<TARGET_FILE:nse_test_boundary_non_reflecting>
    )
    set_tests_properties(
      nse_boundary_non_reflecting
      PROPERTIES
        PASS_REGULAR_EXPRESSION "non-reflecting mixed boundary tests passed"
    )
    add_executable(
      nse_test_boundary_reflective
      tests/test_boundary_reflective.f90
    )
    set_target_properties(
      nse_test_boundary_reflective
      PROPERTIES
        Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}/boundary_reflective_test"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
    )
    target_link_libraries(nse_test_boundary_reflective PRIVATE NSE::core)
    nse_apply_fortran_options(nse_test_boundary_reflective)
    add_test(
      NAME nse_boundary_reflective
      COMMAND mpiexec -n 4 $<TARGET_FILE:nse_test_boundary_reflective>
    )
    set_tests_properties(
      nse_boundary_reflective
      PROPERTIES
        PASS_REGULAR_EXPRESSION
          "reflective boundary face/edge/corner test passed"
    )
    add_executable(
      nse_test_boundary_dirichlet
      tests/test_boundary_dirichlet.f90
    )
    set_target_properties(
      nse_test_boundary_dirichlet
      PROPERTIES
        Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}/boundary_dirichlet_test"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
    )
    target_link_libraries(nse_test_boundary_dirichlet PRIVATE NSE::core)
    nse_apply_fortran_options(nse_test_boundary_dirichlet)
    add_test(
      NAME nse_boundary_dirichlet
      COMMAND mpiexec -n 4 $<TARGET_FILE:nse_test_boundary_dirichlet>
    )
    set_tests_properties(
      nse_boundary_dirichlet
      PROPERTIES
        PASS_REGULAR_EXPRESSION
          "Dirichlet boundary face/edge/corner test passed"
    )
    add_test(
      NAME nse_non_reflecting_smoke
      COMMAND mpiexec -n 4 $<TARGET_FILE:nse_solver>
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_tgv_non_reflecting_small.dat"
    )
    set_tests_properties(
      nse_non_reflecting_smoke
      PROPERTIES
        PASS_REGULAR_EXPRESSION "NSE calculation completed successfully:"
        FAIL_REGULAR_EXPRESSION "floating-point exceptions are signalling"
    )
    add_test(
      NAME nse_reflective_smoke
      COMMAND mpiexec -n 4 $<TARGET_FILE:nse_solver>
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_tgv_reflective_small.dat"
    )
    set_tests_properties(
      nse_reflective_smoke
      PROPERTIES
        PASS_REGULAR_EXPRESSION "NSE calculation completed successfully:"
        FAIL_REGULAR_EXPRESSION "floating-point exceptions are signalling"
    )
    add_test(
      NAME nse_dirichlet_smoke
      COMMAND mpiexec -n 4 $<TARGET_FILE:nse_solver>
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_tgv_dirichlet_small.dat"
    )
    set_tests_properties(
      nse_dirichlet_smoke
      PROPERTIES
        PASS_REGULAR_EXPRESSION "NSE calculation completed successfully:"
        FAIL_REGULAR_EXPRESSION "floating-point exceptions are signalling"
    )
  endif()
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
   NSE_INIT_FFT_BACKEND STREQUAL "2decomp_fftw")
  add_test(
    NAME nse_hit_targets_2decomp_smoke
    COMMAND mpiexec -n 4 $<TARGET_FILE:nse_solver>
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_hit_targets_small.dat"
  )
  set_tests_properties(
    nse_hit_targets_2decomp_smoke
    PROPERTIES
      PASS_REGULAR_EXPRESSION
        "hit_integral_reynolds =;NSE calculation completed successfully:"
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
        "floating-point exceptions are signalling;HIT low-wavenumber isotropy did not converge"
  )
endif()
if(BUILD_TESTING AND NSE_GPU_BACKEND STREQUAL "cuda" AND NOT NSE_USE_MPI)
  if(NSE_FORCING_FFT_BACKEND STREQUAL "cufft")
    add_executable(nse_test_cuda_positivity tests/test_cuda_positivity.f90)
    set_target_properties(nse_test_cuda_positivity PROPERTIES
      Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
    target_link_libraries(nse_test_cuda_positivity PRIVATE NSE::core)
    nse_apply_fortran_options(nse_test_cuda_positivity)
    add_test(NAME nse_cuda_positivity COMMAND $<TARGET_FILE:nse_test_cuda_positivity>)
    add_test(NAME nse_cuda_positivity_fixed_reject
      COMMAND $<TARGET_FILE:nse_test_cuda_positivity> fixed-reject)
    set_tests_properties(nse_cuda_positivity_fixed_reject PROPERTIES WILL_FAIL TRUE)
  endif()
  add_executable(
    nse_test_cuda_boundary_dirichlet
    tests/test_cuda_boundary_dirichlet.f90
  )
  set_target_properties(
    nse_test_cuda_boundary_dirichlet
    PROPERTIES
      Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
  )
  target_link_libraries(nse_test_cuda_boundary_dirichlet PRIVATE NSE::core)
  nse_apply_fortran_options(nse_test_cuda_boundary_dirichlet)
  add_test(
    NAME nse_cuda_boundary_dirichlet_compare
    COMMAND $<TARGET_FILE:nse_test_cuda_boundary_dirichlet>
  )
  set_tests_properties(
    nse_cuda_boundary_dirichlet_compare
    PROPERTIES
      PASS_REGULAR_EXPRESSION "CUDA Dirichlet boundary comparison passed"
      FAIL_REGULAR_EXPRESSION "CUDA backend failure;floating-point exceptions are signalling"
  )
  add_executable(
    nse_test_cuda_boundary_reflective
    tests/test_cuda_boundary_reflective.f90
  )
  set_target_properties(
    nse_test_cuda_boundary_reflective
    PROPERTIES
      Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
  )
  target_link_libraries(nse_test_cuda_boundary_reflective PRIVATE NSE::core)
  nse_apply_fortran_options(nse_test_cuda_boundary_reflective)
  add_test(
    NAME nse_cuda_boundary_reflective_compare
    COMMAND $<TARGET_FILE:nse_test_cuda_boundary_reflective>
  )
  set_tests_properties(
    nse_cuda_boundary_reflective_compare
    PROPERTIES
      PASS_REGULAR_EXPRESSION "CUDA reflective boundary comparison passed"
      FAIL_REGULAR_EXPRESSION "CUDA backend failure;floating-point exceptions are signalling"
  )
  add_executable(
    nse_test_cuda_boundary_non_reflecting
    tests/test_cuda_boundary_non_reflecting.f90
  )
  set_target_properties(
    nse_test_cuda_boundary_non_reflecting
    PROPERTIES
      Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin"
  )
  target_link_libraries(nse_test_cuda_boundary_non_reflecting PRIVATE NSE::core)
  nse_apply_fortran_options(nse_test_cuda_boundary_non_reflecting)
  add_test(
    NAME nse_cuda_boundary_non_reflecting_compare
    COMMAND $<TARGET_FILE:nse_test_cuda_boundary_non_reflecting>
  )
  set_tests_properties(
    nse_cuda_boundary_non_reflecting_compare
    PROPERTIES
      PASS_REGULAR_EXPRESSION "CUDA non-reflecting boundary comparison passed"
      FAIL_REGULAR_EXPRESSION "CUDA backend failure;floating-point exceptions are signalling"
  )
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
    src/numerics/reconstruction/mod_reconstruction_weno5z.f90
    src/numerics/riemann/mod_riemann_roe.f90
    src/numerics/convective/mod_convective_keep.f90
    src/numerics/convective/mod_convective_weno5z_roe.f90
    src/numerics/convective/mod_convective_leaf_registry.f90
    src/numerics/convective/mod_convective_hybrid.f90
    src/numerics/convective/mod_convective_dispatch.f90
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
    PROPERTIES PASS_REGULAR_EXPRESSION "CUDA convective comparison passed"
  )
  add_test(NAME nse_cuda_keep2_compare
    COMMAND $<TARGET_FILE:nse_test_cuda_keep_compare> 2)
  set_tests_properties(
    nse_cuda_keep2_compare
    PROPERTIES PASS_REGULAR_EXPRESSION "CUDA convective comparison passed"
  )
  add_test(NAME nse_cuda_weno5z_roe_compare
    COMMAND $<TARGET_FILE:nse_test_cuda_keep_compare> 5)
  set_tests_properties(
    nse_cuda_weno5z_roe_compare
    PROPERTIES PASS_REGULAR_EXPRESSION "CUDA convective comparison passed"
  )
  add_test(NAME nse_cuda_hybrid_compare
    COMMAND $<TARGET_FILE:nse_test_cuda_keep_compare> 7)
  set_tests_properties(
    nse_cuda_hybrid_compare
    PROPERTIES PASS_REGULAR_EXPRESSION "CUDA convective comparison passed"
  )
  add_executable(
    nse_test_cuda_viscous_compare
    tests/test_cuda_viscous_compare.f90
    src/numerics/reconstruction/mod_reconstruction_weno5z.f90
    src/numerics/riemann/mod_riemann_roe.f90
    src/numerics/convective/mod_convective_keep.f90
    src/numerics/convective/mod_convective_weno5z_roe.f90
    src/numerics/convective/mod_convective_leaf_registry.f90
    src/numerics/convective/mod_convective_hybrid.f90
    src/numerics/convective/mod_convective_dispatch.f90
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
  foreach(mode IN ITEMS fluctuating fh_zero fh_keep2 fh_weno fh_hybrid)
    add_test(NAME nse_cuda_${mode}_compare
      COMMAND $<TARGET_FILE:nse_test_cuda_viscous_compare> ${mode})
    set_tests_properties(nse_cuda_${mode}_compare PROPERTIES
      PASS_REGULAR_EXPRESSION "Resident CUDA LLNS comparison passed" TIMEOUT 60)
  endforeach()
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
  add_test(
    NAME nse_cuda_non_reflecting_smoke
    COMMAND $<TARGET_FILE:nse_solver>
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_tgv_non_reflecting_small.dat"
  )
  set_tests_properties(
    nse_cuda_non_reflecting_smoke
    PROPERTIES
      PASS_REGULAR_EXPRESSION "NSE calculation completed successfully:"
      FAIL_REGULAR_EXPRESSION "CUDA backend failure;floating-point exceptions are signalling"
  )
  add_test(
    NAME nse_cuda_reflective_smoke
    COMMAND $<TARGET_FILE:nse_solver>
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_tgv_reflective_small.dat"
  )
  set_tests_properties(
    nse_cuda_reflective_smoke
    PROPERTIES
      PASS_REGULAR_EXPRESSION "NSE calculation completed successfully:"
      FAIL_REGULAR_EXPRESSION "CUDA backend failure;floating-point exceptions are signalling"
  )
  add_test(
    NAME nse_cuda_dirichlet_smoke
    COMMAND $<TARGET_FILE:nse_solver>
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_tgv_dirichlet_small.dat"
  )
  set_tests_properties(
    nse_cuda_dirichlet_smoke
    PROPERTIES
      PASS_REGULAR_EXPRESSION "NSE calculation completed successfully:"
      FAIL_REGULAR_EXPRESSION "CUDA backend failure;floating-point exceptions are signalling"
  )
  add_test(
    NAME nse_cuda_weno5z_roe_smoke
    COMMAND $<TARGET_FILE:nse_solver>
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_tgv_cuda_weno5z_roe_small.dat"
  )
  set_tests_properties(
    nse_cuda_weno5z_roe_smoke
    PROPERTIES
      PASS_REGULAR_EXPRESSION "NSE calculation completed successfully:"
      FAIL_REGULAR_EXPRESSION "CUDA backend failure;floating-point exceptions are signalling"
  )
  add_test(
    NAME nse_cuda_hybrid_smoke
    COMMAND $<TARGET_FILE:nse_solver>
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_tgv_hybrid_small.dat"
  )
  set_tests_properties(
    nse_cuda_hybrid_smoke
    PROPERTIES
      PASS_REGULAR_EXPRESSION "NSE calculation completed successfully:"
      FAIL_REGULAR_EXPRESSION "CUDA backend failure;floating-point exceptions are signalling"
  )
  if(NSE_INIT_FFT_BACKEND STREQUAL "cufft")
    add_test(
      NAME nse_cuda_hit_smoke
      COMMAND $<TARGET_FILE:nse_solver>
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_hit_cuda_small.dat"
    )
    set_tests_properties(
      nse_cuda_hit_smoke
      PROPERTIES
        PASS_REGULAR_EXPRESSION
          "Single-GPU spectral HIT initialization"
        FAIL_REGULAR_EXPRESSION
          "CUDA HIT initialization error;CUDA backend failure"
    )
    add_test(
      NAME nse_cuda_hit_targets_smoke
      COMMAND $<TARGET_FILE:nse_solver>
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_hit_targets_small.dat"
    )
    set_tests_properties(
      nse_cuda_hit_targets_smoke
      PROPERTIES
        PASS_REGULAR_EXPRESSION
          "hit_integral_reynolds =;NSE calculation completed successfully:"
        FAIL_REGULAR_EXPRESSION
          "CUDA HIT initialization error;CUDA backend failure"
    )
  endif()
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

if(BUILD_TESTING AND NSE_GPU_BACKEND STREQUAL "cuda" AND NSE_USE_MPI AND
   NSE_ENABLE_MULTI_GPU_TESTS)
  add_executable(nse_test_cuda_fluctuating_mpi tests/test_cuda_fluctuating_mpi.f90)
  set_target_properties(nse_test_cuda_fluctuating_mpi PROPERTIES
    Fortran_MODULE_DIRECTORY "${NSE_MODULE_DIR}/fh_mpi_test"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
  target_link_libraries(nse_test_cuda_fluctuating_mpi PRIVATE NSE::core)
  nse_apply_fortran_options(nse_test_cuda_fluctuating_mpi)
  foreach(np IN ITEMS 2 4)
    add_test(NAME nse_cuda_fh_mpi_${np} COMMAND "${MPIEXEC_EXECUTABLE}"
      ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:nse_test_cuda_fluctuating_mpi> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(nse_cuda_fh_mpi_${np} PROPERTIES
      PASS_REGULAR_EXPRESSION "MPI resident LLNS matches" TIMEOUT 120 RUN_SERIAL TRUE)
    if(NSE_ENABLE_CUDA_AWARE_MPI)
      set_tests_properties(nse_cuda_fh_mpi_${np} PROPERTIES ENVIRONMENT "NSE_CUDA_MPI_TRANSPORT=device")
    endif()
  endforeach()
  if(NOT NSE_ENABLE_CUDA_AWARE_MPI)
    foreach(mode IN ITEMS device invalid)
      add_test(NAME nse_cuda_mpi_transport_reject_${mode}
        COMMAND "${CMAKE_COMMAND}" "-DLAUNCHER=${MPIEXEC_EXECUTABLE}"
          "-DNP_FLAG=${MPIEXEC_NUMPROC_FLAG}" "-DSOLVER=$<TARGET_FILE:nse_solver>"
          "-DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/tests/input_fluctuating.dat" "-DMODE=${mode}"
          -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/check_cuda_mpi_transport_reject.cmake")
      set_tests_properties(nse_cuda_mpi_transport_reject_${mode} PROPERTIES TIMEOUT 60 RUN_SERIAL TRUE)
    endforeach()
  endif()
  add_test(NAME nse_cuda_fh_mpi_smoke COMMAND "${MPIEXEC_EXECUTABLE}"
    ${MPIEXEC_NUMPROC_FLAG} ${NSE_MULTI_GPU_TEST_PROCESSES} ${MPIEXEC_PREFLAGS}
    $<TARGET_FILE:nse_solver> ${MPIEXEC_POSTFLAGS}
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_fluctuating.dat")
  set_tests_properties(nse_cuda_fh_mpi_smoke PROPERTIES
    PASS_REGULAR_EXPRESSION "completed successfully" TIMEOUT 60)
  add_test(
    NAME nse_mpi_cuda_tgv_smoke
    COMMAND "${MPIEXEC_EXECUTABLE}" ${MPIEXEC_NUMPROC_FLAG}
      ${NSE_MULTI_GPU_TEST_PROCESSES} ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:nse_solver> ${MPIEXEC_POSTFLAGS}
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_tgv_cuda_small.dat"
  )
  set_tests_properties(
    nse_mpi_cuda_tgv_smoke
    PROPERTIES
      PASS_REGULAR_EXPRESSION
        "NSE MPI\\+CUDA calculation completed successfully:"
      FAIL_REGULAR_EXPRESSION
        "CUDA backend failure;floating-point exceptions are signalling"
  )
  add_test(
    NAME nse_mpi_cuda_non_reflecting_smoke
    COMMAND "${MPIEXEC_EXECUTABLE}" ${MPIEXEC_NUMPROC_FLAG}
      ${NSE_MULTI_GPU_TEST_PROCESSES} ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:nse_solver> ${MPIEXEC_POSTFLAGS}
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_tgv_non_reflecting_all_faces_small.dat"
  )
  set_tests_properties(
    nse_mpi_cuda_non_reflecting_smoke
    PROPERTIES
      PASS_REGULAR_EXPRESSION
        "NSE MPI\\+CUDA calculation completed successfully:"
      FAIL_REGULAR_EXPRESSION
        "CUDA backend failure;floating-point exceptions are signalling"
  )
  add_test(
    NAME nse_mpi_cuda_reflective_smoke
    COMMAND "${MPIEXEC_EXECUTABLE}" ${MPIEXEC_NUMPROC_FLAG}
      ${NSE_MULTI_GPU_TEST_PROCESSES} ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:nse_solver> ${MPIEXEC_POSTFLAGS}
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_tgv_reflective_small.dat"
  )
  set_tests_properties(
    nse_mpi_cuda_reflective_smoke
    PROPERTIES
      PASS_REGULAR_EXPRESSION
        "NSE MPI\\+CUDA calculation completed successfully:"
      FAIL_REGULAR_EXPRESSION
        "CUDA backend failure;floating-point exceptions are signalling"
  )
  add_test(
    NAME nse_mpi_cuda_dirichlet_smoke
    COMMAND "${MPIEXEC_EXECUTABLE}" ${MPIEXEC_NUMPROC_FLAG}
      ${NSE_MULTI_GPU_TEST_PROCESSES} ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:nse_solver> ${MPIEXEC_POSTFLAGS}
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_tgv_dirichlet_small.dat"
  )
  set_tests_properties(
    nse_mpi_cuda_dirichlet_smoke
    PROPERTIES
      PASS_REGULAR_EXPRESSION
        "NSE MPI\\+CUDA calculation completed successfully:"
      FAIL_REGULAR_EXPRESSION
        "CUDA backend failure;floating-point exceptions are signalling"
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
message(STATUS "NSE test suite: ${NSE_TEST_SUITE}")
if(NSE_INIT_FFT_BACKEND STREQUAL "2decomp_fftw")
  message(STATUS "NSE 2DECOMP root: ${NSE_2DECOMP_ROOT}")
endif()
message(STATUS "NSE Fortran modules: ${NSE_MODULE_DIR}")
message(STATUS "NSE executable: ${CMAKE_CURRENT_BINARY_DIR}/bin/${NSE_OUTPUT_NAME}${CMAKE_EXECUTABLE_SUFFIX}")
