# GP3Dのソース選択、外部ライブラリ検出、実行ファイル、CTestを定義する。
# 通常は直接オプションを入力せず、case.yamlから生成されたworkflow.json経由で使用する。
cmake_minimum_required(VERSION 3.20)

project(gp3d LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

# 同じ公開モジュール名を持つ実装から、MPI・FFT・GPUバックエンドを一つずつ選ぶ。
option(USE_MPI "Build with MPI support" OFF)
set(FFT_BACKEND "dft" CACHE STRING "CPU FFT backend: dft or fftw")
set_property(CACHE FFT_BACKEND PROPERTY STRINGS dft fftw)
set(GPU_BACKEND "none" CACHE STRING "GPU backend: none, cuda, or cufftmp")
set_property(CACHE GPU_BACKEND PROPERTY STRINGS none cuda cufftmp)
set(GP3D_CUDA_ARCHITECTURES "" CACHE STRING
  "Optional CUDA architectures, for example 86 for an RTX 3050")
set(CUFFTMP_ROOT "" CACHE PATH "Optional cuFFTMp installation prefix")
set(NVSHMEM_ROOT "" CACHE PATH "Optional NVSHMEM installation prefix")
set(CUFFTMP_API "auto" CACHE STRING "cuFFTMp API generation: auto, modern, or legacy")
set_property(CACHE CUFFTMP_API PROPERTY STRINGS auto modern legacy)

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod")

option(GP3D_ENABLE_WARNINGS "Enable additional Fortran compiler warnings" ON)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(GP3D_DEFAULT_COMMON_FLAGS -fallow-invalid-boz -fallow-argument-mismatch)
  set(GP3D_DEFAULT_WARNING_FLAGS -Wall -Wextra)
  set(GP3D_DEFAULT_DEBUG_FLAGS "")
  set(GP3D_DEFAULT_RELEASE_FLAGS "")
  set(GP3D_DEFAULT_RELWITHDEBINFO_FLAGS "")
else()
  set(GP3D_DEFAULT_COMMON_FLAGS "")
  set(GP3D_DEFAULT_WARNING_FLAGS "")
  set(GP3D_DEFAULT_DEBUG_FLAGS "")
  set(GP3D_DEFAULT_RELEASE_FLAGS "")
  set(GP3D_DEFAULT_RELWITHDEBINFO_FLAGS "")
endif()
set(GP3D_COMMON_FLAGS "${GP3D_DEFAULT_COMMON_FLAGS}" CACHE STRING "Common Fortran flags")
set(GP3D_WARNING_FLAGS "${GP3D_DEFAULT_WARNING_FLAGS}" CACHE STRING "Fortran warning flags")
set(GP3D_DEBUG_FLAGS "${GP3D_DEFAULT_DEBUG_FLAGS}" CACHE STRING "Debug Fortran flags")
set(GP3D_RELEASE_FLAGS "${GP3D_DEFAULT_RELEASE_FLAGS}" CACHE STRING "Release Fortran flags")
set(
  GP3D_RELWITHDEBINFO_FLAGS "${GP3D_DEFAULT_RELWITHDEBINFO_FLAGS}"
  CACHE STRING
  "RelWithDebInfo Fortran flags"
)

foreach(_flag IN LISTS GP3D_COMMON_FLAGS)
  add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:${_flag}>)
endforeach()
if(GP3D_ENABLE_WARNINGS)
  foreach(_flag IN LISTS GP3D_WARNING_FLAGS)
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:${_flag}>)
  endforeach()
endif()
foreach(_flag IN LISTS GP3D_DEBUG_FLAGS)
  add_compile_options($<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:${_flag}>)
endforeach()
foreach(_flag IN LISTS GP3D_RELEASE_FLAGS)
  add_compile_options($<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Release>>:${_flag}>)
endforeach()
foreach(_flag IN LISTS GP3D_RELWITHDEBINFO_FLAGS)
  add_compile_options($<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RelWithDebInfo>>:${_flag}>)
endforeach()

if(NOT GPU_BACKEND STREQUAL "none" AND
   NOT GPU_BACKEND STREQUAL "cuda" AND
   NOT GPU_BACKEND STREQUAL "cufftmp")
  message(FATAL_ERROR
    "Unsupported GPU_BACKEND=${GPU_BACKEND}. Use none, cuda, or cufftmp.")
endif()
if(GPU_BACKEND STREQUAL "cuda" AND USE_MPI)
  message(FATAL_ERROR "The CUDA backend is single-GPU only. Configure with USE_MPI=OFF.")
endif()
if(GPU_BACKEND STREQUAL "cufftmp" AND NOT USE_MPI)
  message(FATAL_ERROR "The cuFFTMp backend requires USE_MPI=ON.")
endif()
if(GPU_BACKEND STREQUAL "cufftmp" AND WIN32)
  message(FATAL_ERROR "cuFFTMp is a Linux backend and cannot be built on Windows.")
endif()
if(NOT CUFFTMP_API STREQUAL "auto" AND
   NOT CUFFTMP_API STREQUAL "modern" AND
   NOT CUFFTMP_API STREQUAL "legacy")
  message(FATAL_ERROR "Unsupported CUFFTMP_API=${CUFFTMP_API}.")
endif()

set(GP3D_GPU_FORTRAN_SOURCE "")
set(GP3D_MAIN_SOURCE src/main/main.f90)
if(GPU_BACKEND STREQUAL "cuda" OR GPU_BACKEND STREQUAL "cufftmp")
  # GNU Fortran and MSVC-hosted nvcc use different executable linkers on Windows.
  # CUDA is built as a bridge library, so compiler checks only need a static archive.
  set(GP3D_FORTRAN_ARCHIVER "${CMAKE_AR}")
  if(GPU_BACKEND STREQUAL "cuda" AND WIN32)
    find_program(GP3D_MSVC_LIBRARIAN NAMES lib.exe REQUIRED)
    get_filename_component(GP3D_MSVC_BIN_DIR
      "${GP3D_MSVC_LIBRARIAN}" DIRECTORY)
    find_program(GP3D_MSVC_LINKER NAMES link.exe
      HINTS "${GP3D_MSVC_BIN_DIR}" NO_DEFAULT_PATH REQUIRED)
    find_program(GP3D_MSVC_RESOURCE_COMPILER NAMES rc.exe REQUIRED)
    set(CMAKE_AR "${GP3D_MSVC_LIBRARIAN}")
  endif()
  set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)
  if(GPU_BACKEND STREQUAL "cufftmp")
    enable_language(CXX)
  endif()
  enable_language(CUDA)
  # The final executable is linked by GNU Fortran.  Keep CUDA runtime linkage
  # explicit so CMake does not append MSVC-style bare library names to it.
  set(CMAKE_CUDA_RUNTIME_LIBRARY None)
  if(GPU_BACKEND STREQUAL "cuda" AND WIN32)
    # CUDA's Windows platform rule expands <CMAKE_LINKER> at generation time.
    # Override the GNU linker inherited from the Fortran-first project.
    set(CMAKE_LINKER "${GP3D_MSVC_LINKER}")
    set(CMAKE_CUDA_HOST_LINK_LAUNCHER "${GP3D_MSVC_LINKER}")
    set(CMAKE_RC_COMPILER "${GP3D_MSVC_RESOURCE_COMPILER}")
    set(CMAKE_SHARED_LIBRARY_CREATE_CUDA_FLAGS "")
    set(CMAKE_CUDA_CREATE_STATIC_LIBRARY
      "\"${GP3D_MSVC_LIBRARIAN}\" /nologo <LINK_FLAGS> /out:<TARGET> <OBJECTS>")
    set(CMAKE_AR "${GP3D_FORTRAN_ARCHIVER}")
  endif()
  find_package(CUDAToolkit REQUIRED)
endif()

if(GPU_BACKEND STREQUAL "cuda")
  set(GP3D_MPI_SOURCE src/mpi/gp3d_mpi_stub.f90)
  set(GP3D_EXECUTABLE_NAME gp3d_cuda)
  set(GP3D_MAIN_SOURCE src/main/main_cuda.f90)
  set(GP3D_GPU_FORTRAN_SOURCE src/gpu/gp3d_gpu_cuda.f90)
elseif(GPU_BACKEND STREQUAL "cufftmp")
  set(GP3D_MPI_SOURCE src/mpi/gp3d_mpi_real.f90)
  set(GP3D_EXECUTABLE_NAME gp3d_cufftmp)
  set(GP3D_MAIN_SOURCE src/main/main_cufftmp.f90)
  set(GP3D_GPU_FORTRAN_SOURCE src/gpu/gp3d_gpu_cufftmp.f90)
elseif(USE_MPI)
  set(GP3D_MPI_SOURCE src/mpi/gp3d_mpi_real.f90)
  set(GP3D_EXECUTABLE_NAME gp3d_mpi)
else()
  set(GP3D_MPI_SOURCE src/mpi/gp3d_mpi_stub.f90)
  set(GP3D_EXECUTABLE_NAME gp3d_sequential)
endif()

if(USE_MPI)
  set(GP3D_FFT_SOURCE src/fft/gp3d_fft_mpi.f90)
  if(FFT_BACKEND STREQUAL "dft")
    set(GP3D_LOCAL_FFT_SOURCE src/fft/gp3d_local_fft_dft.f90)
  elseif(FFT_BACKEND STREQUAL "fftw")
    set(GP3D_LOCAL_FFT_SOURCE src/fft/gp3d_local_fft_fftw.f90)
  else()
    message(FATAL_ERROR "Unsupported FFT_BACKEND=${FFT_BACKEND}. Use dft or fftw.")
  endif()
else()
  set(GP3D_LOCAL_FFT_SOURCE "")
  if(FFT_BACKEND STREQUAL "dft")
    set(GP3D_FFT_SOURCE src/fft/gp3d_fft.f90)
  elseif(FFT_BACKEND STREQUAL "fftw")
    set(GP3D_FFT_SOURCE src/fft/gp3d_fft_fftw.f90)
  else()
    message(FATAL_ERROR "Unsupported FFT_BACKEND=${FFT_BACKEND}. Use dft or fftw.")
  endif()
endif()

set(GP3D_LIBRARY_SOURCES
  src/common/gp3d_types.f90
  ${GP3D_MPI_SOURCE}
  src/grid/gp3d_grid.f90
  ${GP3D_LOCAL_FFT_SOURCE}
  ${GP3D_FFT_SOURCE}
  src/solver/gp3d_solver.f90
  src/init/gp3d_initial_conditions.f90
  src/io/gp3d_restart.f90
  src/io/gp3d_io.f90
  src/io/gp3d_input.f90
  ${GP3D_GPU_FORTRAN_SOURCE}
)

if(USE_MPI)
  if(GPU_BACKEND STREQUAL "cufftmp")
    find_package(MPI REQUIRED COMPONENTS Fortran CXX)
  else()
    find_package(MPI COMPONENTS Fortran)
  endif()
endif()

add_library(gp3d STATIC ${GP3D_LIBRARY_SOURCES})
target_include_directories(gp3d PUBLIC "${CMAKE_Fortran_MODULE_DIRECTORY}")

if(FFT_BACKEND STREQUAL "fftw")
  set(FFTW_ROOT "" CACHE PATH "Optional FFTW installation prefix")
  find_library(FFTW3_LIBRARY
    NAMES fftw3 libfftw3-3
    HINTS "${FFTW_ROOT}"
    PATH_SUFFIXES lib lib64
  )
  if(NOT FFTW3_LIBRARY)
    message(FATAL_ERROR
      "FFTW backend selected, but FFTW3 library was not found. Set FFTW_ROOT or FFTW3_LIBRARY.")
  endif()
  target_link_libraries(gp3d PUBLIC "${FFTW3_LIBRARY}")
endif()

if(GPU_BACKEND STREQUAL "cuda")
  # A DLL keeps the MSVC runtime needed by nvcc on the CUDA side of the C ABI.
  # GNU Fortran links only against its small import library.
  add_library(gp3d_cuda_bridge SHARED src/gpu/gp3d_cuda_bridge.cu)
  set_target_properties(gp3d_cuda_bridge PROPERTIES
    CUDA_STANDARD 17
    CUDA_STANDARD_REQUIRED ON
    # CUDA::cudart below supplies the absolute import-library path.  Disabling
    # CMake's implicit runtime avoids bare MSVC library names in gfortran's link.
    CUDA_RUNTIME_LIBRARY None
  )
  if(GP3D_CUDA_ARCHITECTURES)
    set_property(TARGET gp3d_cuda_bridge PROPERTY
      CUDA_ARCHITECTURES "${GP3D_CUDA_ARCHITECTURES}")
  endif()
  if(WIN32)
    target_compile_options(gp3d_cuda_bridge PRIVATE
      $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/utf-8>)
  endif()
  target_link_libraries(gp3d_cuda_bridge PRIVATE CUDA::cudart CUDA::cufft)
  target_link_libraries(gp3d PUBLIC gp3d_cuda_bridge)
elseif(GPU_BACKEND STREQUAL "cufftmp")
  find_path(CUFFTMP_INCLUDE_DIR
    NAMES cufftMp.h
    HINTS "${CUFFTMP_ROOT}" "$ENV{CUFFTMP_HOME}" "$ENV{NVHPC_ROOT}"
    PATH_SUFFIXES
      include
      include/cufftmp
      math_libs/include
      "math_libs/${CUDAToolkit_VERSION_MAJOR}.${CUDAToolkit_VERSION_MINOR}/include"
      "math_libs/${CUDAToolkit_VERSION_MAJOR}.${CUDAToolkit_VERSION_MINOR}/include/cufftmp"
  )
  find_library(CUFFTMP_LIBRARY
    NAMES cufftMp
    HINTS "${CUFFTMP_ROOT}" "$ENV{CUFFTMP_HOME}" "$ENV{NVHPC_ROOT}"
    PATH_SUFFIXES
      lib lib64
      math_libs/lib64
      "math_libs/${CUDAToolkit_VERSION_MAJOR}.${CUDAToolkit_VERSION_MINOR}/lib64"
  )
  find_library(NVSHMEM_HOST_LIBRARY
    NAMES nvshmem_host
    HINTS "${NVSHMEM_ROOT}" "$ENV{NVSHMEM_HOME}" "$ENV{NVHPC_ROOT}"
    PATH_SUFFIXES
      lib lib64
      comm_libs/nvshmem/lib
      comm_libs/nvshmem/lib64
      "comm_libs/${CUDAToolkit_VERSION_MAJOR}.${CUDAToolkit_VERSION_MINOR}/nvshmem/lib"
      "comm_libs/${CUDAToolkit_VERSION_MAJOR}.${CUDAToolkit_VERSION_MINOR}/nvshmem/lib64"
  )
  if(NOT CUFFTMP_INCLUDE_DIR OR NOT CUFFTMP_LIBRARY OR NOT NVSHMEM_HOST_LIBRARY)
    message(FATAL_ERROR
      "cuFFTMp backend requires cufftMp.h, libcufftMp.so, and "
      "libnvshmem_host.so. Set CUFFTMP_ROOT/CUFFTMP_LIBRARY and "
      "NVSHMEM_ROOT/NVSHMEM_HOST_LIBRARY for the loaded NVIDIA HPC SDK.")
  endif()

  set(GP3D_CUFFTMP_LEGACY_API OFF)
  if(CUFFTMP_API STREQUAL "legacy")
    set(GP3D_CUFFTMP_LEGACY_API ON)
  elseif(CUFFTMP_API STREQUAL "auto")
    file(STRINGS "${CUFFTMP_INCLUDE_DIR}/cufftMp.h"
      GP3D_CUFFTMP_MODERN_DECLARATION REGEX "cufftMpMakePlan3d")
    if(NOT GP3D_CUFFTMP_MODERN_DECLARATION)
      set(GP3D_CUFFTMP_LEGACY_API ON)
    endif()
  endif()

  add_library(gp3d_cufftmp_bridge SHARED src/gpu/gp3d_cufftmp_bridge.cu)
  set_target_properties(gp3d_cufftmp_bridge PROPERTIES
    CUDA_STANDARD 17
    CUDA_STANDARD_REQUIRED ON
    CUDA_RUNTIME_LIBRARY Shared
  )
  if(GP3D_CUDA_ARCHITECTURES)
    set_property(TARGET gp3d_cufftmp_bridge PROPERTY
      CUDA_ARCHITECTURES "${GP3D_CUDA_ARCHITECTURES}")
  endif()
  # cuFFTMp's cufft.h must precede the CUDA Toolkit's regular cuFFT header.
  target_include_directories(gp3d_cufftmp_bridge BEFORE PRIVATE
    "${CUFFTMP_INCLUDE_DIR}")
  if(GP3D_CUFFTMP_LEGACY_API)
    target_compile_definitions(gp3d_cufftmp_bridge PRIVATE
      GP3D_CUFFTMP_LEGACY_API=1)
    message(STATUS "cuFFTMp plan API: legacy cufftMpAttachComm")
  else()
    message(STATUS "cuFFTMp plan API: modern cufftMpMakePlan3d")
  endif()
  target_link_libraries(gp3d_cufftmp_bridge PRIVATE
    "${CUFFTMP_LIBRARY}"
    "${NVSHMEM_HOST_LIBRARY}"
    CUDA::cudart
    CUDA::cuda_driver
    MPI::MPI_CXX
  )
  get_filename_component(GP3D_CUFFTMP_LIBRARY_DIR "${CUFFTMP_LIBRARY}" DIRECTORY)
  get_filename_component(GP3D_NVSHMEM_LIBRARY_DIR "${NVSHMEM_HOST_LIBRARY}" DIRECTORY)
  set_target_properties(gp3d_cufftmp_bridge PROPERTIES
    BUILD_RPATH "${GP3D_CUFFTMP_LIBRARY_DIR};${GP3D_NVSHMEM_LIBRARY_DIR}"
    INSTALL_RPATH "${GP3D_CUFFTMP_LIBRARY_DIR};${GP3D_NVSHMEM_LIBRARY_DIR}"
  )
  target_link_libraries(gp3d PUBLIC gp3d_cufftmp_bridge)
endif()

if(USE_MPI)
  if(MPI_Fortran_FOUND)
    target_link_libraries(gp3d PUBLIC MPI::MPI_Fortran)
  elseif(WIN32)
    set(MPI_ROOT "C:/Program Files (x86)/Microsoft SDKs/MPI" CACHE PATH
      "Microsoft MPI installation prefix")
    target_include_directories(gp3d PUBLIC
      "${MPI_ROOT}/Include"
      "${MPI_ROOT}/Include/x64"
    )
    target_link_directories(gp3d PUBLIC "${MPI_ROOT}/Lib/x64")
    target_link_libraries(gp3d PUBLIC msmpifec msmpi)
  else()
    message(FATAL_ERROR "USE_MPI=ON, but MPI Fortran was not found.")
  endif()
endif()

add_executable(${GP3D_EXECUTABLE_NAME} ${GP3D_MAIN_SOURCE})
target_link_libraries(${GP3D_EXECUTABLE_NAME} PRIVATE gp3d)

include(CTest)
if(BUILD_TESTING)
  add_executable(gp3d_test_restart_slf tests/test_restart_slf.f90)
  target_link_libraries(gp3d_test_restart_slf PRIVATE gp3d)
  add_test(NAME restart_slf_roundtrip
    COMMAND $<TARGET_FILE:gp3d_test_restart_slf>)
  set_tests_properties(restart_slf_roundtrip PROPERTIES
    PASS_REGULAR_EXPRESSION "SLF restart roundtrip passed"
    FIXTURES_SETUP restart_slf_data)

  if(GPU_BACKEND STREQUAL "cuda")
    add_executable(gp3d_test_cuda_splitstep tests/test_cuda_splitstep.f90)
    target_link_libraries(gp3d_test_cuda_splitstep PRIVATE gp3d)
    add_test(NAME cuda_splitstep_compare
      COMMAND $<TARGET_FILE:gp3d_test_cuda_splitstep>)
    set_tests_properties(cuda_splitstep_compare PROPERTIES
      PASS_REGULAR_EXPRESSION "CUDA split-step comparison passed")
    add_test(NAME cuda_smoke
      COMMAND $<TARGET_FILE:${GP3D_EXECUTABLE_NAME}>
        "${CMAKE_SOURCE_DIR}/tests/cuda_smoke.nml")
    set_tests_properties(cuda_smoke PROPERTIES
      PASS_REGULAR_EXPRESSION "split-step timing steps=")
    add_test(NAME cuda_restart_main
      COMMAND $<TARGET_FILE:${GP3D_EXECUTABLE_NAME}>
        "${CMAKE_SOURCE_DIR}/tests/restart_main_smoke.nml")
    set_tests_properties(cuda_restart_main PROPERTIES
      FIXTURES_REQUIRED restart_slf_data
      PASS_REGULAR_EXPRESSION "restart step=7"
      TIMEOUT 20)
  elseif(GPU_BACKEND STREQUAL "cufftmp")
    add_executable(gp3d_test_cufftmp_splitstep
      tests/test_cufftmp_splitstep.f90)
    target_link_libraries(gp3d_test_cufftmp_splitstep PRIVATE gp3d)
    if(MPIEXEC_EXECUTABLE)
      add_test(NAME cufftmp_splitstep_np2
        COMMAND "${MPIEXEC_EXECUTABLE}" ${MPIEXEC_NUMPROC_FLAG} 2
          ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:gp3d_test_cufftmp_splitstep>
          ${MPIEXEC_POSTFLAGS})
      set_tests_properties(cufftmp_splitstep_np2 PROPERTIES
        PASS_REGULAR_EXPRESSION "cuFFTMp split-step comparison passed"
        PROCESSORS 2
        RUN_SERIAL TRUE
        TIMEOUT 120)

      add_test(NAME cufftmp_argle_smoke_np2
        COMMAND "${MPIEXEC_EXECUTABLE}" ${MPIEXEC_NUMPROC_FLAG} 2
          ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:${GP3D_EXECUTABLE_NAME}>
          "${CMAKE_SOURCE_DIR}/tests/cuda_smoke.nml"
          ${MPIEXEC_POSTFLAGS})
      set_tests_properties(cufftmp_argle_smoke_np2 PROPERTIES
        PASS_REGULAR_EXPRESSION "split-step timing steps="
        PROCESSORS 2
        RUN_SERIAL TRUE
        TIMEOUT 120)
    else()
      message(WARNING "MPI launcher was not found; cuFFTMp tests are disabled.")
    endif()
  elseif(USE_MPI)
    add_executable(gp3d_test_distributed_fft tests/test_distributed_fft.f90)
    target_link_libraries(gp3d_test_distributed_fft PRIVATE gp3d)
    find_program(GP3D_MPIEXEC_EXECUTABLE NAMES mpiexec mpirun)
    if(GP3D_MPIEXEC_EXECUTABLE)
      add_test(NAME distributed_fft_np2
        COMMAND "${GP3D_MPIEXEC_EXECUTABLE}" -n 2 $<TARGET_FILE:gp3d_test_distributed_fft>)
      add_test(NAME distributed_fft_np3
        COMMAND "${GP3D_MPIEXEC_EXECUTABLE}" -n 3 $<TARGET_FILE:gp3d_test_distributed_fft>)
      add_test(NAME taylor_green_argle_np2
        COMMAND "${GP3D_MPIEXEC_EXECUTABLE}" -n 2 $<TARGET_FILE:${GP3D_EXECUTABLE_NAME}>
          "${CMAKE_SOURCE_DIR}/tests/taylor_green_smoke.nml")
      set_tests_properties(taylor_green_argle_np2 PROPERTIES
        PASS_REGULAR_EXPRESSION "split-step timing steps=")
      add_test(NAME mpi_restart_main_np2
        COMMAND "${GP3D_MPIEXEC_EXECUTABLE}" -n 2 $<TARGET_FILE:${GP3D_EXECUTABLE_NAME}>
          "${CMAKE_SOURCE_DIR}/tests/restart_main_smoke.nml")
      set_tests_properties(mpi_restart_main_np2 PROPERTIES
        FIXTURES_REQUIRED restart_slf_data
        PASS_REGULAR_EXPRESSION "restart step=7"
        TIMEOUT 20)
    endif()
  else()
    add_test(NAME taylor_green_argle_serial
      COMMAND $<TARGET_FILE:${GP3D_EXECUTABLE_NAME}>
        "${CMAKE_SOURCE_DIR}/tests/taylor_green_smoke.nml")
    set_tests_properties(taylor_green_argle_serial PROPERTIES
      PASS_REGULAR_EXPRESSION "split-step timing steps=")
    add_test(NAME serial_restart_main
      COMMAND $<TARGET_FILE:${GP3D_EXECUTABLE_NAME}>
        "${CMAKE_SOURCE_DIR}/tests/restart_main_smoke.nml")
    set_tests_properties(serial_restart_main PROPERTIES
      FIXTURES_REQUIRED restart_slf_data
      PASS_REGULAR_EXPRESSION "restart step=7"
      TIMEOUT 20)
  endif()
endif()
