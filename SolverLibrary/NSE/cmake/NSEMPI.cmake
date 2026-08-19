include_guard(GLOBAL)

function(nse_configure_mpi target_name)
  string(TOUPPER "${NSE_MPI_PROVIDER}" _provider)
  if(_provider STREQUAL "AUTO")
    if(WIN32)
      set(_provider "MSMPI")
    else()
      set(_provider "SYSTEM")
    endif()
  endif()

  if(NOT _provider MATCHES "^(MSMPI|SYSTEM)$")
    message(FATAL_ERROR "NSE_MPI_PROVIDER must be AUTO, MSMPI, or SYSTEM")
  endif()

  add_library(${target_name} INTERFACE)

  if(_provider STREQUAL "MSMPI")
    set(_msmpi_default "")
    if(DEFINED ENV{MSMPI_ROOT} AND NOT "$ENV{MSMPI_ROOT}" STREQUAL "")
      set(_msmpi_default "$ENV{MSMPI_ROOT}")
    elseif(WIN32)
      set(_msmpi_default "C:/Program Files (x86)/Microsoft SDKs/MPI")
    endif()
    set(
      MSMPI_ROOT "${_msmpi_default}"
      CACHE PATH
      "Microsoft MPI SDK root containing Include and Lib/x64"
    )

    find_path(
      MSMPI_FORTRAN_INCLUDE_DIR
      NAMES mpif.h
      HINTS "${MSMPI_ROOT}/Include" ENV MSMPI_INC
    )
    find_path(
      MSMPI_X64_INCLUDE_DIR
      NAMES mpifptr.h
      HINTS "${MSMPI_ROOT}/Include/x64" ENV MSMPI_INC
      PATH_SUFFIXES x64
    )
    find_file(
      MSMPI_FORTRAN_LIBRARY
      NAMES msmpifec.lib libmsmpifec.a
      HINTS "${MSMPI_ROOT}/Lib/x64" ENV MSMPI_LIB64
    )
    find_file(
      MSMPI_LIBRARY
      NAMES msmpi.lib libmsmpi.a
      HINTS "${MSMPI_ROOT}/Lib/x64" ENV MSMPI_LIB64
    )

    set(_missing "")
    foreach(
      _variable
      MSMPI_FORTRAN_INCLUDE_DIR
      MSMPI_X64_INCLUDE_DIR
      MSMPI_FORTRAN_LIBRARY
      MSMPI_LIBRARY
    )
      if(NOT ${_variable})
        list(APPEND _missing "${_variable}")
      endif()
    endforeach()
    if(_missing)
      list(JOIN _missing ", " _missing_text)
      message(
        FATAL_ERROR
        "Microsoft MPI SDK was not found (${_missing_text}). "
        "Set MSMPI_ROOT to the SDK directory, for example "
        "C:/Program Files (x86)/Microsoft SDKs/MPI."
      )
    endif()

    target_include_directories(
      ${target_name}
      INTERFACE
        "${MSMPI_FORTRAN_INCLUDE_DIR}"
        "${MSMPI_X64_INCLUDE_DIR}"
    )
    target_link_libraries(
      ${target_name}
      INTERFACE
        "${MSMPI_FORTRAN_LIBRARY}"
        "${MSMPI_LIBRARY}"
    )

    # External Fortran MPI libraries commonly export a dependency on the
    # standard CMake target even when Microsoft MPI is configured manually.
    if(NOT TARGET MPI::MPI_Fortran)
      add_library(MPI::MPI_Fortran INTERFACE IMPORTED GLOBAL)
      set_property(
        TARGET MPI::MPI_Fortran
        PROPERTY INTERFACE_INCLUDE_DIRECTORIES
          "${MSMPI_FORTRAN_INCLUDE_DIR};${MSMPI_X64_INCLUDE_DIR}"
      )
      set_property(
        TARGET MPI::MPI_Fortran
        PROPERTY INTERFACE_LINK_LIBRARIES
          "${MSMPI_FORTRAN_LIBRARY};${MSMPI_LIBRARY}"
      )
    endif()

    message(STATUS "Microsoft MPI include: ${MSMPI_FORTRAN_INCLUDE_DIR}")
    message(STATUS "Microsoft MPI x64 include: ${MSMPI_X64_INCLUDE_DIR}")
    message(STATUS "Microsoft MPI Fortran library: ${MSMPI_FORTRAN_LIBRARY}")
  else()
    find_package(MPI REQUIRED COMPONENTS Fortran)
    target_link_libraries(${target_name} INTERFACE MPI::MPI_Fortran)
  endif()

  set(NSE_MPI_PROVIDER_RESOLVED "${_provider}" CACHE INTERNAL "Resolved NSE MPI provider" FORCE)
endfunction()
