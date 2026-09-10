execute_process(COMMAND "${CMAKE_COMMAND}" -E env "NSE_CUDA_MPI_TRANSPORT=${MODE}"
  "${LAUNCHER}" "${NP_FLAG}" 2 "${SOLVER}" "${INPUT}"
  RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE errors TIMEOUT 40)
if(MODE STREQUAL "device")
  set(expected "Device MPI requested but unavailable")
else()
  set(expected "NSE_CUDA_MPI_TRANSPORT must be")
endif()
if(result STREQUAL "0" OR NOT "${output}${errors}" MATCHES "${expected}")
  message(FATAL_ERROR "Expected transport rejection, got ${result}: ${output}${errors}")
endif()
message(STATUS "Unsupported/invalid CUDA MPI transport rejected safely")
