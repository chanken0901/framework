if(NOT DEFINED CONTRACT_FILE)
  message(FATAL_ERROR "CONTRACT_FILE is required")
endif()
file(READ "${CONTRACT_FILE}" _contract)
foreach(_required IN ITEMS
    "backend: cufftmp"
    "status: implemented"
    "host_gather_allowed: false"
    "state_location: distributed_device_resident"
    "forward_distributed_fft"
    "project_local_spectral_pencil"
    "allreduce_forcing_denominators"
    "inverse_distributed_fft"
    "add_local_momentum_rhs"
    "runtime_available: true")
  string(FIND "${_contract}" "${_required}" _position)
  if(_position EQUAL -1)
    message(FATAL_ERROR "cuFFTMp contract is missing: ${_required}")
  endif()
endforeach()
message(STATUS "cuFFTMp forcing contract test passed")
