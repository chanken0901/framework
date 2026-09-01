module mod_model_config
  use mod_precision, only : dp
  implicit none
  private

  public :: gpe_config, nse_config
  public :: init_gpe_config, init_nse_config
  public :: print_gpe_config, print_nse_config
  public :: resolve_nse_flow_parameters
  public :: nse_boundary_face_count
  public :: nse_face_x_min, nse_face_x_max, nse_face_y_min
  public :: nse_face_y_max, nse_face_z_min, nse_face_z_max

  integer, parameter :: nse_boundary_face_count = 6
  integer, parameter :: nse_face_x_min = 1
  integer, parameter :: nse_face_x_max = 2
  integer, parameter :: nse_face_y_min = 3
  integer, parameter :: nse_face_y_max = 4
  integer, parameter :: nse_face_z_min = 5
  integer, parameter :: nse_face_z_max = 6
  character(len=5), parameter :: nse_boundary_face_name(6) = [ &
    character(len=5) :: 'x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max']

  type :: nse_config
    integer :: nv = 5
    real(dp) :: gamma = 1.4_dp
    real(dp) :: cfl = 0.5_dp
    real(dp) :: small_rho = 1.0e-12_dp
    real(dp) :: small_p   = 1.0e-12_dp
    real(dp) :: rho0 = 1.0_dp
    real(dp) :: mach = 0.5_dp
    real(dp) :: reynolds = 0.0_dp
    real(dp) :: prandtl = 0.72_dp
    character(len=32) :: convective_scheme = 'keep6'
    character(len=32) :: hybrid_smooth_scheme = 'keep6'
    character(len=32) :: hybrid_shock_scheme = 'weno5z_roe'
    character(len=32) :: hybrid_sensor = 'ducros_pressure'
    real(dp) :: hybrid_sensor_onset = 0.01_dp
    real(dp) :: hybrid_sensor_full = 0.10_dp
    character(len=32) :: viscous_scheme = 'none'
    character(len=32) :: boundary_condition = 'periodic'
    character(len=32) :: boundary_face_type(nse_boundary_face_count) = 'periodic'
    real(dp) :: boundary_reference_rho(nse_boundary_face_count) = -1.0_dp
    real(dp) :: boundary_reference_velocity(3,nse_boundary_face_count) = 0.0_dp
    real(dp) :: boundary_reference_p(nse_boundary_face_count) = -1.0_dp
    real(dp) :: boundary_relaxation_strength = 0.1_dp
    real(dp) :: boundary_length_scale = -1.0_dp
    character(len=32) :: time_integrator = 'ssprk3'
    character(len=32) :: hit_spectrum = 'johnsen'
    integer :: hit_seed = 13579
    real(dp) :: hit_turbulent_mach = -1.0_dp
    real(dp) :: hit_turbulent_reynolds = -1.0_dp
    real(dp) :: hit_rms_velocity = -1.0_dp
    real(dp) :: hit_peak_wavenumber = 4.0_dp
    real(dp) :: hit_integral_length = 1.0_dp
    real(dp) :: hit_kolmogorov_length = 0.02_dp
    real(dp) :: hit_integral_reynolds = -1.0_dp
    real(dp) :: hit_taylor_microscale = -1.0_dp
    real(dp) :: hit_kinematic_viscosity = -1.0_dp
    real(dp) :: hit_johnsen_length_scale_ratio = 2.0_dp
    real(dp) :: hit_pope_energy_constant = 1.5_dp
    real(dp) :: hit_pope_large_scale_constant = 6.78_dp
    real(dp) :: hit_pope_dissipation_constant = 0.40_dp
    real(dp) :: hit_pope_large_scale_exponent = 2.0_dp
    real(dp) :: hit_pope_dissipation_exponent = 5.2_dp
    real(dp) :: hit_dealias_fraction = 2.0_dp / 3.0_dp
    character(len=32) :: hit_isotropy_mode = 'projected_shell'
    real(dp) :: hit_isotropy_k_cutoff = 2.5_dp
    real(dp) :: hit_isotropy_tolerance = 1.0e-8_dp
    integer :: hit_isotropy_max_iterations = 80
    character(len=512) :: imported_turbulence_file = ''
    character(len=32) :: imported_turbulence_mode = 'embed'
    real(dp) :: imported_turbulence_x_start = -1.0e300_dp
    integer :: imported_turbulence_blend_cells = 0
    real(dp) :: imported_turbulence_velocity_offset_x = 0.0_dp
    real(dp) :: imported_turbulence_velocity_offset_y = 0.0_dp
    real(dp) :: imported_turbulence_velocity_offset_z = 0.0_dp
    real(dp) :: imported_turbulence_background_rho = -1.0_dp
    real(dp) :: imported_turbulence_background_u = 0.0_dp
    real(dp) :: imported_turbulence_background_v = 0.0_dp
    real(dp) :: imported_turbulence_background_w = 0.0_dp
    real(dp) :: imported_turbulence_background_p = -1.0_dp
    real(dp) :: planar_shock_position = -1.0e300_dp
    character(len=32) :: planar_shock_direction = 'positive_x'
    real(dp) :: planar_shock_mach = -1.0_dp
    real(dp) :: planar_shock_upstream_rho = -1.0_dp
    real(dp) :: planar_shock_upstream_u = 0.0_dp
    real(dp) :: planar_shock_upstream_v = 0.0_dp
    real(dp) :: planar_shock_upstream_w = 0.0_dp
    real(dp) :: planar_shock_upstream_p = -1.0_dp
    real(dp) :: planar_shock_downstream_rho = -1.0_dp
    real(dp) :: planar_shock_downstream_u = 0.0_dp
    real(dp) :: planar_shock_downstream_v = 0.0_dp
    real(dp) :: planar_shock_downstream_w = 0.0_dp
    real(dp) :: planar_shock_downstream_p = -1.0_dp
    real(dp) :: shock_tube_diaphragm_position = -1.0e300_dp
    real(dp) :: shock_tube_driver_rho = -1.0_dp
    real(dp) :: shock_tube_driver_u = 0.0_dp
    real(dp) :: shock_tube_driver_v = 0.0_dp
    real(dp) :: shock_tube_driver_w = 0.0_dp
    real(dp) :: shock_tube_driver_p = -1.0_dp
    real(dp) :: shock_tube_driven_rho = -1.0_dp
    real(dp) :: shock_tube_driven_u = 0.0_dp
    real(dp) :: shock_tube_driven_v = 0.0_dp
    real(dp) :: shock_tube_driven_w = 0.0_dp
    real(dp) :: shock_tube_driven_p = -1.0_dp
    character(len=32) :: forcing_scheme = 'none'
    character(len=32) :: forcing_spectrum = 'low_wavenumber'
    character(len=32) :: forcing_fft_backend = 'auto'
    real(dp) :: forcing_k_cutoff = 2.5_dp
    real(dp) :: forcing_target_dissipation = 0.0_dp
    real(dp) :: forcing_dilatational_ratio = 0.0_dp
    real(dp) :: forcing_denominator_floor = 1.0e-14_dp
    real(dp) :: forcing_max_coefficient = 0.0_dp
    integer :: forcing_report_interval = 100
  end type nse_config

  type :: gpe_config
    real(dp) :: g = 50.0_dp
    real(dp) :: sigma0 = 1.0_dp
    real(dp) :: wx = 1.0_dp
    real(dp) :: wy = 1.0_dp
    real(dp) :: wz = 1.0_dp
    real(dp) :: hbar = 1.0_dp
    real(dp) :: mass = 1.0_dp
  end type gpe_config

contains

  subroutine init_gpe_config(cfg)
    type(gpe_config), intent(out) :: cfg
    cfg = gpe_config()
  end subroutine init_gpe_config

  subroutine init_nse_config(cfg)
    type(nse_config), intent(out) :: cfg
    cfg = nse_config()
  end subroutine init_nse_config

  subroutine resolve_nse_flow_parameters(cfg, initial_condition)
    type(nse_config), intent(inout) :: cfg
    character(len=*), intent(in) :: initial_condition
    real(dp) :: component_rms, velocity_scale
    logical :: has_mach_target, has_reynolds_target

    select case (trim(adjustl(initial_condition)))
    case ('hit', 'hit_spectral', 'homogeneous_isotropic_turbulence')
      continue
    case default
      return
    end select

    has_mach_target = cfg%hit_turbulent_mach > 0.0_dp
    has_reynolds_target = cfg%hit_turbulent_reynolds > 0.0_dp
    if (.not. has_mach_target .and. .not. has_reynolds_target) return
    if (has_mach_target .neqv. has_reynolds_target) then
      error stop 'HIT requires both turbulent Mach and Reynolds targets'
    end if
    if (cfg%hit_integral_length <= 0.0_dp) then
      error stop 'HIT characteristic/integral length must be positive'
    end if
    if (cfg%prandtl <= 0.0_dp) then
      error stop 'HIT transport resolution requires prandtl > 0'
    end if

    ! The original code treats Re_ini as Re_lambda and uses
    ! Re_L = 3 Re_lambda^2 / 20. Its integral-scale velocity is
    ! sqrt(3/2) times the one-component RMS velocity.
    component_rms = cfg%hit_turbulent_mach / sqrt(3.0_dp)
    velocity_scale = sqrt(1.5_dp) * component_rms
    cfg%hit_integral_reynolds = 3.0_dp * &
      cfg%hit_turbulent_reynolds**2 / 20.0_dp
    cfg%hit_kinematic_viscosity = velocity_scale * &
      cfg%hit_integral_length / cfg%hit_integral_reynolds
    cfg%hit_taylor_microscale = cfg%hit_integral_length * &
      sqrt(10.0_dp / cfg%hit_integral_reynolds)
    cfg%hit_kolmogorov_length = cfg%hit_integral_length * &
      cfg%hit_integral_reynolds**(-0.75_dp)

    cfg%mach = cfg%hit_turbulent_mach
    cfg%hit_rms_velocity = component_rms
    cfg%reynolds = 1.0_dp / cfg%hit_kinematic_viscosity

    select case (trim(adjustl(cfg%hit_spectrum)))
    case ('johnsen', 'k4_gaussian')
      if (cfg%hit_johnsen_length_scale_ratio <= 0.0_dp) then
        error stop 'Johnsen length-scale ratio must be positive'
      end if
      cfg%hit_peak_wavenumber = 2.0_dp * &
        cfg%hit_johnsen_length_scale_ratio / cfg%hit_integral_length
    case ('pope')
      continue
    case default
      write(*,'(A,A)') 'ERROR: unsupported HIT spectrum: ', &
        trim(cfg%hit_spectrum)
      error stop
    end select
  end subroutine resolve_nse_flow_parameters

  subroutine print_gpe_config(cfg, unit)
    type(gpe_config), intent(in) :: cfg
    integer, intent(in), optional :: unit
    integer :: u
    u = 6
    if (present(unit)) u = unit
    write(u,'(A)') '--- gpe_config ---'
    write(u,'(A,ES16.8)') 'g      = ', cfg%g
    write(u,'(A,ES16.8)') 'sigma0 = ', cfg%sigma0
    write(u,'(A,3ES16.8)') 'wx, wy, wz = ', cfg%wx, cfg%wy, cfg%wz
  end subroutine print_gpe_config

  subroutine print_nse_config(cfg, unit)
    type(nse_config), intent(in) :: cfg
    integer, intent(in), optional :: unit
    integer :: u, face
    u = 6
    if (present(unit)) u = unit
    write(u,'(A)') '--- nse_config ---'
    write(u,'(A,I10)') 'nv     = ', cfg%nv
    write(u,'(A,ES16.8)') 'gamma  = ', cfg%gamma
    write(u,'(A,ES16.8)') 'cfl    = ', cfg%cfl
    write(u,'(A,ES16.8)') 'rho0   = ', cfg%rho0
    write(u,'(A,ES16.8)') 'mach   = ', cfg%mach
    write(u,'(A,ES16.8)') 'reynolds   = ', cfg%reynolds
    write(u,'(A,ES16.8)') 'prandtl   = ', cfg%prandtl
    write(u,'(A,A)') 'convective_scheme = ', trim(cfg%convective_scheme)
    write(u,'(A,A)') 'hybrid_smooth_scheme = ', &
      trim(cfg%hybrid_smooth_scheme)
    write(u,'(A,A)') 'hybrid_shock_scheme  = ', &
      trim(cfg%hybrid_shock_scheme)
    write(u,'(A,A)') 'hybrid_sensor        = ', trim(cfg%hybrid_sensor)
    write(u,'(A,ES16.8)') 'hybrid_sensor_onset = ', &
      cfg%hybrid_sensor_onset
    write(u,'(A,ES16.8)') 'hybrid_sensor_full  = ', &
      cfg%hybrid_sensor_full
    write(u,'(A,A)') 'viscous_scheme    = ', trim(cfg%viscous_scheme)
    write(u,'(A,A)') 'boundary_condition = ', trim(cfg%boundary_condition)
    do face = 1, nse_boundary_face_count
      write(u,'(A,A,A,A)') 'boundary_', &
        trim(nse_boundary_face_name(face)), ' = ', &
        trim(cfg%boundary_face_type(face))
      if (trim(adjustl(cfg%boundary_face_type(face))) == 'non_reflecting') then
        write(u,'(A,5ES16.8)') '  reference rho,u,v,w,p = ', &
          cfg%boundary_reference_rho(face), &
          cfg%boundary_reference_velocity(:,face), &
          cfg%boundary_reference_p(face)
      end if
    end do
    write(u,'(A,ES16.8)') 'boundary_relaxation_strength = ', &
      cfg%boundary_relaxation_strength
    write(u,'(A,ES16.8)') 'boundary_length_scale = ', &
      cfg%boundary_length_scale
    write(u,'(A,A)') 'time_integrator   = ', trim(cfg%time_integrator)
    write(u,'(A,A)') 'hit_spectrum      = ', trim(cfg%hit_spectrum)
    write(u,'(A,I10)') 'hit_seed          = ', cfg%hit_seed
    write(u,'(A,ES16.8)') 'hit_turbulent_mach = ', &
      cfg%hit_turbulent_mach
    write(u,'(A,ES16.8)') 'hit_turbulent_reynolds = ', &
      cfg%hit_turbulent_reynolds
    write(u,'(A,ES16.8)') 'hit_rms_velocity  = ', cfg%hit_rms_velocity
    write(u,'(A,ES16.8)') 'hit_peak_wavenumber = ', cfg%hit_peak_wavenumber
    write(u,'(A,ES16.8)') 'hit_integral_length = ', cfg%hit_integral_length
    write(u,'(A,ES16.8)') 'hit_kolmogorov_length = ', &
      cfg%hit_kolmogorov_length
    write(u,'(A,ES16.8)') 'hit_integral_reynolds = ', &
      cfg%hit_integral_reynolds
    write(u,'(A,ES16.8)') 'hit_taylor_microscale = ', &
      cfg%hit_taylor_microscale
    write(u,'(A,ES16.8)') 'hit_kinematic_viscosity = ', &
      cfg%hit_kinematic_viscosity
    write(u,'(A,ES16.8)') 'hit_johnsen_length_scale_ratio = ', &
      cfg%hit_johnsen_length_scale_ratio
    write(u,'(A,ES16.8)') 'hit_pope_energy_constant = ', &
      cfg%hit_pope_energy_constant
    write(u,'(A,ES16.8)') 'hit_pope_large_scale_constant = ', &
      cfg%hit_pope_large_scale_constant
    write(u,'(A,ES16.8)') 'hit_pope_dissipation_constant = ', &
      cfg%hit_pope_dissipation_constant
    write(u,'(A,ES16.8)') 'hit_pope_large_scale_exponent = ', &
      cfg%hit_pope_large_scale_exponent
    write(u,'(A,ES16.8)') 'hit_pope_dissipation_exponent = ', &
      cfg%hit_pope_dissipation_exponent
    write(u,'(A,ES16.8)') 'hit_dealias_fraction = ', &
      cfg%hit_dealias_fraction
    write(u,'(A,A)') 'hit_isotropy_mode = ', trim(cfg%hit_isotropy_mode)
    write(u,'(A,ES16.8)') 'hit_isotropy_k_cutoff = ', &
      cfg%hit_isotropy_k_cutoff
    write(u,'(A,ES16.8)') 'hit_isotropy_tolerance = ', &
      cfg%hit_isotropy_tolerance
    write(u,'(A,I10)') 'hit_isotropy_max_iterations = ', &
      cfg%hit_isotropy_max_iterations
    if (len_trim(cfg%imported_turbulence_file) > 0) then
      write(u,'(A,A)') 'imported_turbulence_file = ', &
        trim(cfg%imported_turbulence_file)
      write(u,'(A,A)') 'imported_turbulence_mode = ', &
        trim(cfg%imported_turbulence_mode)
      write(u,'(A,ES16.8)') 'imported_turbulence_x_start = ', &
        cfg%imported_turbulence_x_start
      write(u,'(A,I10)') 'imported_turbulence_blend_cells = ', &
        cfg%imported_turbulence_blend_cells
      write(u,'(A,3ES16.8)') 'imported_turbulence_velocity_offset = ', &
        cfg%imported_turbulence_velocity_offset_x, &
        cfg%imported_turbulence_velocity_offset_y, &
        cfg%imported_turbulence_velocity_offset_z
      write(u,'(A,5ES16.8)') 'imported_turbulence_background = ', &
        cfg%imported_turbulence_background_rho, &
        cfg%imported_turbulence_background_u, &
        cfg%imported_turbulence_background_v, &
        cfg%imported_turbulence_background_w, &
        cfg%imported_turbulence_background_p
    end if
    if (cfg%planar_shock_position > -1.0e250_dp) then
      write(u,'(A,ES16.8)') 'planar_shock_position = ', &
        cfg%planar_shock_position
      write(u,'(A,A)') 'planar_shock_direction = ', &
        trim(cfg%planar_shock_direction)
      write(u,'(A,ES16.8)') 'planar_shock_mach = ', cfg%planar_shock_mach
      write(u,'(A,5ES16.8)') 'planar_shock_upstream = ', &
        cfg%planar_shock_upstream_rho, cfg%planar_shock_upstream_u, &
        cfg%planar_shock_upstream_v, cfg%planar_shock_upstream_w, &
        cfg%planar_shock_upstream_p
      write(u,'(A,5ES16.8)') 'planar_shock_downstream = ', &
        cfg%planar_shock_downstream_rho, cfg%planar_shock_downstream_u, &
        cfg%planar_shock_downstream_v, cfg%planar_shock_downstream_w, &
        cfg%planar_shock_downstream_p
    end if
    if (cfg%shock_tube_diaphragm_position > -1.0e250_dp) then
      write(u,'(A,ES16.8)') 'shock_tube_diaphragm_position = ', &
        cfg%shock_tube_diaphragm_position
      write(u,'(A,5ES16.8)') 'shock_tube_driver = ', &
        cfg%shock_tube_driver_rho, cfg%shock_tube_driver_u, &
        cfg%shock_tube_driver_v, cfg%shock_tube_driver_w, &
        cfg%shock_tube_driver_p
      write(u,'(A,5ES16.8)') 'shock_tube_driven = ', &
        cfg%shock_tube_driven_rho, cfg%shock_tube_driven_u, &
        cfg%shock_tube_driven_v, cfg%shock_tube_driven_w, &
        cfg%shock_tube_driven_p
    end if
    write(u,'(A,A)') 'forcing_scheme    = ', trim(cfg%forcing_scheme)
    write(u,'(A,A)') 'forcing_spectrum  = ', trim(cfg%forcing_spectrum)
    write(u,'(A,A)') 'forcing_fft_backend = ', trim(cfg%forcing_fft_backend)
    write(u,'(A,ES16.8)') 'forcing_k_cutoff = ', cfg%forcing_k_cutoff
    write(u,'(A,ES16.8)') 'forcing_target_dissipation = ', &
      cfg%forcing_target_dissipation
    write(u,'(A,ES16.8)') 'forcing_dilatational_ratio = ', &
      cfg%forcing_dilatational_ratio
    write(u,'(A,ES16.8)') 'forcing_denominator_floor = ', &
      cfg%forcing_denominator_floor
    write(u,'(A,ES16.8)') 'forcing_max_coefficient = ', &
      cfg%forcing_max_coefficient
    write(u,'(A,I10)') 'forcing_report_interval = ', &
      cfg%forcing_report_interval
  end subroutine print_nse_config

end module mod_model_config
