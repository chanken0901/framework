#pragma once

// Included inside the CUDA bridge's anonymous namespace after the KEEP and
// WENO face-flux functions. New leaf schemes only need one branch in
// leaf_face_flux_cuda; sensor implementations remain independent.

__device__ inline void hybrid_velocity_pressure_at_cuda(
    const double* q,
    const GridView& grid,
    int i,
    int j,
    int k,
    double gamma,
    double small_rho,
    double small_p,
    double velocity[3],
    double& pressure) {
  double density;
  primitive_state(
      q, grid, cell_index(grid, i, j, k), gamma, small_rho, small_p,
      density, velocity[0], velocity[1], velocity[2], pressure);
}

__device__ inline double hybrid_pressure_at_cuda(
    const double* q,
    const GridView& grid,
    int i,
    int j,
    int k,
    double gamma,
    double small_rho,
    double small_p) {
  double velocity[3];
  double pressure;
  hybrid_velocity_pressure_at_cuda(
      q, grid, i, j, k, gamma, small_rho, small_p, velocity, pressure);
  return pressure;
}

__device__ inline double ducros_pressure_cell_sensor_cuda(
    const double* q,
    const GridView& grid,
    int i,
    int j,
    int k,
    int direction,
    double gamma,
    double small_rho,
    double small_p,
    double inverse_dx,
    double inverse_dy,
    double inverse_dz) {
  constexpr double sensor_epsilon = 1.0e-30;
  double velocity_xm[3], velocity_xp[3];
  double velocity_ym[3], velocity_yp[3];
  double velocity_zm[3], velocity_zp[3];
  double unused_pressure;
  hybrid_velocity_pressure_at_cuda(
      q, grid, i - 1, j, k, gamma, small_rho, small_p,
      velocity_xm, unused_pressure);
  hybrid_velocity_pressure_at_cuda(
      q, grid, i + 1, j, k, gamma, small_rho, small_p,
      velocity_xp, unused_pressure);
  hybrid_velocity_pressure_at_cuda(
      q, grid, i, j - 1, k, gamma, small_rho, small_p,
      velocity_ym, unused_pressure);
  hybrid_velocity_pressure_at_cuda(
      q, grid, i, j + 1, k, gamma, small_rho, small_p,
      velocity_yp, unused_pressure);
  hybrid_velocity_pressure_at_cuda(
      q, grid, i, j, k - 1, gamma, small_rho, small_p,
      velocity_zm, unused_pressure);
  hybrid_velocity_pressure_at_cuda(
      q, grid, i, j, k + 1, gamma, small_rho, small_p,
      velocity_zp, unused_pressure);

  const double divergence =
      0.5 * inverse_dx * (velocity_xp[0] - velocity_xm[0])
      + 0.5 * inverse_dy * (velocity_yp[1] - velocity_ym[1])
      + 0.5 * inverse_dz * (velocity_zp[2] - velocity_zm[2]);
  const double vorticity_x =
      0.5 * inverse_dy * (velocity_yp[2] - velocity_ym[2])
      - 0.5 * inverse_dz * (velocity_zp[1] - velocity_zm[1]);
  const double vorticity_y =
      0.5 * inverse_dz * (velocity_zp[0] - velocity_zm[0])
      - 0.5 * inverse_dx * (velocity_xp[2] - velocity_xm[2]);
  const double vorticity_z =
      0.5 * inverse_dx * (velocity_xp[1] - velocity_xm[1])
      - 0.5 * inverse_dy * (velocity_yp[0] - velocity_ym[0]);
  const double rate_scale = fmax(
      fabs(divergence),
      fmax(fabs(vorticity_x), fmax(fabs(vorticity_y), fabs(vorticity_z))));
  double ducros_factor = 0.0;
  if (rate_scale > 1.0e-15) {
    const double scaled_divergence = divergence / rate_scale;
    const double scaled_vorticity_x = vorticity_x / rate_scale;
    const double scaled_vorticity_y = vorticity_y / rate_scale;
    const double scaled_vorticity_z = vorticity_z / rate_scale;
    const double vorticity_squared =
        scaled_vorticity_x * scaled_vorticity_x
        + scaled_vorticity_y * scaled_vorticity_y
        + scaled_vorticity_z * scaled_vorticity_z;
    const double compression = fmin(scaled_divergence, 0.0);
    ducros_factor = compression * compression / (
        scaled_divergence * scaled_divergence
        + vorticity_squared + sensor_epsilon);
  }

  const int di = direction == 0 ? 1 : 0;
  const int dj = direction == 1 ? 1 : 0;
  const int dk = direction == 2 ? 1 : 0;
  const double pressure_minus = hybrid_pressure_at_cuda(
      q, grid, i - di, j - dj, k - dk,
      gamma, small_rho, small_p);
  const double pressure_center = hybrid_pressure_at_cuda(
      q, grid, i, j, k, gamma, small_rho, small_p);
  const double pressure_plus = hybrid_pressure_at_cuda(
      q, grid, i + di, j + dj, k + dk,
      gamma, small_rho, small_p);
  const double pressure_curvature = fabs(
      pressure_plus - 2.0 * pressure_center + pressure_minus) / (
      pressure_plus + 2.0 * pressure_center + pressure_minus
      + sensor_epsilon);
  return fmax(0.0, fmin(1.0, pressure_curvature * ducros_factor));
}

__device__ inline double hybrid_face_weight_cuda(
    const double* q,
    const GridView& grid,
    int i,
    int j,
    int k,
    int direction,
    int sensor,
    double sensor_onset,
    double sensor_full,
    double gamma,
    double small_rho,
    double small_p,
    double inverse_dx,
    double inverse_dy,
    double inverse_dz) {
  if (sensor != sensor_ducros_pressure) {
    return 0.0;
  }
  const int di = direction == 0 ? 1 : 0;
  const int dj = direction == 1 ? 1 : 0;
  const int dk = direction == 2 ? 1 : 0;
  const double left_sensor = ducros_pressure_cell_sensor_cuda(
      q, grid, i, j, k, direction, gamma, small_rho, small_p,
      inverse_dx, inverse_dy, inverse_dz);
  const double right_sensor = ducros_pressure_cell_sensor_cuda(
      q, grid, i + di, j + dj, k + dk, direction,
      gamma, small_rho, small_p, inverse_dx, inverse_dy, inverse_dz);
  double scaled = (fmax(left_sensor, right_sensor) - sensor_onset)
      / (sensor_full - sensor_onset);
  scaled = fmax(0.0, fmin(1.0, scaled));
  return scaled * scaled * (3.0 - 2.0 * scaled);
}

__device__ inline void leaf_face_flux_cuda(
    int scheme,
    const double* q,
    const GridView& grid,
    int i,
    int j,
    int k,
    int direction,
    double gamma,
    double small_rho,
    double small_p,
    double flux[5]) {
  if (scheme == convective_weno5z_roe) {
    weno5z_roe_face_flux_cuda(
        q, grid, i, j, k, direction, gamma, small_rho, small_p, flux);
  } else {
    const int keep_order = scheme == convective_keep2 ? 2 : 6;
    keep_face_flux_cuda(
        q, grid, i, j, k, direction, keep_order,
        gamma, small_rho, small_p, flux);
  }
}

__device__ inline void hybrid_face_flux_cuda(
    const double* q,
    const GridView& grid,
    int i,
    int j,
    int k,
    int direction,
    int smooth_scheme,
    int shock_scheme,
    int sensor,
    double sensor_onset,
    double sensor_full,
    double gamma,
    double small_rho,
    double small_p,
    double inverse_dx,
    double inverse_dy,
    double inverse_dz,
    double flux[5]) {
  const double alpha = hybrid_face_weight_cuda(
      q, grid, i, j, k, direction, sensor, sensor_onset, sensor_full,
      gamma, small_rho, small_p, inverse_dx, inverse_dy, inverse_dz);
  if (alpha <= 0.0) {
    leaf_face_flux_cuda(
        smooth_scheme, q, grid, i, j, k, direction,
        gamma, small_rho, small_p, flux);
    return;
  }
  if (alpha >= 1.0) {
    leaf_face_flux_cuda(
        shock_scheme, q, grid, i, j, k, direction,
        gamma, small_rho, small_p, flux);
    return;
  }
  double smooth_flux[5], shock_flux[5];
  leaf_face_flux_cuda(
      smooth_scheme, q, grid, i, j, k, direction,
      gamma, small_rho, small_p, smooth_flux);
  leaf_face_flux_cuda(
      shock_scheme, q, grid, i, j, k, direction,
      gamma, small_rho, small_p, shock_flux);
  for (int variable = 0; variable < 5; ++variable) {
    flux[variable] =
        (1.0 - alpha) * smooth_flux[variable] + alpha * shock_flux[variable];
  }
}

__global__ void rhs_hybrid_kernel(
    const double* q,
    double* rhs,
    GridView grid,
    int smooth_scheme,
    int shock_scheme,
    int sensor,
    double sensor_onset,
    double sensor_full,
    double gamma,
    double small_rho,
    double small_p,
    double inverse_dx,
    double inverse_dy,
    double inverse_dz,
    std::size_t physical_count) {
  const std::size_t linear =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (linear >= physical_count) {
    return;
  }
  const int i = static_cast<int>(linear % grid.nx) + grid.nghost;
  const std::size_t jk = linear / grid.nx;
  const int j = static_cast<int>(jk % grid.ny) + grid.nghost;
  const int k = static_cast<int>(jk / grid.ny) + grid.nghost;
  const std::size_t center = cell_index(grid, i, j, k);
  const double inverse_spacing[3] = {inverse_dx, inverse_dy, inverse_dz};
  double minus_flux[5], plus_flux[5];
  double result[5] = {0.0, 0.0, 0.0, 0.0, 0.0};

  for (int direction = 0; direction < 3; ++direction) {
    const int di = direction == 0 ? 1 : 0;
    const int dj = direction == 1 ? 1 : 0;
    const int dk = direction == 2 ? 1 : 0;
    hybrid_face_flux_cuda(
        q, grid, i - di, j - dj, k - dk, direction,
        smooth_scheme, shock_scheme, sensor, sensor_onset, sensor_full,
        gamma, small_rho, small_p, inverse_dx, inverse_dy, inverse_dz,
        minus_flux);
    hybrid_face_flux_cuda(
        q, grid, i, j, k, direction,
        smooth_scheme, shock_scheme, sensor, sensor_onset, sensor_full,
        gamma, small_rho, small_p, inverse_dx, inverse_dy, inverse_dz,
        plus_flux);
    for (int variable = 0; variable < 5; ++variable) {
      result[variable] -= inverse_spacing[direction]
          * (plus_flux[variable] - minus_flux[variable]);
    }
  }
  for (int variable = 0; variable < 5; ++variable) {
    rhs[center + grid.cell_count * static_cast<std::size_t>(variable)] =
        result[variable];
  }
}
