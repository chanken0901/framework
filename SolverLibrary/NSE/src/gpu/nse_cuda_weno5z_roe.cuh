#pragma once

// This header is included inside the CUDA bridge's anonymous namespace after
// GridView and the state indexing helpers have been declared.

__device__ inline void weno5z_smoothness(
    const double value[5], double beta[3]) {
  const double d20 = value[0] - 2.0 * value[1] + value[2];
  const double d21 = value[1] - 2.0 * value[2] + value[3];
  const double d22 = value[2] - 2.0 * value[3] + value[4];
  const double d10 = value[0] - 4.0 * value[1] + 3.0 * value[2];
  const double d11 = value[1] - value[3];
  const double d12 = 3.0 * value[2] - 4.0 * value[3] + value[4];
  beta[0] = (13.0 / 12.0) * d20 * d20 + 0.25 * d10 * d10;
  beta[1] = (13.0 / 12.0) * d21 * d21 + 0.25 * d11 * d11;
  beta[2] = (13.0 / 12.0) * d22 * d22 + 0.25 * d12 * d12;
}

__device__ inline double reconstruct_weno5z_left_cuda(
    const double value[5]) {
  constexpr double epsilon = 1.0e-20;
  const double candidate[3] = {
      (2.0 * value[0] - 7.0 * value[1] + 11.0 * value[2]) / 6.0,
      (-value[1] + 5.0 * value[2] + 2.0 * value[3]) / 6.0,
      (2.0 * value[2] + 5.0 * value[3] - value[4]) / 6.0};
  double beta[3];
  weno5z_smoothness(value, beta);
  const double tau5 = fabs(beta[0] - beta[2]);
  const double ratio0 = tau5 / (beta[0] + epsilon);
  const double ratio1 = tau5 / (beta[1] + epsilon);
  const double ratio2 = tau5 / (beta[2] + epsilon);
  const double alpha[3] = {
      0.1 * (1.0 + ratio0 * ratio0),
      0.6 * (1.0 + ratio1 * ratio1),
      0.3 * (1.0 + ratio2 * ratio2)};
  const double weight_sum = alpha[0] + alpha[1] + alpha[2];
  return (alpha[0] * candidate[0] + alpha[1] * candidate[1]
          + alpha[2] * candidate[2]) / weight_sum;
}

__device__ inline double reconstruct_weno5z_right_cuda(
    const double value[5]) {
  constexpr double epsilon = 1.0e-20;
  const double candidate[3] = {
      (-value[0] + 5.0 * value[1] + 2.0 * value[2]) / 6.0,
      (2.0 * value[1] + 5.0 * value[2] - value[3]) / 6.0,
      (11.0 * value[2] - 7.0 * value[3] + 2.0 * value[4]) / 6.0};
  double beta[3];
  weno5z_smoothness(value, beta);
  const double tau5 = fabs(beta[0] - beta[2]);
  const double ratio0 = tau5 / (beta[0] + epsilon);
  const double ratio1 = tau5 / (beta[1] + epsilon);
  const double ratio2 = tau5 / (beta[2] + epsilon);
  const double alpha[3] = {
      0.3 * (1.0 + ratio0 * ratio0),
      0.6 * (1.0 + ratio1 * ratio1),
      0.1 * (1.0 + ratio2 * ratio2)};
  const double weight_sum = alpha[0] + alpha[1] + alpha[2];
  return (alpha[0] * candidate[0] + alpha[1] * candidate[1]
          + alpha[2] * candidate[2]) / weight_sum;
}

__device__ inline void load_normal_state_cuda(
    const double* q,
    const GridView& grid,
    int i,
    int j,
    int k,
    int direction,
    int offset,
    double state[5]) {
  const int x = i + (direction == 0 ? offset : 0);
  const int y = j + (direction == 1 ? offset : 0);
  const int z = k + (direction == 2 ? offset : 0);
  const std::size_t cell = cell_index(grid, x, y, z);
  state[0] = q[cell];
  state[4] = q[cell + 4 * grid.cell_count];
  if (direction == 0) {
    state[1] = q[cell + grid.cell_count];
    state[2] = q[cell + 2 * grid.cell_count];
    state[3] = q[cell + 3 * grid.cell_count];
  } else if (direction == 1) {
    state[1] = q[cell + 2 * grid.cell_count];
    state[2] = q[cell + grid.cell_count];
    state[3] = q[cell + 3 * grid.cell_count];
  } else {
    state[1] = q[cell + 3 * grid.cell_count];
    state[2] = q[cell + grid.cell_count];
    state[3] = q[cell + 2 * grid.cell_count];
  }
}

__device__ inline void roe_primitive_state_cuda(
    const double state[5],
    double gamma,
    double small_rho,
    double small_p,
    double& rho,
    double& u,
    double& v,
    double& w,
    double& pressure,
    double& enthalpy) {
  rho = fmax(state[0], small_rho);
  u = state[1] / rho;
  v = state[2] / rho;
  w = state[3] / rho;
  pressure = fmax(
      (gamma - 1.0)
          * (state[4] - 0.5 * rho * (u * u + v * v + w * w)),
      small_p);
  enthalpy = (state[4] + pressure) / rho;
}

__device__ inline void roe_eigensystem_cuda(
    const double left_state[5],
    const double right_state[5],
    double gamma,
    double small_rho,
    double small_p,
    double right_matrix[5][5],
    double left_matrix[5][5],
    double eigenvalue[5]) {
  double rho_l, u_l, v_l, w_l, pressure_l, enthalpy_l;
  double rho_r, u_r, v_r, w_r, pressure_r, enthalpy_r;
  roe_primitive_state_cuda(
      left_state, gamma, small_rho, small_p,
      rho_l, u_l, v_l, w_l, pressure_l, enthalpy_l);
  roe_primitive_state_cuda(
      right_state, gamma, small_rho, small_p,
      rho_r, u_r, v_r, w_r, pressure_r, enthalpy_r);

  const double sqrt_l = sqrt(rho_l);
  const double sqrt_r = sqrt(rho_r);
  const double denominator = sqrt_l + sqrt_r;
  const double roe_density = fmax(sqrt_l * sqrt_r, small_rho);
  const double u = (sqrt_l * u_l + sqrt_r * u_r) / denominator;
  const double v = (sqrt_l * v_l + sqrt_r * v_r) / denominator;
  const double w = (sqrt_l * w_l + sqrt_r * w_r) / denominator;
  const double enthalpy =
      (sqrt_l * enthalpy_l + sqrt_r * enthalpy_r) / denominator;
  const double velocity_squared = u * u + v * v + w * w;
  const double sound_speed_squared = fmax(
      (gamma - 1.0) * (enthalpy - 0.5 * velocity_squared),
      small_p / roe_density);
  const double sound_speed = sqrt(sound_speed_squared);

  eigenvalue[0] = u - sound_speed;
  eigenvalue[1] = u;
  eigenvalue[2] = u;
  eigenvalue[3] = u;
  eigenvalue[4] = u + sound_speed;
  for (int row = 0; row < 5; ++row) {
    for (int column = 0; column < 5; ++column) {
      right_matrix[row][column] = 0.0;
      left_matrix[row][column] = 0.0;
    }
  }
  const double right_columns[5][5] = {
      {1.0, u - sound_speed, v, w, enthalpy - u * sound_speed},
      {1.0, u, v, w, 0.5 * velocity_squared},
      {0.0, 0.0, 1.0, 0.0, v},
      {0.0, 0.0, 0.0, 1.0, w},
      {1.0, u + sound_speed, v, w, enthalpy + u * sound_speed}};
  for (int column = 0; column < 5; ++column) {
    for (int row = 0; row < 5; ++row) {
      right_matrix[row][column] = right_columns[column][row];
    }
  }

  const double gamma_minus_one = gamma - 1.0;
  const double b1 = gamma_minus_one * velocity_squared
      / (2.0 * sound_speed_squared);
  const double b2 = gamma_minus_one / sound_speed_squared;
  left_matrix[0][0] = 0.5 * (b1 + u / sound_speed);
  left_matrix[0][1] = 0.5 * (-b2 * u - 1.0 / sound_speed);
  left_matrix[0][2] = -0.5 * b2 * v;
  left_matrix[0][3] = -0.5 * b2 * w;
  left_matrix[0][4] = 0.5 * b2;
  left_matrix[1][0] = 1.0 - b1;
  left_matrix[1][1] = b2 * u;
  left_matrix[1][2] = b2 * v;
  left_matrix[1][3] = b2 * w;
  left_matrix[1][4] = -b2;
  left_matrix[2][0] = -v;
  left_matrix[2][2] = 1.0;
  left_matrix[3][0] = -w;
  left_matrix[3][3] = 1.0;
  left_matrix[4][0] = 0.5 * (b1 - u / sound_speed);
  left_matrix[4][1] = 0.5 * (-b2 * u + 1.0 / sound_speed);
  left_matrix[4][2] = -0.5 * b2 * v;
  left_matrix[4][3] = -0.5 * b2 * w;
  left_matrix[4][4] = 0.5 * b2;
}

__device__ inline void euler_physical_flux_cuda(
    const double state[5],
    double gamma,
    double small_rho,
    double small_p,
    double flux[5]) {
  double rho, u, v, w, pressure, enthalpy;
  roe_primitive_state_cuda(
      state, gamma, small_rho, small_p,
      rho, u, v, w, pressure, enthalpy);
  flux[0] = rho * u;
  flux[1] = rho * u * u + pressure;
  flux[2] = rho * u * v;
  flux[3] = rho * u * w;
  flux[4] = u * (state[4] + pressure);
}

__device__ inline double harten_hyman_fix_cuda(
    double eigenvalue, double delta) {
  if (delta > 0.0 && fabs(eigenvalue) < delta) {
    return 0.5 * (eigenvalue * eigenvalue / delta + delta);
  }
  return fabs(eigenvalue);
}

__device__ inline void roe_numerical_flux_cuda(
    const double left_state[5],
    const double right_state[5],
    double gamma,
    double small_rho,
    double small_p,
    double flux[5]) {
  double left_flux[5], right_flux[5];
  double right_matrix[5][5], left_matrix[5][5], eigenvalue[5];
  double rho_l, u_l, v_l, w_l, pressure_l, enthalpy_l;
  double rho_r, u_r, v_r, w_r, pressure_r, enthalpy_r;
  euler_physical_flux_cuda(
      left_state, gamma, small_rho, small_p, left_flux);
  euler_physical_flux_cuda(
      right_state, gamma, small_rho, small_p, right_flux);
  roe_eigensystem_cuda(
      left_state, right_state, gamma, small_rho, small_p,
      right_matrix, left_matrix, eigenvalue);
  roe_primitive_state_cuda(
      left_state, gamma, small_rho, small_p,
      rho_l, u_l, v_l, w_l, pressure_l, enthalpy_l);
  roe_primitive_state_cuda(
      right_state, gamma, small_rho, small_p,
      rho_r, u_r, v_r, w_r, pressure_r, enthalpy_r);
  const double sound_l = sqrt(gamma * pressure_l / rho_l);
  const double sound_r = sqrt(gamma * pressure_r / rho_r);
  const double entropy_delta[5] = {
      fmax(0.0, (u_r - sound_r) - (u_l - sound_l)),
      fmax(0.0, u_r - u_l),
      fmax(0.0, u_r - u_l),
      fmax(0.0, u_r - u_l),
      fmax(0.0, (u_r + sound_r) - (u_l + sound_l))};
  double wave_strength[5];
  for (int wave = 0; wave < 5; ++wave) {
    wave_strength[wave] = 0.0;
    for (int variable = 0; variable < 5; ++variable) {
      wave_strength[wave] += left_matrix[wave][variable]
          * (right_state[variable] - left_state[variable]);
    }
    wave_strength[wave] *=
        harten_hyman_fix_cuda(eigenvalue[wave], entropy_delta[wave]);
  }
  for (int variable = 0; variable < 5; ++variable) {
    double dissipation = 0.0;
    for (int wave = 0; wave < 5; ++wave) {
      dissipation += right_matrix[variable][wave] * wave_strength[wave];
    }
    flux[variable] =
        0.5 * (left_flux[variable] + right_flux[variable] - dissipation);
  }
}

__device__ inline void rotate_flux_to_global_cuda(
    const double normal_flux[5], int direction, double global_flux[5]) {
  global_flux[0] = normal_flux[0];
  global_flux[4] = normal_flux[4];
  if (direction == 0) {
    global_flux[1] = normal_flux[1];
    global_flux[2] = normal_flux[2];
    global_flux[3] = normal_flux[3];
  } else if (direction == 1) {
    global_flux[1] = normal_flux[2];
    global_flux[2] = normal_flux[1];
    global_flux[3] = normal_flux[3];
  } else {
    global_flux[1] = normal_flux[2];
    global_flux[2] = normal_flux[3];
    global_flux[3] = normal_flux[1];
  }
}

// Test raw conservative variables, not the pressure floored for flux evaluation.
__device__ inline bool admissible_roe_state_cuda(
    const double s[5], double gamma, double small_rho, double small_p) {
  for (int v = 0; v < 5; ++v) if (!isfinite(s[v])) return false;
  if (s[0] < small_rho) return false;
  const double p = (gamma-1.0)*(s[4]-0.5*(s[1]*s[1]+s[2]*s[2]+s[3]*s[3])/s[0]);
  return isfinite(p) && p >= small_p;
}

__device__ inline bool limit_roe_state_cuda(const double center[5],
    double state[5], double gamma, double small_rho, double small_p) {
  if (admissible_roe_state_cuda(state, gamma, small_rho, small_p)) return false;
  double original[5];
  for (int v = 0; v < 5; ++v) {
    if (!isfinite(state[v])) {
      for (int c = 0; c < 5; ++c) state[c] = center[c];
      return true;
    }
    original[v] = state[v];
  }
  double low = 0.0, high = 1.0;
  for (int iteration = 0; iteration < 50; ++iteration) {
    const double theta = 0.5*(low+high);
    for (int v = 0; v < 5; ++v) state[v] = center[v]+theta*(original[v]-center[v]);
    if (admissible_roe_state_cuda(state, gamma, small_rho, small_p)) low = theta;
    else high = theta;
  }
  for (int v = 0; v < 5; ++v) state[v] = center[v]+(0.99*low)*(original[v]-center[v]);
  return true;
}

__device__ inline void weno5z_roe_face_flux_cuda(
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
  double center_left[5], center_right[5];
  double right_matrix[5][5], left_matrix[5][5], eigenvalue[5];
  double stencil_state[5], characteristic_stencil[5];
  double left_characteristic[5], right_characteristic[5];
  double left_state[5], right_state[5], normal_flux[5];
  load_normal_state_cuda(
      q, grid, i, j, k, direction, 0, center_left);
  load_normal_state_cuda(
      q, grid, i, j, k, direction, 1, center_right);
  roe_eigensystem_cuda(
      center_left, center_right, gamma, small_rho, small_p,
      right_matrix, left_matrix, eigenvalue);

  for (int characteristic = 0; characteristic < 5; ++characteristic) {
    for (int point = 0; point < 5; ++point) {
      load_normal_state_cuda(
          q, grid, i, j, k, direction, point - 2, stencil_state);
      characteristic_stencil[point] = 0.0;
      for (int variable = 0; variable < 5; ++variable) {
        characteristic_stencil[point] +=
            left_matrix[characteristic][variable] * stencil_state[variable];
      }
    }
    left_characteristic[characteristic] =
        reconstruct_weno5z_left_cuda(characteristic_stencil);
    for (int point = 0; point < 5; ++point) {
      load_normal_state_cuda(
          q, grid, i, j, k, direction, point - 1, stencil_state);
      characteristic_stencil[point] = 0.0;
      for (int variable = 0; variable < 5; ++variable) {
        characteristic_stencil[point] +=
            left_matrix[characteristic][variable] * stencil_state[variable];
      }
    }
    right_characteristic[characteristic] =
        reconstruct_weno5z_right_cuda(characteristic_stencil);
  }

  for (int variable = 0; variable < 5; ++variable) {
    left_state[variable] = 0.0;
    right_state[variable] = 0.0;
    for (int characteristic = 0; characteristic < 5; ++characteristic) {
      left_state[variable] += right_matrix[variable][characteristic]
          * left_characteristic[characteristic];
      right_state[variable] += right_matrix[variable][characteristic]
          * right_characteristic[characteristic];
    }
  }
  const bool limited_left = limit_roe_state_cuda(center_left, left_state, gamma, small_rho, small_p);
  const bool limited_right = limit_roe_state_cuda(center_right, right_state, gamma, small_rho, small_p);
  double rl, ul, vl, wl, pl, hl, rr, ur, vr, wr, pr, hr;
  roe_primitive_state_cuda(left_state, gamma, small_rho, small_p, rl, ul, vl, wl, pl, hl);
  roe_primitive_state_cuda(right_state, gamma, small_rho, small_p, rr, ur, vr, wr, pr, hr);
  const double cl = sqrt(gamma*pl/rl), cr = sqrt(gamma*pr/rr);
  if (limited_left || limited_right || ur-ul > 2.0*fmin(cl, cr)) {
    // Shared local Lax-Friedrichs face flux retains conservation.
    double fl[5], fr[5];
    euler_physical_flux_cuda(left_state, gamma, small_rho, small_p, fl);
    euler_physical_flux_cuda(right_state, gamma, small_rho, small_p, fr);
    const double speed = fmax(fabs(ul)+cl, fabs(ur)+cr);
    for (int v = 0; v < 5; ++v)
      normal_flux[v] = 0.5*(fl[v]+fr[v]-speed*(right_state[v]-left_state[v]));
  } else {
    roe_numerical_flux_cuda(left_state, right_state, gamma, small_rho, small_p, normal_flux);
  }
  rotate_flux_to_global_cuda(normal_flux, direction, flux);
}

__global__ void rhs_weno5z_roe_kernel(
    const double* q,
    double* rhs,
    GridView grid,
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
    weno5z_roe_face_flux_cuda(
        q, grid, i - di, j - dj, k - dk, direction,
        gamma, small_rho, small_p, minus_flux);
    weno5z_roe_face_flux_cuda(
        q, grid, i, j, k, direction,
        gamma, small_rho, small_p, plus_flux);
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
