"""Offline SciPy reference reactor. Production reactor uses Fortran DVODE."""
from dataclasses import dataclass
import math

from .mechanism import MechanismError, positive
from .thermo import R


@dataclass(frozen=True)
class ReactorResult:
    mode: str
    species: tuple
    time: tuple
    temperature: tuple
    pressure: tuple
    density: tuple
    mass_fractions: tuple
    ignition_delay: float | None
    diagnostics: dict


def integrate(mechanism, *, temperature, pressure, mass_fractions, end_time,
              mode='constant_volume', rtol=1.e-7, atol_species=1.e-14,
              atol_temperature=1.e-6, max_step=None, max_steps=100000,
              ignition_temperature_rise=400.0):
    """Return every accepted BDF step (no output-state clipping/normalization).

    State is [T, Y_1, ..., Y_N]. Numerical-Jacobian/Newton trial compositions
    use a nonnegative, normalized extension of the RHS. Accepted states must
    satisfy the original strict composition contract or integration fails.
    Ignition delay is the first upward T0+rise crossing, not peak dT/dt.
    """
    import numpy as np
    from scipy.integrate import BDF
    from scipy.optimize import brentq

    if mode not in ('constant_volume', 'constant_pressure'):
        raise MechanismError('mode must be constant_volume or constant_pressure')
    for label, value in [('temperature', temperature), ('pressure', pressure),
                         ('end_time', end_time), ('rtol', rtol),
                         ('atol_species', atol_species), ('atol_temperature', atol_temperature),
                         ('ignition_temperature_rise', ignition_temperature_rise)]:
        positive(value, label)
    if rtol < 1.e-12 or rtol > 1.e-2:
        raise MechanismError('rtol must be between 1e-12 and 1e-2')
    if isinstance(max_steps, bool) or not isinstance(max_steps, int) or max_steps < 1:
        raise MechanismError('max_steps must be a positive integer')
    step_limit = end_time if max_step is None else positive(max_step, 'max_step')
    gas = mechanism.gas
    y0 = gas.fractions(mass_fractions)
    initial = gas.properties(temperature, y0, pressure)
    rho0 = initial['density']
    cv_mode = mode == 'constant_volume'
    energy_key = 'e' if cv_mode else 'h'
    energy0 = initial[energy_key]
    energy_scale = max(1., abs(energy0), initial['cp']*temperature)
    extended_calls = 0

    def rhs(time, state):
        nonlocal extended_calls
        t = float(state[0])
        raw = state[1:]
        if not np.all(np.isfinite(state)):
            raise MechanismError('Nonfinite BDF trial state')
        # This extension is for unconstrained Newton/Jacobian trials only.
        # Never write these values back into the solver's solution vector.
        if np.any(raw < 0):
            extended_calls += 1
        y = np.maximum(raw, 0.)
        total = float(y.sum())
        if total <= 0:
            raise MechanismError('BDF trial has no positive species')
        y = (y/total).tolist()
        props = gas.properties(t, y, pressure)
        rho = rho0 if cv_mode else props['density']
        rates = mechanism.kinetics.evaluate(t, rho, y)
        dy = np.asarray(rates.mass_production)/rho
        energies = [(n.molar(t)[1] - (R*t if cv_mode else 0.))/m
                    for n, m in zip(gas.thermo, gas.molar_masses)]
        dt = -math.fsum(e*v for e, v in zip(energies, dy))/props['cv' if cv_mode else 'cp']
        return np.r_[dt, dy]

    def elements(y):
        return tuple(math.fsum(v/m*atoms[e] for v, m, atoms in
                     zip(y, gas.molar_masses, mechanism.topology.compositions))
                     for e in range(len(mechanism.topology.elements)))

    elem0 = elements(y0)
    times, temperatures, pressures, densities, fractions = [], [], [], [], []
    errors = dict(mass_sum_error=0., element_relative_error=0., energy_relative_error=0.)

    def record(time, state):
        t = float(state[0])
        y = tuple(float(v) for v in state[1:])
        try:
            props = gas.properties(t, y, pressure)
        except MechanismError as exc:
            raise MechanismError(f'Invalid accepted state at t={time:.17g}: {exc}; '
                                 'reduce tolerances/max_step; state was not clipped') from exc
        rho = rho0 if cv_mode else props['density']
        p = rho*props['gas_constant']*t if cv_mode else pressure
        errors['mass_sum_error'] = max(errors['mass_sum_error'], abs(math.fsum(y)-1))
        errors['element_relative_error'] = max(errors['element_relative_error'],
            max((abs(a-b)/max(1., abs(b)) for a,b in zip(elements(y), elem0)), default=0.))
        errors['energy_relative_error'] = max(errors['energy_relative_error'],
                                             abs(props[energy_key]-energy0)/energy_scale)
        if errors['element_relative_error'] > 1.e-8 or errors['energy_relative_error'] > max(100*rtol, 1.e-7):
            raise MechanismError(f'Conservation check failed at t={time:.17g}: {errors}')
        times.append(float(time)); temperatures.append(t); pressures.append(float(p))
        densities.append(float(rho)); fractions.append(y)

    state0 = np.asarray([temperature, *y0], dtype=float)
    record(0., state0)
    solver = BDF(rhs, 0., state0, end_time, rtol=rtol,
                 atol=np.asarray([atol_temperature]+[atol_species]*len(y0)),
                 max_step=step_limit, jac=None)
    ignition = None
    threshold = temperature+ignition_temperature_rise
    while solver.status == 'running':
        if len(times)-1 >= max_steps:
            raise MechanismError('BDF exceeded max_steps; no successful result returned')
        previous_time, previous_temperature = solver.t, float(solver.y[0])
        message = solver.step()
        if solver.status == 'failed':
            raise MechanismError(f'BDF failed at t={solver.t}: {message}')
        record(solver.t, solver.y)
        if ignition is None and previous_temperature < threshold <= solver.y[0]:
            dense = solver.dense_output()
            ignition = float(brentq(lambda t: float(dense(t)[0])-threshold,
                previous_time, solver.t, xtol=max(1.e-30, end_time*1.e-13)))
    errors.update(steps=len(times)-1, nfev=solver.nfev, njev=solver.njev, nlu=solver.nlu,
                  negative_trial_rhs_calls=extended_calls, energy_quantity=energy_key,
                  energy_scale=energy_scale, ignition_temperature=threshold)
    return ReactorResult(mode, gas.names, tuple(times), tuple(temperatures), tuple(pressures),
                         tuple(densities), tuple(fractions), ignition, errors)
