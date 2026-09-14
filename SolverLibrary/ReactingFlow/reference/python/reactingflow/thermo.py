"""Offline reference thermodynamics in SI mol units; runtime implementation is Fortran."""
from dataclasses import dataclass
from bisect import bisect_left
import math
from .mechanism import MechanismError, positive

R = 8.31446261815324  # J/(mol K)


@dataclass(frozen=True)
class NASA:
    model: str
    bounds: tuple[float, ...]
    coefficients: tuple[tuple[float, ...], ...]
    reference_pressure: float = 101325.0

    def __post_init__(self):
        object.__setattr__(self, 'bounds', tuple(self.bounds))
        object.__setattr__(self, 'coefficients', tuple(tuple(a) for a in self.coefficients))
        width = {'NASA7': 7, 'NASA9': 9}.get(self.model)
        if width is None or len(self.bounds) < 2 or len(self.coefficients) != len(self.bounds)-1:
            raise MechanismError('Invalid NASA model/region count')
        for t in self.bounds:
            positive(t, 'temperature bound')
        if any(b <= a for a,b in zip(self.bounds, self.bounds[1:])):
            raise MechanismError('Temperature bounds must strictly increase')
        positive(self.reference_pressure, 'reference pressure')
        for row in self.coefficients:
            if len(row) != width or any(isinstance(a,bool) or not isinstance(a,(int,float)) or not math.isfinite(a) for a in row):
                raise MechanismError('Invalid NASA coefficients')

    def molar(self, temperature):
        """cp[J/mol/K], h[J/mol], s[J/mol/K], g[J/mol] at reference pressure."""
        t = positive(temperature, 'temperature')
        if not self.bounds[0] <= t <= self.bounds[-1]:
            raise MechanismError('Temperature outside NASA validity interval; no extrapolation')
        a = self.coefficients[max(0, bisect_left(self.bounds, t)-1)]
        if self.model == 'NASA7':
            cp = a[0]+a[1]*t+a[2]*t**2+a[3]*t**3+a[4]*t**4
            h = a[0]+a[1]*t/2+a[2]*t**2/3+a[3]*t**3/4+a[4]*t**4/5+a[5]/t
            s = a[0]*math.log(t)+a[1]*t+a[2]*t**2/2+a[3]*t**3/3+a[4]*t**4/4+a[6]
        else:
            cp = a[0]/t**2+a[1]/t+a[2]+a[3]*t+a[4]*t**2+a[5]*t**3+a[6]*t**4
            h = -a[0]/t**2+a[1]*math.log(t)/t+a[2]+a[3]*t/2+a[4]*t**2/3+a[5]*t**3/4+a[6]*t**4/5+a[7]/t
            s = -a[0]/(2*t**2)-a[1]/t+a[2]*math.log(t)+a[3]*t+a[4]*t**2/2+a[5]*t**3/3+a[6]*t**4/4+a[8]
        if not all(math.isfinite(v) for v in (cp,h,s)) or cp <= 1:
            raise MechanismError('Nonfinite properties or nonpositive ideal-gas cv')
        return cp*R, h*R*t, s*R, (h-s)*R*t


@dataclass(frozen=True)
class IdealGas:
    names: tuple[str, ...]
    molar_masses: tuple[float, ...]
    thermo: tuple[NASA, ...]

    def __post_init__(self):
        for key in ('names','molar_masses','thermo'):
            object.__setattr__(self, key, tuple(getattr(self,key)))
        if not self.names or len(set(self.names)) != len(self.names) or not (len(self.names)==len(self.molar_masses)==len(self.thermo)):
            raise MechanismError('Invalid mixture species layout')
        for m in self.molar_masses:
            positive(m,'molar mass')
        if self.temperature_bounds[0] >= self.temperature_bounds[1]:
            raise MechanismError('Species have no common temperature interval')

    @property
    def temperature_bounds(self):
        return max(n.bounds[0] for n in self.thermo), min(n.bounds[-1] for n in self.thermo)

    def fractions(self, y):
        y = tuple(y)
        if len(y) != len(self.names) or any(isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) or v<0 for v in y):
            raise MechanismError('Invalid mass fractions')
        if not math.isclose(sum(y),1,rel_tol=0,abs_tol=1.e-12):
            raise MechanismError('Mass fractions must sum to one; no silent normalization')
        return y

    def properties(self, temperature, y, pressure=101325.0):
        y = self.fractions(y)
        t = positive(temperature,'temperature')
        p = positive(pressure,'pressure')
        species = [n.molar(t) for n in self.thermo]
        amounts = [v/m for v,m in zip(y,self.molar_masses)]
        total = sum(amounts)
        gas_constant = R*total
        cp = sum(v*s[0] for v,s in zip(amounts,species))
        h = sum(v*s[1] for v,s in zip(amounts,species))
        entropy = sum(v*(s[2]-R*math.log((v/total)*p/n.reference_pressure))
                      for v,s,n in zip(amounts,species,self.thermo) if v>0)
        cv = cp-gas_constant
        return dict(cp=cp, cv=cv, h=h, e=h-gas_constant*t, s=entropy, g=h-t*entropy,
                    gas_constant=gas_constant, gamma=cp/cv, density=p/(gas_constant*t),
                    sound_speed=math.sqrt(cp/cv*gas_constant*t))

    def temperature_from_energy(self, energy, y):
        """Bracketed fixed-composition inversion, e includes formation energy."""
        y = self.fractions(y)
        if isinstance(energy,bool) or not isinstance(energy,(int,float)) or not math.isfinite(energy):
            raise MechanismError('Invalid specific internal energy')
        lo,hi = self.temperature_bounds
        low = self.properties(lo,y)['e']
        high = self.properties(hi,y)['e']
        if not low <= energy <= high:
            raise MechanismError('Energy outside common temperature range')
        for _ in range(100):
            mid = (lo+hi)/2
            residual = self.properties(mid,y)['e']-energy
            if abs(residual) <= 1.e-11*max(1,abs(energy)):
                return mid
            if residual>0:
                hi=mid
            else:
                lo=mid
        raise MechanismError('Temperature inversion did not converge; check polynomial continuity')
