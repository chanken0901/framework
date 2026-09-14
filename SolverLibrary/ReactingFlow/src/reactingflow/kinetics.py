"""Independent neutral ideal-gas rate evaluation in mol, m, s, K.

Cantera objects are consumed only by compile_kinetics, never by evaluate.
No chemical integration, density clipping or additional CFD energy source here.
"""
from dataclasses import dataclass
from bisect import bisect_right
import math
from .mechanism import MechanismError, positive
from .thermo import R, IdealGas

NEG_INF = float('-inf')


def finite(value, label):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise MechanismError(f'{label}: finite number required')
    return float(value)


def exponential(log_value):
    if log_value == NEG_INF:
        return 0.0
    try:
        value = math.exp(log_value)
    except OverflowError as exc:
        raise MechanismError('Reaction rate overflow; state is outside numerical range') from exc
    if not math.isfinite(value):
        raise MechanismError('Nonfinite reaction rate')
    return value


@dataclass(frozen=True)
class Arrhenius:
    log_a: float
    b: float
    ea_over_r: float
    sign: int = 1

    def log_abs(self, t):
        return self.log_a + self.b*math.log(t) - self.ea_over_r/t


def arrhenius(rate, order, signed=False):
    a=finite(rate.pre_exponential_factor, 'Arrhenius A')
    if a < 0 and not signed:
        raise MechanismError('Negative A is supported only inside positive-sum PLOG groups')
    # Cantera: kmol/m3. Internal: mol/m3. k_mol = k_kmol * 1000**(1-order).
    return Arrhenius(math.log(abs(a))+(1-order)*math.log(1000) if a else NEG_INF,
                     finite(rate.temperature_exponent,'Arrhenius b'),
                     finite(rate.activation_energy,'Arrhenius Ea')/(1000*R), -1 if a<0 else 1)


def chebyshev_sequence(x, n):
    result=[1.0]
    if n>1:
        result.append(x)
    for _ in range(2,n):
        result.append(2*x*result[-1]-result[-2])
    return result


@dataclass(frozen=True)
class Rate:
    kind: str
    high: Arrhenius | None = None
    low: Arrhenius | None = None
    efficiencies: tuple[float, ...] = ()
    parameters: tuple[float, ...] = ()
    plog: tuple[tuple[float, tuple[Arrhenius, ...]], ...] = ()
    coefficients: tuple[tuple[float, ...], ...] = ()
    bounds: tuple[float, ...] = ()

    def log_coefficient(self,t,p,c):
        if self.kind == 'Arrhenius':
            return self.high.log_abs(t)
        if self.kind == 'PLOG':
            def group(i):
                terms=[a.log_abs(t) for a in self.plog[i][1]]
                peak=max(terms)
                if peak==NEG_INF:
                    raise MechanismError('PLOG group has zero total rate')
                total=math.fsum(a.sign*math.exp(v-peak) for a,v in zip(self.plog[i][1],terms))
                if total<=0:
                    raise MechanismError('PLOG group sum must be positive at this temperature')
                return peak+math.log(total)
            pressures=[row[0] for row in self.plog]
            if p<=pressures[0]: return group(0)
            if p>=pressures[-1]: return group(len(pressures)-1)
            i=bisect_right(pressures,p)-1
            f=math.log(p/pressures[i])/math.log(pressures[i+1]/pressures[i])
            return (1-f)*group(i)+f*group(i+1)
        if self.kind == 'Chebyshev':
            t0,t1,p0,p1=self.bounds
            # EOS roundoff can place an exact boundary a few ulps outside.
            for endpoint in (t0,t1):
                if math.isclose(t,endpoint,rel_tol=2e-14): t=endpoint
            for endpoint in (p0,p1):
                if math.isclose(p,endpoint,rel_tol=2e-14): p=endpoint
            if not t0<=t<=t1 or not p0<=p<=p1:
                raise MechanismError('Chebyshev rate outside fitted T/P domain; no extrapolation')
            tx=(2/t-1/t0-1/t1)/(1/t1-1/t0)
            px=(2*math.log(p)-math.log(p0)-math.log(p1))/math.log(p1/p0)
            ts=chebyshev_sequence(tx,len(self.coefficients))
            ps=chebyshev_sequence(px,len(self.coefficients[0]))
            return math.log(10)*math.fsum(a*ts[i]*ps[j] for i,row in enumerate(self.coefficients) for j,a in enumerate(row))
        collider=math.fsum(e*v for e,v in zip(self.efficiencies,c))
        if collider==0:
            return NEG_INF
        log_high=self.high.log_abs(t)
        if self.kind == 'three-body':
            return log_high+math.log(collider)
        log_low=self.low.log_abs(t)
        if log_low==NEG_INF or log_high==NEG_INF:
            return NEG_INF
        log_pr=log_low+math.log(collider)-log_high
        # Stable log(Pr/(1+Pr)) for both pressure limits.
        blending=-math.log1p(math.exp(-log_pr)) if log_pr>=0 else log_pr-math.log1p(math.exp(log_pr))
        correction=0.0
        if self.kind=='Troe':
            a,t3,t1,*optional=self.parameters
            fc=(1-a)*math.exp(-t/t3)+a*math.exp(-t/t1)
            if optional: fc+=math.exp(-optional[0]/t)
            if fc<=0: raise MechanismError('Nonpositive Troe Fcent')
            fc_log=math.log10(fc)
            x=log_pr/math.log(10)-.4-.67*fc_log
            denominator=.75-1.27*fc_log-.14*x
            # Equivalent to 1/(1+(x/denominator)^2), including zero denominator.
            weight=denominator**2/(denominator**2+x*x)
            correction=math.log(10)*fc_log*weight
        elif self.kind=='SRI':
            a,b,c_sri,d,e=self.parameters
            base=a*math.exp(-b/t)+math.exp(-t/c_sri)
            if base<=0: raise MechanismError('Nonpositive SRI base')
            correction=math.log(d)+math.log(base)/(1+(log_pr/math.log(10))**2)+e*math.log(t)
        return log_high+blending+correction


@dataclass(frozen=True)
class Reaction:
    rate: Rate
    reactants: tuple[float, ...]
    products: tuple[float, ...]
    orders: tuple[float, ...]
    reversible: bool


@dataclass(frozen=True)
class Rates:
    forward: tuple[float, ...]  # mol/(m3 s), includes collider concentration
    reverse: tuple[float, ...]
    net: tuple[float, ...]
    molar_production: tuple[float, ...]
    mass_production: tuple[float, ...]  # kg/(m3 s)
    heat_release: float  # -sum(h_molar * omega_molar), W/m3, diagnostic only


@dataclass(frozen=True)
class Kinetics:
    gas: IdealGas
    reactions: tuple[Reaction, ...]

    def evaluate(self,temperature,density,y):
        t=positive(temperature,'temperature')
        rho=positive(density,'density')
        y=self.gas.fractions(y)
        c=tuple(rho*v/m for v,m in zip(y,self.gas.molar_masses))
        p=R*t*sum(c)
        positive(p,'pressure')
        properties=[n.molar(t) for n in self.gas.thermo]
        chemical=[-s[3]/(R*t)+math.log(n.reference_pressure/(R*t))
                  for s,n in zip(properties,self.gas.thermo)]
        def product(orders):
            if any(v==0 and o>0 for v,o in zip(c,orders)): return NEG_INF
            return math.fsum(o*math.log(v) for v,o in zip(c,orders) if o)
        forward,reverse,net=[],[],[]
        for r in self.reactions:
            log_k=r.rate.log_coefficient(t,p,c)
            lf=log_k+product(r.orders)
            lr=NEG_INF
            if r.reversible:
                log_kc=math.fsum((b-a)*mu for a,b,mu in zip(r.reactants,r.products,chemical))
                lr=log_k-log_kc+product(r.products)
            qf,qr=exponential(lf),exponential(lr)
            if lf==lr: q=0.0
            elif lf>lr: q=qf*(-math.expm1(lr-lf))
            else: q=qr*math.expm1(lf-lr)
            forward.append(qf); reverse.append(qr); net.append(q)
        omega=tuple(math.fsum((r.products[i]-r.reactants[i])*q for r,q in zip(self.reactions,net))
                    for i in range(len(y)))
        mass=tuple(v*m for v,m in zip(omega,self.gas.molar_masses))
        heat=-math.fsum(s[1]*v for s,v in zip(properties,omega))
        if not all(math.isfinite(v) for v in (*omega,*mass,heat)):
            raise MechanismError('Nonfinite production rate')
        return Rates(tuple(forward),tuple(reverse),tuple(net),omega,mass,heat)


def compile_kinetics(gas, solution):
    """Convert supported Cantera 3.2 rate objects to immutable SI data."""
    compiled=[]
    for i,r in enumerate(solution.reactions()):
        try:
            reactants=tuple(float(r.reactants.get(s,0)) for s in gas.names)
            products=tuple(float(r.products.get(s,0)) for s in gas.names)
            overrides=dict(r.orders)
            if set(overrides)-set(gas.names): raise MechanismError('Unknown order species')
            if r.reversible and overrides: raise MechanismError('Custom orders on reversible reactions are unsupported')
            orders=tuple(finite(overrides.get(s,a),'reaction order') for s,a in zip(gas.names,reactants))
            if any(o<0 for o in orders): raise MechanismError('Negative reaction orders are unsupported')
            n=sum(orders)
            rate=r.rate
            body=r.third_body
            eff=tuple(finite(body.efficiencies.get(s,body.default_efficiency),'third-body efficiency') for s in gas.names) if body else ()
            if any(v<0 for v in eff): raise MechanismError('Negative collider efficiency')
            if rate.type=='Arrhenius':
                own=Rate('three-body' if body else 'Arrhenius',high=arrhenius(rate,n+bool(body)),efficiencies=eff)
            elif rate.type=='falloff' and not rate.chemically_activated:
                if not body or rate.sub_type not in ('Lindemann','Troe','SRI'):
                    raise MechanismError('Unsupported falloff subtype')
                params=tuple(finite(float(v),'falloff parameter') for v in rate.falloff_coeffs)
                if rate.sub_type=='Troe':
                    if len(params) not in (3,4) or not 0<=params[0]<=1 or min(params[1:3])<=0 or (len(params)==4 and params[3]<0):
                        raise MechanismError('Unsupported Troe parameters (positive T1/T3 required)')
                if rate.sub_type=='SRI':
                    if len(params)==3: params+= (1.,0.)
                    if len(params)!=5 or params[0]<0 or params[2]<=0 or params[3]<=0:
                        raise MechanismError('Invalid SRI parameters')
                own=Rate(rate.sub_type,high=arrhenius(rate.high_rate,n),low=arrhenius(rate.low_rate,n+1),
                         efficiencies=eff,parameters=params)
            elif rate.type=='pressure-dependent-Arrhenius' and body is None:
                groups={}
                for pressure,a in rate.rates:
                    pressure=float(positive(float(pressure),'PLOG pressure'))
                    groups.setdefault(pressure,[]).append(arrhenius(a,n,signed=True))
                if not groups: raise MechanismError('Empty PLOG')
                own=Rate('PLOG',plog=tuple((p,tuple(v)) for p,v in sorted(groups.items())))
            elif rate.type=='Chebyshev' and body is None:
                data=[[finite(float(v),'Chebyshev coefficient') for v in row] for row in rate.data]
                data[0][0]+=(1-n)*math.log10(1000)
                bounds=tuple(float(v) for v in (*rate.temperature_range,*rate.pressure_range))
                if any(not math.isfinite(v) or v<=0 for v in bounds) or bounds[1]<=bounds[0] or bounds[3]<=bounds[2]:
                    raise MechanismError('Invalid Chebyshev bounds')
                own=Rate('Chebyshev',coefficients=tuple(tuple(row) for row in data),bounds=bounds)
            else:
                raise MechanismError(f'Unsupported reaction type: {r.reaction_type}')
            compiled.append(Reaction(own,reactants,products,orders,bool(r.reversible)))
        except (ValueError,OverflowError) as exc:
            raise MechanismError(f'Reaction {i+1} ({r.equation}): {exc}') from exc
    return Kinetics(gas,tuple(compiled))
