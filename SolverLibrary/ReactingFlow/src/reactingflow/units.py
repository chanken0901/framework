"""Single boundary between dimensionless CFD variables and SI thermo inputs."""
from dataclasses import dataclass
from .mechanism import positive, MechanismError
import math


@dataclass(frozen=True)
class ReferenceScales:
    density: float  # kg/m3
    velocity: float  # m/s (not implicitly sound speed)
    length: float  # m
    temperature: float  # K

    def __post_init__(self):
        for key in ('density','velocity','length','temperature'):
            positive(getattr(self,key),key)

    def scale(self, quantity):
        scales=dict(density=self.density,velocity=self.velocity,length=self.length,
                    temperature=self.temperature,time=self.length/self.velocity,
                    pressure=self.density*self.velocity**2,energy=self.velocity**2,
                    energy_density=self.density*self.velocity**2)
        if quantity not in scales:
            raise MechanismError(f'Unknown reference quantity: {quantity}')
        positive(scales[quantity],quantity+' scale')
        return scales[quantity]

    def to_si(self, value, quantity):
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):
            raise MechanismError('Invalid dimensionless value')
        return value*self.scale(quantity)

    def from_si(self, value, quantity):
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):
            raise MechanismError('Invalid SI value')
        return value/self.scale(quantity)
