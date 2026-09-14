"""Independent reacting-flow foundation. No NSE imports or runtime dependency."""
from .mechanism import Mechanism, MechanismError, validate_topology

__all__ = ['Mechanism', 'MechanismError', 'validate_topology']
