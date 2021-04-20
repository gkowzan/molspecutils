from sqlalchemy.orm import registry, relationship
from sqlalchemy import Column, Integer, String, ForeignKey, Float
from spectroscopy.alchemy.meta import RovibMixin, TransitionMixin, BaseMixin, LineMixin
from spectroscopy.happier import CH3Cl_gq_to_dict, CH3Cl_lq_to_dict

mapper_registry = registry()
Base = mapper_registry.generate_base(cls=BaseMixin)

def local_state_convert(llq: str, luq: str):
    return CH3Cl_lq_to_dict(llq), CH3Cl_lq_to_dict(luq)


def global_state_convert(glq: str, guq: str):
    return CH3Cl_gq_to_dict(glq), CH3Cl_gq_to_dict(guq) 


class RotState(Base):
    """Corresponds to local upper/lower state in HITRAN."""
    __tablename__ = 'rot_state'
    
    j = Column(Integer)
    k = Column(Integer)
    l = Column(Integer)
    f = Column(Float)
    sym = Column(String(2))

    def __repr__(self):
        return f"RotState(j={self.j!r}, k={self.k!r}, l={self.l!r}, f={self.f!r}, sym={self.sym!r})"

class VibState(Base):
    """Corresponds to global upper/lower state in HITRAN."""
    __tablename__ = 'vib_state'

    nu1 = Column(Integer)
    nu2 = Column(Integer)
    nu3 = Column(Integer)
    nu4 = Column(Integer)
    nu5 = Column(Integer)
    nu6 = Column(Integer)

    def __repr__(self):
        return (f"VibState(nu1={self.nu1!r}, nu2={self.nu2!r}, nu3={self.nu3!r}, "
                f"nu4={self.nu4!r}, nu5={self.nu5!r}, nu6={self.nu6!r})")

# All the structure is in mix-ins. These classes are defined here to inherit
# from CH3Cl's Base.
class RovibState(RovibMixin, Base):
    """"""

class TransitionPair(TransitionMixin, Base):
    """"""

class LineParameters(LineMixin, Base):
    """"""
