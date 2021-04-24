from sqlalchemy.orm import registry, relationship
from sqlalchemy import Column, Integer, String, ForeignKey, Float
from molspecutils.alchemy.meta import RovibMixin, TransitionMixin, BaseMixin, LineMixin
from molspecutils.happier import CO_llq_to_pair

mapper_registry = registry()
Base = mapper_registry.generate_base(cls=BaseMixin)

def local_state_convert(llq: str, luq: str):
    j, jp = CO_llq_to_pair(llq)

    return dict(j=j), dict(j=jp)


def global_state_convert(glq: str, guq: str):
    return dict(nu=int(glq.strip())), dict(nu=int(guq.strip())) 


class RotState(Base):
    """Corresponds to local upper/lower state in HITRAN."""
    __tablename__ = 'rot_state'
    j = Column(Integer, unique=True)

    def __repr__(self):
        return f"RotState(j={self.j!r})"
   
class VibState(Base):
    """Corresponds to global upper/lower state in HITRAN."""
    __tablename__ = 'vib_state'
    nu = Column(Integer, unique=True)

    def __repr__(self):
        return f"VibState(nu={self.nu!r})"

# All the structure is in mix-ins. These classes are defined here to inherit
# from CO's Base.
class RovibState(RovibMixin, Base):
    """"""

class TransitionPair(TransitionMixin, Base):
    """"""

class LineParameters(LineMixin, Base):
    """"""
