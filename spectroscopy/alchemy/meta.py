"""Take data from HAPI and put it in an sqlite database with sqlalchemy."""
from pathlib import Path
from sqlalchemy.orm import registry, relationship, declarative_mixin, declared_attr
from sqlalchemy import Column, Integer, String, ForeignKey, Float, UniqueConstraint
import appdirs

dirs = appdirs.AppDirs('happier', 'gkowzan')
hitran_cache = str(Path(dirs.user_cache_dir) / 'db')

class BaseMixin:
    id = Column(Integer, primary_key=True)

@declarative_mixin
class RovibMixin:
    __tablename__ = 'rovib_state'

    @declared_attr
    def nu_id(cls):
        return Column(Integer, ForeignKey('vib_state.id'))

    @declared_attr
    def j_id(cls):
        return Column(Integer, ForeignKey('rot_state.id'))
    
    energy = Column(Float)
    g = Column(Float)

    @declared_attr
    def j(cls):
        return relationship("RotState")

    @declared_attr
    def nu(cls):
        return relationship("VibState")

    __table_args__ = (UniqueConstraint('nu_id', 'j_id'),)

    def __repr__(self):
        return f"RovibState(nu={self.nu!r}, j={self.j!r}, energy={self.energy!r}, g={self.g!r})"

@declarative_mixin
class TransitionMixin:
    __tablename__ = 'transition_pair'

    @declared_attr
    def statepp_id(cls):
        return Column(Integer, ForeignKey('rovib_state.id'))

    @declared_attr
    def statep_id(cls):
        return Column(Integer, ForeignKey('rovib_state.id'))

    @declared_attr
    def statepp(cls):
        return relationship("RovibState", foreign_keys="TransitionPair.statepp_id")

    @declared_attr
    def statep(cls):
        return relationship("RovibState", foreign_keys="TransitionPair.statep_id")

    __table_args__ = (UniqueConstraint('statepp_id', 'statep_id'),)

    def __repr__(self):
        return f"TransitionPair(statepp={self.statepp!r}, statep={self.statep!r})"

@declarative_mixin
class LineMixin:
    __tablename__ = 'line_parameters'

    nu = Column(Float)
    sw = Column(Float)
    A = Column(Float)
    gamma_air = Column(Float)
    gamma_self = Column(Float)
    n_air = Column(Float)
    delta_air = Column(Float)

    @declared_attr
    def transition_id(cls):
        return Column(Integer, ForeignKey('transition_pair.id'), unique=True)

    @declared_attr
    def transition(cls):
        return relationship("TransitionPair")

    def __repr__(self):
        return (f"LineParameters(nu={self.nu!r}, sw={self.sw!r}, A={self.A!r}, gamma_air={self.gamma_air!r}, "
                f"gamma_self={self.gamma_self!r}, n_air={self.n_air!r}, delta_air={self.delta_air!r}, "
                f"transition={self.transition!r})")
    
if __name__ == '__main__':
    from sqlalchemy import create_engine, select
    from sqlalchemy.orm import Session

    engine = create_engine("sqlite+pysqlite:///:memory:", echo=False, future=True)
    Base.metadata.create_all(engine)
    with Session(bind=engine) as session:
        for j in range(20):
            session.add(CORotState(j=j))
        for nu in range(3):
            session.add(COVibState(nu=nu))

        # add radiative states
        statep = RovibState(nu=session.get(COVibState, 1),
                            j=session.get(CORotState, 2))
        statepp = RovibState(nu=session.get(COVibState, 0),
                             j=session.get(CORotState, 1))
        session.add(statep)
        session.add(statepp)

        # add transition pair
        transition = TransitionPair(statepp=statepp, statep=statep)
        session.add(transition)
        print(session.execute(select(TransitionPair)).one())
