import pytest
import molspecutils.molecule as mol

def test_rotstate_inheritance():
    state_dia = mol.DiatomState(nu=0, j=10)
    assert isinstance(state_dia, mol.RotState)

    state_top = mol.SymTopState(nu=0, j=10, k=5)
    assert isinstance(state_top, mol.RotState)


def test_rotstate_inequality():
    state_dia = mol.DiatomState(nu=0, j=10)
    state_top = mol.SymTopState(nu=0, j=10, k=5)
    assert state_dia != state_top


@pytest.fixture(scope='module')
def ch3cl_mode():
    return mol.CH3ClAlchemyMode()

@pytest.fixture
def state_pair():
    return (mol.SymTopState(nu=0, j=5, k=1), mol.SymTopState(nu=1, j=6, k=1))

def test_ch3cl_nu(ch3cl_mode, state_pair):
    assert pytest.approx(ch3cl_mode.nu(state_pair)) == 2.212456947926596e+13

def test_ch3cl_mu(ch3cl_mode, state_pair):
    assert pytest.approx(ch3cl_mode.mu(state_pair)) == 1.165889466675292e-30

def test_ch3cl_gamma(ch3cl_mode, state_pair):
    assert pytest.approx(ch3cl_mode.gamma(state_pair)) == 3627488741.7999997

def test_ch3cl_delta(ch3cl_mode, state_pair):
    assert pytest.approx(ch3cl_mode.delta(state_pair)) == 0.0

def test_ch3cl_equilibrium_pop(ch3cl_mode, state_pair):
    assert pytest.approx(ch3cl_mode.equilibrium_pop(state_pair[0], 296.0)) == 0.00034793245147648864


@pytest.fixture(scope='module')
def co_mode():
    return mol.COAlchemyMode()

@pytest.fixture
def co_state_pair():
    return (mol.DiatomState(nu=0, j=1), mol.DiatomState(nu=1, j=2))

def test_co_nu(co_mode, co_state_pair):
    assert pytest.approx(co_mode.nu(co_state_pair)) == 64481039954923.66

def test_co_mu(co_mode, co_state_pair):
    assert pytest.approx(co_mode.mu(co_state_pair)) == 2.8977216594262325e-31

def test_co_gamma(co_mode, co_state_pair):
    assert pytest.approx(co_mode.gamma(co_state_pair)) == 2269428907.06

def test_co_delta(co_mode, co_state_pair):
    assert pytest.approx(co_mode.delta(co_state_pair)) == -71950189.92

def test_co_equilibrium_pop(co_mode, co_state_pair):
    assert pytest.approx(co_mode.equilibrium_pop(co_state_pair[0], 296.0)) == 0.02745169596805075
