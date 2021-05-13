import scipy.constants as C

def wn2nu(wn):
    return wn*C.c*1e2


def nu2wn(nu):
    return nu/C.c*1e-2


def wn2joule(wn):
    return C.h*wn2nu(wn)


def joule2wn(E):
    return nu2wn(E/C.h)
