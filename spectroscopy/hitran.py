"""Deprecated HITRAN class and helper functions."""
import os
import io
from glob import glob
from pathlib import Path
from functools import lru_cache
import numpy as np
from scipy.special import erfcx
import scipy.constants as C
from shed.units import wn2nu, nu2wn

DATA_DIR = str(Path.home() / 'documents/nauka/cs/hitran2012')
DATA_DIR = '/mnt/d/cs/hitran2012'

c = C.c
N_A = C.N_A
k = C.k
amu = C.value('atomic mass constant')
opj = os.path.join


def cef(z):
    return erfcx(-z*1.0j)


def complex_voigt(v, v0, doppler, homo):
    deltav = v-v0
    a = np.sqrt(np.log(2))/doppler*(deltav + homo*1.0j)
    voigt = cef(a)
    peak_value = c*100*np.sqrt(np.log(2)/np.pi)/doppler

    return peak_value*voigt


def simplecount(filehandle):
    lines = 0
    for line in filehandle:
        lines += 1
    return lines


class HITRAN:
    trans = {'1': 'H2O', '2': 'CO2', '3': 'O3', '4': 'N2O', '5': 'CO',
             '6': 'CH4', '7': 'O2', '8': 'NO', '9': 'SO2', '10': 'NO2',
             '11': 'NH3', '12': 'HNO3', '13': 'OH', '14': 'HF', '15':
             'HCl', '16': 'HBr', '17': 'HI', '18': 'ClO', '19': 'OCS',
             '20': 'H2CO', '21': 'HOCl', '22': 'N2', '23': 'HCN', '24':
             'CH3Cl', '25': 'H2O2', '26': 'C2H2', '27': 'C2H6', '28':
             'PH3', '29': 'COF2', '30': 'SF6', '31': 'H2S', '32':
             'HCOOH', '33': 'HO2', '34': 'O', '35': 'ClONO2', '36':
             'NO+', '37': 'HOBr', '38': 'C2H4', '39': 'CH3OH', '40':
             'CH3Br', '41': 'CH3CN', '42': 'CF4', '43': 'C4H2', '44':
             'HC3N', '45': 'H2', '46': 'CS', '47': 'SO3'}

    recip_volumes = {'CO2': 1.977e-3/44.01*N_A,
                     'CO': 1.145e-3/28.010*N_A}

    mass = {'CO2': 44.01*amu,
            'CO': 28.01*amu}

    def __init__(self, molecule, dataDir=DATA_DIR):
        self.molecule = molecule.upper()
        self.num = self.name2num(self.molecule)
        self.dataDir = dataDir

        fpathPar = glob(opj(dataDir, '%02d_hit*.par' % self.num))
        fpathZip = glob(opj(dataDir, '%02d_hit*.zip' % self.num))

        if fpathPar:
            self.fpath = fpathPar[0]
        elif fpathZip:
            self.fpath = fpathZip[0]
        else:
            raise ValueError("Data for molecule '%s' is unavailable" %
                             self.molecule)

        self.get_data = lru_cache()(self.get_data)

    @classmethod
    def num2name(cls, M):
        '''For a given input molecule identifier number, return the
        corresponding molecular formula.

        Parameters
        ----------
        M : int
            The HITRAN molecule identifier number.

        Returns
        -------
        molecular_formula : str
            The string describing the molecule.
        '''
        return cls.trans[str(M)]

    @classmethod
    def name2num(cls, molecule_name):
        """For a given input molecular formula, return the corresponding
        HITRAN molecule identifier number.

        Parameters
        ----------
        molecular_formula : str
            The string describing the molecule.

        Returns
        -------
        M : int
            The HITRAN molecular identified number.
        """
        # trans = {v:k for k,v in trans.items()}
        return(int({v: k for k, v in cls.trans.items()}[molecule_name]))

    @staticmethod
    def parread(filename):
        if filename.endswith('.zip'):
            import zipfile
            zip = zipfile.ZipFile(filename, 'r')
            insideName = os.path.split(filename[:-3]+'par')[-1]
            filehandle = io.TextIOWrapper(zip.open(insideName, 'r'))
        else:
            filehandle = open(filename, 'r')

        lines = simplecount(filehandle)
        if filename.endswith('.zip'):
            insideName = os.path.split(filename[:-3]+'par')[-1]
            filehandle = io.TextIOWrapper(zip.open(insideName, 'r'))
        else:
            filehandle.seek(0)

        data = np.empty((lines, ), [
            ('M', np.uint8),
            ('I', np.uint8),
            ('linecenter', np.float64),
            ('S', np.float64),
            ('Acoeff', np.float64),
            ('gamma-air', np.float64),
            ('gamma-self', np.float64),
            ('Epp', np.float64),
            ('N', np.float64),
            ('delta', np.float64),
            ('Vp', 'a15'),
            ('Vpp', 'a15'),
            ('Qp', 'a15'),
            ('Qpp', 'a15'),
            ('Ierr', 'a6'),
            ('Iref', 'a12'),
            ('flag', 'a1'),
            ('gp', np.float64),
            ('gpp', np.float64)])

        for i in range(lines):
            line = filehandle.readline()
            data[i] = (np.uint8(line[0:2]),
                       np.uint8(line[2]),
                       np.float64(line[3:15]),
                       np.float64(line[15:25]),
                       np.float64(line[25:35]),
                       np.float64(line[35:40]),
                       np.float64(line[40:45]),
                       np.float64(line[45:55]),
                       np.float64(line[55:59]),
                       np.float64(line[59:67]),
                       line[67:82],
                       line[82:97],
                       line[97:112],
                       line[112:127],
                       line[127:133],
                       line[133:145],
                       line[145],
                       np.float64(line[146:153]),
                       np.float64(line[153:160]))

        if filename.endswith('.zip'):
            zip.close()
        else:
            filehandle.close()

        return data

    def get_data(self, wn_min, wn_max, Smin=None):
        """Returns HITRAN data for specific molecule from wn_min to
        wn_max and reciprocal of its molar volume.
        """
        data = self.parread(self.fpath)
        data = data[(data['linecenter'] > wn_min) &
                    (data['linecenter'] < wn_max)]
        if Smin:
            data = data[data['S'] >= Smin]

        return data

    @staticmethod
    def pshift(data, p):
        """
        Line position pressure shift at pressure 'p'.
        """
        p = p*1/760
        return data['delta']*p

    @staticmethod
    def pressureHWHM(data, p):
        """
        Pressure-induced HWHM at presssure 'p'.
        """
        p = p*1/760
        return data['gamma-air']*p

    def dopplerHWHM(self, data, T):
        """
        Doppler line width HWHM at temperature 'T'.
        """
        m = self.mass[self.molecule]
        v = np.sqrt(2*k*T/m)
        return nu2wn(np.sqrt(np.log(2))*data['linecenter']*1e2*v)

    def cross_section(self, wns, p, T, disp=False, Smin=None,
                      overrride_data=None, homocoeff=1.0, shiftcoeff=1.0,
                      vipa=0.0, sd=None):
        """Returns cross_sections for wavenumbers 'wns' at pressure 'p'
        with line strength cut off 'Smin', based on reciprocal molar
        volume 'recip_vol' and HITRAN data in 'data'

        homocoeff and shiftcoeff is for adjusting broadening and
        shifting coefficients when, e.g. measurements from HITRAN were
        performed in nitrogen, but ours were in argon.

        vipa is HWHM of etalon resolution in Hz.
        """
        nA = self.recip_volumes[self.molecule]
        m = self.mass[self.molecule]
        v = np.sqrt(2*k*T/m)
        p *= 1/760
        freqs = wn2nu(wns)
        lims = sorted((wns[0], wns[-1]))
        if overrride_data is not None:
            data = overrride_data
        else:
            data = self.get_data(lims[0], lims[1], Smin)
        lc_freqs = wn2nu(data['linecenter'])
        pshift = shiftcoeff*wn2nu(data['delta'])*p
        homo = homocoeff*wn2nu(data['gamma-air'])*p + vipa
        doppler = np.sqrt(np.log(2))*data['linecenter']*1e2*v
        if disp:
            xs = np.zeros(wns.size, dtype=np.complex)
        else:
            xs = np.zeros(wns.size)
        for i in range(lc_freqs.size):
            width = 0.5346*homo[i]+np.sqrt(0.2166*homo[i]**2+doppler[i]**2)
            center = lc_freqs[i]+pshift[i]
            span = (freqs > center-50*width) & (freqs < center+50*width)
            if sd is None:
                line = complex_voigt(freqs[span], center, doppler[i],
                                     homo[i]).astype(xs.dtype)
                xs[span] += data['S'][i]*273.15/T*nA*p*line
            else:
                line = sd[0](lc_freqs[i], doppler[i], homo[i], sd[1],
                             pshift[i], sd[2],
                             freqs[span]).astype(xs.dtype)
                xs[span] += data['S'][i]*273.15/T*nA*p*line

        return xs
