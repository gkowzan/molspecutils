"""Parse MIRS output files."""
from collections import namedtuple

MIRSRow = namedtuple("MIRSRow", "nu sw Ppp Jpp rovibpp eigenpp Pp Jp rovibp eigenp Kpp Kp Vpp Vp nupp vibsympp nup vibsymp elower")

def parse_prd_row(row):
    return MIRSRow(
        nu=float(row[:11].strip()),   # line frequency
        sw=float(row[11:21].strip()), # line strength
        Ppp=int(row[21:24].strip()), # lower polyad number
        Jpp=int(row[24:27].strip()), # lower J number
        rovibpp=row[27:31].strip(),  # rovib symmetry
        eigenpp=int(row[31:35].strip()), # eigenvalue number within sym block
        Pp=int(row[35:38].strip()),
        Jp=int(row[38:41].strip()),
        rovibp = row[41:45].strip(),
        eigenp = row[45:49].strip(),
        Kpp = int(row[51:53].strip()),
        Kp = int(row[54:56].strip()),
        Vpp = int(row[58:61].strip()), # no. of vib levels in polyad
        Vp = int(row[61:65].strip()),
        nupp = row[86:92],
        vibsympp = row[93:95].strip(),
        nup = row[96:102],
        vibsymp = row[103:105].strip(),
        elower = float(row[106:117].strip())
    )


def rot_kwargs(mirs_row):
    """Lower and upper rot state kwargs."""
    return (
        dict(j=mirs_row.Jpp, k=mirs_row.Kpp,
             l=0, sym=mirs_row.rovibpp, f=0),
        dict(j=mirs_row.Jp, k=mirs_row.Kp,
             l=0, sym=mirs_row.rovibp, f=0)
    )


def vib_kwargs(mirs_row):
    vibs = ['nu'+str(i) for i in range(1, 7)]
    vibpp = dict(zip(vibs, [int(c) for c in mirs_row.nupp]))
    vibp = dict(zip(vibs, [int(c) for c in mirs_row.nup]))

    return vibpp, vibp
