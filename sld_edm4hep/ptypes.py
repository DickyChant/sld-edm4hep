"""SLD Monte Carlo particle-type codes -> PDG.

The SLD ``MCPART.ptype`` word uses the JETSET 6.3 particle numbering (with SLD
extensions) for codes below 1000 and ``1000000 + KF`` (JETSET 7.3 numbering) above.
A positive code is the particle, a negative code its antiparticle.

The table was established empirically from the released mini-DST sample (median
mass from stored energy, charge, decay products, GEANT provenance) and checked
against the JETSET 6.3 numbering where it is known (7-12 leptons, 17-38 mesons,
41-75 baryons, 500+ flavours with u=1, d=2 as in JETSET 6.x).

Fields: name, PDG code of the positive SLD code, mass in GeV, charge of the
positive code in units of e, self-conjugate flag.  Codes whose PDG identity could
not be fixed from the corpus are kept with PDG 0 and a descriptive name; their raw
SLD code is always preserved in the ``MCParticlesSLDType`` collection.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class PType:
    name: str
    pdg: int
    mass: float
    charge: float
    self_conjugate: bool = False
    note: str = ""


_T = PType

# fmt: off
TABLE: dict[int, PType] = {
    # gauge bosons and leptons (JETSET 6.3: 1 gamma, 7 e-, 8 nu_e, 9 mu-, 10 nu_mu, 11 tau-, 12 nu_tau)
    1:  _T("gamma",  22, 0.0,       0.0, True),
    7:  _T("e-",     11, 0.000511, -1.0),
    8:  _T("nu_e",   12, 0.0,       0.0),
    9:  _T("mu-",    13, 0.10566,  -1.0),
    10: _T("nu_mu",  14, 0.0,       0.0),
    11: _T("tau-",   15, 1.77686,  -1.0),
    12: _T("nu_tau", 16, 0.0,       0.0),
    # light and heavy mesons (JETSET 6.3 numbering)
    17: _T("pi+",    211, 0.13957,  1.0),
    18: _T("K+",     321, 0.49368,  1.0),
    19: _T("K0",     311, 0.49761,  0.0),
    20: _T("D0",     421, 1.86484,  0.0),
    21: _T("D+",     411, 1.86966,  1.0),
    22: _T("D_s+",   431, 1.96835,  1.0),
    23: _T("pi0",    111, 0.13498,  0.0, True),
    24: _T("eta",    221, 0.54786,  0.0, True),
    25: _T("eta'",   331, 0.95778,  0.0, True),
    26: _T("eta_c",  441, 2.9839,   0.0, True),
    27: _T("rho+",   213, 0.77526,  1.0),
    28: _T("K*+",    323, 0.89167,  1.0),
    29: _T("K*0",    313, 0.89555,  0.0),
    30: _T("D*0",    423, 2.00685,  0.0),
    31: _T("D*+",    413, 2.01026,  1.0),
    32: _T("D_s*+",  433, 2.1122,   1.0),
    33: _T("rho0",   113, 0.77526,  0.0, True),
    34: _T("omega",  223, 0.78266,  0.0, True),
    35: _T("phi",    333, 1.01946,  0.0, True),
    36: _T("J/psi",  443, 3.09690,  0.0, True),
    37: _T("K_S0",   310, 0.49761,  0.0, True),
    38: _T("K_L0",   130, 0.49761,  0.0, True),
    # SLD extensions identified by mass and decay
    6:  _T("a_1+",    20213, 1.23,   1.0, note="m~1.22, ->rho pi"),
    16: _T("a_1(0)",  20113, 1.23,   0.0, True, note="m~1.20, ->rho pi"),
    94: _T("K_1(1270)+", 10323, 1.272, 1.0, note="m~1.29, ->K pi, K rho, K* pi"),
    97: _T("K_1(1400)0", 20313, 1.403, 0.0, note="m~1.41, ->K* pi"),
    # baryons (JETSET 6.3 numbering)
    41: _T("p",       2212, 0.93827,  1.0),
    42: _T("n",       2112, 0.93957,  0.0),
    43: _T("Sigma+",  3222, 1.18937,  1.0),
    44: _T("Sigma0",  3212, 1.19264,  0.0),
    45: _T("Sigma-",  3112, 1.19745, -1.0),
    46: _T("Xi0",     3322, 1.31486,  0.0),
    47: _T("Xi-",     3312, 1.32171, -1.0),
    48: _T("Sigma_c++", 4222, 2.4530, 2.0),
    49: _T("Sigma_c+",  4212, 2.4529, 1.0),
    50: _T("Sigma_c0",  4112, 2.4538, 0.0),
    51: _T("Xi_c'+",    4322, 2.5784, 1.0, note="m=2.55 in JETSET, ->Xi_c+ gamma"),
    57: _T("Lambda",    3122, 1.11568, 0.0),
    58: _T("Lambda_c+", 4122, 2.28646, 1.0),
    59: _T("Xi_c+",     4232, 2.4679,  1.0),
    60: _T("Xi_c0",     4132, 2.4710,  0.0),
    61: _T("Delta++",   2224, 1.232,   2.0),
    62: _T("Delta+",    2214, 1.232,   1.0),
    63: _T("Delta0",    2114, 1.232,   0.0),
    64: _T("Delta-",    1114, 1.232,  -1.0),
    65: _T("Sigma*+",   3224, 1.3828,  1.0),
    66: _T("Sigma*0",   3214, 1.3837,  0.0),
    67: _T("Sigma*-",   3114, 1.3872, -1.0),
    68: _T("Xi*0",      3324, 1.5318,  0.0),
    69: _T("Xi*-",      3314, 1.5350, -1.0),
    70: _T("Omega-",    3334, 1.67245, -1.0),
    71: _T("Sigma_c*++", 4224, 2.5184, 2.0),
    72: _T("Sigma_c*+",  4214, 2.5175, 1.0),
    73: _T("Sigma_c*0",  4114, 2.5180, 0.0),
    74: _T("Xi_c*+",     4324, 2.6459, 1.0),
    75: _T("Xi_c*0",     4314, 2.6461, 0.0),
    # b hadrons (SLD extension of the 6.3 table)
    101: _T("B+",      521, 5.27934, 1.0),
    102: _T("B0",      511, 5.27965, 0.0),
    103: _T("B_s0",    531, 5.36688, 0.0),
    123: _T("B*+",     523, 5.3247,  1.0),
    124: _T("B*0",     513, 5.3247,  0.0),
    125: _T("B_s*0",   533, 5.4154,  0.0),
    145: _T("Sigma_b+",  5222, 5.81056, 1.0),
    146: _T("Sigma_b0",  5212, 5.8135,  0.0),
    147: _T("Sigma_b-",  5112, 5.81564, -1.0),
    241: _T("Lambda_b0", 5122, 5.61960, 0.0),
    242: _T("Xi_b0",     5232, 5.7919,  0.0),
    243: _T("Xi_b-",     5132, 5.7970, -1.0),
    293: _T("Sigma_b*+", 5224, 5.83032, 1.0),
    294: _T("Sigma_b*0", 5214, 5.8327,  0.0),
    295: _T("Sigma_b*-", 5114, 5.83474, -1.0),
    296: _T("Xi_b*0",    5324, 5.9520,  0.0),
    297: _T("Xi_b*-",    5314, 5.9555, -1.0),
    # excited charm mesons: m 2.38-2.63 GeV, decays to D(*) pi / D(*) K; the JETSET 6.3
    # extension table that fixes which state is which is not in the corpus.
    109: _T("D**(2420)0?", 0, 2.4228, 0.0, note="->D* pi, D pi; likely D_1(2420)0"),
    110: _T("D**(2380)+?", 0, 2.384,  1.0, note="->D pi"),
    111: _T("D**(2450)0?", 0, 2.45,   0.0, note="->D* pi"),
    112: _T("D**(2450)+?", 0, 2.45,   1.0, note="->D* pi"),
    113: _T("D**(2420)0?", 0, 2.4233, 0.0, note="->D* pi"),
    114: _T("D**(2420)+?", 0, 2.4224, 1.0, note="->D* pi"),
    115: _T("D**(2460)0?", 0, 2.458,  0.0, note="->D(*) pi"),
    116: _T("D**(2460)+?", 0, 2.46,   1.0, note="->D(*) pi"),
    117: _T("D_s**(2535)+?", 0, 2.535, 1.0, note="->D K"),
    118: _T("D_s**(2620)+?", 0, 2.62,  1.0, note="->D* K"),
    119: _T("D_s**(2535)+?", 0, 2.535, 1.0, note="->D* K"),
    120: _T("D_s**(2570)+?", 0, 2.572, 1.0, note="->D K"),
    # partons: JETSET 6.x flavour order u=1, d=2, s=3, c=4, b=5; 500 = gluon
    500: _T("g", 21, 0.0,  0.0, True),
    501: _T("u",  2, 0.0,  2.0 / 3.0),
    502: _T("d",  1, 0.0, -1.0 / 3.0),
    503: _T("s",  3, 0.0, -1.0 / 3.0),
    504: _T("c",  4, 1.35, 2.0 / 3.0),
    505: _T("b",  5, 5.0, -1.0 / 3.0),
    # GEANT-produced nuclei (detector secondaries), identified by mass
    91: _T("deuteron", 1000010020, 1.87561, 1.0),
    92: _T("triton",   1000010030, 2.80892, 1.0),
    93: _T("alpha",    1000020040, 3.72738, 2.0),
    # 1000000 + JETSET 7.3 KF
    1000091: _T("cluster",   91, 0.0, 0.0, True, note="JETSET 7.3 KF 91; mass from stored E"),
    1000092: _T("string",    92, 0.0, 0.0, True, note="JETSET 7.3 KF 92; mass from stored E"),
    1000094: _T("CMshower",  94, 0.0, 0.0, True, note="JETSET 7.3 KF 94; mass from stored E"),
    1000445: _T("chi_c2",    445, 3.55617, 0.0, True),
    1010221: _T("f_0(980)",  9010221, 0.990, 0.0, True),
    1010321: _T("K_0*(1430)+", 10321, 1.425, 1.0),
    1010441: _T("chi_c0",    10441, 3.41471, 0.0, True),
    1020443: _T("chi_c1",    20443, 3.51067, 0.0, True),
    1030443: _T("psi(2S)",   100443, 3.68610, 0.0, True),
    # unidentified: m 3.4-3.8 GeV, decays Lambda_c+ + (e, gamma, Lambda); a few per 10^5 events
    1000073: _T("unknown(1000073)", 0, 3.47, 0.0, note="->Lambda_c e"),
    1000074: _T("unknown(1000074)", 0, 3.47, 0.0, note="->Lambda_c- e+"),
    1000075: _T("unknown(1000075)", 0, 3.50, 0.0, note="->Lambda_c gamma"),
    1000076: _T("unknown(1000076)", 0, 3.40, 0.0, note="->Lambda_c- gamma"),
    1000077: _T("unknown(1000077)", 0, 3.84, 0.0, note="->Lambda_c Lambda"),
    1000078: _T("unknown(1000078)", 0, 3.82, 0.0, note="->Lambda_c- Lambda"),
}
# fmt: on

# Codes for which the mass must be taken from the stored energy (virtual partons and
# colour-singlet systems); for every other known code the tabulated mass is used.
MASS_FROM_ENERGY = frozenset({500, 501, 502, 503, 504, 505, 1000091, 1000092, 1000094})

# Mass table used by jazelle_reader when the energy word was omitted in the bank.
# It is wrong for this corpus (codes 7-10 are leptons, not K/p/n), and the value it
# inserted must be recognised so that the energy can be recomputed.
READER_MASS = {
    1: 0.000511, 2: 0.0, 3: 0.0, 4: 0.10566, 5: 0.13957, 6: 0.13498,
    7: 0.49368, 8: 0.49761, 9: 0.93827, 10: 0.93957,
}


def lookup(code: int) -> tuple[int, float, float, str]:
    """Return (pdg, mass, charge, name) for a signed SLD code.

    Unknown codes give (0, 0.0, 0.0, "unknown(code)").
    """
    entry = TABLE.get(abs(code))
    if entry is None:
        return 0, 0.0, 0.0, f"unknown({code})"
    if code < 0 and not entry.self_conjugate:
        return -entry.pdg, entry.mass, -entry.charge, entry.name + "bar"
    return entry.pdg, entry.mass, entry.charge, entry.name
