"""Multiscale stiffness homogenization of native and transparent wood.

Procedural port of 02_stiffness/hom_TranspWood_model.m (lumped cell wall,
no ML/S2 split). Kelvin-Mandel notation, order 11,22,33,23,13,12, GPa.
RVE chain: pn -> cel -> cw -> EWuc/LWuc -> vesselwood -> ringwood -> ray
-> clearwood.
"""
import itertools
import json
import os
from math import gamma, sqrt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
I6 = np.eye(6)
IVOL = np.zeros((6, 6)); IVOL[:3, :3] = 1 / 3
IDEV = I6 - IVOL


# ---------------------------------------------------------------- tensors
def C_kmu(k, mu):
    return 3 * k * IVOL + 2 * mu * IDEV


def C_Enu(E, nu):
    return C_kmu(E / (3 * (1 - 2 * nu)), E / (2 * (1 + nu)))


def Q4_bp(azi, zeni):
    """Rotation of Mandel 4th-order tensors, Pichler convention (fun_Q4_bp)."""
    ca, sa, cz, sz = np.cos(azi), np.sin(azi), np.cos(zeni), np.sin(zeni)
    Q = np.array([[ca * cz, sa * cz, -sz], [-sa, ca, 0.0], [ca * sz, sa * sz, cz]])
    p = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
    r2 = sqrt(2)
    M = np.empty((6, 6))
    for a, (i, j) in enumerate(p):
        for b, (k, l) in enumerate(p):
            if a < 3 and b < 3:
                M[a, b] = Q[i, k] ** 2
            elif a < 3:
                M[a, b] = r2 * Q[i, k] * Q[i, l]
            elif b < 3:
                M[a, b] = r2 * Q[i, k] * Q[j, k]
            else:
                M[a, b] = Q[i, k] * Q[j, l] + Q[i, l] * Q[j, k]
    return M


def isotropy(C, tol=1e-8):
    """'iso', 'transiso' (axis e3) or 'aniso' - same test as fun_check_isotropy."""
    la, mu = C[0, 1], C[3, 3] / 2
    Ci = la * np.ones((6, 6)) * (np.arange(6) < 3)[:, None] * (np.arange(6) < 3) + 2 * mu * I6
    if np.abs(Ci - C).sum() < tol:
        return 'iso'
    Ct = np.zeros((6, 6))
    Ct[0, 0] = Ct[1, 1] = C[0, 0]
    Ct[0, 1] = Ct[1, 0] = C[0, 0] - C[5, 5]
    Ct[0, 2] = Ct[1, 2] = Ct[2, 0] = Ct[2, 1] = C[0, 2]
    Ct[2, 2], Ct[3, 3], Ct[4, 4], Ct[5, 5] = C[2, 2], C[3, 3], C[3, 3], C[5, 5]
    return 'transiso' if np.abs(Ct - C).sum() < tol else 'aniso'


# ------------------------------------------------------------ Hill tensors
_PAIRS = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
_W = np.array([1, 1, 1, sqrt(2), sqrt(2), sqrt(2)])
_GL = np.polynomial.legendre.leggauss(16)


def _to_tensor(M):
    T = np.zeros((3, 3, 3, 3))
    for a, (i, j) in enumerate(_PAIRS):
        for b, (k, l) in enumerate(_PAIRS):
            v = M[a, b] / (_W[a] * _W[b])
            T[i, j, k, l] = T[j, i, k, l] = T[i, j, l, k] = T[j, i, l, k] = v
    return T


def _panels(ratio):
    """Gauss-Legendre nodes on [0, pi/2], refined where atan(ratio*10^k)."""
    k = np.arange(-2, 3) if not 0.2 < ratio < 5 else np.array([0])
    br = np.arctan(ratio * 10.0 ** k)
    br = np.unique(np.r_[0, br[(br > 1e-12) & (br < np.pi / 2 - 1e-12)], np.pi / 2])
    x, w = _GL
    a, b = br[:-1, None], br[1:, None]
    return ((b - a) / 2 * x + (a + b) / 2).ravel(), ((b - a) / 2 * w).ravel()


def hill_numeric(C0, asp, slend):
    """Hill tensor of an ellipsoid (a1/a2 = asp, a1/a3 = slend, axis e3) in
    an arbitrary anisotropic matrix, same integrand as fun_P_ellipsoid_aniso
    but one vectorised fixed quadrature for all 21 components."""
    th, wt = _panels(slend)                       # hemisphere, integrand even in n
    u, wu = _panels(1 / asp)
    ph = np.r_[u, np.pi - u, np.pi + u, 2 * np.pi - u]
    wp = np.r_[wu, wu, wu, wu]
    TH, PH = np.meshgrid(th, ph, indexing='ij')
    w = (np.outer(wt, wp) * np.sin(TH)).ravel()
    n = np.stack([np.cos(PH) * np.sin(TH), asp * np.sin(PH) * np.sin(TH),
                  slend * np.cos(TH)], -1).reshape(-1, 3)
    K = np.einsum('ijkl,pj,pl->pik', _to_tensor(C0), n, n, optimize=True)
    a, b, c, d, e, f = K[:, 0, 0], K[:, 1, 1], K[:, 2, 2], K[:, 1, 2], K[:, 0, 2], K[:, 0, 1]
    A11, A22, A33 = b * c - d * d, a * c - e * e, a * b - f * f        # adjugate of the
    A23, A13, A12 = e * f - a * d, f * d - b * e, d * e - c * f        # symmetric K
    s = w / (a * A11 + f * A12 + e * A13)
    Kinv = np.stack([A11, A12, A13, A12, A22, A23, A13, A23, A33], -1).reshape(-1, 3, 3) * s[:, None, None]
    T = np.einsum('pik,pj,pl->ijkl', Kinv, n, n, optimize=True) * (2 / (4 * np.pi))
    P = np.empty((6, 6))
    for a, (i, j) in enumerate(_PAIRS):
        for b, (k, l) in enumerate(_PAIRS):
            P[a, b] = _W[a] * _W[b] * (T[i, j, k, l] + T[j, i, k, l]
                                       + T[i, j, l, k] + T[j, i, l, k]) / 4
    return P


def hill_iso(C0, slend):
    """Analytic Hill tensor of sphere / spheroid in isotropic matrix."""
    D0 = np.linalg.inv(C0)
    E0, mu0 = 1 / D0[0, 0], C0[3, 3] / 2
    if slend == 1:                                # fun_P_sphere_iso
        k0 = E0 * mu0 / (3 * (3 * mu0 - E0))
        S = 3 * k0 / (3 * k0 + 4 * mu0) * IVOL + 6 * (k0 + 2 * mu0) / (5 * (3 * k0 + 4 * mu0)) * IDEV
        return S @ np.linalg.inv(C0)
    nu, s = E0 / (2 * mu0) - 1, 1 / slend          # fun_P_spheroid_iso
    if s > 1:
        g = s / (s * s - 1) ** 1.5 * (s * sqrt(s * s - 1) - np.arccosh(s))
    else:
        g = s / (1 - s * s) ** 1.5 * (np.arccos(s) - s * sqrt(1 - s * s))
    q, c = s * s / (s * s - 1), 1 / (s * s - 1)
    S11 = 3 / (8 * (1 - nu)) * q + 1 / (4 * (1 - nu)) * (1 - 2 * nu - 9 / 4 * c) * g
    S33 = 1 / (2 * (1 - nu)) * (1 - 2 * nu + (3 * s * s - 1) * c - (1 - 2 * nu + 3 * q) * g)
    S12 = 1 / (4 * (1 - nu)) * (q / 2 - (1 - 2 * nu + 3 / 4 * c) * g)
    S13 = -1 / (2 * (1 - nu)) * q + 1 / (4 * (1 - nu)) * (3 * q - (1 - 2 * nu)) * g
    S31 = 1 / (2 * (1 - nu)) * (-1 + 2 * nu + 1 / (1 - s * s)) + 1 / (4 * (1 - nu)) * (2 * (1 - 2 * nu) - 3 / (1 - s * s)) * g
    S66 = 1 / (8 * (1 - nu)) * q + 1 / (4 * (1 - nu)) * (1 - 2 * nu - 3 / 4 * c) * g
    S44 = 1 / (4 * (1 - nu)) * (1 - 2 * nu - (s * s + 1) * c) - 1 / (8 * (1 - nu)) * (1 - 2 * nu - 3 * (s * s + 1) * c) * g
    S = np.array([[S11, S12, S13, 0, 0, 0], [S12, S11, S13, 0, 0, 0], [S31, S31, S33, 0, 0, 0],
                  [0, 0, 0, 2 * S44, 0, 0], [0, 0, 0, 0, 2 * S44, 0], [0, 0, 0, 0, 0, 2 * S66]])
    return S @ np.linalg.inv(C0)


def R_spheroid(sr, alpha, beta):
    """Imperfect-interface tensor (fun_R_spheroid, Rao 2017)."""
    rho = 1 / sr
    if rho == 1:
        return alpha * IVOL + (3 * alpha + 2 * beta) / 5 * IDEV
    r2 = complex(rho * rho - 1)                   # oblate branch is complex, as in MATLAB
    at = np.arcsin(np.sqrt(r2) / rho) if rho > 1 else np.arcsinh(np.sqrt(r2) / rho)
    P33 = 1.5 * (rho / r2 ** 1.5 * at - 1 / (rho * r2))
    P11 = 0.75 * (rho * (rho ** 2 - 2) / r2 ** 1.5 * at + rho / r2)
    Q33 = 1.5 * ((2 * rho ** 2 + 1) / (rho * r2 ** 2) - 3 * rho / r2 ** 2.5 * at)
    Q13 = 0.75 * (rho * (rho ** 2 + 2) / r2 ** 2.5 * at - 3 * rho / r2 ** 2)
    Q11 = 9 / 16 * (rho * (2 + rho ** 2) / r2 ** 2 + rho ** 3 * (rho ** 2 - 4) / r2 ** 2.5 * at)
    P11, P33, Q11, Q13, Q33 = (np.real(v) for v in (P11, P33, Q11, Q13, Q33))
    P = np.diag([P11, P11, P33, P33 / 2, P33 / 2, P11])
    Q = np.zeros((6, 6))
    Q[:3, :3] = [[Q11, Q11 / 3, Q13], [Q11 / 3, Q11, Q13], [Q13, Q13, Q33]]
    return alpha * P + (beta - alpha) * Q


# ---------------------------------------------------------- homogenization
def phase(C, vol, asp=1.0, slend=1.0, ori=None, IF=None, matrix=False):
    return dict(C=C, vol=vol, asp=asp, slend=slend, ori=ori, IF=IF, matrix=matrix)


def _hill(C0, ph):
    if isotropy(C0) == 'iso' and ph['asp'] == 1:
        return hill_iso(C0, ph['slend'])
    return hill_numeric(C0, ph['asp'], ph['slend'])


def _concentration(ph, C0):
    """Dilute strain concentration A_inf, Chom and E_inf^-1 contributions."""
    if ph['matrix']:
        return ph['C'], I6
    Q = Q4_bp(*ph['ori']) if ph['ori'] is not None else I6
    P = Q.T @ _hill(Q @ C0 @ Q.T, ph) @ Q
    R = np.zeros((6, 6)) if ph['IF'] is None else Q.T @ R_spheroid(ph['slend'], *ph['IF']) @ Q
    Ci = Q.T @ ph['C'] @ Q
    Ainf = np.linalg.inv(I6 + P @ (Ci - C0) + (I6 - P @ C0) @ R @ Ci)
    return Ci @ Ainf, (I6 + R @ Ci) @ Ainf


def hom(phases, scheme, tol=8):
    """Mori-Tanaka ('MT', one phase flagged matrix) or self-consistent ('SCS')."""
    V = sum(p['vol'] for p in phases)
    f = [p['vol'] / V for p in phases]

    def step(C0):
        CA, EA = np.zeros((6, 6)), np.zeros((6, 6))
        for fi, p in zip(f, phases):
            if fi == 0:                           # e.g. vessels in softwood
                continue
            a, e = _concentration(p, C0)
            CA += fi * a
            EA += fi * e
        return CA @ np.linalg.inv(EA)

    if scheme == 'MT':
        C = step(next(p['C'] for p in phases if p['matrix']))
    else:
        C0 = sum(fi * p['C'] for fi, p in zip(f, phases))       # Voigt start
        err, it, tolerance = 1.0, 0, 10.0 ** -tol / 100
        while err > tolerance and it < tol ** 1.5:
            it += 1
            C = step(C0)
            err = np.abs(C - C0).sum() / C[:3, :3].mean()
            C0 = C
    return C, V                                   # raw; callers pass on sym(C)


# ------------------------------------------- data (all in this folder, JSON)
def _load(fn):
    with open(os.path.join(HERE, fn), encoding='utf8') as f:
        return json.load(f)


def _iso(d):
    return C_Enu(d['E_GPa'], d['nu']) if 'E_GPa' in d else C_kmu(d['k_GPa'], d['mu_GPa'])


PH = _load('phases.json')                         # cell-wall constituents
WOODS = _load('woods.json')                       # species data (from data_Wood.xlsx)
POLYMERS = {k: (v['E_GPa'], v['nu']) for k, v in _load('polymers.json').items() if k[0] != '_'}
CRYCEL = {k: np.array(v['C_GPa']) for k, v in PH['crystalline_cellulose'].items() if isinstance(v, dict)}
AMCEL, HEMCEL, LIGNIN = (_iso(PH[k]) for k in ('amorphous_cellulose', 'hemicellulose', 'lignin'))
PORE = np.zeros((6, 6))
RHO = dict(crycel=PH['crystalline_cellulose']['density'], amcel=PH['amorphous_cellulose']['density'],
           hemcel=PH['hemicellulose']['density'], lignin=PH['lignin']['density'],
           extr=PH['extractives']['density'])


def defaults(name):
    """Model inputs prefilled from one species (editable afterwards)."""
    w = WOODS[name]
    nan = lambda v: np.nan if v is None else float(v)
    ewlw, tol_ew = [nan(v) for v in w['EWLW_um'].values()], 1e-10
    if 'EWLW_surrogate' in w and not np.all(np.isfinite(ewlw)):   # e.g. Cedar -> Douglas Fir
        ewlw, tol_ew = [nan(v) for v in WOODS[w['EWLW_surrogate']]['EWLW_um'].values()], 0.05
    c = w['chemistry_dry_mass']
    return dict(wood=name, density=float(w['density_kg_m3']['mean']), moisture=12.0, CI=float(w['CI']),
                MFA=[float(w['MFA_deg']['mean'])], vessel=w['cells_pct']['vessels'] / 100,
                ray=w['cells_pct']['rays'] / 100,
                mass_fr=[float(c[k]) for k in ('cellulose', 'hemicellulose', 'lignin', 'extractives')],
                EWLW=ewlw, EWLW_tol=tol_ew, cell_ar=1 / 100, n_super=4.0, state=1, polymer='none',
                E_poly=np.nan, nu_poly=np.nan, IF=(0.0, 0.0), nfam=20, slend_mf=1e-20,
                cellulose='Dri2014_TI', tol=8, swelling=False, swelling_pct=0.0)


def volfrac(p):
    """Lumped-cell-wall volume fractions (comp_volfrac_TranspWood), plus wall swelling."""
    m = dict(zip(['cel', 'hemcel', 'lignin', 'extr'], p['mass_fr']))
    m['hemcel'] += 1 - sum(m.values())
    m['crycel'], m['amcel'] = m['cel'] * p['CI'], m['cel'] * (1 - p['CI'])
    H = p['moisture'] / 100
    w = {k: m[k] / (1 + H) for k in ('crycel', 'amcel', 'hemcel', 'lignin', 'extr')}
    wH = H / (1 + H)
    rhoH = -1.17 * H ** 3 + 2.35 * H ** 2 - 1.22 * H + 1.3
    rbcw = 1 / (sum(w[k] / RHO[k] for k in w) + wH / rhoH)
    fcw = {k: w[k] * rbcw / RHO[k] for k in ('crycel', 'amcel', 'hemcel', 'lignin')}
    fcw['pore'] = wH * rbcw / rhoH + w['extr'] * rbcw / RHO['extr']
    s = sum(fcw.values())
    fcw = {k: v / s for k, v in fcw.items()}
    bcw = p['density'] / 1000 / rbcw
    fib = 1 - p['vessel'] - p['ray']
    nonv = p['ray'] + fib
    lum = 1 - bcw - p['vessel']
    fclw = dict(bcw=bcw, vessel=p['vessel'], fib=fib, bcw_ray=bcw * p['ray'] / nonv,
                bcw_fib=bcw * fib / nonv, lum_ray=lum * p['ray'] / nonv, lum_fib=lum * fib / nonv)
    por = fclw['lum_fib'] / (fclw['bcw_fib'] + fclw['lum_fib'])
    # Swelling after Schindler et al. (Compos. Sci. Technol., Supplementary Information):
    # the lumped wall (ML+S2 together) thickens by s at fixed cell-centre
    # distance, i.e. inwards into the lumen (as comp_volfrac_TranspWood). Square
    # cells: wall/half-width = 1-sqrt(porosity); EW/LW cells: walls 2W -> 2W(1+s).
    s = p['swelling_pct'] / 100 if p['swelling'] else 0.0
    square = lambda phi: (1 - (1 + s) * (1 - sqrt(phi))) ** 2 if (1 + s) * (1 - sqrt(phi)) < 1 else -1
    e = p['EWLW']
    if not np.isfinite(e[0]):
        EW, fEW, fLW = 1.0, square(por), 0.0
    else:
        sf = gamma(1 + 1 / p['n_super']) ** 2 / gamma(1 + 2 / p['n_super'])
        lt = e[0] - e[3] + e[1]
        cells = lambda s: (max(e[0] - e[1] * s, 0) ** 2 * sf / (e[0] + e[1]) ** 2,
                           max(lt - e[3] * s, 0) * max(e[2] - e[3] * s, 0) * sf / ((lt + e[3]) * (e[2] + e[3])))
        fEW, fLW = cells(0)
        EW = (por - fLW) / (fEW - fLW)
        if not -p['EWLW_tol'] <= EW <= 1 + p['EWLW_tol']:
            raise ValueError(f'earlywood fraction {EW:.4g} outside [0,1]')
        EW = min(max(EW, 0.0), 1.0)
        fEW, fLW = cells(s)                       # anatomy (EW share) fixed, lumens shrink
    sw = dict(pct=100 * s, porosity_fib=por, k_wall=1.0)
    if s:
        por_f, por_r = EW * fEW + (1 - EW) * fLW, square(por)
        if min(por_f, por_r, fEW) <= 0:
            raise ValueError(f'{100 * s:g} % wall swelling closes the lumens')
        nf, nr = fclw['bcw_fib'] + fclw['lum_fib'], fclw['bcw_ray'] + fclw['lum_ray']
        fclw.update(bcw_fib=(1 - por_f) * nf, lum_fib=por_f * nf, bcw_ray=(1 - por_r) * nr, lum_ray=por_r * nr)
        k = (fclw['bcw_fib'] + fclw['bcw_ray']) / bcw          # wall volume ratio swollen/native
        fclw['bcw'] *= k
        fcw = {n: (x + (k - 1) * (n == 'pore')) / k for n, x in fcw.items()}   # new volume = pore
        sw.update(porosity_fib_swollen=por_f, k_wall=k)
    return dict(m=m, rho_bcw=rbcw, fcw=fcw, fclw=fclw, EW=EW, LW=1 - EW, fEW=fEW, fLW=fLW, swelling=sw)


# ------------------------------------------------------------- the model
def run(p):
    """Homogenize one parameter set p (see defaults()) for every MFA in p['MFA']."""
    v = volfrac(p)
    fcw, fclw = v['fcw'], v['fclw']
    poly = None
    if p['state'] != 1:
        E, nu = POLYMERS.get(p['polymer'], (p['E_poly'], p['nu_poly']))
        if not (np.isfinite(E) and np.isfinite(nu)):
            raise ValueError('polymer E and nu required for states 2 and 3')
        poly = C_Enu(E, nu)
    pore = PORE if poly is None else poly
    lig = poly if p['state'] == 3 else LIGNIN
    tol, IF = p['tol'], tuple(p['IF'])
    L = {}

    def H(name, phases, scheme):                  # store raw Chom, pass on symmetric part
        L[name], V = hom(phases, scheme, tol)
        return (L[name] + L[name].T) / 2, V

    # RVE 1: polymer network, RVE 2: cellulose fibril
    C_pn, V_pn = H('pn', [phase(HEMCEL, fcw['hemcel']), phase(lig, fcw['lignin']),
                      phase(pore, fcw['pore'])], 'SCS')
    C_cel, V_cel = H('cel', [phase(AMCEL, fcw['amcel'], matrix=True),
                        phase(CRYCEL[p['cellulose']], fcw['crycel'], slend=1e-20, ori=(0, 0))], 'MT')
    out = dict(MFA=list(p['MFA']), C_clearwood=[], C_cellwall=[], levels=[], vol=v)
    for mfa in p['MFA']:
        a = mfa * np.pi / 180
        L = dict(pn=L['pn'], cel=L['cel'])
        # RVE 3: cell wall, fibrils in nfam orientation families
        if a == 0:
            cel = [phase(C_cel, V_cel, slend=p['slend_mf'], ori=(0, 0), IF=IF)]
        else:
            n = p['nfam']
            cel = [phase(C_cel, V_cel / n, slend=p['slend_mf'], ori=(i / n * 2 * np.pi, a), IF=IF)
                   for i in range(n)]
        C_cw, _ = H('cw', cel + [phase(C_pn, V_pn, matrix=True)], 'MT')
        # RVE 4a/b: earlywood / latewood unit cells
        ew = v['EW'] * fclw['fib']
        C_ew, V_ew = H('EWuc', [phase(pore, v['fEW'] * ew, slend=p['cell_ar'], ori=(0, 0)),
                          phase(C_cw, (1 - v['fEW']) * ew, matrix=True)], 'MT')
        LW = v['LW'] > 0
        if LW:
            lw = v['LW'] * fclw['fib']
            C_lw, V_lw = H('LWuc', [phase(pore, v['fLW'] * lw, asp=p['EWLW'][2] / p['EWLW'][0],
                                    slend=p['cell_ar'], ori=(0, 0)),
                              phase(C_cw, (1 - v['fLW']) * lw, matrix=True)], 'MT')
        # RVE 5: vesselwood, RVE 6: annual ring
        C_vw, V_vw = H('vesselwood', [phase(pore, fclw['vessel'], slend=1 / 1000, ori=(0, 0)),
                          phase(C_ew, V_ew, matrix=True)], 'MT')
        if LW:
            Q = Q4_bp(0, np.pi / 2)
            rot = lambda C: (Q @ C @ Q.T + (Q @ C @ Q.T).T) / 2
            C_rw, V_rw = H('ringwood', [phase(rot(C_vw), V_vw, slend=1e20, ori=(0, np.pi / 2)),
                              phase(rot(C_lw), V_lw, slend=1e20, ori=(0, np.pi / 2))], 'SCS')
        else:
            C_rw, V_rw = C_vw, V_vw
            L['ringwood'] = L['vesselwood']
        # RVE 4c: ray cell (rays run radially, i.e. 90° from L)
        C_ray, V_ray = H('ray', [phase(pore, fclw['lum_ray'], slend=1 / 5, ori=(0, 0)),   # ray lumen 1/5
                            phase(C_cw, fclw['bcw_ray'], matrix=True)], 'MT')
        # RVE 7: clearwood
        C, _ = H('clearwood', [phase(C_ray, V_ray, asp=10, slend=1e-20, ori=(0, np.pi / 2)),
                    phase(C_rw, V_rw, matrix=True)], 'MT')
        out['C_clearwood'].append(L['clearwood'])
        out['C_cellwall'].append(L['cw'])
        out['levels'].append(dict(L))
    out['C_clearwood'], out['C_cellwall'] = np.array(out['C_clearwood']), np.array(out['C_cellwall'])
    out.update(engineering(out['C_clearwood']))
    out['M_cw'] = np.array([indentation_modulus(C) for C in out['C_cellwall']])
    return out


def engineering(Cs):
    """E_R/T/L, G_TL/RL/RT, nu_TR/LR/TL/LT as in hom_TranspWood_model."""
    r = {k: [] for k in ('E_R', 'E_T', 'E_L', 'G_TL', 'G_RL', 'G_RT', 'nu_TR', 'nu_LR', 'nu_TL', 'nu_LT')}
    for C in Cs:
        D = np.linalg.inv(C)
        ER, ET, EL = 1 / D[0, 0], 1 / D[1, 1], 1 / D[2, 2]
        vals = (ER, ET, EL, 1 / (2 * D[3, 3]), 1 / (2 * D[4, 4]), 1 / (2 * D[5, 5]),
                -D[0, 1] * ET, -D[0, 2] * EL, -D[1, 2] * ET,
                (C[0, 0] * C[1, 2] - C[0, 1] * C[0, 2]) / (C[0, 0] * C[1, 1] - C[0, 1] ** 2))
        for k, x in zip(r, vals):
            r[k].append(x)
    return {k: np.array(x) for k, x in r.items()}


def indentation_modulus(C):
    """Delafargue & Ulm (2004), transversely isotropic about e3."""
    C11, C13, C44, C31 = C[0, 0], C[0, 2], C[3, 3] / 2, sqrt(C[0, 0] * C[2, 2])
    return 2 * sqrt((C31 ** 2 - C13 ** 2) / C11 / (1 / C44 + 2 / (C31 + C13)))


# ---------------------------------------------------- input files and batches
# Flat input keys (GUI fields = input-file keys). Numeric values may be a number,
# a list [a, b, c] or a range string 'start:step:stop' -> all combinations are run.
INPUTS = {'density': 'Density [kg/m³]', 'moisture': 'Moisture [% dry mass]',
          'E_poly': 'E polymer [GPa]', 'nu_poly': 'nu polymer [-]',
          'IF_alpha': 'IF alpha tangential [1/GPa]', 'IF_beta': 'IF beta normal [1/GPa]',
          'swelling_pct': 'Wall thickness increase [%]', 'CI': 'Crystallinity CI [-]',
          'cellulose': 'Cellulose [dry mass fr.]', 'hemicellulose': 'Hemicellulose [dry mass fr.]',
          'lignin': 'Lignin [dry mass fr.]', 'extractives': 'Extractives [dry mass fr.]',
          'vessel_pct': 'Vessels [%]', 'ray_pct': 'Rays [%]', 'L_EW': 'L EW [µm]', 'W2_EW': '2W EW [µm]',
          'L_LW': 'L LW [µm]', 'W2_LW': '2W LW [µm]', 'cell_aspect_ratio': 'Cell aspect ratio',
          'superellipse_n': 'Superellipse n', 'n_families': 'Fibril orientation families',
          'microfibril_slenderness': 'Microfibril slenderness', 'tolerance': 'Tolerance (1-16)'}
STATES = {1: 'native', 2: 'infiltrated', 3: 'delignified + infiltrated'}
OUTPUTS = {'E_R': 'E_R [GPa]', 'E_T': 'E_T [GPa]', 'E_L': 'E_L [GPa]', 'G_TL': 'G_TL [GPa]',
           'G_RL': 'G_RL [GPa]', 'G_RT': 'G_RT [GPa]', 'nu_TR': 'nu_TR', 'nu_LR': 'nu_LR',
           'nu_TL': 'nu_TL', 'nu_LT': 'nu_LT', 'M_cw': 'M_cellwall [GPa]'}
POLY = ('E_poly', 'nu_poly')                      # irrelevant for native wood


def flat_defaults(wood):
    """Complete flat input set for one species (what the GUI shows)."""
    p = defaults(wood)
    E, nu = POLYMERS['HEMA']
    f = dict(wood=wood, states=[1, 2, 3], polymer='HEMA', E_poly=E, nu_poly=nu,
             density=p['density'], moisture=p['moisture'], CI=p['CI'], MFA=p['MFA'],
             vessel_pct=100 * p['vessel'], ray_pct=100 * p['ray'], IF_alpha=0.0, IF_beta=0.0,
             swelling=False, swelling_pct=10.0, cell_aspect_ratio=p['cell_ar'],
             superellipse_n=p['n_super'], n_families=p['nfam'],
             microfibril_slenderness=p['slend_mf'], tolerance=p['tol'], cellulose_material=p['cellulose'])
    f.update(zip(('cellulose', 'hemicellulose', 'lignin', 'extractives'), p['mass_fr']))
    f.update(zip(('L_EW', 'W2_EW', 'L_LW', 'W2_LW'), p['EWLW']))
    return f


def values(v):
    """5 -> [5.]; [1, 2] -> [1., 2.]; '1, 2' -> [1., 2.]; '0:0.5:2' -> [0., 0.5, ..., 2.]."""
    if isinstance(v, (list, tuple)):
        return [np.nan if x is None else float(x) for x in v]
    if not isinstance(v, str):
        return [np.nan if v is None else float(v)]
    out = []
    for tok in v.replace(';', ',').replace(' ', ',').split(','):
        if ':' in tok:
            a, s, b = (float(x) for x in tok.split(':'))
            out += list(np.round(np.arange(a, b + s * 1e-9, s), 12))
        elif tok:
            out.append(float(tok))
    if not out:
        raise ValueError(f'empty input "{v}"')
    return out


def run_batch(inp, progress=None):
    """Run all combinations of an input dict (keys as in flat_defaults/INPUTS).

    Returns dict(inputs, rows, varied, errors, n); rows are (inputs d, result r,
    MFA index i) with one row per state, combination and MFA."""
    f = flat_defaults(inp.get('wood', 'Birch'))
    if inp.get('polymer') in POLYMERS and 'E_poly' not in inp:
        f['E_poly'], f['nu_poly'] = POLYMERS[inp['polymer']]
    f.update({k: v for k, v in inp.items() if not k.startswith('_')})
    lists = {k: values(f[k]) for k in INPUTS}
    mfa, states = values(f['MFA']), [int(s) for s in values(f['states'])]
    varied = [k for k in INPUTS if len(lists[k]) > 1] + (['MFA'] if len(mfa) > 1 else [])
    combos = [(st, dict(zip(INPUTS, (lists[k][j] for k, j in zip(INPUTS, idx)))))
              for st in states for idx in itertools.product(*(range(len(v)) for v in lists.values()))
              if st != 1 or all(idx[list(INPUTS).index(k)] == 0 for k in POLY)]   # native: polymer once
    rows, errors = [], []
    for n, (st, d) in enumerate(combos):
        p = defaults(f['wood'])
        p.update(state=st, polymer='custom', MFA=mfa, cellulose=f['cellulose_material'],
                 swelling=bool(f['swelling']), E_poly=d['E_poly'], nu_poly=d['nu_poly'],
                 density=d['density'], moisture=d['moisture'], CI=d['CI'],
                 vessel=d['vessel_pct'] / 100, ray=d['ray_pct'] / 100, IF=(d['IF_alpha'], d['IF_beta']),
                 swelling_pct=d['swelling_pct'], cell_ar=d['cell_aspect_ratio'],
                 n_super=d['superellipse_n'], nfam=int(d['n_families']),
                 slend_mf=d['microfibril_slenderness'], tol=int(d['tolerance']),
                 mass_fr=[d[k] for k in ('cellulose', 'hemicellulose', 'lignin', 'extractives')],
                 EWLW=[d[k] for k in ('L_EW', 'W2_EW', 'L_LW', 'W2_LW')])
        try:
            r = run(p)
        except Exception as e:                    # keep the batch running, report at the end
            errors.append(f'{STATES[st]} {d}: {e}')
            continue
        rows += [(dict(d, state=st, MFA=r['MFA'][i]), r, i) for i in range(len(r['MFA']))]
        if progress:
            progress(n + 1, len(combos))
    return dict(inputs=f, rows=rows, varied=varied, errors=errors, n=len(combos))


def write_csv(b, fn):
    """One line per row of run_batch: all inputs, outputs, upper triangles of C."""
    import csv
    iu, f = np.triu_indices(6), b['inputs']
    with open(fn, 'w', newline='', encoding='utf8') as fh:
        w = csv.writer(fh)
        w.writerow(['wood', 'polymer', 'cellulose_material', 'swelling', 'state'] + list(INPUTS)
                   + ['MFA'] + list(OUTPUTS) + [f'Cclw_{a + 1}{c + 1}' for a, c in zip(*iu)]
                   + [f'Ccellwall_{a + 1}{c + 1}' for a, c in zip(*iu)])
        for d, r, i in b['rows']:
            w.writerow([f['wood'], f['polymer'], f['cellulose_material'], f['swelling'], STATES[d['state']]]
                       + [d[k] for k in INPUTS] + [d['MFA']] + [r[k][i] for k in OUTPUTS]
                       + list(r['C_clearwood'][i][iu]) + list(r['C_cellwall'][i][iu]))


if __name__ == '__main__':
    # python hom_tw.py input.json [results.csv]
    import sys
    if len(sys.argv) < 2:
        sys.exit('usage: python hom_tw.py input.json [results.csv]')
    with open(sys.argv[1], encoding='utf8') as fh:
        inp = json.load(fh)
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
        os.path.dirname(os.path.abspath(sys.argv[1])),
        inp.get('output', os.path.splitext(os.path.basename(sys.argv[1]))[0] + '_results.csv'))
    b = run_batch(inp, lambda n, N: print(f'\r{n}/{N} runs', end='', flush=True))
    write_csv(b, out)
    print(f"\n{len(b['rows'])} rows -> {out}" + (f", {len(b['errors'])} failed:" if b['errors'] else ''))
    for e in b['errors']:
        print('  ', e)
