"""Compare hom_tw.py with MATLAB reference results (matlab/export_benchmark.m).

Prints the relative Frobenius error per RVE level and case, plus timing.
"""
import os
import time

import numpy as np
from scipy.io import loadmat

import wood_stiffness_homogenization.homogenization_tw as h

HERE = os.path.dirname(os.path.abspath(__file__))
B = dict(density=646, moisture=10, CI=0.58, MFA=[10])
CASES = {  # same inputs as the user structs in export_benchmark.m
    'birch_s1': ('Birch', dict(B, state=1)),
    'birch_s2_HEMA': ('Birch', dict(B, state=2, polymer='HEMA')),
    'birch_s3_HEMA': ('Birch', dict(B, state=3, polymer='HEMA')),
    'birch_sweep': ('Birch', dict(MFA=[0, 5, 10, 20, 30])),
    'birch_IF': ('Birch', dict(MFA=[10], IF=(0.1, 0.05))),
    'birch_s2_PMMA': ('Birch', dict(MFA=[10], state=2, polymer='PMMA')),
    'spruce_s1': ('Spruce', {}),
    'spruce_s2_HEMA': ('Spruce', dict(state=2, polymer='HEMA')),
    'spruce_sweep': ('Spruce', dict(MFA=[0, 10, 20, 30])),
    'pine_s1': ('Pine', {}),
}
LEVELS = ['pn', 'cel', 'cw', 'EWuc', 'LWuc', 'vesselwood', 'ringwood', 'ray', 'clearwood']


def rel(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


def test_vs_matlab():
    tol = 1e-6
    M = loadmat(os.path.join(HERE, 'bench_matlab.mat'), squeeze_me=True, struct_as_record=False)
    ok = True

    # 1) raw Hill tensors versus fun_P_ellipsoid_aniso
    hl = M['hill']
    err = max(rel(h.hill_numeric(hl.C0[:, :, k], *hl.shapes[s]), hl.P[:, :, s, k])
              for k in range(3) for s in range(len(hl.shapes)))
    t = time.perf_counter()
    for k in range(3):
        for s in range(len(hl.shapes)):
            h.hill_numeric(hl.C0[:, :, k], *hl.shapes[s])
    tp = (time.perf_counter() - t) / (3 * len(hl.shapes))
    print(f'Hill tensor: max rel. error {err:.1e}, {tp * 1e3:.1f} ms/call '
          f'(MATLAB {hl.time_per_call:.2f} s/call)\n')

    # 2) full model, every RVE level
    print(f"{'case':16s}" + ''.join(f'{l:>11s}' for l in LEVELS) + '   t_py   t_mat')
    for name, (wood, upd) in CASES.items():
        p = h.defaults(wood)
        p.update(upd)
        t = time.perf_counter()
        r = h.run(p)
        tpy = time.perf_counter() - t
        m = getattr(M['bench'], name)
        row = f'{name:16s}'
        for lev in LEVELS:
            key = 'L_' + lev
            if not hasattr(m, key):
                row += f"{'-':>11s}"
                continue
            ref = np.reshape(getattr(m, key), (6, 6, -1), order='F')
            e = max(rel(L[lev], ref[:, :, i]) for i, L in enumerate(r['levels']))
            ok &= e < tol
            row += f'{e:11.1e}'
        print(row + f'{tpy:7.2f}{m.time:8.2f}')
        eng = np.r_[r['E_R'], r['E_T'], r['E_L'], r['G_TL'], r['G_RL'], r['G_RT'],
                    r['nu_TR'], r['nu_LR'], r['nu_TL'], r['nu_LT'], r['M_cw']]
        ref = np.r_[np.reshape(m.E, (3, -1)).ravel(), np.reshape(m.G, (3, -1)).ravel(),
                    np.reshape(m.nu, (4, -1)).ravel(), np.ravel(m.M_cw)]
        e_eng = np.max(np.abs(eng - ref) / np.abs(ref))
        ok &= e_eng < tol
        if e_eng >= tol:
            print(f'   engineering constants max rel. error {e_eng:.1e}')

    # 3) hom_TranspWood_final.m (state 1) against Python birch_s1
    p = h.defaults('Birch')
    p.update(B)
    e = rel(h.run(p)['C_clearwood'][0], M['final'].C_clearwood)
    ok &= e < tol
    print(f'\nhom_TranspWood_final.m vs Python: rel. error {e:.1e} '
          f"(MATLAB {M['final'].time:.1f} s)")
    print('\nALL PASSED' if ok else '\nFAILED (tol %.0e)' % tol)
    
    assert ok, f"Some tests failed (tol {tol:.0e}). See printed output for details."
