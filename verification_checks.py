#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=============================================================================
 VERIFICATION SUITE
 Companion to "Three-Dimensional Structural Analysis with OpenSeesPy"
 Dr. Wahab, Revision 1.0

 Reproduces the nine independent checks of Part IV of the tutorial.
 Run AFTER portico_3d_3niveles.py is in the same directory:

     python verification_checks.py

 Each check prints the computed value, the independently derived reference
 value, and the discrepancy. A model that passes all nine is one you can
 defend.
=============================================================================
"""

import numpy as np
import openseespy.opensees as ops
import portico_3d_3niveles as M

BAR = '-' * 74


def header(n, title):
    print(f'\n{BAR}\n CHECK {n}.  {title}\n{BAR}')


# =============================================================================
# Build once: gravity state
# =============================================================================
mod  = M.build_model()
grav = M.apply_gravity(mod, factored=True)
M.run_gravity()
ops.reactions()

# -----------------------------------------------------------------------------
header(1, 'Global vertical equilibrium')
Rz = [ops.nodeReaction(M.nid(i, 0), 3) for i in M.PLAN]
W_beams = sum(grav['w'][k] for k in grav['w']) * M.Lper
W_cols  = 4 * M.FACT_D * M.A_col * M.gamma_c * M.Hs * M.NS
print(f'  beams          : {W_beams:12.3f} kN')
print(f'  columns (self) : {W_cols:12.3f} kN')
print(f'  applied total  : {grav["W_applied"]:12.3f} kN')
print(f'  sum(Rz)        : {sum(Rz):12.3f} kN')
print(f'  discrepancy    : {abs(sum(Rz)-grav["W_applied"])/grav["W_applied"]*100:.2e} %')

# -----------------------------------------------------------------------------
header(2, 'Symmetry of base reactions')
print(f'  Rz per column  : {[round(r, 4) for r in Rz]}')
print(f'  W/4            : {grav["W_applied"]/4:.4f} kN')
print(f'  max deviation  : {max(abs(r - grav["W_applied"]/4) for r in Rz):.2e} kN')

# -----------------------------------------------------------------------------
header(3, 'Beam internal forces vs closed-form solutions')
fb = ops.eleResponse(20000 + 100 * 1 + 1, 'localForce')
w1 = grav['w'][1]
Mi, Mj, V = abs(fb[4]), abs(fb[10]), abs(fb[2])
Mfix, Msim = w1 * M.Lx**2 / 12, w1 * M.Lx**2 / 8
Mmid = Msim - 0.5 * (Mi + Mj)
print(f'  span {M.Lx:.2f} m, w = {w1:.3f} kN/m')
print(f'  end shear   FE = {V:8.3f} kN     wL/2  = {w1*M.Lx/2:8.3f} kN'
      f'   ({abs(V-w1*M.Lx/2):.2e} kN)')
print(f'  end moment  FE = {Mi:8.3f} kNm    bounds: 0 < M < wL^2/12 = {Mfix:.3f}')
print(f'  fixity ratio   = {Mi/Mfix:.4f}   (1.0 = rigid joints, 0 = pinned)')
print(f'  M_end + M_mid  = {Mi+Mmid:8.3f} kNm    wL^2/8 = {Msim:8.3f} kNm'
      f'   ({abs(Mi+Mmid-Msim):.2e} kNm)')

# -----------------------------------------------------------------------------
w, T = M.run_modal(min(3 * M.NS, 9))
F, V_b, kexp, d, drift, Rx, ecol, ebm = M.run_elf(T[0], 1)

header(4, 'Base-shear equilibrium')
print(f'  applied V      : {V_b:12.3f} kN')
print(f'  sum(Rx)        : {abs(Rx):12.3f} kN')
print(f'  discrepancy    : {abs(abs(Rx)-V_b)/V_b*100:.2e} %')

# -----------------------------------------------------------------------------
header(5, 'Overturning-moment equilibrium')
ops.reactions()
My   = sum(ops.nodeReaction(M.nid(i, 0), 5) for i in M.PLAN)
Rz2  = {i: ops.nodeReaction(M.nid(i, 0), 3) for i in M.PLAN}
Wt   = sum(Rz2.values())
coup = sum(Rz2[i] * M.PLAN[i][0] for i in M.PLAN) - Wt * M.Lx / 2
OTM  = sum(F[k] * k * M.Hs for k in F)
print(f'  applied  sum(F*h)        : {OTM:10.3f} kNm')
print(f'  resisting column moments : {abs(My):10.3f} kNm  ({abs(My)/OTM*100:5.1f} %)')
print(f'  resisting axial couple   : {coup:10.3f} kNm  ({coup/OTM*100:5.1f} %)')
print(f'  residual                 : {OTM-(abs(My)+coup):10.4f} kNm'
      f'  ({(OTM-(abs(My)+coup))/OTM*100:+.3f} %)')
print(f'  [with TRANSF_COL = "{M.TRANSF_COL}"; the residual is the P-Delta')
print( '   second-order moment and vanishes for "Linear"]')

# -----------------------------------------------------------------------------
header(6, 'Storey stiffness bounds (storey 1 only)')
Vst = [sum(F[j] for j in range(k, M.NS + 1)) for k in range(1, M.NS + 1)]
kst = [Vst[i] / (drift[i] * M.Hs) for i in range(M.NS)]
kmin = 4 * 3 * M.Ec * M.Iy_col / M.Hs**3
kmax = 4 * 12 * M.Ec * M.Iy_col / M.Hs**3
c = M.Ec * M.Iy_col / M.Hs
b = M.Ec * M.Iy_bm / M.Lx
kpor = 2 * 12 * c * (c + 6 * b) / (M.Hs**2 * (2 * c + 3 * b))
print(f'  lower bound  4*3EI/h^3   : {kmin:12.0f} kN/m')
print(f'  computed k1  = V1/D1     : {kst[0]:12.0f} kN/m')
print(f'  upper bound  4*12EI/h^3  : {kmax:12.0f} kN/m')
print(f'  isolated one-bay portal  : {kpor:12.0f} kN/m  (two frames)')
print(f'  within bounds            : {kmin < kst[0] < kmax}')
print(f'  storey stiffnesses       : {[round(v) for v in kst]} kN/m')

# -----------------------------------------------------------------------------
header(7, 'Second-order (P-Delta) effects, isolated')
saved = M.TRANSF_COL
res = {}
for tr in ('Linear', 'PDelta'):
    M.TRANSF_COL = tr
    m2 = M.build_model(); M.apply_gravity(m2, True); M.run_gravity()
    w2, T2 = M.run_modal(2)
    F2, Vb2, _, d2, dr2, Rx2, _, _ = M.run_elf(T2[0], 1)
    ops.reactions()
    My2 = abs(sum(ops.nodeReaction(M.nid(i, 0), 5) for i in M.PLAN))
    Rz3 = {i: ops.nodeReaction(M.nid(i, 0), 3) for i in M.PLAN}
    Wt2 = sum(Rz3.values())
    cp2 = sum(Rz3[i] * M.PLAN[i][0] for i in M.PLAN) - Wt2 * M.Lx / 2
    O2  = sum(F2[k] * k * M.Hs for k in F2)
    res[tr] = (T2[0], My2, cp2, O2 - (My2 + cp2))
M.TRANSF_COL = saved
print('  transf      T1 [s]   sum|My|     couple    residual')
for tr, (t1, mm, cc, rr) in res.items():
    print(f'  {tr:<10s} {t1:7.4f} {mm:10.3f} {cc:10.3f} {rr:11.4f}')

# -----------------------------------------------------------------------------
header(8, 'Stability coefficient  theta = P*drift_e / V')
Pfl = []
wcs = M.A_col * M.gamma_c * M.Hs * 4
for k in range(1, M.NS + 1):
    q = (M.qD_roof if k == M.NS else M.qD_floor) * M.Aplan \
        + M.A_bm * M.gamma_c * M.Lper
    Pfl.append(q + wcs * (0.5 if k == M.NS else 1.0))
print('  storey      P [kN]   drift_e        V [kN]     theta   limit 0.10')
for k in range(M.NS):
    P = sum(Pfl[k:])
    th = P * drift[k] / Vst[k]
    print(f'    {k+1:<5d} {P:10.2f}   {drift[k]:.6f}   {Vst[k]:9.2f}   {th:.5f}'
          f'   {"OK" if th < 0.10 else "P-Delta required"}')

# -----------------------------------------------------------------------------
header(9, 'Time history vs response-spectrum prediction')
mod3 = M.build_model(); M.apply_gravity(mod3, factored=False); M.run_gravity()
w3, T3 = M.run_modal(min(3 * M.NS, 9))
phi  = np.array([ops.nodeEigenvector(M.mid(k), 1, 1) for k in range(1, M.NS + 1)])
mv   = np.array([mod3['Wseis'][k] / M.g for k in range(1, M.NS + 1)])
G1   = (phi * mv).sum() / (phi**2 * mv).sum()


def sdof_peak(Tn, z, ag, dt):
    """Peak relative displacement of a linear SDOF, Newmark average acceleration."""
    wn = 2 * np.pi / Tn
    k, c_, m_ = wn**2, 2 * z * wn, 1.0
    gm, bt = 0.5, 0.25
    kh = k + gm * c_ / (bt * dt) + m_ / (bt * dt**2)
    A1 = m_ / (bt * dt) + gm * c_ / bt
    A2 = m_ / (2 * bt) + dt * c_ * (gm / (2 * bt) - 1)
    u = v = a = 0.0
    um = 0.0
    for i in range(1, len(ag)):
        du = (-m_ * (ag[i] - ag[i-1]) + A1 * v + A2 * a) / kh
        dv = gm * du / (bt * dt) - gm * v / bt + dt * a * (1 - gm / (2 * bt))
        da = du / (bt * dt**2) - v / (bt * dt) - a / (2 * bt)
        u += du; v += dv; a += da
        um = max(um, abs(u))
    return um


t, ag, dt = M.load_record()
Sd  = sdof_peak(T3[0], 0.05, ag, dt)
pred = G1 * phi[-1] * Sd
R = M.run_time_history(zeta=0.05, direction=1)
fe = R['umax'][-1]
print(f'  T1                        : {T3[0]:.4f} s')
print(f'  mode 1 shape (roof = 1)   : {np.round(phi/phi[-1], 4)}')
print(f'  Gamma1 * phi_roof         : {G1*phi[-1]:.4f}')
print(f'  Sd(T1, 5%)                : {Sd*1000:.3f} mm  '
      f'(Sa = {Sd*(2*np.pi/T3[0])**2/M.g:.3f} g)')
print(f'  predicted roof peak       : {pred*1000:.3f} mm')
print(f'  time-history roof peak    : {fe*1000:.3f} mm')
print(f'  difference                : {abs(pred-fe)/fe*100:.2f} %')
print(f'\n{BAR}\n All nine checks complete.\n{BAR}')
