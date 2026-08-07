#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=============================================================================
 FIGURE GENERATION
 Companion to "Three-Dimensional Structural Analysis with OpenSeesPy"
 Dr. Wahab, Revision 1.0

 Produces every figure used in the tutorial, directly from the model.
 Run after portico_3d_3niveles.py is in the same directory:

     python make_figures.py

 Figures are written to ./figures/ as PNG at 150 dpi.
=============================================================================
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')                      # headless backend
import matplotlib.pyplot as plt
import openseespy.opensees as ops
import opsvis as opsv
import portico_3d_3niveles as M

OUT = os.environ.get('FIGDIR', 'figures')
os.makedirs(OUT, exist_ok=True)
DPI = 150

plt.rcParams.update({
    'font.size': 9, 'axes.grid': True, 'grid.alpha': 0.3,
    'axes.titlesize': 10, 'axes.titleweight': 'bold',
    'figure.facecolor': 'white', 'savefig.bbox': 'tight',
})

BLUE, RED, GREY = '#1f4e79', '#c00000', '#808080'


def save(name):
    plt.savefig(f'{OUT}/{name}.png', dpi=DPI)
    plt.close('all')
    print(f'  wrote {OUT}/{name}.png')


def shapes(model):
    d = {t: ['rect', [M.ac, M.ac]] for t in model['cols']}
    d.update({t: ['rect', [M.hb, M.bb]] for t in model['beams']})
    return d


print('Generating figures...')

# =============================================================================
# Fig 1 - model geometry (wireframe + extruded)
# =============================================================================
mod = M.build_model()
opsv.plot_model(node_labels=1, element_labels=0, fig_wi_he=(20., 16.))
plt.title('Model geometry and node numbering')
save('fig01_model')

mod = M.build_model()
opsv.plot_extruded_shapes_3d(shapes(mod), fig_wi_he=(20., 16.))
plt.title('Extruded sections: columns 45x45, beams 30x60')
save('fig02_extruded')

# =============================================================================
# Fig 3-4 - gravity: deformed shape and bending moment
# =============================================================================
mod = M.build_model()
grav = M.apply_gravity(mod, factored=True)
M.run_gravity()

opsv.plot_defo(sfac=300., fig_wi_he=(20., 16.), node_supports=False)
plt.title('Deformed shape under 1.2D + 1.6L  (magnified 300x)')
save('fig03_gravity_defo')

plt.figure(figsize=(8.5, 6.5))
opsv.section_force_diagram_3d('My', 0.02, node_supports=False,
                              end_max_values=False)
plt.title('Bending moment $M_y$ under 1.2D + 1.6L')
save('fig04_gravity_My')

# --- companion 2D diagram for one beam, with the hand check overlaid --------
fb = ops.eleResponse(20000 + 100 * 1 + 1, 'localForce')
wb = grav['w'][1]
Mend, Vend = abs(fb[4]), abs(fb[2])
xb = np.linspace(0, M.Lx, 201)
Mx = wb * M.Lx * xb / 2 - wb * xb**2 / 2 - Mend          # sagging positive
Mmid = wb * M.Lx**2 / 8 - Mend

fig, ax = plt.subplots(figsize=(7.4, 3.9))
ax.fill_between(xb, 0, Mx, where=Mx >= 0, color=BLUE, alpha=0.25)
ax.fill_between(xb, 0, Mx, where=Mx < 0, color=RED, alpha=0.25)
ax.plot(xb, Mx, color=BLUE, lw=2.0)
ax.axhline(0, color='k', lw=1.2)
ax.axhline(wb * M.Lx**2 / 12, color=GREY, ls='--', lw=1.2,
           label='$wL^2/12$ = %.2f (fully fixed)' % (wb * M.Lx**2 / 12))
ax.axhline(-wb * M.Lx**2 / 12, color=GREY, ls='--', lw=1.2)
ax.plot([0, M.Lx], [-Mend, -Mend], 'o', color=RED, ms=6)
ax.plot([M.Lx / 2], [Mmid], 'o', color=BLUE, ms=6)
ax.annotate(f'$M_{{end}}$ = {Mend:.2f}', (0, -Mend), xytext=(8, -16),
            textcoords='offset points', color=RED, fontsize=9)
ax.annotate(f'$M_{{mid}}$ = {Mmid:.2f}', (M.Lx / 2, Mmid), xytext=(0, 10),
            textcoords='offset points', ha='center', color=BLUE, fontsize=9)
ax.set(xlabel='Position along beam 1-2, level 1 [m]',
       ylabel='Bending moment [kN$\\cdot$m]',
       title='Beam 1-2, level 1:  $w$ = %.3f kN/m,  $L$ = %.2f m\n'
             '$M_{end} + M_{mid}$ = %.3f = $wL^2/8$ = %.3f  '
             % (wb, M.Lx, Mend + Mmid, wb * M.Lx**2 / 8),
       xlim=(0, M.Lx))
ax.invert_yaxis()
ax.legend(fontsize=8, loc='lower right')
plt.tight_layout()
save('fig04b_beam_bmd')

plt.figure(figsize=(8.5, 6.5))
opsv.section_force_diagram_3d('N', 0.02, node_supports=False, end_max_values=False)
plt.title('Axial force $N$ under 1.2D + 1.6L')
save('fig04c_gravity_N')

plt.figure(figsize=(8.5, 6.5))
opsv.section_force_diagram_3d('Vz', 0.02, node_supports=False, end_max_values=False)
plt.title('Shear force $V_z$ under 1.2D + 1.6L')
save('fig04d_gravity_Vz')

# =============================================================================
# Fig 5 - mode shapes 1 to 3
# =============================================================================
w, T = M.run_modal(min(3 * M.NS, 9))
labels = ['Mode 1 - translation X', 'Mode 2 - translation Y', 'Mode 3 - torsion']
for i in range(3):
    opsv.plot_mode_shape(i + 1, endDispFlag=0, node_supports=False,
                         fig_wi_he=(18., 15.))
    plt.title(f'{labels[i]}\n$T_{{{i+1}}}$ = {T[i]:.4f} s     '
              f'$f_{{{i+1}}}$ = {1/T[i]:.3f} Hz')
    save(f'fig05_mode{i+1}')

# =============================================================================
# Fig 6 - modal periods and participating mass
# =============================================================================
try:
    mp = ops.modalProperties('-return')
    mx = np.array(mp['partiMassRatiosMX'])
    my = np.array(mp['partiMassRatiosMY'])
    mz = np.array(mp['partiMassRatiosRMZ'])
except Exception:
    mx = my = mz = np.zeros(len(T))

n = len(T)
idx = np.arange(1, n + 1)
fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.4))
ax[0].bar(idx, T, color=BLUE, width=0.6)
ax[0].axhline(0.0466 * M.Htot ** 0.9, color=RED, ls='--', lw=1.2,
              label='$T_a$ = 0.0466$h^{0.9}$ = %.3f s' % (0.0466 * M.Htot ** 0.9))
ax[0].set(xlabel='Mode', ylabel='Period $T_n$ [s]', title='Modal periods')
ax[0].set_xticks(idx); ax[0].legend(fontsize=8)

wid = 0.27
ax[1].bar(idx - wid, mx, wid, label='$U_X$', color=BLUE)
ax[1].bar(idx,       my, wid, label='$U_Y$', color=RED)
ax[1].bar(idx + wid, mz, wid, label='$R_Z$', color=GREY)
ax[1].set(xlabel='Mode', ylabel='Participating mass [%]',
          title='Participating mass by mode')
ax[1].set_xticks(idx); ax[1].legend(fontsize=8)
plt.tight_layout()
save('fig06_modal_summary')

# =============================================================================
# Fig 7 - ELF: force, displacement and drift profiles
# =============================================================================
F, V, kexp, d, drift, Rx, ecol, ebm = M.run_elf(T[0], direction=1)
lev = np.arange(0, M.NS + 1)
z = lev * M.Hs
u = np.array(d) * 1000.0
Fi = np.array([0.0] + [F[k] for k in range(1, M.NS + 1)])
Vst = np.array([sum(F[j] for j in range(k, M.NS + 1)) for k in range(1, M.NS + 1)])
dr = np.array(drift)
dra = dr * M.Cd_fac / M.Ie_fac

fig, ax = plt.subplots(1, 4, figsize=(12.5, 3.8), sharey=True)

ax[0].barh(z[1:], Fi[1:], height=0.5, color=BLUE)
for k in range(1, M.NS + 1):
    ax[0].annotate(f'{Fi[k]:.1f}', (Fi[k], z[k]), xytext=(3, 0),
                   textcoords='offset points', va='center', fontsize=8)
ax[0].set(xlabel='Lateral force $F_i$ [kN]', ylabel='Elevation [m]',
          title='ELF distribution ($k$ = %.1f)' % kexp)
ax[0].set_xlim(0, Fi.max() * 1.3)

def stair(vals, zz):
    """Storey-constant quantity as a vertical staircase: x piecewise const in y."""
    xs, ys = [], []
    for i, v in enumerate(vals):
        xs += [v, v]
        ys += [zz[i], zz[i + 1]]
    return np.array(xs), np.array(ys)

xs, ys = stair(Vst, z)
ax[1].plot(xs, ys, color=RED, lw=1.8)
ax[1].fill_betweenx(ys, 0, xs, color=RED, alpha=0.15)
for k in range(M.NS):
    ax[1].annotate(f'{Vst[k]:.1f}', (Vst[k], 0.5*(z[k]+z[k+1])), xytext=(4, 0),
                   textcoords='offset points', va='center', fontsize=8)
ax[1].set(xlabel='Storey shear $V_i$ [kN]', title='Storey shear',
          xlim=(0, Vst.max() * 1.3))

ax[2].plot(u, z, 'o-', color=BLUE, lw=1.8, ms=4)
for k in range(1, M.NS + 1):
    ax[2].annotate(f'{u[k]:.2f}', (u[k], z[k]), xytext=(-4, 7),
                   textcoords='offset points', ha='right', fontsize=8)
ax[2].set(xlabel='Displacement $u_x$ [mm]', title='Lateral displacement',
          xlim=(0, u.max() * 1.35))

xa, ya = stair(dra * 100, z)
xe, ye = stair(dr * 100, z)
ax[3].plot(xa, ya, color=RED, lw=1.8, label='design $C_d\\delta_{xe}/I_e$')
ax[3].fill_betweenx(ya, 0, xa, color=RED, alpha=0.15)
ax[3].plot(xe, ye, color=GREY, lw=1.4, ls='--', label='elastic')
ax[3].axvline(M.DRIFT_LIMIT * 100, color='k', ls=':', lw=1.6,
              label='limit 2.0 %')
for k in range(M.NS):
    ax[3].annotate(f'{dra[k]*100:.2f}', (dra[k]*100, 0.5*(z[k]+z[k+1])),
                   xytext=(4, 0), textcoords='offset points', va='center',
                   fontsize=8)
ax[3].set(xlabel='Storey drift ratio [%]', title='Drift check',
          xlim=(0, M.DRIFT_LIMIT * 100 * 1.15))
ax[3].legend(fontsize=7, loc='upper right')

for a in ax:
    a.set_ylim(-0.2, M.Htot + 0.3)
    a.set_yticks(z)
plt.tight_layout()
save('fig07_elf_profiles')

# =============================================================================
# Fig 8 - deformed shape under ELF
# =============================================================================
opsv.plot_defo(sfac=300., fig_wi_he=(20., 16.), node_supports=False)
plt.title('Deformed shape under equivalent lateral force, X direction\n'
          '(magnified 300x)')
save('fig08_elf_defo')

# =============================================================================
# Fig 9 - overturning resistance breakdown
# =============================================================================
ops.reactions()
My_b = abs(sum(ops.nodeReaction(M.nid(i, 0), 5) for i in M.PLAN))
Rz2 = {i: ops.nodeReaction(M.nid(i, 0), 3) for i in M.PLAN}
Wt = sum(Rz2.values())
couple = sum(Rz2[i] * M.PLAN[i][0] for i in M.PLAN) - Wt * M.Lx / 2
OTM = sum(F[k] * k * M.Hs for k in F)

fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.4))
ax[0].bar(['Applied\n$\\Sigma F_i h_i$'], [OTM], color=GREY, width=0.5)
ax[0].bar(['Resisting'], [couple], color=BLUE, width=0.5, label='axial couple')
ax[0].bar(['Resisting'], [My_b], bottom=[couple], color=RED, width=0.5,
          label='column base moments')
ax[0].set(ylabel='Moment [kN$\\cdot$m]', title='Overturning equilibrium')
ax[0].legend(fontsize=8)
ax[0].annotate(f'{OTM:.1f}', (0, OTM), ha='center', va='bottom', fontsize=8)
ax[0].annotate(f'{couple+My_b:.1f}', (1, couple + My_b), ha='center',
               va='bottom', fontsize=8)

cols = [1, 2, 3, 4]
grav_col = Wt / 4
ax[1].bar([str(c) for c in cols], [Rz2[c] for c in cols], color=BLUE, width=0.55)
ax[1].axhline(grav_col, color=RED, ls='--', lw=1.3,
              label=f'gravity only = {grav_col:.1f} kN')
ax[1].set(xlabel='Column line', ylabel='Axial reaction $R_z$ [kN]',
          title='Axial force redistribution', ylim=(0, max(Rz2.values()) * 1.25))
for c in cols:
    ax[1].annotate(f'{Rz2[c]:.1f}', (str(c), Rz2[c]), ha='center',
                   va='bottom', fontsize=8)
ax[1].legend(fontsize=8)
plt.tight_layout()
save('fig09_overturning')

# =============================================================================
# Fig 10 - the artificial ground motion: a, v, d
# =============================================================================
R = M.run_time_history(zeta=0.05, direction=1)
acc, dt_r, tg = R['ag'], R['dt'], R['time_in']
vg, dg = M.integrate_motion(acc, dt_r)

fig, ax = plt.subplots(3, 1, figsize=(10.5, 6.0), sharex=True)
ax[0].plot(tg, acc / M.g, color=RED, lw=0.6)
ax[0].set(ylabel='$a_g$ [g]')
ax[0].set_title('Artificial ground motion matched to the ASCE 7-22 design spectrum '
                f'($S_{{DS}}$ = {M.SDS:.2f} g, $S_{{D1}}$ = {M.SD1:.2f} g)')
ax[0].annotate(f'PGA = {np.max(np.abs(acc))/M.g:.3f} g', (0.99, 0.90),
               xycoords='axes fraction', ha='right', fontsize=8)
ax[1].plot(tg, vg, color=BLUE, lw=0.7)
ax[1].set(ylabel='$v_g$ [m/s]')
ax[1].annotate(f'PGV = {np.max(np.abs(vg)):.3f} m/s', (0.99, 0.90),
               xycoords='axes fraction', ha='right', fontsize=8)
ax[2].plot(tg, dg, color='#3a7d44', lw=0.9)
ax[2].set(ylabel='$d_g$ [m]', xlabel='Time [s]')
ax[2].annotate(f'PGD = {np.max(np.abs(dg)):.3f} m', (0.99, 0.90),
               xycoords='axes fraction', ha='right', fontsize=8)
for a_ in ax:
    a_.set_xlim(0, tg[-1])
plt.tight_layout()
save('fig10_artificial_motion')

# =============================================================================
# Fig 11 - spectral match and Fourier content
# =============================================================================
Tv = np.geomspace(0.05, 2.0, 120)
Sa_a = M.response_spectrum(acc, dt_r, Tv, 0.05)
Sa_t = M.target_spectrum(Tv)
_, akt, dkt = M.synthetic_record()
Sa_kt = M.response_spectrum(akt, dkt, Tv, 0.05)
ratio = Sa_a / Sa_t

fig, ax = plt.subplots(1, 3, figsize=(13.0, 3.7))
ax[0].plot(Tv, Sa_t / M.g, color='k', lw=2.0, label='target (ASCE 7-22)')
ax[0].plot(Tv, Sa_a / M.g, color=RED, lw=1.4, label='artificial, matched')
ax[0].plot(Tv, Sa_kt / M.g, color=GREY, lw=1.0, ls='--',
           label='unmatched filtered noise')
for i, Ti in enumerate(R['T'][:3]):
    ax[0].axvline(Ti, color=BLUE, ls=':', lw=1.0)
ax[0].annotate('$T_1$', (R['T'][0], 0.06), color=BLUE, fontsize=9, ha='center')
ax[0].set(xscale='log', xlabel='Period $T$ [s]', ylabel='$S_a$ [g]',
          title='Spectral matching ($\\zeta$ = 5 %)')
ax[0].legend(fontsize=7.5)

ax[1].plot(Tv, ratio, color=RED, lw=1.4)
ax[1].axhline(1.0, color='k', lw=1.2)
ax[1].axhspan(0.90, 1.10, color=BLUE, alpha=0.12, label='$\\pm$10 % band')
ax[1].axvline(R['T'][0], color=BLUE, ls=':', lw=1.0)
bnd = (Tv >= 0.1) & (Tv <= 2.0)
ax[1].set(xscale='log', xlabel='Period $T$ [s]', ylabel='$S_a$ / target',
          title='Match ratio  (mean |misfit| = %.2f %% over 0.1-2 s)'
                % (np.mean(np.abs(ratio[bnd] - 1)) * 100),
          ylim=(0.75, 1.25))
ax[1].legend(fontsize=7.5)

freq = np.fft.rfftfreq(len(acc), dt_r)
FA = np.abs(np.fft.rfft(acc)) * dt_r
ax[2].plot(freq, FA, color=BLUE, lw=0.6)
ax[2].axvline(1 / R['T'][0], color=RED, ls='--', lw=1.2,
              label='$f_1$ = %.2f Hz' % (1 / R['T'][0]))
ax[2].set(xscale='log', xlabel='Frequency [Hz]', ylabel='|FFT($a_g$)| [m/s]',
          title='Fourier amplitude spectrum', xlim=(0.1, 25))
ax[2].legend(fontsize=8)
plt.tight_layout()
save('fig11_spectral_match')

# =============================================================================
# Fig 12 - structural response to the artificial wave
# =============================================================================
t_r, U, Vb = R['t'], R['U'], R['Vb']
fig, ax = plt.subplots(2, 1, figsize=(10.5, 5.2), sharex=True)
cmap = [BLUE, '#e08000', RED]
for k in range(M.NS):
    ax[0].plot(t_r, U[:, k] * 1000, lw=0.7, color=cmap[k % 3],
               label=f'level {k+1}')
ipk = int(np.argmax(np.abs(U[:, -1])))
ax[0].plot(t_r[ipk], U[ipk, -1] * 1000, 'ko', ms=5)
ax[0].annotate(f'peak roof {abs(U[ipk,-1])*1000:.1f} mm',
               (t_r[ipk], U[ipk, -1] * 1000), xytext=(14, -18),
               textcoords='offset points', fontsize=8)
ax[0].set(ylabel='$u_x$ [mm]', title='Floor displacement response, X direction')
ax[0].legend(fontsize=8, ncol=3, loc='upper right')

ax[1].plot(t_r, Vb, color='#3a7d44', lw=0.7)
ax[1].axhline(R['Vb_max'], color=RED, ls='--', lw=1.0,
              label=f'peak {R["Vb_max"]:.0f} kN')
ax[1].axhline(-R['Vb_max'], color=RED, ls='--', lw=1.0)
ax[1].axhline(V, color='k', ls=':', lw=1.2, label=f'ELF design V = {V:.0f} kN')
ax[1].axhline(-V, color='k', ls=':', lw=1.2)
ax[1].set(xlabel='Time [s]', ylabel='Base shear [kN]',
          title='Base shear response')
ax[1].legend(fontsize=8, loc='upper right')
for a_ in ax:
    a_.set_xlim(0, t_r[-1])
plt.tight_layout()
save('fig12_response_history')

# =============================================================================
# Fig 13 - response envelopes: time history vs ELF
# =============================================================================
zz = np.arange(0, M.NS + 1) * M.Hs
u_tha = np.r_[0, R['umax']] * 1000
u_elf = np.array(d) * 1000
u_elf_amp = u_elf * M.Cd_fac / M.Ie_fac
dr_tha = R['drift_max']


def stair2(vals, zv):
    xs, ys = [], []
    for i, v in enumerate(vals):
        xs += [v, v]; ys += [zv[i], zv[i + 1]]
    return np.array(xs), np.array(ys)


fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.9))
ax[0].plot(u_tha, zz, 'o-', color=RED, lw=1.8, ms=4, label='time history (elastic)')
ax[0].plot(u_elf_amp, zz, 's--', color=BLUE, lw=1.6, ms=4,
           label='ELF $\\times C_d/I_e$')
ax[0].plot(u_elf, zz, '^:', color=GREY, lw=1.3, ms=4, label='ELF (reduced)')
ax[0].set(xlabel='Peak displacement [mm]', ylabel='Elevation [m]',
          title='Displacement envelope')
ax[0].legend(fontsize=7.5)

xa, ya = stair2(dr_tha * 100, zz)
xb, yb = stair2(np.array(drift) * M.Cd_fac / M.Ie_fac * 100, zz)
ax[1].plot(xa, ya, color=RED, lw=1.8, label='time history')
ax[1].fill_betweenx(ya, 0, xa, color=RED, alpha=0.15)
ax[1].plot(xb, yb, color=BLUE, lw=1.6, ls='--', label='ELF $\\times C_d/I_e$')
ax[1].axvline(M.DRIFT_LIMIT * 100, color='k', ls=':', lw=1.5, label='limit 2 %')
ax[1].set(xlabel='Peak storey drift ratio [%]', title='Drift envelope',
          xlim=(0, max(dr_tha.max() * 100, M.DRIFT_LIMIT * 100) * 1.2))
ax[1].legend(fontsize=7.5)

names = ['ELF design\n$C_sW$', 'ELF $\\times R$', 'Time history\n(elastic)',
         '$S_a(T_1)WM_1^*$']
try:
    mx1 = np.array(ops.modalProperties('-return')['partiMassRatiosMX'])[0] / 100
except Exception:
    mx1 = 0.8457
Vpred = np.interp(R['T'][0], Tv, Sa_a) / M.g * sum(mod['Wseis'].values()) * mx1
vals = [V, V * M.R_fac, R['Vb_max'], Vpred]
ax[2].bar(names, vals, color=[GREY, '#9fb8cd', RED, '#3a7d44'], width=0.6)
for i, v_ in enumerate(vals):
    ax[2].annotate(f'{v_:.0f}', (i, v_), ha='center', va='bottom', fontsize=8)
ax[2].set(ylabel='Base shear [kN]', title='Base shear comparison',
          ylim=(0, max(vals) * 1.25))
ax[2].tick_params(axis='x', labelsize=7.5)

for a_ in ax[:2]:
    a_.set_ylim(-0.2, M.Htot + 0.3); a_.set_yticks(zz)
plt.tight_layout()
save('fig13_response_envelopes')

# =============================================================================
# Fig 14 - single-mode spectrum prediction vs full time history
# =============================================================================
mod3 = M.build_model(); M.apply_gravity(mod3, factored=False); M.run_gravity()
w3, T3 = M.run_modal(3)
phi = np.array([ops.nodeEigenvector(M.mid(k), 1, 1) for k in range(1, M.NS + 1)])
mv = np.array([mod3['Wseis'][k] / M.g for k in range(1, M.NS + 1)])
G1 = (phi * mv).sum() / (phi**2 * mv).sum()
Sd1 = np.interp(T3[0], Tv, Sa_a) / (2 * np.pi / T3[0])**2
pred = G1 * phi[-1] * Sd1
fe = R['umax'][-1]

fig, ax = plt.subplots(1, 2, figsize=(10.0, 3.6))
ax[0].plot(phi / phi[-1], np.arange(1, M.NS + 1) * M.Hs, 'o-', color=BLUE,
           lw=1.8, ms=5, label='mode 1 shape')
ax[0].plot(np.r_[0, R['umax']] / R['umax'][-1], zz, 's--', color=RED, lw=1.5,
           ms=4, label='THA envelope (normalised)')
ax[0].set(xlabel='Normalised amplitude', ylabel='Elevation [m]',
          title='First mode vs response envelope', ylim=(-0.2, M.Htot + 0.3))
ax[0].set_yticks(zz); ax[0].legend(fontsize=8)

ax[1].bar(['Spectrum\n$\\Gamma_1\\phi_{1,roof}S_d$', 'Time history\n(3D FE)'],
          [pred * 1000, fe * 1000], color=[GREY, BLUE], width=0.5)
for i, v_ in enumerate([pred * 1000, fe * 1000]):
    ax[1].annotate(f'{v_:.2f} mm', (i, v_), ha='center', va='bottom', fontsize=9)
ax[1].set(ylabel='Peak roof displacement [mm]',
          title=f'Independent verification: {abs(pred-fe)/fe*100:.2f} % difference',
          ylim=(0, max(pred, fe) * 1000 * 1.25))
plt.tight_layout()
save('fig14_spectrum_check')

# =============================================================================
# Fig 12 - parametric study: T1 vs number of storeys
# =============================================================================
NS_saved = M.NS
ns_list = list(range(1, 11))
T1_list = []
for ns in ns_list:
    M.NS = ns
    m_ = M.build_model()
    M.apply_gravity(m_, factored=False)
    M.run_gravity()
    _, Tn = M.run_modal(min(3 * ns, 6))
    T1_list.append(Tn[0])
M.NS = NS_saved

h = np.array(ns_list) * M.Hs
Ta = 0.0466 * h ** 0.9
fig, ax = plt.subplots(figsize=(6.2, 3.8))
ax.plot(h, T1_list, 'o-', color=BLUE, lw=1.8, ms=5, label='computed $T_1$ (OpenSees)')
ax.plot(h, Ta, 's--', color=RED, lw=1.5, ms=4,
        label='$T_a = 0.0466\\,h^{0.9}$ (ASCE 7-22)')
ax.fill_between(h, Ta, 1.4 * Ta, color=RED, alpha=0.10,
                label='$T_a$ to $C_uT_a$ ($C_u$ = 1.4)')
ax.set(xlabel='Building height $h$ [m]', ylabel='Fundamental period $T_1$ [s]',
       title='Parametric study: fundamental period vs height')
ax.legend(fontsize=8)
plt.tight_layout()
save('fig15_period_vs_height')

print('\nAll figures written to', OUT)
