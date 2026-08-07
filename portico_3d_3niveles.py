#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=============================================================================
 THREE-DIMENSIONAL STRUCTURAL ANALYSIS WITH OpenSeesPy
 Reinforced concrete space frame - N storeys (default 3)

 Companion script to the tutorial:
   "Three-Dimensional Structural Analysis with OpenSeesPy:
    A Complete Tutorial - From Model Construction to Verified Results"
   Dr. Wahab, Revision 1.0, August 2026

 The script builds a single-bay, NS-storey reinforced concrete space frame and
 runs four analyses in sequence:

   [1] section property report          (Tutorial Sec. 9)
   [2] gravity analysis, 1.2D + 1.6L    (Tutorial Sec. 15)
   [3] modal analysis + participating mass (Tutorial Sec. 16)
   [4] equivalent lateral force + drift (Tutorial Sec. 17)
   [5] linear time-history analysis     (Tutorial Sec. 18)

 Built-in verification (Tutorial Sec. 20):
   - global vertical equilibrium: sum of reactions vs applied load
   - symmetry of base reactions
   - beam shear and moment against closed-form solutions
   - base-shear equilibrium
   - computed T1 against the ASCE 7-22 empirical period

 Modelling features:
   - fully parametric in the number of storeys
   - Saint-Venant torsional constant computed from the section aspect ratio
   - rigid floor diaphragms with polar mass moment of inertia
   - optional ACI 318-19 effective (cracked) stiffness
   - P-Delta geometric transformation for columns
   - gravity actions and seismic mass derived from the same floor pressures
   - full Rayleigh damping pair anchored at modes 1 and 2
   - self-contained synthetic ground motion if no record file is supplied

 UNITS: kN, m, s  =>  mass in kN*s^2/m = tonne;  stress in kPa
=============================================================================
"""

import os
import numpy as np
import openseespy.opensees as ops

# =============================================================================
# 0.  ANALYSIS SWITCHES
# =============================================================================
NS          = 3          # number of storeys  (set to 1 to recover the original)
CRACKED     = False      # True -> ACI 318-19 effective stiffness (0.35Ig / 0.70Ig)
DIAPHRAGM   = True       # True -> rigid floor diaphragms + rotational inertia
TRANSF_COL  = 'PDelta'   # 'Linear' | 'PDelta' | 'Corotational'  (columns)
TRANSF_BEAM = 'Linear'
RECORD_FILE = 'registro.txt'   # if absent, a record is generated (see RECORD_TYPE)
RECORD_TYPE = 'artificial'     # 'artificial' = spectrum-matched | 'kt' = filtered noise
PLOT        = False      # True -> opsvis / matplotlib figures (needs display)

# =============================================================================
# 1.  UNITS
# =============================================================================
m_   = 1.0
kN   = 1.0
sec  = 1.0
cm   = 0.01 * m_
kPa  = kN / m_**2
MPa  = 1.0e3 * kPa
g    = 9.80665 * m_ / sec**2
ton  = kN * sec**2 / m_          # 1 kN*s^2/m == 1 t == 1000 kg

# =============================================================================
# 2.  GEOMETRY
# =============================================================================
Lx  = 6.50 * m_      # bay length, global X
Ly  = 4.80 * m_      # bay length, global Y   (B in the original notebook)
Hs  = 2.70 * m_      # storey height
Htot = NS * Hs

# Plan position of the four column lines  (i = 1..4)
PLAN = {1: (0.0, 0.0), 2: (Lx, 0.0), 3: (Lx, Ly), 4: (0.0, Ly)}

# =============================================================================
# 3.  MATERIAL
# =============================================================================
fc    = 21.0 * MPa                       # f'c = 21 MPa (210 kgf/cm2)
Ec    = 4700.0 * np.sqrt(fc / MPa) * MPa # ACI 318-19, Eq. 19.2.2.1.b -> 2.154e7 kPa
nu    = 0.20                             # ACI 318-19 R19.2.2
Gc    = Ec / (2.0 * (1.0 + nu))
gamma_c = 24.0 * kN / m_**3              # reinforced-concrete unit weight

# =============================================================================
# 4.  SECTIONS
# =============================================================================
ac = 0.45 * m_                  # square column  45 x 45 cm
bb, hb = 0.30 * m_, 0.60 * m_   # beam 30 x 60 cm  (bb = width, hb = depth)


def beta_torsion(ratio):
    """Saint-Venant coefficient beta for a solid rectangle, J = beta*a*b^3,
    with a = long side, b = short side.  (Timoshenko & Goodier, Table.)"""
    r  = np.array([1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0, 1e6])
    bt = np.array([0.1406, 0.1661, 0.1958, 0.2287, 0.2494, 0.2633,
                   0.2808, 0.2913, 0.3123, 0.3333])
    return float(np.interp(ratio, r, bt))


def rect_props(b, h):
    """Gross properties of a rectangle of width b and depth h.
    Iy is taken about the axis parallel to b (bending in the depth direction)."""
    A  = b * h
    Iy = b * h**3 / 12.0
    Iz = h * b**3 / 12.0
    a_, b_ = max(b, h), min(b, h)
    J  = beta_torsion(a_ / b_) * a_ * b_**3
    return A, Iy, Iz, J


A_col, Iy_col, Iz_col, J_col = rect_props(ac, ac)
A_bm,  Iy_bm,  Iz_bm,  J_bm  = rect_props(bb, hb)

# Effective (cracked) stiffness - ACI 318-19 Table 6.6.3.1.1(a)
if CRACKED:
    fI_col, fI_bm, fJ = 0.70, 0.35, 0.20   # fJ per ASCE 41-17 recommendation
else:
    fI_col = fI_bm = fJ = 1.00

Iy_col_e, Iz_col_e, J_col_e = fI_col*Iy_col, fI_col*Iz_col, fJ*J_col
Iy_bm_e,  Iz_bm_e,  J_bm_e  = fI_bm*Iy_bm,  fI_bm*Iz_bm,  fJ*J_bm

# =============================================================================
# 5.  LOADS AND MASS SOURCE
# =============================================================================
qD_floor = 6.00 * kN / m_**2    # slab + finishes + partitions (typical floor)
qL_floor = 2.00 * kN / m_**2    # live load, residential/office
qD_roof  = 5.00 * kN / m_**2
qL_roof  = 1.00 * kN / m_**2

FACT_D, FACT_L = 1.2, 1.6       # ACI 318-19 / ASCE 7-22 combination 1.2D + 1.6L
PSI_L          = 0.25           # live-load fraction in the seismic mass

Aplan = Lx * Ly
Lper  = 2.0 * (Lx + Ly)

# Seismic-design parameters (illustrative, ASCE 7-22 special RC moment frame)
R_fac, Cd_fac, Ie_fac = 8.0, 5.5, 1.0
Cs_base = 0.10                  # V/W, assumed base-shear coefficient
DRIFT_LIMIT = 0.020             # ASCE 7-22 Table 12.12-1, Risk Category II

# =============================================================================
# 6.  MODEL BUILDER
# =============================================================================
def nid(i, k):
    """Node tag of column line i (1..4) at level k (0 = foundation)."""
    return 100 * k + i


def mid(k):
    """Diaphragm master-node tag at level k.

    Tag bands (see tutorial Sec. 7): real nodes occupy 100*k + i, i.e. up to
    100*NS + 4.  Master nodes therefore start at 9000 so that the two bands
    cannot collide for any practical NS.  Using 1000 + k here would clash with
    node 1001 (= level 10, line 1) as soon as NS reaches 10.
    """
    return 9000 + k


def build_model():
    """Create the finite-element model and return a dictionary of tags/weights."""
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)

    # ---------- 6.1 Nodes ----------
    for k in range(NS + 1):
        z = k * Hs
        for i, (x, y) in PLAN.items():
            ops.node(nid(i, k), x, y, z)

    # ---------- 6.2 Boundary conditions ----------
    for i in PLAN:
        ops.fix(nid(i, 0), 1, 1, 1, 1, 1, 1)     # fixed bases

    # ---------- 6.3 Geometric transformations ----------
    ops.geomTransf(TRANSF_COL,  1, 1.0, 0.0, 0.0)   # columns  (local z -> global X)
    ops.geomTransf(TRANSF_BEAM, 2, 0.0, 0.0, 1.0)   # beams    (local z -> global Z)

    # ---------- 6.4 Elements ----------
    cols, beams = [], []
    for k in range(1, NS + 1):
        for i in PLAN:                                   # columns
            tag = 10000 + 100 * k + i          # columns: 10000 + 100*storey + line
            ops.element('elasticBeamColumn', tag, nid(i, k-1), nid(i, k),
                        A_col, Ec, Gc, J_col_e, Iy_col_e, Iz_col_e, 1)
            cols.append(tag)
        for j, (a, b) in enumerate([(1, 2), (2, 3), (3, 4), (4, 1)], start=1):
            tag = 20000 + 100 * k + j         # beams:   20000 + 100*storey + bay
            ops.element('elasticBeamColumn', tag, nid(a, k), nid(b, k),
                        A_bm, Ec, Gc, J_bm_e, Iy_bm_e, Iz_bm_e, 2)
            beams.append(tag)

    # ---------- 6.5 Storey weights ----------
    Wbeam  = A_bm  * gamma_c * Lper                 # beams of one floor
    Wcol_s = A_col * gamma_c * Hs * 4.0             # columns of one storey
    Wfloor, Wseis = {}, {}
    for k in range(1, NS + 1):
        roof = (k == NS)
        qD = qD_roof if roof else qD_floor
        qL = qL_roof if roof else qL_floor
        # half of the storey below + half of the storey above (none at roof)
        Wstruct = Wbeam + 0.5 * Wcol_s * (1.0 if roof else 2.0)
        Wfloor[k] = {'D': qD*Aplan + Wstruct, 'L': qL*Aplan}
        Wseis[k]  = Wfloor[k]['D'] + PSI_L * Wfloor[k]['L']

    # ---------- 6.6 Rigid diaphragms and mass ----------
    if DIAPHRAGM:
        for k in range(1, NS + 1):
            ops.node(mid(k), Lx/2.0, Ly/2.0, k*Hs)
            ops.fix(mid(k), 0, 0, 1, 1, 1, 0)            # only UX, UY, RZ active
            ops.rigidDiaphragm(3, mid(k), *[nid(i, k) for i in PLAN])
            mk  = Wseis[k] / g                            # translational mass
            Jmk = mk * (Lx**2 + Ly**2) / 12.0             # polar mass inertia
            ops.mass(mid(k), mk, mk, 0.0, 0.0, 0.0, Jmk)
            for i in PLAN:                                # small vertical mass
                ops.mass(nid(i, k), 0.0, 0.0, mk/4.0, 0.0, 0.0, 0.0)
    else:
        for k in range(1, NS + 1):
            mk = Wseis[k] / g / 4.0
            for i in PLAN:
                ops.mass(nid(i, k), mk, mk, mk, 0.0, 0.0, 0.0)

    return {'cols': cols, 'beams': beams, 'Wfloor': Wfloor, 'Wseis': Wseis}


# =============================================================================
# 7.  GRAVITY ANALYSIS
# =============================================================================
def apply_gravity(model, factored=True):
    """Apply the floor pressure to the perimeter beams.

    SIMPLIFICATION: the total floor action q*Lx*Ly is distributed among the
    four perimeter beams in proportion to their length, i.e.
        w = q*Lx*Ly / (2*(Lx+Ly))
    This preserves the total load and the centroid (the plan is doubly
    symmetric) but not the exact tributary (trapezoidal/triangular) pattern.
    Beam self-weight is added explicitly.
    """
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    w_sw = A_bm * gamma_c
    wk, Wtot = {}, 0.0
    for k in range(1, NS + 1):
        qD = qD_roof if k == NS else qD_floor
        qL = qL_roof if k == NS else qL_floor
        if factored:
            w = (FACT_D * (qD * Aplan / Lper + w_sw) + FACT_L * qL * Aplan / Lper)
        else:
            w = (qD * Aplan / Lper + w_sw) + qL * Aplan / Lper
        for j in range(1, 5):
            ops.eleLoad('-ele', 20000 + 100 * k + j, '-type', '-beamUniform', 0.0, -w, 0.0)
        # column self-weight, applied as an axial uniform load along local x
        wc = (FACT_D if factored else 1.0) * A_col * gamma_c
        for i in PLAN:
            ops.eleLoad('-ele', 10000 + 100 * k + i, '-type', '-beamUniform',
                        0.0, 0.0, -wc)
        wk[k] = w
        Wtot += w * Lper + 4.0 * wc * Hs
    return {'w': wk, 'W_applied': Wtot}


def run_gravity(nsteps=10):
    ops.wipeAnalysis()
    ops.system('BandGeneral')
    ops.numberer('RCM')
    ops.constraints('Transformation')
    ops.test('NormDispIncr', 1.0e-10, 20)
    ops.algorithm('Linear')
    ops.integrator('LoadControl', 1.0 / nsteps)
    ops.analysis('Static')
    ok = ops.analyze(nsteps)
    ops.loadConst('-time', 0.0)          # freeze gravity, reset pseudo-time
    return ok


# =============================================================================
# 8.  MODAL ANALYSIS
# =============================================================================
def run_modal(nmodes):
    ops.wipeAnalysis()
    ops.system('BandGeneral')
    ops.numberer('RCM')
    ops.constraints('Transformation')
    try:
        lam = ops.eigen('-genBandArpack', nmodes)
    except Exception:
        lam = ops.eigen('-fullGenLapack', nmodes)
    lam = np.array(lam)
    w   = np.sqrt(lam)
    T   = 2.0 * np.pi / w
    return w, T


# =============================================================================
# 9.  EQUIVALENT LATERAL FORCE (ELF)
# =============================================================================
def elf_forces(Wseis, T1, direction=1):
    """ASCE 7-22 Sec. 12.8.3 vertical distribution:  Fi = V * wi*hi^k / sum."""
    if T1 <= 0.5:
        kexp = 1.0
    elif T1 >= 2.5:
        kexp = 2.0
    else:
        kexp = 1.0 + (T1 - 0.5) / 2.0
    W = sum(Wseis.values())
    V = Cs_base * W
    num = {k: Wseis[k] * (k * Hs) ** kexp for k in Wseis}
    den = sum(num.values())
    return {k: V * num[k] / den for k in Wseis}, V, kexp


def run_elf(T1, direction=1):
    """Rebuild -> gravity -> ELF, and return elastic storey drifts."""
    model = build_model()
    apply_gravity(model, factored=True)
    run_gravity()
    F, V, kexp = elf_forces(model['Wseis'], T1)
    ops.timeSeries('Linear', 2)
    ops.pattern('Plain', 2, 2)
    for k, Fk in F.items():
        node = mid(k) if DIAPHRAGM else nid(1, k)
        if DIAPHRAGM:
            p = [0.0] * 6
            p[direction - 1] = Fk
            ops.load(node, *p)
        else:
            for i in PLAN:
                p = [0.0] * 6
                p[direction - 1] = Fk / 4.0
                ops.load(nid(i, k), *p)
    ops.wipeAnalysis()
    ops.system('BandGeneral')
    ops.numberer('RCM')
    ops.constraints('Transformation')
    ops.test('NormDispIncr', 1.0e-10, 20)
    ops.algorithm('Linear')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')
    ops.analyze(1)

    ops.reactions()
    Rx = sum(ops.nodeReaction(nid(i, 0), direction) for i in PLAN)

    d = [0.0]
    for k in range(1, NS + 1):
        node = mid(k) if DIAPHRAGM else nid(1, k)
        d.append(ops.nodeDisp(node, direction))
    drift = [(d[k] - d[k-1]) / Hs for k in range(1, NS + 1)]

    # representative design actions
    ecol = ops.eleResponse(10000 + 100 * 1 + 1, 'localForce')   # storey-1 column
    ebm  = ops.eleResponse(20000 + 100 * 1 + 1, 'localForce')   # level-1 beam, X
    return F, V, kexp, d, drift, Rx, ecol, ebm


# =============================================================================
# 10.  GROUND MOTION
# =============================================================================
# Design spectrum parameters (ASCE 7-22 Sec. 11.4.6).  These define the TARGET
# that the artificial accelerogram is matched to.
SDS, SD1, TL = 1.00, 0.60, 8.0          # in units of g


def target_spectrum(T, SDS=SDS, SD1=SD1, TL=TL):
    """ASCE 7-22 design response spectrum, Sa in m/s^2 (input T in s)."""
    T  = np.atleast_1d(np.asarray(T, dtype=float))
    Ts = SD1 / SDS
    T0 = 0.2 * Ts
    Sa = np.where(T < T0,  SDS * (0.4 + 0.6 * T / T0),
         np.where(T <= Ts, SDS,
         np.where(T <= TL, SD1 / T, SD1 * TL / T**2)))
    return Sa * g


def response_spectrum(ag, dt, Tv, zeta=0.05):
    """5%-damped pseudo-acceleration spectrum by the Nigam-Jennings method.

    Exact for a piecewise-linear excitation, and vectorised over the period
    vector, so the whole spectrum costs one loop over time rather than one
    loop per period.  Returns Sa in the same units as ag.
    """
    Tv = np.asarray(Tv, dtype=float)
    wn = 2.0 * np.pi / Tv
    wd = wn * np.sqrt(1.0 - zeta**2)
    e  = np.exp(-zeta * wn * dt)
    s  = np.sin(wd * dt)
    c  = np.cos(wd * dt)
    z  = zeta / np.sqrt(1.0 - zeta**2)
    q  = zeta / (wn * dt)

    A = e * (z * s + c)
    B = e * s / wd
    C = (1.0/wn**2) * (2.0*q + e*(((1.0-2.0*zeta**2)/(wd*dt) - z)*s
                                  - (1.0 + 2.0*q)*c))
    D = (1.0/wn**2) * (1.0 - 2.0*q + e*((2.0*zeta**2-1.0)/(wd*dt)*s + 2.0*q*c))
    Ap = -e * wn * s / np.sqrt(1.0 - zeta**2)
    Bp =  e * (c - z * s)
    Cp = (1.0/wn**2) * (-1.0/dt + e*((wn/np.sqrt(1.0-zeta**2) + z/dt)*s + c/dt))
    Dp = (1.0/(wn**2 * dt)) * (1.0 - e*(z*s + c))

    p  = -np.asarray(ag, dtype=float)
    u  = np.zeros_like(wn)
    v  = np.zeros_like(wn)
    umax = np.zeros_like(wn)
    for i in range(len(p) - 1):
        un = A*u + B*v + C*p[i] + D*p[i+1]
        vn = Ap*u + Bp*v + Cp*p[i] + Dp*p[i+1]
        u, v = un, vn
        umax = np.maximum(umax, np.abs(u))
    return umax * wn**2                       # pseudo-acceleration


def envelope(t, t1=2.0, t2=14.0, decay=0.30):
    """Jennings-Housner-Tsai intensity envelope: parabolic rise, plateau,
    exponential decay."""
    return np.where(t < t1, (t / t1)**2,
           np.where(t < t2, 1.0, np.exp(-decay * (t - t2))))


def baseline_correct(a, dt, order=2):
    """Remove a low-order polynomial trend so that integrated velocity and
    displacement do not drift.  Applied to velocity, then differentiated back."""
    t = np.arange(len(a)) * dt
    v = np.concatenate([[0.0], np.cumsum(0.5*(a[1:] + a[:-1])) * dt])
    cf = np.polyfit(t, v, order)
    v_corr = v - np.polyval(cf, t)
    a_corr = np.gradient(v_corr, dt)
    return a_corr - a_corr.mean()


def integrate_motion(a, dt):
    """Acceleration -> velocity, displacement, with polynomial detrending."""
    t = np.arange(len(a)) * dt
    v = np.concatenate([[0.0], np.cumsum(0.5*(a[1:] + a[:-1])) * dt])
    v = v - np.polyval(np.polyfit(t, v, 2), t)
    d = np.concatenate([[0.0], np.cumsum(0.5*(v[1:] + v[:-1])) * dt])
    d = d - np.polyval(np.polyfit(t, d, 2), t)
    return v, d


def highpass_taper(freq, f1=0.20, f2=0.40):
    """Cosine taper that removes Fourier content below f1 Hz and passes it
    fully above f2 Hz.  Without this, spectral matching injects spurious
    very-long-period energy and the integrated displacement drifts to
    physically impossible values."""
    w = np.ones_like(freq)
    w[freq <= f1] = 0.0
    m = (freq > f1) & (freq < f2)
    w[m] = 0.5 * (1.0 - np.cos(np.pi * (freq[m] - f1) / (f2 - f1)))
    return w


def artificial_record(dt=0.01, tdur=30.0, seed=2026, niter=30,
                      Tmin=0.05, Tmax=2.0, nT=120, zeta=0.05, verbose=False):
    """Generate a SPECTRUM-COMPATIBLE artificial accelerogram.

    Method (standard practice, e.g. ASCE 7-22 Ch. 16 / Gasparini & Vanmarcke):
      1. start from Gaussian white noise shaped by an intensity envelope;
      2. compute its 5%-damped response spectrum;
      3. scale the Fourier amplitude at each frequency by the ratio
         Sa_target / Sa_current, keeping the phase unchanged;
      4. re-apply the envelope and repeat until the misfit is acceptable.

    Phase is never altered, so the non-stationary character of the signal
    survives the matching.  Returns (t, a, dt) with a in m/s^2.
    """
    rng = np.random.default_rng(seed)
    n   = int(round(tdur / dt))
    t   = np.arange(n) * dt
    env = envelope(t)

    a = rng.standard_normal(n) * env
    Tv = np.geomspace(Tmin, Tmax, nT)
    St = target_spectrum(Tv)
    fT = 1.0 / Tv[::-1]                       # ascending frequency
    freq = np.fft.rfftfreq(n, dt)

    for it in range(niter):
        Sc = response_spectrum(a, dt, Tv, zeta)
        ratio = St / np.maximum(Sc, 1e-12)
        # clamp to 1.0 outside the matched band: extrapolating the ratio there
        # would amplify unconstrained long-period content on every iteration
        scale = np.interp(freq, fT, ratio[::-1], left=1.0, right=1.0)
        A = np.fft.rfft(a) * scale
        a = np.fft.irfft(A, n)
        if it < niter - 3:            # keep the signal non-stationary; the
            a *= env                  # last few passes refine the match only
        if verbose:
            mis = np.mean(np.abs(Sc - St) / St) * 100
            print(f'    iteration {it+1:2d}: mean spectral misfit {mis:6.2f} %')

    a = baseline_correct(a, dt)
    return t, a, dt


def synthetic_record(dt=0.01, tdur=30.0, pga=0.35 * 9.80665, seed=11):
    """Kanai-Tajimi filtered noise with an intensity envelope, m/s^2.
    Kept for comparison with the spectrum-matched record above."""
    rng = np.random.default_rng(seed)
    n   = int(tdur / dt)
    t   = np.arange(n) * dt
    wg, zg = 2.0 * np.pi * 2.5, 0.60            # firm-soil filter parameters
    f   = np.fft.rfftfreq(n, dt) * 2.0 * np.pi
    H   = (wg**2 + 2j*zg*wg*f) / (wg**2 - f**2 + 2j*zg*wg*f)
    a   = np.fft.irfft(np.fft.rfft(rng.standard_normal(n)) * np.abs(H), n)
    a  *= envelope(t)
    a  *= pga / np.max(np.abs(a))
    return t, a, dt


def load_record(path=RECORD_FILE, scale=1.0, kind=None):
    """Return (t, a, dt).  Priority: real file -> RECORD_TYPE."""
    if os.path.exists(path):
        data = np.loadtxt(path, skiprows=1)
        t, a = data[:, 0], data[:, 1] * scale
        dt = float(np.round(t[1] - t[0], 6))
        return t, a, dt
    kind = kind or RECORD_TYPE
    return artificial_record() if kind == 'artificial' else synthetic_record()


# =============================================================================
# 11.  LINEAR TIME-HISTORY ANALYSIS
# =============================================================================
def run_time_history(zeta=0.05, direction=1, record=None):
    """Run the THA and return the full response history of every floor.

    Returns a dict with the input motion, the floor displacement histories,
    the base-shear history, and the peak envelopes over the height.
    """
    model = build_model()
    apply_gravity(model, factored=False)     # service gravity for the THA
    run_gravity()
    w, T = run_modal(min(3 * NS, 9))

    # --- full Rayleigh pair anchored at modes 1 and 2 of the excited direction
    w1, w2 = w[0], w[1]
    a0 = 2.0 * zeta * w1 * w2 / (w1 + w2)
    a1 = 2.0 * zeta / (w1 + w2)
    ops.rayleigh(a0, 0.0, 0.0, a1)           # betaKcomm for elastic elements

    t, acc, dt = record if record is not None else load_record()
    ops.timeSeries('Path', 3, '-dt', dt, '-values', *acc.tolist(), '-factor', 1.0)
    ops.pattern('UniformExcitation', 3, direction, '-accel', 3)

    ops.wipeAnalysis()
    ops.system('BandGeneral')
    ops.numberer('RCM')
    ops.constraints('Transformation')
    ops.test('NormDispIncr', 1.0e-10, 20)
    ops.algorithm('Linear')
    ops.integrator('Newmark', 0.50, 0.25)    # average acceleration: stable
    ops.analysis('Transient')

    nodes = [(mid(k) if DIAPHRAGM else nid(1, k)) for k in range(1, NS + 1)]
    nsteps = len(acc)
    tt = np.zeros(nsteps)
    U  = np.zeros((nsteps, NS))              # floor displacements
    Vb = np.zeros(nsteps)                    # base shear
    u0 = np.array([ops.nodeDisp(nd, direction) for nd in nodes])

    for i in range(nsteps):
        if ops.analyze(1, dt) != 0:
            print(f'  !! transient analysis failed at step {i}')
            tt, U, Vb = tt[:i], U[:i], Vb[:i]
            break
        tt[i] = ops.getTime()
        U[i]  = np.array([ops.nodeDisp(nd, direction) for nd in nodes]) - u0
        ops.reactions()
        Vb[i] = -sum(ops.nodeReaction(nid(j, 0), direction) for j in PLAN)

    # --- peak envelopes over the height
    umax = np.max(np.abs(U), axis=0)
    Uall = np.column_stack([np.zeros(len(U)), U])          # add ground level
    drift_hist = np.diff(Uall, axis=1) / Hs
    dmax = np.max(np.abs(drift_hist), axis=0)

    return {'t': tt, 'ag': acc, 'dt': dt, 'time_in': t,
            'U': U, 'Vb': Vb, 'T': T, 'rayleigh': (a0, a1),
            'umax': umax, 'drift_max': dmax,
            'Vb_max': np.max(np.abs(Vb)) if len(Vb) else np.nan,
            'nodes': nodes}


# =============================================================================
# 12.  MAIN
# =============================================================================
if __name__ == '__main__':

    line = '=' * 74
    print(line)
    print(f' 3-D RC SPACE FRAME  |  {NS} storeys  |  {Lx:.2f} x {Ly:.2f} x '
          f'{Htot:.2f} m  |  units: kN, m, s')
    print(line)

    print('\n[1] SECTION PROPERTIES (gross)')
    print(f'  Ec = {Ec:,.0f} kPa   Gc = {Gc:,.0f} kPa   nu = {nu:.2f}')
    print(f'  Column {ac*100:.0f}x{ac*100:.0f}: A={A_col:.4f} m2  '
          f'Iy=Iz={Iy_col:.5e} m4  J={J_col:.5e} m4')
    print(f'  Beam   {bb*100:.0f}x{hb*100:.0f}: A={A_bm:.4f} m2  '
          f'Iy={Iy_bm:.5e} m4  Iz={Iz_bm:.5e} m4  J={J_bm:.5e} m4')
    print(f'  Effective stiffness (ACI 318-19): {"ON" if CRACKED else "OFF"}')

    # ---------------- gravity ----------------
    model = build_model()
    grav = apply_gravity(model, factored=True)
    ok = run_gravity()
    print('\n[2] GRAVITY ANALYSIS (1.2D + 1.6L)')
    print(f'  status = {ok} (0 = converged)')
    for k, wv in grav['w'].items():
        print(f'    level {k}: uniform beam load w = {wv:7.3f} kN/m')

    ops.reactions()
    Rz = sum(ops.nodeReaction(nid(i, 0), 3) for i in PLAN)
    W_applied = grav['W_applied']
    print(f'  Sum of vertical reactions      = {Rz:12.3f} kN')
    print(f'  Total applied vertical load    = {W_applied:12.3f} kN')
    print(f'  Equilibrium error              = {abs(Rz - W_applied)/W_applied*100:8.2e} %')

    # mid-span deflection of a roof beam, node-based check
    dz = min(ops.nodeDisp(nid(i, NS), 3) for i in PLAN)
    print(f'  Max. vertical joint displacement (axial shortening) = {dz*1000:.3f} mm')

    fb   = ops.eleResponse(20000 + 100 * 1 + 1, 'localForce')
    w1   = grav['w'][1]
    Mfix = w1 * Lx**2 / 12.0
    Msim = w1 * Lx**2 / 8.0
    Mi, Mj = abs(fb[4]), abs(fb[10])
    Mmid = Msim - 0.5 * (Mi + Mj)
    print('  Verification, beam 1-2 at level 1 (span 6.50 m, w = '
          f'{w1:.2f} kN/m):')
    print(f'    end moments from FE analysis   |M| = {Mi:.2f} / {Mj:.2f} kNm')
    print(f'    ideal fixed-end value  wL^2/12     = {Mfix:.2f} kNm')
    print(f'    simply-supported value wL^2/8      = {Msim:.2f} kNm')
    print(f'    mid-span moment (equilibrium)      = {Mmid:.2f} kNm')
    print(f'    end shear from FE  V = {abs(fb[2]):.2f} kN   vs  wL/2 = '
          f'{w1*Lx/2:.2f} kN')

    # ---------------- modal ----------------
    nmodes = min(3 * NS, 9)
    w, T = run_modal(nmodes)
    print(f'\n[3] MODAL ANALYSIS ({nmodes} modes)')
    print('   Mode      T [s]      f [Hz]     w [rad/s]')
    for i in range(nmodes):
        print(f'   {i+1:>4d}   {T[i]:9.4f}   {1/T[i]:9.4f}   {w[i]:10.3f}')

    try:
        mp = ops.modalProperties('-return')
        mx = np.array(mp['partiMassRatiosMX'])
        my = np.array(mp['partiMassRatiosMY'])
        mz = np.array(mp['partiMassRatiosRMZ'])
        print('\n   Participating mass ratios [%]')
        print('   Mode      UX        UY       RZ        (values in %)')
        for i in range(nmodes):
            print(f'   {i+1:>4d}  {mx[i]:7.2f}  {my[i]:7.2f}  {mz[i]:7.2f}')
        print(f'   Sum:  {mx.sum():7.2f}  {my.sum():7.2f}  {mz.sum():7.2f}'
              '   (>= 90 % required by ASCE 7-22 Sec. 12.9.1.1)')
    except Exception as e:
        print(f'   [modalProperties unavailable: {e}]')

    # empirical period for comparison (ASCE 7-22 Eq. 12.8-7, RC moment frames)
    Ta = 0.0466 * Htot ** 0.9
    print(f'\n   Empirical period  Ta = 0.0466*h^0.90 = {Ta:.4f} s')
    print(f'   Computed T1 / Ta = {T[0]/Ta:.2f}   '
          f'(Cu*Ta limit with Cu=1.4: {1.4*Ta:.4f} s)')

    # ---------------- ELF ----------------
    F, V, kexp, d, drift, Rx, ecol, ebm = run_elf(T[0], direction=1)
    Wtot = sum(model['Wseis'].values())
    print('\n[4] EQUIVALENT LATERAL FORCE, X DIRECTION')
    print(f'  Seismic weight W = {Wtot:.2f} kN;  Cs = {Cs_base:.3f};  '
          f'V = {V:.2f} kN;  k = {kexp:.2f}')
    print('  Level    h [m]     F [kN]     ux [mm]   drift_e   drift_a=Cd/Ie   limit')
    for k in range(1, NS + 1):
        da = drift[k-1] * Cd_fac / Ie_fac
        flag = 'OK ' if da <= DRIFT_LIMIT else 'NG!'
        print(f'   {k:>3d}   {k*Hs:7.2f}  {F[k]:9.2f}  {d[k]*1000:9.3f}  '
              f'{drift[k-1]:8.5f}   {da:8.5f}  {DRIFT_LIMIT:6.3f}  {flag}')
    print(f'  Base-shear check: sum(Rx) = {-Rx:.3f} kN vs V = {V:.3f} kN  '
          f'(error {abs(abs(Rx)-V)/V*100:.2e} %)')
    print('\n  Local stress resultants under 1.2D + 1.6L + E_x')
    print('  Column 1, storey 1 (base end i / top end j):')
    print(f'    N   = {ecol[0]:9.2f} / {ecol[6]:9.2f} kN')
    print(f'    Vy  = {ecol[1]:9.2f} / {ecol[7]:9.2f} kN     Vz = {ecol[2]:9.2f} /'
          f' {ecol[8]:9.2f} kN')
    print(f'    T   = {ecol[3]:9.2f} / {ecol[9]:9.2f} kNm')
    print(f'    My  = {ecol[4]:9.2f} / {ecol[10]:9.2f} kNm    Mz = {ecol[5]:9.2f} /'
          f' {ecol[11]:9.2f} kNm')
    print('  Beam 1-2, level 1 (local z vertical -> My is the flexural action):')
    print(f'    Vz  = {ebm[2]:9.2f} / {ebm[8]:9.2f} kN')
    print(f'    My  = {ebm[4]:9.2f} / {ebm[10]:9.2f} kNm')

    # ---------------- time history ----------------
    print('\n[5] LINEAR TIME-HISTORY ANALYSIS')
    if os.path.exists(RECORD_FILE):
        print(f'  Record: {RECORD_FILE} (from file)')
    elif RECORD_TYPE == 'artificial':
        print(f'  Record: artificial, matched to the ASCE 7-22 design spectrum '
              f'(SDS = {SDS:.2f} g, SD1 = {SD1:.2f} g)')
    else:
        print('  Record: synthetic Kanai-Tajimi filtered noise')

    R = run_time_history(zeta=0.05, direction=1)
    acc, dt_r = R['ag'], R['dt']
    v_g, d_g = integrate_motion(acc, dt_r)
    print(f'  n = {len(acc)};  dt = {dt_r:.4f} s;  duration = {len(acc)*dt_r:.1f} s')
    print(f'  PGA = {np.max(np.abs(acc))/g:.3f} g   PGV = {np.max(np.abs(v_g)):.3f} m/s'
          f'   PGD = {np.max(np.abs(d_g)):.3f} m')

    # spectral match quality
    Tv = np.geomspace(0.05, 2.0, 120)
    Sa_a = response_spectrum(acc, dt_r, Tv, 0.05)
    Sa_t = target_spectrum(Tv)
    rat  = Sa_a / Sa_t
    band = (Tv >= 0.1) & (Tv <= 2.0)
    print(f'  Spectral match over 0.10-2.00 s: mean |misfit| = '
          f'{np.mean(np.abs(rat[band]-1))*100:.2f} %,  '
          f'ratio range {rat[band].min():.3f} - {rat[band].max():.3f}')
    print(f'  Sa at T1 = {R["T"][0]:.4f} s : record {np.interp(R["T"][0], Tv, Sa_a)/g:.3f} g'
          f'  vs target {target_spectrum(R["T"][0])[0]/g:.3f} g')

    a0, a1 = R['rayleigh']
    print(f'  Rayleigh: alphaM = {a0:.5f} 1/s, betaK = {a1:.6f} s  '
          f'(zeta = 5 % at modes 1-2)')

    print('\n  Peak response envelopes (time history vs ELF)')
    print('  Level    u_max [mm]   drift_max    ELF drift    THA/ELF')
    for k in range(NS):
        print(f'   {k+1:>3d}   {R["umax"][k]*1000:10.3f}   {R["drift_max"][k]:9.5f}'
              f'    {drift[k]:9.5f}   {R["drift_max"][k]/drift[k]:7.2f}')
    print(f'  Peak roof displacement = {R["umax"][-1]*1000:.3f} mm '
          f'= H/{Htot/R["umax"][-1]:,.0f}')
    print(f'  Peak base shear        = {R["Vb_max"]:.2f} kN '
          f'= {R["Vb_max"]/sum(model["Wseis"].values()):.4f} W')
    print(f'  ELF base shear         = {V:.2f} kN = {Cs_base:.4f} W')
    print(f'  ratio THA / ELF        = {R["Vb_max"]/V:.2f}   '
          f'(compare R = {R_fac:.1f}: the ELF forces are already reduced by R,')
    print( '                                the time history is fully elastic)')

    # independent check: elastic base shear from the first mode alone
    try:
        mp2 = ops.modalProperties('-return')
        mx1 = np.array(mp2['partiMassRatiosMX'])[0] / 100.0
    except Exception:
        mx1 = 0.8457
    Sa1 = np.interp(R['T'][0], Tv, Sa_a)
    Vel = Sa1 / g * sum(model['Wseis'].values()) * mx1
    print(f'  Check: Sa(T1)*W*M1* = {Vel:.2f} kN vs THA peak {R["Vb_max"]:.2f} kN'
          f'  ({abs(Vel-R["Vb_max"])/R["Vb_max"]*100:.1f} %)')
    print('\n' + line)

    # ---------------- optional plots ----------------
    if PLOT:
        import matplotlib.pyplot as plt
        import opsvis as opsv
        build_model(); apply_gravity(model); run_gravity()
        ele_shapes = {t_: ['rect', [ac, ac]] for t_ in model['cols']}
        ele_shapes.update({t_: ['rect', [hb, bb]] for t_ in model['beams']})
        opsv.plot_extruded_shapes_3d(ele_shapes); plt.show()
        opsv.plot_defo(1e3); plt.show()
        opsv.section_force_diagram_3d('My', 0.5e-1); plt.show()
        w, T = run_modal(nmodes)
        for i in range(1, min(nmodes, 3) + 1):
            opsv.plot_mode_shape(i, endDispFlag=0, node_supports=False)
            plt.title(f'Mode {i}  -  T = {T[i-1]:.4f} s'); plt.show()
        if len(hist):
            fig, ax = plt.subplots(2, 1, figsize=(10, 5), sharex=True)
            ax[0].plot(t, acc / g, 'r', lw=0.6); ax[0].set_ylabel('a_g [g]')
            ax[1].plot(hist[:, 0], hist[:, 1] * 1000, 'k', lw=0.6)
            ax[1].set_ylabel('u_roof [mm]'); ax[1].set_xlabel('t [s]')
            for a in ax: a.grid(alpha=.3)
            plt.tight_layout(); plt.show()
