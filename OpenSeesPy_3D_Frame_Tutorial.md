# Three-Dimensional Structural Analysis with OpenSeesPy

## A Complete Tutorial — From Model Construction to Verified Results

**Dr. Wahab**
Revision 1.0 — August 2026

---

## Contents

**Part I — Foundations**
1. [Purpose and scope](#1-purpose-and-scope)
2. [Installation and how to run the model](#2-installation-and-how-to-run-the-model)
3. [How OpenSees thinks: the domain model](#3-how-opensees-thinks-the-domain-model)
4. [Units: the single most common source of error](#4-units-the-single-most-common-source-of-error)

**Part II — Building the model**
5. [The worked example](#5-the-worked-example)
6. [Step 1 — Model space and degrees of freedom](#6-step-1--model-space-and-degrees-of-freedom)
7. [Step 2 — Nodes and a numbering scheme that scales](#7-step-2--nodes-and-a-numbering-scheme-that-scales)
8. [Step 3 — Boundary conditions](#8-step-3--boundary-conditions)
9. [Step 4 — Section properties: A, I and J](#9-step-4--section-properties-a-i-and-j)
10. [Step 5 — Geometric transformations and local axes](#10-step-5--geometric-transformations-and-local-axes)
11. [Step 6 — Elements](#11-step-6--elements)
12. [Step 7 — Diaphragms and mass](#12-step-7--diaphragms-and-mass)

**Part III — Loads and analysis**
13. [Step 8 — Loads, time series and patterns](#13-step-8--loads-time-series-and-patterns)
14. [Step 9 — The analysis stack](#14-step-9--the-analysis-stack)
15. [Step 10 — Gravity analysis](#15-step-10--gravity-analysis)
16. [Step 11 — Modal analysis](#16-step-11--modal-analysis)
17. [Step 12 — Equivalent lateral force and drift](#17-step-12--equivalent-lateral-force-and-drift)
18. [Step 13 — Linear time-history analysis](#18-step-13--linear-time-history-analysis)

**Part IV — Interpretation and verification**
19. [Reading and interpreting the results](#19-reading-and-interpreting-the-results)
20. [Verification: nine independent checks](#20-verification-nine-independent-checks)
21. [Ten pitfalls that silently corrupt a model](#21-ten-pitfalls-that-silently-corrupt-a-model)
22. [Exercises](#22-exercises)

**Appendices**
- [A — Notation](#appendix-a--notation)
- [B — Reference console output](#appendix-b--reference-console-output)
- [C — Formula summary](#appendix-c--formula-summary)

---

# Part I — Foundations

## 1. Purpose and scope

This tutorial teaches three-dimensional structural analysis with **OpenSeesPy**, the Python
interface to the *Open System for Earthquake Engineering Simulation*. It assumes you have
never opened OpenSees before, but that you have taken a first course in structural analysis:
you should be comfortable with the words *stiffness*, *degree of freedom*, *moment of
inertia*, *mode shape* and *fixed-end moment*.

We build one structure and take it all the way through:

- a three-storey, single-bay reinforced concrete space frame;
- gravity analysis under a factored load combination;
- modal analysis, including participating mass;
- an equivalent lateral force (ELF) seismic analysis with storey-drift verification;
- a linear time-history analysis with Rayleigh damping.

Every result in this document was produced by running the accompanying script,
`portico_3d_3niveles.py`, and every result is then checked — by hand where a hand check is
possible, and against an independent calculation where it is not. That last part is the
point of the tutorial. Producing colourful plots from a finite element program is easy.
Knowing that the numbers behind them are right is the skill worth having.

A deliberate emphasis throughout: **OpenSees will not warn you when your model is wrong.**
It will happily analyse a structure whose mass is ten times too large, whose torsional
constant is nonsense, or whose gravity load keeps increasing during the earthquake. It
converges, it prints numbers, and the numbers are wrong. Section 20 and Section 21 exist
because of this.

### What this tutorial does not cover

The model is **linear elastic** throughout. It therefore says nothing about ductility
demand, plastic hinge formation, strength degradation or collapse capacity. Nonlinear
modelling with fibre sections is the natural next step and is outlined in the exercises.
Soil–structure interaction, joint panel-zone flexibility, and slab in-plane flexibility are
also outside the scope.

---

## 2. Installation and how to run the model

OpenSeesPy is distributed as a compiled Python package. Two commands are enough:

```bash
pip install openseespy       # the solver (required)
pip install opsvis           # plotting helpers (optional)
```

`opsvis` depends on `matplotlib` and `numpy`; both are pulled in automatically. If you only
want numbers and not pictures, you can skip `opsvis` entirely — the script guards all
plotting behind a `PLOT` switch that defaults to `False`.

Verify the installation:

```python
import openseespy.opensees as ops
print(ops.version() if hasattr(ops, "version") else "OpenSeesPy imported successfully")
```

Then run the tutorial model:

```bash
python portico_3d_3niveles.py       # the model and all five analyses
python verification_checks.py       # the nine independent checks of Part IV
```

The full expected output is reproduced in [Appendix B](#appendix-b--reference-console-output).
If your numbers match, your installation is sound and you can follow along with confidence.

### The switches at the top of the script

```python
NS          = 3          # number of storeys — set to any integer
CRACKED     = False      # True -> ACI 318-19 effective (cracked) stiffness
DIAPHRAGM   = True       # True -> rigid floor diaphragms + rotational inertia
TRANSF_COL  = 'PDelta'   # 'Linear' | 'PDelta' | 'Corotational'
RECORD_FILE = 'registro.txt'
PLOT        = False
```

The model is fully parametric in `NS`. Setting `NS = 1` recovers a single-storey frame;
setting `NS = 10` builds a ten-storey frame with no other changes. Building parametric
models from the very first line is a habit worth forming immediately: it is the difference
between a script you use once and a script you keep.

---

## 3. How OpenSees thinks: the domain model

Most commercial software presents you with a drawing canvas. OpenSees presents you with an
**object database**, called the *domain*, which you populate by issuing commands. Nothing is
graphical, and nothing is inferred. If you do not tell OpenSees that a node has mass, it has
no mass.

The objects fall into four groups, and you must create them **in this order**, because later
objects refer to earlier ones by integer tag.

### 3.1 Model objects — the structure itself

| Object | Command | What it is |
|---|---|---|
| Node | `ops.node` | A point with coordinates and DOFs |
| Constraint | `ops.fix`, `ops.rigidDiaphragm` | Restraint or kinematic tie |
| Transformation | `ops.geomTransf` | Local↔global axis mapping for a member |
| Element | `ops.element` | A member connecting nodes |
| Mass | `ops.mass` | Inertia attached to a node |

### 3.2 Load objects — what acts on the structure

| Object | Command | What it is |
|---|---|---|
| Time series | `ops.timeSeries` | A scalar function λ(t) |
| Pattern | `ops.pattern` | A container of loads, scaled by one time series |
| Nodal load | `ops.load` | Point force/moment at a node |
| Element load | `ops.eleLoad` | Distributed load along a member |

The distinction between *pattern* and *time series* is the one beginners find least
intuitive. A pattern holds the **spatial distribution** of a load — which nodes, which
members, how much. The time series holds the **scalar multiplier** applied to that whole
distribution. The load actually applied at pseudo-time `t` is

$$\mathbf{P}(t) = \lambda(t)\,\mathbf{P}_{\text{ref}}$$

where **P**<sub>ref</sub> is what you typed into `ops.load` and `ops.eleLoad`, and λ(t) comes
from the time series. With `timeSeries('Linear', 1)`, λ(t) = t. This is exactly why gravity
must be frozen with `loadConst` before a dynamic analysis begins — see §15.3.

### 3.3 Analysis objects — how the equations get solved

These seven objects are recreated for each analysis stage. They are not properties of the
structure; they are properties of the *solution procedure*.

| Object | Command | Role |
|---|---|---|
| Constraint handler | `ops.constraints` | How constraints enter the equations |
| DOF numberer | `ops.numberer` | Equation ordering (bandwidth) |
| System | `ops.system` | Matrix storage and solver |
| Convergence test | `ops.test` | When to stop iterating |
| Algorithm | `ops.algorithm` | How to iterate within a step |
| Integrator | `ops.integrator` | How to advance from step to step |
| Analysis | `ops.analysis` | Assembles the six above |

### 3.4 Recorder and query objects — getting results out

`ops.nodeDisp`, `ops.nodeReaction`, `ops.eleResponse` and `ops.recorder`. We use the query
functions rather than file recorders here, because it keeps everything visible in one place.

### 3.5 The mental model in one sentence

> You *declare* a structure, you *declare* what acts on it, you *declare* how to solve, and
> then you *ask* for answers — in that order, every time.

---

## 4. Units: the single most common source of error

OpenSees is **dimensionless**. It stores numbers, not units. Consistency is entirely your
responsibility, and the price of inconsistency is a model that runs perfectly and answers
the wrong question.

This tutorial uses the **kN – m – s** system:

| Quantity | Unit | Note |
|---|---|---|
| Force | kN | |
| Length | m | |
| Time | s | |
| Stress, modulus | kN/m² = kPa | 1 MPa = 10³ kPa |
| Moment | kN·m | |
| **Mass** | **kN·s²/m** | **= 1 tonne = 1000 kg** |
| Acceleration | m/s² | g = 9.80665 |
| Unit weight | kN/m³ | concrete ≈ 24 |

The derived mass unit is where most people come unstuck. From F = ma, a mass of 1 kN·s²/m
accelerated at 1 m/s² generates 1 kN of inertia force. Since 1 kN ≈ 102 kgf, that mass is
1000 kg — one tonne.

The script declares units explicitly so that they can be read and checked:

```python
m_   = 1.0
kN   = 1.0
sec  = 1.0
kPa  = kN / m_**2
MPa  = 1.0e3 * kPa
g    = 9.80665 * m_ / sec**2
ton  = kN * sec**2 / m_          # 1 kN·s²/m == 1 t == 1000 kg
```

Writing `fc = 21.0 * MPa` instead of `fc = 21000` costs nothing and documents intent.

> **Rule.** Before you run anything, write your unit system at the top of the file as a
> comment, and convert *every* input to it at the point of entry. Never convert halfway
> through.

---

# Part II — Building the model

## 5. The worked example

### 5.1 Geometry

A single-bay, three-storey reinforced concrete space frame — four columns, four beams per
floor, no interior members.

```
        z
        |
        |        7 ——————————— 6          Level 3   z = 8.10 m
        |       /|            /|
        |      / |           / |
        |     8 ——————————— 5  |          (plan repeats at every level)
        |     |  |          |  |
        |     |  3 —————————|— 2          Level 2   z = 5.40 m
        |     | /           | /
        |     4 ——————————— 1              Level 1   z = 2.70 m
        |
        +——————————————————— x            Level 0   z = 0     (fixed)
       /
      y
```

| Parameter | Symbol | Value |
|---|---|---|
| Bay length, global X | L<sub>x</sub> | 6.50 m |
| Bay length, global Y | L<sub>y</sub> | 4.80 m |
| Storey height | H<sub>s</sub> | 2.70 m |
| Number of storeys | NS | 3 |
| Total height | H | 8.10 m |
| Plan area | A<sub>plan</sub> | 31.20 m² |
| Beam perimeter | L<sub>per</sub> | 22.60 m |

Column lines in plan, numbered 1 to 4 counter-clockwise:

| Line | (x, y) |
|---|---|
| 1 | (0, 0) |
| 2 | (6.50, 0) |
| 3 | (6.50, 4.80) |
| 4 | (0, 4.80) |

### 5.2 Materials

| Property | Symbol | Value | Source |
|---|---|---|---|
| Compressive strength | f′<sub>c</sub> | 21 MPa | assumed |
| Elastic modulus | E<sub>c</sub> | 21 538 106 kPa | ACI 318-19 Eq. 19.2.2.1.b, E<sub>c</sub> = 4700√f′<sub>c</sub> (MPa) |
| Poisson's ratio | ν | 0.20 | ACI 318-19 R19.2.2 |
| Shear modulus | G<sub>c</sub> | 8 974 211 kPa | G = E/[2(1+ν)] |
| Unit weight | γ<sub>c</sub> | 24 kN/m³ | reinforced concrete |

Check the modulus by hand: 4700 × √21 = 4700 × 4.5826 = 21 538 MPa = 2.1538 × 10⁷ kPa. ✓

### 5.3 Sections

| Member | Dimensions | b (width) | h (depth) |
|---|---|---|---|
| Columns | 45 × 45 cm | 0.45 m | 0.45 m |
| Beams | 30 × 60 cm | 0.30 m | 0.60 m |

### 5.4 Loads

| Load | Floors 1–2 | Roof |
|---|---|---|
| Dead, q<sub>D</sub> (slab, finishes, partitions) | 6.00 kN/m² | 5.00 kN/m² |
| Live, q<sub>L</sub> | 2.00 kN/m² | 1.00 kN/m² |

Self-weight of beams and columns is computed from γ<sub>c</sub> and applied separately.

- **Gravity combination:** 1.2D + 1.6L (ACI 318-19 / ASCE 7-22)
- **Seismic mass source:** D + 0.25L

> **Design principle.** Gravity actions and seismic mass must come from the *same* physical
> description of the building. If you tell the program the floors weigh 6 kN/m² for load
> purposes and 10 kN/m² for mass purposes, you have analysed two different structures and
> combined their answers.

---

## 6. Step 1 — Model space and degrees of freedom

```python
ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 6)
```

- `ops.wipe()` destroys any previous domain. Always call it first. In an interactive
  notebook, forgetting it is the classic cause of *"element with tag 1 already exists"*.
- `-ndm 3` — three-dimensional space (coordinates x, y, z).
- `-ndf 6` — six degrees of freedom per node: three translations
  (U<sub>x</sub>, U<sub>y</sub>, U<sub>z</sub>) and three rotations
  (R<sub>x</sub>, R<sub>y</sub>, R<sub>z</sub>), in that order.

The DOF order matters constantly. It is the order of the arguments in `ops.fix`, the order in
`ops.load`, the order in `ops.mass`, and the index you pass to `ops.nodeDisp`. Memorise it:

| Index | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| DOF | U<sub>x</sub> | U<sub>y</sub> | U<sub>z</sub> | R<sub>x</sub> | R<sub>y</sub> | R<sub>z</sub> |

Common alternatives: `-ndm 2 -ndf 3` for plane frames, `-ndm 3 -ndf 3` for space trusses.

---

## 7. Step 2 — Nodes and a numbering scheme that scales

A node is a tag and a set of coordinates:

```python
ops.node(nodeTag, x, y, z)
```

For a single-storey frame you can write eight `ops.node` calls by hand. For a thirty-storey
frame you cannot, and you should not want to. Use a **numbering function**:

```python
def nid(i, k):
    """Node tag of column line i (1..4) at level k (0 = foundation)."""
    return 100 * k + i
```

This maps (line, level) → tag as follows:

| Level | Line 1 | Line 2 | Line 3 | Line 4 |
|---|---|---|---|---|
| 0 (base) | 1 | 2 | 3 | 4 |
| 1 | 101 | 102 | 103 | 104 |
| 2 | 201 | 202 | 203 | 204 |
| 3 | 301 | 302 | 303 | 304 |

The tag now *tells you where you are*. Reading `203` you immediately know: level 2, column
line 3. When an error message mentions node 203, you do not have to consult a sketch.

Generating all nodes is then a double loop:

```python
PLAN = {1: (0.0, 0.0), 2: (Lx, 0.0), 3: (Lx, Ly), 4: (0.0, Ly)}

for k in range(NS + 1):          # k = 0 (base) .. NS (roof)
    z = k * Hs
    for i, (x, y) in PLAN.items():
        ops.node(nid(i, k), x, y, z)
```

For `NS = 3` this creates 16 nodes. Change `NS` to 10 and it creates 44. Nothing else in the
script needs editing.

> **Practice.** Reserve numbering *bands* for different object classes. In this model, nodes
> occupy 1–999, diaphragm master nodes occupy 1001–1999, column elements 10000+, beam
> elements 20000+. Bands prevent the tag collisions that appear the moment a model grows.

---

## 8. Step 3 — Boundary conditions

```python
ops.fix(nodeTag, Ux, Uy, Uz, Rx, Ry, Rz)     # 1 = restrained, 0 = free
```

The four base nodes are fully fixed:

```python
for i in PLAN:
    ops.fix(nid(i, 0), 1, 1, 1, 1, 1, 1)
```

Some common alternatives:

| Support | Command |
|---|---|
| Fixed | `ops.fix(n, 1,1,1,1,1,1)` |
| Pinned (3D) | `ops.fix(n, 1,1,1,0,0,0)` |
| Roller on z | `ops.fix(n, 0,0,1,0,0,0)` |

A full fixity is a modelling *assumption*, not a fact. Real foundations rotate. Fixed bases
shorten the computed fundamental period and understate drift. For a preliminary analysis it
is acceptable and it is what we adopt here; for a final design of a drift-critical structure
it should be revisited with rotational springs at the base.

---

## 9. Step 4 — Section properties: A, I and J

The `elasticBeamColumn` element needs five section constants: A, I<sub>y</sub>,
I<sub>z</sub>, J, and (through E and G) the material. Three of these are routine; the fourth,
J, is where most models are quietly wrong.

### 9.1 Area and moments of inertia

For a rectangle of width *b* and depth *h*:

$$A = bh, \qquad I_y = \frac{bh^3}{12}, \qquad I_z = \frac{hb^3}{12}$$

Here I<sub>y</sub> is taken about the axis parallel to *b*, i.e. bending occurs in the
*depth* direction. Which one OpenSees calls "strong" depends entirely on the geometric
transformation — see §10. This is the crux of 3D frame modelling.

### 9.2 The torsional constant J is not the polar moment

For a **circular** section, and only for a circular section,
J = I<sub>p</sub> = I<sub>y</sub> + I<sub>z</sub>. For any other shape, the cross-section
warps out of plane under torsion, and the Saint-Venant torsional constant is smaller.
For a solid rectangle:

$$J = \beta\, a\, b^3$$

where *a* is the **long** side, *b* is the **short** side, and β depends on a/b:

| a/b | 1.0 | 1.2 | 1.5 | 2.0 | 2.5 | 3.0 | 4.0 | 5.0 | 10.0 | ∞ |
|---|---|---|---|---|---|---|---|---|---|---|
| β | 0.1406 | 0.1661 | 0.1958 | 0.2287 | 0.2494 | 0.2633 | 0.2808 | 0.2913 | 0.3123 | 0.3333 |

*(Saint-Venant torsion of a solid rectangle; classical elasticity solution.)*

```python
def beta_torsion(ratio):
    r  = np.array([1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0, 1e6])
    bt = np.array([0.1406, 0.1661, 0.1958, 0.2287, 0.2494, 0.2633,
                   0.2808, 0.2913, 0.3123, 0.3333])
    return float(np.interp(ratio, r, bt))

def rect_props(b, h):
    A  = b * h
    Iy = b * h**3 / 12.0
    Iz = h * b**3 / 12.0
    a_, b_ = max(b, h), min(b, h)
    J  = beta_torsion(a_ / b_) * a_ * b_**3
    return A, Iy, Iz, J
```

### 9.3 Computed values, checked by hand

**Column 0.45 × 0.45 m:**

- A = 0.45² = **0.2025 m²**
- I<sub>y</sub> = I<sub>z</sub> = 0.45⁴/12 = 0.04100625/12 = **3.41719 × 10⁻³ m⁴**
- a/b = 1.0 → β = 0.1406 → J = 0.1406 × 0.45 × 0.45³ = 0.1406 × 0.04100625 =
  **5.76548 × 10⁻³ m⁴**

**Beam 0.30 × 0.60 m:**

- A = 0.30 × 0.60 = **0.1800 m²**
- I<sub>y</sub> = 0.30 × 0.60³/12 = 0.30 × 0.216/12 = **5.40000 × 10⁻³ m⁴**
- I<sub>z</sub> = 0.60 × 0.30³/12 = 0.60 × 0.027/12 = **1.35000 × 10⁻³ m⁴**
- a/b = 0.60/0.30 = 2.0 → β = 0.2287 → J = 0.2287 × 0.60 × 0.027 =
  **3.70494 × 10⁻³ m⁴**

Compare the beam's J with its polar moment:
I<sub>p</sub> = 5.400 + 1.350 = 6.750 × 10⁻³ m⁴. Using I<sub>p</sub> in place of J would
overstate torsional stiffness by 82 %. Using a round number such as 0.05 — which is what
happens when a placeholder is never replaced — overstates it by a factor of **13.5**.

### 9.4 Effective (cracked) stiffness

Reinforced concrete cracks under service loads. For seismic drift verification, codes
require reduced stiffnesses. ACI 318-19 Table 6.6.3.1.1(a) gives:

| Member | I<sub>eff</sub>/I<sub>g</sub> |
|---|---|
| Beams | 0.35 |
| Columns | 0.70 |
| Walls (uncracked / cracked) | 0.70 / 0.35 |

Areas are not reduced. The code is silent on J; ASCE 41-17 suggests a substantially larger
reduction for torsion, and 0.20 J<sub>g</sub> is used here.

```python
if CRACKED:
    fI_col, fI_bm, fJ = 0.70, 0.35, 0.20
else:
    fI_col = fI_bm = fJ = 1.00
```

This tutorial reports gross-section results (`CRACKED = False`) so that the hand checks in
§20 are clean. Expect periods to lengthen by roughly 30–40 % and drifts to increase
accordingly when you switch it on. **For a real drift check, switch it on.**

---

## 10. Step 5 — Geometric transformations and local axes

This is the section to read twice.

Every frame element has its own **local** coordinate system. Element stiffness is formulated
in local axes and then rotated into global axes. `ops.geomTransf` defines that rotation.

```python
ops.geomTransf(transfType, transfTag, vecxzX, vecxzY, vecxzZ)
```

### 10.1 How OpenSees constructs the local triad

1. **Local x** points from node *i* to node *j*. You control it by the node order in the
   element definition.
2. **`vecxz`** is any vector you supply that lies in the local x–z plane and is *not*
   parallel to local x.
3. The remaining axes follow from cross products:

$$\hat{\mathbf{y}}_{\text{loc}} = \hat{\mathbf{v}}_{xz} \times \hat{\mathbf{x}}_{\text{loc}},
\qquad
\hat{\mathbf{z}}_{\text{loc}} = \hat{\mathbf{x}}_{\text{loc}} \times \hat{\mathbf{y}}_{\text{loc}}$$

### 10.2 Worked out for this model

```python
ops.geomTransf(TRANSF_COL,  1, 1.0, 0.0, 0.0)   # columns
ops.geomTransf(TRANSF_BEAM, 2, 0.0, 0.0, 1.0)   # beams
```

**Beams (tag 2, vecxz = (0,0,1)).** Take beam 1→2, which runs along global X:

- x̂<sub>loc</sub> = (1, 0, 0)
- ŷ<sub>loc</sub> = (0,0,1) × (1,0,0) = (0, 1, 0) → along global **Y**
- ẑ<sub>loc</sub> = (1,0,0) × (0,1,0) = (0, 0, 1) → along global **Z (vertical)**

Now repeat for beam 2→3, which runs along global Y:

- x̂<sub>loc</sub> = (0, 1, 0)
- ŷ<sub>loc</sub> = (0,0,1) × (0,1,0) = (−1, 0, 0)
- ẑ<sub>loc</sub> = (0,1,0) × (−1,0,0) = (0, 0, 1) → **again vertical**

That is the whole point of choosing vecxz = (0,0,1) for beams: **local z is vertical for
every horizontal member, regardless of its plan direction.** Gravity therefore acts along
local −z for all of them, and vertical bending is bending about local **y** for all of them.

Consequence: the strong-axis inertia of a beam must be supplied as **I<sub>y</sub>**, and it
equals *b·h³/12* with *h* the depth. For our 30 × 60 beam, I<sub>y</sub> = 5.40 × 10⁻³ m⁴.
Swapping I<sub>y</sub> and I<sub>z</sub> here would lay every beam on its side and make the
floor system four times more flexible — with no error message.

**Columns (tag 1, vecxz = (1,0,0)).** For column 1→101, running along global Z:

- x̂<sub>loc</sub> = (0, 0, 1)
- ŷ<sub>loc</sub> = (1,0,0) × (0,0,1) = (0, −1, 0)
- ẑ<sub>loc</sub> = (0,0,1) × (0,−1,0) = (1, 0, 0) → along global X

For a square column I<sub>y</sub> = I<sub>z</sub>, so the orientation is immaterial here.
For a rectangular column it is critical, and you would need to verify it exactly as above.

> **Why not vecxz = (0,0,1) for columns?** Because local x is already (0,0,1) and the two
> vectors would be parallel. OpenSees cannot form the triad and the analysis fails.

### 10.3 Choosing the transformation type

| Type | Kinematics | Use for |
|---|---|---|
| `Linear` | Small displacement, no geometric stiffness | Gravity, preliminary work, exact hand-check comparisons |
| `PDelta` | Includes P-Δ geometric stiffness | Seismic analysis of buildings — the standard choice |
| `Corotational` | Full large displacement | Very flexible structures, buckling, large rotations |

We use `PDelta` for columns. Section 20.7 demonstrates the numerical consequence, and
Section 20.8 quantifies how small it is for this particular building.

---

## 11. Step 6 — Elements

```python
ops.element('elasticBeamColumn', eleTag, *eleNodes,
            Area, E_mod, G_mod, Jxx, Iy, Iz, transfTag)
```

Note the argument order carefully: **J comes before I<sub>y</sub> and I<sub>z</sub>**. It is
easy to transpose them.

```python
cols, beams = [], []
for k in range(1, NS + 1):
    for i in PLAN:                                    # 4 columns per storey
        tag = 10000 + 100 * k + i
        ops.element('elasticBeamColumn', tag, nid(i, k-1), nid(i, k),
                    A_col, Ec, Gc, J_col_e, Iy_col_e, Iz_col_e, 1)
        cols.append(tag)
    for j, (a, b) in enumerate([(1, 2), (2, 3), (3, 4), (4, 1)], start=1):
        tag = 20000 + 100 * k + j                     # 4 beams per floor
        ops.element('elasticBeamColumn', tag, nid(a, k), nid(b, k),
                    A_bm, Ec, Gc, J_bm_e, Iy_bm_e, Iz_bm_e, 2)
        beams.append(tag)
```

For `NS = 3` this creates 12 columns and 12 beams.

### 11.1 Traverse the perimeter consistently

The beam list `[(1,2), (2,3), (3,4), (4,1)]` walks the perimeter in one continuous
counter-clockwise loop. It would be geometrically equivalent to write `(1,4)` for the last
beam instead of `(4,1)` — the same two nodes are connected, and the analysis gives the same
displacements. But the local x axis of that member would point the opposite way, and so
would its local y and z. Its bending moment diagram would then be plotted with the opposite
sign from its three neighbours.

Nothing fails. The plot simply misleads you. **Keep node ordering consistent around a loop.**

### 11.2 One element per member

Each beam and column is a single element here. For an elastic analysis this is exact for
end forces, because the elastic beam-column element carries the exact solution of the
Euler–Bernoulli equation for a prismatic member. It does mean that:

- deflected shapes plot as cubics between nodes (fine), and
- mid-span moments must be recovered by equilibrium rather than read from a node (see §20.3).

Subdivide members when you need intermediate output points, when you have intermediate
loads or supports, or when you move to distributed-plasticity nonlinear elements.

---

## 12. Step 7 — Diaphragms and mass

### 12.1 Why a rigid diaphragm

A concrete floor slab is very stiff in its own plane. Modelling it as rigid means all points
on a floor share the same three in-plane motions: two translations and one rotation about
the vertical axis. This:

- collapses the in-plane DOFs of a floor from 12 to 3, which speeds up the eigen-solution;
- lets us attach the floor mass and its **rotational inertia** at one point;
- makes torsional modes appear correctly.

```python
for k in range(1, NS + 1):
    ops.node(mid(k), Lx/2.0, Ly/2.0, k*Hs)          # master node at the centre
    ops.fix(mid(k), 0, 0, 1, 1, 1, 0)               # only Ux, Uy, Rz remain active
    ops.rigidDiaphragm(3, mid(k), *[nid(i, k) for i in PLAN])
```

- `rigidDiaphragm(perpDirn, masterNode, *slaveNodes)`: `perpDirn = 3` means the diaphragm
  plane is perpendicular to global Z, i.e. horizontal.
- The master node's out-of-plane DOFs (U<sub>z</sub>, R<sub>x</sub>, R<sub>y</sub>) are
  restrained. They carry no mass and no stiffness, so leaving them free would make the
  stiffness matrix singular.
- The slave nodes keep their vertical and out-of-plane rotational freedom, so the beams still
  bend under gravity. Only in-plane motion is tied.

A rigid diaphragm is a **multi-point constraint**. It therefore requires
`ops.constraints('Transformation')` — `'Plain'` cannot handle it. This is a frequent cause
of confusing failures.

### 12.2 Mass is weight divided by g

This deserves its own heading because it is the error I see most often.

$$m = \frac{W}{g}$$

In the kN–m–s system, a floor weighing 352.92 kN has a mass of
352.92 / 9.80665 = **35.99 kN·s²/m** (i.e. 35.99 tonnes). If you enter 352.92 as the mass,
every period in your model is √9.81 = **3.13 times too long**, and every spectral
acceleration you subsequently read off a design spectrum is wrong.

### 12.3 Rotational mass inertia

For a rectangular floor of mass *m*, rotating about its own centroid:

$$J_m = \frac{m\,(L_x^2 + L_y^2)}{12}$$

For floor 1: J<sub>m</sub> = 35.99 × (6.50² + 4.80²)/12 = 35.99 × 65.29/12 =
**195.8 t·m²**. Omit this and torsional modes vanish from your model entirely.

```python
mk  = Wseis[k] / g                       # translational mass
Jmk = mk * (Lx**2 + Ly**2) / 12.0        # polar mass moment of inertia
ops.mass(mid(k), mk, mk, 0.0, 0.0, 0.0, Jmk)
for i in PLAN:                            # vertical mass stays on the real nodes
    ops.mass(nid(i, k), 0.0, 0.0, mk/4.0, 0.0, 0.0, 0.0)
```

Note the argument order matches the DOF order:
(m<sub>x</sub>, m<sub>y</sub>, m<sub>z</sub>, J<sub>x</sub>, J<sub>y</sub>, J<sub>z</sub>).

### 12.4 Seismic weight, computed and checked

```python
Wbeam  = A_bm  * gamma_c * Lper                  # all beams of one floor
Wcol_s = A_col * gamma_c * Hs * 4.0              # all columns of one storey
Wstruct = Wbeam + 0.5 * Wcol_s * (1.0 if roof else 2.0)
Wfloor[k] = {'D': qD*Aplan + Wstruct, 'L': qL*Aplan}
Wseis[k]  = Wfloor[k]['D'] + 0.25 * Wfloor[k]['L']
```

Each floor takes **half the column height below it and half above** — the roof gets only the
half below.

Hand check for floor 1:

- Beams: 0.18 × 24 × 22.60 = 97.632 kN
- Columns (one full storey height): 0.2025 × 24 × 2.70 × 4 = 52.488 kN
- Slab dead: 6.00 × 31.20 = 187.20 kN
- **D total = 337.32 kN**; L = 2.00 × 31.20 = 62.40 kN
- W<sub>seis,1</sub> = 337.32 + 0.25 × 62.40 = **352.92 kN** ✓

Roof:

- 5.00 × 31.20 + 97.632 + 0.5 × 52.488 = 156.00 + 97.632 + 26.244 = **279.876 kN**
- W<sub>seis,3</sub> = 279.876 + 0.25 × 31.20 = **287.68 kN** ✓

| Level | W<sub>seis</sub> (kN) | m (t) | J<sub>m</sub> (t·m²) |
|---|---|---|---|
| 1 | 352.92 | 35.99 | 195.8 |
| 2 | 352.92 | 35.99 | 195.8 |
| 3 (roof) | 287.68 | 29.34 | 159.6 |
| **Total** | **993.52** | **101.31** | |

---

# Part III — Loads and analysis

## 13. Step 8 — Loads, time series and patterns

### 13.1 Creating the pattern

```python
ops.timeSeries('Linear', 1)      # λ(t) = t
ops.pattern('Plain', 1, 1)       # pattern tag 1, driven by time series 1
```

Everything issued *after* this call, until the next `pattern`, belongs to pattern 1.

### 13.2 Distributed loads on beams

```python
ops.eleLoad('-ele', eleTag, '-type', '-beamUniform', Wy, Wz, Wx)
```

For a 3D element the three components are in **local** axes, in the order
(W<sub>y</sub>, W<sub>z</sub>, W<sub>x</sub>). Since we established in §10.2 that local z is
vertical for every beam, gravity is a negative W<sub>z</sub>:

```python
ops.eleLoad('-ele', tag, '-type', '-beamUniform', 0.0, -w, 0.0)
```

### 13.3 Converting floor pressure to beam line load

There are no interior beams, so the slab spans between the perimeter beams. The exact
tributary distribution for a two-way slab is trapezoidal on the long beams and triangular on
the short ones. For a teaching model we adopt a simpler equivalent that preserves both the
**total load** and its **centroid** (legitimate here because the plan is doubly symmetric):

$$w = \frac{q\,L_x L_y}{2(L_x + L_y)} + \gamma_c A_{\text{beam}}$$

State this assumption explicitly whenever you use it. It is accurate for global response
(reactions, periods, drift) but it redistributes local beam moments between the long and
short beams relative to a true tributary analysis.

Hand check for floor 1, factored:

- Slab: 1.2 × (6.00 × 31.20 / 22.60) = 1.2 × 8.2832 = 9.9398 kN/m
- Beam self-weight: 1.2 × (0.18 × 24) = 1.2 × 4.32 = 5.1840 kN/m
- Live: 1.6 × (2.00 × 31.20 / 22.60) = 1.6 × 2.7611 = 4.4177 kN/m
- **w = 19.542 kN/m** ✓ (script prints 19.542)

Roof: 1.2 × (6.9027 + 4.32) + 1.6 × 1.3805 = 13.4672 + 2.2088 = **15.676 kN/m** ✓

### 13.4 Column self-weight

Applied as a uniform **axial** load along the column's local x axis, which points upward
from node *i* to node *j*, so gravity is negative:

```python
wc = FACT_D * A_col * gamma_c              # 1.2 × 0.2025 × 24 = 5.832 kN/m
ops.eleLoad('-ele', coltag, '-type', '-beamUniform', 0.0, 0.0, -wc)
```

Including this matters: it is what makes the equilibrium check in §20.1 close exactly. If
column weight appears in the mass but not in the load, your reactions will not match your
applied load and you will not know why.

### 13.5 Nodal loads

```python
ops.load(nodeTag, Fx, Fy, Fz, Mx, My, Mz)
```

Used in §17 for the equivalent lateral forces.

---

## 14. Step 9 — The analysis stack

Seven objects, always in the same conceptual order. This is the configuration used for the
static analyses:

```python
ops.wipeAnalysis()                        # discard the previous configuration
ops.system('BandGeneral')
ops.numberer('RCM')
ops.constraints('Transformation')
ops.test('NormDispIncr', 1.0e-10, 20)
ops.algorithm('Linear')
ops.integrator('LoadControl', 1.0/nsteps)
ops.analysis('Static')
ops.analyze(nsteps)
```

### 14.1 What each one does

**`constraints`** — how restraints and multi-point constraints enter the equations.

| Handler | Handles MP constraints? | Note |
|---|---|---|
| `Plain` | **No** | Homogeneous single-point constraints only |
| `Transformation` | Yes | Condenses out slave DOFs; the general choice |
| `Penalty` | Yes | Needs a well-chosen penalty number; can spoil conditioning |
| `Lagrange` | Yes | Exact, but enlarges the system |

Because we use rigid diaphragms, `Transformation` is mandatory throughout.

**`numberer`** — assigns equation numbers. `RCM` (Reverse Cuthill–McKee) minimises
bandwidth. Use it by default.

**`system`** — matrix storage and solver. `BandGeneral` suits frames. `UmfPack` or
`SparseGeneral` are better for large models; `ProfileSPD` exploits symmetry and positive
definiteness.

**`test`** — the convergence criterion, `test(type, tol, maxIter)`. `NormDispIncr` tests the
norm of the displacement increment. For a linear analysis the tolerance is essentially
irrelevant, because convergence is reached in one iteration.

**`algorithm`** — the iteration strategy within a step. **For a linear elastic model, use
`Linear`.** It forms and factorises the stiffness matrix once and solves directly, which is
exact and fastest. `Newton` will also work and give the same answer, but it iterates for no
reason. `KrylovNewton` and `NewtonLineSearch` belong to nonlinear analysis.

**`integrator`** — how load or time advances.

| Integrator | Use |
|---|---|
| `LoadControl(dλ)` | Static, load-controlled |
| `DisplacementControl` | Static pushover |
| `Newmark(γ, β)` | Transient (dynamic) |

**`analysis`** — `'Static'` or `'Transient'`; assembles the six objects above.

### 14.2 Reading the return value

`ops.analyze()` returns **0 for success** and a negative number for failure. This is the
opposite of the usual Python truth convention, and it must be checked:

```python
if ops.analyze(nsteps) != 0:
    raise RuntimeError('gravity analysis failed to converge')
```

---

## 15. Step 10 — Gravity analysis

### 15.1 Why apply gravity in steps

```python
ops.integrator('LoadControl', 1.0/10)
ops.analyze(10)
```

For a linear model, one step of size 1.0 gives the identical answer. Ten steps of 0.1 cost
almost nothing and mean that the same script still works when you later replace the elastic
elements with nonlinear ones, where incremental application is essential. Build the habit
now.

### 15.2 Results

```
Sum of vertical reactions      =     1426.512 kN
Total applied vertical load    =     1426.512 kN
Equilibrium error              = 3.19e-14 %
Max. vertical joint displacement (axial shortening) = -0.413 mm
```

The 3 × 10⁻¹⁴ % is floating-point round-off — the check is exact. See §20.1.

### 15.3 Freezing gravity: `loadConst`

```python
ops.loadConst('-time', 0.0)
```

This is the most consequential single line in the script, and the easiest to omit.

Recall from §3.2 that the applied load is λ(t)·**P**<sub>ref</sub>, with λ(t) = t for a
`Linear` series. After the static analysis, pseudo-time is 1.0 and gravity is fully applied.
If you now start a transient analysis and let the clock run to t = 30 s, gravity becomes
**thirty times** its correct value — smoothly, without warning, and with the structure
crushing itself as the earthquake proceeds.

`loadConst` does two things: it makes the existing pattern's load constant regardless of
time, and `-time 0.0` resets the clock so the ground motion starts at t = 0.

> **Rule.** Every dynamic analysis that follows a gravity analysis must have `loadConst`
> between them. No exceptions.

---

## 16. Step 11 — Modal analysis

### 16.1 The eigenproblem

OpenSees solves

$$\left(\mathbf{K} - \omega_n^2 \mathbf{M}\right)\boldsymbol{\phi}_n = \mathbf{0}$$

returning the eigenvalues λ<sub>n</sub> = ω<sub>n</sub>².

```python
ops.wipeAnalysis()
ops.system('BandGeneral')
ops.numberer('RCM')
ops.constraints('Transformation')
lam = ops.eigen('-genBandArpack', nmodes)
w = np.sqrt(np.array(lam))
T = 2.0 * np.pi / w
```

Two solver options: `-genBandArpack` (default, iterative, efficient, and tolerant of DOFs
with zero mass) and `-fullGenLapack` (direct, robust, but slow and requiring all modes). If
ARPACK fails on a small model, fall back to LAPACK.

Because our columns use the `PDelta` transformation, the eigen-analysis runs on the
gravity-loaded state and includes the geometric stiffness. Gravity *softens* the structure,
so periods are marginally longer than a purely linear analysis would give (§20.7).

### 16.2 Computed periods

| Mode | T (s) | f (Hz) | ω (rad/s) | Character |
|---|---|---|---|---|
| 1 | 0.3162 | 3.163 | 19.873 | Translation X |
| 2 | 0.2939 | 3.402 | 21.377 | Translation Y |
| 3 | 0.1633 | 6.125 | 38.482 | Torsion |
| 4 | 0.0962 | 10.397 | 65.324 | 2nd translation X |
| 5 | 0.0917 | 10.903 | 68.504 | 2nd translation Y |
| 6 | 0.0539 | 18.547 | 116.535 | 3rd translation X |
| 7 | 0.0532 | 18.806 | 118.163 | 3rd translation Y |
| 8 | 0.0516 | 19.398 | 121.883 | 2nd torsion |
| 9 | 0.0316 | 31.623 | 198.690 | 3rd torsion |

**Read the physics from the numbers.** T₁ (X) is longer than T₂ (Y) because the bay in X is
6.50 m against 4.80 m in Y: the longer beam is more flexible in rotation, so it restrains
the column tops less, and the frame sways more easily in X. The building is therefore softer
in its long plan direction — which is the correct and expected result for this geometry.

### 16.3 Participating mass

```python
mp = ops.modalProperties('-return')     # returns percentages, not fractions
```

| Mode | UX (%) | UY (%) | RZ (%) |
|---|---|---|---|
| 1 | 84.57 | — | — |
| 2 | — | 85.36 | — |
| 3 | — | — | 85.89 |
| 4 | 12.05 | — | — |
| 5 | — | 11.62 | — |
| 6 | 3.38 | — | — |
| 7 | — | 3.02 | — |
| 8 | — | — | 11.15 |
| **Σ** | **100.00** | **100.00** | **97.04** |

Three observations worth making explicit:

1. Modes are **cleanly separated** — mode 1 is pure X, mode 2 pure Y, mode 3 pure torsion,
   with no cross-coupling. That is the signature of a doubly symmetric plan with the centre
   of mass at the centre of rigidity. In an asymmetric building these would be coupled, and
   the participating mass table is how you detect it.
2. The first mode captures only 84.6 % of the X mass. A single-mode approximation would
   underestimate base shear by about 15 %. ASCE 7-22 §12.9.1.1 requires ≥ 90 %, so at least
   two modes per direction are needed here.
3. `modalProperties` returns **percentages**. Multiplying by 100 again is an easy mistake
   and produces the absurd "8457 %" that should immediately tell you something is wrong.

### 16.4 Sanity check against the code period

ASCE 7-22 Eq. 12.8-7 for concrete moment-resisting frames (SI):

$$T_a = C_t h_n^{x} = 0.0466\,(8.10)^{0.90} = 0.3062~\text{s}$$

Computed T₁ / T<sub>a</sub> = 0.3162 / 0.3062 = **1.03**, and the upper limit
C<sub>u</sub>T<sub>a</sub> = 1.4 × 0.3062 = 0.4287 s is not exceeded.

This agreement is the single best global indicator that the mass–stiffness balance is
physically sensible. Empirical period formulas are calibrated on measured buildings; if your
computed period differs from T<sub>a</sub> by a factor of two or three, suspect your mass
units before you suspect the code.

---

## 17. Step 12 — Equivalent lateral force and drift

### 17.1 Vertical distribution

ASCE 7-22 §12.8.3 distributes the base shear V over the height as

$$F_i = V\,\frac{w_i h_i^{k}}{\sum_j w_j h_j^{k}}, \qquad
k = \begin{cases}
1.0 & T \le 0.5~\text{s}\\
1.0 + \dfrac{T - 0.5}{2} & 0.5 < T < 2.5~\text{s}\\
2.0 & T \ge 2.5~\text{s}
\end{cases}$$

The exponent *k* reflects higher-mode participation: tall, long-period buildings carry more
of their shear near the top. Here T₁ = 0.316 s < 0.5 s, so **k = 1.0** — a linear
(inverted-triangular) distribution.

With an assumed seismic coefficient C<sub>s</sub> = 0.10 (a placeholder; a real project needs
a site-specific spectrum):

$$V = 0.10 \times 993.52 = 99.35~\text{kN}$$

Hand check of the distribution, with Σ w<sub>i</sub>h<sub>i</sub> =
352.92 × 2.70 + 352.92 × 5.40 + 287.68 × 8.10 = 952.88 + 1905.77 + 2330.21 = 5188.86 kN·m:

| Level | h (m) | w<sub>i</sub>h<sub>i</sub> | F<sub>i</sub> = V·w<sub>i</sub>h<sub>i</sub>/Σ |
|---|---|---|---|
| 1 | 2.70 | 952.88 | 18.24 kN |
| 2 | 5.40 | 1905.77 | 36.49 kN |
| 3 | 8.10 | 2330.21 | 44.62 kN |
| | | **5188.86** | **99.35 kN** ✓ |

The script prints 18.25 / 36.49 / 44.62 kN. ✓

### 17.2 Applying the forces

Because a rigid diaphragm exists, the storey force is applied **at the master node**, which
distributes it to the columns according to their relative stiffness — exactly as a real slab
would.

```python
p = [0.0] * 6
p[direction - 1] = Fk
ops.load(mid(k), *p)
```

### 17.3 Drift results

```
Level    h [m]     F [kN]     ux [mm]   drift_e   drift_a=Cd/Ie   limit
   1      2.70      18.25      1.159   0.00043    0.00236   0.020  OK
   2      5.40      36.49      2.727   0.00058    0.00319   0.020  OK
   3      8.10      44.62      3.754   0.00038    0.00209   0.020  OK
Base-shear check: sum(Rx) = 99.352 kN vs V = 99.352 kN  (error 1.00e-13 %)
```

Storey drift ratio is computed from the difference of consecutive floor displacements:

$$\theta_{\text{drift},i} = \frac{u_i - u_{i-1}}{H_s}$$

Level 2, by hand: (2.727 − 1.159)/2700 = 1.568/2700 = **0.000581** ✓

### 17.4 Elastic versus design drift

The displacements above come from the *reduced* design forces (V = C<sub>s</sub>W, where
C<sub>s</sub> already contains the division by R). The expected inelastic displacement is
recovered by amplifying:

$$\delta_x = \frac{C_d\,\delta_{xe}}{I_e}$$

With C<sub>d</sub> = 5.5 and I<sub>e</sub> = 1.0 (special reinforced concrete moment frame,
Risk Category II), the largest amplified drift ratio is
0.000581 × 5.5 = **0.00319**, against a limit of 0.020 (ASCE 7-22 Table 12.12-1). The frame
passes with a large margin — unsurprising, since 45 × 45 columns over an 8.1 m height are
generously proportioned.

Forgetting the C<sub>d</sub> amplification is one of the most common errors in seismic
practice: it makes every building appear to pass by a factor of five.

### 17.5 Storey stiffness

Dividing storey shear by storey drift displacement:

| Storey | V<sub>i</sub> (kN) | Δ<sub>i</sub> (mm) | k<sub>i</sub> = V/Δ (kN/m) |
|---|---|---|---|
| 1 | 99.35 | 1.159 | 85 693 |
| 2 | 81.11 | 1.568 | 51 726 |
| 3 | 44.62 | 1.027 | 43 454 |

Storey 1 is the stiffest because its columns are fixed at the base. See §20.6 for the bound
check on this number.

---

## 18. Step 13 — Linear time-history analysis

### 18.1 Rayleigh damping

Structural damping is represented as a linear combination of mass and stiffness:

$$\mathbf{C} = a_0\mathbf{M} + a_1\mathbf{K}$$

The damping ratio at circular frequency ω is then

$$\zeta(\omega) = \frac{a_0}{2\omega} + \frac{a_1\omega}{2}$$

Requiring ζ = ζ<sub>target</sub> at two chosen frequencies ω₁ and ω₂ gives

$$a_0 = \frac{2\zeta\,\omega_1\omega_2}{\omega_1 + \omega_2}, \qquad
a_1 = \frac{2\zeta}{\omega_1 + \omega_2}$$

```python
w1, w2 = w[0], w[1]
a0 = 2.0 * zeta * w1 * w2 / (w1 + w2)
a1 = 2.0 * zeta / (w1 + w2)
ops.rayleigh(a0, 0.0, 0.0, a1)     # (alphaM, betaKcurrent, betaKinit, betaKcommitted)
```

**Both terms are needed.** Setting a₀ = 0 and keeping only the stiffness term is a common
shortcut, but it gives

$$\zeta(\omega_1) = \frac{a_1\omega_1}{2} = \zeta\,\frac{\omega_1}{\omega_1+\omega_2}
\approx 0.48\,\zeta$$

— less than half the intended damping in the fundamental mode, which is precisely the mode
that dominates the response. The script prints a₀ = 1.03043 s⁻¹, a₁ = 0.002423 s.

Note also that ζ(ω) rises linearly at high frequency, so the stiffness term heavily damps
high modes. That is usually desirable: it suppresses spurious high-frequency numerical
noise.

### 18.2 The ground motion

```python
ops.timeSeries('Path', 3, '-dt', dt, '-values', *acc, '-factor', 1.0)
ops.pattern('UniformExcitation', 3, direction, '-accel', 3)
```

`UniformExcitation` applies a rigid base acceleration in the given global direction,
generating the effective force −**M ι** a<sub>g</sub>(t).

**The `-factor` argument is a unit conversion, not a scale factor** — or at least, it must
contain one. If your record file is in cm/s² and your model is in metres, the factor must
include 0.01. If it is in units of g, the factor must include 9.80665. Getting this wrong is
undetectable from the output: the response is simply too big or too small by a constant.
Always confirm the peak ground acceleration in g after loading:

```python
print(f"PGA = {np.max(np.abs(acc))/g:.3f} g")
```

If no record file is present, the script generates one — filtered noise through a
Kanai–Tajimi filter with a Jennings-type intensity envelope, scaled to 0.35 g. This keeps
the tutorial self-contained and reproducible.

### 18.3 Newmark integration

```python
ops.integrator('Newmark', 0.50, 0.25)
ops.analysis('Transient')
```

γ = 0.50, β = 0.25 is the **average acceleration** method: second-order accurate,
unconditionally stable, and with no algorithmic damping. It is the right default.
γ = 0.5, β = 1/6 (linear acceleration) is conditionally stable and requires
Δt ≤ 0.551 T<sub>min</sub>. Choosing γ > 0.5 introduces artificial damping.

**Time step.** A useful rule is Δt ≤ T<sub>min</sub>/10 for the highest mode you care to
capture, and also Δt ≤ the record's own sampling interval. With Δt = 0.01 s and T₁ = 0.316 s
we have about 32 steps per fundamental cycle, and mode 3 (T = 0.163 s) still gets 16 — ample.

### 18.4 The stepping loop

```python
top = mid(NS)
d0 = ops.nodeDisp(top, direction)
for i in range(nsteps):
    if ops.analyze(1, dt) != 0:
        print(f'  !! transient analysis failed at step {i}')
        break
    hist[i, 0] = ops.getTime()
    hist[i, 1] = ops.nodeDisp(top, direction) - d0
```

Subtracting `d0` reports displacement **relative to the gravity-deformed state**, which is
what you want. The return value of `analyze` is checked every step, so a divergence is
caught where it happens rather than discovered later in a nonsensical plot.

### 18.5 Result

```
Record: synthetic (Kanai-Tajimi);  n = 3000;  dt = 0.0100 s;  PGA = 0.350 g
Rayleigh: alphaM = 1.03043 1/s, betaK = 0.002423 s  (zeta = 5% at modes 1-2)
Peak roof displacement (X) = 32.157 mm = H/252
Peak roof drift ratio      = 0.00397
```

Section 20.9 verifies this 32.157 mm against an independent response-spectrum calculation.

---

# Part IV — Interpretation and verification

## 19. Reading and interpreting the results

### 19.1 Global versus local force output

Two different functions return element forces, and confusing them produces plausible-looking
nonsense:

| Call | Returns |
|---|---|
| `ops.eleForce(tag)` | Resisting forces in **global** axes |
| `ops.eleResponse(tag, 'localForce')` | Stress resultants in **local** axes |

For design you want the local ones: axial force, two shears, torsion and two bending moments
at each end. The 12 values are ordered

$$[N_i,\;V_{y,i},\;V_{z,i},\;T_i,\;M_{y,i},\;M_{z,i},\;
   N_j,\;V_{y,j},\;V_{z,j},\;T_j,\;M_{y,j},\;M_{z,j}]$$

so index 0 is N at end *i*, index 4 is M<sub>y</sub> at end *i*, index 10 is M<sub>y</sub> at
end *j*, and so on.

### 19.2 Design actions from this model

```
Column 1, storey 1 (base end i / top end j):
  N   =    323.75 /   -308.00 kN
  Vy  =     -7.51 /      7.51 kN     Vz =    -10.52 /     10.52 kN
  T   =     -0.00 /      0.00 kNm
  My  =     32.88 /     -4.48 kNm    Mz =     -6.76 /    -13.52 kNm

Beam 1-2, level 1:
  Vz  =     48.97 /     78.05 kN
  My  =    -13.01 /    107.51 kNm
```

Interpretation:

- **Column 1 is on the windward side** (x = 0) for lateral load in +X, so overturning
  *relieves* its axial force: 323.75 kN against a gravity-only value of 356.63 kN. Column 2
  at x = 6.50 m carries 389.51 kN. The difference, ±32.88 kN, is the overturning couple.
- The sign change in N between ends (+323.75 / −308.00) is the sign convention of the local
  force vector, not a change from compression to tension. Both ends are in compression. The
  difference of 15.75 kN is the column's own factored self-weight
  (5.832 × 2.70 = 15.75 kN ✓ — a small but satisfying internal check).
- **The beam's end moments are strongly unequal** (13.01 vs 107.51 kN·m) because gravity and
  seismic moments add at one end and oppose at the other. This is exactly why seismic codes
  require top *and* bottom continuous reinforcement in beams, and why the reversal must be
  designed for both directions of loading.
- Column torsion is essentially zero, as it must be for a symmetric structure under a load
  through the centre of rigidity.

### 19.3 What the mode shapes mean

- **Mode 1 (0.3162 s)** — the whole frame sways in X, all floors moving in the same
  direction, largest at the roof. The normalised shape is [0.309, 0.729, 1.000] from level 1
  to roof: close to linear, slightly convex, which is characteristic of a shear-dominated
  frame rather than a flexural cantilever.
- **Mode 2 (0.2939 s)** — the same in Y.
- **Mode 3 (0.1633 s)** — floors rotate about the vertical axis, in phase. That
  T₃/T₁ = 0.52 is a useful diagnostic: when the first torsional period *approaches or
  exceeds* the first translational period, the building is torsionally flexible and codes
  impose penalties. At 0.52 this frame is comfortably torsionally stiff.
- **Modes 4–5** — second translational modes, with a sign reversal up the height.

---

## 20. Verification: nine independent checks

A model is not finished when it runs. It is finished when you have convinced yourself it is
right. Here are nine checks, ordered from cheapest to most demanding. Perform at least the
first four on **every** model you ever build.

All nine are automated in the companion script `verification_checks.py`, which prints the
computed value, the independently derived reference value, and the discrepancy for each.
Work through the hand derivations below first, then run the script and confirm that your
arithmetic and the program agree.

### 20.1 Global vertical equilibrium

The sum of vertical reactions must equal the total applied vertical load.

Applied:

- Beams: (19.542 + 19.542 + 15.676) kN/m × 22.60 m = 54.760 × 22.60 = 1237.55 kN
- Columns: 4 × 5.832 kN/m × 2.70 m × 3 storeys = 4 × 15.7464 × 3 = 188.96 kN
- **Total = 1426.51 kN**

Computed: Σ R<sub>z</sub> = **1426.512 kN**, error 3 × 10⁻¹⁴ %. ✓

This single check catches: wrong load signs, loads applied to the wrong elements, missing
self-weight, and unit errors in the load definition. It costs three lines of code.

### 20.2 Symmetry

A doubly symmetric structure under symmetric load must have four identical base reactions:

```
base Rz per column = [356.628, 356.628, 356.628, 356.628]
W/4                =  356.628
```

Exact to all digits. ✓ Any asymmetry in the output of a symmetric model means an asymmetry
in the input — a mistyped coordinate, a load on three beams instead of four.

### 20.3 Beam internal forces against closed-form solutions

Beam 1–2 at level 1: span L = 6.50 m, w = 19.542 kN/m.

**Shear.** Statics is unavoidable — whatever the end fixity, the end shear of a symmetric
uniformly loaded beam is wL/2:

$$V = \frac{19.542 \times 6.50}{2} = 63.51~\text{kN}$$

Computed: **63.51 kN**. ✓

**Moments.** The end moment must lie between the two classical bounds:

$$\frac{wL^2}{12} = 68.80~\text{kN·m (fully fixed ends)} \quad\text{and}\quad
0 \text{ (pinned ends)}$$

Computed end moment: **60.256 kN·m** — a fixity ratio of 60.256/68.802 = **0.876**, i.e.
12.4 % below the fully fixed value. That shortfall is the rotational flexibility supplied by
the columns, and it is physically sensible: stiffer columns would push the ratio toward 1.0,
more slender ones toward 0.

**Equilibrium of the moment diagram.** For any uniformly loaded prismatic beam, whatever the
end conditions:

$$M_{\text{end}} + M_{\text{mid}} = \frac{wL^2}{8}$$

From the model: 60.256 + 42.948 = **103.204 kN·m**, and wL²/8 = 19.542 × 42.25/8 =
**103.204 kN·m** — agreement to the last printed digit. ✓

This is the check that recovers mid-span moment when you have modelled each beam with a
single element.

### 20.4 Base shear equilibrium

```
sum(Rx) = 99.352 kN vs V = 99.352 kN  (error 1.00e-13 %)
```

Trivial to perform, and it catches lateral loads applied to the wrong DOF — the classic
"I meant X and typed Y" error, which produces a perfectly convergent analysis of a
completely different load case.

### 20.5 Overturning moment equilibrium

Applied overturning moment about the base:

$$M_{OT} = \sum F_i h_i = 18.245(2.70) + 36.490(5.40) + 44.616(8.10) = 607.70~\text{kN·m}$$

This must be resisted by two mechanisms — column base moments and the axial couple:

| Mechanism | Value |
|---|---|
| Σ column base moments M<sub>y</sub> | 182.51 kN·m |
| Axial couple: ΔN × L<sub>x</sub> = 2 × 32.88 × 6.50 | 425.19 kN·m |
| **Total resisting** | **607.70 kN·m** ✓ |

With `TRANSF_COL = 'Linear'` this closes to **0.0000 kN·m residual — exact**.

This check is more informative than it looks: it tells you *how* the frame resists
overturning. Here 70 % comes from the axial couple and 30 % from column bending, which is
typical of a low-rise frame with a wide bay.

### 20.6 Storey stiffness bounds

The first-storey lateral stiffness must lie between two classical limits, because the base
is fixed and only the beam restraint at the top is uncertain:

- **Lower bound** — beams infinitely flexible, columns act as four cantilevers:

$$k_{\min} = 4 \times \frac{3EI}{H_s^3} = 4 \times \frac{3(21\,538\,106)(3.41719\times10^{-3})}{2.70^3} = 44\,871~\text{kN/m}$$

- **Upper bound** — beams infinitely rigid, columns fixed–fixed:

$$k_{\max} = 4 \times \frac{12EI}{H_s^3} = 179\,484~\text{kN/m}$$

Computed: k₁ = **85 693 kN/m**, i.e. 48 % of the rigid-beam bound. Comfortably inside the
bracket. ✓

We can go further. For an isolated one-bay portal, slope-deflection gives the exact sway
stiffness. Let c = EI<sub>c</sub>/H<sub>s</sub> and b = EI<sub>b</sub>/L<sub>x</sub>:

$$k_{\text{portal}} = \frac{12c\,(c + 6b)}{H_s^2\,(2c + 3b)}$$

With c = 73 600.2/2.70 = 27 259.3 kN·m and b = 116 305.8/6.50 = 17 893.2 kN·m:

$$k_{\text{portal}} = \frac{12(27\,259.3)(27\,259.3 + 107\,359.2)}{7.29\,(54\,518.6 + 53\,679.6)}
= 55\,828~\text{kN/m per frame} \;\Rightarrow\; 111\,656~\text{kN/m for two frames}$$

The three-storey frame gives 85 693 kN/m, which is lower — as it should be. In an isolated
portal the level-1 beam restrains only the column below it; in a multi-storey frame it must
share its rotational restraint with the column continuing above. Approximately half the
restraint is diverted upward, and the storey is correspondingly softer.

> Note that this bound applies to **storey 1 only**. For upper storeys, V/Δ is not bounded
> this way, because the drift of an upper storey includes chord rotation accumulated from
> below. The printed k₃ = 43 454 kN/m is *below* the "cantilever lower bound" for exactly
> this reason, and that is not an error.

### 20.7 Second-order effects, isolated and quantified

Re-running the overturning check with both column transformations:

| Transformation | T₁ (s) | ΣM<sub>base</sub> | Axial couple | Residual |
|---|---|---|---|---|
| `Linear` | 0.3153 | 182.510 | 425.191 | **0.0000 kN·m** |
| `PDelta` | 0.3162 | 183.585 | 427.408 | **−3.2914 kN·m** (0.54 %) |

The residual under P-Δ is not an error — it is the second-order moment, which a first-order
equilibrium statement cannot contain. Estimating it by hand as Σ P<sub>i</sub>δ<sub>i</sub>
over the three storeys gives approximately 3.5 kN·m, of the same order as the 3.29 kN·m
observed. (The remaining difference reflects the specific geometric-stiffness formulation of
the `PDelta` transformation.)

The period also lengthens from 0.3153 s to 0.3162 s — gravity softens the structure. Both
effects point the same way, and both are small, which brings us to the next check.

### 20.8 Stability coefficient

ASCE 7-22 Eq. 12.8-16 defines

$$\theta = \frac{P_x \Delta I_e}{V_x h_{sx} C_d}$$

which, since Δ already contains the C<sub>d</sub>/I<sub>e</sub> amplification, reduces to the
convenient form θ = P<sub>x</sub> θ<sub>drift,e</sub> / V<sub>x</sub>.

| Storey | P<sub>x</sub> (kN) | θ<sub>drift,e</sub> | V<sub>x</sub> (kN) | θ |
|---|---|---|---|---|
| 1 | 954.52 | 0.000429 | 99.35 | 0.0041 |
| 2 | 617.20 | 0.000581 | 81.11 | 0.0044 |
| 3 | 279.88 | 0.000380 | 44.62 | 0.0024 |

All values are far below the 0.10 threshold at which P-Δ effects must be explicitly
included, and far below the 0.25 stability limit. This *independently confirms* the 0.54 %
residual found in §20.7: a structure with θ ≈ 0.004 should show second-order effects of
roughly half a percent, and it does. Two unrelated calculations agreeing is the strongest
form of verification available.

### 20.9 Dynamic response predicted from the response spectrum

This is the check that closes the loop on the time-history analysis.

The peak modal displacement of a multi-degree-of-freedom system responding predominantly in
its first mode is

$$u_{\text{roof}}^{\max} \approx \Gamma_1\,\phi_{1,\text{roof}}\;S_d(T_1,\zeta)$$

where the modal participation factor is

$$\Gamma_1 = \frac{\boldsymbol{\phi}_1^{T}\mathbf{M}\boldsymbol{\iota}}
{\boldsymbol{\phi}_1^{T}\mathbf{M}\boldsymbol{\phi}_1}$$

From the eigen-solution: mode shape normalised to the roof is [0.3087, 0.7290, 1.0000], and
Γ₁·φ₁,roof = **1.2850**. Computing the 5 %-damped spectral displacement of the *same* record
at T₁ = 0.3160 s by independent SDOF Newmark integration in NumPy gives
S<sub>d</sub> = **25.207 mm** (S<sub>a</sub> = 1.016 g).

$$u_{\text{roof}}^{\max} \approx 1.2850 \times 25.207 = 32.39~\text{mm}$$

The full time-history analysis gave **32.157 mm** — a difference of **0.7 %**.

Two entirely independent computational paths — a 3D finite element transient analysis, and a
single-degree-of-freedom spectrum combined with a modal participation factor — agree to
within one percent. The residual difference is the contribution of modes 2 and above, which
the single-mode estimate omits. This is exactly the level of agreement you should expect and
should insist on.

### 20.10 Summary table

| # | Check | Result |
|---|---|---|
| 1 | Vertical equilibrium | error 3 × 10⁻¹⁴ % |
| 2 | Symmetry of reactions | exact |
| 3 | Beam V, M against closed form | exact (V), 103.204 vs 103.204 (M) |
| 4 | Base shear | error 1 × 10⁻¹³ % |
| 5 | Overturning equilibrium | exact (linear transf.) |
| 6 | Storey stiffness bounds | 44 871 < 85 693 < 179 484 ✓ |
| 7 | P-Δ residual | 0.54 %, sign and magnitude as expected |
| 8 | Stability coefficient | θ ≈ 0.004 ≪ 0.10, consistent with #7 |
| 9 | Time history vs spectrum | 32.16 vs 32.39 mm, 0.7 % |

---

## 21. Ten pitfalls that silently corrupt a model

Each of these produces a model that runs, converges and prints results. None produces an
error message. This is why §20 exists.

**1. Entering weight where mass is required.**
Writing `mass = W * 9.81` rather than `W / 9.81`. Periods come out √9.81 = 3.13 times too
long; the building appears far more flexible than it is; spectral accelerations read off the
descending branch of the design spectrum are far too low.
*Detect it:* compare T₁ with the code empirical period T<sub>a</sub> (§16.4).

**2. Reusing a variable defined for a different member.**
The archetype is defining `J` when setting up beams and then using that same `J` for
columns. Torsional stiffness is silently wrong for every column in the model.
*Detect it:* print every section constant before the analysis and check the magnitudes by
hand (§9.3).

**3. Using the polar moment as the torsional constant.**
J = I<sub>y</sub> + I<sub>z</sub> is true only for circular sections. For the 30 × 60 beam it
overstates torsional stiffness by 82 %.
*Detect it:* J should always be *smaller* than I<sub>y</sub> + I<sub>z</sub> for a solid
non-circular section.

**4. Omitting `loadConst` before a transient analysis.**
Gravity, driven by a `Linear` time series, grows in proportion to elapsed time. Over a 30 s
record it reaches thirty times its correct value.
*Detect it:* plot vertical reaction against time — it should be flat.

**5. Swapping I<sub>y</sub> and I<sub>z</sub> on beams.**
With vecxz = (0,0,1), vertical bending is governed by I<sub>y</sub>. Supplying the weak-axis
value lays every beam on its side and makes the floor four times more flexible.
*Detect it:* mid-span deflections and periods will be far too large; the storey-stiffness
bound check (§20.6) will fail.

**6. Mass and load derived from different assumptions.**
Applying 1.2 kN/m to the beams while assigning mass from a 10 kN/m² floor pressure means
your gravity design and your seismic analysis describe different buildings.
*Detect it:* compute total applied gravity load and total seismic weight and compare them —
they should differ only by the load factors and the live-load participation factor.

**7. Applying lateral load to the wrong DOF.**
`ops.load(n, 0, 10, 0, 0, 0, 0)` puts the force in Y, not X. Everything converges and you
interpret the answer as an X-direction result.
*Detect it:* the base shear check of §20.4, performed in the direction you *intended*.

**8. Stiffness-proportional damping alone.**
`rayleigh(0, 0, 0, a1)` gives ζ ≈ 0.48 ζ<sub>target</sub> in the fundamental mode.
*Detect it:* evaluate ζ(ω₁) = a₀/(2ω₁) + a₁ω₁/2 by hand and confirm it equals your target.

**9. Unknown or unconverted record units.**
A `-factor` chosen without knowing whether the record is in g, cm/s² or m/s².
*Detect it:* always print PGA in g immediately after loading the record (§18.2).

**10. Inconsistent node ordering around a closed loop.**
Defining the last beam of a perimeter as (1,4) rather than (4,1) reverses its local axes and
therefore the sign of its plotted moment diagram.
*Detect it:* under symmetric gravity load, the four beams of a floor must produce
symmetric-looking diagrams.

**Bonus — mislabelled plots.** Plotting `nodeDisp(...)` and labelling the axis
"Acceleration". The analysis is right; the conclusion drawn from it may not be. Label axes
from the variable actually plotted, not from what you expected to plot.

---

## 22. Exercises

**E1 — Parametric study.** Set `NS = 6` and re-run. Plot T₁ against building height for
NS = 1…10 and compare the trend with T<sub>a</sub> = 0.0466 h<sup>0.9</sup>. At what height
does the model start to diverge from the empirical curve, and why?

**E2 — Cracked sections.** Set `CRACKED = True`. Report the change in T₁, in maximum drift,
and in the stability coefficient θ. Does the frame still satisfy the drift limit? Which of
the three quantities is most sensitive?

**E3 — Torsional irregularity.** Move the diaphragm master nodes to
(0.6 L<sub>x</sub>, 0.5 L<sub>y</sub>) to introduce a mass eccentricity. Examine the
participating-mass table. What happens to the clean separation between modes 1, 2 and 3?
Quantify the increase in the corner column's demand.

**E4 — Torsional constant sensitivity.** Multiply J for all members by 10 and re-run. Which
periods change and which do not? Explain the pattern in terms of mode shapes.

**E5 — Damping model.** Repeat the time-history analysis with a₀ = 0 (stiffness-proportional
damping only). Predict the change in peak roof displacement before running it, using
ζ(ω₁) = a₁ω₁/2, then check your prediction.

**E6 — Time step convergence.** Rerun the time history with Δt = 0.02, 0.01, 0.005 and
0.0025 s. Plot peak roof displacement against Δt and identify the point at which the answer
has converged.

**E7 — Verification by hand.** Without running anything, estimate T₁ for the `NS = 1` case
using the portal formula of §20.6 and the floor mass of §12.4. Then run it and compare.

**E8 — Toward nonlinearity.** Replace the `elasticBeamColumn` elements with
`forceBeamColumn` using fibre sections (`Concrete02` and `Steel02`), and perform a
displacement-controlled pushover. Compare the initial stiffness of the pushover curve with
the elastic k₁ = 85 693 kN/m found in §17.5. They should agree — and if they do not, the
fibre section discretisation is the first place to look.

---

## Appendix A — Notation

| Symbol | Meaning | Unit |
|---|---|---|
| A | Cross-sectional area | m² |
| a₀, a₁ | Rayleigh mass- and stiffness-proportional coefficients | s⁻¹, s |
| C<sub>d</sub> | Deflection amplification factor | – |
| C<sub>s</sub> | Seismic response coefficient, V/W | – |
| E<sub>c</sub> | Concrete elastic modulus | kPa |
| F<sub>i</sub> | Lateral force at level *i* | kN |
| G<sub>c</sub> | Shear modulus | kPa |
| g | Gravitational acceleration, 9.80665 | m/s² |
| H<sub>s</sub>, H | Storey height, total height | m |
| I<sub>e</sub> | Importance factor | – |
| I<sub>y</sub>, I<sub>z</sub> | Second moments of area, local axes | m⁴ |
| J | Saint-Venant torsional constant | m⁴ |
| J<sub>m</sub> | Polar mass moment of inertia | t·m² |
| k<sub>i</sub> | Storey lateral stiffness | kN/m |
| L<sub>x</sub>, L<sub>y</sub> | Bay lengths | m |
| M<sub>OT</sub> | Overturning moment | kN·m |
| m | Mass | kN·s²/m = t |
| R | Response modification factor | – |
| S<sub>d</sub>, S<sub>a</sub> | Spectral displacement, acceleration | m, m/s² |
| T<sub>n</sub> | Period of mode *n* | s |
| T<sub>a</sub> | Empirical (code) period | s |
| V | Base shear | kN |
| W | Seismic weight | kN |
| w | Uniform line load | kN/m |
| Γ<sub>n</sub> | Modal participation factor | – |
| γ<sub>c</sub> | Unit weight of concrete | kN/m³ |
| ζ | Damping ratio | – |
| θ | Stability coefficient | – |
| ν | Poisson's ratio | – |
| φ<sub>n</sub> | Mode shape *n* | – |
| ω<sub>n</sub> | Circular frequency of mode *n* | rad/s |

---

## Appendix B — Reference console output

Running `python portico_3d_3niveles.py` with the default switches produces:

```
==========================================================================
 3-D RC SPACE FRAME  |  3 storeys  |  6.50 x 4.80 x 8.10 m  |  units: kN, m, s
==========================================================================

[1] SECTION PROPERTIES (gross)
  Ec = 21,538,106 kPa   Gc = 8,974,211 kPa   nu = 0.20
  Column 45x45: A=0.2025 m2  Iy=Iz=3.41719e-03 m4  J=5.76548e-03 m4
  Beam   30x60: A=0.1800 m2  Iy=5.40000e-03 m4  Iz=1.35000e-03 m4  J=3.70494e-03 m4
  Effective stiffness (ACI 318-19): OFF

[2] GRAVITY ANALYSIS (1.2D + 1.6L)
  status = 0 (0 = converged)
    level 1: uniform beam load w =  19.542 kN/m
    level 2: uniform beam load w =  19.542 kN/m
    level 3: uniform beam load w =  15.676 kN/m
  Sum of vertical reactions      =     1426.512 kN
  Total applied vertical load    =     1426.512 kN
  Equilibrium error              = 3.19e-14 %
  Max. vertical joint displacement (axial shortening) = -0.413 mm
  Verification, beam 1-2 at level 1 (span 6.50 m, w = 19.54 kN/m):
    end moments from FE analysis   |M| = 60.26 / 60.26 kNm
    ideal fixed-end value  wL^2/12     = 68.80 kNm
    simply-supported value wL^2/8      = 103.20 kNm
    mid-span moment (equilibrium)      = 42.95 kNm
    end shear from FE  V = 63.51 kN   vs  wL/2 = 63.51 kN

[3] MODAL ANALYSIS (9 modes)
   Mode      T [s]      f [Hz]     w [rad/s]
      1      0.3162      3.1629       19.873
      2      0.2939      3.4023       21.377
      3      0.1633      6.1247       38.482
      4      0.0962     10.3966       65.324
      5      0.0917     10.9027       68.504
      6      0.0539     18.5471      116.535
      7      0.0532     18.8062      118.163
      8      0.0516     19.3983      121.883
      9      0.0316     31.6225      198.690

   Participating mass ratios [%]
   Mode      UX        UY       RZ        (values in %)
      1    84.57     0.00     0.00
      2     0.00    85.36     0.00
      3     0.00     0.00    85.89
      4    12.05     0.00     0.00
      5     0.00    11.62     0.00
      6     3.38     0.00     0.00
      7     0.00     3.02     0.00
      8     0.00     0.00    11.15
      9     0.00     0.00     0.00
   Sum:   100.00   100.00    97.04   (>= 90 % required by ASCE 7-22 Sec. 12.9.1.1)

   Empirical period  Ta = 0.0466*h^0.90 = 0.3062 s
   Computed T1 / Ta = 1.03   (Cu*Ta limit with Cu=1.4: 0.4287 s)

[4] EQUIVALENT LATERAL FORCE, X DIRECTION
  Seismic weight W = 993.52 kN;  Cs = 0.100;  V = 99.35 kN;  k = 1.00
  Level    h [m]     F [kN]     ux [mm]   drift_e   drift_a=Cd/Ie   limit
     1      2.70      18.25      1.159   0.00043    0.00236   0.020  OK
     2      5.40      36.49      2.727   0.00058    0.00319   0.020  OK
     3      8.10      44.62      3.754   0.00038    0.00209   0.020  OK
  Base-shear check: sum(Rx) = 99.352 kN vs V = 99.352 kN  (error 1.00e-13 %)

  Local stress resultants under 1.2D + 1.6L + E_x
  Column 1, storey 1 (base end i / top end j):
    N   =    323.75 /   -308.00 kN
    Vy  =     -7.51 /      7.51 kN     Vz =    -10.52 /     10.52 kN
    T   =     -0.00 /      0.00 kNm
    My  =     32.88 /     -4.48 kNm    Mz =     -6.76 /    -13.52 kNm
  Beam 1-2, level 1 (local z vertical -> My is the flexural action):
    Vz  =     48.97 /     78.05 kN
    My  =    -13.01 /    107.51 kNm

[5] LINEAR TIME-HISTORY ANALYSIS
  Record: synthetic (Kanai-Tajimi);  n = 3000;  dt = 0.0100 s;  PGA = 0.350 g
  Rayleigh: alphaM = 1.03043 1/s, betaK = 0.002423 s  (zeta = 5% at modes 1-2)
  Peak roof displacement (X) = 32.157 mm = H/252
  Peak roof drift ratio      = 0.00397
==========================================================================
```

Small differences in the last decimal place are possible between platforms and library
versions. Differences in the first or second significant figure are not, and indicate a
problem worth tracking down.

---

## Appendix C — Formula summary

**Section properties (rectangle b × h):**

$$A = bh \qquad I_y = \frac{bh^3}{12} \qquad I_z = \frac{hb^3}{12} \qquad J = \beta a b^3$$

**Material:**

$$E_c = 4700\sqrt{f'_c}~[\text{MPa}] \qquad G = \frac{E}{2(1+\nu)}$$

**Mass:**

$$m = \frac{W}{g} \qquad J_m = \frac{m(L_x^2 + L_y^2)}{12}$$

**Local axis construction:**

$$\hat{\mathbf{y}}_{\text{loc}} = \hat{\mathbf{v}}_{xz} \times \hat{\mathbf{x}}_{\text{loc}}
\qquad
\hat{\mathbf{z}}_{\text{loc}} = \hat{\mathbf{x}}_{\text{loc}} \times \hat{\mathbf{y}}_{\text{loc}}$$

**Eigenproblem:**

$$(\mathbf{K} - \omega_n^2\mathbf{M})\boldsymbol{\phi}_n = \mathbf{0}
\qquad T_n = \frac{2\pi}{\omega_n}$$

**Rayleigh damping:**

$$a_0 = \frac{2\zeta\omega_1\omega_2}{\omega_1+\omega_2} \qquad
a_1 = \frac{2\zeta}{\omega_1+\omega_2} \qquad
\zeta(\omega) = \frac{a_0}{2\omega} + \frac{a_1\omega}{2}$$

**Equivalent lateral force:**

$$F_i = V\frac{w_ih_i^k}{\sum_j w_jh_j^k} \qquad
\delta_x = \frac{C_d\delta_{xe}}{I_e} \qquad
\theta = \frac{P_x\,\theta_{\text{drift},e}}{V_x}$$

**Verification identities:**

$$\sum R_z = \sum P_z \qquad
V_{\text{beam}} = \frac{wL}{2} \qquad
M_{\text{end}} + M_{\text{mid}} = \frac{wL^2}{8}$$

$$4\frac{3EI}{H_s^3} \le k_1 \le 4\frac{12EI}{H_s^3} \qquad
k_{\text{portal}} = \frac{12c(c+6b)}{H_s^2(2c+3b)}$$

$$T_{\text{Rayleigh}} = 2\pi\sqrt{\frac{\sum W_iu_i^2}{g\sum F_iu_i}} \qquad
u_{\text{roof}}^{\max} \approx \Gamma_1\phi_{1,\text{roof}}S_d(T_1,\zeta)$$

---

## References

1. ACI Committee 318. *Building Code Requirements for Structural Concrete (ACI 318-19) and
   Commentary.* American Concrete Institute, 2019.
2. ASCE/SEI 7-22. *Minimum Design Loads and Associated Criteria for Buildings and Other
   Structures.* American Society of Civil Engineers, 2022.
3. ASCE/SEI 41-17. *Seismic Evaluation and Retrofit of Existing Buildings.* American Society
   of Civil Engineers, 2017.
4. Chopra, A. K. *Dynamics of Structures: Theory and Applications to Earthquake
   Engineering.* 5th ed., Pearson.
5. Timoshenko, S. P. and Goodier, J. N. *Theory of Elasticity.* 3rd ed., McGraw-Hill —
   Saint-Venant torsion of rectangular sections.
6. OpenSeesPy documentation: https://openseespydoc.readthedocs.io/en/latest/

---

*Accompanying scripts:*
- `portico_3d_3niveles.py` — the model and the five analyses
- `verification_checks.py` — the nine verification checks of Part IV

*Prepared by Dr. Wahab.*
