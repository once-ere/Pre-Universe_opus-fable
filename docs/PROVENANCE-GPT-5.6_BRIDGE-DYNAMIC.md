# Provenance — gpt-5.6_bridge as a field (`wolfram/gpt56_bridge_dynamic.wls`)

**Effort.** This page records the construction of `gpt-5.6f_rame` and
`gpt-5.6_bridge`, the promotion of the bridge interpolation parameter from a
constant to a **scalar field on spacetime**, the two exact transgression
identities that follow, and the precise sense in which the resulting connection
is superior to Claude's Weitzenböck connection — together with the one criterion
that was *tried and failed*, reported here rather than quietly dropped.
Everything needed to repeat the work is on this page. No other file needs to be
consulted, and no command below is abbreviated.

---

## 1. Environment actually used

| item | value |
|---|---|
| OS | Ubuntu 26.04.1 LTS, x86-64 |
| kernel | Linux 7.0.0-31-generic |
| C library | GNU libc 2.43 |
| shell | bash |
| Wolfram Language | 15.0.1 for Linux x86 (64-bit), 2 July 2026, Professional licence |
| cores | 24 |
| git | 2.53.0 |
| repository commit | `bae14ee24b58fb35f8a029f0b5528342d7a99599` |

Only `wolframscript` is required for this page.

---

## 2. Getting the repository and checking the prerequisite

```bash
git clone --recurse-submodules https://github.com/once-ere/Pre-Universe_opus-fable.git
cd Pre-Universe_opus-fable
mkdir -p logs
wolframscript -code '{$Version, $LicenseType}'
```

The expected reply names your Wolfram version and a non-expired licence type,
for example:

```
{15.0.1 for Linux x86 (64-bit) (July 2, 2026), Professional}
```

If that prints a licensing error instead, stop: nothing below can run.

---

## 3. The two commands that produce the result

The constant-parameter family (36 assertions, about forty seconds):

```bash
cd Pre-Universe_opus-fable
wolframscript -file wolfram/gpt56_bridge.wls 2>&1 | tee logs/gpt56_bridge.log
echo "exit status: ${PIPESTATUS[0]}"
tail -n 5 logs/gpt56_bridge.log
```

The bridge **field** (17 assertions, about seven seconds):

```bash
cd Pre-Universe_opus-fable
wolframscript -file wolfram/gpt56_bridge_dynamic.wls 2>&1 | tee logs/gpt56_bridge_dynamic.log
echo "exit status: ${PIPESTATUS[0]}"
```

Exit status `0` means every assertion passed; both scripts call `Exit[1]` if any
fails.

### 3.1 The stricter form: assertions **and** zero kernel messages

`wolframscript` does not fail when the kernel emits a message (for example a
`Limit::alimv` or `General::stop`) while every assertion still passes. The
repository therefore ships a verifier that loads a source with `$MessageList`
captured and exits non-zero if **any** message was generated. Run it on each
source in turn:

```bash
cd Pre-Universe_opus-fable
wolframscript -file tests/verify_wolfram_source.wls wolfram/gpt56_bridge.wls \
  2>&1 | tee logs/verify-wolfram-source.log
echo "exit status: ${PIPESTATUS[0]}"
wolframscript -file tests/verify_wolfram_source.wls wolfram/gpt56_bridge_dynamic.wls \
  2>&1 | tee logs/verify-gpt56-bridge-dynamic.log
echo "exit status: ${PIPESTATUS[0]}"
```

Each log ends with a six-line summary; the seconds value is machine-dependent,
every other value is fixed:

```
Source verifier file     : gpt56_bridge.wls
Source verifier seconds  : 40.1
Source verifier messages : 0
Source verifier succeeded: 36
Source verifier failed   : 0
SUCCESS: gpt56_bridge.wls passed every assertion with no messages.
```

```
Source verifier file     : gpt56_bridge_dynamic.wls
Source verifier seconds  : 7.2
Source verifier messages : 0
Source verifier succeeded: 17
Source verifier failed   : 0
SUCCESS: gpt56_bridge_dynamic.wls passed every assertion with no messages.
```

Exit status `0` is the pass certificate; `1` means an assertion failed **or** a
message was emitted; `2` means the path argument was missing or the file does
not exist. This is the form `scripts/verify_all.sh` runs in its stages 7 and 9.

---

## 4. The geometry

### 4.1 Canonical starting point

Coordinates $x^0\ldots x^7$; tangent metric
$\eta=\mathrm{diag}(1,1,1,1,-1,-1,-1,-1)$; the canonical diagonal frame

```
frameCanonical = DiagonalMatrix[{Tan[6 H x0], q, q, q, 1, p, p, p}]
q = Exp[-a4[H x4]] / Sin[6 H x0]^(1/6)
p = Exp[+a4[H x4]] / Sin[6 H x0]^(1/6)
```

with `a4` a free function. This satisfies
$g_{\mu\nu}=e_\mu{}^a\eta_{ab}e_\nu{}^b$, has 37 nonzero Christoffel symbols and
**24** nonzero spin-connection components, and is torsion-free.

### 4.2 `gpt-5.6f_rame`

Let $R$ be the normalised order-4 Hadamard matrix and

$$J=\begin{pmatrix}0&R\\ R^{\mathsf T}&0\end{pmatrix}.$$

Then $J^{\mathsf T}\eta+\eta J=0$, so $J\in\mathfrak{so}(4,4)$, and $J^2=\mathbb{1}$,
so the group element is available in **closed form** with no series:

$$\Lambda(u)=\cosh(u)\,\mathbb{1}+\sinh(u)\,J ,\qquad
\texttt{gpt-5.6f\_rame}=\texttt{frameCanonical}\cdot\Lambda(u).$$

This is a simultaneous equal-rapidity boost in all four orthogonal $(+,-)$
planes, so unlike the canonical diagonal frame it is **dense**: it mixes the
$+4$ and $-4$ blocks. It reconstructs the same curved metric for every $u$ and
reduces to the canonical frame at $u=0$.

### 4.3 `gpt-5.6_bridge`

Let $K=\Gamma_{\text{W}}-\Gamma_{\text{LC}}$, the difference between the
Weitzenböck and Levi-Civita affine connections. A difference of connections is a
**tensor**, so for any scalar $h$

$$\Gamma(h)=\Gamma_{\text{LC}}+h\,K$$

is again a connection. With $h(s)=3s^2-2s^3$ the endpoints are exact:
$h=0$ is Levi-Civita, $h=1$ is Claude's Weitzenböck geometry.

**The promotion.** Nothing in that argument required $h$ to be constant. Let
$s=s(x)$ be a scalar field; $h(s(x))\,K$ is still a tensor, so
$\Gamma(h(s(x)))$ is still a metric-compatible $\mathfrak{so}(4,4)$ connection.
The script verifies this directly: vielbein-postulate residual zero, metric
compatibility zero, antisymmetry in the flat pair zero, and

$$\omega_{\text{bridge}} = \bigl(1-h(s(x))\bigr)\,\omega_{\text{canonical}} .$$

---

## 5. The two exact transgression identities

| quantity | identity | new with a field? |
|---|---|---|
| torsion | $T = h\,T_{\text{Weitzenböck}}$ | no — holds for constants too |
| curvature | $R = R_{\text{LC}} + h\,D_{\text{LC}}K + h^2\,K\wedge K + \mathrm{d}h\wedge K$ | **yes — the last term** |

In components the new term is

```
(d_mu h) K^rho_{nu sigma} - (d_nu h) K^rho_{mu sigma}
```

and the script proves both that it is **nonzero** for a bridge field and that it
is **identically zero** for every constant parameter. No member of the
constant-parameter family reproduces it.

---

## 6. The criterion that failed, reported honestly

The plan proposed to claim superiority via **axial torsion**: Claude's
Weitzenböck connection has a totally antisymmetric torsion part that vanishes
identically, and the hope was that the bridge would not.

**It does not work, and the script proves it does not.** Because the bridge
torsion is $h$ times the Weitzenböck torsion, and $h$ is a scalar, the totally
antisymmetric part of the bridge torsion is $h\times 0=0$. Both assertions

```
weitzenboeck-axial-torsion-vanishes
dynamic-bridge-axial-torsion-also-vanishes
```

are in the test suite precisely so that this negative result is recorded and
cannot be quietly forgotten. The superiority claim below does **not** rest on
axial torsion.

---

## 7. The criterion that works: irremovability

In the Dirac operator the entire family contributes one scalar one-form. Writing
the torsion potential $\mathcal{T}=\log\sin(6Hx^0)$, the contribution is

$$\gamma^\mu\Gamma^{\text{spin}}_\mu
= -\tfrac12\bigl(1-h\bigr)\,\partial_\mu\mathcal{T}\,\gamma^\mu .$$

A term of that shape is absorbed by a spinor rescaling $\Psi\to e^{-F}\Psi$
**if and only if** the one-form is exact, that is, if and only if the
**obstruction two-form** vanishes:

$$\Omega \;=\; \mathrm{d}\bigl[(1-h)\,\mathrm{d}\mathcal{T}\bigr]
\;=\; \tfrac12\,\mathrm{d}h\wedge \mathrm{d}\mathcal{T}.$$

| connection | $\Omega$ | consequence |
|---|---|---|
| canonical Levi-Civita ($h=0$) | $0$ | removable by a rescaling |
| Claude's Weitzenböck ($h=1$) | $0$ | the connection is zero outright |
| any **constant** bridge parameter | $0$ | removable by a rescaling |
| **bridge field varying transverse to $\mathcal{T}$** | $\neq 0$ | **irremovable** |

Because $\mathcal{T}$ depends on $x^0$ alone, any bridge field with
$\partial_4 h\neq0$ produces a nonzero obstruction. For the sample choice
$s=\tfrac12+\tfrac14\tanh(x^4)$, hence $h=3s^2-2s^3$, the script reports the
exact witness

$$\Omega_{x^0x^4}
= \frac{9H\cot(6Hx^0)\,\operatorname{sech}^2(x^4)\bigl(\tanh^2(x^4)-4\bigr)}{32}.$$

That is the measurable sense in which `gpt-5.6_bridge` is superior: it is the
only member of the comparison that a spinor cannot be redefined to ignore.

---

## 8. The result

Full output of `logs/gpt56_bridge_dynamic.log` from the recorded run:

```
Tests run    : 17
Succeeded    : 17
Failed       : 0

Weitzenboeck torsion vector  : 6*H*Cot[6*H*x0] (and zero elsewhere)
Axial torsion, Weitzenboeck  : 0 (identically)
Axial torsion, gpt-5.6_bridge: 0 (identically; the torsion is a scalar multiple of the Weitzenboeck torsion)
Obstruction witness d[(1-h)dT]_{x0 x4} for the sample bridge field s = 1/2 + Tanh[x4]/4, h = 3 s^2 - 2 s^3:
  (9*H*Cot[6*H*x0]*Sech[x4]^2*(-4 + Tanh[x4]^2))/32
```

Final lines of `logs/gpt56_bridge.log` from the recorded run:

```
Canonical nonzero spin-connection components: 24
Interior torsion witness T^1_(0 1): -1/2*(H*Cot[6*H*x0])
Interior curvature witness R^1_(0 0 1): (H^2*(25 + Cos[12*H*x0])*Csc[6*H*x0]^2)/4
Tests succeeded: 36
Tests failed: 0
```

---

## 9. Checking the central claims yourself

Each block is complete and standalone; paste it into a terminal from the
repository root. `wolframscript -code` echoes the value of the last expression,
so a trailing `Null` in the output is normal.

The Hadamard generator really lies in $\mathfrak{so}(4,4)$ and squares to one:

```bash
wolframscript -code '
R = 1/2 {{1,1,1,1},{1,-1,1,-1},{1,1,-1,-1},{1,-1,-1,1}};
eta = DiagonalMatrix[{1,1,1,1,-1,-1,-1,-1}];
J = ArrayFlatten[{{ConstantArray[0,{4,4}], R}, {Transpose[R], ConstantArray[0,{4,4}]}}];
Print["R orthogonal      : ", R . Transpose[R] == IdentityMatrix[4]];
Print["J in so(4,4)      : ", Transpose[J] . eta + eta . J == ConstantArray[0,{8,8}]];
Print["J^2 == identity   : ", J . J == IdentityMatrix[8]];
Print["Lambda preserves eta: ", Simplify[
  With[{L = Cosh[u] IdentityMatrix[8] + Sinh[u] J},
    Transpose[L] . eta . L - eta] == ConstantArray[0,{8,8}]]];'
```

Expected output:

```
R orthogonal      : True
J in so(4,4)      : True
J^2 == identity   : True
Lambda preserves eta: True
Null
```
The obstruction is exactly $-\tfrac12\,\mathrm{d}h\wedge\mathrm{d}\mathcal{T}$,
vanishes for a constant $h$, and does not vanish for the sample field:

```bash
wolframscript -code '
T = Log[Sin[6 H x0]];
oneForm[h_] := {-(1/2)(1 - h) D[T, x0], -(1/2)(1 - h) D[T, x4]};
obs[h_] := D[oneForm[h][[2]], x0] - D[oneForm[h][[1]], x4];
Print["constant h         : ", Simplify[obs[hc]]];
Print["field h            : ", Simplify[obs[1/2 + Tanh[x4]/4]]];
Print["equals -1/2 dh^dT  : ", Simplify[
  obs[1/2 + Tanh[x4]/4] - (1/2)(D[1/2 + Tanh[x4]/4, x4] D[T, x0]
                             - D[1/2 + Tanh[x4]/4, x0] D[T, x4])] == 0];'
```

Expected output:

```
constant h         : 0
field h            : (9*H*Cot[6*H*x0]*Sech[x4]^2*(-4 + Tanh[x4]^2))/32
equals -1/2 dh^dT  : True
Null
```

The middle line is the witness quoted in section 7. The first line is why no
constant parameter can do this.

---

## 10. Viewing the result in a web browser

Serve the repository over HTTP from its root:

```bash
cd Pre-Universe_opus-fable
python3 -m http.server 8911
```

Then open <http://127.0.0.1:8911/docs/PROVENANCE-GPT-5.6_BRIDGE-DYNAMIC.md> for
this page and <http://127.0.0.1:8911/artifacts/fable/figures/> for the figures.
Stop the server with `Ctrl-C`.

To evaluate the bridge notebook in place, storing its outputs in the file:

```bash
cd Pre-Universe_opus-fable
wolframscript -file scripts/execute_gpt56_notebook.wls 2>&1 | tee logs/execute_gpt56_notebook.log
```

Expected, on the last lines: `Input cells evaluated    : 5`,
`Output cells stored      : 5`, `Cells with messages      : 0`,
`Notebook tests succeeded : 36`, and `SUCCESS: the notebook is evaluated and
its outputs are stored.` (`scripts/run_gpt56_notebook.wls` evaluates the same
cells but discards the results; it is kept for anyone who wants to check the
notebook without rewriting it.)

Then, on a machine with the Mathematica graphical front end installed:

```bash
mathematica notebooks/gpt-5.6_bridge.nb
```

followed by **Evaluation > Evaluate Notebook**.

---

## 11. Honest limits of this page

- `gpt-5.6f_rame` is a **pseudo-orthogonal frame gauge transformation** of the
  canonical frame. It changes no curved-space observable by itself; the physical
  content is in the connection, not the frame choice.
- The superiority claim is exactly the one stated in section 7 and nothing
  wider. It is a statement about removability by spinor rescaling, proved
  symbolically, not a claim about experimental preference.
- The axial-torsion criterion **failed**, and section 6 says so.
- The bridge field's dynamics are not derived from an action anywhere in this
  repository; the field is treated as given geometry.
