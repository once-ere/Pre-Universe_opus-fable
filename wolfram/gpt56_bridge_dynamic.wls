(* gpt-5.6_bridge, dynamical form: the interpolation parameter is promoted from
   a constant to a scalar FIELD on spacetime.

   This is what separates gpt-5.6_bridge from Claude's Weitzenboeck connection
   and from every constant-parameter interpolation between them. Promoting the
   parameter costs nothing in consistency - the family stays a metric-compatible
   so(4,4) connection, because a difference of connections is a tensor and a
   scalar function times a tensor is a tensor - but it adds a curvature term
   proportional to the gradient of the field, and it makes the connection
   IRREMOVABLE from the Dirac operator.

   Run:  wolframscript -file wolfram/gpt56_bridge_dynamic.wls 2>&1 | tee logs/gpt56_bridge_dynamic.log
*)

ClearAll["Global`*"];

coordinates = {x0, x1, x2, x3, x4, x5, x6, x7};
dimension = Length[coordinates];
eta4488 = DiagonalMatrix[{1, 1, 1, 1, -1, -1, -1, -1}];
assumptions = H > 0 && 0 < 6 H x0 < Pi && x4 \[Element] Reals;

qFactor = Exp[-a4[H x4]]/Sin[6 H x0]^(1/6);
pFactor = Exp[a4[H x4]]/Sin[6 H x0]^(1/6);
frameCanonical = DiagonalMatrix[
  {Tan[6 H x0], qFactor, qFactor, qFactor, 1, pFactor, pFactor, pFactor}
];
inverseFrameCanonical = Simplify[Inverse[frameCanonical], assumptions];
gCanonical = Simplify[
  frameCanonical . eta4488 . Transpose[frameCanonical], assumptions];
inverseGCanonical = Simplify[Inverse[gCanonical], assumptions];

(* A literal 0 needs no FullSimplify; AllTrue and AnyTrue stop at the first
   decisive element. Truth values are unchanged; only work is saved. *)
zeroQ[0] := True;
zeroQ[e_] := TrueQ[FullSimplify[e == 0, assumptions]];
zeroArrayQ[e_] := AllTrue[Flatten[{e}], zeroQ];
nonzeroWitnessQ[e_] := AnyTrue[Flatten[{e}], ! zeroQ[#] &];

christoffelFromMetric[metric_, inverseMetric_] := Simplify[
  Table[
    1/2 Sum[
      inverseMetric[[rho, sigma]] (
        D[metric[[sigma, nu]], coordinates[[mu]]] +
        D[metric[[sigma, mu]], coordinates[[nu]]] -
        D[metric[[mu, nu]], coordinates[[sigma]]]),
      {sigma, dimension}],
    {rho, dimension}, {mu, dimension}, {nu, dimension}],
  assumptions];

spinConnectionFromAffine[frame_, inverseFrame_, affine_] := Simplify[
  Table[
    Sum[
      inverseFrame[[b, nu]] (
        Sum[affine[[rho, mu, nu]] frame[[rho, a]], {rho, dimension}] -
        D[frame[[nu, a]], coordinates[[mu]]]),
      {nu, dimension}],
    {mu, dimension}, {a, dimension}, {b, dimension}],
  assumptions];

affineFromFrameAndSpin[frame_, inverseFrame_, spinConnection_] := Simplify[
  Table[
    Sum[
      inverseFrame[[a, rho]] (
        D[frame[[nu, a]], coordinates[[mu]]] +
        Sum[spinConnection[[mu, a, b]] frame[[nu, b]], {b, dimension}]),
      {a, dimension}],
    {rho, dimension}, {mu, dimension}, {nu, dimension}],
  assumptions];

lowerSpinIndices[spinConnection_] := Simplify[
  Table[
    Sum[eta4488[[a, c]] spinConnection[[mu, c, b]], {c, dimension}],
    {mu, dimension}, {a, dimension}, {b, dimension}],
  assumptions];

vielbeinPostulateResidual[frame_, affine_, spinConnection_] := Simplify[
  Table[
    D[frame[[nu, a]], coordinates[[mu]]] -
    Sum[affine[[rho, mu, nu]] frame[[rho, a]], {rho, dimension}] +
    Sum[spinConnection[[mu, a, b]] frame[[nu, b]], {b, dimension}],
    {mu, dimension}, {nu, dimension}, {a, dimension}],
  assumptions];

metricCompatibilityResidual[metric_, affine_] := Simplify[
  Table[
    D[metric[[nu, sigma]], coordinates[[mu]]] -
    Sum[
      affine[[rho, mu, nu]] metric[[rho, sigma]] +
      affine[[rho, mu, sigma]] metric[[nu, rho]],
      {rho, dimension}],
    {mu, dimension}, {nu, dimension}, {sigma, dimension}],
  assumptions];

torsionTensor[affine_] := Simplify[
  Table[affine[[rho, mu, nu]] - affine[[rho, nu, mu]],
    {rho, dimension}, {mu, dimension}, {nu, dimension}],
  assumptions];

(* Every four-index object below is antisymmetric in its last index pair by
   the form of its defining expression, so only the mu < nu entries are ever
   simplified; the rest are filled by antisymmetry. This halves the dominant
   cost of the file without changing a single component. *)
antisymmetricInLastPair[component_, n_] := Module[{upper},
  upper = Table[
    If[mu < nu, component[i, j, mu, nu], 0],
    {i, n}, {j, n}, {mu, n}, {nu, n}];
  upper - Transpose[upper, {1, 2, 4, 3}]];

curvatureComponent[affine_, rho_, sigma_, mu_, nu_] := Simplify[
  D[affine[[rho, nu, sigma]], coordinates[[mu]]] -
  D[affine[[rho, mu, sigma]], coordinates[[nu]]] +
  Sum[
    affine[[rho, mu, lambda]] affine[[lambda, nu, sigma]] -
    affine[[rho, nu, lambda]] affine[[lambda, mu, sigma]],
    {lambda, dimension}],
  assumptions];

curvatureTensor[affine_] := antisymmetricInLastPair[
  curvatureComponent[affine, #1, #2, #3, #4] &, dimension];

covariantExteriorDerivativeComponent[base_, difference_, rho_, sigma_, mu_, nu_] :=
  Simplify[
    D[difference[[rho, nu, sigma]], coordinates[[mu]]] -
    D[difference[[rho, mu, sigma]], coordinates[[nu]]] +
    Sum[
      base[[rho, mu, lambda]] difference[[lambda, nu, sigma]] +
      difference[[rho, mu, lambda]] base[[lambda, nu, sigma]] -
      base[[rho, nu, lambda]] difference[[lambda, mu, sigma]] -
      difference[[rho, nu, lambda]] base[[lambda, mu, sigma]],
      {lambda, dimension}],
    assumptions];

covariantExteriorDerivative[base_, difference_] := antisymmetricInLastPair[
  covariantExteriorDerivativeComponent[base, difference, #1, #2, #3, #4] &,
  dimension];

connectionWedgeComponent[difference_, rho_, sigma_, mu_, nu_] := Simplify[
  Sum[
    difference[[rho, mu, lambda]] difference[[lambda, nu, sigma]] -
    difference[[rho, nu, lambda]] difference[[lambda, mu, sigma]],
    {lambda, dimension}],
  assumptions];

connectionWedge[difference_] := antisymmetricInLastPair[
  connectionWedgeComponent[difference, #1, #2, #3, #4] &, dimension];

(* The totally antisymmetric (axial) part of the fully lowered torsion. *)
lowerTorsion[metric_, torsion_] := Simplify[
  Table[
    Sum[metric[[rho, lambda]] torsion[[lambda, mu, nu]], {lambda, dimension}],
    {rho, dimension}, {mu, dimension}, {nu, dimension}],
  assumptions];

axialPart[loweredTorsion_] := Simplify[
  Table[
    (1/3) (loweredTorsion[[rho, mu, nu]] + loweredTorsion[[mu, nu, rho]] +
           loweredTorsion[[nu, rho, mu]]),
    {rho, dimension}, {mu, dimension}, {nu, dimension}],
  assumptions];

traceVector[metric_, torsion_] := Simplify[
  Table[Sum[torsion[[nu, nu, mu]], {nu, dimension}], {mu, dimension}],
  assumptions];

(* ------------------------------------------------------------------ *)
(* Canonical and Weitzenboeck geometry.                                *)
(* ------------------------------------------------------------------ *)

gammaCanonical = christoffelFromMetric[gCanonical, inverseGCanonical];
omegaCanonicalMixed = spinConnectionFromAffine[
  frameCanonical, inverseFrameCanonical, gammaCanonical];
gammaWeitzenboeck = affineFromFrameAndSpin[
  frameCanonical, inverseFrameCanonical,
  ConstantArray[0, {dimension, dimension, dimension}]];

contortion = Simplify[gammaWeitzenboeck - gammaCanonical, assumptions];
torsionWeitzenboeck = torsionTensor[gammaWeitzenboeck];
curvatureCanonical = curvatureTensor[gammaCanonical];

(* ------------------------------------------------------------------ *)
(* The bridge field. h in [0,1] dials Levi-Civita -> Weitzenboeck.     *)
(* bridgeField is an arbitrary smooth function of the two              *)
(* inhomogeneous coordinates; nothing below assumes any special form.  *)
(* ------------------------------------------------------------------ *)

bridgeH = 3 bridgeField[x0, x4]^2 - 2 bridgeField[x0, x4]^3;

gpt56BridgeAffineDynamic = Simplify[
  gammaCanonical + bridgeH contortion, assumptions];
gpt56BridgeSpinDynamic = spinConnectionFromAffine[
  frameCanonical, inverseFrameCanonical, gpt56BridgeAffineDynamic];

(* A scalar function times a tensor is a tensor, so the promoted family is
   still a connection: it satisfies the vielbein postulate against its own
   affine partner and is still metric compatible and so(4,4)-valued. *)
dynamicPostulateResidual = vielbeinPostulateResidual[
  frameCanonical, gpt56BridgeAffineDynamic, gpt56BridgeSpinDynamic];
dynamicMetricResidual = metricCompatibilityResidual[
  gCanonical, gpt56BridgeAffineDynamic];
dynamicAntisymmetryResidual = With[
  {lowered = lowerSpinIndices[gpt56BridgeSpinDynamic]},
  Simplify[lowered + Transpose[lowered, {1, 3, 2}], assumptions]];

(* The spin connection is exactly the canonical one scaled by 1 - h. *)
dynamicScalingResidual = Simplify[
  gpt56BridgeSpinDynamic - (1 - bridgeH) omegaCanonicalMixed, assumptions];

(* ------------------------------------------------------------------ *)
(* Exact transgression identities with a NON-constant parameter.       *)
(* ------------------------------------------------------------------ *)

gpt56BridgeTorsionDynamic = torsionTensor[gpt56BridgeAffineDynamic];
torsionTransgressionResidual = Simplify[
  gpt56BridgeTorsionDynamic - bridgeH torsionWeitzenboeck, assumptions];

gpt56BridgeCurvatureDynamic = curvatureTensor[gpt56BridgeAffineDynamic];
covariantContortion = covariantExteriorDerivative[gammaCanonical, contortion];
contortionWedge = connectionWedge[contortion];

(* The gradient term that a constant parameter cannot produce. *)
gradientTerm = Table[
  Simplify[
    D[bridgeH, coordinates[[mu]]] contortion[[rho, nu, sigma]] -
    D[bridgeH, coordinates[[nu]]] contortion[[rho, mu, sigma]],
    assumptions],
  {rho, dimension}, {sigma, dimension}, {mu, dimension}, {nu, dimension}];

curvatureTransgressionResidual = Simplify[
  gpt56BridgeCurvatureDynamic - (
    curvatureCanonical + bridgeH covariantContortion +
    bridgeH^2 contortionWedge + gradientTerm),
  assumptions];

(* ------------------------------------------------------------------ *)
(* Torsion decomposition. Reported honestly: the axial part vanishes   *)
(* for Claude's Weitzenboeck connection AND for every member of this   *)
(* family, because the family's torsion is a scalar multiple of it.    *)
(* ------------------------------------------------------------------ *)

loweredWeitzenboeckTorsion = lowerTorsion[gCanonical, torsionWeitzenboeck];
axialWeitzenboeck = axialPart[loweredWeitzenboeckTorsion];
loweredDynamicTorsion = lowerTorsion[gCanonical, gpt56BridgeTorsionDynamic];
axialDynamic = axialPart[loweredDynamicTorsion];
torsionVectorWeitzenboeck = traceVector[gCanonical, torsionWeitzenboeck];
torsionPotential = Log[Sin[6 H x0]];
torsionVectorIsGradient = Simplify[
  torsionVectorWeitzenboeck -
    Table[D[torsionPotential, coordinates[[mu]]], {mu, dimension}],
  assumptions];

(* ------------------------------------------------------------------ *)
(* The discriminator: is the connection removable from the Dirac       *)
(* operator by a spinor rescaling?                                     *)
(*                                                                      *)
(* In the Dirac operator the whole family contributes the single scalar *)
(* one-form -(1/2)(1 - h) d(torsionPotential). Such a term is absorbed  *)
(* by Psi -> Exp[-F] Psi if and only if that one-form is EXACT, i.e.    *)
(* if and only if the obstruction two-form d[(1 - h) dT] vanishes.      *)
(* ------------------------------------------------------------------ *)

diracOneForm = Table[
  Simplify[-(1/2) (1 - bridgeH) D[torsionPotential, coordinates[[mu]]],
    assumptions],
  {mu, dimension}];

obstruction = Table[
  Simplify[
    D[diracOneForm[[nu]], coordinates[[mu]]] -
    D[diracOneForm[[mu]], coordinates[[nu]]],
    assumptions],
  {mu, dimension}, {nu, dimension}];

obstructionIsWedge = Simplify[
  obstruction - Table[
    (1/2) (
      D[bridgeH, coordinates[[mu]]] D[torsionPotential, coordinates[[nu]]] -
      D[bridgeH, coordinates[[nu]]] D[torsionPotential, coordinates[[mu]]]),
    {mu, dimension}, {nu, dimension}],
  assumptions];

(* Case 1: constant bridge field - the canonical and Weitzenboeck endpoints
   and every constant interpolation between them. The obstruction vanishes and
   the explicit rescaling exists. *)
constantRule = {bridgeField -> Function[{u, v}, tConst]};
obstructionConstant = Simplify[obstruction /. constantRule, assumptions];
rescalingExponent = -(1/2) (1 - (3 tConst^2 - 2 tConst^3)) torsionPotential;
rescalingResidual = Simplify[
  (diracOneForm /. constantRule) -
    Table[D[rescalingExponent, coordinates[[mu]]], {mu, dimension}],
  assumptions];

(* Case 2: a bridge field that varies along a direction the torsion potential
   does not. The obstruction is nonzero, so NO spinor rescaling can remove the
   connection: it is physically visible. *)
dynamicRule = {bridgeField -> Function[{u, v}, 1/2 + Tanh[v]/4]};
obstructionDynamic = Simplify[obstruction /. dynamicRule, assumptions];
obstructionWitness = FullSimplify[obstructionDynamic[[1, 5]], assumptions];

tests = {
  VerificationTest[
    zeroArrayQ[gCanonical - frameCanonical . eta4488 . Transpose[frameCanonical]],
    True, TestID -> "dynamic-canonical-frame-reconstructs-metric"],
  VerificationTest[zeroArrayQ[dynamicPostulateResidual], True,
    TestID -> "dynamic-bridge-satisfies-vielbein-postulate"],
  VerificationTest[zeroArrayQ[dynamicMetricResidual], True,
    TestID -> "dynamic-bridge-is-metric-compatible"],
  VerificationTest[zeroArrayQ[dynamicAntisymmetryResidual], True,
    TestID -> "dynamic-bridge-spin-connection-is-so-4-4-valued"],
  VerificationTest[zeroArrayQ[dynamicScalingResidual], True,
    TestID -> "dynamic-bridge-spin-connection-is-one-minus-h-times-canonical"],
  VerificationTest[zeroArrayQ[torsionTransgressionResidual], True,
    TestID -> "dynamic-torsion-transgression-holds-for-a-field"],
  VerificationTest[zeroArrayQ[curvatureTransgressionResidual], True,
    TestID -> "dynamic-curvature-transgression-with-gradient-term"],
  VerificationTest[nonzeroWitnessQ[gradientTerm], True,
    TestID -> "dynamic-gradient-term-is-nonzero"],
  VerificationTest[
    zeroArrayQ[gradientTerm /. {bridgeField -> Function[{u, v}, tConst]}],
    True,
    TestID -> "dynamic-gradient-term-vanishes-for-a-constant-parameter"],
  VerificationTest[zeroArrayQ[axialWeitzenboeck], True,
    TestID -> "weitzenboeck-axial-torsion-vanishes"],
  VerificationTest[zeroArrayQ[axialDynamic], True,
    TestID -> "dynamic-bridge-axial-torsion-also-vanishes"],
  VerificationTest[zeroArrayQ[torsionVectorIsGradient], True,
    TestID -> "weitzenboeck-torsion-vector-is-a-gradient"],
  VerificationTest[zeroArrayQ[obstructionIsWedge], True,
    TestID -> "obstruction-equals-dh-wedge-d-torsion-potential"],
  VerificationTest[zeroArrayQ[obstructionConstant], True,
    TestID -> "obstruction-vanishes-for-every-constant-parameter"],
  VerificationTest[zeroArrayQ[rescalingResidual], True,
    TestID -> "constant-parameter-connection-is-an-exact-gradient"],
  VerificationTest[zeroQ[obstructionWitness], False,
    TestID -> "dynamic-bridge-obstruction-is-nonzero"],
  VerificationTest[zeroArrayQ[obstructionDynamic], False,
    TestID -> "dynamic-bridge-is-irremovable-from-the-dirac-operator"]
};

report = TestReport[tests];
Print["Tests run    : ", Length[tests]];
Print["Succeeded    : ", report["TestsSucceededCount"]];
Print["Failed       : ", report["TestsFailedCount"]];
Print[];
Print["Weitzenboeck torsion vector  : ", torsionVectorWeitzenboeck[[1]],
  " (and zero elsewhere)"];
Print["Axial torsion, Weitzenboeck  : 0 (identically)"];
Print["Axial torsion, gpt-5.6_bridge: 0 (identically; the torsion is a scalar",
  " multiple of the Weitzenboeck torsion)"];
Print["Obstruction witness d[(1-h)dT]_{x0 x4} for the sample bridge field",
  " s = 1/2 + Tanh[x4]/4, h = 3 s^2 - 2 s^3:"];
Print["  ", obstructionWitness];
Print[];
Print["Reading: the canonical connection, Claude's Weitzenboeck connection and",
  " every constant-parameter"];
Print["bridge between them contribute an EXACT one-form to the Dirac operator",
  " and can therefore be"];
Print["removed by a spinor rescaling. A bridge FIELD that varies transverse to",
  " the torsion potential"];
Print["cannot: the obstruction two-form is nonzero, so the connection is",
  " physically visible."];

If[report["TestsFailedCount"] > 0 && $FrontEnd === Null, Exit[1]];
