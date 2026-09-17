(* gpt-5.6_bridge: a verified interpolation between Levi-Civita and
   Weitzenboeck geometry in signature (4,4). *)

(* ::NotebookSection:: Canonical 4+4 geometry *)
ClearAll["Global`*"];

coordinates = {x0, x1, x2, x3, x4, x5, x6, x7};
dimension = Length[coordinates];
eta4488 = DiagonalMatrix[{1, 1, 1, 1, -1, -1, -1, -1}];
assumptions = H > 0 && kappa \[Element] Reals &&
  t \[Element] Reals && 0 < 6 H x0 < Pi;

qFactor = Exp[-a4[H x4]]/Sin[6 H x0]^(1/6);
pFactor = Exp[a4[H x4]]/Sin[6 H x0]^(1/6);
frameCanonical = DiagonalMatrix[
  {Tan[6 H x0], qFactor, qFactor, qFactor, 1, pFactor, pFactor, pFactor}
];
inverseFrameCanonical = Simplify[Inverse[frameCanonical], assumptions];
gCanonical = Simplify[
  frameCanonical . eta4488 . Transpose[frameCanonical], assumptions
];
inverseGCanonical = Simplify[Inverse[gCanonical], assumptions];

zeroQ[expression_] := TrueQ[FullSimplify[expression == 0, assumptions]];
zeroArrayQ[expression_] := And @@ (zeroQ /@ Flatten[expression]);
nonzeroCount[expression_] := Count[
  FullSimplify[#, assumptions] & /@ Flatten[expression], Except[0]
];

christoffelFromMetric[metric_, inverseMetric_] := Simplify[
  Table[
    1/2 Sum[
      inverseMetric[[rho, sigma]] (
        D[metric[[sigma, nu]], coordinates[[mu]]] +
        D[metric[[sigma, mu]], coordinates[[nu]]] -
        D[metric[[mu, nu]], coordinates[[sigma]]]
      ),
      {sigma, dimension}
    ],
    {rho, dimension}, {mu, dimension}, {nu, dimension}
  ],
  assumptions
];

spinConnectionFromAffine[frame_, inverseFrame_, affine_] := Simplify[
  Table[
    Sum[
      inverseFrame[[b, nu]] (
        Sum[affine[[rho, mu, nu]] frame[[rho, a]], {rho, dimension}] -
        D[frame[[nu, a]], coordinates[[mu]]]
      ),
      {nu, dimension}
    ],
    {mu, dimension}, {a, dimension}, {b, dimension}
  ],
  assumptions
];

lowerSpinIndices[spinConnection_] := Simplify[
  Table[
    Sum[eta4488[[a, c]] spinConnection[[mu, c, b]], {c, dimension}],
    {mu, dimension}, {a, dimension}, {b, dimension}
  ],
  assumptions
];

affineFromFrameAndSpin[frame_, inverseFrame_, spinConnection_] := Simplify[
  Table[
    Sum[
      inverseFrame[[a, rho]] (
        D[frame[[nu, a]], coordinates[[mu]]] +
        Sum[
          spinConnection[[mu, a, b]] frame[[nu, b]],
          {b, dimension}
        ]
      ),
      {a, dimension}
    ],
    {rho, dimension}, {mu, dimension}, {nu, dimension}
  ],
  assumptions
];

vielbeinPostulateResidual[frame_, affine_, spinConnection_] := Simplify[
  Table[
    D[frame[[nu, a]], coordinates[[mu]]] -
    Sum[affine[[rho, mu, nu]] frame[[rho, a]], {rho, dimension}] +
    Sum[
      spinConnection[[mu, a, b]] frame[[nu, b]],
      {b, dimension}
    ],
    {mu, dimension}, {nu, dimension}, {a, dimension}
  ],
  assumptions
];

torsionTensor[affine_] := Simplify[
  Table[
    affine[[rho, mu, nu]] - affine[[rho, nu, mu]],
    {rho, dimension}, {mu, dimension}, {nu, dimension}
  ],
  assumptions
];

curvatureComponent[affine_, rho_, sigma_, mu_, nu_] := Simplify[
  D[affine[[rho, nu, sigma]], coordinates[[mu]]] -
  D[affine[[rho, mu, sigma]], coordinates[[nu]]] +
  Sum[
    affine[[rho, mu, lambda]] affine[[lambda, nu, sigma]] -
    affine[[rho, nu, lambda]] affine[[lambda, mu, sigma]],
    {lambda, dimension}
  ],
  assumptions
];

curvatureTensor[affine_] := Table[
  curvatureComponent[affine, rho, sigma, mu, nu],
  {rho, dimension}, {sigma, dimension}, {mu, dimension}, {nu, dimension}
];

metricCompatibilityResidual[metric_, affine_] := Simplify[
  Table[
    D[metric[[nu, sigma]], coordinates[[mu]]] -
    Sum[
      affine[[rho, mu, nu]] metric[[rho, sigma]] +
      affine[[rho, mu, sigma]] metric[[nu, rho]],
      {rho, dimension}
    ],
    {mu, dimension}, {nu, dimension}, {sigma, dimension}
  ],
  assumptions
];

covariantExteriorDerivative[base_, difference_] := Table[
  Simplify[
    D[difference[[rho, nu, sigma]], coordinates[[mu]]] -
    D[difference[[rho, mu, sigma]], coordinates[[nu]]] +
    Sum[
      base[[rho, mu, lambda]] difference[[lambda, nu, sigma]] +
      difference[[rho, mu, lambda]] base[[lambda, nu, sigma]] -
      base[[rho, nu, lambda]] difference[[lambda, mu, sigma]] -
      difference[[rho, nu, lambda]] base[[lambda, mu, sigma]],
      {lambda, dimension}
    ],
    assumptions
  ],
  {rho, dimension}, {sigma, dimension}, {mu, dimension}, {nu, dimension}
];

connectionWedge[difference_] := Table[
  Simplify[
    Sum[
      difference[[rho, mu, lambda]] difference[[lambda, nu, sigma]] -
      difference[[rho, nu, lambda]] difference[[lambda, mu, sigma]],
      {lambda, dimension}
    ],
    assumptions
  ],
  {rho, dimension}, {sigma, dimension}, {mu, dimension}, {nu, dimension}
];

torsionForm[frame_, spinConnection_] := Table[
  Simplify[
    D[frame[[nu, a]], coordinates[[mu]]] -
    D[frame[[mu, a]], coordinates[[nu]]] +
    Sum[
      spinConnection[[mu, a, b]] frame[[nu, b]] -
      spinConnection[[nu, a, b]] frame[[mu, b]],
      {b, dimension}
    ],
    assumptions
  ],
  {a, dimension}, {mu, dimension}, {nu, dimension}
];

spinCurvature[spinConnection_] := Table[
  Simplify[
    D[spinConnection[[nu, a, b]], coordinates[[mu]]] -
    D[spinConnection[[mu, a, b]], coordinates[[nu]]] +
    Sum[
      spinConnection[[mu, a, c]] spinConnection[[nu, c, b]] -
      spinConnection[[nu, a, c]] spinConnection[[mu, c, b]],
      {c, dimension}
    ],
    assumptions
  ],
  {mu, dimension}, {nu, dimension}, {a, dimension}, {b, dimension}
];

coordinateTriples = Subsets[Range[dimension], {3}];

firstCartanResidual[frame_, spinConnection_, torsion_, curvature_] :=
  Table[
    With[
      {mu = triple[[1]], nu = triple[[2]], rho = triple[[3]]},
      Simplify[
        D[torsion[[a, nu, rho]], coordinates[[mu]]] +
        D[torsion[[a, rho, mu]], coordinates[[nu]]] +
        D[torsion[[a, mu, nu]], coordinates[[rho]]] +
        Sum[
          spinConnection[[mu, a, b]] torsion[[b, nu, rho]] +
          spinConnection[[nu, a, b]] torsion[[b, rho, mu]] +
          spinConnection[[rho, a, b]] torsion[[b, mu, nu]] -
          curvature[[mu, nu, a, b]] frame[[rho, b]] -
          curvature[[nu, rho, a, b]] frame[[mu, b]] -
          curvature[[rho, mu, a, b]] frame[[nu, b]],
          {b, dimension}
        ],
        assumptions
      ]
    ],
    {a, dimension}, {triple, coordinateTriples}
  ];

secondCartanResidual[spinConnection_, curvature_] := Table[
  With[
    {mu = triple[[1]], nu = triple[[2]], rho = triple[[3]]},
    Simplify[
      D[curvature[[nu, rho, a, b]], coordinates[[mu]]] +
      D[curvature[[rho, mu, a, b]], coordinates[[nu]]] +
      D[curvature[[mu, nu, a, b]], coordinates[[rho]]] +
      Sum[
        spinConnection[[mu, a, c]] curvature[[nu, rho, c, b]] -
        curvature[[nu, rho, a, c]] spinConnection[[mu, c, b]] +
        spinConnection[[nu, a, c]] curvature[[rho, mu, c, b]] -
        curvature[[rho, mu, a, c]] spinConnection[[nu, c, b]] +
        spinConnection[[rho, a, c]] curvature[[mu, nu, c, b]] -
        curvature[[mu, nu, a, c]] spinConnection[[rho, c, b]],
        {c, dimension}
      ],
      assumptions
    ]
  ],
  {a, dimension}, {b, dimension}, {triple, coordinateTriples}
];

gammaCanonical = christoffelFromMetric[gCanonical, inverseGCanonical];
omegaCanonicalMixed = spinConnectionFromAffine[
  frameCanonical, inverseFrameCanonical, gammaCanonical
];
omegaCanonical = lowerSpinIndices[omegaCanonicalMixed];
torsionCanonical = torsionTensor[gammaCanonical];

gammaWeitzenboeck = affineFromFrameAndSpin[
  frameCanonical,
  inverseFrameCanonical,
  ConstantArray[0, {dimension, dimension, dimension}]
];
torsionWeitzenboeck = torsionTensor[gammaWeitzenboeck];
curvatureWeitzenboeck = curvatureTensor[gammaWeitzenboeck];

Print["[1] canonical 4+4 geometry: dimension = ", dimension,
  ", tangent signature = ", Diagonal[eta4488],
  ", Det g = ", Simplify[Det[gCanonical], assumptions]];

(* ::NotebookSection:: Split-signature Clifford algebra and type-one spinors *)
(* A real irreducible Cl(4,4) representation on the exterior algebra of R^4.
   Even forms occupy components 1..8 (type 1); odd forms occupy 9..16. *)
exteriorBasis = Join[
  Select[Range[0, 15], EvenQ[DigitCount[#, 2, 1]] &],
  Select[Range[0, 15], OddQ[DigitCount[#, 2, 1]] &]
];
basisPosition[state_] := First[FirstPosition[exteriorBasis, state]];
occupiedQ[state_, index_] := BitAnd[state, 2^(index - 1)] != 0;
orientationSign[state_, index_] := (-1)^Total[
  Boole[occupiedQ[state, #]] & /@ Range[index - 1]
];

creationOperator[index_] := SparseArray[
  Table[
    With[{state = exteriorBasis[[column]]},
      If[
        occupiedQ[state, index],
        Nothing,
        {basisPosition[state + 2^(index - 1)], column} ->
          orientationSign[state, index]
      ]
    ],
    {column, 16}
  ],
  {16, 16}
];

annihilationOperator[index_] := SparseArray[
  Table[
    With[{state = exteriorBasis[[column]]},
      If[
        occupiedQ[state, index],
        {basisPosition[state - 2^(index - 1)], column} ->
          orientationSign[state, index],
        Nothing
      ]
    ],
    {column, 16}
  ],
  {16, 16}
];

gamma4488 = Join[
  Table[Normal[creationOperator[index] + annihilationOperator[index]],
    {index, 4}],
  Table[Normal[creationOperator[index] - annihilationOperator[index]],
    {index, 4}]
];
gammaCommutators = Table[
  gamma4488[[a]] . gamma4488[[b]] -
  gamma4488[[b]] . gamma4488[[a]],
  {a, dimension}, {b, dimension}
];
cliffordResidual = Table[
  gamma4488[[a]] . gamma4488[[b]] +
  gamma4488[[b]] . gamma4488[[a]] -
  2 eta4488[[a, b]] IdentityMatrix[16],
  {a, dimension}, {b, dimension}
];

spinorGaugeConnection[loweredSpinConnection_] := Simplify[
  Table[
    1/8 Sum[
      loweredSpinConnection[[mu, a, b]] gammaCommutators[[a, b]],
      {a, dimension}, {b, dimension}
    ],
    {mu, dimension}
  ],
  assumptions
];

canonicalSpinorGaugeConnection = spinorGaugeConnection[omegaCanonical];

Print["[2] real Clifford algebra: gamma dimensions = ", Dimensions[gamma4488],
  ", all entries integer = ", And @@ (IntegerQ /@ Flatten[gamma4488]),
  ", max |Clifford residual| = ", Max[Abs[Flatten[cliffordResidual]]]];

(* ::NotebookSection:: The gpt-5.6f_rame frame and its spin lift *)
hadamard4 = 1/2 {
  {1, 1, 1, 1},
  {1, -1, 1, -1},
  {1, 1, -1, -1},
  {1, -1, -1, 1}
};
bridgeGenerator = ArrayFlatten[{
  {ConstantArray[0, {4, 4}], hadamard4},
  {Transpose[hadamard4], ConstantArray[0, {4, 4}]}
}];
bridgeRapidity = kappa x0 x4;
bridgeLambda = Simplify[
  Cosh[bridgeRapidity] IdentityMatrix[dimension] +
  Sinh[bridgeRapidity] bridgeGenerator,
  assumptions
];
inverseBridgeLambda = Simplify[
  Cosh[bridgeRapidity] IdentityMatrix[dimension] -
  Sinh[bridgeRapidity] bridgeGenerator,
  assumptions
];

rotatedNegativeGamma = Table[
  Sum[hadamard4[[index, j]] gamma4488[[4 + j]], {j, 4}],
  {index, 4}
];
bridgeSpinFactors = Table[
  Cosh[bridgeRapidity/2] IdentityMatrix[16] -
  Sinh[bridgeRapidity/2]
    gamma4488[[index]] . rotatedNegativeGamma[[index]],
  {index, 4}
];
inverseBridgeSpinFactors = Table[
  Cosh[bridgeRapidity/2] IdentityMatrix[16] +
  Sinh[bridgeRapidity/2]
    gamma4488[[index]] . rotatedNegativeGamma[[index]],
  {index, 4}
];
bridgeSpinLift = Simplify[Dot @@ bridgeSpinFactors, assumptions];
inverseBridgeSpinLift = Simplify[
  Dot @@ Reverse[inverseBridgeSpinFactors], assumptions
];
spinLiftCliffordResidual = Simplify[
  Table[
    bridgeSpinLift . gamma4488[[a]] . inverseBridgeSpinLift -
    Sum[bridgeLambda[[a, b]] gamma4488[[b]], {b, dimension}],
    {a, dimension}
  ],
  assumptions
];

(* Exact visible label: gpt-5.6f_rame. Valid Wolfram symbol: gpt56Frame. *)
gpt56Frame = Simplify[frameCanonical . bridgeLambda, assumptions];
inverseGpt56Frame = Simplify[
  inverseBridgeLambda . inverseFrameCanonical, assumptions
];

omegaLeviCivitaInBridgeFrame = spinConnectionFromAffine[
  gpt56Frame, inverseGpt56Frame, gammaCanonical
];
omegaWeitzenboeckInBridgeFrame = spinConnectionFromAffine[
  gpt56Frame, inverseGpt56Frame, gammaWeitzenboeck
];

Print["[3] gpt-5.6f_rame: Hadamard generator squares to the identity = ",
  bridgeGenerator . bridgeGenerator === IdentityMatrix[dimension],
  ", Det Lambda = ", Simplify[Det[bridgeLambda], assumptions],
  ", rapidity = ", bridgeRapidity];

(* ::NotebookSection:: The gpt-5.6_bridge connection family *)
(* Exact visible label: gpt-5.6_bridge. Valid Wolfram symbol: gpt56Bridge. *)
gpt56Contortion = Simplify[
  omegaWeitzenboeckInBridgeFrame - omegaLeviCivitaInBridgeFrame,
  assumptions
];
gpt56BridgeGeneral = Simplify[
  omegaLeviCivitaInBridgeFrame +
  t gpt56Contortion,
  assumptions
];
gpt56Bridge[parameter_] := gpt56BridgeGeneral /. t -> parameter;

gpt56BridgeAffineGeneral = Simplify[
  gammaCanonical + t (gammaWeitzenboeck - gammaCanonical),
  assumptions
];
gpt56BridgeAffine[parameter_] :=
  gpt56BridgeAffineGeneral /. t -> parameter;

affineContortionGPT = Simplify[
  gammaWeitzenboeck - gammaCanonical, assumptions
];
curvatureCanonical = curvatureTensor[gammaCanonical];
covariantContortionGPT = covariantExteriorDerivative[
  gammaCanonical, affineContortionGPT
];
contortionWedgeGPT = connectionWedge[affineContortionGPT];
gpt56BridgeCurvatureGeneral = curvatureTensor[gpt56BridgeAffineGeneral];
gpt56BridgeCurvature[parameter_] :=
  gpt56BridgeCurvatureGeneral /. t -> parameter;
curvatureTransgressionGPT = Simplify[
  curvatureCanonical + t covariantContortionGPT +
  t^2 contortionWedgeGPT,
  assumptions
];

gpt56BridgeTorsionGeneral = torsionTensor[gpt56BridgeAffineGeneral];
gpt56BridgeTorsion[parameter_] :=
  gpt56BridgeTorsionGeneral /. t -> parameter;

gpt56BridgeSpinorConnectionGeneral = spinorGaugeConnection[
  lowerSpinIndices[gpt56BridgeGeneral]
];
gpt56BridgeSpinorConnection[parameter_] :=
  gpt56BridgeSpinorConnectionGeneral /. t -> parameter;

omegaGPT56 = gpt56Bridge[1/2];

spinorConnectionLeviCivitaGPT = spinorGaugeConnection[
  lowerSpinIndices[omegaLeviCivitaInBridgeFrame]
];
spinorConnectionWeitzenboeckGPT = spinorGaugeConnection[
  lowerSpinIndices[omegaWeitzenboeckInBridgeFrame]
];
spinorConnectionGaugeTransform[connection_] := Simplify[
  Table[
    inverseBridgeSpinLift . connection[[mu]] . bridgeSpinLift +
    inverseBridgeSpinLift . D[bridgeSpinLift, coordinates[[mu]]],
    {mu, dimension}
  ],
  assumptions
];

gammaCurvedCanonical = Simplify[
  Table[
    Sum[
      inverseFrameCanonical[[a, mu]] gamma4488[[a]],
      {a, dimension}
    ],
    {mu, dimension}
  ],
  assumptions
];
gammaCurvedGPT = Simplify[
  Table[
    Sum[inverseGpt56Frame[[a, mu]] gamma4488[[a]], {a, dimension}],
    {mu, dimension}
  ],
  assumptions
];

omegaBridgeCanonicalGeneral = Simplify[
  (1 - t) omegaCanonicalMixed, assumptions
];
spinorConnectionBridgeCanonical = spinorGaugeConnection[
  lowerSpinIndices[omegaBridgeCanonicalGeneral]
];
torsionFormBridgeCanonical = torsionForm[
  frameCanonical, omegaBridgeCanonicalGeneral
];
spinCurvatureBridgeCanonical = spinCurvature[
  omegaBridgeCanonicalGeneral
];
firstCartanBridgeResidual = firstCartanResidual[
  frameCanonical,
  omegaBridgeCanonicalGeneral,
  torsionFormBridgeCanonical,
  spinCurvatureBridgeCanonical
];
secondCartanBridgeResidual = secondCartanResidual[
  omegaBridgeCanonicalGeneral, spinCurvatureBridgeCanonical
];
torsionPotentialGPT = Log[Sin[6 H x0]];
diracConnectionCanonical = Simplify[
  Sum[
    gammaCurvedCanonical[[mu]] . canonicalSpinorGaugeConnection[[mu]],
    {mu, dimension}
  ],
  assumptions
];
diracBridgeConjugacyResidual = Simplify[
  Sum[
    gammaCurvedCanonical[[mu]] . (
      canonicalSpinorGaugeConnection[[mu]] -
      spinorConnectionBridgeCanonical[[mu]] +
      (t/2) D[torsionPotentialGPT, coordinates[[mu]]]
        IdentityMatrix[16]
    ),
    {mu, dimension}
  ],
  assumptions
];

interiorTorsionWitness = FullSimplify[
  gpt56BridgeTorsion[1/2][[2, 1, 2]],
  assumptions
];
interiorCurvatureWitness = FullSimplify[
  curvatureComponent[gpt56BridgeAffine[1/2], 2, 1, 1, 2],
  assumptions
];

Print["[4] gpt-5.6_bridge family: interior torsion witness T^1_(0 1) = ",
  interiorTorsionWitness,
  ", interior curvature witness R^1_(0 0 1) = ", interiorCurvatureWitness];

(* ::NotebookSection:: Exact verification report *)
tests = {
  VerificationTest[Dimensions[frameCanonical], {8, 8},
    TestID -> "canonical-frame-dimensions"],
  VerificationTest[
    zeroArrayQ[
      gCanonical - frameCanonical . eta4488 . Transpose[frameCanonical]
    ],
    True,
    TestID -> "canonical-frame-reconstructs-metric"
  ],
  VerificationTest[
    zeroArrayQ[
      vielbeinPostulateResidual[
        frameCanonical, gammaCanonical, omegaCanonicalMixed
      ]
    ],
    True,
    TestID -> "canonical-vielbein-postulate"
  ],
  VerificationTest[zeroArrayQ[torsionCanonical], True,
    TestID -> "canonical-torsion-free"],
  VerificationTest[
    zeroArrayQ[omegaCanonical + Transpose[omegaCanonical, {1, 3, 2}]],
    True,
    TestID -> "canonical-spin-connection-antisymmetric"
  ],
  VerificationTest[nonzeroCount[omegaCanonical], 24,
    TestID -> "canonical-spin-connection-has-24-components"],
  VerificationTest[zeroArrayQ[curvatureWeitzenboeck], True,
    TestID -> "weitzenboeck-curvature-free"],
  VerificationTest[zeroArrayQ[cliffordResidual], True,
    TestID -> "clifford-4-plus-4"],
  VerificationTest[
    zeroArrayQ[gammaCommutators[[All, All, 1 ;; 8, 9 ;; 16]]] &&
    zeroArrayQ[gammaCommutators[[All, All, 9 ;; 16, 1 ;; 8]]],
    True,
    TestID -> "type-one-spinors-preserved"
  ],
  VerificationTest[
    zeroArrayQ[
      Transpose[bridgeGenerator] . eta4488 + eta4488 . bridgeGenerator
    ],
    True,
    TestID -> "bridge-generator-is-so-4-4"
  ],
  VerificationTest[
    zeroArrayQ[
      bridgeLambda . eta4488 . Transpose[bridgeLambda] - eta4488
    ],
    True,
    TestID -> "bridge-lambda-is-pseudo-orthogonal"
  ],
  VerificationTest[
    zeroArrayQ[
      bridgeSpinLift . inverseBridgeSpinLift - IdentityMatrix[16]
    ],
    True,
    TestID -> "bridge-spin-lift-inverse"
  ],
  VerificationTest[zeroArrayQ[spinLiftCliffordResidual], True,
    TestID -> "bridge-spin-lift-implements-lambda"],
  VerificationTest[
    zeroArrayQ[bridgeSpinLift[[1 ;; 8, 9 ;; 16]]] &&
    zeroArrayQ[bridgeSpinLift[[9 ;; 16, 1 ;; 8]]],
    True,
    TestID -> "bridge-spin-lift-preserves-type-one"
  ],
  VerificationTest[
    zeroArrayQ[
      gammaCurvedGPT - Table[
        inverseBridgeSpinLift . gammaCurvedCanonical[[mu]] .
          bridgeSpinLift,
        {mu, dimension}
      ]
    ],
    True,
    TestID -> "curved-gamma-frame-covariance"
  ],
  VerificationTest[
    zeroArrayQ[
      gpt56Frame . eta4488 . Transpose[gpt56Frame] - gCanonical
    ],
    True,
    TestID -> "gpt56-frame-preserves-metric"
  ],
  VerificationTest[
    zeroArrayQ[
      vielbeinPostulateResidual[
        gpt56Frame, gpt56BridgeAffine[t], gpt56Bridge[t]
      ]
    ],
    True,
    TestID -> "bridge-vielbein-postulate"
  ],
  VerificationTest[
    zeroArrayQ[
      lowerSpinIndices[gpt56Bridge[t]] +
      Transpose[lowerSpinIndices[gpt56Bridge[t]], {1, 3, 2}]
    ],
    True,
    TestID -> "bridge-spin-connection-antisymmetric"
  ],
  VerificationTest[
    zeroArrayQ[
      metricCompatibilityResidual[gCanonical, gpt56BridgeAffine[t]]
    ],
    True,
    TestID -> "bridge-metric-compatible"
  ],
  VerificationTest[
    zeroArrayQ[gpt56BridgeAffine[0] - gammaCanonical],
    True,
    TestID -> "bridge-canonical-endpoint"
  ],
  VerificationTest[
    zeroArrayQ[gpt56BridgeAffine[1] - gammaWeitzenboeck],
    True,
    TestID -> "bridge-weitzenboeck-endpoint"
  ],
  VerificationTest[
    zeroArrayQ[
      gpt56BridgeTorsion[t] - t torsionWeitzenboeck
    ],
    True,
    TestID -> "bridge-torsion-transgression"
  ],
  VerificationTest[
    zeroArrayQ[
      gpt56BridgeCurvature[t] - curvatureTransgressionGPT
    ],
    True,
    TestID -> "bridge-curvature-transgression"
  ],
  VerificationTest[
    zeroArrayQ[gpt56BridgeCurvature[0] - curvatureCanonical],
    True,
    TestID -> "bridge-curvature-canonical-endpoint"
  ],
  VerificationTest[zeroArrayQ[gpt56BridgeCurvature[1]], True,
    TestID -> "bridge-curvature-weitzenboeck-endpoint"],
  VerificationTest[zeroQ[interiorTorsionWitness], False,
    TestID -> "bridge-interior-torsion-nonzero"],
  VerificationTest[zeroQ[interiorCurvatureWitness], False,
    TestID -> "bridge-interior-curvature-nonzero"],
  VerificationTest[
    zeroArrayQ[
      gpt56BridgeSpinorConnection[t][[All, 1 ;; 8, 9 ;; 16]]
    ] &&
    zeroArrayQ[
      gpt56BridgeSpinorConnection[t][[All, 9 ;; 16, 1 ;; 8]]
    ],
    True,
    TestID -> "bridge-preserves-type-one-spinors"
  ],
  VerificationTest[
    zeroArrayQ[
      spinorConnectionLeviCivitaGPT -
      spinorConnectionGaugeTransform[canonicalSpinorGaugeConnection]
    ],
    True,
    TestID -> "levi-civita-spinor-gauge-law"
  ],
  VerificationTest[
    zeroArrayQ[
      spinorConnectionWeitzenboeckGPT -
      spinorConnectionGaugeTransform[
        ConstantArray[0, {dimension, 16, 16}]
      ]
    ],
    True,
    TestID -> "weitzenboeck-spinor-pure-gauge"
  ],
  VerificationTest[
    zeroArrayQ[
      gpt56BridgeSpinorConnection[t] -
      spinorConnectionGaugeTransform[spinorConnectionBridgeCanonical]
    ],
    True,
    TestID -> "bridge-spinor-gauge-law"
  ],
  VerificationTest[
    zeroArrayQ[
      diracConnectionCanonical +
      (1/2) Sum[
        D[torsionPotentialGPT, coordinates[[mu]]]
          gammaCurvedCanonical[[mu]],
        {mu, dimension}
      ]
    ],
    True,
    TestID -> "canonical-dirac-connection-is-torsion-gradient"
  ],
  VerificationTest[zeroArrayQ[diracBridgeConjugacyResidual], True,
    TestID -> "bridge-dirac-scalar-conjugacy"
  ],
  VerificationTest[
    zeroArrayQ[
      torsionFormBridgeCanonical - Table[
        Sum[
          frameCanonical[[rho, a]]
            gpt56BridgeTorsion[t][[rho, mu, nu]],
          {rho, dimension}
        ],
        {a, dimension}, {mu, dimension}, {nu, dimension}
      ]
    ],
    True,
    TestID -> "bridge-frame-and-affine-torsion-agree"
  ],
  VerificationTest[zeroArrayQ[firstCartanBridgeResidual], True,
    TestID -> "bridge-first-cartan-bianchi-identity"],
  VerificationTest[zeroArrayQ[secondCartanBridgeResidual], True,
    TestID -> "bridge-second-cartan-bianchi-identity"
  ]
};

report = TestReport[tests];
(* The raw TestReportObject carries timestamps and a CreationID, so printing it
   would make every captured notebook output differ between runs. *)
Print["Tests run: ", Length[tests]];
Print["Canonical nonzero spin-connection components: ",
  nonzeroCount[omegaCanonical]];
Print["Interior torsion witness T^1_(0 1): ", interiorTorsionWitness];
Print["Interior curvature witness R^1_(0 0 1): ",
  interiorCurvatureWitness];
Print["Tests succeeded: ", report["TestsSucceededCount"]];
Print["Tests failed: ", report["TestsFailedCount"]];

If[report["TestsFailedCount"] > 0 && $FrontEnd === Null, Exit[1]];