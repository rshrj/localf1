(* ::Package:: *)

(*
  SWModularWorkbench.wl

  Lightweight Mathematica/Wolfram Language utilities for plotting modular
  fundamental domains, congruence-subgroup domains, cusps, and images under
  Seiberg-Witten/modular maps.

  Author: ChatGPT for Rishi Raj
  Version: 0.1.0

  Design notes:
  - Matrices are 2 x 2 integer representatives of PSL(2,Z).
  - Gamma and -Gamma are identified by PSLNormalize.
  - The standard SL(2,Z) tile is represented by a polygonal approximation
    to the usual boundary |Re[tau]| <= 1/2, |tau| >= 1.
  - For plotting, cusps at Infinity are recorded symbolically but not drawn.
*)

BeginPackage["SWModularWorkbench`"];

(* ---------------------------------------------------------------------- *)
(* Public usage messages                                                   *)
(* ---------------------------------------------------------------------- *)

SMatrix::usage = "SMatrix is the standard modular generator {{0,-1},{1,0}}.";
TMatrix::usage = "TMatrix[n] is the modular translation matrix {{1,n},{0,1}}. TMatrix[] is TMatrix[1].";
IdMatrix::usage = "IdMatrix is the 2 x 2 identity matrix.";

Mobius::usage = "Mobius[gamma][tau] applies gamma={{a,b},{c,d}} to tau by (a tau+b)/(c tau+d).";
TransformPoint::usage = "TransformPoint[gamma,z] applies gamma to z, allowing z === Infinity.";
TransformBoundary::usage = "TransformBoundary[gamma,pts] applies gamma to a list of complex boundary points.";
CPoint::usage = "CPoint[z] returns {Re[z],Im[z]} for a finite complex number z.";

PSLNormalize::usage = "PSLNormalize[gamma] chooses a canonical representative of {gamma,-gamma}.";
PSLEqualQ::usage = "PSLEqualQ[gamma1,gamma2] tests equality in PSL(2,Z), i.e. up to sign.";
SL2ZMatrixQ::usage = "SL2ZMatrixQ[gamma] tests whether gamma is a 2 x 2 integer matrix of determinant 1.";

FundamentalBoundary::usage = "FundamentalBoundary[ymax,n] returns complex points approximating the boundary of the standard SL(2,Z) fundamental tile truncated at Im[tau]=ymax.";
FundamentalPolygon::usage = "FundamentalPolygon[ymax,n] returns a Polygon object for the standard tile.";

PlotDomains::usage = "PlotDomains[gammas,opts] plots the images gamma.F of the standard SL(2,Z) tile.";
PlotModularDomain::usage = "PlotModularDomain[gammas,opts] plots a union of modular tiles, optionally marks finite cusps, and optionally draws the full SL(2,Z) tessellation as a background (ShowBackground->True).";
PlotMappedDomain::usage = "PlotMappedDomain[f,gammas,opts] plots the image of tile boundaries under f after modular transformations.";

Gamma0Q::usage = "Gamma0Q[N][gamma] tests whether gamma lies in Gamma0(N).";
Gamma1Q::usage = "Gamma1Q[N][gamma] tests whether gamma lies in Gamma1(N).";
GammaNQ::usage = "GammaNQ[N][gamma] tests whether gamma lies in the principal congruence subgroup Gamma(N).";

IndexGamma0::usage = "IndexGamma0[N] gives [SL(2,Z):Gamma0(N)].";
IndexGammaN::usage = "IndexGammaN[N] gives [SL(2,Z):Gamma(N)] for N>=1, with the standard PSL convention for N=1,2 handled separately.";

SL2ZWords::usage = "SL2ZWords[maxLen] returns PSL-distinct words in S,T,T^-1 up to length maxLen.";
CosetReps::usage = "CosetReps[subgroupQ,targetIndex,maxLen] finds left coset representatives gamma such that gamma.F tiles a subgroup fundamental domain. subgroupQ is a predicate on matrices.";
CosetRepsGamma0::usage = "CosetRepsGamma0[N,maxLen] returns representatives for Gamma0(N)\\SL(2,Z), suitable for plotting a Gamma0(N) domain as union of gamma.F.";
CosetRepsGamma1::usage = "CosetRepsGamma1[N,maxLen] returns representatives for Gamma1(N)\\SL(2,Z), using a brute-force search.";
CosetRepsGammaN::usage = "CosetRepsGammaN[N,maxLen] returns representatives for Gamma(N)\\SL(2,Z), using a brute-force search.";
SameLeftCosetQ::usage = "SameLeftCosetQ[subgroupQ][g1,g2] tests whether g1 and g2 are in the same left coset of a subgroup, i.e. g1.g2^-1 lies in the subgroup.";

CuspsFromReps::usage = "CuspsFromReps[reps] gives DeleteDuplicates[gamma.Infinity] for representatives reps.";
FiniteCuspsFromReps::usage = "FiniteCuspsFromReps[reps] gives finite cusps from CuspsFromReps[reps].";

LambdaModular::usage = "LambdaModular[tau] is the modular lambda function theta2[tau]^4/theta3[tau]^4.";
Theta2C::usage = "Theta2C[tau] is theta_2(tau)=EllipticTheta[2,0,Exp[Pi I tau]].";
Theta3C::usage = "Theta3C[tau] is theta_3(tau)=EllipticTheta[3,0,Exp[Pi I tau]].";
Theta4C::usage = "Theta4C[tau] is theta_4(tau)=EllipticTheta[4,0,Exp[Pi I tau]].";
JPureSU2::usage = "JPureSU2[tau] is currently an alias for UPureSU2[tau]. Check conventions for your paper before using it physically.";
UPureSU2::usage = "UPureSU2[tau] is a common normalized pure SU(2)-style Hauptmodul: (theta3^4+theta4^4)/(2 theta2^4). Check conventions for your paper.";

KnownCosetReps::usage = "KnownCosetReps[name] returns small hard-coded representatives for names like \"SL2Z\", \"Gamma0[2]\", \"Gamma0[4]\", \"Gamma[2]\".";

Begin["`Private`"];

(* ---------------------------------------------------------------------- *)
(* Basic matrices and PSL(2,Z) action                                      *)
(* ---------------------------------------------------------------------- *)

IdMatrix = IdentityMatrix[2];
SMatrix = {{0, -1}, {1, 0}};
TMatrix[n_Integer : 1] := {{1, n}, {0, 1}};

SL2ZMatrixQ[gamma_] := MatrixQ[gamma, IntegerQ] && Dimensions[gamma] === {2, 2} && Det[gamma] === 1;

FirstNonzeroEntry[m_] := Module[{v = Select[Flatten[m], # =!= 0 &]}, If[v === {}, 0, First[v]]];

PSLNormalize[gamma_?MatrixQ] := Module[{g = gamma}, If[FirstNonzeroEntry[g] < 0, -g, g]];

PSLEqualQ[g1_?MatrixQ, g2_?MatrixQ] := PSLNormalize[g1] === PSLNormalize[g2];

Mobius[gamma_?MatrixQ][tau_] := Module[{a, b, c, d},
  {{a, b}, {c, d}} = gamma;
  (a tau + b)/(c tau + d)
];

TransformPoint[gamma_?MatrixQ, z_] := Module[{a, b, c, d},
  {{a, b}, {c, d}} = gamma;
  If[z === Infinity,
    If[c === 0, Infinity, a/c],
    (a z + b)/(c z + d)
  ]
];

TransformBoundary[gamma_?MatrixQ, pts_List] := TransformPoint[gamma, #] & /@ pts;

CPoint[z_] /; z =!= Infinity := {Re[N[z]], Im[N[z]]};

(* ---------------------------------------------------------------------- *)
(* Fundamental tile                                                        *)
(* ---------------------------------------------------------------------- *)

FundamentalBoundary[ymax_: 3, n_Integer : 200] /; ymax > Sqrt[3]/2 && n >= 8 := Module[
  {nSide, nArc, left, arc, right},
  nSide = Max[4, Floor[n/3]];
  nArc = Max[8, n - 2 nSide];
  left = Table[-1/2 + I y, {y, Sqrt[3]/2, ymax, (ymax - Sqrt[3]/2)/nSide}];
  arc = Table[Exp[I theta], {theta, 2 Pi/3, Pi/3, -(Pi/3)/nArc}];
  right = Table[1/2 + I y, {y, Sqrt[3]/2, ymax, (ymax - Sqrt[3]/2)/nSide}];
  Join[left, Reverse[right], Reverse[arc]]
];

FundamentalPolygon[ymax_: 3, n_Integer : 200] := Polygon[CPoint /@ FundamentalBoundary[ymax, n]];

(* ---------------------------------------------------------------------- *)
(* Congruence subgroup tests                                               *)
(* ---------------------------------------------------------------------- *)

Gamma0Q[N_Integer?Positive][gamma_?MatrixQ] := SL2ZMatrixQ[gamma] && Mod[gamma[[2, 1]], N] === 0;

Gamma1Q[N_Integer?Positive][gamma_?MatrixQ] := Module[{a, b, c, d},
  If[! SL2ZMatrixQ[gamma], Return[False]];
  {{a, b}, {c, d}} = gamma;
  Mod[c, N] === 0 && Mod[a - 1, N] === 0 && Mod[d - 1, N] === 0
];

GammaNQ[N_Integer?Positive][gamma_?MatrixQ] := SL2ZMatrixQ[gamma] && Mod[gamma - IdMatrix, N] === ConstantArray[0, {2, 2}];

IndexGamma0[N_Integer?Positive] := N Times @@ ((1 + 1/#) & /@ FactorInteger[N][[All, 1]]);

IndexGammaN[1] := 1;
IndexGammaN[2] := 6;
IndexGammaN[N_Integer?Positive] /; N >= 3 := (N^3/2) Times @@ ((1 - 1/#^2) & /@ FactorInteger[N][[All, 1]]);

(* ---------------------------------------------------------------------- *)
(* Word generation and cosets                                              *)
(* ---------------------------------------------------------------------- *)

DeletePSLDuplicates[list_List] := DeleteDuplicatesBy[list, PSLNormalize];

SL2ZWords[maxLen_Integer?NonNegative] := Module[{gens, layers, current, all},
  gens = {SMatrix, TMatrix[1], TMatrix[-1]};
  current = {IdMatrix};
  all = current;
  Do[
    current = DeletePSLDuplicates[Flatten[Outer[Dot, current, gens, 1], 1]];
    all = DeletePSLDuplicates[Join[all, current]],
    {maxLen}
  ];
  all
];

SameLeftCosetQ[subgroupQ_][g1_?MatrixQ, g2_?MatrixQ] := Module[{delta},
  (* g1 and g2 are in the same left coset Gamma\SL2Z iff g1.g2^-1 in Gamma. *)
  delta = g1 . Inverse[g2];
  MatrixQ[delta, IntegerQ] && Det[delta] === 1 && TrueQ[subgroupQ[delta]]
];

CosetReps[subgroupQ_, targetIndex_Integer?Positive, maxLen_Integer : 8] := Module[
  {candidates, reps = {}},
  candidates = SL2ZWords[maxLen];
  Do[
    If[! AnyTrue[reps, SameLeftCosetQ[subgroupQ][gamma, #] &], AppendTo[reps, PSLNormalize[gamma]]];
    If[Length[reps] >= targetIndex, Break[]],
    {gamma, candidates}
  ];
  If[Length[reps] < targetIndex,
    Message[CosetReps::few, Length[reps], targetIndex, maxLen]
  ];
  reps
];

CosetReps::few = "Found only `1` representatives, expected `2`. Increase maxLen beyond `3`.";

CosetRepsGamma0[N_Integer?Positive, maxLen_Integer : 8] := CosetReps[Gamma0Q[N], IndexGamma0[N], maxLen];

(* Geometric/PSL index for Gamma1(N). The SL index is N^2 prod(1-1/p^2); divide by 2 for N>=3 because -I is not in Gamma1(N) but acts trivially. *)
IndexGamma1[1] := 1;
IndexGamma1[N_Integer?Positive] /; N >= 2 := Module[{slIndex},
  slIndex = N^2 Times @@ ((1 - 1/#^2) & /@ FactorInteger[N][[All, 1]]);
  If[N <= 2, slIndex, slIndex/2]
];

CosetRepsGamma1[N_Integer?Positive, maxLen_Integer : 9] := CosetReps[Gamma1Q[N], IndexGamma1[N], maxLen];
CosetRepsGammaN[N_Integer?Positive, maxLen_Integer : 10] := CosetReps[GammaNQ[N], IndexGammaN[N], maxLen];

KnownCosetReps["SL2Z"] := {IdMatrix};
KnownCosetReps["Gamma0[2]"] := {IdMatrix, SMatrix, SMatrix . TMatrix[1]};
KnownCosetReps["Gamma[2]"] := {IdMatrix, TMatrix[1], SMatrix, SMatrix . TMatrix[1], TMatrix[1] . SMatrix, SMatrix . TMatrix[1] . SMatrix};
KnownCosetReps["Gamma0[4]"] := CosetRepsGamma0[4, 8];
KnownCosetReps[name_String] := (Message[KnownCosetReps::unk, name]; $Failed);
KnownCosetReps::unk = "Unknown coset-representative set `1`.";

(* ---------------------------------------------------------------------- *)
(* Cusps                                                                   *)
(* ---------------------------------------------------------------------- *)

NormalizeRationalCusp[x_Rational] := x;
NormalizeRationalCusp[x_Integer] := x;
NormalizeRationalCusp[x_] := x;

CuspsFromReps[reps_List] := DeleteDuplicates[NormalizeRationalCusp /@ (TransformPoint[#, Infinity] & /@ reps)];
FiniteCuspsFromReps[reps_List] := DeleteCases[CuspsFromReps[reps], Infinity];

(* ---------------------------------------------------------------------- *)
(* Plotting                                                                *)
(* ---------------------------------------------------------------------- *)

Options[PlotDomains] = {
  YMax -> 3,
  BoundaryPoints -> 220,
  DomainOpacity -> 0.35,
  PlotRange -> Automatic,
  Frame -> True,
  Axes -> False,
  AspectRatio -> 1,
  ImageSize -> Large
};

PlotDomains[gammas_List, OptionsPattern[]] := Module[
  {ymax, n, base, transformed, pr},
  ymax = OptionValue[YMax];
  n = OptionValue[BoundaryPoints];
  base = FundamentalBoundary[ymax, n];
  transformed = TransformBoundary[#, base] & /@ gammas;
  pr = Replace[OptionValue[PlotRange], Automatic :> {{-3, 3}, {0, ymax}}];
  Graphics[
    {
      Opacity[OptionValue[DomainOpacity]],
      Table[
        {ColorData[97][1 + Mod[i - 1, 97]], EdgeForm[Directive[Black, Thin]], Polygon[CPoint /@ DeleteCases[transformed[[i]], Infinity]]},
        {i, Length[transformed]}
      ]
    },
    Frame -> OptionValue[Frame],
    Axes -> OptionValue[Axes],
    AspectRatio -> OptionValue[AspectRatio],
    PlotRange -> pr,
    ImageSize -> OptionValue[ImageSize]
  ]
];

Options[PlotModularDomain] = Join[
  Options[PlotDomains],
  {
    ShowCusps -> True,
    CuspStyle -> {Red, PointSize[Large]},
    CuspLabelStyle -> {Black, 12},
    CuspLabelOffset -> {0.08, 0.08},
    ShowBackground -> False,
    BackgroundStyle -> {GrayLevel[0.9], EdgeForm[Directive[GrayLevel[0.5], Thin]]},
    BackgroundMaxLen -> 6
  }
];

PlotModularDomain[gammas_List, OptionsPattern[]] := Module[
  {domain, finiteCusps, cuspGraphics, pr, ymax, n},
  ymax = OptionValue[YMax];
  n    = OptionValue[BoundaryPoints];
  pr   = Replace[OptionValue[PlotRange], Automatic :> {{-3, 3}, {0, ymax}}];
  domain = PlotDomains[
    gammas,
    YMax -> ymax, BoundaryPoints -> n,
    DomainOpacity -> OptionValue[DomainOpacity],
    PlotRange -> pr,
    Frame -> OptionValue[Frame], Axes -> OptionValue[Axes],
    AspectRatio -> OptionValue[AspectRatio], ImageSize -> OptionValue[ImageSize]
  ];
  finiteCusps = FiniteCuspsFromReps[gammas];
  cuspGraphics = If[TrueQ[OptionValue[ShowCusps]] && finiteCusps =!= {},
    Graphics[{
      OptionValue[CuspStyle], Point[CPoint /@ finiteCusps],
      Table[Text[Style[TraditionalForm[finiteCusps[[i]]], Sequence @@ OptionValue[CuspLabelStyle]],
                 CPoint[finiteCusps[[i]]] + OptionValue[CuspLabelOffset]],
            {i, Length[finiteCusps]}]
    }],
    Graphics[{}]
  ];
  If[TrueQ[OptionValue[ShowBackground]],
    Module[{x1, x2, y1, y2, base30, cands, bgReps, bgBase, bgTiles},
      {{x1, x2}, {y1, y2}} = pr;
      base30  = FundamentalBoundary[ymax, 30];
      cands   = SL2ZWords[OptionValue[BackgroundMaxLen]];
      bgReps  = Select[cands, Function[g,
        AnyTrue[
          CPoint /@ DeleteCases[TransformBoundary[g, base30], Infinity],
          x1 <= #[[1]] <= x2 && y1 <= #[[2]] <= y2 &
        ]
      ]];
      bgBase  = FundamentalBoundary[ymax, n];
      bgTiles = TransformBoundary[#, bgBase] & /@ bgReps;
      Show[
        Graphics[
          {Sequence @@ OptionValue[BackgroundStyle],
           Table[Polygon[CPoint /@ DeleteCases[bgTiles[[i]], Infinity]], {i, Length[bgTiles]}]},
          Frame -> OptionValue[Frame], Axes -> OptionValue[Axes],
          AspectRatio -> OptionValue[AspectRatio], PlotRange -> pr,
          ImageSize -> OptionValue[ImageSize]
        ],
        domain, cuspGraphics
      ]
    ],
    Show[domain, cuspGraphics]
  ]
];

Options[PlotMappedDomain] = {
  YMax -> 3,
  BoundaryPoints -> 220,
  PlotRange -> All,
  Frame -> True,
  Axes -> False,
  AspectRatio -> 1,
  ImageSize -> Large,
  BoundaryStyle -> Thick
};

PlotMappedDomain[f_, gammas_List, OptionsPattern[]] := Module[
  {ymax, n, base, curves},
  ymax = OptionValue[YMax];
  n = OptionValue[BoundaryPoints];
  base = FundamentalBoundary[ymax, n];
  curves = Table[
    Quiet[Check[f /@ TransformBoundary[gammas[[i]], base], {}]],
    {i, Length[gammas]}
  ];
  Graphics[
    Table[
      {ColorData[97][1 + Mod[i - 1, 97]], OptionValue[BoundaryStyle], Line[CPoint /@ DeleteCases[curves[[i]], Infinity]]},
      {i, Length[curves]}
    ],
    Frame -> OptionValue[Frame],
    Axes -> OptionValue[Axes],
    AspectRatio -> OptionValue[AspectRatio],
    PlotRange -> OptionValue[PlotRange],
    ImageSize -> OptionValue[ImageSize]
  ]
];

(* ---------------------------------------------------------------------- *)
(* Theta/modular functions useful in SW notebooks                          *)
(* ---------------------------------------------------------------------- *)

Theta2C[tau_] := EllipticTheta[2, 0, Exp[Pi I tau]];
Theta3C[tau_] := EllipticTheta[3, 0, Exp[Pi I tau]];
Theta4C[tau_] := EllipticTheta[4, 0, Exp[Pi I tau]];

LambdaModular[tau_] := Theta2C[tau]^4/Theta3C[tau]^4;

UPureSU2[tau_] := (Theta3C[tau]^4 + Theta4C[tau]^4)/(2 Theta2C[tau]^4);
JPureSU2[tau_] := UPureSU2[tau];

End[];

EndPackage[];
