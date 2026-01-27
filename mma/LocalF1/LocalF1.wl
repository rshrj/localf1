(* ::Package:: *)

(* :Title:    LocalF1 Package
   :Context:  LocalF1`
   :Author:   Rishi Raj
   :Summary:  Utilities for local F1 geometry, periods, and BPS indices.
   :Package Version: 0.1.0
   :Mathematica Version: 14.3
   :History:
     0.1.0 - 2025-12-01 - Initial version.
*)



(* ::Package:: *)
(**)


BeginPackage["LocalF1`"];


(* ===================== *)
(*   Public symbols      *)
(* ===================== *)

(* Usage messages for exported functions *)

RichardsonResum::usage = "TODO";

MirrorCurve::usage = "MirrorCurve[x, y, m, u] is the affine mirror-curve equation of local F1 \
in the variables x, y and parameters m, u.";

MirrorCurvef3::usage = "TODO";

MirrorCurveg2::usage = "TODO";

MirrorCurveg3::usage = "TODO";

MirrorCurveDelta::usage = "TODO";

MirrorCurveJ::usage = "TODO";

URoot::usage = "URoot[la, k] is the kth root of the order 4 part of the discriminant.";

URoots::usage = "TODO";

PlotURoots::usage = "TODO";

NodalCurveParam::usage = "TODO";

NodalCurveX0::usage = "TODO";

NodalCurveX0tot1::usage = "TODO";

WeierstrassChangeOfVariables::usage = "TODO";

WeierstrassChangeOfVariablesInverse::usage = "TODO";

LiEval::usage = "Replacement rule that numerically evaluates expressions like Li2[x], BW[x] etc";

RationalData::usage = "RationalData[expr, var] gives the roots and exponents of rational expression expr in var.";

NodalCurvebFromt1::usage = "TODO";

NodalCurvebFrommu::usage = "TODO";

KerrDoranP::usage = "TODO";

KerrDoranQ::usage = "TODO";

BW::usage = "BW[z] Bloch Wigner Dilog";

Li2::usage = "Ordinary dilogarithm";


(* Large volume stuff *)

Ch::usage = "Ch[p, q] represents the Chern vector of the sheaf O(p F + q B).";
GV::usage = "Ch[p, q] represents the Chern vector of the sheaf O(p F + q B).";
y::usage = "Refinement parameter y";

repCh::usage = "Expands Ch[p, q] Ch[p, q][1] and GV[p, q, n]";

Euler::usage = "TODO";

DSZ::usage = "TODO";

ToFB::usage = "TODO";

ToHC::usage = "TODO";

GiesekerSlope::usage = "GiesekerSlope[{r, dF, dB, ch2}, M]";
GiesekerDim::usage = "GiesekerDim[{r, dF, dB, ch2}]";

BogomolovDiscriminant::usage = "TODO";
SecondChernClass::usage = "SecondChernClass[gam] computes the second chern class."

ZLV::usage = "ZLV[{r, dF, dB, ch2}, T, m] is the large volume central charge for local F1 given by Z = -ch2 + dF m + (m^2 r)/2 + 2 dF T - m r T - 4 r T^2 + dB (-m + T)";

SpectralFlow::usage = "TODO";

Sigma::usage = "Sigma such that ZLV = gam.Sigma.Pi\[Transpose]";

MLV::usage = "LV monodromy";

MSF::usage = "MSF[mF, mB] -> Spectral flow monodromy by object Ch[mF, mB]";

MonodromyFromSphericalTwist::usage = "MonodromyFromSphericalTwist[mF, mB] -> Spherical twist monodromy by spherical object Ch[mF, mB]";

M1p::usage = "F1Global.pdf eq. (5.59): Conifold monodromy at special lambda";

M2p::usage = "F1Global.pdf eq. (5.59): Conifold monodromy at special lambda";

MonodromyOnCharge::usage = "MonodromyOnCharge[M, {r, dF, dB, ch2}] -> Action of M on charge gam = Sigma.Inverse[M].Inverse[Sigma]";

MonodromyOnTau::usage = "MonodromyOnCharge[M, tau] is the action on tau";


(* Large volume scattering *)

Rays::usage = "Rays[{r, dF, dB, ch2}, t, psi, m] computes s = Re(T) for the ray Re(e^{-i psi} Z) = 0 at a given value of t = Im(T)";

Wall::usage = "Wall[{r, dF, dB, ch2}, {rr, ddF, ddB, cch2}, {s, t}, m] computes the equation of the wall of marginal stability between the two objects. Only works for real m.";

WallRadius::usage = "WallRadius[{r, dF, dB, ch2}, {rr, ddF, ddB, cch2}, {s, t}, m] computes the radius of the wall circle. Only works for real m.";

WallCircle::usage = "WallCircle[{r, dF, dB, ch2}, {rr, ddF, ddB, cch2}, {s, t}, m] gives the Circle object for the wall circle. Only works for real m.";

DiscF1::usage = "DiscF1[gam, m] computes the ray discriminant of gam. Only works for real m.";

InitialPosition::usage = "InitialPosition[{r, dF, dB, ch2}, m] calculates the initial s value of the ray when it intersects the real T axis. Only works for real m.";

CollDBridgeland::usage = "";
CollDArcaraII::usage = "";
CollDArcaraI::usage = "";

CollBridgeland::usage = "";
CollArcaraII::usage = "";
CollArcaraI::usage = "";

(* simple captioning mechanism *)
GamToString::usage = "";
GamCaption::usage = "";

ListSpecFlows::usage = "ListSpecFlows[Coll, psi, m, L, xmin, xmax, tmax] lists spectral flow parameters (mF, mB) in the range -L <= m_i <= L that can be applied to the collection Coll such that the resulting quiver domain has a non-empty intersection with the region xmin < Re[T] < xmax, 0 < Im[T] < tmax."

(* XY plane scattering *)
XYRay::usage = "XYRay[{r, dF, dB, ch2}, {x, y}, mu] gives the ray equation in the XY coordinates.";
XYParabolaY::usage = "XYParabola[]";
InitialPositionXY::usage = "InitialPositionXY[gam, m] gives the initial position of ray gam in the x y plane (point with 0 cost function). Only works for real m.";
XYTost::usage = "XYTost[{x, y}, psi, m] gives the coordinates in s, t plane for a given point {x, y} in the XY coordinates. Only works for real m.";
XYRaySegment::usage = "XYRaySegment[Ray, mu, L, styleMap, defaultStyle] draws the line corresponding to Ray of length L with style given by a depth-style map styleMap.";
CollPlotData::usage = "Plot data for ConstructDomain";
RayStyleMap::usage = "Default depth -- style map for rays.";
ConstructDomain::usage = "ConstructDomain[PlotAssociation, mu, L, {x, y}] plots the quiver domain associated with the collection in PlotAssociation in the XY plane.";

NaiveIntersection::usage = "NaiveIntersection[gam, gamp, mu] gives the intersection points of the lines corresponding to gam and gamp without reference to the initial points and the future directions.";

CollInitial::usage = "The three sets of initial charges identified for fixed m slice.";
InitialRaysFromColl::usage = "InitialRaysFromColl[Coll, xmin, xmax, m] produces initial rays for ConstructLVDiagram by taking all translates of Coll that start between x = xmin and x = xmax. Only works for real m.";

CollInequalities::usage="CollInequalities[coll] gives the list of inequalities Re(e^(-i \\psi) Z_gam) < 0 for all gam in coll.";

(* from F0Scattering.m *)

GCD1::usage = "";
SameHalfPlaneQ::usage = "SameHalfPlaneQ[Zlist] gives True if all elements of Zlist are in a common half plane";

QuiverDomain::usage = "QuiverDomain[Coll, psi, m] plots the region where the LV central charges Z of Coll have Re[e^{-I psi} Z] < 0."; 

ChernToCh::usage = "ChernToCh[gam] replaces {r,dF,dB,ch2} by r Ch[dF/r,dB/r] whenever the discriminant vanishes";

ExtFromStrong::usage = "ExtFromStrong[Coll] computes the Chern vectors of the objects in the Ext collection dual to the given strong collection Coll"
StrongFromExt::usage = "StrongFromExt[Coll] computes the Chern vectors of the objects in the strong collection dual to the given Ext collection Coll";

MutateCollection::usage = "MutateCollection[Coll, klist] acts on the list of Chern vectors Coll by the successive mutations in klist, which is a list of {node number,sign}, with sign=1 for right mutation, -1 for left mutation";

ScattCheck::usage = "ScattCheck[Tree, m] returns {charge,{x,y}} of the root vertex if Tree is consistent, otherwise {total charge,{}}";
ScattSort::usage = "ScattSort[LiTree, m] sorts trees in LiTree by growing radius";
ScattGraph::usage = "ScattGraph[Tree, m] extracts the list of vertices and adjacency matrix of Tree";

FOmbToOm::usage = "FOmbToOm[OmbList] computes integer index from list of rational indices, used internally by FScattIndex";
ScattIndexImproved::usage = "ScattIndexImproved[TreeList, opt] computes the index for each tree in TreeList, taking care of non-primitive internal states";
(* ScattIndexImprovedInternal::usage = "ScattIndexImprovedInternal[Tree, opt] computes the index for Tree, taking care of non-primitive internal states"; *)

(* TreeFromListRays::usage = "TreeFormListRays[ListRays, k] "; *)
LVTreesFromListRays::usage = "LVTreesFromListRays[ListRays, gam, m] extract the trees with given charge in the List of rays, constructed by ConstructLVDiagram";

KroneckerDims::usage = "KroneckerDims[m, Nn] gives the list of populated dimension vectors {n1,n2} for Kronecker quiver with m arrows, with (n1,n2) coprime and 0<=n1,n2<=Nn"; 
IntersectRaysNoTest::usage = "IntersectRays[{r, dF, dB, ch2}, {rr, ddB, ddF, cch2}, z, zz, m] returns intersection point (x,y) of two rays if the intersection point lies strictly upward from z and z', or {} otherwise, without testing non-vanishing of DSZ product";  
CostPhi::usage = "CostPhi[{r, dF, dB, ch2}, s, mu] gives the cost function \\phi_s(\\gamma) = dB + 2 dF - r (mu + 8 s)";
ConstructLVDiagram::usage = "ConstructLVDiagram[smin, smax, phimax, Nm, m, ListRays] constructs the LV scattering diagram for F1 with initial rays in the interval [smin,smax], cost function up to phimax, scattering products with n1 + n2 <= Nn at each intersection; m is assumed to be real; The output consists of a list of  { charge, {x,y}, parent1, parent2, n1, n2 }; If ListRays is not empty, then uses it as initial rays.";
ConstructLVDiagramOpt::usage = "ConstructLVDiagram[phimax, Nm, m, InitialRays]";

(* Period expansions in Pi-Stability slice *)
F1Series::usage = "F1Series[la, u, Nn, Nr] computes large volume period expansion of F1 periods, its first and second derivatives using Richardson acceleration.";
PicardFuchsP::usage = "PicardFuchsP[la, u] computes the coefficient of d^2 t/d u^2 in the normalized 3rd order Picard Fuchs equations for the CY3 periods.";
PicardFuchsQ::usage = "PicardFuchsQ[la, u] computes the coefficient of d t/d u in the normalized 3rd order Picard Fuchs equations for the CY3 periods.";
PicardFuchs::usage = "PicardFuchs[f, la, u]";
PicardFuchsA::usage = "PicardFuchsA[la, u]";
SystemMatrix::usage = "SystemMatrix[{{t, td}, {dt, dtd}, {ddt, ddtd}}, la, u]";
ComputeTransition::usage = "ComputeTransition[la, upath, Nn, Nr]";

G10FundamentalDomain::usage = "G10FundamentalDomain plots the funamental domain of the noncongruence subgroup Gamma(10;0,3,0,1;8).";
G10RepeatDomain::usage = "G10RepeatDomain[{kmin, kmax}] repeats G10FundamentalDomain from x = 8*kmin to 8*kmax.";

MirrorCurveJ2::usage = "MirrorCurveJ2[el, ur]";
PicardFuchsP2::usage = "PicardFuchsP2[el, ur]";
PicardFuchsQ2::usage = "PicardFuchsQ2[el, ur]";
PicardFuchs2::usage = "PicardFuchs2[f, el, ur]";
PicardFuchsA2::usage = "PicardFuchsA2[el, ur]";
F1Series2::usage = "F1Series2[el, ur, Nn, Nr]";
SystemMatrix2::usage = "SystemMatrix2[dat, el, ur]";

Begin["`Private`"];

(* ===================== *)
(*   Internal helpers    *)
(* ===================== *)



(* ===================== *)
(*   Public definitions  *)
(* ===================== *)


RichardsonResum[PartSum_, n_] := 
  With[{d = Length[PartSum] - n - 1}, 
   Sum[PartSum[[Length[PartSum] + k - n]] (d + k)^
      n (-1)^(k + n)/k!/(n - k)!, {k, 0, n}]];




MirrorCurve[x_,y_,m_,u_] := x + y + 1/x/y + m/x - 1/u;

MirrorCurvef3[X_, m_, u_] := 5832 u^6 - 1728 m^3 u^6 - 
 432 m^2 u^4 (-3 + X) + (3 - 2 X)^2 (3 + X) + 324 u^3 (-3 + 2 X) + 
 108 m u^2 (-3 + 36 u^3 + 2 X);

MirrorCurveg2[m_, u_] := 27 (1 - 8 m u^2 - 24 u^3 + 16 m^2 u^4);

MirrorCurveg3[m_, u_] := 27 (-1 + 12 m u^2 + 36 u^3 - 48 m^2 u^4 - 144 m u^5 - 216 u^6 + 
   64 m^3 u^6);

MirrorCurveDelta[m_, u_] := m + u - 8 m^2 u^2 - 36 m u^3 - 27 u^4 + 16 m^3 u^4;

MirrorCurveJ[m_,u_]:=(1-8 m u^2-24 u^3+16 m^2 u^4)^3/(u^8 (m+u-8 m^2 u^2-36 m u^3-27 u^4+16 m^3 u^4));

URoot[la_, k_] := Root[la + #1 - 8 la^2 #1^2 - 
   36 la #1^3 + (-27 + 16 la^3) #1^4 &, k];

URoots[m_] := u /. NSolve[MirrorCurveDelta[m, u] == 0, u];

PlotURoots[m_] := ComplexListPlot[Thread[{URoots[m]}],
   PlotRange -> {{-1, 1}, {-1, 1}},
   PlotStyle -> {PointSize[.02]},
   AspectRatio -> 1,
   PlotLegends -> N[Abs /@ URoots[m]]
   ];

NodalCurveParam[X_,Y_,t_,X0_] := {X -> 1/2 (t^2 - 4 X0), Y -> (t (t^2 - 6 X0))/Sqrt[2]};

NodalCurveX0[m_,u_] := (3 + 12 u^2 (-3 m - 9 u + 12 m^2 u^2 + 36 m u^3 + 
    2 (27 - 8 m^3) u^4))/(2 + 16 u^2 (-3 u + m (-1 + 2 m u^2)));

NodalCurveX0tot1 = {X0->t1^2/2+t1};


WeierstrassChangeOfVariables[{x_, y_}->{X_,Y_}, m_, u_] := {x -> (-9 + 36 u^2 (m + 3 u) + 6 X + Sqrt[2] Y)/(
  6 u (-3 + 12 m u^2 + 2 X)), y -> (36 u^2)/(3 - 12 m u^2 - 2 X)};


WeierstrassChangeOfVariablesInverse[{X_,Y_}->{x_,y_},m_,u_]:={X -> 3/2 - (6 u^2 (3 + m y))/y, 
 Y -> -((54 Sqrt[2] u^2 (-1 + u (2 x + y)))/y)};


LiEval = {
  BW[x_] :> Im[PolyLog[2, x]] + Arg[1 - x] Log[Abs[x]], 
  Li2[z_] :> PolyLog[2, z]
};


RationalData[expr_, var_] := 
  Module[{t, num, den, fn, fd, all, 
    linear},(*make sure it's really a single rational expression*)
   t = Together[expr];
   num = Numerator[t];
   den = Denominator[t];
   (*factor numerator and denominator in var*)
   fn = Rest@FactorList[num];
   (*{{factor,exp},...}*)
   fd = Rest@FactorList[den];
   (*denominator factors get negative exponents*)
   (* all = Join[fn, {{#[[1]], -#[[2]]} & /@ fd} // Flatten[#, 1] &]; *)
   all = Join[fn, ({#[[1]], -#[[2]]} & /@ fd)];
   (*keep only true linear polynomials in var*)
   linear = 
    Select[all, 
     PolynomialQ[First[#], var] && Exponent[First[#], var] == 1 &];
   (*extract exponents and roots*)
   {
    (*(a) exponents*)
    linear[[All, 2]],
    (*(b) roots:for a z+b==0\[RightArrow]z=-b/a*)
    Simplify[-Coefficient[#[[1]], var, 0]/
        Coefficient[#[[1]], var, 1]] & /@ linear
    }];


NodalCurvebFromt1[t1_] := (-3 - 2 t1 + Sqrt[3] Sqrt[t1 (2 + t1)])/(3 + t1);

NodalCurvebFrommu[m_, u_] := (1 - 2 Sqrt[1 + 12 m u^2] + Sqrt[
 3 + 36 m u^2 - 6 Sqrt[1 + 12 m u^2]])/(1 + Sqrt[1 + 12 m u^2]);

KerrDoranP[b_] := (5 BW[b] - 4 BW[b^2] + BW[b^3])/(2 \[Pi]);

KerrDoranQ[b_] := 1/(6 \[Pi]) (3 (Log[1/b] + Log[b]) Log[1/(2 + 1/b + b)^(
     2/3)] + (Log[1/b^2] - Log[1/b] + Log[b]) Log[2 + 1/b + b] - 
   3 (Li2[1/b^3] - Li2[1/b^2] - Li2[b] + Li2[b^2] - 
      2 (Li2[1/b^2] + Log[1 - 1/b^2] Log[1/b^2]) + 
      Log[1 - b] Log[1/b] + (Log[1/b^2] - Log[1/b]) Log[(-1 + b)/b] + 
      2 (Li2[1/b] + Log[1/b] Log[(-1 + b)/b]) + 
      Log[1 - 1/b^3] (Log[1/b^2] - Log[b]) - 
      Log[1 - 1/b^2] (Log[1/b] - Log[b]) + Log[(-1 + b)/b] Log[b] - 
      2 (Li2[b] + Log[1 - b] Log[b]) - (Log[1/b] - Log[b]) Log[
        1 - b^2]));

(* add KerrDoran direct from m, u *)


repCh = {Ch[mF_, mB_][1] :> -{1, mF, mB, -(mB^2/2) + mB mF}, 
   Ch[mF_, mB_] :> {1, mF, mB, -(mB^2/2) + mB mF}, 
   GV[mF_, mB_, n_] :> {0, mF, mB, n}};

Euler[{r_, dF_, dB_, ch2_}, {rp_, dFp_, dBp_, ch2p_}] := dB dBp - dBp dF - dB dFp + ch2p r + (dBp r)/2 + dFp r + ch2 rp - (dB rp)/2 - dF rp + r rp;

DSZ[{r_, dF_, dB_, ch2_}, {rp_, dFp_, dBp_, ch2p_}] := dBp r + 2 dFp r - dB rp - 2 dF rp;

ToFB[{r_, dH_, dC_, ch2_}] := {r, dH, dH + dC, ch2};

ToHC[{r_, dF_, dB_, ch2_}] := {r, dF, dB - dF, ch2};

GiesekerSlope[{r_, dF_, dB_, ch2_}, M_]:= (dB - dB M + dF (2 + M))/r;

GiesekerDim[{r_, df_, db_, ch2_}] := 1 - db^2 + 2 db df - 2 ch2 r - r^2;

BogomolovDiscriminant[{r_, dF_, dB_, ch2_}] := -(dB^2/(2 r^2)) + (dB dF)/r^2 - ch2/r;
SecondChernClass[{r_, dF_, dB_, ch2_}] := -ch2 - dB^2/2 + dB dF;

ZLV[{r_, dF_, dB_, ch2_}, T_, m_] := -ch2 + dF m + (m^2 r)/2 + 2 dF T - m r T - 4 r T^2 + dB (-m + T);
ZLV[{r_, dF_, dB_, ch2_}, {s_,t_}, m_] := -ch2 + dB (-m + s + I t) + 
 1/2 (2 dF + r (m - 4 s - 4 I t)) (m + 2 s + 2 I t);


SpectralFlow[charges : {{_, _, _, _} ..}, sf : {_, _}] := SpectralFlow[#, sf] & /@ charges;

SpectralFlow[{r_, dF_, dB_, ch2_}, {mF_, mB_}] :=
  {r, dF + r mF, dB + r mB,
   ch2 - dB mB + dF mB + dB mF - (mB^2 r)/2 + mB mF r};

Sigma = {{0, 0, 0, -1}, {0, 1, 2, 0}, {0, -1, 1, 0}, {-1, 0, 0, 0}};

MLV = {{1, 0, 0, 0}, {0, 1, 0, 0}, {1, 0, 1, 0}, {4, 1, 8, 1}};

MSF[mF_, mB_] := {{1, mF, mB, -(mB^2/2) + mB mF}, {0, 1, 0, mB}, {0, 0, 
    1, -mB + mF}, {0, 0, 0, 1}};

(* MonodromyFromSphericalTwist[mF_, mB_] := {{1, 0, 0, 0}, {0, 1, 0, 0}, {1/2 mB (mB - 2 mF), -mB + mF,
     1 + mB + 2 mF, -1}, {1/2 mB (mB^2 - 4 mF^2), -mB^2 - mB mF + 
     2 mF^2, (mB + 2 mF)^2, 1 - mB - 2 mF}}; *)
MonodromyFromSphericalTwist[{rr_,ddF_,ddB_,cch2_}] := {{1, 0, 0, 0}, {0, 1, 0, 0}, {-cch2 rr, (-ddB + ddF) rr, 
  1 + ddB rr + 
   2 ddF rr, -rr^2}, {-cch2 (ddB + 2 ddF), (-ddB + ddF) (ddB + 
     2 ddF), (ddB + 2 ddF)^2, 1 - (ddB + 2 ddF) rr}};

M1p = {{1, 0, 0, 0}, {0, 1, 0, 0}, {-1/2, 1, 0, 1}, {-1/2, 1, -1, 1}};
M2p = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 1, 3, 1}, {0, -2, -4, -1}};

MonodromyOnCharge[M_, gam_] :=
   If[Length[M] == 0, gam,
    If[ArrayDepth[M] == 2,
     (*a single matrix*)
     gam . Sigma . Inverse[M] . Inverse[Sigma],
     (*a list of matrices*)
     MonodromyOnCharge[Drop[M, 1], 
      gam . Sigma . Inverse[First[M]] . Inverse[Sigma]]
     ]
    ];

MonodromyOnTau[M_, tau_] := (M[[4, 3]] + tau M[[4, 4]])/(M[[3, 3]] + tau M[[3, 4]]);


Rays[{r_, dF_, dB_, ch2_}, t_, psi_, m_ : 0] := (
   ch2 + (dB - dF) Re[m] - (dB + 2 dF) t Tan[psi] + (dB - dF) Im[
      m] Tan[psi])/(dB + 2 dF) /; r == 0;
Rays[{r_, dF_, dB_, ch2_}, t_, psi_, m_ : 0] := -(1/(
    8 r)) (-dB - 2 dF + r Re[m] + r (8 t + Im[m]) Tan[psi] + 
     Sec[psi] \[Sqrt](Cos[
          psi]^2 ((dB + 2 dF)^2 + 
           9 r^2 (-Im[m] + Re[m]) (Im[m] + Re[m]) - 
           2 r (8 ch2 + (9 dB - 6 dF) Re[m]) + 
           r^2 (8 t + Im[m])^2 Sec[psi]^2 + 
           6 r Im[m] (-3 dB + 2 dF + 3 r Re[m]) Tan[psi])));

Wall[{r_, dF_, dB_, ch2_}, {rr_, ddF_, ddB_, cch2_}, {s_, t_}, 
   m_ : 0] := (cch2 (dB + 2 dF - m r) - ch2 (ddB + 2 ddF - m rr) + 
     1/2 m (-6 dB ddF + 6 ddB dF - ddB m r + 
        4 ddF m r + (dB - 4 dF) m rr) + 
     8 (-((cch2 + (ddB - ddF) m) r) + (ch2 + (dB - dF) m) rr) s + 
     4 (ddB r + 2 ddF r - (dB + 2 dF) rr) (s^2 + t^2))/(4 (ddB r + 
       2 ddF r - dB rr - 2 dF rr));

WallRadius[{r_, dF_, dB_, ch2_}, {rr_, ddF_, ddB_, cch2_}, 
   m_ : 0] := (r (2 (ddB + 2 ddF) (ch2 ddB + 2 ch2 ddF - 
           cch2 (dB + 2 dF) + 3 dB ddF m - 3 ddB dF m) + (2 cch2 + 
           3 ddB m) (4 cch2 + 3 (ddB - 2 ddF) m) r) + 
     2 ((dB + 2 dF) (-ch2 (ddB + 2 ddF) + cch2 (dB + 2 dF) - 
           3 dB ddF m + 3 ddB dF m) + (-8 cch2 ch2 + 
           3 (-3 cch2 dB - 3 ch2 ddB + 2 ch2 ddF + 2 cch2 dF) m + 
           9 (dB (-ddB + ddF) + ddB dF) m^2) r) rr + (2 ch2 + 
        3 dB m) (4 ch2 + 3 (dB - 2 dF) m) rr^2)/(8 (ddB r + 
       2 ddF r - (dB + 2 dF) rr)^2);

WallCircle[{r_, dF_, dB_, ch2_}, {rr_, ddF_, ddB_, cch2_}, m_ : 0] := 
  Module[{R},
   R = WallRadius[{r, dF, dB, ch2}, {rr, ddF, ddB, cch2}, m];
   If[R > 0, 
    Circle[{(cch2 r + ddB m r - ddF m r - ch2 rr - dB m rr + 
       dF m rr)/(ddB r + 2 ddF r - dB rr - 2 dF rr), 0}, 
     Sqrt[R], {0, Pi}], {}]];

DiscF1[{r_, dF_, dB_, ch2_}, m_] := 
  16 BogomolovDiscriminant[{r, dF, dB, ch2}] + (-3 dB + 2 dF + 
     3 m r)^2/r^2;

InitialPosition[{r_, dF_, dB_, ch2_}, m_ : 0] :=
  Module[{}, 
   If[! IntegerQ[SecondChernClass[{r,dF,dB,ch2}]], 
    Print["Non integer second Chern class !"]];
   If[r == 0,
    (ch2 + dB m - dF m)/(dB + 2 dF),
    (dB + 2 dF - m r)/8/r - 
     Sign[r]/8 Sqrt[Max[DiscF1[{r, dF, dB, ch2}, m], 0]]
    ]
   ];

CollDBridgeland = {{1, 0, 0, 0}, {1, 1, 0, 0}, {1, 1, 1, 1/2}, {1, 2, 1, 3/2}};

CollDArcaraII = {{1, 0, 0, 0}, {1, 1, 0, 0}, {1, 1, 1, 1/2}, {2, 1, 1, -(1/2)}};

CollDArcaraI = {{1, 0, 0, 0}, {1, 0, 1, -(1/2)}, {1, 1, 1, 1/2}, {2, 1, 1, -(1/2)}};

CollBridgeland = {{1, 0, 0, 0}, {-1, 1, 0, 0}, {-1, 0, 1, 1/2}, {1, -1, -1, 1/2}};

CollArcaraII = {{1, 0, 0, 0}, {0, 0, -1, 1/2}, {1, -2, -1, 3/2}, {-1, 1, 1, -(1/2)}};

CollArcaraI = {{1, 0, 0, 0}, {0, 0, 1, -(1/2)}, {1, -2, -1, 3/2}, {-1, 1, 0, 0}};

(* simple captioning mechanism *)
GamToString = {
   -Ch[p_, q_] :> "(" <> ToString[p] <> ", " <> ToString[q] <> ")[1]",
   Ch[p_, q_] :> "(" <> ToString[p] <> ", " <> ToString[q] <> ")"
   };

GamCaption[gam_, t_, psi_, m_] := Module[{s, tt, pos, dir, normal},
  s = Sign[gam[[1]]];
  pos[tt_] = {Rays[gam /. repCh, tt, psi, m], tt};
  dir = Normalize[D[pos[tt], tt] /. {tt -> t}];
  normal = {dir[[2]], -dir[[1]]};
  Text[gam /. GamToString, pos[t], -s*1.5*normal, -s*dir]
  ];

ListSpecFlows[Coll_, psi_, m_, L_, xmin_, xmax_, tmax_] := Module[{s, t}, 
Select[Flatten[Table[{mF, mB}, {mF, -L, L}, {mB, -L, L}], 1], 
 With[{z = ZLV[#, {s, t}, m] & /@ SpectralFlow[Coll, #]}, 
    Reduce[FullSimplify@
      ComplexExpand[
       AllTrue[z, 
         Re[Exp[-I psi] #] < 0 &] && (xmin < s < xmax && 
          0 < t < tmax)], {s, t}]] =!= False &]
];


(* XY plane scattering *)

XYRay[{r_, dF_, dB_, ch2_}, {x_, y_}, mu_] := 
  r y + (2 dF + dB) x + (dF - dB) mu - ch2;

XYParabolaY[x_,mu_,m2_,psi_] := 1/4 (18 m2^2 + mu^2 - 2 mu x - 
    8 x^2 + (mu^2 - 2 mu x - 8 x^2) Cos[2 psi]) Sec[psi]^2;


XYTost[{x_, y_}, psi_, m_] := {x - (Sqrt[-(m^2/2) + m x + 4 x^2 + y] Tan[psi])/(
    2 Sqrt[Sec[psi]^2]), Sqrt[-(m^2/2) + m x + 4 x^2 + y]/(
   2 Sqrt[Sec[psi]^2])};

RayStyleMap = <|
   0 -> Directive[Black, Thickness[.002]],
   1 -> Directive[Opacity[.6, Red], Thickness[.0015]],
   2 -> Directive[Opacity[.3, Blue], Thickness[.001]],
   3 -> Directive[Opacity[.1, Green], Thickness[.0007]]
   |>;

XYRaySegment[{c : {r_, dF_, dB_, ch2_}, z : {x0_, y0_}, ___, depth_}, mu_, 
   L_ : 10, styleMap_ : <||>, 
   defaultStyle_ : Directive[Black, Thickness[.0002]]] := 
  Module[{v, sty},
   v = {-r, dB + 2 dF};
   sty = Lookup[styleMap, depth, defaultStyle];
   {sty, Line[{z, z + L Normalize[v]}]}];

CollInequalities[Coll_, {x_, y_}, mu_] := 
  Table[gam[[
       1]] y + (gam[[3]] + 2 gam[[2]]) x + (gam[[2]] - gam[[3]]) mu - 
     gam[[4]] < 0, {gam, Coll}];

CollPlotData = {
   <|"coll" -> CollBridgeland, 
    "mB" -> Function[{m, mu}, (2 m + 3 mu - 1)/3], 
    "color" -> RGBColor[0, 0.8260000000000001, 1]|>,
   <|"coll" -> CollArcaraI, 
    "mB" -> Function[{m, mu}, (2 m + 3 mu - 3)/3], 
    "color" -> RGBColor[0.46900000000000003`, 1, 0]|>,
   <|"coll" -> CollArcaraII, 
    "mB" -> Function[{m, mu}, (2 m + 3 mu - 1)/3], 
    "color" -> RGBColor[1, 0, 0.192]|>
   };

ConstructDomain[asc_, mu_, L_] := 
  Module[{ineq, x, y}, 
   ineq = Table[
     CollInequalities[
      SpectralFlow[
       asc["coll"] /. repCh, {mF, Ceiling[asc["mB"][mF, mu]]}], {x, y},
       mu], {mF, -L, L}];
   {EdgeForm[Opacity[.25, asc["color"]]], 
    DiscretizeRegion[ImplicitRegion[And @@ #, {x, y}]] & /@ ineq}
   ];


(* from F0Scattering.m *)
SameHalfPlaneQ[{}] := True;
SameHalfPlaneQ[Zlist_List] := 
  If[AnyTrue[Zlist, # == 0 &], 
   False, -Subtract @@ MinMax[Arg[Zlist/Zlist[[1]]]] < Pi];

Options[QuiverDomain] = {PlotStyle -> LightBlue, 
   PlotRange -> {{-1.5, 1.5}, {0, 1}}, PlotPoints -> 100};
QuiverDomain[coll_, psi_, m_, OptionsPattern[]] := 
  Module[{sty = OptionValue[PlotStyle], pr = OptionValue[PlotRange], 
    pp = OptionValue[PlotPoints], z = ZLV[#, {s, t}, m] & /@ coll}, 
   RegionPlot[
    t > 0 && AllTrue[z, Re[Exp[-I psi] #] < 0 &], {s, pr[[1, 1]], 
     pr[[1, 2]]}, {t, pr[[2, 1]], pr[[2, 2]]},
    PlotPoints -> pp,
    AspectRatio -> 1,
    PlotStyle -> {sty, Opacity[.5]},
    BoundaryStyle -> None
    ]];

ExtFromStrong[Coll_]:=Module[{S,Si},
   S=Table[Euler[Coll[[i]],Coll[[j]]],{i,Length[Coll]},{j,Length[Coll]}];
   Si=Inverse[Transpose[S]];
   Si . Coll
   ];

StrongFromExt[Coll_]:=Module[{S,Si},
   S=Table[Euler[Coll[[j]],Coll[[i]]],{i,Length[Coll]},{j,Length[Coll]}];
   Si=Inverse[Transpose[S]];
   Si . Coll
   ];

(* Coll is a list of Chern vectors, klist a list of {node,\pm 1} *)
MutateCollection[Coll_, klist_] := Module[{Coll0, k, eps},
    Coll0 = 
   If[Length[klist] > 1, MutateCollection[Coll, Drop[klist, -1]], 
    Coll];
    k = Last[klist][[1]]; eps = Last[klist][[2]];
    Table[
   If[i == k, -Coll0[[k]], 
    Coll0[[i]] + 
     Max[0, eps DSZ[Coll0[[i]], Coll0[[k]]]] Coll0[[k]]], {i, 
    Length[Coll0]}]];

TreeFromListRays[ListRays_,k_]:=If[ListRays[[k,3]]==0,ListRays[[k,1]],{ListRays[[k,5]]TreeFromListRays[ListRays,ListRays[[k,3]]],ListRays[[k,6]]TreeFromListRays[ListRays,ListRays[[k,4]]]}];

GCD1[{r_,dF_,dB_,ch2_}]:=Module[{d},d=GCD[r,dF,dB];
If[EvenQ[(dB-2ch2)/d],d,If[EvenQ[d],d/2,1]]];

(*Check consistency of single tree,returns {charge,{xf,\
yf}} if tree is consistent,otherwise {charge,{}};\
ignore whether leaves have Delta=0 or not*)
ScattCheck[Tree_, m_] :=
  Module[{S1, S2, z, r, dH, dC, ch2},
   If[! ListQ[Tree] || Length[Tree] > 2,
    (*tree consists of a single node*)
    {r, dF, dB, ch2} = Tree /. repCh;
    z = {(dB + 2 dF - 
         m r)/(8 r), -((dB^2 + 4 dB dF + 4 dF^2 - 8 ch2 r - 
           9 dB m r + 6 dF m r)/(8 r^2))};
    (*Print["Initial pt:",z];*)
    {Tree /. repCh, z},
    (*otherwise,check each of the two branches*)
    S1 = ScattCheck[Tree[[1]], m];
    S2 = ScattCheck[Tree[[2]], m];
    If[Length[S1[[2]]] > 0 && Length[S2[[2]]] > 0,
     z = 
      IntersectRays[S1[[1]] /. repCh, S2[[1]] /. repCh, S1[[2]], 
       S2[[2]], m];
     (*Print[{S1[[1]],S2[[1]],S1[[2]],S2[[2]],z}];*)
     {S1[[1]] + S2[[1]], z}, {S1[[1]] + S2[[1]], {}}
     ]
    ]
   ];

ScattSort[LiTree_, m_ : 0] :=(*sort trees by decreasing radius*)
  Reverse[Map[#[[2]] &, 
    SortBy[Table[{XYTost[ScattCheck[LiTree[[i]], m][[2]], 0, m][[2]], 
       LiTree[[i]]}, {i, Length[LiTree]}], N[First[#]] &]]];


(*construct total charge,\
coordinate of root and list of line segments in (s,t) coordinates,\
{min(x),max(x)}*)
ScattGraphInternal[Tree_, m_ : 0] := 
  Module[{S1, S2, TreeNum, sInit, z, Li},
   If[! ListQ[Tree] || Length[Tree] > 2,
    (* branch 1 *)
    TreeNum = Tree /. repCh;
    sInit = InitialPosition[TreeNum, m];
    {Tree, {sInit, m^2/2 - m sInit - 4 sInit^2}, {}},
    (* branch 2 *)
    S1 = ScattGraphInternal[Tree[[1]], m];
    S2 = ScattGraphInternal[Tree[[2]], m];
    z = IntersectRays[S1[[1]] /. repCh, S2[[1]] /. repCh, m];
    If[Length[z] == 0,
     Print["Illegal tree"],
     Li = {S1[[3]], S2[[3]], Arrow[{S1[[2]], z}], Arrow[{S2[[2]], z}]};
     {S1[[1]] + S2[[1]], z, Li}
     ]
    ]
   ];

(*extracts list of vertices in (x,y) plane and adjacency matrix*)
ScattGraph[Tree_, m_ : 0] := Module[{T, LiArrows, LiVertex},
   T = ScattGraphInternal[Tree, m];
   LiArrows = Cases[Flatten[T[[3]]], x_Arrow] /. Arrow[x_] :> x;
   LiVertex = Union[Flatten[LiArrows, 1]];
   {LiVertex, 
    Table[If[i != j, Sign[Count[LiArrows, {LiVertex[[i]], LiVertex[[j]]}]],
       0], {i, Length[LiVertex]}, {j, Length[LiVertex]}]}
   ];

LVTreesFromListRays[ListRays_,{r_,dF_,dB_,ch2_},m_]:=Module[{Lipos,div,LiTrees},
   div=Divisors[GCD1[{r,dF,dB,ch2}]];
   Lipos=Flatten[Join[Table[Position[ListRays,{r,dF,dB,ch2}/k],{k,div}]],1];
   If[Lipos=={},
   Print["No such dimension vector in the list"],
   LiTrees=(GCD1[{r,d1,d2,ch2}]/GCD1[ListRays[[#,1]]])TreeFromListRays[ListRays,#]&/@First[Transpose[Lipos]];
   ScattSort[DeleteDuplicatesBy[SortBy[LiTrees,Length[TreeConstituents[#]]&],ScattGraph[#,m]&],m]
]];


FOmbToOm[OmbList_, y_] := Module[{n},
   If[Length[OmbList] < 2,
    First@OmbList,
    n = Length[OmbList];
    DivisorSum[
     n, (MoebiusMu[#] (y - 
           y^-1)/(# (y^# - y^-#)) (OmbList[[n/#]] /. {y -> y^#})) &]]
   ];

(*compute index for each tree in the list*)
ScattIndexImproved[TreeList_, y_, opt : OptionsPattern[]] := 
  Table[Simplify[
    FOmbToOm[
     Last@ScattIndexImprovedInternal[TreeList[[i]], y, opt][[2]], 
     y]], {i, Length[TreeList]}];

Options[ScattIndexImprovedInternal] = {"Debug" -> False};
ScattIndexImprovedInternal[Tree_, y_, opt : OptionsPattern[]] := 
  Module[{S1, S2, g1, g2, gFinal, kappa, Li, tem, repOmAttb, 
    rrr},(*compute {total charge,
   list of Kronecker indices associated to each vertex*)
   If[! ListQ[Tree] || 
     Length[Tree] > 
      2, {Tree, {Join[{1}, 
       Table[(y - y^-1)/(j (y^j - y^-j)), {j, 2, GCD1@(Tree /. repCh)}]]}},
     If[OptionValue["Debug"], 
     Print["Calling with args: ", Tree[[1]], "  |  ", Tree[[2]]]];
    S1 = ScattIndexImprovedInternal[Tree[[1]], y, opt] /. repCh;
    S2 = ScattIndexImprovedInternal[Tree[[2]], y, opt] /. repCh;
    If[OptionValue["Debug"], Print["S1 is: ", S1, "   S2 is: ", S2]];
    g1 = GCD1@S1[[1]]; g2 = GCD1@S2[[1]];
    gFinal = GCD1@(S1[[1]] + S2[[1]]);
    kappa = Abs[DSZ[S1[[1]], S2[[1]]]]/g1/g2;
    Li = Join[S1[[2]], S2[[2]]];
    If[OptionValue["Debug"], 
     Print["Li is: ", Li, "  g1 is: ", g1, "  g2 is: ", g2, 
      "  gFinal is: ", gFinal]];
    AppendTo[Li, 
     repOmAttb = 
      Join[Table[
        CoulombHiggs`OmAttb[{P, 0}, _] -> Last[S1[[2]]][[P]], {P, 1, 
         g1}], Table[
        CoulombHiggs`OmAttb[{0, Q}, _] -> Last[S2[[2]]][[Q]], {Q, 1, 
         g2}]];
     If[OptionValue["Debug"], Print["repOmAttb is: ", repOmAttb]];
     tem = 
      Table[rrr = 
        If[And @@ (IntegerQ /@ {P g1/gFinal, P g2/gFinal}), 
         CoulombHiggs`FlowTreeFormulaRat[{{0, kappa}, {-kappa, 
            0}}, {g2, -g1}, {P g1/gFinal, P g2/gFinal}, y], 0];
       Simplify[
        rrr /. repOmAttb /. {CoulombHiggs`OmAttb[{p_, q_}, ___] :> 
           0 /; p > 1 || q > 1 || p q != 0}], {P, 1, gFinal}];
     If[OptionValue["Debug"], Print["tem is: ", tem]]; tem];
    (*If[GCD@@(S1[[1]]+S2[[1]])!=1,Print[
    "Beware, non-primitive state"]];*){S1[[1]] + S2[[1]], Li}]];

ChernToCh[f_]:=f/.{{r_Integer,dF_Integer,dB_Integer,ch2_}:>If[r==0,GV[dF,dB,ch2],If[BogomolovDiscriminant[{r,dF,dB,ch2}]==0,If[r>0,r Ch[dF/r,dB/r],-r Ch[dF/r,dB/r][1]],{r,dF,dB,ch2}]]};	

KroneckerDims[m_Integer?NonNegative, Nn_Integer?NonNegative] := KroneckerDims[m,Nn] = Module[{Ta={}},
   Do[If[m n1 n2-n1^2-n2^2+1>=0&&GCD[n1,n2]==1,AppendTo[Ta,{n1,n2}]],{n1,0,Nn},{n2,0,Nn-n1}];Drop[Ta,2]];

(* IntersectRaysNoTest[{r_, dF_, dB_, ch2_}, {rr_, ddF_, ddB_, cch2_}, 
   z_, zz_, mu_ : 0] :=(*returns (x,
  y) coordinate of intersection point of two rays,
  or {} if they don't intersect*)(*here do not test if DSZ<>0,
  and require strictly in future of z and zz*)
  Module[{zi}, 
   zi = {(cch2 r + ddB mu r - ddF mu r - ch2 rr - dB mu rr + 
      dF mu rr)/(ddB r + 2 ddF r - dB rr - 2 dF rr), (
     ch2 ddB + 2 ch2 ddF - cch2 (dB + 2 dF) + 3 dB ddF mu - 
      3 ddB dF mu)/(ddB r + 2 ddF r - (dB + 2 dF) rr)};
   If[CostPhi[{r, dF, dB, ch2}, zi[[1]], mu] > 
      CostPhi[{r, dF, dB, ch2}, z[[1]], mu] && 
     CostPhi[{rr, ddF, ddB, cch2}, zi[[1]], mu] > 
      CostPhi[{rr, ddF, ddB, cch2}, zz[[1]], mu], zi, {}]]; *)

IntersectRaysNoTest[{r_, dF_, dB_, ch2_}, {rr_, ddF_, ddB_, cch2_}, 
   z_, zz_, mu_ : 0] :=
  (*returns (x,y) coordinate of intersection point of two rays,
  or {} if they don't intersect*)(*here do not test if DSZ<>0,
  and require strictly in future of z and zz*)Module[{zi, eta, eeta},
   zi = NaiveIntersection[{r, dF, dB, ch2}, {rr, ddF, ddB, cch2}, mu];
   eta = {-r, 2 dF + dB};
   eeta = {-rr, 2 ddF + ddB};
   If[eta . (zi - z) > 0 && eeta . (zi - zz) > 0, zi, {}]];

IntersectRays[{r_, dF_, dB_, ch2_}, {rr_, ddF_, ddB_, cch2_}, 
   mu_ : 0] := 
  If[DSZ[{r, dF, dB, ch2}, {rr, ddF, ddB, cch2}] == 0, {}, 
   NaiveIntersection[{r, dF, dB, ch2}, {rr, ddF, ddB, cch2}, mu]];

IntersectRays[{r_, dF_, dB_, ch2_}, {rr_, ddF_, ddB_, cch2_}, z_, zz_,
    mu_ : 0] :=
  Module[{zi, eta, eeta},
   If[DSZ[{r, dF, dB, ch2}, {rr, ddF, ddB, cch2}] == 0, Return[{}]];
   zi = NaiveIntersection[{r, dF, dB, ch2}, {rr, ddF, ddB, cch2}, mu];
   eta = {-r, 2 dF + dB};
   eeta = {-rr, 2 ddF + ddB};
   If[eta . (zi - z) >= 0 && eeta . (zi - zz) >= 0, zi, {}]];
   

CostPhi[{r_,dF_,dB_,ch2_},s_,m_]:=dB + 2 dF - r (m + 8 s);

InitialPositionXY[{r_, dF_, dB_, ch2_}, m_] := {(dB + 2 dF - m r)/(
   8 r), -((dB^2 + 4 dB dF + 4 dF^2 - 8 ch2 r - 9 dB m r + 6 dF m r)/(
    8 r^2))};

NaiveIntersection[{r_, dF_, dB_, ch2_}, {rr_, ddF_, ddB_, cch2_}, mu_] := {(cch2 r + ddB mu r - ddF mu r - ch2 rr - dB mu rr + dF mu rr)/(
 ddB r + 2 ddF r - dB rr - 2 dF rr), (
 ch2 ddB + 2 ch2 ddF - cch2 (dB + 2 dF) + 3 dB ddF mu - 3 ddB dF mu)/(
 ddB r + 2 ddF r - (dB + 2 dF) rr)};

(* InitialRaysFromColl[Coll_, xmin_, xmax_, m_] := Module[{x, y},
   Flatten[
    Table[
     {x, y} = InitialPositionXY[gam, m];
     Table[
      With[{newgam = SpectralFlow[gam, {3 k, 2 k}]},
       {sign*newgam, InitialPositionXY[newgam, m], 0, 0, 0, 0, 0}
       ], {k, 
       Range[Ceiling[xmin - x], Floor[xmax - x]]}, {sign, {-1, 
        1}}], {gam, Coll}], 2]
   ]; *)

CollInitial = {
   {Ch[0, 0], Ch[1, 1], Ch[2, 1], Ch[2, 2]},
   {Ch[0, 0], Ch[0, 1], Ch[1, 1], Ch[2, 2]},
   {Ch[0, 0], Ch[1, 0], Ch[1, 1], Ch[2, 1]}
};

InitialRaysFromColl[Coll_, xmin_, xmax_, m_] := 
  Module[{x, y, first, last, rest, inter, init},
   first = 
    NaiveIntersection[First[Coll], SpectralFlow[Last[Coll], {-3, -2}],
      m];
   last = 
    NaiveIntersection[Last[Coll], SpectralFlow[First[Coll], {3, 2}], 
     m];
   rest = 
    Table[NaiveIntersection[Coll[[i]], Coll[[i + 1]], m], {i, 
      Length[Coll] - 1}];
   inter = Join[{first}, rest, {last}];
   init = 
    Simplify[
     Table[(inter[[i]] + inter[[i + 1]])/2, {i, Length[inter] - 1}]];
   Flatten[
    Table[
     {x, y} = init[[i]];
     Table[With[{gamp = SpectralFlow[Coll[[i]], {3 k, 2 k}]},
       {sign*gamp,
        {k + x, -4 k^2 - k (m + 8 x) + y}, 0, 0, 0, 0, 0}], {k, 
       Range[Ceiling[xmin - x], Floor[xmax - x]]}, {sign, {-1, 1}}
      ], {i, Length[Coll]}], 2]];


ConstructLVDiagram[phimax_, Nm_, m_, ListRays0_] := 
  Module[{Inter, ListInter, ListRays, ListNewRays, kappa, KTab},
   (*initial rays {charge,{x,y},parent1,parent2,n1,n2}*)
   
   (* The list of rays IS already provided *)
   ListRays = ListRays0;
   ListInter = Select[ListRays[[All, {3, 4}]], First[#] > 0 &];
   
   CheckAbort[
    While[True, ListNewRays = {};
     Monitor[
      Do[If[! MemberQ[ListInter, {i, j}], AppendTo[ListInter, {i, j}];
        kappa = DSZ[ListRays[[i, 1]], ListRays[[j, 1]]];
        If[kappa != 0, 
         Inter = IntersectRaysNoTest[ListRays[[i, 1]], ListRays[[j, 1]],
            ListRays[[i, 2]], ListRays[[j, 2]], m];
         If[Inter != {}, KTab = KroneckerDims[Abs[kappa], Nm];
          
          Do[If[CostPhi[
              KTab[[k, 1]] ListRays[[i, 1]] + 
               KTab[[k, 2]] ListRays[[j, 1]], Inter[[1]], m] <= 
             phimax, AppendTo[
             ListNewRays, {KTab[[k, 1]] ListRays[[i, 1]] + 
               KTab[[k, 2]] ListRays[[j, 1]], Inter, i, j, KTab[[k, 1]],
               KTab[[k, 2]]}]], {k, Length[KTab]}]]]], {i, 
        Length[ListRays]}, {j, i + 1, Length[ListRays]}], {i, j}];
     If[ListNewRays == {}, Break[], 
      Print["Adding ", Length[ListNewRays], " rays, "];
      ListRays = Flatten[{ListRays, ListNewRays}, 1];]];
    Print[Length[ListRays], " in total."];
    ListRays,
    Print["Aborted \[LongDash] returning partial ListRays (", 
     Length[ListRays], ")."];
    ListRays
    ]];


ConstructLVDiagramOpt[phimax_, Nm_, m_, ListRays0_] :=
 Module[{ListRays = ListRays0, seen, n, new, kappa, inter, ktab},
   seen = AssociationThread[
     Select[ListRays[[All, {3, 4}]], First[#] > 0 &] -> True
  ];
  
  CheckAbort[
   While[True,
    n = Length[ListRays];
    
    new = Reap[
       Do[
        If[! KeyExistsQ[seen, {i, j}],
         seen[{i, j}] = True;
         
         kappa = DSZ[ListRays[[i, 1]], ListRays[[j, 1]]];
         If[kappa =!= 0,
          inter = IntersectRaysNoTest[
            ListRays[[i, 1]], ListRays[[j, 1]],
            ListRays[[i, 2]], ListRays[[j, 2]], m
          ];
          
          If[inter =!= {},
           ktab = KroneckerDims[Abs[kappa], Nm];
           Do[
            If[
             CostPhi[ktab[[k, 1]] ListRays[[i, 1]] + ktab[[k, 2]] ListRays[[j, 1]],
               inter[[1]], m] <= phimax,
             Sow[{
               ktab[[k, 1]] ListRays[[i, 1]] + ktab[[k, 2]] ListRays[[j, 1]],
               inter, i, j, ktab[[k, 1]], ktab[[k, 2]],
               1 + Max[ListRays[[i, 7]], ListRays[[j, 7]]]
               }]
            ],
            {k, Length[ktab]}
           ];
          ]];
         ],
        {i, n}, {j, i + 1, n}
       ]
      ][[2]];
    
    new = If[new === {}, {}, First[new]];
    If[new === {} || new === Null, Break[],
     Print["Adding ", Length[new], " rays, "];
     ListRays = Join[ListRays, new];
    ];
   ];
   Print[Length[ListRays], " in total."];
   ListRays,
   
   Print["Aborted -- returning partial ListRays (", Length[ListRays], ")."];
   ListRays
  ]
 ];



(* Period expansions in Pi-Stability slice *)

F1Series[la_, u_, Nn_, Nr_] := 
  Module[{varpi2tab, varpi3tab, varpi2dtab, varpi3dtab, varpi2ddtab, 
    varpi3ddtab, k, l},
   varpi2tab = Table[u^n Sum[k = n - 2 l; 
       (k + 2 l)!/((l)! (k!)^2 (l - k)!) la^(l - k) /(k + 2 l),
       {l, 0, Floor[n/2]}], {n, 1, Nn}];
   varpi3tab = Table[u^n Sum[k = n - 2 l; 
       la^(l - k) (If[
           k > l, ((-1)^(k + l) Gamma[k - l] Gamma[1 + k + 2 l])/(
           4 (k + 2 l) l! Gamma[1 + k]^2),
           (( (k + 2 l)! (4 HarmonicNumber[k] + 3 HarmonicNumber[l] + 
                 HarmonicNumber[-k + l] - 
                 8 HarmonicNumber[k + 2 l]))/(4 (k + 
                 2 l) (k!)^2 l! Gamma[1 - k + l]))] + (
          2  (k + 2 l)!)/((k + 2 l)^2 (k!)^2 l! (-k + l)!)),
       {l, 0, Floor[n/2]}], {n, 1, Nn}];
   varpi2dtab = Table[n u^(n - 1) Sum[k = n - 2 l; 
       (k + 2 l)!/((l)! (k!)^2 (l - k)!) la^(l - k) /(k + 2 l),
       {l, 0, Floor[n/2]}], {n, 1, Nn}];
   varpi3dtab = Table[n u^(n - 1) Sum[k = n - 2 l; 
       la^(l - k) (If[
           k > l, ((-1)^(k + l) Gamma[k - l] Gamma[1 + k + 2 l])/(
           4 (k + 2 l) l! Gamma[1 + k]^2),
           (( (k + 2 l)! (4 HarmonicNumber[k] + 3 HarmonicNumber[l] + 
                 HarmonicNumber[-k + l] - 
                 8 HarmonicNumber[k + 2 l]))/(4 (k + 
                 2 l) (k!)^2 l! Gamma[1 - k + l]))] + (
          2  (k + 2 l)!)/((k + 2 l)^2 (k!)^2 l! (-k + l)!)),
       {l, 0, Floor[n/2]}], {n, 1, Nn}];
   varpi2ddtab = Table[n (n - 1) u^(n - 2) Sum[k = n - 2 l; 
       (k + 2 l)!/((l)! (k!)^2 (l - k)!) la^(l - k) /(k + 2 l),
       {l, 0, Floor[n/2]}], {n, 2, Nn}];
   varpi3ddtab = Table[n (n - 1) u^(n - 2) Sum[k = n - 2 l; 
       la^(l - k) (If[
           k > l, ((-1)^(k + l) Gamma[k - l] Gamma[1 + k + 2 l])/(
           4 (k + 2 l) l! Gamma[1 + k]^2),
           (( (k + 2 l)! (4 HarmonicNumber[k] + 3 HarmonicNumber[l] + 
                 HarmonicNumber[-k + l] - 
                 8 HarmonicNumber[k + 2 l]))/(4 (k + 
                 2 l) (k!)^2 l! Gamma[1 - k + l]))] + (
          2  (k + 2 l)!)/((k + 2 l)^2 (k!)^2 l! (-k + l)!)),
       {l, 0, Floor[n/2]}], {n, 2, Nn}];
   {{RichardsonResum[Accumulate[varpi2tab], Nr], 
     RichardsonResum[Accumulate[varpi3tab], Nr]}, {RichardsonResum[
      Accumulate[varpi2dtab], Nr], 
     RichardsonResum[Accumulate[varpi3dtab], Nr]}, {RichardsonResum[
      Accumulate[varpi2ddtab], Nr], 
     RichardsonResum[Accumulate[varpi3ddtab], Nr]}}
    ];

PicardFuchsP[la_, u_] := (
  50 u la + 24 la^2 - 2016 u^3 la^2 + 
   u^2 (27 - 320 la^3) + 54 u^5 (-27 + 16 la^3) + 
   4 u^4 la (-783 + 224 la^3))/(
  u (9 u + 8 la) (u + la - 36 u^3 la - 
     8 u^2 la^2 + u^4 (-27 + 16 la^3)));

PicardFuchsQ[la_, u_] := (
  16 u la + 8 la^2 - 1860 u^3 la^2 + 
   u^2 (9 - 256 la^3) + 54 u^5 (-27 + 16 la^3) + 
   16 u^4 la (-189 + 64 la^3))/(
  u^2 (9 u + 8 la) (u + la - 36 u^3 la - 
     8 u^2 la^2 + u^4 (-27 + 16 la^3)));

PicardFuchs[f_, la_, u_] := 
  D[f, {u, 3}] + D[f, {u, 2}] PicardFuchsP[la, u] + 
   D[f, u] PicardFuchsQ[la, u];

PicardFuchsA[la_, u_] = {{0, 1, 0}, {0, 0, 
    1}, {0, -PicardFuchsQ[la, u], -PicardFuchsP[la, u]}};

SystemMatrix[{{t_, td_}, {dt_, dtd_}, {ddt_, ddtd_}}, la_, 
   u_] := Transpose@
   Prepend[Transpose[{{-((I (t + Log[u]))/(2 \[Pi])), 
       1/6 + (td + Log[u]^2 + Log[la]^2/8 - 
         1/4 (t + Log[u]) (8 Log[u] + Log[la]))/\[Pi]^2}, {-((
        I (dt + 1/u))/(2 \[Pi])), -((
        8 t - 4 dtd u + 8 (1 + dt u) Log[u] + Log[la] + 
         dt u Log[la])/(4 \[Pi]^2 u))}, {-((I (ddt - 1/u^2))/(
        2 \[Pi])), (-8 + 8 t - 16 dt u + 
        4 ddtd u^2 + (8 - 8 ddt u^2) Log[u] + Log[la] - 
        ddt u^2 Log[la])/(4 \[Pi]^2 u^2)}}], {1, 0, 0}];

Options[ComputeTransition] = {WorkingPrecision -> 60};
ComputeTransition[la_, upath_, Nn_ : 512, Nr_ : 5, 
   OptionsPattern[]] := Module[{sol, Psi, mon, SystMat, msg, tt = 0},
   Monitor[
    msg = "Integrating the connection...";
    sol = NDSolve[{
       Psi'[t] == upath'[t]*PicardFuchsA[la, upath[t]] . Psi[t],
       Psi[0] == IdentityMatrix[3]
       }, Psi, {t, 0, 1}, 
      WorkingPrecision -> OptionValue[WorkingPrecision], 
      MaxSteps -> 10^6, StepMonitor :> (tt = t)];
    mon = Psi[1] /. First[sol];
    msg = "Solving for the system matrix...";
    SystMat = 
     SystemMatrix[F1Series[la, upath[0], Nn, Nr], la, 
      upath[0]];
    LinearSolve[SystMat, mon . SystMat], 
    msg <> " " <> ToString[Floor[100 tt]] <> "%."]
   ];

G10FundamentalDomain = Graphics[{Table[{Circle[{n, 0}, 1, {Pi/3, 2 Pi/3}], 
     Line[{{n, 1}, {n, 2}}], Thin, 
     Line[{{n + 1/2, Sqrt[3]/2}, {n + 1/2, 2}}]}, {n, 1, 
     7}], {Circle[{2 + 1/2, 0}, 1/2, {0, Pi/2}], 
    Circle[{5 + 1/2, 0}, 1/2, {0, Pi/2}], 
    Line[{{2 + 1/2, 1/2}, {2 + 1/2, Sqrt[3]/2}}], Line[{{3, 0}, {3, 1}}], 
    Line[{{5 + 1/2, 1/2}, {5 + 1/2, Sqrt[3]/2}}], Line[{{6, 0}, {6, 1}}],
     Thin, Circle[{2, 0}, 1, {0, Pi/3}], 
    Circle[{5, 0}, 1, {0, Pi/3}]}, {Circle[{0, 0}, 1, {Pi/3, Pi/2}], 
    Line[{{0, 0}, {0, 2}}], Circle[{8, 0}, 1, {Pi/2, 2 Pi/3}], 
    Line[{{8, 0}, {8, 2}}], Thin, Line[{{1/2, Sqrt[3]/2}, {1/2, 2}}]}}];

G10RepeatDomain[{kmin_, kmax_}] := 
  Table[G10FundamentalDomain /. {l_Circle :> Translate[l, {8*k, 0}], 
     l_Line :> Translate[l, {8*k, 0}]}, {k, kmin, kmax}];


MirrorCurveJ2[el_, ur_] := (1 + 16 el^2 ur^4 - 
  8 el ur^2 (1 + 3 ur))^3/(el^3 ur^8 (1 + ur - 8 el ur^2 - 
   36 el ur^3 + el (-27 + 16 el) ur^4));
PicardFuchsP2[el_, ur_] := (24 + 50 ur + (27 - 320 el) ur^2 - 2016 el ur^3 + 
 4 el (-783 + 224 el) ur^4 + 
 54 el (-27 + 16 el) ur^5)/(ur (8 + 9 ur) (1 + ur - 8 el ur^2 - 
   36 el ur^3 + el (-27 + 16 el) ur^4));
PicardFuchsQ2[el_, ur_] := (8 + 16 ur + (9 - 256 el) ur^2 - 1860 el ur^3 + 
 16 el (-189 + 64 el) ur^4 + 
 54 el (-27 + 16 el) ur^5)/(ur^2 (8 + 9 ur) (1 + ur - 8 el ur^2 - 
   36 el ur^3 + el (-27 + 16 el) ur^4));
PicardFuchs2[f_, el_, ur_] := D[f, {ur, 3}] + PicardFuchsP2[el, ur] * D[f, {ur, 2}] + PicardFuchsQ2[el, ur] * D[f, {ur, 1}];
PicardFuchsA2[el_, ur_] := {{0, 1, 0}, {0, 0, 1}, {0, -PicardFuchsQ2[el, ur], -PicardFuchsP2[el, ur]}};

F1Series2[el_, ur_, Nn_ : 512, Nr_ : 10] := 
  Module[{varpi2tab, varpi3tab, varpi2dtab, varpi3dtab, varpi2ddtab, 
    varpi3ddtab, k, l}, varpi2tab = Table[ur^n  Sum[k = n - 2 l;
       n!/((l)! (k!)^2 (l - k)!) el^l/n, {l, 0, Floor[n/2]}], {n, 1, 
      Nn}];
   varpi3tab = Table[ur^n  Sum[k = n - 2 l;
       el^
         l (If[k > 
            l, ((-1)^(k + l) Gamma[k - l] Gamma[
               1 + k + 2 l])/(4 (k + 2 l) l! Gamma[
                1 + k]^2), (((k + 2 l)! (4 HarmonicNumber[k] + 
                 3 HarmonicNumber[l] + HarmonicNumber[-k + l] - 
                 8 HarmonicNumber[k + 2 l]))/(4 (k + 
                 2 l) (k!)^2 l! Gamma[1 - k + l]))] + (2 (k + 
                2 l)!)/((k + 2 l)^2 (k!)^2 l! (-k + l)!)), {l, 0, 
        Floor[n/2]}], {n, 1, Nn}];
   varpi2dtab = Table[n ur^(n - 1) Sum[k = n - 2 l;
       n!/((l)! (k!)^2 (l - k)!) el^l/n, {l, 0, Floor[n/2]}], {n, 1, 
      Nn}];
   varpi3dtab = Table[n ur^(n - 1) Sum[k = n - 2 l;
       el^
         l (If[k > 
            l, ((-1)^(k + l) Gamma[k - l] Gamma[
               1 + k + 2 l])/(4 (k + 2 l) l! Gamma[
                1 + k]^2), (((k + 2 l)! (4 HarmonicNumber[k] + 
                 3 HarmonicNumber[l] + HarmonicNumber[-k + l] - 
                 8 HarmonicNumber[k + 2 l]))/(4 (k + 
                 2 l) (k!)^2 l! Gamma[1 - k + l]))] + (2 (k + 
                2 l)!)/((k + 2 l)^2 (k!)^2 l! (-k + l)!)), {l, 0, 
        Floor[n/2]}], {n, 1, Nn}];
   varpi2ddtab = Table[n (n - 1) ur^(n - 2) Sum[k = n - 2 l;
       n!/((l)! (k!)^2 (l - k)!) el^l/n, {l, 0, Floor[n/2]}], {n, 2, 
      Nn}];
   varpi3ddtab = Table[n (n - 1) ur^(n - 2) Sum[k = n - 2 l;
       el^
         l (If[k > 
            l, ((-1)^(k + l) Gamma[k - l] Gamma[
               1 + k + 2 l])/(4 (k + 2 l) l! Gamma[
                1 + k]^2), (((k + 2 l)! (4 HarmonicNumber[k] + 
                 3 HarmonicNumber[l] + HarmonicNumber[-k + l] - 
                 8 HarmonicNumber[k + 2 l]))/(4 (k + 
                 2 l) (k!)^2 l! Gamma[1 - k + l]))] + (2 (k + 
                2 l)!)/((k + 2 l)^2 (k!)^2 l! (-k + l)!)), {l, 0, 
        Floor[n/2]}], {n, 2, Nn}];
   {{RichardsonResum[Accumulate[varpi2tab], Nr], 
     RichardsonResum[Accumulate[varpi3tab], Nr]}, {RichardsonResum[
      Accumulate[varpi2dtab], Nr], 
     RichardsonResum[Accumulate[varpi3dtab], Nr]}, {RichardsonResum[
      Accumulate[varpi2ddtab], Nr], 
     RichardsonResum[Accumulate[varpi3ddtab], Nr]}}];

SystemMatrix2[{{v1_, v2_}, {v1p_, v2p_}, {v1pp_, v2pp_}}, el_, 
   ur_] := {{1, (Log[ur] + 1/3 Log[el] + v1)/(2 Pi I), (
    4 (\[Pi]^2 + 6 v2) - 
     3 (Log[el]^2 + 6 Log[el] (v1 + Log[ur]) + 
        8 Log[ur] (2 v1 + Log[ur])))/(
    24 \[Pi]^2)}, {0, -((I (1/ur + v1p))/(2 \[Pi])), -((
     8 v1 - 4 ur v2p + 3 Log[el] + 8 Log[ur] + 
      ur v1p (3 Log[el] + 8 Log[ur]))/(4 \[Pi]^2 ur))}, {0, -((
     I (-(1/ur^2) + v1pp))/(2 \[Pi])), (-8 + 8 v1 - 16 ur v1p + 
     4 ur^2 v2pp + 3 Log[el] - 3 ur^2 v1pp Log[el] + 8 Log[ur] - 
     8 ur^2 v1pp Log[ur])/(4 \[Pi]^2 ur^2)}};


End[]; (* `Private` *)

EndPackage[];
