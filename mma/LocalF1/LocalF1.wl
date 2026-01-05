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

repCh::usage = "Expands Ch[p, q] Ch[p, q][1] and GV[p, q, n]";

Euler::usage = "TODO";

DSZ::usage = "TODO";

ToFB::usage = "TODO";

ToHC::usage = "TODO";

GiesekerSlope::usage = "TODO";

BogomolovDiscriminant::usage = "TODO";

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

InitialPosition::usage = "InitialPosition[{r, dF, dB, ch2}, m] calculates the initial s value of the ray when it intersects the real axis. Only works for real m.";

CollDBridgeland::usage = "";
CollDArcaraII::usage = "";
CollDArcaraI::usage = "";

CollBridgeland::usage = "";
CollArcaraII::usage = "";
CollArcaraI::usage = "";

(* XY plane scattering *)
XYRay::usage = "";
XYParabolaY::usage = "";
InitialPositionXY::usage = "InitialPositionXY[gam, m] gives the initial position of ray gam in the x y plane (point with 0 cost function). Only works for real m.";

InitialRaysFromColl::usage = "InitialRaysFromColl[Coll, xmin, xmax, m] produces initial rays for ConstructLVDiagram by taking all translates of Coll that start between x = xmin and x = xmax. Only works for real m.";

(* from F0Scattering.m *)

GCD1::usage = "";
SameHalfPlaneQ::usage = "SameHalfPlaneQ[Zlist] gives True if all elements of Zlist are in a common half plane";

QuiverDomain::usage = "QuiverDomain[Coll, psi, m] plots the region where the LV central charges Z of Coll have Re[e^{-I psi} Z] < 0, and the region where they are in same half-plane."; 

ChernToCh::usage = "";

ExtFromStrong::usage = "ExtFromStrong[Coll] computes the Chern vectors of the objects in the Ext collection dual to the given strong collection Coll"
StrongFromExt::usage = "StrongFromExt[Coll] computes the Chern vectors of the objects in the strong collection dual to the given Ext collection Coll";

TreeFromListRays::usage = "Remains to be implemented";
LVTreesFromListRays::usage = "Remains to be implemented";

KroneckerDims::usage = "KroneckerDims[m, Nn] gives the list of populated dimension vectors {n1,n2} for Kronecker quiver with m arrows, with (n1,n2) coprime and 0<=n1,n2<=Nn"; 
IntersectRaysNoTest::usage = "IntersectRays[{r, dF, dB, ch2}, {rr, ddB, ddF, cch2}, z, zz, m] returns intersection point (x,y) of two rays if the intersection point lies strictly upward from z and z', or {} otherwise, without testing non-vanishing of DSZ product";  
CostPhi::usage = "CostPhi[{r, dF, dB, ch2}, s, mu] gives the cost function \\phi_s(\\gamma) = dB + 2 dF - r (mu + 8 s)";
ConstructLVDiagram::usage = "ConstructLVDiagram[smin, smax, phimax, Nm, m, ListRays] constructs the LV scattering diagram for F1 with initial rays in the interval [smin,smax], cost function up to phimax, scattering products with n1 + n2 <= Nn at each intersection; m is assumed to be real; The output consists of a list of  { charge, {x,y}, parent1, parent2, n1, n2 }; If ListRays is not empty, then uses it as initial rays.";
ConstructLVDiagramOpt::usage = "";

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

BogomolovDiscriminant[{r_, dF_, dB_, ch2_}] := -(dB^2/(2 r^2)) + (dB dF)/r^2 - ch2/r;

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

MonodromyFromSphericalTwist[mF_, mB_] := {{1, 0, 0, 0}, {0, 1, 0, 0}, {1/2 mB (mB - 2 mF), -mB + mF,
     1 + mB + 2 mF, -1}, {1/2 mB (mB^2 - 4 mF^2), -mB^2 - mB mF + 
     2 mF^2, (mB + 2 mF)^2, 1 - mB - 2 mF}};

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

MonodromyOnTau[M_, tau_] := (M[[4, 3]] + 8 tau M[[4, 4]])/(8 M[[3, 3]] + 
     64 tau M[[3, 4]]);  


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

InitialPosition[{r_, dF_, dB_, ch2_}, m_: 0] := (dB + 2 dF - m r - Sqrt[
 dB^2 - 16 ch2 r + 2 dB (2 dF - 9 m r) + (2 dF + 3 m r)^2])/(8 r);

CollDBridgeland = {{1, 0, 0, 0}, {1, 1, 0, 0}, {1, 1, 1, 1/2}, {1, 2, 1, 3/2}};

CollDArcaraII = {{1, 0, 0, 0}, {1, 1, 0, 0}, {1, 1, 1, 1/2}, {2, 1, 1, -(1/2)}};

CollDArcaraI = {{1, 0, 0, 0}, {1, 0, 1, -(1/2)}, {1, 1, 1, 1/2}, {2, 1, 1, -(1/2)}};

CollBridgeland = {{1, 0, 0, 0}, {-1, 1, 0, 0}, {-1, 0, 1, 1/2}, {1, -1, -1, 1/2}};

CollArcaraII = {{1, 0, 0, 0}, {0, 0, -1, 1/2}, {1, -2, -1, 3/2}, {-1, 1, 1, -(1/2)}};

CollArcaraI = {{1, 0, 0, 0}, {0, 0, 1, -(1/2)}, {1, -2, -1, 3/2}, {-1, 1, 0, 0}};


(* XY plane scattering *)

XYRay[{r_, dF_, dB_, ch2_}, {x_, y_}, mu_] := 
  r y + (2 dF + dB) x + (dF - dB) mu - ch2;

XYParabolaY[x_,mu_,m2_,psi_] := 1/4 (18 m2^2 + mu^2 - 2 mu x - 
    8 x^2 + (mu^2 - 2 mu x - 8 x^2) Cos[2 psi]) Sec[psi]^2;


(* from F0Scattering.m *)
SameHalfPlaneQ[{}] := True;
SameHalfPlaneQ[Zlist_List] := 
  If[AnyTrue[Zlist, # == 0 &], 
   False, -Subtract @@ MinMax[Arg[Zlist/Zlist[[1]]]] < Pi];

Options[QuiverDomain] = {PlotStyle -> LightBlue};
QuiverDomain[Coll_, psi_, m_, 
   OptionsPattern[]] := {RegionPlot[(And @@ 
       Table[Re[Exp[-I psi] ZLV[Coll[[i]], {s, t}, m]] < 0, {i, 
         Length[Coll]}]) && t > 0, {s, -1.5, 1.5}, {t, 0, 1}, 
    PlotPoints -> 100, AspectRatio -> 1, 
    PlotStyle -> Flatten[{OptionValue[PlotStyle], Opacity[.5]}]], 
   RegionPlot[
    SameHalfPlaneQ[
     Table[ZLV[Coll[[i]], {s, t}, m], {i, Length[Coll]}]], {s, -1.5, 
     1.5}, {t, 0, 1}, PlotPoints -> 100, AspectRatio -> 1, 
    BoundaryStyle -> Directive[Dashed], 
    PlotStyle -> Flatten[{OptionValue[PlotStyle], Opacity[.3]}]]};

ExtFromStrong[Coll_]:=Module[{S,Si},
   S=Table[Euler[Coll[[i]],Coll[[j]]],{i,Length[Coll]},{j,Length[Coll]}];
   Si=Inverse[Transpose[S]];
   Si . Coll
   ]

StrongFromExt[Coll_]:=Module[{S,Si},
   S=Table[Euler[Coll[[j]],Coll[[i]]],{i,Length[Coll]},{j,Length[Coll]}];
   Si=Inverse[Transpose[S]];
   Si . Coll
   ]

TreeFromListRays[ListRays_,k_]:=If[ListRays[[k,3]]==0,ListRays[[k,1]],{ListRays[[k,5]]TreeFromListRays[ListRays,ListRays[[k,3]]],ListRays[[k,6]]TreeFromListRays[ListRays,ListRays[[k,4]]]}];

GCD1[{r_,dF_,dB_,ch2_}]:=Module[{d},d=GCD[r,dF,dB];
If[EvenQ[(dB-2ch2)/d],d,If[EvenQ[d],d/2,1]]];

LVTreesFromListRays[ListRays_,{r_,dF_,dB_,ch2_},m_]:=Module[{Lipos,div,LiTrees},
   div=Divisors[GCD1[{r,dF,dB,ch2}]];
   Lipos=Flatten[Join[Table[Position[ListRays,{r,dF,dB,ch2}/k],{k,div}]],1];
   If[Lipos=={},
   Print["No such dimension vector in the list"],
   LiTrees=(GCD1[{r,d1,d2,ch2}]/GCD1[ListRays[[#,1]]])TreeFromListRays[ListRays,#]&/@First[Transpose[Lipos]]
   (* ScattSort[DeleteDuplicatesBy[SortBy[LiTrees,Length[TreeConstituents[#]]&],ScattGraph[#,m]&],m] *)
]];

ChernToCh[f_]:=f/.{{r_Integer,dF_Integer,dB_Integer,ch2_}:>If[r==0,GV[dF,dB,ch2],If[BogomolovDiscriminant[{r,dF,dB,ch2}]==0,If[r>0,r Ch[dF/r,dB/r],-r Ch[dF/r,dB/r][1]],{r,dF,dB,ch2}]]};	

KroneckerDims[m_Integer?NonNegative, Nn_Integer?NonNegative] := KroneckerDims[m,Nn] = Module[{Ta={}},
   Do[If[m n1 n2-n1^2-n2^2+1>=0&&GCD[n1,n2]==1,AppendTo[Ta,{n1,n2}]],{n1,0,Nn},{n2,0,Nn-n1}];Drop[Ta,2]];

(* TODO: can be simplified *)
IntersectRaysNoTest[{r_, dF_, dB_, ch2_}, {rr_, ddF_, ddB_, cch2_}, 
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
      CostPhi[{rr, ddF, ddB, cch2}, zz[[1]], mu], zi, {}]];

CostPhi[{r_,dF_,dB_,ch2_},s_,m_]:=dB + 2 dF - r (m + 8 s);

InitialPositionXY[{r_, dF_, dB_, ch2_}, m_] := {(dB + 2 dF - m r)/(
   8 r), -((dB^2 + 4 dB dF + 4 dF^2 - 8 ch2 r - 9 dB m r + 6 dF m r)/(
    8 r^2))};

InitialRaysFromColl[Coll_, xmin_, xmax_, m_] := Module[{x, y},
   Flatten[
    Table[
     {x, y} = InitialPositionXY[gam, m];
     Table[
      With[{newgam = SpectralFlow[gam, {3 k, 2 k}]},
       {sign*newgam, InitialPositionXY[newgam, m], 0, 0, 0, 0, 0}
       ], {k, 
       Range[Ceiling[xmin - x], Floor[xmax - x]]}, {sign, {-1, 
        1}}], {gam, Coll}], 2]
   ];


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
        If[! TrueQ @ Lookup[seen, {i, j}, False],
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


End[]; (* `Private` *)

EndPackage[];
