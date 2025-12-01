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

BeginPackage["LocalF1`"];


(* ===================== *)
(*   Public symbols      *)
(* ===================== *)

MirrorCurve::usage = "MirrorCurveF1[x, y, m, u] is the affine mirror-curve equation of local F1 \
in the variables x, y and parameters m, u.";

MirrorCurvef3::usage = "TODO";

MirrorCurveg2::usage = "TODO";

MirrorCurveg3::usage = "TODO";

MirrorCurveDelta::usage = "TODO";

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

(* Usage messages for exported functions *)

Begin["`Private`"];

(* ===================== *)
(*   Internal helpers    *)
(* ===================== *)



(* ===================== *)
(*   Public definitions  *)
(* ===================== *)

MirrorCurve[x_,y_,m_,u_] := x + y + 1/x/y + m/x - 1/u;

MirrorCurvef3[X_, m_, u_] := 5832 u^6 - 1728 m^3 u^6 - 
 432 m^2 u^4 (-3 + X) + (3 - 2 X)^2 (3 + X) + 324 u^3 (-3 + 2 X) + 
 108 m u^2 (-3 + 36 u^3 + 2 X);

MirrorCurveg2[m_, u_] := 27 (1 - 8 m u^2 - 24 u^3 + 16 m^2 u^4);

MirrorCurveg3[m_, u_] := 27 (-1 + 12 m u^2 + 36 u^3 - 48 m^2 u^4 - 144 m u^5 - 216 u^6 + 
   64 m^3 u^6);

MirrorCurveDelta[m_, u_] := m + u - 8 m^2 u^2 - 36 m u^3 - 27 u^4 + 16 m^3 u^4;

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
   all = Join[fn, {{#[[1]], -#[[2]]} & /@ fd} // Flatten[#, 1] &];
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



End[]; (* `Private` *)

EndPackage[];