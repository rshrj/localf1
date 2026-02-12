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

RichardsonResum::usage =
"RichardsonResum[partSums, n] returns the n-th order Richardson transform of the list of partial sums partSums. \
Typical use: RichardsonResum[Accumulate[data], n].";

MirrorCurve::usage =
"MirrorCurve[x, y, m, u] gives the affine mirror-curve equation of local F1 in variables x, y and parameters m, u \
(i.e. an expression whose vanishing defines the curve).";

MirrorCurvef3::usage =
"MirrorCurvef3[X, m, u] gives the cubic f3(X) on the right-hand side of the Weierstrass form \
Y^2 = f3(X) = 4 X^3 - g2 X - g3 for the local F1 mirror curve.";

MirrorCurveg2::usage =
"MirrorCurveg2[m, u] gives the Weierstrass invariant g2 (weight-4 Eisenstein invariant) of the local F1 mirror curve.";

MirrorCurveg3::usage =
"MirrorCurveg3[m, u] gives the Weierstrass invariant g3 (weight-6 Eisenstein invariant) of the local F1 mirror curve.";

MirrorCurveDelta::usage =
"MirrorCurveDelta[m, u] gives the order-4 factor in the discriminant of the local F1 mirror curve (as a polynomial in u). \
Its roots locate the discriminant locus in the (m,u)-slice.";

MirrorCurveJ::usage =
"MirrorCurveJ[m, u] gives the j-invariant of the local F1 mirror curve as a rational function of (m,u).";

URoot::usage =
"URoot[la, k] gives the k-th Root object u solving MirrorCurveDelta[la, u] == 0.";

URoots::usage =
"URoots[m] returns the list of numerical roots u of MirrorCurveDelta[m, u] == 0.";

PlotURoots::usage =
"PlotURoots[m] plots the roots URoots[m] in the complex u-plane and labels them by their absolute values.";

NodalCurveParam::usage =
"NodalCurveParam[X, Y, t, X0] returns replacement rules {X->..., Y->...} giving a nodal (rational) parametrization \
of the Weierstrass curve in terms of parameter t and node location X0.";

NodalCurveX0::usage =
"NodalCurveX0[m, u] gives the X0-parameter used in the nodal parametrization at the discriminant locus.";

NodalCurveX0tot1::usage =
"NodalCurveX0tot1 is a replacement rule {X0 -> ...} expressing X0 in terms of t1.";

WeierstrassChangeOfVariables::usage =
"WeierstrassChangeOfVariables[{x,y}->{X,Y}, m, u] gives replacement rules for {x,y} in terms of {X,Y} \
implementing the change of variables to Weierstrass coordinates.";

WeierstrassChangeOfVariablesInverse::usage =
"WeierstrassChangeOfVariablesInverse[{X,Y}->{x,y}, m, u] gives replacement rules for {X,Y} in terms of {x,y} \
implementing the inverse change of variables from Weierstrass coordinates.";

LiEval::usage =
"LiEval is a list of replacement rules that numerically evaluates special functions such as Li2[z] and BW[z] \
(e.g. expr /. LiEval).";

RationalData::usage =
"RationalData[expr, var] returns {exponents, roots} for the linear factors in var appearing in the rational expression expr. \
Exponents are positive for numerator factors and negative for denominator factors.";

NodalCurvebFromt1::usage =
"NodalCurvebFromt1[t1] gives the parameter b as a function of t1.";

NodalCurvebFrommu::usage =
"NodalCurvebFrommu[m, u] gives the parameter b as a function of (m,u).";

KerrDoranP::usage =
"KerrDoranP[b] evaluates the Kerr--Doran quantity P(b) built from Bloch--Wigner dilogarithms and the correction term.";

KerrDoranQ::usage =
"KerrDoranQ[b] evaluates the Kerr--Doran quantity Q(b) built from logarithms and dilogarithms.";

(* Boris *)

ToFund::usage =
"ToFund[tau] maps a point tau in the upper half-plane to the standard SL(2,Z) fundamental domain \
(by repeated application of tau -> tau - n to enforce -1/2 <= Re[tau] < 1/2 and, if |tau|<1, the inversion tau -> -1/tau). \
The result satisfies -1/2 <= Re[tau] < 1/2 and Abs[tau] >= 1 (up to boundary conventions).";

InverseJ::usage =
"InverseJ[j, prec] returns a numerical value of tau in the SL(2,Z) fundamental domain such that j(tau) = j, \
computed with numerical precision prec (default: MachinePrecision).\n\n\
The inversion is performed via the classical hypergeometric uniformization: an auxiliary parameter a is obtained from \
4 a (1 - a) == 1728/j, and tau is evaluated as\n\
  tau = i * 2F1(1/6,5/6;1;1-a) / 2F1(1/6,5/6;1;a),\n\
followed by reduction to the fundamental domain using ToFund.\n\n\
Note: The inverse of the j-invariant is globally multivalued; this routine returns a single branch determined by the \
chosen solution for a and the fundamental-domain reduction.";

Mon3ToMon4::usage =
"Mon3ToMon4[m0, M] embeds a 3x3 matrix M into a 4x4 matrix with a fixed 2x2 identity block in the upper-left, \
using the complex parameter m0 to separate real and imaginary parts of the first-column data.\n\n\
If Im[m0] != 0, the first-column entries are adjusted so that the dependence on m0 is made explicit and the remaining \
coefficients are extracted from imaginary parts. The result is rounded to 10^-3.\n\
If Im[m0] == 0, a fallback embedding is used (with a placeholder 999 in the (3,2) and (4,2) entries), and the result is \
again rounded to 10^-3.";

tauB::usage =
"tauB[\\[ScriptL], prec] returns the modular parameter tau (in the SL(2,Z) fundamental domain) associated to the parameter \
\\[ScriptL], computed with numerical precision prec (default: MachinePrecision).\n\n\
It evaluates the j-invariant\n\
  j(\\[ScriptL]) = (27 + 256 \\[ScriptL]) (243 + 256 \\[ScriptL])^3 / (2^24 \\[ScriptL]^3)\n\
and then applies InverseJ[j(\\[ScriptL]), prec].";

Mon4ToMon3::usage =
"Mon4ToMon3[m0, M] converts a 4x4 monodromy matrix M to a 3x3 monodromy matrix by restricting to the period subspace \
and substituting the mass parameter m = m0 via the combination M[[i,1]] + m0 M[[i,2]] in the first column (rows 3 and 4).";

FindPeriodCombination::usage =
"FindPeriodCombination[m0, P0, x, m, Q] finds an integer linear relation among {1, m0, P0, x}. \
If it finds integers {v1,v2,v3,v4} with v1 + v2 m0 + v3 P0 + v4 x == 0 and v4 != 0, it returns the corresponding \
symbolic expression -(v1 + v2 m + v3 Q)/v4, i.e. x written as a linear combination of the symbols m and Q with the \
same integer coefficients determined at (m0,P0,x). If v4 == 0, it returns {}.";

KerrDoranQ1::usage =
"KerrDoranQ1[b] For Re[b]>0 TODO";

KerrDoranFix0::usage =
"KerrDoranFix0[b] returns a piecewise-constant correction term (a multiple of \\[Pi]) used to fix branch choices in the Kerr--Doran prescription. \
It is implemented via HeavisideTheta conditions in the complex b-plane and is nonzero only in the left half-plane (Re[b]<0), \
with additional subregion structure depending on |b| and angular wedges.";

KerrDoranFix1::usage =
"KerrDoranFix1[b] returns a piecewise-constant integer multiple of I used as the coefficient of Log[b] in the Kerr--Doran branch-fixing prescription. \
It is supported in the region Re[b] < -1/2 and depends on the sign of Im[b], with wedge exclusions implemented via HeavisideTheta conditions.";

KerrDoranQ1Fix::usage =
"KerrDoranQ1Fix[b] gives the branch-corrected Kerr--Doran quantity Q1(b), defined as KerrDoranQ1[b] plus the constant correction KerrDoranFix0[b] \
and the logarithmic correction KerrDoranFix1[b] Log[b].";

KerrDoranFormula::usage =
"KerrDoranFormula[avec, bvec, dvec, evec, A, B] evaluates the Kerr--Doran dilogarithmic combination built from \
the lists avec, bvec (complex parameters) with weights dvec, evec and scalars A, B. \
It returns a linear combination of Li2[avec[[j]]/bvec[[k]]] plus logarithmic correction terms, with diagonal terms (avec[[j]]==bvec[[k]]) omitted.";

KerrDoranFormulaExact::usage = "KerrDoranFormulaExact[avec, bvec, dvec, evec, A, B] TODO";

KerrDoranElFromb::usage =
"KerrDoranElFromb[b] gives the parameter el as a rational function of the Kerr--Doran parameter b.";

KerrDoranUrFromb::usage =
"KerrDoranUrFromb[b] gives the parameter ur as a rational function of the Kerr--Doran parameter b.";

KerrDoranElUrFromb::usage = "KerrDoranElUrFromb[b] returns {el, ur} as rational functions of the Kerr--Doran parameter b.";

KerrDoranbFromEl::usage =
"KerrDoranbFromEl[el] solves KerrDoranElFromb[b] == el for b, selecting the branch with Abs[b] > 1 (or Abs[b] == 1 and Im[b] >= 0).";

AnalyticContinuationOfRoots::usage =
"AnalyticContinuationOfRoots[poly, x, t, t0, t1] numerically tracks the roots in x of the polynomial equation poly==0 \
as the parameter t is continuously deformed from t0 to t1.\n\n\
The function first solves poly /. t->t0 for x to obtain the initial root set, then follows these roots along a discrete \
path in t by repeatedly solving poly at intermediate t-values and matching roots between consecutive steps by proximity.\n\n\
Output:\n\
Returns an Association with keys:\n\
  \"tValues\"      -> list of sampled parameter values along the path,\n\
  \"rootsPaths\"   -> list of root trajectories; each element is the sequence of one root tracked across \"tValues\",\n\
  \"initialRoots\" -> roots at t=t0,\n\
  \"finalRoots\"   -> roots at t=t1.\n\n\
Options:\n\
  Steps -> n (default 100)\n\
    Number of sample points used along the t-path (including endpoints). Larger values typically improve tracking \
    when roots move quickly or approach each other.\n\
  Path -> \"Straight\" (default) or an explicit list\n\
    If \"Straight\", uses a linear interpolation path t(s)=t0+(t1-t0)s sampled at Steps points.\n\
    If a list is given, it is interpreted as an explicit sequence of t-values to traverse; the first element should \
    equal t0 and the last should equal t1.\n\n\
Algorithm notes:\n\
At each step i, roots are computed by NSolve[poly /. t->tValues[[i]], x]. When the number of roots matches the previous \
step, a distance matrix |x_j^{(i-1)}-x_k^{(i)}| is formed and a greedy nearest-neighbor assignment is used to permute \
the new roots so that each tracked trajectory remains continuous.\n\n\
Typical use cases:\n\
  • Following algebraic branches across parameter space.\n\
  • Detecting monodromy/permutations of roots after traversing loops (compare \"initialRoots\" and \"finalRoots\").\n\
  • Visualizing root motion by plotting the entries of \"rootsPaths\" against \"tValues\".\n\n\
Caveats:\n\
  • If two roots become very close between steps (near a discriminant point), the greedy matching can swap branches; \
    increase Steps or choose a path that avoids the collision.\n\
  • If NSolve returns a different number/order of solutions between steps (e.g. due to numerical issues), tracking may fail \
    or produce discontinuities.\n\
  • This routine assumes poly is a polynomial in x (possibly with parameter t) and relies on numerical root finding.";

BW::usage =
"BW[z] is the Bloch--Wigner dilogarithm of z.";

Li2::usage =
"Li2[z] is the classical (Spence) dilogarithm Li_2(z).";


(* Large volume stuff *)

Ch::usage =
"Ch[p, q] represents the Chern vector of the sheaf O(p F + q B) (as a symbolic object used by repCh and related routines).";

GV::usage =
"GV[p, q, n] represents the Chern vector of a rank-zero object with (dF,dB,ch2) = (p,q,n) (as a symbolic object used by repCh).";

y::usage =
"y is the refinement parameter used in refined indices (e.g. in Om / Coulomb-branch formulas).";

repCh::usage =
"repCh is a list of replacement rules expanding symbolic charges such as Ch[p,q], Ch[p,q][1], and GV[p,q,n] into vectors {r,dF,dB,ch2}.";

Euler::usage =
"Euler[gam, gamp] gives the Euler pairing \\[Chi](gam, gamp) for charges gam and gamp written as {r,dF,dB,ch2}.";

DSZ::usage =
"DSZ[gam, gamp] gives the Dirac--Schwinger--Zwanziger pairing <gam, gamp> for charges gam and gamp written as {r,dF,dB,ch2}.";

ToFB::usage =
"ToFB[gam] converts charges from (H,C) coordinates to (F,B) coordinates.";

ToHC::usage =
"ToHC[gam] converts charges from (F,B) coordinates to (H,C) coordinates.";

GiesekerSlope::usage =
"GiesekerSlope[{r, dF, dB, ch2}, M] returns the Gieseker slope \\[Mu]_M for the given charge at polarization parameter M.";

GiesekerDim::usage =
"GiesekerDim[{r, dF, dB, ch2}] returns the expected dimension (virtual dimension) associated to the given charge.";

BogomolovDiscriminant::usage =
"BogomolovDiscriminant[{r, dF, dB, ch2}] gives the Bogomolov discriminant \\[Delta] of the charge.";

SecondChernClass::usage =
"SecondChernClass[gam] computes the second Chern class c2 from the charge vector gam = {r,dF,dB,ch2}.";

ZLV::usage =
"ZLV[gam, T, m] gives the large-volume central charge for local F1, with gam = {r,dF,dB,ch2}, T complex and m a parameter. \
ZLV[gam, {s,t}, m] uses T = s + I t.";

SpectralFlow::usage =
"SpectralFlow[gam, {mF, mB}] applies the spectral-flow transform by integers (mF,mB) to the charge gam. \
SpectralFlow[charges, {mF,mB}] maps over a list of charges.";

Sigma::usage =
"Sigma is the matrix relating the charge basis to the period basis used in monodromy actions (so that ZLV can be written via Sigma and periods).";

MLV::usage =
"MLV is the large-volume monodromy matrix.";

MSF::usage =
"MSF[mF, mB] gives the monodromy matrix implementing spectral flow by Ch[mF, mB].";

MonodromyFromSphericalTwist::usage =
"MonodromyFromSphericalTwist[gam] returns the monodromy matrix of the spherical twist by a spherical object of charge gam.";

M1p::usage =
"M1p is the conifold monodromy matrix at the special \\[Lambda] point (F1Global.pdf, eq. (5.59)).";

M2p::usage =
"M2p is the conifold monodromy matrix at the special \\[Lambda] point (F1Global.pdf, eq. (5.59)).";

MonodromyOnCharge::usage =
"MonodromyOnCharge[M, gam] applies the monodromy matrix (or a list of matrices) M to a charge vector gam using Sigma.";

MonodromyOnTau::usage =
"MonodromyOnTau[M, tau] gives the fractional-linear action of monodromy matrix M on \\[Tau].";


(* Large volume scattering *)

Rays::usage =
"Rays[gam, t, psi, m] computes s = Re[T] for the ray defined by Re(Exp[-I psi] ZLV[gam, T, m]) == 0 at fixed t = Im[T].";

Wall::usage =
"Wall[gam, gamp, {s,t}, m] gives the wall equation of marginal stability between charges gam and gamp in the (s,t)-plane. \
Only works for real m.";

WallRadius::usage =
"WallRadius[gam, gamp, m] returns the squared radius of the wall circle in the (s,t)-plane. Only works for real m.";

WallCircle::usage =
"WallCircle[gam, gamp, m] returns a Circle[...] object for the wall in the (s,t)-plane (upper semicircle when present). \
Only works for real m.";

DiscF1::usage =
"DiscF1[gam, m] computes the ray discriminant of gam used to determine whether the corresponding ray intersects the real axis. \
Only works for real m.";

InitialPosition::usage =
"InitialPosition[gam, m] gives the initial s-position where the ray of charge gam meets the real T-axis (t=0). \
Only works for real m.";

CollDBridgeland::usage =
"CollDBridgeland is a list of Chern vectors for the dual (Ext) collection associated with the Bridgeland--Stern--Perling collection.";

CollDArcaraII::usage =
"CollDArcaraII is a list of Chern vectors for the dual (Ext) collection associated with the Arcara--Miles II collection.";

CollDArcaraI::usage =
"CollDArcaraI is a list of Chern vectors for the dual (Ext) collection associated with the Arcara--Miles I collection.";

CollBridgeland::usage =
"CollBridgeland is a list of Chern vectors for the Bridgeland--Stern--Perling strong exceptional collection.";

CollArcaraII::usage =
"CollArcaraII is a list of Chern vectors for the Arcara--Miles II strong exceptional collection.";

CollArcaraI::usage =
"CollArcaraI is a list of Chern vectors for the Arcara--Miles I strong exceptional collection.";

GamToString::usage =
"GamToString is a list of replacement rules mapping charges (e.g. Ch[p,q]) to compact strings used for plot labels.";

GamCaption::usage =
"GamCaption[gam, t, psi, m] produces a Text[...] label for charge gam placed near the ray at parameter value t.";

ListSpecFlows::usage =
"ListSpecFlows[coll, psi, m, L, xmin, xmax, tmax] lists spectral-flow parameters {mF,mB} with -L <= mF,mB <= L \
such that the translated collection has a non-empty quiver domain intersecting xmin < Re[T] < xmax and 0 < Im[T] < tmax.";


(* XY plane scattering *)

XYRay::usage =
"XYRay[gam, {x,y}, mu] gives the ray equation in XY-coordinates for the charge gam.";

XYParabolaY::usage =
"XYParabolaY[x, mu, m2, psi] gives the parabola y(x) used in the XY-to-(s,t) mapping at parameters (mu,m2,psi).";

InitialPositionXY::usage =
"InitialPositionXY[gam, m] gives the initial point {x,y} of the ray for charge gam in the XY-plane, the point where the cost function vanishes. Only works for real m.";

XYTost::usage =
"XYTost[{x,y}, psi, m] converts a point in XY-coordinates to {s,t} coordinates. Only works for real m.";

XYRaySegment::usage =
"XYRaySegment[rayData, mu, L, styleMap, defaultStyle] draws the line segment for a ray of length L with styling determined by styleMap.";

CollPlotData::usage =
"CollPlotData is a list of associations encoding collections and plotting data used by ConstructDomain.";

RayStyleMap::usage =
"RayStyleMap is the default depth-to-style map used for ray plotting.";

ConstructDomain::usage =
"ConstructDomain[plotAssociation, mu, L] constructs the quiver domain region(s) in the XY-plane for the collection specified by plotAssociation.";

NaiveIntersection::usage =
"NaiveIntersection[gam, gamp, mu] returns the intersection point {x,y} of the two XY-ray lines associated with gam and gamp (no future-direction checks).";

CollInitial::usage =
"CollInitial is a list of three initial charge collections identified for a fixed m-slice. CollInitial[[1]] for 0 < m < 1/3, CollInitial[[2]] for 1/3 < m < 2/3, CollInitial[[3]] for -1/3 < m < 0.";

InitialRaysFromColl::usage =
"InitialRaysFromColl[coll, xmin, xmax, m] generates initial rays for ConstructLVDiagram by translating coll so rays start between x=xmin and x=xmax. \
Only works for real m.";

CollInequalities::usage =
"CollInequalities[coll, {x,y}, mu] returns the list of inequalities defining the quiver domain in XY-coordinates for all gam in coll.";


(* from F0Scattering.m *)

GCD1::usage =
"GCD1[gam] returns the modified gcd used in the scattering algorithm for the charge gam.";

SameHalfPlaneQ::usage =
"SameHalfPlaneQ[zlist] returns True iff all complex numbers in zlist lie in a common open half-plane.";

QuiverDomain::usage =
"QuiverDomain[coll, psi, m] plots the region in the (s,t)-plane where all central charges ZLV for coll satisfy Re(Exp[-I psi] Z) < 0.";

ReZLV::usage = "ReZLV[{r, df, db, ch2}, {s, t}, m] in the real part of ZLV[{r, df, db, ch2}, {s, t}, m].";
QConditions::usage = "QConditions[Coll, m, {s, t}] define the domain in (s, t) plane where the collection Coll is valid.";
QDomain::usage = "QDomain[Coll, m, tooltip_label] creates a region plot in (s, t) plane where the collection Coll is valid. It generates a tooltip with the argument 'tooltip_label' activated by hovering on the region. It takes any option that RegionPlot takes.";

ChernToCh::usage =
"ChernToCh[expr] rewrites charge vectors {r,dF,dB,ch2} into symbolic Ch[...] or GV[...] when applicable (e.g. at vanishing discriminant).";

ExtFromStrong::usage =
"ExtFromStrong[coll] computes the Chern vectors of the Ext (dual) collection associated to the strong collection coll.";

StrongFromExt::usage =
"StrongFromExt[coll] computes the Chern vectors of the strong collection dual to the Ext collection coll.";

MutateCollection::usage =
"MutateCollection[coll, klist] applies successive mutations to the collection coll. \
klist is a list of {node, sign}, with sign = +1 for right mutation and sign = -1 for left mutation.";

ScattCheck::usage =
"ScattCheck[tree, m] checks consistency of a scattering tree. \
Returns {totalCharge, {x,y}} if consistent, otherwise {totalCharge, {}}.";

ScattSort::usage =
"ScattSort[treeList, m] sorts scattering trees by decreasing radius (using the root position).";

ScattGraph::usage =
"ScattGraph[tree, m] returns {vertices, adjacencyMatrix} for the scattering tree.";

FOmbToOm::usage =
"FOmbToOm[ombList, y] converts a list of rational indices (primitive data) into an integer/refined index via the Möbius inversion step.";

ScattIndexImproved::usage =
"ScattIndexImproved[treeList, y, opt] computes the refined index for each tree in treeList, taking care of non-primitive internal states.";

LVTreesFromListRays::usage =
"LVTreesFromListRays[listRays, gam, m] extracts scattering trees of total charge gam from listRays (as produced by ConstructLVDiagram).";

KroneckerDims::usage =
"KroneckerDims[m, Nn] lists populated coprime dimension vectors {n1,n2} for the m-Kronecker quiver with 0 <= n1,n2 <= Nn.";

IntersectRaysNoTest::usage =
"IntersectRaysNoTest[gam, gamp, z, zp, mu] returns the intersection point {x,y} of two rays if it lies strictly in the future of z and z', \
or {} otherwise, without testing vanishing of the DSZ pairing.";

CostPhi::usage =
"CostPhi[gam, s, mu] gives the cost function \\[Phi]_s(\\[Gamma]) = dB + 2 dF - r (mu + 8 s) for gam = {r,dF,dB,ch2}.";

ConstructLVDiagram::usage =
"ConstructLVDiagram[phimax, Nm, m, listRays] constructs the large-volume scattering diagram for local F1. \
Rays are iteratively added by intersecting existing rays and applying Kronecker quiver products with n1+n2 <= Nm and cost <= phimax. \
Assumes m is real. Output is a list of {charge, {x,y}, parent1, parent2, n1, n2, ...}.";

ConstructLVDiagramOpt::usage =
"ConstructLVDiagramOpt[phimax, Nm, m, listRays] is an optimized version of ConstructLVDiagram with the same output conventions.";


(* Period expansions in Pi-Stability slice *)

F1Series::usage =
"F1Series[la, u, Nn, Nr] computes large-volume series expansions of the F1 periods and their first and second derivatives, \
with Richardson acceleration of order Nr.";

PicardFuchsP::usage =
"PicardFuchsP[la, u] gives the coefficient P(la,u) of f'' in the normalized third-order Picard--Fuchs equation f''' + P f'' + Q f' == 0.";

PicardFuchsQ::usage =
"PicardFuchsQ[la, u] gives the coefficient Q(la,u) of f' in the normalized third-order Picard--Fuchs equation f''' + P f'' + Q f' == 0.";

PicardFuchs::usage =
"PicardFuchs[f, la, u] evaluates the normalized third-order Picard--Fuchs operator on f as f''' + P f'' + Q f'.";

PicardFuchsA::usage =
"PicardFuchsA[la, u] returns the 3x3 companion matrix A(la,u) for the first-order system Psi' = A Psi equivalent to the Picard--Fuchs equation.";

SystemMatrix::usage =
"SystemMatrix[dat, la, u] constructs the system matrix used to convert between period bases, given dat = {{t,td},{dt,dtd},{ddt,ddtd}}.";

ComputeTransition::usage =
"ComputeTransition[la, upath, Nn, Nr] numerically integrates the Picard--Fuchs connection along the path upath[t] (t in [0,1]) \
and returns the transition/monodromy matrix in the chosen basis.";

TTDFromSystemMatrix::usage = "TTDFromSystemMatrix[SystemMatrix] TODO";

TauFromSystemMatrix::usage = "TauFromSystemMatrix[SystemMatrix] TODO";

G10FundamentalDomain::usage =
"G10FundamentalDomain is a Graphics object plotting the fundamental domain of the noncongruence subgroup Gamma(10;0,3,0,1;8).";

G10RepeatDomain::usage =
"G10RepeatDomain[{kmin, kmax}] repeats G10FundamentalDomain translated by x -> x + 8 k for kmin <= k <= kmax.";

MirrorCurve2::usage =
"MirrorCurve2[x, y, el, ur] gives the affine mirror-curve equation in variables x, y and parameters el, ur \
(i.e. an expression whose vanishing defines the curve). \
It corresponds to the alternative (el,ur) parameterization used elsewhere in the package.";

MirrorCurveg22::usage =
"MirrorCurveg22[el, ur] gives the Weierstrass invariant g2 (weight-4 Eisenstein invariant) of the mirror curve \
in the (el,ur) parameterization.";

MirrorCurveg32::usage =
"MirrorCurveg32[el, ur] gives the Weierstrass invariant g3 (weight-6 Eisenstein invariant) of the mirror curve \
in the (el,ur) parameterization.";

MirrorCurveJ2::usage =
"MirrorCurveJ2[el, ur] gives the j-invariant j(el,ur) for the alternative parameterization used in the degree-12 map.";

MirrorCurveDelta2::usage =
"MirrorCurveDelta2[el, ur] gives the order-4 factor in the discriminant of the local F1 mirror curve in variables (el,ur). \
Its roots locate the discriminant locus in the (el,ur)-slice.";

PicardFuchsP2::usage =
"PicardFuchsP2[el, ur] gives the coefficient P(el,ur) of f'' in the normalized third-order Picard--Fuchs equation in variables (el,ur).";

PicardFuchsQ2::usage =
"PicardFuchsQ2[el, ur] gives the coefficient Q(el,ur) of f' in the normalized third-order Picard--Fuchs equation in variables (el,ur).";

PicardFuchs2::usage =
"PicardFuchs2[f, el, ur] evaluates the normalized third-order Picard--Fuchs operator in variables (el,ur) on f.";

PicardFuchsA2::usage =
"PicardFuchsA2[el, ur] returns the 3x3 companion matrix A(el,ur) for the first-order system equivalent to the Picard--Fuchs equation.";

PicardFuchsD1::usage =
"PicardFuchsD1[f, el, ur] applies the first original Picard--Fuchs differential operator for local F1 to the function f(el,ur). \
This is a PDE operator in (el,ur), not the reduced third-order ODE in ur.";

PicardFuchsD2::usage =
"PicardFuchsD2[f, el, ur] applies the second original Picard--Fuchs differential operator for local F1 to the function f(el,ur). \
This is a PDE operator in (el,ur), not the reduced third-order ODE in ur.";

PicardFuchsD1BDelta::usage =
"PicardFuchsD1BDelta[g, b, du] applies the second-order Picard--Fuchs differential operator D_1 (in b and du) to g, returning the resulting expression.";

PicardFuchsD2BDelta::usage =
"PicardFuchsD2BDelta[g, b, du] applies the second-order Picard--Fuchs differential operator D_2 (in b and du) to g, returning the resulting expression.";

F1Series2::usage =
"F1Series2[el, ur, Nn, Nr] computes large-volume series expansions in variables (el,ur), including first and second derivatives, \
with Richardson acceleration of order Nr.";

SystemMatrix2::usage =
"SystemMatrix2[dat, el, ur] constructs the system matrix used to convert between period bases for the (el,ur) parameterization.";

F1PeriodsLV::usage = "F1PeriodsLV[el, ur] or F1PeriodsLV[el, ur, Nn, Nr] TODO";

TreeCharge::usage = "TreeCharge[tree] computes the total charge of a tree.";

PrecomputedLVTrees::usage = "PrecomputedLVTrees[gam, m=-1/10] gives a few precomputed trees for m = -1/10.";

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

BW[x_] := Im[PolyLog[2, x]] + Arg[1 - x] Log[Abs[x]];


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

(* old routines *)

(* KerrDoranP[b_] := (5 BW[b] - 4 BW[b^2] + BW[b^3])/(2 \[Pi]); *)

(* KerrDoranQ[b_] := 1/(6 \[Pi]) (3 (Log[1/b] + Log[b]) Log[1/(2 + 1/b + b)^(
     2/3)] + (Log[1/b^2] - Log[1/b] + Log[b]) Log[2 + 1/b + b] - 
   3 (Li2[1/b^3] - Li2[1/b^2] - Li2[b] + Li2[b^2] - 
      2 (Li2[1/b^2] + Log[1 - 1/b^2] Log[1/b^2]) + 
      Log[1 - b] Log[1/b] + (Log[1/b^2] - Log[1/b]) Log[(-1 + b)/b] + 
      2 (Li2[1/b] + Log[1/b] Log[(-1 + b)/b]) + 
      Log[1 - 1/b^3] (Log[1/b^2] - Log[b]) - 
      Log[1 - 1/b^2] (Log[1/b] - Log[b]) + Log[(-1 + b)/b] Log[b] - 
      2 (Li2[b] + Log[1 - b] Log[b]) - (Log[1/b] - Log[b]) Log[
        1 - b^2])); *)

(* Boris *)

ToFund[tau_] :=
  Which[Re[tau] < -1/2 || Re[tau] >= 1/2,
   ToFund[tau - Floor[Re[tau] + 1/2]], Abs[tau] < 1, -1/tau, True,
   tau];

InverseJ[j_, Pre__:MachinePrecision] :=
  Module[{so, tau0, a}, so = Solve[4 a (1 - a) == 1728/j, a];
   tau0 =
    N[I Hypergeometric2F1[1/6, 5/6, 1, 1 - a /. Last[so]]/
       Hypergeometric2F1[1/6, 5/6, 1, a /. Last[so]], Pre];
   ToFund[tau0]];

Mon3ToMon4[m0_, M_] :=
  If[Im[m0] != 0,
   Round[
     1000 {{1, 0, 0, 0}, {0, 1, 0,
        0}, {M[[2, 1]] - Im[M[[2, 1]]]/Im[m0] m0,
        Im[M[[2, 1]]/Im[m0]], M[[2, 2]],
        M[[2, 3]]}, {M[[3, 1]] - Im[M[[3, 1]]]/Im[m0] m0,
        Im[M[[3, 1]]]/Im[m0], M[[3, 2]], M[[3, 3]]}}]/1000,
   Round[
     1000 {{1, 0, 0, 0}, {0, 1, 0, 0}, {M[[2, 1]], 999, M[[2, 2]],
        M[[2, 3]]}, {M[[3, 1]], 999, M[[3, 2]], M[[3, 3]]}}]/1000];

tauB[\[ScriptL]_, Pre__:MachinePrecision] :=
  InverseJ[(27 + 256 \[ScriptL]) (243 + 256 \[ScriptL])^3/
      2^24/\[ScriptL]^3, Pre];

Mon4ToMon3[m0_, M_] := {{1, 0, 0}, {M[[3, 1]] + m0 M[[3, 2]], M[[3, 3]],
    M[[3, 4]]}, {M[[4, 1]] + m0 M[[4, 2]], M[[4, 3]], M[[4, 4]]}}; 

FindPeriodCombination[m0_, P0_, x_, m_, Q_] :=
 Module[{V}, V = FindIntegerNullVector[{1, m0, P0, x}];
  If[V[[4]] == 0, {}, -(V[[1]] + m  V[[2]] + Q V[[3]])/V[[4]]]];


(* For Re[b]>0 *)
KerrDoranQ1[b_] :=
  1/(2 Pi) (1/(2 Pi) (PolyLog[2, b^3] - 4 PolyLog[2, b^2] + 5 PolyLog[2, b]) +
   1/(4 Pi) Log[
     b]*(6 Log[b^2 + b + 1] + Log[b] - 16  Log[b + 1] - 0 4 I Pi) - \[Pi]/6);

(* For Re[b]<0, Im[b]>0
KerrDoranQ2[b_] :=
  1/(12 \[Pi]) (-4 I \[Pi] (Log[b] - 8 Log[1 + b] +
        3 Log[1 + b + b^2]) +
     3 Log[b] (Log[b] - 16 Log[1 + b] + 6 Log[1 + b + b^2]) +
     6 (5 PolyLog[2, b] - 4 PolyLog[2, b^2] + PolyLog[2, b^3]));

(* For Re[b]<0, Im[b]<0 *)
KerrDoranQ3[b_] :=
 1/(12 \[Pi]) (4 I \[Pi] (Log[b] - 8 Log[1 + b] +
       3 Log[1 + b + b^2]) +
    3 Log[b] (Log[b] - 16 Log[1 + b] + 6 Log[1 + b + b^2]) +
    6 (5 PolyLog[2, b] - 4 PolyLog[2, b^2] + PolyLog[2, b^3])); *)

KerrDoranQ[b_] := 
  1/(2 \[Pi])(1/(2 \[Pi]) (-Log[1 - b] Log[
       1/b] - (Log[1/b^2] - Log[1/b]) Log[(-1 + b)/b] - 
     Log[1 - 1/b^3] (Log[1/b^2] - Log[b]) + 
     Log[1 - 1/b^2] (Log[1/b] - Log[b]) - 
     Log[(-1 + b)/b] Log[b] + (Log[1/b] - Log[b]) Log[
       1 - b^2] + (Log[1/b] + 
        Log[b]) Log[-((b (1 + b + b^2))/(1 + b)^4)] + (Log[1/b^2] - 
        Log[1/b] + Log[b]) Log[-((1 + b + b^2)/(1 + b)^2)] - 
     PolyLog[2, 1/b^3] + PolyLog[2, 1/b^2] + 
     2 (Log[1 - 1/b^2] Log[1/b^2] + PolyLog[2, 1/b^2]) - 
     2 (Log[1/b] Log[(-1 + b)/b] + PolyLog[2, 1/b]) + PolyLog[2, b] + 
     2 (Log[1 - b] Log[b] + PolyLog[2, b]) - PolyLog[2, b^2]) - \[Pi]/6);

KerrDoranFix0[b_] := -Pi HeavisideTheta[-Re[b]] (Sign[1 - Abs[b]] + 
   2 HeavisideTheta[Abs[b] - 1] HeavisideTheta[
     Im[b] - Sqrt[3] Re[b]] HeavisideTheta[-Im[b] - Sqrt[3] Re[b]]);

KerrDoranFix1[b_] := 3 I Sign[
  Im[b]] HeavisideTheta[-(1/2) - Re[b]] (1 - 
   HeavisideTheta[-Im[b] - Sqrt[3] Re[b]] HeavisideTheta[
     Im[b] - Sqrt[3] Re[b]]);

KerrDoranQ1Fix[b_] := KerrDoranQ1[b] + KerrDoranFix0[b] + KerrDoranFix1[b] Log[b];

KerrDoranP[b_] := 1/(2 \[Pi])((5 BW[b] - 4 BW[b^2] + BW[b^3])/(2 \[Pi]) + ((Arg[b] Log[Abs[(b (1 + b + b^2)^3)/(1 + b)^8]])/(2 \[Pi])));

KerrDoranFormula[avec_, bvec_, dvec_, evec_, A_, 
   B_] := -Sum[
     If[avec[[j]] === bvec[[k]], 0, 
      dvec[[j]]*evec[[k]] (PolyLog[2, 
          avec[[j]]/bvec[[k]]] + (Log[avec[[j]]] - 
            Log[bvec[[k]]]) Log[1 - avec[[j]]/bvec[[k]]])], {j, 
      Length[dvec]}, {k, Length[evec]}] - 
   Log[B] Sum[dvec[[j]] Log[avec[[j]]], {j, Length[dvec]}] + 
   Log[A] Sum[evec[[k]] Log[bvec[[k]]], {k, Length[evec]}];

KerrDoranFormulaExact[avec_, bvec_, dvec_, evec_, A_, 
   B_] := -Sum[
     If[avec[[j]] == bvec[[k]], 0, 
      dvec[[j]]*evec[[k]] (PolyLog[2, 
          avec[[j]]/bvec[[k]]] + (Log[avec[[j]]] - 
            Log[bvec[[k]]]) Log[1 - avec[[j]]/bvec[[k]]])], {j, 
      Length[dvec]}, {k, Length[evec]}] - 
   Log[B] Sum[dvec[[j]] Log[avec[[j]]], {j, Length[dvec]}] + 
   Log[A] Sum[evec[[k]] Log[bvec[[k]]], {k, Length[evec]}];

KerrDoranElFromb[b_] := -((b (1 + b + b^2)^3)/(1 + b)^8);

KerrDoranUrFromb[b_] := -((1 + b)^4/((1 + b + b^2) (1 + 4 b + b^2)));

KerrDoranElUrFromb[b_] := {KerrDoranElFromb[b], KerrDoranUrFromb[b]};

KerrDoranbFromEl[el_] := b/.Solve[{KerrDoranElFromb[b]==el, (Abs[b]>1)||(Abs[b]==1&&Im[b]>=0)}, b];

(*Options*)
Options[AnalyticContinuationOfRoots] = {"Steps" -> 100, 
   "Path" -> "Straight"};
(**)(*Method:\
Track roots along a path in parameter space by continuous deformation*)
AnalyticContinuationOfRoots[poly_, x_, t_, t0_, t1_, 
   OptionsPattern[]] := 
  Module[{roots0, nSteps, tValues, rootsPath, currentRoots, nextRoots,
     distances, permutation, minPos, usedIndices},(*Step 1:
   Find initial roots at t=t0*)
   roots0 = x /. NSolve[poly /. t -> t0, x];
   (*Step 2:Create path from t0 to t1*)
   nSteps = OptionValue["Steps"];
   If[OptionValue["Path"] === "Straight",
    tValues = Table[t0 + (t1 - t0) s, {s, 0, 1, 1/(nSteps - 1)}],
    tValues = OptionValue["Path"]];
   (*Step 3:Track roots along the path*)
   rootsPath = {roots0};
   currentRoots = roots0;
   Do[(*Solve polynomial at next t value*)
    nextRoots = x /. NSolve[poly /. t -> tValues[[i]], x];
    (*Match roots by proximity using Hungarian-style greedy matching*)
    If[Length[currentRoots] == Length[nextRoots],
     (*Compute distance matrix*)
     distances = 
      Table[Abs[currentRoots[[j]] - nextRoots[[k]]], {j, 
        Length[currentRoots]}, {k, Length[nextRoots]}];
     (*Greedy assignment:for each current root,find closest next root*)
     permutation = Table[0, {Length[currentRoots]}];
     usedIndices = {};
     Do[(*Find minimum in row j excluding already used columns*)
      minPos = 0;
      minVal = Infinity;
      Do[
       If[! MemberQ[usedIndices, k] && distances[[j, k]] < minVal, 
        minVal = distances[[j, k]];
        minPos = k;], {k, Length[nextRoots]}];
      permutation[[j]] = minPos;
      AppendTo[usedIndices, minPos];, {j, Length[currentRoots]}];
     (*Reorder nextRoots according to permutation*)
     nextRoots = nextRoots[[permutation]];];
    currentRoots = nextRoots;
    AppendTo[rootsPath, currentRoots];, {i, 2, Length[tValues]}
    ];
   (*Return association with path and roots*)
   <|"tValues" -> tValues, "rootsPaths" -> Transpose[rootsPath], 
    "finalRoots" -> Last[rootsPath], 
    "initialRoots" -> First[rootsPath]|>];


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
   If[! IntegerQ[SecondChernClass[{r,dF,dB,ch2}]], 
    Print["Non integer second Chern class !"]];
   If[r == 0,
    (ch2 + dB m - dF m)/(dB + 2 dF),
    (dB + 2 dF - m r)/8/r - 
     Sign[r]/8 Sqrt[Max[DiscF1[{r, dF, dB, ch2}, m], 0]]
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
   1 -> Directive[Opacity[1, Red], Thickness[.0017]],
   2 -> Directive[Opacity[1, Blue], Thickness[.0015]],
   3 -> Directive[Opacity[1, Green], Thickness[.0012]]
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

ReZLV[{r_, df_, db_, ch2_}, {s_, t_}, m_] := -ch2 - db m + df m + (
   m^2 r)/2 + (db + 2 df - m r) s - 4 r s^2 + 4 r t^2;
QConditions[Coll_, m_, {s_, t_}] := 
  And @@ (ReZLV[#, {s, t}, m] < 0 & /@ Coll);
ClearAll[QDomain];
Options[QDomain] = Options[RegionPlot];
QDomain[Coll_, m_, tooltip_, opt___ : OptionsPattern[]] := 
  RegionPlot[Tooltip[QConditions[Coll, m, {s, t}], tooltip], {s, -1, 2}, {t, 0, 1}, opt];

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

TTDFromSystemMatrix[SystMat_] := SystMat[[1, {2, 3}]];

TauFromSystemMatrix[SystMat_] := SystMat[[2, 3]]/SystMat[[2, 2]];

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

MirrorCurve2[x_, y_, el_, ur_] := -(1/ur) + x + y + (el (1 + y))/(x y);
MirrorCurveg22[el_, ur_] := 108 (1 + 16 el^2 ur^4 - 8 el ur^2 (1 + 3 ur));
MirrorCurveg32[el_, ur_] := -216 (1 - 64 el^3 ur^6 - 12 el ur^2 (1 + 3 ur) + 
   24 el^2 ur^4 (2 + 6 ur + 9 ur^2));
MirrorCurveJ2[el_, ur_] := (1 + 16 el^2 ur^4 - 
  8 el ur^2 (1 + 3 ur))^3/(el^3 ur^8 (1 + ur - 8 el ur^2 - 
   36 el ur^3 + el (-27 + 16 el) ur^4));
MirrorCurveDelta2[el_, ur_] := 1 + ur - 8 el ur^2 - 36 el ur^3 + el (-27 + 16 el) ur^4;
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

PicardFuchsD1[f_, el_, ur_] := (-2 ur^3 D[f, ur] - ur^4 D[f, {ur, 2}] + 3 D[f, el] - 
  ur D[f, el, ur] + 3 el D[f, {el, 2}]);

PicardFuchsD2[f_, el_, ur_] := 
 ur (1 + ur) D[f, ur] + ur^2 (1 + ur) D[f, {ur, 2}] + 
  el (4 D[f, el] - ur (4 + 3 ur) D[f, ur, el] + 4 el D[f, {el, 2}]);

PicardFuchsD1BDelta[g_, b_, du_] := (
  (
    (-1 + b)^5 (
        -5 + 6 du - 6 du^2 + 2 du^3
        + b^16 (-5 + 6 du - 6 du^2 + 2 du^3)
        + 6 b^2  (-66 + 61 du - 69 du^2 + 26 du^3)
        + 6 b^14 (-66 + 61 du - 69 du^2 + 26 du^3)
        + b    (-64 + 72 du - 78 du^2 + 28 du^3)
        + b^15 (-64 + 72 du - 78 du^2 + 28 du^3)
        - 2 b^6  (8954 + 87 du - 111 du^2 + 50 du^3)
        - 2 b^10 (8954 + 87 du - 111 du^2 + 50 du^3)
        + 12 b^5  (-852 + 135 du - 150 du^2 + 56 du^3)
        + 12 b^11 (-852 + 135 du - 150 du^2 + 56 du^3)
        + 4 b^3  (-396 + 261 du - 297 du^2 + 113 du^3)
        + 4 b^13 (-396 + 261 du - 297 du^2 + 113 du^3)
        - 6 b^8 (4653 + 662 du - 738 du^2 + 276 du^3)
        + 2 b^4  (-2298 + 894 du - 1008 du^2 + 377 du^3)
        + 2 b^12 (-2298 + 894 du - 1008 du^2 + 377 du^3)
        - 2 b^7 (12496 + 1368 du - 1533 du^2 + 568 du^3)
        - 2 b^9 (12496 + 1368 du - 1533 du^2 + 568 du^3)
      ) D[g, du]
    +
    (-1 + b)^7 du (
        -3 + 6 du - 4 du^2 + du^3
        + b^14 (-3 + 6 du - 4 du^2 + du^3)
        + 4 b   (-9 + 21 du - 15 du^2 + 4 du^3)
        + 4 b^13 (-9 + 21 du - 15 du^2 + 4 du^3)
        + 4 b^3  (-174 + 504 du - 379 du^2 + 107 du^3)
        + 4 b^11 (-174 + 504 du - 379 du^2 + 107 du^3)
        + b^2  (-201 + 528 du - 392 du^2 + 109 du^3)
        + b^12 (-201 + 528 du - 392 du^2 + 109 du^3)
        + 24 b^7 (-198 + 704 du - 527 du^2 + 148 du^3)
        + 4 b^5 (-759 + 2547 du - 1913 du^2 + 539 du^3)
        + 4 b^9 (-759 + 2547 du - 1913 du^2 + 539 du^3)
        + 3 b^6 (-1419 + 4970 du - 3724 du^2 + 1046 du^3)
        + 3 b^8 (-1419 + 4970 du - 3724 du^2 + 1046 du^3)
        + b^4  (-1683 + 5292 du - 3984 du^2 + 1124 du^3)
        + b^10 (-1683 + 5292 du - 3984 du^2 + 1124 du^3)
      ) D[g, {du, 2}]
    +
    (1 + b)^9 (1 + 5 b + 6 b^2 + 5 b^3 + b^4) (
        3 (1 + 4 b + b^2)^2 (-1 - 7 b - 15 b^2 - 18 b^3 - 11 b^4 - 3 b^5 + b^6)
          D[g, b]
        +
        (-1 + b) (
            (-1 + b)^3 (
                -1 + b^6 (-1 + du) + du
                + b   (-2 + 9 du)
                + b^5 (-2 + 9 du)
                + b^2 (1 + 27 du)
                + b^4 (1 + 27 du)
                + b^3 (4 + 34 du)
              ) D[g, b, du]
            +
            3 b (1 + b) (1 + 5 b + 6 b^2 + 5 b^3 + b^4)^2
              D[g, {b, 2}]
          )
      )
  ) / (3 (-1 + b) b (1 + b)^10 (1 + b + b^2)^3 (1 + 4 b + b^2)^3)
);

PicardFuchsD2BDelta[g_, b_, du_] :=
 Module[
  {
   den,
   a01, a02, a10, a11, a20,
   pref
  },

  (* Common denominator *)
  den = 4 (-1 + b) b^2 (1 + b)^2 (1 + b + b^2)^4 (1 + 4 b + b^2)^3;

  (* Coefficients multiplying derivatives *)
  a01 = (-1 + b)^5 (
      (-1 + du) du
      + b^12 (-1 + du) du
      + b^5  (-454 + 82 du - 90 du^2)
      + b^7  (-454 + 82 du - 90 du^2)
      + b    (-5 - 10 du + 12 du^2)
      + b^11 (-5 - 10 du + 12 du^2)
      + b^4  (-272 - 35 du + 27 du^2)
      + b^8  (-272 - 35 du + 27 du^2)
      + b^2  (-34 - 37 du + 51 du^2)
      + b^10 (-34 - 37 du + 51 du^2)
      + b^3  (-117 - 64 du + 86 du^2)
      + b^9  (-117 - 64 du + 86 du^2)
      - 2 b^6 (270 - 65 du + 87 du^2)
    );

  a02 = (-1 + b)^7 du (
      (-1 + du)^2
      + b^10 (-1 + du)^2
      + 2 b   (4 - 11 du + 7 du^2)
      + 2 b^9 (4 - 11 du + 7 du^2)
      + 3 b^2 (11 - 34 du + 26 du^2)
      + 3 b^8 (11 - 34 du + 26 du^2)
      + 12 b^5 (16 - 49 du + 41 du^2)
      + 4 b^3 (22 - 69 du + 57 du^2)
      + 4 b^7 (22 - 69 du + 57 du^2)
      + b^4 (158 - 492 du + 405 du^2)
      + b^6 (158 - 492 du + 405 du^2)
    );

  pref = b (1 + 6 b + 11 b^2 + 11 b^3 + 6 b^4 + b^5);

  a10 = 4 (-1 - 6 b - 8 b^2 - 4 b^3 + b^4) (1 + 5 b + 6 b^2 + 5 b^3 + b^4)^2;

  a11 = (-1 + b)^4 (
      -1 - 2 du + 3 du^2
      + b^8 (-1 - 2 du + 3 du^2)
      + 2 b   (-2 - 7 du + 15 du^2)
      + 2 b^7 (-2 - 7 du + 15 du^2)
      + b^2 (-4 - 44 du + 111 du^2)
      + b^6 (-4 - 44 du + 111 du^2)
      + 2 b^4 (5 - 58 du + 132 du^2)
      + b^3 (4 - 98 du + 210 du^2)
      + b^5 (4 - 98 du + 210 du^2)
    );

  a20 = 4 b (-1 + b^2) (1 + b + b^2)^3 (1 + 4 b + b^2)^2;

  (* Assemble operator *)
  (
    a01 D[g, du]
    + a02 D[g, {du, 2}]
    + pref (
        a10 D[g, b]
        + a11 D[g, b, du]
        + a20 D[g, {b, 2}]
      )
  )/den
 ];


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

F1PeriodsLV[el_, ur_, Nn_ : 512, Nr_ : 10] := 
  SystemMatrix2[F1Series2[el, ur, Nn, Nr], el, ur];

TreeCharge[arg : {r_, d1_, d2_, ch2_} /; FreeQ[arg, Ch]] := {r, d1, d2, ch2};
TreeCharge[trees_List] := Total[TreeCharge /@ trees];
TreeCharge[arg : Except[_List]] := arg /. repCh;

(* precomputed trees for m = -1/10 *)
PrecomputedLVTrees[{0, 0, 1, 1/2}, -1/10] = {{{-1, -1, 0, 0}, {1, 1, \
1, 1/2}}};
PrecomputedLVTrees[{0, 1, 0, 0}, -1/10] = {{{-1, 0, 0, 0}, {1, 1, 0, \
0}}};
PrecomputedLVTrees[{0, 1, 1, 1/2}, -1/10] = {{{-1, 0, 0, 0}, {1, 1, 1, \
1/2}}};
PrecomputedLVTrees[{0, 2, 1, 1/2}, -1/10] = {{{1, 1, 1, 1/2}, {{-2, 0, \
0, 0}, {1, 1, 0, 0}}}};
PrecomputedLVTrees[{0, 2, 2, 0}, -1/10] = {{{1, 1, 1, 1/2}, {-1, 1, 1, \
-1/2}}};
PrecomputedLVTrees[{0, 3, 2, 0}, -1/10] = {{{-1, 1, 1, -1/2}, {{1, 1, \
1, 1/2}, {{-1, 0, 0, 0}, {1, 1, 0, 0}}}}};
PrecomputedLVTrees[{0, 3, 3, 1/2}, -1/10] = {{{-1, 1, 1, -1/2}, {{-1, \
0, 0, 0}, {2, 2, 2, 1}}}};
PrecomputedLVTrees[{0, 4, 2, 0}, -1/10] = {{{-1, 2, 1, -3/2}, {1, 2, \
1, 3/2}}, {{-1, 1, 1, -1/2}, {{1, 1, 1, 1/2}, {{-2, 0, 0, 0}, {2, 2, \
0, 0}}}}};
PrecomputedLVTrees[{1, 2, 2, 2}, -1/10] = {{{1, 2, 1, 3/2}, {{-1, -1, \
0, 0}, {1, 1, 1, 1/2}}}};
PrecomputedLVTrees[{1, 0, 1, -1/2}, -1/10] = {{{1, 0, 0, 0}, {{-1, 2, \
2, -2}, {1, -2, -1, 3/2}}}};
PrecomputedLVTrees[{1, 1, 2, 0}, -1/10] = {{{1, 1, 1, 1/2}, {{-1, 2, \
2, -2}, {1, -2, -1, 3/2}}}};
PrecomputedLVTrees[{1, 2, 3, 3/2}, -1/10] = {{{{-1, 2, 2, -2}, {1, -2, \
-1, 3/2}}, {{1, 2, 1, 3/2}, {{-1, -1, 0, 0}, {1, 1, 1, 1/2}}}}};
PrecomputedLVTrees[{-1, 0, 1, 1/2}, -1/10] = {{{-1, 0, 0, 0}, {{-1, \
-1, 0, 0}, {1, 1, 1, 1/2}}}};
PrecomputedLVTrees[{-1, -2, -1, -3/2}, -1/10] = {{-1, -2, -1, -3/2}};
PrecomputedLVTrees[{1, 0, 0, -1}, -1/10] = {{{-1, 2, 2, -2}, {2, -2, \
-2, 1}}, {{{-1, 2, 2, -2}, {1, -2, -1, 3/2}}, {{-1, 2, 1, -3/2}, {2, \
-2, -2, 1}}}};
PrecomputedLVTrees[{1, 0, 0, -2}, -1/10] = {{{1, -1, -1, 1/2}, {{-1, \
3, 2, -4}, {1, -2, -1, 3/2}}}, {{{-1, 5, 4, -12}, {1, -5, -3, 21/2}}, \
{{-1, 2, 1, -3/2}, {2, -2, -2, 1}}}, {{{-1, 3, 2, -4}, {1, -2, -2, 2}}, \
{{1, -1, -1, 1/2}, {{-1, 2, 2, -2}, {1, -2, -1, 3/2}}}}};
PrecomputedLVTrees[{2, 1, 1, -1/2}, -1/10] = {{{3, 0, 0, 0}, {-1, 1, \
1, -1/2}}};
PrecomputedLVTrees[{2, 1, 0, -1}, -1/10] = {{{2, 0, 0, 0}, {{-1, 2, 1, \
-3/2}, {1, -1, -1, 1/2}}}};


End[]; (* `Private` *)

EndPackage[];
