(* ::Package:: *)

(* ::Title:: *)
(*Scattering diagram for local F0,F1 - routines*)


Print["F0Scattering 1.18, June 22, 2025 - A package for evaluating DT invariants on F0,F1"];


BeginPackage["IndexVars`"];
Protect[tau, y];
EndPackage[];

BeginPackage["F0Scattering`"]
Needs["IndexVars`"]
Needs["CoulombHiggs`"]


(* symbols *)
Ch::usage = "Ch[m1,m2] denotes the charge vector for O(m1,m2) on F0";
GV::usage = "GV[m1,m2,ch2] denotes the charge vector for [0,m1,m2,ch2]";
Kr::usage = "\!\(\*SubscriptBox[\(Kr\), \(m_\)]\)[p_,q_] denotes the index of the Kronecker quiver with m arrows, dimension vector (p,q)";
\[Gamma]::usage = "\!\(\*SubscriptBox[\(\[Gamma]\), \(i\)]\) denotes the charge vector for the i-th node of the McKay quiver";
a::usage = "Parameter running from 0 to 1";
y::usage = "Refinement parameter";
tau::usage = "Kahler modulus"; 
tau1::usage = "Real part of tau";
tau2::usage = "Imaginary part of tau";
gam1::usage = "Charge vector [r,d1,d2,ch2] for E1[-1] for orbifold quiver";
gam2::usage = "Charge vector [r,d1,d2,ch2] for E2[-1] for orbifold quiver";
gam3::usage = "Charge vector [r,d1,d2,ch2] for E3[-1] for orbifold quiver";
gam4::usage = "Charge vector [r,d1,d2,ch2] for E4[-1] for orbifold quiver";
gam1p::usage = "Charge vector [r,d1,d2,ch2] for E1[-1] for phase I quiver";
gam2p::usage = "Charge vector [r,d1,d2,ch2] for E2[-1] for phase I quiver";
gam3p::usage = "Charge vector [r,d1,d2,ch2] for E3[-1] for phase I quiver";
gam4p::usage = "Charge vector [r,d1,d2,ch2] for E4[-1] for phase I quiver";
gamd::usage= "gamd[i] is the Chern vector attached to i-th node of quiver";
Ch1::usage = "Ch1[mH,mC] denotes the charge vector for O(mH,mC) on F1";
GV1::usage = "GV[mH,mC,ch2] denotes the charge vector for [0,mH,mC,ch2]";

(* global variables *)
LVTreesF0::usage = "LVTreesF0[{r_,d1_,d2_,ch2_] gives the list of precomputed trees at m=1/2, when available";
LVTreesF1::usage = "LVTreesF1[{r_,dH_,dC_,ch2_] gives the list of precomputed trees, when available";
McKayTrees::usage = "McKayTrees[{n1_,n2_,n3_,n4_}] gives the list of precomputed McKay trees, when available";
ListSubsetsAny::usage = "Precomputed list of all binary splits of Range[n] for n=2..10, used by ListStableTrees";
FourierC::usage = "FourierC[m_] gives the list of precomputed Fourier coefficients of C[m]";
FourierCp::usage = "FourierCp[m_] gives the list of precomputed Fourier coefficients of C'[m]";
FourierCpp::usage = "FourierCpp[m_] gives the list of precomputed Fourier coefficients of C''[m]";
Initialtau1::usage = "Initial value of tau1 for FindRoot";
Initialtau2::usage = "Initial value of tau2 for FindRoot";
MLV::usage = "Monodromy matrix around cusp at infinity";
MLR1::usage = "Monodromy matrix under m->m+1,T->T-1";
MLR2::usage = "Monodromy matrix under m->m-1,T->T";
MCon::usage = "Monodromy matrix around cusp at 0";
MOrb::usage = "Monodromy matrix under tau->(tau-1)/(4tau-3)";
MFB::usage = "Fiber-base duality matrix";
MB::usage ="MB[k_] gives the Z2-monodromy around tauB";
MBt::usage ="MBt[k_] gives the tilde Z2-monodromy around tauB";
MZ4::usage ="MZ4[k_] gives the Z4-monodromy around tauB";
MZ4t::usage ="MZ4t[k_] gives the tilde Z4-monodromy around tauB";
Sigma::usage = "Matrix such that Z=gamma.Sigma.Transpose[{1,m,T,TD}]";
Mat::usage = "DSZ matrix for McKay scattering diagram";
InitialRaysOrigin::usage ="InitialRaysOrigin[m_] is a List of initial points (u_i,v_i) for initial rays in McKay scattering diagram";

(* Manipulating Chern vectors *)
Euler::usage = "Euler[{r_,d1_,d2_,ch2_},{rr_,dd1_,dd2_,cch2_}] computes the Euler form on D(F_0)";
DSZ::usage = "DSZ[{r_,d1_,d2_,ch2_},{rr_,dd1_,dd2_,cch2_}] computes the antisymmetrized Euler form 2rk (cc11+cc12)-2 rrk (c11+c12) on D(F_0)";
Disc::usage = "Disc[{r_,d1_,d2_,ch2_}] computes the discriminant -ch2/r+d1 d2/r^2";
DiscF0::usage = "DiscF0[{r_,d1_,d2_,ch2_},m_:0] computes the discriminant ((d1-d2-m r)^2+4 d1 d2-4 r ch2)/(8 r^2)";
DimGieseker::usage = "DimGieseker[{r_,d1_,d2_,ch2_}]Computes expected dimension of moduli space of Gieseker-semi stable sheaves";
repChO::usage = "Replaces Ch[m1_,m2_] by string O(m1,m2)";
repCh::usage = "Replacement Ch[m1,m2] by the charge vector {1,m1,m2,m1 m2}";
ChernToCh::usage="ChernToCh[f_] replaces {r,d1,d2,ch2} by r Ch[d1/r,d2/r] whenever the discriminant vanishes";
GenSlope::usage = "GenSlope[{r_,d1_,d2_,ch2_},m_] computes the slope (d1+d2)/(2r) if r<>0, or (ch2-m d2)/d1+d2) if r=0";
MutateCollection::usage="MutateCollection[Coll_,klist_] acts on the list of Chern vectors Coll by the successive mutations in klist, which is a list of {node number,sign}, with sign=1 for right mutation, -1 for left mutation";
ExtFromStrong::usage="ExtFromStrong[Coll_] computes the Chern vectors of the objects in the Ext collection dual to a strong collection"; 
StrongFromExt::usage="ExtFromStrong[Coll_] computes the Chern vectors of the objects in the strong collection dual to an Ext collection"; 
SameHalfPlaneQ::usage="SameHalfPlaneQ[Zlist_] gives True if all elements of Zlist are in a common half plane";
DimFromCh::usage="DimFromCh[Coll_,gam_] computes the dimension vector from Chern vector gam with respect to Ext-collection Coll";
ChFromDim::usage="ChFromDim[Coll_,Nvec_] computes the Chern vector from dimension vector Nvec with respect to Ext-collection Coll";


repKr::usage = "replaces \!\(\*SubscriptBox[\(Kr\), \(m_\)]\)[p_,q_]\[RuleDelayed]1";
SpecFlow::usage = "SpecFlow[{r_,d1_,d2_,ch2_},{m1_,m2_}] gives {r,d1+r m1,d2+r m2,d2 m1+d1 m2+m1 m2 r}";

(* Large volume central charge, walls of marginal stability *)
ZLV::usage = "ZLV[{r_,d1_,d2_,ch2_},{s_,t_},m_] computes the standard central charge Z=-r tau(tau+m)+d1 tau+d2(tau+m)-ch2 with tau=s+I t";
Wall::usage = "Wall[{r_,d1_,d2_,ch2_},{rr_,dd1_,dd2_,cch2_},{s_,t_},m_] computes Im[Z[gamma] Conjugate[Z[ggamma]]/t";
QuiverDomain::usage="QuiverDomain[Coll_,psi_,m_,OptionsPattern[]] plots the Region where the LV central charges of the vectors in Coll have Re[e^{-I psi}Z]<0, and the region where they are in same half-plane."; 
StabilityRegion::usage="StabilityRegion[gam1_,gam2_,m_] plots the region where Im[Z[gam1]Conjugate[Z[gam2]] DSZ[gam1,gam2]>0";


(* basic operations on trees and lists of trees *)
TreeCharge::usage = "TreeCharge[tree_] is the total charge {r,d1,d2,ch2} carried by tree (or list of trees).";
TreeHue::usage= "TreeHue[i_,n_] specifies the color for the i-th tree among a list of n";
TreeConstituents::usage = "TreeConstituents[Tree_] gives the flattened list of constituents of Tree";
FlipTree::usage = "FlipTree[Tree_] constructs the reflected tree under Ch[m1,m2]:>-Ch[-m1,-m2]";
ShiftTree::usage = "ShiftTree[{k1_,k2},Tree_] constructs the shifted tree under Ch[m1,m2]:>Ch[m1+k1,m2+k2]";
ScattCheck::usage = "ScattCheck[Tree_,m_]returns {charge,{x,y}} of the root vertex if Tree is consistent, otherwise {total charge,{}}";
ScattSort::usage = "ScattSort[LiTree_,m_] sorts trees in LiTree by growing radius";
ScattGraph::usage = "ScattGraph[Tree_,m_] extracts the list of vertices and adjacency matrix of Tree";
ScattGraphInternal::usage = "ScattGraphInternal[Tree_,m_] is a light version of ScattDiagInternal used to construct the scattering graph";

(* constructing scattering trees *)
ScanConstituents::usage = "ScanConstituents[gam_,{smin_,smax_},{n_,np_},m_,phimax_] searches possible list of constituents \[PlusMinus]O(k,k), \[PlusMinus]O(k+1,k) with slope in [smin,smax], number of D/Dbar constituents<=(n,np), cost function less than phimax and charges adding up to gam";
ListStableTrees::usage = "ListStableTrees[LiCh_,{s0_,t0_},m_] constructs consistent trees from constituents list LiCh={k_i Ch[m_i,n_i]}, which are stable at (s0,t0)";
ScanAllTrees::usage = "ScanAllTrees[{r_,d1_,d2_,ch2_},{s0_,t0_},n_,m_] constructs all possible trees with up to (n,n) left- and right- moving charges adding up to [r,d1,d2,ch2) leading to an outgoing ray through the point (s,t)";
ContractInitialRays::usage="ContractInitialRays[Trees_,m_] replaces primary scatterings in LV+eps diagram by their total charge";
DecontractInitialRays::usage="DecontractInitialRays[Trees_,m_] replaces initial rays in LV tree by initial rays in LV+eps scattering diagram";
KroneckerDims::usage="KroneckerDims[m_,Nn_] gives the list of populated dimension vectors {n1,n2} for Kronecker quiver with m arrows, with (n1,n2) coprime and 0<=n1,n2<=Nn"; 
ConstructLVDiagram::usage="ConstructLVDiagram[smin_,smax_,phimax_,Nm_,m_,ListRays_] constructs the LV scattering diagram with initial rays in the interval [smin,smax], cost function up to phi, scattering products up to (Nm,Nm) at each intersection; m is assumed to be real; The output consists of a list of  {charge, {x,y}, parent1, parent2,n1,n2 }; If ListRays is not empty, then uses it as initial rays."
LVTreesFromListRays::usage="lVTreesFromListRays[ListRays_,{r_,d1_,d2_,ch2_},m_] extract the trees with given charge in the List of rays, constructed by ConstructLVDiagram";
CostPhi::usage="CostPhi[{r_,d1_,d2_,ch2_},s_,m_]:=d1+d2-r(2s+m)";


McKayrep::usage = "replaces {n1_,n2_,n3_,n4_} by n1 \!\(\*SubscriptBox[\(\[Gamma]\), \(1\)]\)+n2  \!\(\*SubscriptBox[\(\[Gamma]\), \(\(2\)\(\\\ \)\)]\)+n3 \!\(\*SubscriptBox[\(\[Gamma]\), \(3\)]\)+n4 \!\(\*SubscriptBox[\(\[Gamma]\), \(4\)]\)";
McKayrepi::usage = "replaces  n1 \!\(\*SubscriptBox[\(\[Gamma]\), \(1\)]\)+n2  \!\(\*SubscriptBox[\(\[Gamma]\), \(\(2\)\(\\\ \)\)]\)+n3 \!\(\*SubscriptBox[\(\[Gamma]\), \(3\)]\)+n4 \!\(\*SubscriptBox[\(\[Gamma]\), \(4\)]\) by {n1_,n2_,n3_,n4_}";
(* routines for phase II scattering diagram *)
(*McKayDim::usage = "McKayDim[{n1_,n2_,n3_,n4_}] computes the dimension of quiver moduli space in anti-attractor chamber";*)
nFromCh::usage = "nFromCh[{rk_,c11_,c12_,ch2_}] produces the corresponding dimension vector {n1,n2,n3,n4} for phase II quiver";
ChFromn::usage = "ChFromn[{n1_,n2_,n3_,n4_}]:= produces the Chern character associated to dimension vector {n1,n2,n3,n4} for phase II quiver";
repChn::usage = "repChn transforms Ch[m1_,m2_] into corresponding dimension vector";
(* routines for phase I scattering diagram *)
ChFromnI::usage = "ChFromnI[{n1_,n2_,n3_,n4_}]:= produces the Chern character associated to dimension vector {n1,n2,n3,n4} for phase I quiver";
nIFromCh::usage = "nIFromCh[{rk_,c11_,c12_,ch2_}] produces the corresponding dimension vector {n1,n2,n3,n4} for phase I quiver";
repChnI::usage = "repChnI transforms Ch[m1_,m2_] into corresponding dimension vector";
(*McKayIDim::usage = "McKayIDim:[{n1_,n2_,n3_,n4_}] computes the dimension of quiver moduli space in anti-attractor chamber";*)

(* common routines for quiver scattering diagram *)
TreeToGamma::usage="TreeToGamma[Tree_] replaces dimension vectors by linear combinations of Subscript[\[Gamma], i]";
McKayDSZ::usage = "McKayDSZ[{n1_,n2_,n3_,n4_},{nn1_,nn2_,nn3_,nn4_}] computes the antisymmetrized Euler form";
McKayRayEq::usage = "McKayRayEq[{n1_,n2_,n3_,n4_},{u_,v_}] gives the linear form vanishing on the scattering ray";
McKayVec::usage = "McKayVec[{n1_,n2_,n3_,n4_}] computes the positive vector along the ray";
McKayRay::usage = "McKayRay[{n1_,n2_,n3_,n4_},{u_,v_},{k1_,k2_},tx_] produces from {u,v}+k1 vec to (u,v)+k2 vec, decorated with text at the target";
McKayScattIndex::usage = "McKayScattIndex[TreeList_] computes the index for each tree in TreeList; do not trust the result if internal lines have non-primitive charges";
McKayScattIndexInternal::usage = "McKayScattIndexInternal[Tree_] computes {total charge, list of Kronecker indices associated to each vertex in Tree}";
McKayScanAllTrees::usage = "McKayScanIAllTrees[Nvec_,m_] generates consistent scattering trees with leaves carrying charge {p,0,0,0}, {0,p,0,0},{0,0,p,0},{0,0,0,p}, with non-zero DSZ pairing at each vertex, with distinct support";
McKayListAllTrees::usage = "McKayIListAllTrees[Nvec_] generates all trees with leaves carrying charge {p,0,0,0}, {0,p,0,0},{0,0,p,0},{0,0,0,p} and with non-zero DSZ pairing at each vertex ";
McKayScattCheck::usage = "McKayScattCheck[Tree_,m_]returns {charge,{u,v}coordinate} of the root vertex if Tree is consistent, otherwise {total charge,{}}";
McKayScattGraph::usage = "McKayScattGraph[Tree_,m_]extracts the list of vertices and adjacency matrix of Tree";
McKayScattDiag::usage= "McKayScattDiag[TreeList_,m_] draws McKay scattering diagram in (u,v) plane for each tree in Treelist"; 
McKayScattDiagInternal::usage= "McKayScattDiagInternal[Tree_,m_] constructs total charge, coordinate of root and list of line segments in (u,v) coordinates; used by McKayScattDiag"; 
McKayIntersectRays::usage = "McKayIntersectRays[{n1_,n2_,n3_,n4_},{nn1_,nn2_,nn3_,nn4_},m_] returns intersection point (u,v) of two rays, or {} if they are collinear;
  McKayIintersectRays[{n1_,n2_,n3_,n4_},{nn1_,nn2_,nn3_,nn4_},z_,zz,m_]returns intersection point (u,v) of two rays if the intersection point lies upward from z and z', or {} otherwise";
McKayIntersectRaysNoTest::usage = "McKayIntersectRays[{n1_,n2_,n3_,n4_},{nn1_,nn2_,nn3_,nn4_},,z_,zz,m_] returns intersection point (u,v) of two rays if the intersection point lies strictly upward from z and z', or {} otherwise, without testing non-vanishing of DSZ product";
McKayInitialRays::usage="McKayInitialRays[L_,m_] draws the initial rays in (u',v') plane, rescaling each arrow by a factor of L";
McKayScattIndexImproved::usage="McKayScattIndexImproved[TreeList_, opt_] computes the index for each tree in TreeList, taking care of non-primitive internal states";
McKayScattIndexImprovedInternal::usage="McKayScattIndexImprovedInternal[Tree_, opt_] computes the index for Tree, taking care of non-primitive internal states";
FOmbToOm::usage="FOmbToOm[OmbList_] computes integer index from list of rational indices, used internally by FScattIndex";
ConstructMcKayDiagram::usage="ConstructMcKayDiagram[phimax_,m_,ListRays_] constructs the quiver scattering diagram with height up to phimax, scattering products up to n1+n2<=Nm at each intersection; m is assumed to be real; The output consists of a list of  {charge, {x,y}, parent1, parent2,n1,n2 }; If ListRays is not empty, then uses it as initial rays.";
McKayTreesFromListRays::usage="McKayTreesFromListRays[ListRays_,{n1_,n2_,n3_,n4_},m_] extracts the list of distinct trees with given dimension vector  ";


(* from CoulombHiggs *)
F0HiggsBranchFormula::usage = "F0HiggsBranchFormula[Mat_,Cvec_,Nvec_] computes the refined index of a quiver with DSZ matrix Mat, FI parameters Cvec, dimension vector Nvec using Reineke's formula. Assumes that the quiver has no oriented loops.";
F0RationalInvariant::usage = "F0RationalInvariant[Mat_,Cvec_,Nvec_,y_] computes the rational invariant of a quiver with DSZ matrix Mat, dimension vector Nvec and FI parameters Cvec computed using Reineke's formula.";
F0StackInvariant::usage = "F0StackInvariant[Mat_,Cvec_,Nvec_,y_] gives the stack invariant of a quiver with DSZ matrix Mat, dimension vector Nvec and FI parameters Cvec computed using Reineke's formula.";
F0QDeformedFactorial::usage = "F0QDeformedFactorial[n_,y_] gives the q-deformed factorial ";
F0ListAllPartitions::usage = "F0ListAllPartitions[gam_] returns the list of unordered partitions of the positive integer vector gam as a sum of positive integer vectors "; 
F0BinarySplits::usage="F0BinarySplits[Nvec_] gives the list of dimension vectors which are smaller than Nvec/2";
EvaluateKronecker::usage = "EvaluateKronecker[f_] replaces each \!\(\*SubscriptBox[\(Kr\), \(m_\)]\)[p_,q_] with the index of the Kronecker quiver with m arrows and dimension vector {p,q}";


(* Large volume scattering diagram *)
xytost::usage = "xytost[{x_,y_},m:0]:={x,Sqrt[x^2+m x+2y]}";
sttoxy::usage = "sttoxy[{s_,t_},m:0]:={s,-1/2(m s+s^2-t^2)}";
Rayxy::usage = "Rayxy[{r_,d1_,d2_,ch2_},x_,L_,m_] gives Graphics directive for ray starting from (x,?) if r<>0, or (?,x) if r=0, extending by distance L"; 
IntersectRays::usage = "IntersectRays[{r_,d1_,d2_,ch2_},{rr_,dd1_,dd2_,cch2_},m_] returns intersection point (x,y) of two rays, or {} if they are collinear;
  IntersectRays[{r_,d1_,d2_,ch2_},{rr_,dd1_,dd2_,cch2_},z_,zz_,m_]returns intersection point (x,y) of two rays if the intersection point lies upward from z and z', or {} otherwise";  
IntersectRaysNoTest::usage = "IntersectRays[{r_,d1_,d2_,ch2_},{rr_,dd1_,dd2_,cch2_},z_,zz_,m_]returns intersection point (x,y) of two rays if the intersection point lies strictly upward from z and z', or {} otherwise, without testing non-vanishing of DSZ product";  
IntersectRaysSt::usage = "IntersectRaysSt[{r_,d1_,d2_,ch2_},{rr_,dd1_,dd2_,cch2_},psi_,m_] returns intersection point (s,t) for rotated scattering rays, or {} if they are collinear";
InitialPosition::usage="InitialPosition[{r_,d1_,d2_,ch2_},m_] computes the s-value where the ray intersects the real axis"; 
TestBranch::usage = "TestBranch[{r_,d1_,d2_,ch2_},s_]  tests if (s,.) is on the branch with ImZ>0";
Rays::usage = "Rays[{r_,d1_,d2_,ch2_},t_,psi_,m_] computes  s(t) for the ray Re[E^(-I psi) Z_gam]\[Equal]0 with LV central charge"; 
WallCircle::usage = "WallCircle[{r_,d1_,d2_,ch2_},{rr_,dd1_,dd2_,cch2_},m_] constructs the graphics directive for the wall of marginal stability in (s,t) plane";
WallRadius::usage = "WallRadius[{r_,d1_,d2_,ch2_},{rr_,dd1_,dd2_,cch2_},m_] computes the radius SQUARE of the wall of marginal stability in (s,t) plane";

ScattDiag::usage= "ScattDiag[TreeList_,m_] draws scattering diagram in (x,y) plane for each tree in Treelist"; 
ScattDiagInternal::usage= "ScattDiagInternal[Tree_,m_] constructs total charge, coordinate of root and list of line segments in (x,y) coordinates, {min(x), max(x)}; used by ScattDiag"; 
ScattDiagSt::usage = "ScattDiagSt[TreeList_] overlays scattering diagrams in (s,t) plane for all trees in TreeList;
  ScattDiagSt[TreeList_,psi_] does the same for rotated scattering diagrams";
ScattDiagInternalSt::usage= "ScattDiagInternalSt[Tree_,m_] constructs total charge, coordinate of root and list of line segments in (s,t) coordinates, {min(s), max(s)}; used by ScattDiag;
  ScattDiagInternalSt[Tree_,psi_] does the same for rotated scattering diagram"; 
ScattDiagLV::usage = "ScattDiagLV[TreeList_,psi_,m_] overlays scattering diagrams in (s,t) plane for all trees in TreeList, using exact hyperbolaes for rays";
ScattDiagInternalLV::usage= "ScattDiagInternalSt[Tree,psi_,Styl_,m_] constructs total charge, coordinate of root and list of line segments in (s,t) coordinates, {min(s), max(s)}, using PlotStyle->Styl for plotting rays ; used by ScattDiagLV;"; 
InitialLabelPosition::usage= "InitialLabelPosition[m_] returns a position (s,t) for the label for an initial ray with slope m; this position is updated on each call using variables LiSlopes and LiHeights"; 
ScattIndexImproved::usage="ScattIndexImproved[TreeList_, opt_] computes the index for each tree in TreeList, taking care of non-primitive internal states";
ScattIndexImprovedInternal::usage="ScattIndexImprovedInternal[Tree_, opt_] computes the index for Tree, taking care of non-primitive internal states";
(* computing the index *)
ScattIndex::usage = "ScattIndex[TreeList_] computes the index for each tree in TreeList; do not trust the result if internal lines have non-primitive charges";
ScattIndexInternal::usage = "ScattIndexInternal[Tree_] computes {total charge, list of Kronecker indices associated to each vertex in Tree}";
BinaryTreeIndex::usage="BinaryTreeIndex[TreeList_] computes the index of each binary tree in TreeList";
BinaryTreeIndexInternal::usage="BinaryTreeIndexInternal[Tree_] produces a list of wall-crossing factors from each vertex in a binary tree";



MonodromyOnCharge::usage = "MonodromyOnCharge[M_,{r_,d1_,d2_,ch2_}] computes the image of charge vector under sequence of monodromies";
MonodromyOnTau::usage = "MonodromyOnTau[M_,tau_] computes the image of tau under sequence of monodromies";
MonodromyFromSphericalTwist::usage = "MonodromyFromSphericalTwist[{r_,d1_,d2_,ch2_}] computes the action of spherical twist in period basis";
(* Pi stability *)
FundDomainC::usage = "FundDomainC[k_] produces the Graphics directives for the fundamental domain of Gamma_1(4) on the interval [k,k+1]";
ConvergenceDomain::usage ="ConvergenceDomain[m_] plots the convergence domain of the expansions around tau=Infinity, 0, 1/2, 1";
ToFundDomain::usage = "ToFundDomain[tau_] produces {tau',M'} such that M'.tau'=tau and tau' lies in fundamendal domain of Gamma_1(4)";
ToFundDomainApprox::usage = "ToFundDomainApprox[tau_,precision_] produces {tau',M'} such that M'.tau'=tau and tau' lies in fundamendal domain of Gamma_1(4) centered around conifold";
UHPGeodesic::usage="UHPGeodesic[tau,tau'] draws the upper-half-plane geodesic between these two points (which may lie on the boundary)";

QuantumVolume::usage = "QuantumVolume[m_] gives the period T(0)/I";
QuantumVolumep::usage = "QuantumVolume[m_] gives the period T(1/2)/I";
Eichlerf1::usage = "Eichlerf1[tau_,m_]:=Sum[FourierC[m][n+1]]/n Exp[2Pi I n tau],n>0]";
Eichlerf2::usage ="Eichlerf2b[tau_,m_]:=-Sum[FourierC[m][[n+1]]/n^2 Exp[2Pi I n tau],n>0]";
Eichlerf1p::usage = "Eichlerf1p[taup_,m_]:=Sum[FourierCp[m][[n+1]]/n Exp[2Pi I n taup],n>0]";
Eichlerf2p::usage = "Eichlerf2p[taup_,m_]:=2Pi I(taup)Eichlerf1p[taup,m]-Sum[FourierCp[m][[n+1]]/n^2 Exp[2Pi I n taup],n>0]";
Eichlerf1pp::usage = "Eichlerf1pp[taup_,m_]:=Sum[FourierCpp[m][[n+1]]/n Exp[2Pi I n taup],n>0]";
Eichlerf2pp::usage = "Eichlerf2pp[taup_,m_]:=2Pi I(taup)Eichlerf1pp[taup,m]-Sum[FourierCpp[m][[n+1]]/n^2 Exp[2Pi I n taup],n>0]";
EichlerTLV::usage= "EichlerTLV[tau_,m_] gives numerical value of T[tau] using Eichler integral near LV point";
EichlerTDLV::usage= "EichlerTDLV[tau_,m_] gives numerical value of T[tau] using Eichler integral near LV point";
EichlerTp::usage = "EichlerTp[tau_,m_] gives numerical value of T[tau] using Eichler integral near tau=0";
EichlerTDp::usage = "EichlerTDp[tau_,m_] gives numerical value of T[tau] using Eichler integral near tau=0";
EichlerTpp::usage = "EichlerTpp[tau_,m_] gives numerical value of T[tau] using Eichler integral near tau=1/2";
EichlerTDpp::usage = "EichlerTDpp[tau_,m_] gives numerical value of T[tau] using Eichler integral near tau=1/2";
om2B::usage = "om2B[tau_,m_,NMax_:20] computes om2B using the 1/z expansion to order NMax";
om3B::usage = "om3B[tau_,m_,NMax_:20] computes om3B using the 1/z expansion to order NMax";
EichlerTB::usage= "EichlerTB[tau_,m_] gives numerical value of T[tau] using expansion near branch point";
EichlerTDB::usage= "EichlerTDB[tau_,m_] gives numerical value of T[tau] using expansion near branch point";
EichlerTBtauR::usage = "EichlerTBtauR[tb_] provides the radius of convergence of the expansion in tau/(1-2tau) near the branch point tau=tb"

EichlerTTDLVcut::usage= "EichlerTTDDLVcut[tau_,m_] gives numerical value of {T[tau],TD[tau]} using Eichler integral near LV point";
EichlerTTDpcut::usage = "EichlerTTDpcut[tau_,m_] gives numerical value of {T[tau],TD[tau]} using Eichler integral near tau=0";
EichlerTTDppcut::usage = "EichlerTTDppcut[tau_,m_] gives numerical value of {T[tau],TD[tau]} using Eichler integral near tau=1/2";
EichlerTTDpppcut::usage = "EichlerTTDpppcut[tau_,m_] gives numerical value of {T[tau],TD[tau]} using Eichler integral near tau=1";
EichlerTTDBcut::usage = "EichlerTTDBcut[tau_,m_] gives numerical value of {T[tau],TD[tau]} using Eichler integral near tau=tauB";
SwitchMin::usage="SwitchMin[{\!\(\*SubscriptBox[\(x\), \(1\)]\),\[Ellipsis],\!\(\*SubscriptBox[\(x\), \(n\)]\)},\!\(\*SubscriptBox[\(code\), \(1\)]\),\[Ellipsis],\!\(\*SubscriptBox[\(code\), \(n\)]\)] runs the \!\(\*SubscriptBox[\(code\), \(i\)]\) corresponding to the \!\(\*SubscriptBox[\(x\), \(i\)]\) that is minimum";

EichlerZ::usage = "EichlerZ[{r_,d1_,d2_,ch2_},tau_,m_] gives numerical value of central charge at tau=t1+I t2, by mapping back to fundamental domain F_C";
EichlerZp::usage = "EichlerZp[{r_,d1_,d2_,ch2_},taup_,m_] gives numerical value of central charge at tau', using expansion around conifold";
EichlerZpp::usage = "EichlerZpp[{r_,d1_,d2_,ch2_},taupp_,m_] gives numerical value of central charge at tau'' sing expansion around dual conifold";
EichlerT::usage= "EichlerT[tau_,m_] gives numerical value of T[tau] using Eichler integral, by mapping back to fundamental domain F_C";
EichlerTD::usage= "EichlerTD[tau_,m_] gives numerical value of TD[tau] using Eichler integral, by mapping back to fundamental domain F_C";

EichlerZ0::usage = "EichlerZ0[{r_,d1_,d2_,ch2_},tau_,m_] gives numerical value of central charge at tau=t1+I t2, by mapping back to fundamental domain F_C, not using expansion around tauB";
EichlerT0::usage= "EichlerT0[tau_,m_] gives numerical value of T[tau] using Eichler integral, by mapping back to fundamental domain F_C, not using expansion around tauB";
EichlerTD0::usage= "EichlerTD0[tau_,m_] gives numerical value of TD[tau] using Eichler integral, by mapping back to fundamental domain F_C, not using expansion around tauB";

EichlerInt::usage= "EichlerInt[tau0_,tau1_,m_] computes the integral int {1,tau} C(tau) dtau along a straight line from tau0 to tau1";

MeijerT0::usage ="MeijerT0[tau_] computes T[tau] at m=0 using MeijerG representation";
MeijerTD0::usage ="MeijerTD0[tau_] computes T[tau] at m=0 using MeijerG representation";
MeijerT::usage ="MeijerT0[tau_] computes T[tau] at m=1/2 using MeijerG representation";
MeijerTD::usage ="MeijerTD0[tau_] computes T[tau] at m=1/2 using MeijerG representation";



tauB::usage="tauB[m_] computes the location of branch point in tau plane";
InversetauB::usage="InversetauB[tau_] computes the value of m leading to desired tauB[m]";

InverseJ::usage= "InverseJ[j_] computes tau such that J[tau]=j";
InverseJ4::usage= "InverseJ4[j_]Computes tau such that J4[tau]=j";
ModularJ4::usage= "ModularJ4[tau_] computes 8+(DedekindEta[tau]/DedekindEta[4tau])^8";

DtauT::usage="DtauT[\[Tau]_,m_] is T'[\[Tau]], the derivative of T with respect to \[Tau], calculated efficiently.";
DtauZ::usage="DtauZ[\[Gamma]_,\[Tau]_,m_] is \!\(\*SubscriptBox[\(Z\), \(\[Gamma]\)]\)'[\[Tau]], namely D[EichlerZ[\[Gamma],\[Tau]],\[Tau]], calculated efficiently.";
UnitDtauZ::usage="UnitDtauZ[\[Gamma]_,\[Tau]_,m_] is Normalize[\!\(\*SubscriptBox[\(Z\), \(\[Gamma]\)]\)'[\[Tau]]].";

ArgDtauT::usage="ArgDtauT[\[Tau]_,m_] is the argument of T'[\[Tau]], between -Pi and Pi.";
ArgDtauTD::usage="ArgDtauTD[\[Tau]_,m_] is the argument of \!\(\*SubscriptBox[\(T\), \(D\)]\)'[\[Tau]], between -Pi and Pi.";
UnitDtauT::usage="UnitDtauT[\[Tau]_,m_] is Normalize[T'[\[Tau]]], namely T'[\[Tau]]/Abs[T'[\[Tau]]].";
UnitDtauTD::usage="UnitDtauTD[\[Tau]_,m_] is Normalize[\!\(\*SubscriptBox[\(T\), \(D\)]\)'[\[Tau]]].";

IntegralCurve::usage="IntegralCurve[tauinit_,tangent_,{ainit_,amin_,amax_},boundaries_] produces a function on [0,1] (it normalizes the domain).  Starting at tauinit at a=ainit, following the tangent direction (an expression in tau) and stopping at the boundaries (default: {Im[tau]==0.01}).  The range of integration parameters {amin,amax} can be infinite provided the actual rays remain finite by hitting the boundaries.";
NormalizeFunctionDomain::usage="NormalizeFunctionDomain[fun_InterpolatingFunction] rescales the argument of fun to interval [0,1]";
NormalizeApprox::usage="NormalizeApprox[z_,eps_] normalizes z\[Element]\!\(\*TemplateBox[{},\n\"Complexes\"]\) (complex) to roughly unit length for large z, but behaves smoothly near zero.";
RayCh::usage="RayCh[\[Psi]_,m_] gives an initial ray starting at 0, namely a function [0,1]\[RightArrow]\[DoubleStruckCapitalH] starting (close to) 0 and following a ray where \!\(\*SubscriptBox[\(Z\), \([1, 0, 0]\)]\)=-\!\(\*SubscriptBox[\(T\), \(D\)]\) has phase (\[Psi]+\[Pi]/2 mod \[Pi]).  Shifting \[Psi] by 2\[Pi] gives a different ray, corresponding to a homological shift by 2.  Values are cached.";
RayGeneral::usage="RayGeneral[\[Gamma]_,\[Psi]_,m_,tauexpr_,start_List] gives a ray where \!\(\*SubscriptBox[\(Z\), \(\[Gamma]\)]\) has phase (\[Psi]+\[Pi]/2 mod \[Pi]).  The starting point is obtained using FindRoot[\!\(\*SubscriptBox[\(\[Ellipsis]Z\), \(\[Gamma]\)]\)[tauexpr]\[Ellipsis],start], see documentation of FindRoot.";

IntersectExactRaysLV::usage="IntersectExactRaysLV[{r_,d1_,d2_,ch2_},{rr_,dd1_,dd2_,cch2_},psi_,m_]returns value of tau of intersection point of two exact rays using EichlerTLV to evaluate the periods, or {} if they are collinear";
IntersectExactRaysC::usage="IntersectExactRaysLV[{r_,d1_,d2_,ch2_},{rr_,dd1_,dd2_,cch2_},psi_,m_]returns value of tau of intersection point of two exact rays using EichlerTp to evaluate the periods, or {} if they are collinear";
CriticalPsi::usage = "CriticalPsi[m_]:=Values of psi such that {V_psi=0,tildeV_psi=1-mu,mu=0}";
XY::usage = "XY[tau_,psi_,m_] computes the coordinates {x,y} such that scattering  rays are straight lines ry+(d1+d2)x+(d2 m-ch2)=0";



(* Mapping local F1 to local F0 *)
F1ToF0::usage = "F1ToF0[{r_,dH_,dC_,ch2_}]:={r,1/2(dH-dC),dH+dC,ch2}";
F0ToF1::usage = "F0ToF1[{r_,d1_,d2_,ch2_}]:={r,d1+1/2 d2,-d1+1/2 d2,ch2}";
repCh1::usage = "Replaces Ch1[mH_,mC_] by {1,mH,mC,1/2 (mH^2-mC^2)},GV[mH_,mC_,n_]:>{0,mH,mC,n}";
ChernToCh1::usage="ChernToCh1[f_] replaces {r,dH,dC,ch2} by r Ch1[dH/r,dC/r] whenever the discriminant vanishes";
DSZ1::usage = "DSZ1[gam1_,gam2_] computes the DSZ pairing between F1 Chern vectors";
ZLV1::usage = "ZLV1[{r_,dH_,dC_,ch2_},{s_,t_},m_] gives the large volume central charge on local F1";
Euler1::usage = "Euler1[{r_,dH_,dC_,ch2_},{rr_,ddH_,ddC_,cch2_}] computes the Euler form on D(F_1)";
Disc1::usage = "Disc1[{r_,dH_,dC_,ch2_}] computes the discriminant -ch2/r+(dH^2-dC^2)/(2r^2)";
DimGieseker1::usage = "DimGieseker1[{r_,dH_,dC_,ch2_}] computes expected dimension of moduli space of Gieseker-semi stable sheaves";


MutateCollection1::usage="MutateCollection1[Coll_,klist_] acts on the list of Chern vectors Coll by the successive mutations in klist, which is a list of {node number,sign}, with sign=1 for right mutation, -1 for left mutation";
ExtFromStrong1::usage="ExtFromStrong1[Coll_] computes the Chern vectors of the objects in the Ext collection dual to a strong collection"; 
StrongFromExt1::usage="ExtFromStrong1[Coll_] computes the Chern vectors of the objects in the strong collection dual to an Ext collection"; 
SpecFlow1::usage = "SpecFlow1[{r_,dH_,dC_,ch2_},{mH_,mC_}] gives CHern vector of tensor product by O(mH,mC)";
ZLV1::usage = "ZLV1[{r_,dH_,dC_,ch2_},{s_,t_},m_] computes the large volume central charge Z for F1";
InitialPosition1::usage="InitialPosition1[{r_,dH_,dC_,ch2_},m_] computes the s-value where the ray intersects the real axis"; 
QuiverDomain1::usage="QuiverDomain[Coll_,psi_,m_,OptionsPattern[]] plots the Region where the LV central charges of the vectors in Coll have Re[e^{-I psi}Z]<0, and the region where they are in same half-plane.";
StabilityRegion1::usage="StabilityRegion1[gam1_,gam2_,m_] plots the region where Im[Z[gam1]Conjugate[Z[gam2]] DSZ1[gam1,gam2]>0"; 
Wall1::usage = "Wall1[{r_,dH_,dC_,ch2_},{rr_,ddH_,ddC_,cch2_},{s_,t_},m_] computes Im[Z[gamma] Conjugate[Z[ggamma]]/t";
GenSlope1::usage = "GenSlope1[{r_,dH_,dC_,ch2_},m_] computes the slope (3dH+dC)/(4r) if r<>0, or s-value if r=0";
IntersectRays1::usage = "IntersectRays1[{r_,dH_,dC_,ch2_},{rr_,ddH_,ddC_,cch2_},m_] returns intersection point (x,y) of two rays, or {} if they are collinear;
  IntersectRays1[{r_,dH_,dC_,ch2_},{rr_,ddH_,ddC_,cch2_},z_,zz_,m_]returns intersection point (x,y) of two rays if the intersection point lies upward from z and z', or {} otherwise";  
IntersectRaysNoTest1::usage = "IntersectRays1[{r_,dH_,dC_,ch2_},{rr_,ddH_,ddC_,cch2_},z_,zz_,m_]returns intersection point (x,y) of two rays if the intersection point lies strictly upward from z and z', or {} otherwise, without testing non-vanishing of DSZ product";  
IntersectRaysSt1::usage = "IntersectRaysSt1[{r_,dH_,dC_,ch2_},{rr_,ddH_,ddC_,cch2_},psi_,m_] returns intersection point (s,t) for rotated scattering rays, or {} if they are collinear";

ScanConstituents1::usage = "ScanConstituents[gam_,{smin_,smax_},{n_,np_},m_,phimax_] searches possible list of constituents \[PlusMinus]O(k,k), \[PlusMinus]O(k+1,k) with slope in [smin,smax], number of D/Dbar constituents<=(n,np), cost function less than phimax and charges adding up to gam";
ListStableTrees1::usage = "ListStableTrees[LiCh_,{s0_,t0_},m_] constructs consistent trees from constituents list LiCh={k_i Ch[m_i,n_i]}, which are stable at (s0,t0)";
ScanAllTrees1::usage = "ScanAllTrees[{r_,dH_,dC_,ch2_},{s0_,t0_},n_,m_] constructs all possible trees with up to (n,n) left- and right- moving charges adding up to [r,d1,d2,ch2) leading to an outgoing ray through the point (s,t)";
ContractInitialRays1::usage="ContractInitialRays1[Trees_,m_] replaces primary scatterings in LV+eps diagram by their total charge";
ScattGraph1::usage = "ScattGraph1[Tree_,m_:1/2] extracts the list of vertices and adjacency matrix of Tree";
ScattGraphInternal1::usage = "ScattGraphInternal1[Tree_,m_:1/2] is a light version of ScattDiagInternal1 used to construct the scattering graph";
ScattIndexImproved1::usage="ScattIndexImproved1[TreeList_, opt_] computes the index for each tree in TreeList, taking care of non-primitive internal states";
ScattIndexImprovedInternal1::usage="ScattIndexImprovedInternal1[Tree_, opt_] computes the index for Tree, taking care of non-primitive internal states";
ScattCheck1::usage = "ScattCheck1[Tree_,m_]returns {charge,{x,y}} of the root vertex if Tree is consistent, otherwise {total charge,{}}";
ScattSort1::usage = "ScattSort1[LiTree_,m_] sorts trees in LiTree by growing radius";
(* computing the index *)
ScattIndex1::usage = "ScattIndex[TreeList_] computes the index for each tree in TreeList; do not trust the result if internal lines have non-primitive charges";
ScattIndexInternal1::usage = "ScattIndexInternal[Tree_] computes {total charge, list of Kronecker indices associated to each vertex in Tree}";
GCD1::usage="GCD1[gam_] gives the largest d such that gam/d belongs to the charge lattice";
CostPhi1::usage="CostPhi1[{r_,dH_,dC_,ch2_},s_,m_]:=1/2(3dH+dC)-r(2s+m)";
ConstructLVDiagram1::usage="ConstructLVDiagram1[smin_,smax_,phimax_,Nm_,m_,ListRays_] constructs the LV scattering diagram for F1 with initial rays in the interval [smin,smax], cost function up to phi, scattering products with n1+n2<=Nn at each intersection; m is assumed to be real; The output consists of a list of  {charge, {x,y}, parent1, parent2,n1,n2 }; If ListRays is not empty, then uses it as initial rays.";
LVTreesFromListRays1::usage="lVTreesFromListRays1[ListRays_,{r_,dH_,dC_,ch2_},m_] extract the trees with given charge in the List of rays, constructed by ConstructLVDiagram1";



Begin["`Private`"]

General::argerr="Wrong arguments in the expression `1`";
General::notimpl="`1` not implemented";

Euler[{rk_,c11_,c12_,ch2_},{rrk_,cc11_,cc12_,cch2_}]:=-c11 cc12-cc11 c12+rk (cc11+cc12) +cch2 rk-rrk(c11+c12)  +ch2 rrk+rk rrk;
SpecFlow[{rk_,c11_,c12_,ch2_},{m1_,m2_}]:={rk,c11+rk m1,c12+rk m2,c12 m1+c11 m2+m1 m2 rk+ch2};
Disc[{rk_,c11_,c12_,ch2_}]:=-ch2/rk+c11 c12/rk^2;
DSZ[{rk_,c11_,c12_,ch2_},{rrk_,cc11_,cc12_,cch2_}]:=2rk (cc11+cc12)-2 rrk (c11+c12);
DimGieseker[{rk_,c11_,c12_,ch2_}]:=1+2 c11 c12-2 ch2 rk-rk^2;
repCh={Ch[m1_,m2_][1]:>-{1,m1,m2,m1 m2},Ch[m1_,m2_]:>{1,m1,m2,m1 m2},GV[m1_,m2_,n_]:>{0,m1,m2,n}};
ChernToCh[f_]:=f/.{{r_Integer,d1_Integer,d2_Integer,ch2_}:>If[r==0,GV[d1,d2,ch2],If[Disc[{r,d1,d2,ch2}]==0,If[r>0,r Ch[d1/r,d2/r],-r Ch[d1/r,d2/r][1]],{r,d1,d2,ch2}]]};

MutateCollection[Coll_,klist_]:=Module[{Coll0,k,eps},
(* Coll is a list of Chern vectors, klist a list of {node,\pm 1} *)
Coll0=If[Length[klist]>1, MutateCollection[Coll,Drop[klist,-1]], Coll];
k=Last[klist][[1]]; eps=Last[klist][[2]];
Table[If[i==k,-Coll0[[k]],Coll0[[i]]+Max[0,eps DSZ[Coll0[[i]],Coll0[[k]]]]Coll0[[k]]],{i,Length[Coll0]}]]

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

DimFromCh[Coll_,gam_]:=(gam/.repCh/.repCh1) . Inverse[Coll];
ChFromDim[Coll_,Nvec_]:=Flatten[{Nvec} . Coll];




ZLV[{r_,d1_,d2_,ch2_},{s_,t_},m_]:=-r (s+I t)(s+I t+m)+d1 (s+I t)+d2(s+I t+m)-ch2;

MLV=\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{"1", "0", "0", "0"},
{"0", "1", "0", "0"},
{"1", "0", "1", "0"},
{
FractionBox["1", "2"], 
FractionBox["1", "2"], "1", "1"}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\);MLR1=\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{"1", "0", "0", "0"},
{"1", "1", "0", "0"},
{
RowBox[{"-", "1"}], "0", "1", "0"},
{"0", 
RowBox[{"-", 
FractionBox["1", "2"]}], 
RowBox[{"-", 
FractionBox["1", "2"]}], "1"}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\);MLR2=\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{"1", "0", "0", "0"},
{
RowBox[{"-", "1"}], "1", "0", "0"},
{"0", "0", "1", "0"},
{"0", "0", 
RowBox[{"-", 
FractionBox["1", "2"]}], "1"}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\); MCon=\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{"1", "0", "0", "0"},
{"0", "1", "0", "0"},
{"0", "0", "1", 
RowBox[{"-", "4"}]},
{"0", "0", "0", "1"}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\);MOrb=\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{"1", "0", "0", "0"},
{"0", "1", "0", "0"},
{"1", 
RowBox[{"-", "2"}], 
RowBox[{"-", "3"}], "4"},
{
FractionBox["1", "2"], 
RowBox[{"-", 
FractionBox["1", "2"]}], 
RowBox[{"-", "1"}], "1"}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\);Sigma=\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{"0", "0", "0", 
RowBox[{"-", "2"}]},
{"0", "0", "1", "0"},
{"0", "1", "1", "0"},
{
RowBox[{"-", "1"}], "0", "0", "0"}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\);
MB[k_]:=\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{"1", "0", "0", "0"},
{"0", "1", "0", "0"},
{
RowBox[{"1", "+", 
RowBox[{"2", " ", "k"}]}], "0", 
RowBox[{"-", "1"}], "0"},
{
RowBox[{
FractionBox["1", "2"], "+", "k", "+", 
SuperscriptBox["k", "2"]}], 
RowBox[{
FractionBox["1", "2"], "+", "k"}], "0", 
RowBox[{"-", "1"}]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\); MBt[k_]:=\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{"1", "0", "0", "0"},
{"0", "1", "0", "0"},
{
RowBox[{"1", "+", 
RowBox[{"2", " ", "k"}]}], 
RowBox[{"-", "2"}], 
RowBox[{"-", "1"}], "0"},
{
RowBox[{
FractionBox["1", "2"], "+", "k", "+", 
SuperscriptBox["k", "2"]}], 
RowBox[{
RowBox[{"-", 
FractionBox["1", "2"]}], "-", "k"}], "0", 
RowBox[{"-", "1"}]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\);MFB={{1,0,0,0},{0,-1,0,0},{0,1,1,0},{0,0,0,1}};
MZ4[k_]:=\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{"1", "0", "0", "0"},
{"1", 
RowBox[{"-", "1"}], "0", "0"},
{
RowBox[{
RowBox[{"-", "2"}], " ", "k", " ", 
RowBox[{"(", 
RowBox[{"1", "+", "k"}], ")"}]}], 
RowBox[{"1", "+", 
RowBox[{"2", " ", "k"}]}], 
RowBox[{"3", "+", 
RowBox[{"4", " ", "k"}]}], 
RowBox[{"-", "4"}]},
{
RowBox[{
RowBox[{"-", "2"}], " ", "k", " ", 
SuperscriptBox[
RowBox[{"(", 
RowBox[{"1", "+", "k"}], ")"}], "2"]}], 
RowBox[{
FractionBox["1", "2"], " ", 
SuperscriptBox[
RowBox[{"(", 
RowBox[{"1", "+", 
RowBox[{"2", " ", "k"}]}], ")"}], "2"]}], 
RowBox[{
FractionBox["5", "2"], "+", 
RowBox[{"6", " ", "k"}], "+", 
RowBox[{"4", " ", 
SuperscriptBox["k", "2"]}]}], 
RowBox[{
RowBox[{"-", "3"}], "-", 
RowBox[{"4", " ", "k"}]}]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\);
(* old definition, fixed point at m=-1/2: 
MZ4t[k_]:={{1,0,0,0},{-1,-1,0,0},{1-2 k (1+k),3+2 k,3+4 k,-4},{-2 k (1+k)^2,2 (1+k)^2,5/2+6 k+4 k^2,-3-4 k}}; *)

MZ4t[k_]:=\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{"1", "0", "0", "0"},
{"1", 
RowBox[{"-", "1"}], "0", "0"},
{
RowBox[{
RowBox[{"-", "2"}], " ", 
SuperscriptBox["k", "2"]}], 
RowBox[{"1", "+", 
RowBox[{"2", " ", "k"}]}], 
RowBox[{"1", "+", 
RowBox[{"4", " ", "k"}]}], 
RowBox[{"-", "4"}]},
{
RowBox[{
RowBox[{"-", 
SuperscriptBox["k", "2"]}], " ", 
RowBox[{"(", 
RowBox[{"1", "+", 
RowBox[{"2", " ", "k"}]}], ")"}]}], 
RowBox[{
FractionBox["1", "2"], "+", "k", "+", 
RowBox[{"2", " ", 
SuperscriptBox["k", "2"]}]}], 
RowBox[{
FractionBox["1", "2"], "+", 
RowBox[{"2", " ", "k"}], "+", 
RowBox[{"4", " ", 
SuperscriptBox["k", "2"]}]}], 
RowBox[{
RowBox[{"-", "1"}], "-", 
RowBox[{"4", " ", "k"}]}]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\);


SameHalfPlaneQ[{}]:=True;
SameHalfPlaneQ[Zlist_List]:=If[AnyTrue[Zlist,#==0&],False,-Subtract@@MinMax[Arg[Zlist/Zlist[[1]]]]<Pi];



(* Fundamental domain centered around k-1/2 *)
FundDomainC[k_]:=Graphics[{Thick,Blue,Line[{{k,0},{k,1}}],Line[{{k+1,0},{k+1,1}}],
Circle[{k+1/4,0},1/4,{0,Pi}],
Circle[{k+3/4,0},1/4,{0,Pi}],Hue[.9],
Circle[{k+1/8,0},1/8,{0,Pi}],
Circle[{k+7/24,0},1/24,{0,Pi}],
Circle[{k+5/12,0},1/12,{0,Pi}],
Circle[{k+7/8,0},1/8,{0,Pi}],
Circle[{k+17/24,0},1/24,{0,Pi}],
Circle[{k+7/12,0},1/12,{0,Pi}],
Black,
Text["LV",{k+1/2,.9}],Text["C",{k+.03,0.005}],Text["C'",{k+1/2+.03,0.005}],
Text["C",{k+1-.03,0.005}],Text["LV",{k+1/4-.03,.005}],Text["LV",{k+3/4+.03,.005}],
Text["C",{k+1/3+.03,0.005}],Text["C",{k+2/3-0.03,0.005}](*,
Dotted,Circle[{k-1,0},3/2,{0,Pi/9}],
Text["B",{k+.4,.55}]*)
}];

ConvergenceDomain[m_]:=Module[{tb},
tb=tauB[m];Show[
RegionPlot[Im[(t1+I t2)]>Im[tb],{t1,0,1},{t2,0,1},PlotStyle->{Blue,Opacity[.1]}],RegionPlot[Im[-1/4/(t1+I t2)]>Max[Im[-1/4/(tb-1)],Im[-1/4/(tb)]],{t1,0,1},{t2,0,1},PlotStyle->{Red,Opacity[.1]}],
RegionPlot[Im[-1/4/(t1+I t2-1)]>Max[Im[-1/4/(tb-1)],Im[-1/4/(tb)]],{t1,0,1},{t2,0,1},PlotStyle->{Red,Opacity[.1]}],
RegionPlot[Im[(t1+I t2)/(1-2(t1+I t2))]>Im[tb/(1-2tb)],{t1,0,1},{t2,0,1},PlotStyle->{Green,Opacity[.1]}],
Graphics[{PointSize[.02],Point[ReIm[tb]],Point[ReIm[(tb-1)/(4tb-3)]],Point[ReIm[(3tb-1)/(4tb-1)]]}],Graphics[{Dashed,UHPGeodesic[1/2,tb],UHPGeodesic[1/2,(tb-1)/(4tb-3)],UHPGeodesic[1/2,(3tb-1)/(4tb-1)],Text["B",ReIm[tb+.03 I]]}],PlotRange->{{0,1},{0,1}}]];

UHPGeodesic[DirectedInfinity[I],tau_]:=HalfLine[ReIm[tau],{0,1}];
UHPGeodesic[tau_,DirectedInfinity[I]]:=HalfLine[ReIm[tau],{0,1}];
UHPGeodesic[tau_/;Im[tau]>=0,taup_/;Im[taup]>=0]:=
If[Re[tau]==Re[taup],Line[{ReIm[tau],ReIm[taup]}],
With[{x0=(Abs[tau]^2-Abs[taup]^2)/(2 (Re[tau]-Re[taup]))},
Circle[{x0,0},Abs[tau-x0],{Arg[tau-x0],Arg[taup-x0]}]]];
UHPGeodesic[tau_,taup:{_,_}]:=UHPGeodesic[tau,{1,I} . taup];
UHPGeodesic[tau:{_,_},taup_]:=UHPGeodesic[{1,I} . tau,taup];
call:UHPGeodesic[tau_,tau_]:=With[{m=Message[General::argerr,HoldForm[call]]},Null/;m];
call:UHPGeodesic[___]:=With[{m=Message[General::argerr,HoldForm[call]]},Null/;m];


MonodromyOnCharge[M_,{r_,d1_,d2_,ch2_}]:=Module[{Sig={{0,0,0,-2},{0,0,1,0},{0,1,1,0},{-1,0,0,0}}},If[Length[M]==0,{r,d1,d2,ch2},If[
     ArrayDepth[M]==2,
(* a single matrix *)
{r,d1,d2,ch2} . Sig . Inverse[M] . Inverse[Sig],
(* a list of matrices *)
MonodromyOnCharge[Drop[M,1],{r,d1,d2,ch2} . Sig . Inverse[First[M]] . Inverse[Sig]]]]
];

MonodromyOnTau[M_,tau_]:=If[Length[M]==0,tau,If[
ArrayDepth[M]==2,
(* a single matrix *)
(M[[4,4]]tau+M[[4,3]])/(M[[3,4]]tau+M [[3,3]])
,(* a list of matrices *)
MonodromyOnTau[Drop[M,1],(M[[1,4,4]]tau+M[[1,4,3]])/(M[[1,3,4]]tau+M [[1,3,3]])]]
];

MonodromyFromSphericalTwist[{r_,d1_,d2_,ch2_}]:=Inverse[Sigma] . Inverse[Transpose[IdentityMatrix[4]-Transpose[{{r,d1,d2,ch2}}] . {{-2d1 -2d2,2r,2r,0}}]] . Sigma;

ToFundDomainApprox[tau_]:={#[[1,1]],Rationalize[#[[2;;]]]}&[ToFundDomainAux[tau]];
ToFundDomain = ToFundDomainApprox;
ToFundDomainAux=Compile[{{ttau,_Complex}},
  Module[{tau=ttau,M=1.{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}},fl,taupp,flpp},
    Do[
      fl=Floor[Re[tau]];
      tau=tau-fl;
      M=M . {{1,0,0,0},{0,1,0,0},{fl,0,1,0},{fl^2/2,fl/2,fl,1}};
      taupp=.25(.5-tau)^-1+.5;
      If[0<=Re[taupp]<=1,Break[]];
      flpp=Floor[Re[taupp]];
      tau=.5-.25(taupp-flpp-.5)^-1;
      M=M . ({{1,0,0,0},{0,1,0,0},{1/2,-1/2,0,0},{1/4,0,0,0}}+(-1)^flpp {{0,0,0,0},{0,0,0,0},{-1/2,1/2+flpp,1+2flpp,-4flpp},{-1/4,flpp/2,flpp,1-2flpp}});
    ,100];
    Join[{{tau,0,0,0}},M]]];



tauB[0]=1/2;
tauB[m_?NumericQ]:=tauB[m]=If[Im[m]==0&&Re[#]<1/2,(3#-1)/(4#-1),#]&[ToFundDomainApprox[InverseJ4[-4(Exp[I Pi m]+Exp[-I Pi m])]][[1]]];

ModularJ4[Rational[p_,q_]] := If[Mod[q,4]==0,ComplexInfinity,If[Mod[q,2]==0,-8,8]];
ModularJ4[_Integer] := 8;
ModularJ4[_DirectedInfinity]=Infinity;
ModularJ4[tau_]:=8+(DedekindEta[tau]/DedekindEta[4tau])^8;
ModularJ4[tau_] := With[{taurat=Rationalize[tau]},ModularJ4[taurat]/;Head[taurat]==Rational||Head[taurat]==Integer];



InversetauB[tau_]:=If[Im[#]<=0,-#,#]&[1/Pi ArcCos[-1/8 ModularJ4[tau]]];

InverseJ[j_]:=Module[{so,tau0,a},so=Solve[4a(1-a)==1728/j,a];
tau0=N[I Hypergeometric2F1[1/6,5/6,1,1-a/.Last[so ]]/Hypergeometric2F1[1/6,5/6,1,a/.Last[so]]];
Rationalize[Re[tau0]]+I Sqrt[Rationalize[Im[tau0]^2]
]];

InverseJ4[-8]=-1/2;
InverseJ4[-8.]=-.5;
InverseJ4[-8.+0.I]=-.5+0.I;

InverseJ4[j4_]:=Module[{tau0,Li,LiJ4},tau0=InverseJ[(j4^2+240j4+2112)^3/(j4-8)^4/(j4+8)];
Li={tau0,-1/tau0,-1/(tau0+1),-1/(tau0+2),-1/(tau0+3),-1/(2+-1/tau0)};
MinimalBy[Li,N[Abs[ModularJ4[#]-j4]]&][[1]]
];




(* Pi stability *)
ruleSeriesReducedEta[nmax_]:={DedekindEta[(p_Integer:1)tau/;p>0]:>q^(p/24) ReducedDedekindEta[p tau], ReducedDedekindEta[(p_.) tau]:>Sum[q^(p(6k^2+k))-q^(p(1+5k+6k^2))-q^(p(2+7k+6k^2))+q^(p(5+11k+6k^2)),{k,0,Floor[Sqrt[nmax/(6.p)]]}]+O[q]^(nmax+1)};
FourierCn[m_,n_Integer?NonNegative]:=(FourierCCompute[m,n];FourierCn[m,n]);
FourierCCompute[m_,nmax_Integer?NonNegative]:=With[{c=Cos[Pi m],qB2m=Exp[2.Pi Im[tauB[m]]]}, With[{coefs=CoefficientList[(DedekindEta[tau]^4 DedekindEta[2tau]^6)/DedekindEta[4tau]^4 (1+(8-8c)/(ModularJ4[tau]+8c))^(1/2)//.ruleSeriesReducedEta[nmax],q]}, Do[FourierCn[m,n]=coefs[[n+1]]/qB2m^n,{n,1,nmax}]]];
FourierCpn[m_,n_Integer?NonNegative]:=(FourierCpCompute[m,n];FourierCpn[m,n]);
FourierCpCompute[m_,nmax_Integer?NonNegative]:=With[{c=Cos[Pi m/2],qB2m=Exp[2.Pi Im[-.25/(tauB[m]-Boole[Re[tauB[m]]>=1/2])]]}, With[{coefs=CoefficientList[(DedekindEta[4tau]^4 DedekindEta[2tau]^6)/DedekindEta[tau]^4 1/c (1+(16(c^2-1))/(c^2 (ModularJ4[tau]-8)+16))^(1/2)//.ruleSeriesReducedEta[nmax],q]}, Do[FourierCpn[m,n]=coefs[[n+1]]/qB2m^n,{n,1,nmax}]]];
FourierCppn[m_,n_Integer?NonNegative]:=(FourierCppCompute[m,n];FourierCppn[m,n]);
FourierCppCompute[m_,nmax_Integer?NonNegative]:=With[{s=Sin[Pi m/2],qB2m=Exp[2.Pi Im[tauB[m]/(1-2tauB[m])]]}, With[{coefs=CoefficientList[q/(I s) ((DedekindEta[4tau]^8 DedekindEta[tau]^16)/(DedekindEta[2tau]^12 q) 16/(q(ModularJ4[tau]+8-16/s^2)))^(1/2)//.ruleSeriesReducedEta[nmax],q]}, Do[FourierCppn[m,n]=coefs[[n+1]]/qB2m^n,{n,1,nmax}]]];

(* old definitions, kept for reference: *)
FourierC[m_]:=FourierC[m]={1,-2 E^(-I m \[Pi])-2 E^(I m \[Pi]),-16+6 E^(-2 I m \[Pi])+6 E^(2 I m \[Pi]),36 E^(-I m \[Pi])+36 E^(I m \[Pi])-20 E^(-3 I m \[Pi])-20 E^(3 I m \[Pi]),112-128 E^(-2 I m \[Pi])-128 E^(2 I m \[Pi])+70 E^(-4 I m \[Pi])+70 E^(4 I m \[Pi]),-300 E^(-I m \[Pi])-300 E^(I m \[Pi])+500 E^(-3 I m \[Pi])+500 E^(3 I m \[Pi])-252 E^(-5 I m \[Pi])-252 E^(5 I m \[Pi]),-448+1332 E^(-2 I m \[Pi])+1332 E^(2 I m \[Pi])-2016 E^(-4 I m \[Pi])-2016 E^(4 I m \[Pi])+924 E^(-6 I m \[Pi])+924 E^(6 I m \[Pi]),1568 E^(-I m \[Pi])+1568 E^(I m \[Pi])-6272 E^(-3 I m \[Pi])-6272 E^(3 I m \[Pi])+8232 E^(-5 I m \[Pi])+8232 E^(5 I m \[Pi])-3432 E^(-7 I m \[Pi])-3432 E^(7 I m \[Pi]),1136-9216 E^(-2 I m \[Pi])-9216 E^(2 I m \[Pi])+29568 E^(-4 I m \[Pi])+29568 E^(4 I m \[Pi])-33792 E^(-6 I m \[Pi])-33792 E^(6 I m \[Pi])+12870 E^(-8 I m \[Pi])+12870 E^(8 I m \[Pi]),-5994 E^(-I m \[Pi])-5994 E^(I m \[Pi])+53496 E^(-3 I m \[Pi])+53496 E^(3 I m \[Pi])-138024 E^(-5 I m \[Pi])-138024 E^(5 I m \[Pi])+138996 E^(-7 I m \[Pi])+138996 E^(7 I m \[Pi])-48620 E^(-9 I m \[Pi])-48620 E^(9 I m \[Pi]),-2016+48900 E^(-2 I m \[Pi])+48900 E^(2 I m \[Pi])-297600 E^(-4 I m \[Pi])-297600 E^(4 I m \[Pi])+636900 E^(-6 I m \[Pi])+636900 E^(6 I m \[Pi])-572000 E^(-8 I m \[Pi])-572000 E^(8 I m \[Pi])+184756 E^(-10 I m \[Pi])+184756 E^(10 I m \[Pi]),18876 E^(-I m \[Pi])+18876 E^(I m \[Pi])-353804 E^(-3 I m \[Pi])-353804 E^(3 I m \[Pi])+1594296 E^(-5 I m \[Pi])+1594296 E^(5 I m \[Pi])-2906904 E^(-7 I m \[Pi])-2906904 E^(7 I m \[Pi])+2353208 E^(-9 I m \[Pi])+2353208 E^(9 I m \[Pi])-705432 E^(-11 I m \[Pi])-705432 E^(11 I m \[Pi]),3136-216576 E^(-2 I m \[Pi])-216576 E^(2 I m \[Pi])+2329524 E^(-4 I m \[Pi])+2329524 E^(4 I m \[Pi])-8281856 E^(-6 I m \[Pi])-8281856 E^(6 I m \[Pi])+13137696 E^(-8 I m \[Pi])+13137696 E^(8 I m \[Pi])-9674496 E^(-10 I m \[Pi])-9674496 E^(10 I m \[Pi])+2704156 E^(-12 I m \[Pi])+2704156 E^(12 I m \[Pi]),-53404 E^(-I m \[Pi])-53404 E^(I m \[Pi])+1951612 E^(-3 I m \[Pi])+1951612 E^(3 I m \[Pi])-14329172 E^(-5 I m \[Pi])-14329172 E^(5 I m \[Pi])+41953912 E^(-7 I m \[Pi])+41953912 E^(7 I m \[Pi])-58862024 E^(-9 I m \[Pi])-58862024 E^(9 I m \[Pi])+39739336 E^(-11 I m \[Pi])+39739336 E^(11 I m \[Pi])-10400600 E^(-13 I m \[Pi])-10400600 E^(13 I m \[Pi]),-5504+842016 E^(-2 I m \[Pi])+842016 E^(2 I m \[Pi])-15184512 E^(-4 I m \[Pi])-15184512 E^(4 I m \[Pi])+83756288 E^(-6 I m \[Pi])+83756288 E^(6 I m \[Pi])-208161408 E^(-8 I m \[Pi])-208161408 E^(8 I m \[Pi])+261715272 E^(-10 I m \[Pi])+261715272 E^(10 I m \[Pi])-163081408 E^(-12 I m \[Pi])-163081408 E^(12 I m \[Pi])+40116600 E^(-14 I m \[Pi])+40116600 E^(14 I m \[Pi]),140400 E^(-I m \[Pi])+140400 E^(I m \[Pi])-9382800 E^(-3 I m \[Pi])-9382800 E^(3 I m \[Pi])+107068536 E^(-5 I m \[Pi])+107068536 E^(5 I m \[Pi])-470503800 E^(-7 I m \[Pi])-470503800 E^(7 I m \[Pi])+1015008800 E^(-9 I m \[Pi])+1015008800 E^(9 I m \[Pi])-1155823200 E^(-11 I m \[Pi])-1155823200 E^(11 I m \[Pi])+668610000 E^(-13 I m \[Pi])+668610000 E^(13 I m \[Pi])-155117520 E^(-15 I m \[Pi])-155117520 E^(15 I m \[Pi]),9328-2957312 E^(-2 I m \[Pi])-2957312 E^(2 I m \[Pi])+85949440 E^(-4 I m \[Pi])+85949440 E^(4 I m \[Pi])-703078400 E^(-6 I m \[Pi])-703078400 E^(6 I m \[Pi])+2560168832 E^(-8 I m \[Pi])+2560168832 E^(8 I m \[Pi])-4876607488 E^(-10 I m \[Pi])-4876607488 E^(10 I m \[Pi])+5074066432 E^(-12 I m \[Pi])+5074066432 E^(12 I m \[Pi])-2738626560 E^(-14 I m \[Pi])-2738626560 E^(14 I m \[Pi])+601080390 E^(-16 I m \[Pi])+601080390 E^(16 I m \[Pi]),-343332 E^(-I m \[Pi])-343332 E^(I m \[Pi])+40443816 E^(-3 I m \[Pi])+40443816 E^(3 I m \[Pi])-692913336 E^(-5 I m \[Pi])-692913336 E^(5 I m \[Pi])+4373102916 E^(-7 I m \[Pi])+4373102916 E^(7 I m \[Pi])-13570708252 E^(-9 I m \[Pi])-13570708252 E^(9 I m \[Pi])+23133626928 E^(-11 I m \[Pi])+23133626928 E^(11 I m \[Pi])-22156843920 E^(-13 I m \[Pi])-22156843920 E^(13 I m \[Pi])+11207240820 E^(-15 I m \[Pi])+11207240820 E^(15 I m \[Pi])-2333606220 E^(-17 I m \[Pi])-2333606220 E^(17 I m \[Pi]),-12112+9556542 E^(-2 I m \[Pi])+9556542 E^(2 I m \[Pi])-434377728 E^(-4 I m \[Pi])-434377728 E^(4 I m \[Pi])+5100935832 E^(-6 I m \[Pi])+5100935832 E^(6 I m \[Pi])-26058459456 E^(-8 I m \[Pi])-26058459456 E^(8 I m \[Pi])+70374074616 E^(-10 I m \[Pi])+70374074616 E^(10 I m \[Pi])-108535193856 E^(-12 I m \[Pi])-108535193856 E^(12 I m \[Pi])+96291874980 E^(-14 I m \[Pi])+96291874980 E^(14 I m \[Pi])-45823540320 E^(-16 I m \[Pi])-45823540320 E^(16 I m \[Pi])+9075135300 E^(-18 I m \[Pi])+9075135300 E^(18 I m \[Pi]),786980 E^(-I m \[Pi])+786980 E^(I m \[Pi])-159390164 E^(-3 I m \[Pi])-159390164 E^(3 I m \[Pi])+3992189256 E^(-5 I m \[Pi])+3992189256 E^(5 I m \[Pi])-35066249064 E^(-7 I m \[Pi])-35066249064 E^(7 I m \[Pi])+149954943816 E^(-9 I m \[Pi])+149954943816 E^(9 I m \[Pi])-358196020584 E^(-11 I m \[Pi])-358196020584 E^(11 I m \[Pi])+504298646760 E^(-13 I m \[Pi])+504298646760 E^(13 I m \[Pi])-416686719240 E^(-15 I m \[Pi])-416686719240 E^(15 I m \[Pi])+187207076760 E^(-17 I m \[Pi])+187207076760 E^(17 I m \[Pi])-35345263800 E^(-19 I m \[Pi])-35345263800 E^(19 I m \[Pi]),14112-28819200 E^(-2 I m \[Pi])-28819200 E^(2 I m \[Pi])+1998932100 E^(-4 I m \[Pi])+1998932100 E^(4 I m \[Pi])-32871142400 E^(-6 I m \[Pi])-32871142400 E^(6 I m \[Pi])+228458121600 E^(-8 I m \[Pi])+228458121600 E^(8 I m \[Pi])-838270768128 E^(-10 I m \[Pi])-838270768128 E^(10 I m \[Pi])+1794095328100 E^(-12 I m \[Pi])+1794095328100 E^(12 I m \[Pi])-2323175904000 E^(-14 I m \[Pi])-2323175904000 E^(14 I m \[Pi])+1796169636000 E^(-16 I m \[Pi])+1796169636000 E^(16 I m \[Pi])-764221920000 E^(-18 I m \[Pi])-764221920000 E^(18 I m \[Pi])+137846528820 E^(-20 I m \[Pi])+137846528820 E^(20 I m \[Pi]),-1721664 E^(-I m \[Pi])-1721664 E^(I m \[Pi])+582696632 E^(-3 I m \[Pi])+582696632 E^(3 I m \[Pi])-20881318248 E^(-5 I m \[Pi])-20881318248 E^(5 I m \[Pi])+249274348632 E^(-7 I m \[Pi])+249274348632 E^(7 I m \[Pi])-1425187330280 E^(-9 I m \[Pi])-1425187330280 E^(9 I m \[Pi])+4572579836160 E^(-11 I m \[Pi])+4572579836160 E^(11 I m \[Pi])-8861035155840 E^(-13 I m \[Pi])-8861035155840 E^(13 I m \[Pi])+10620867699960 E^(-15 I m \[Pi])+10620867699960 E^(15 I m \[Pi])-7715393448840 E^(-17 I m \[Pi])-7715393448840 E^(17 I m \[Pi])+3117452267160 E^(-19 I m \[Pi])+3117452267160 E^(19 I m \[Pi])-538257874440 E^(-21 I m \[Pi])-538257874440 E^(21 I m \[Pi]),-21312+82010412 E^(-2 I m \[Pi])+82010412 E^(2 I m \[Pi])-8498486304 E^(-4 I m \[Pi])-8498486304 E^(4 I m \[Pi])+191859579076 E^(-6 I m \[Pi])+191859579076 E^(6 I m \[Pi])-1772707940608 E^(-8 I m \[Pi])-1772707940608 E^(8 I m \[Pi])+8577116722328 E^(-10 I m \[Pi])+8577116722328 E^(10 I m \[Pi])-24423289001408 E^(-12 I m \[Pi])-24423289001408 E^(12 I m \[Pi])+43228413250440 E^(-14 I m \[Pi])+43228413250440 E^(14 I m \[Pi])-48224149944960 E^(-16 I m \[Pi])-48224149944960 E^(16 I m \[Pi])+33035211991320 E^(-18 I m \[Pi])+33035211991320 E^(18 I m \[Pi])-12708137133120 E^(-20 I m \[Pi])-12708137133120 E^(20 I m \[Pi])+2104098963720 E^(-22 I m \[Pi])+2104098963720 E^(22 I m \[Pi]),3631056 E^(-I m \[Pi])+3631056 E^(I m \[Pi])-1997901808 E^(-3 I m \[Pi])-1997901808 E^(3 I m \[Pi])+100609870968 E^(-5 I m \[Pi])+100609870968 E^(5 I m \[Pi])-1601903745528 E^(-7 I m \[Pi])-1601903745528 E^(7 I m \[Pi])+11970840291184 E^(-9 I m \[Pi])+11970840291184 E^(9 I m \[Pi])-50080411878576 E^(-11 I m \[Pi])-50080411878576 E^(11 I m \[Pi])+128090240705520 E^(-13 I m \[Pi])+128090240705520 E^(13 I m \[Pi])-208595503068720 E^(-15 I m \[Pi])-208595503068720 E^(15 I m \[Pi])+217614593560320 E^(-17 I m \[Pi])+217614593560320 E^(17 I m \[Pi])-141033661750080 E^(-19 I m \[Pi])-141033661750080 E^(19 I m \[Pi])+51770621014320 E^(-21 I m \[Pi])+51770621014320 E^(21 I m \[Pi])-8233430727600 E^(-23 I m \[Pi])-8233430727600 E^(23 I m \[Pi]),31808-222031872 E^(-2 I m \[Pi])-222031872 E^(2 I m \[Pi])+33753752064 E^(-4 I m \[Pi])+33753752064 E^(4 I m \[Pi])-1029180911616 E^(-6 I m \[Pi])-1029180911616 E^(6 I m \[Pi])+12414452113332 E^(-8 I m \[Pi])+12414452113332 E^(8 I m \[Pi])-77459531102208 E^(-10 I m \[Pi])-77459531102208 E^(10 I m \[Pi])+284938362239232 E^(-12 I m \[Pi])+284938362239232 E^(12 I m \[Pi])-661106293585920 E^(-14 I m \[Pi])-661106293585920 E^(14 I m \[Pi])+996776069525280 E^(-16 I m \[Pi])+996776069525280 E^(16 I m \[Pi])-976520149186560 E^(-18 I m \[Pi])-976520149186560 E^(18 I m \[Pi])+600480962115840 E^(-20 I m \[Pi])+600480962115840 E^(20 I m \[Pi])-210775826626560 E^(-22 I m \[Pi])-210775826626560 E^(22 I m \[Pi])+32247603683100 E^(-24 I m \[Pi])+32247603683100 E^(24 I m \[Pi]),-7388750 E^(-I m \[Pi])-7388750 E^(I m \[Pi])+6479820000 E^(-3 I m \[Pi])+6479820000 E^(3 I m \[Pi])-451563312800 E^(-5 I m \[Pi])-451563312800 E^(5 I m \[Pi])+9443173120000 E^(-7 I m \[Pi])+9443173120000 E^(7 I m \[Pi])-90618384320000 E^(-9 I m \[Pi])-90618384320000 E^(9 I m \[Pi])+483548212480000 E^(-11 I m \[Pi])+483548212480000 E^(11 I m \[Pi])-1585236120800000 E^(-13 I m \[Pi])-1585236120800000 E^(13 I m \[Pi])+3364115718250500 E^(-15 I m \[Pi])+3364115718250500 E^(15 I m \[Pi])-4721450020237500 E^(-17 I m \[Pi])-4721450020237500 E^(17 I m \[Pi])+4359746483850000 E^(-19 I m \[Pi])+4359746483850000 E^(19 I m \[Pi])-2550342399150000 E^(-21 I m \[Pi])-2550342399150000 E^(21 I m \[Pi])+857649034125000 E^(-23 I m \[Pi])+857649034125000 E^(23 I m \[Pi])-126410606437752 E^(-25 I m \[Pi])-126410606437752 E^(25 I m \[Pi]),-35168+575463252 E^(-2 I m \[Pi])+575463252 E^(2 I m \[Pi])-126335574144 E^(-4 I m \[Pi])-126335574144 E^(4 I m \[Pi])+5131307204556 E^(-6 I m \[Pi])+5131307204556 E^(6 I m \[Pi])-79624685496096 E^(-8 I m \[Pi])-79624685496096 E^(8 I m \[Pi])+629645371993404 E^(-10 I m \[Pi])+629645371993404 E^(10 I m \[Pi])-2927497388964096 E^(-12 I m \[Pi])-2927497388964096 E^(12 I m \[Pi])+8647879881415320 E^(-14 I m \[Pi])+8647879881415320 E^(14 I m \[Pi])-16903714660073280 E^(-16 I m \[Pi])-16903714660073280 E^(16 I m \[Pi])+22187267535402840 E^(-18 I m \[Pi])+22187267535402840 E^(18 I m \[Pi])-19373906045295360 E^(-20 I m \[Pi])-19373906045295360 E^(20 I m \[Pi])+10806926725356840 E^(-22 I m \[Pi])+10806926725356840 E^(22 I m \[Pi])-3487900814364096 E^(-24 I m \[Pi])-3487900814364096 E^(24 I m \[Pi])+495918532948104 E^(-26 I m \[Pi])+495918532948104 E^(26 I m \[Pi]),14550840 E^(-I m \[Pi])+14550840 E^(I m \[Pi])-20015031960 E^(-3 I m \[Pi])-20015031960 E^(3 I m \[Pi])+1904702069160 E^(-5 I m \[Pi])+1904702069160 E^(5 I m \[Pi])-51645974352840 E^(-7 I m \[Pi])-51645974352840 E^(7 I m \[Pi])+627387270118632 E^(-9 I m \[Pi])+627387270118632 E^(9 I m \[Pi])-4197622654224648 E^(-11 I m \[Pi])-4197622654224648 E^(11 I m \[Pi])+17260059513187800 E^(-13 I m \[Pi])+17260059513187800 E^(13 I m \[Pi])-46365141858636600 E^(-15 I m \[Pi])-46365141858636600 E^(15 I m \[Pi])+83978297123279400 E^(-17 I m \[Pi])+83978297123279400 E^(17 I m \[Pi])-103513628292821640 E^(-19 I m \[Pi])-103513628292821640 E^(19 I m \[Pi])+85726442440744560 E^(-21 I m \[Pi])+85726442440744560 E^(21 I m \[Pi])-45696528549867312 E^(-23 I m \[Pi])-45696528549867312 E^(23 I m \[Pi])+14177435706634032 E^(-25 I m \[Pi])+14177435706634032 E^(25 I m \[Pi])-1946939425648112 E^(-27 I m \[Pi])-1946939425648112 E^(27 I m \[Pi]),38528-1435184128 E^(-2 I m \[Pi])-1435184128 E^(2 I m \[Pi])+448736049440 E^(-4 I m \[Pi])+448736049440 E^(4 I m \[Pi])-23991839298560 E^(-6 I m \[Pi])-23991839298560 E^(6 I m \[Pi])+473081359196544 E^(-8 I m \[Pi])+473081359196544 E^(8 I m \[Pi])-4675331608679424 E^(-10 I m \[Pi])-4675331608679424 E^(10 I m \[Pi])+27013669185927424 E^(-12 I m \[Pi])+27013669185927424 E^(12 I m \[Pi])-99431247904035840 E^(-14 I m \[Pi])-99431247904035840 E^(14 I m \[Pi])+244773857899186560 E^(-16 I m \[Pi])+244773857899186560 E^(16 I m \[Pi])-412956158084305920 E^(-18 I m \[Pi])-412956158084305920 E^(18 I m \[Pi])+479765994137768520 E^(-20 I m \[Pi])+479765994137768520 E^(20 I m \[Pi])-377833233880327680 E^(-22 I m \[Pi])-377833233880327680 E^(22 I m \[Pi])+192844242067192128 E^(-24 I m \[Pi])+192844242067192128 E^(24 I m \[Pi])-57600019234268672 E^(-26 I m \[Pi])-57600019234268672 E^(26 I m \[Pi])+7648690600760440 E^(-28 I m \[Pi])+7648690600760440 E^(28 I m \[Pi]),-27904380 E^(-I m \[Pi])-27904380 E^(I m \[Pi])+59204317684 E^(-3 I m \[Pi])+59204317684 E^(3 I m \[Pi])-7604280148476 E^(-5 I m \[Pi])-7604280148476 E^(5 I m \[Pi])+264419303687184 E^(-7 I m \[Pi])+264419303687184 E^(7 I m \[Pi])-4018309042087792 E^(-9 I m \[Pi])-4018309042087792 E^(9 I m \[Pi])+33249499873467528 E^(-11 I m \[Pi])+33249499873467528 E^(11 I m \[Pi])-168624301769396760 E^(-13 I m \[Pi])-168624301769396760 E^(13 I m \[Pi])+561203482312881240 E^(-15 I m \[Pi])+561203482312881240 E^(15 I m \[Pi])-1274452047651038760 E^(-17 I m \[Pi])-1274452047651038760 E^(17 I m \[Pi])+2011886133468644520 E^(-19 I m \[Pi])+2011886133468644520 E^(19 I m \[Pi])-2210232479074459320 E^(-21 I m \[Pi])-2210232479074459320 E^(21 I m \[Pi])+1659215071059537456 E^(-23 I m \[Pi])+1659215071059537456 E^(23 I m \[Pi])-812327522159398224 E^(-25 I m \[Pi])-812327522159398224 E^(25 I m \[Pi])+233910865281437456 E^(-27 I m \[Pi])+233910865281437456 E^(27 I m \[Pi])-30067266499541040 E^(-29 I m \[Pi])-30067266499541040 E^(29 I m \[Pi]),-56448+3459034800 E^(-2 I m \[Pi])+3459034800 E^(2 I m \[Pi])-1521256550400 E^(-4 I m \[Pi])-1521256550400 E^(4 I m \[Pi])+105956973447600 E^(-6 I m \[Pi])+105956973447600 E^(6 I m \[Pi])-2627358387062400 E^(-8 I m \[Pi])-2627358387062400 E^(8 I m \[Pi])+32077886522351832 E^(-10 I m \[Pi])+32077886522351832 E^(10 I m \[Pi])-227225343250478400 E^(-12 I m \[Pi])-227225343250478400 E^(12 I m \[Pi])+1024911836619849000 E^(-14 I m \[Pi])+1024911836619849000 E^(14 I m \[Pi])-3110388188013216000 E^(-16 I m \[Pi])-3110388188013216000 E^(16 I m \[Pi])+6553237194852612000 E^(-18 I m \[Pi])+6553237194852612000 E^(18 I m \[Pi])-9718980134526614016 E^(-20 I m \[Pi])-9718980134526614016 E^(20 I m \[Pi])+10125872113522356000 E^(-22 I m \[Pi])+10125872113522356000 E^(22 I m \[Pi])-7261688487062278400 E^(-24 I m \[Pi])-7261688487062278400 E^(24 I m \[Pi])+3415934085283011600 E^(-26 I m \[Pi])+3415934085283011600 E^(26 I m \[Pi])-949492626301296000 E^(-28 I m \[Pi])-949492626301296000 E^(28 I m \[Pi])+118264581564861424 E^(-30 I m \[Pi])+118264581564861424 E^(30 I m \[Pi]),52293776 E^(-I m \[Pi])+52293776 E^(I m \[Pi])-168475462416 E^(-3 I m \[Pi])-168475462416 E^(3 I m \[Pi])+28903490729824 E^(-5 I m \[Pi])+28903490729824 E^(5 I m \[Pi])-1276587551424416 E^(-7 I m \[Pi])-1276587551424416 E^(7 I m \[Pi])+24026489663674768 E^(-9 I m \[Pi])+24026489663674768 E^(9 I m \[Pi])-243100896444018512 E^(-11 I m \[Pi])-243100896444018512 E^(11 I m \[Pi])+1500340664454749200 E^(-13 I m \[Pi])+1500340664454749200 E^(13 I m \[Pi])-6084833511181141200 E^(-15 I m \[Pi])-6084833511181141200 E^(15 I m \[Pi])+16960075301885815200 E^(-17 I m \[Pi])+16960075301885815200 E^(17 I m \[Pi])-33317072386801016928 E^(-19 I m \[Pi])-33317072386801016928 E^(19 I m \[Pi])+46587074684615927952 E^(-21 I m \[Pi])+46587074684615927952 E^(21 I m \[Pi])-46152788726235038864 E^(-23 I m \[Pi])-46152788726235038864 E^(23 I m \[Pi])+31681663949319272896 E^(-25 I m \[Pi])+31681663949319272896 E^(25 I m \[Pi])-14341328444346956864 E^(-27 I m \[Pi])-14341328444346956864 E^(27 I m \[Pi])+3852619080807858592 E^(-29 I m \[Pi])+3852619080807858592 E^(29 I m \[Pi])-465428353255261088 E^(-31 I m \[Pi])-465428353255261088 E^(31 I m \[Pi]),74864-8085110784 E^(-2 I m \[Pi])-8085110784 E^(2 I m \[Pi])+4945643298816 E^(-4 I m \[Pi])+4945643298816 E^(4 I m \[Pi])-444635092680704 E^(-6 I m \[Pi])-444635092680704 E^(6 I m \[Pi])+13740026274020352 E^(-8 I m \[Pi])+13740026274020352 E^(8 I m \[Pi])-205233736325332992 E^(-10 I m \[Pi])-205233736325332992 E^(10 I m \[Pi])+1762529905005617152 E^(-12 I m \[Pi])+1762529905005617152 E^(12 I m \[Pi])-9613459157191557120 E^(-14 I m \[Pi])-9613459157191557120 E^(14 I m \[Pi])+35378366988324574080 E^(-16 I m \[Pi])+35378366988324574080 E^(16 I m \[Pi])-91129094507654086656 E^(-18 I m \[Pi])-91129094507654086656 E^(18 I m \[Pi])+167645653727708307456 E^(-20 I m \[Pi])+167645653727708307456 E^(20 I m \[Pi])-221722195714286616576 E^(-22 I m \[Pi])-221722195714286616576 E^(22 I m \[Pi])+209362464757308529664 E^(-24 I m \[Pi])+209362464757308529664 E^(24 I m \[Pi])-137818113137062772736 E^(-26 I m \[Pi])-137818113137062772736 E^(26 I m \[Pi])+60119341117061259264 E^(-28 I m \[Pi])+60119341117061259264 E^(28 I m \[Pi])-15626184712570077184 E^(-30 I m \[Pi])-15626184712570077184 E^(30 I m \[Pi])+1832624140942590534 E^(-32 I m \[Pi])+1832624140942590534 E^(32 I m \[Pi]),-95866848 E^(-I m \[Pi])-95866848 E^(I m \[Pi])+462980767656 E^(-3 I m \[Pi])+462980767656 E^(3 I m \[Pi])-105103942046136 E^(-5 I m \[Pi])-105103942046136 E^(5 I m \[Pi])+5846793588080676 E^(-7 I m \[Pi])+5846793588080676 E^(7 I m \[Pi])-135110681063399228 E^(-9 I m \[Pi])-135110681063399228 E^(9 I m \[Pi])+1655771588758093392 E^(-11 I m \[Pi])+1655771588758093392 E^(11 I m \[Pi])-12299558790832481520 E^(-13 I m \[Pi])-12299558790832481520 E^(13 I m \[Pi])+59989037184362511720 E^(-15 I m \[Pi])+59989037184362511720 E^(15 I m \[Pi])-201883904425328641560 E^(-17 I m \[Pi])-201883904425328641560 E^(17 I m \[Pi])+483169045488611297904 E^(-19 I m \[Pi])+483169045488611297904 E^(19 I m \[Pi])-835625498440097796816 E^(-21 I m \[Pi])-835625498440097796816 E^(21 I m \[Pi])+1048314083198916339432 E^(-23 I m \[Pi])+1048314083198916339432 E^(23 I m \[Pi])-945543900921061021272 E^(-25 I m \[Pi])-945543900921061021272 E^(25 I m \[Pi])+597882635881880673632 E^(-27 I m \[Pi])+597882635881880673632 E^(27 I m \[Pi])-251665348389532662816 E^(-29 I m \[Pi])-251665348389532662816 E^(29 I m \[Pi])+63356434586872415604 E^(-31 I m \[Pi])+63356434586872415604 E^(31 I m \[Pi])-7219428434016265740 E^(-33 I m \[Pi])-7219428434016265740 E^(33 I m \[Pi]),-78624+18380875116 E^(-2 I m \[Pi])+18380875116 E^(2 I m \[Pi])-15480558500352 E^(-4 I m \[Pi])-15480558500352 E^(4 I m \[Pi])+1781721793466248 E^(-6 I m \[Pi])+1781721793466248 E^(6 I m \[Pi])-68073142578260928 E^(-8 I m \[Pi])-68073142578260928 E^(8 I m \[Pi])+1233603767516803752 E^(-10 I m \[Pi])+1233603767516803752 E^(10 I m \[Pi])-12724300887373852928 E^(-12 I m \[Pi])-12724300887373852928 E^(12 I m \[Pi])+83015111355045907380 E^(-14 I m \[Pi])+83015111355045907380 E^(14 I m \[Pi])-365632798986116122080 E^(-16 I m \[Pi])-365632798986116122080 E^(16 I m \[Pi])+1132774460270815292628 E^(-18 I m \[Pi])+1132774460270815292628 E^(18 I m \[Pi])-2530864240359071788032 E^(-20 I m \[Pi])-2530864240359071788032 E^(20 I m \[Pi])+4129118408681738902128 E^(-22 I m \[Pi])+4129118408681738902128 E^(22 I m \[Pi])-4926342736233601450112 E^(-24 I m \[Pi])-4926342736233601450112 E^(24 I m \[Pi])+4252829605904134432688 E^(-26 I m \[Pi])+4252829605904134432688 E^(26 I m \[Pi])-2587094633434532311552 E^(-28 I m \[Pi])-2587094633434532311552 E^(28 I m \[Pi])+1052090301320636935588 E^(-30 I m \[Pi])+1052090301320636935588 E^(30 I m \[Pi])-256789515991470867552 E^(-32 I m \[Pi])-256789515991470867552 E^(32 I m \[Pi])+28453041475240576740 E^(-34 I m \[Pi])+28453041475240576740 E^(34 I m \[Pi]),172225200 E^(-I m \[Pi])+172225200 E^(I m \[Pi])-1232633788400 E^(-3 I m \[Pi])-1232633788400 E^(3 I m \[Pi])+367155305777568 E^(-5 I m \[Pi])+367155305777568 E^(5 I m \[Pi])-25531640014106400 E^(-7 I m \[Pi])-25531640014106400 E^(7 I m \[Pi])+718951398657458400 E^(-9 I m \[Pi])+718951398657458400 E^(9 I m \[Pi])-10584742330502037600 E^(-11 I m \[Pi])-10584742330502037600 E^(11 I m \[Pi])+93766445999468193000 E^(-13 I m \[Pi])+93766445999468193000 E^(13 I m \[Pi])-544076923588981243272 E^(-15 I m \[Pi])-544076923588981243272 E^(15 I m \[Pi])+2182078027863476127000 E^(-17 I m \[Pi])+2182078027863476127000 E^(17 I m \[Pi])-6259645432780209199800 E^(-19 I m \[Pi])-6259645432780209199800 E^(19 I m \[Pi])+13110341118350814829200 E^(-21 I m \[Pi])+13110341118350814829200 E^(21 I m \[Pi])-20240542000762209646800 E^(-23 I m \[Pi])-20240542000762209646800 E^(23 I m \[Pi])+23019617206052379137232 E^(-25 I m \[Pi])+23019617206052379137232 E^(25 I m \[Pi])-19054866651372337046800 E^(-27 I m \[Pi])-19054866651372337046800 E^(27 I m \[Pi])+11167702937457797392200 E^(-29 I m \[Pi])+11167702937457797392200 E^(29 I m \[Pi])-4392744531774512461800 E^(-31 I m \[Pi])-4392744531774512461800 E^(31 I m \[Pi])+1040447039019991239000 E^(-33 I m \[Pi])+1040447039019991239000 E^(33 I m \[Pi])-112186277816662845432 E^(-35 I m \[Pi])-112186277816662845432 E^(35 I m \[Pi]),84784-40745794176 E^(-2 I m \[Pi])-40745794176 E^(2 I m \[Pi])+46813333282686 E^(-4 I m \[Pi])+46813333282686 E^(4 I m \[Pi])-6846372755841024 E^(-6 I m \[Pi])-6846372755841024 E^(6 I m \[Pi])+321139250710612992 E^(-8 I m \[Pi])+321139250710612992 E^(8 I m \[Pi])-7009145053489963008 E^(-10 I m \[Pi])-7009145053489963008 E^(10 I m \[Pi])+86143073947499319192 E^(-12 I m \[Pi])+86143073947499319192 E^(12 I m \[Pi])-666136492303597226496 E^(-14 I m \[Pi])-666136492303597226496 E^(14 I m \[Pi])+3473935549505492140224 E^(-16 I m \[Pi])+3473935549505492140224 E^(16 I m \[Pi])-12777731958521929452032 E^(-18 I m \[Pi])-12777731958521929452032 E^(18 I m \[Pi])+34112435177858413400568 E^(-20 I m \[Pi])+34112435177858413400568 E^(20 I m \[Pi])-67224315590143108319232 E^(-22 I m \[Pi])-67224315590143108319232 E^(22 I m \[Pi])+98483979360779467515648 E^(-24 I m \[Pi])+98483979360779467515648 E^(24 I m \[Pi])-106999393290507509041152 E^(-26 I m \[Pi])-106999393290507509041152 E^(26 I m \[Pi])+85069184920932046413348 E^(-28 I m \[Pi])+85069184920932046413348 E^(28 I m \[Pi])-48098645086327745962752 E^(-30 I m \[Pi])-48098645086327745962752 E^(30 I m \[Pi])+18319028514411286374048 E^(-32 I m \[Pi])+18319028514411286374048 E^(32 I m \[Pi])-4214301914504204280576 E^(-34 I m \[Pi])-4214301914504204280576 E^(34 I m \[Pi])+442512540276836779204 E^(-36 I m \[Pi])+442512540276836779204 E^(36 I m \[Pi]),-303956332 E^(-I m \[Pi])-303956332 E^(I m \[Pi])+3188230888060 E^(-3 I m \[Pi])+3188230888060 E^(3 I m \[Pi])-1236422051528020 E^(-5 I m \[Pi])-1236422051528020 E^(5 I m \[Pi])+106754207980557800 E^(-7 I m \[Pi])+106754207980557800 E^(7 I m \[Pi])-3638740371699959512 E^(-9 I m \[Pi])-3638740371699959512 E^(9 I m \[Pi])+63903262880178506168 E^(-11 I m \[Pi])+63903262880178506168 E^(11 I m \[Pi])-669809443488809010472 E^(-13 I m \[Pi])-669809443488809010472 E^(13 I m \[Pi])+4582190620094875901952 E^(-15 I m \[Pi])+4582190620094875901952 E^(15 I m \[Pi])-21668677101731926372352 E^(-17 I m \[Pi])-21668677101731926372352 E^(17 I m \[Pi])+73547913927843365204272 E^(-19 I m \[Pi])+73547913927843365204272 E^(19 I m \[Pi])-183546849608345769097488 E^(-21 I m \[Pi])-183546849608345769097488 E^(21 I m \[Pi])+341468867992288512487152 E^(-23 I m \[Pi])+341468867992288512487152 E^(23 I m \[Pi])-475903353665982554147472 E^(-25 I m \[Pi])-475903353665982554147472 E^(25 I m \[Pi])+494910436842213634248552 E^(-27 I m \[Pi])+494910436842213634248552 E^(27 I m \[Pi])-378506920250642133701688 E^(-29 I m \[Pi])-378506920250642133701688 E^(29 I m \[Pi])+206716977191797634247672 E^(-31 I m \[Pi])+206716977191797634247672 E^(31 I m \[Pi])-76309795353305919771912 E^(-33 I m \[Pi])-76309795353305919771912 E^(33 I m \[Pi])+17064779370112381710712 E^(-35 I m \[Pi])+17064779370112381710712 E^(35 I m \[Pi])-1746130564335626209832 E^(-37 I m \[Pi])-1746130564335626209832 E^(37 I m \[Pi]),-109760+88262832180 E^(-2 I m \[Pi])+88262832180 E^(2 I m \[Pi])-137163782089440 E^(-4 I m \[Pi])-137163782089440 E^(4 I m \[Pi])+25317739520861660 E^(-6 I m \[Pi])+25317739520861660 E^(6 I m \[Pi])-1448824956232270080 E^(-8 I m \[Pi])-1448824956232270080 E^(8 I m \[Pi])+37841411551211301480 E^(-10 I m \[Pi])+37841411551211301480 E^(10 I m \[Pi])-550312569148552866880 E^(-12 I m \[Pi])-550312569148552866880 E^(12 I m \[Pi])+5004960182043472595448 E^(-14 I m \[Pi])+5004960182043472595448 E^(14 I m \[Pi])-30629558478969516968832 E^(-16 I m \[Pi])-30629558478969516968832 E^(16 I m \[Pi])+132342668589667097805416 E^(-18 I m \[Pi])+132342668589667097805416 E^(18 I m \[Pi])-416761663742440114408896 E^(-20 I m \[Pi])-416761663742440114408896 E^(20 I m \[Pi])+976118294875687508455416 E^(-22 I m \[Pi])+976118294875687508455416 E^(22 I m \[Pi])-1719466513517759886305792 E^(-24 I m \[Pi])-1719466513517759886305792 E^(24 I m \[Pi])+2285011418183993045430408 E^(-26 I m \[Pi])+2285011418183993045430408 E^(26 I m \[Pi])-2278617593692374014305344 E^(-28 I m \[Pi])-2278617593692374014305344 E^(28 I m \[Pi])+1678798409642347968725464 E^(-30 I m \[Pi])+1678798409642347968725464 E^(30 I m \[Pi])-886636921164038406601344 E^(-32 I m \[Pi])-886636921164038406601344 E^(32 I m \[Pi])+317537568751769135879928 E^(-34 I m \[Pi])+317537568751769135879928 E^(34 I m \[Pi])-69079795476729979369792 E^(-36 I m \[Pi])-69079795476729979369792 E^(36 I m \[Pi])+6892620648693261354600 E^(-38 I m \[Pi])+6892620648693261354600 E^(38 I m \[Pi]),527896512 E^(-I m \[Pi])+527896512 E^(I m \[Pi])-8030564664928 E^(-3 I m \[Pi])-8030564664928 E^(3 I m \[Pi])+4026088424373912 E^(-5 I m \[Pi])+4026088424373912 E^(5 I m \[Pi])-428966557670311128 E^(-7 I m \[Pi])-428966557670311128 E^(7 I m \[Pi])+17593024933382028208 E^(-9 I m \[Pi])+17593024933382028208 E^(9 I m \[Pi])-366265887233332827312 E^(-11 I m \[Pi])-366265887233332827312 E^(11 I m \[Pi])+4511593513177210378656 E^(-13 I m \[Pi])+4511593513177210378656 E^(13 I m \[Pi])-36109974637249434987104 E^(-15 I m \[Pi])-36109974637249434987104 E^(15 I m \[Pi])+199560323821668996322656 E^(-17 I m \[Pi])+199560323821668996322656 E^(17 I m \[Pi])-793021642425472821113376 E^(-19 I m \[Pi])-793021642425472821113376 E^(19 I m \[Pi])+2328015568512074228039824 E^(-21 I m \[Pi])+2328015568512074228039824 E^(21 I m \[Pi])-5135413523357816715334032 E^(-23 I m \[Pi])-5135413523357816715334032 E^(23 I m \[Pi])+8588674488897322890387408 E^(-25 I m \[Pi])+8588674488897322890387408 E^(25 I m \[Pi])-10905864648319573375566032 E^(-27 I m \[Pi])-10905864648319573375566032 E^(27 I m \[Pi])+10445765915856969265093344 E^(-29 I m \[Pi])+10445765915856969265093344 E^(29 I m \[Pi])-7423807031338193294381856 E^(-31 I m \[Pi])-7423807031338193294381856 E^(31 I m \[Pi])+3795680157560072900512576 E^(-33 I m \[Pi])+3795680157560072900512576 E^(33 I m \[Pi])-1319989808913980264371584 E^(-35 I m \[Pi])-1319989808913980264371584 E^(35 I m \[Pi])+279564693510998680542576 E^(-37 I m \[Pi])+279564693510998680542576 E^(37 I m \[Pi])-27217014869199032015600 E^(-39 I m \[Pi])-27217014869199032015600 E^(39 I m \[Pi]),143136-187181414400 E^(-2 I m \[Pi])-187181414400 E^(2 I m \[Pi])+390386842425600 E^(-4 I m \[Pi])+390386842425600 E^(4 I m \[Pi])-90381572518246400 E^(-6 I m \[Pi])-90381572518246400 E^(6 I m \[Pi])+6274121276211770500 E^(-8 I m \[Pi])+6274121276211770500 E^(8 I m \[Pi])-194983070102168969216 E^(-10 I m \[Pi])-194983070102168969216 E^(10 I m \[Pi])+3334914167424804736000 E^(-12 I m \[Pi])+3334914167424804736000 E^(12 I m \[Pi])-35433375091367222681600 E^(-14 I m \[Pi])-35433375091367222681600 E^(14 I m \[Pi])+252540998919812132694400 E^(-16 I m \[Pi])+252540998919812132694400 E^(16 I m \[Pi])-1270515940382533610086400 E^(-18 I m \[Pi])-1270515940382533610086400 E^(18 I m \[Pi])+4670135798039231956077568 E^(-20 I m \[Pi])+4670135798039231956077568 E^(20 I m \[Pi])-12834287825418808838758400 E^(-22 I m \[Pi])-12834287825418808838758400 E^(22 I m \[Pi])+26749548511739260506105700 E^(-24 I m \[Pi])+26749548511739260506105700 E^(24 I m \[Pi])-42578612953138655689676800 E^(-26 I m \[Pi])-42578612953138655689676800 E^(26 I m \[Pi])+51760851005093544399545600 E^(-28 I m \[Pi])+51760851005093544399545600 E^(28 I m \[Pi])-47692076233469494447265792 E^(-30 I m \[Pi])-47692076233469494447265792 E^(30 I m \[Pi])+32736448984509109733437600 E^(-32 I m \[Pi])+32736448984509109733437600 E^(32 I m \[Pi])-16220023582543898058905600 E^(-34 I m \[Pi])-16220023582543898058905600 E^(34 I m \[Pi])+5481868009620079861638400 E^(-36 I m \[Pi])+5481868009620079861638400 E^(36 I m \[Pi])-1131096721836842888960000 E^(-38 I m \[Pi])-1131096721836842888960000 E^(38 I m \[Pi])+107507208733336176461620 E^(-40 I m \[Pi])+107507208733336176461620 E^(40 I m \[Pi]),-903039924 E^(-I m \[Pi])-903039924 E^(I m \[Pi])+19738967891168 E^(-3 I m \[Pi])+19738967891168 E^(3 I m \[Pi])-12709875270129312 E^(-5 I m \[Pi])-12709875270129312 E^(5 I m \[Pi])+1661754656776334568 E^(-7 I m \[Pi])+1661754656776334568 E^(7 I m \[Pi])-81563681598208880856 E^(-9 I m \[Pi])-81563681598208880856 E^(9 I m \[Pi])+2001848504153424982224 E^(-11 I m \[Pi])+2001848504153424982224 E^(11 I m \[Pi])-28805899709636159858544 E^(-13 I m \[Pi])-28805899709636159858544 E^(13 I m \[Pi])+267965372809366635375596 E^(-15 I m \[Pi])+267965372809366635375596 E^(15 I m \[Pi])-1717661830208546832629844 E^(-17 I m \[Pi])-1717661830208546832629844 E^(17 I m \[Pi])+7921428932800764489822720 E^(-19 I m \[Pi])+7921428932800764489822720 E^(19 I m \[Pi])-27069309275392905685212160 E^(-21 I m \[Pi])-27069309275392905685212160 E^(21 I m \[Pi])+69902317696789359250057224 E^(-23 I m \[Pi])+69902317696789359250057224 E^(23 I m \[Pi])-138050337274532836422707256 E^(-25 I m \[Pi])-138050337274532836422707256 E^(25 I m \[Pi])+209608135554675702606986224 E^(-27 I m \[Pi])+209608135554675702606986224 E^(27 I m \[Pi])-244378981692756247696996560 E^(-29 I m \[Pi])-244378981692756247696996560 E^(29 I m \[Pi])+216916370831878934281615620 E^(-31 I m \[Pi])+216916370831878934281615620 E^(31 I m \[Pi])-143972817341564036220180604 E^(-33 I m \[Pi])-143972817341564036220180604 E^(33 I m \[Pi])+69194558935949446178861136 E^(-35 I m \[Pi])+69194558935949446178861136 E^(35 I m \[Pi])-22745181563290004772488304 E^(-37 I m \[Pi])-22745181563290004772488304 E^(37 I m \[Pi])+4575180199512357281822360 E^(-39 I m \[Pi])+4575180199512357281822360 E^(39 I m \[Pi])-424784580848791721628840 E^(-41 I m \[Pi])-424784580848791721628840 E^(41 I m \[Pi]),-154112+389261541312 E^(-2 I m \[Pi])+389261541312 E^(2 I m \[Pi])-1081674478596096 E^(-4 I m \[Pi])-1081674478596096 E^(4 I m \[Pi])+312319193163489624 E^(-6 I m \[Pi])+312319193163489624 E^(6 I m \[Pi])-26163669835946164032 E^(-8 I m \[Pi])-26163669835946164032 E^(8 I m \[Pi])+962497757485589750712 E^(-10 I m \[Pi])+962497757485589750712 E^(10 I m \[Pi])-19257119726747355032832 E^(-12 I m \[Pi])-19257119726747355032832 E^(12 I m \[Pi])+237634789776470840022264 E^(-14 I m \[Pi])+237634789776470840022264 E^(14 I m \[Pi])-1959583085400396925094976 E^(-16 I m \[Pi])-1959583085400396925094976 E^(16 I m \[Pi])+11393255174956471511194808 E^(-18 I m \[Pi])+11393255174956471511194808 E^(18 I m \[Pi])-48457214009697432139699968 E^(-20 I m \[Pi])-48457214009697432139699968 E^(20 I m \[Pi])+154628975304373048510929408 E^(-22 I m \[Pi])+154628975304373048510929408 E^(22 I m \[Pi])-376478356692538620338439936 E^(-24 I m \[Pi])-376478356692538620338439936 E^(24 I m \[Pi])+706345606634352305055329664 E^(-26 I m \[Pi])+706345606634352305055329664 E^(26 I m \[Pi])-1025115552241079082321629952 E^(-28 I m \[Pi])-1025115552241079082321629952 E^(28 I m \[Pi])+1148114154472654609298976792 E^(-30 I m \[Pi])+1148114154472654609298976792 E^(30 I m \[Pi])-983043618914729935228421952 E^(-32 I m \[Pi])-983043618914729935228421952 E^(32 I m \[Pi])+631590104720975516725825944 E^(-34 I m \[Pi])+631590104720975516725825944 E^(34 I m \[Pi])-294705636350227559937308416 E^(-36 I m \[Pi])-294705636350227559937308416 E^(36 I m \[Pi])+94291368595076831942990520 E^(-38 I m \[Pi])+94291368595076831942990520 E^(38 I m \[Pi])-18501728410302928319833920 E^(-40 I m \[Pi])-18501728410302928319833920 E^(40 I m \[Pi])+1678910486211891090247320 E^(-42 I m \[Pi])+1678910486211891090247320 E^(42 I m \[Pi]),1523065676 E^(-I m \[Pi])+1523065676 E^(I m \[Pi])-47432377940508 E^(-3 I m \[Pi])-47432377940508 E^(3 I m \[Pi])+38988879445661728 E^(-5 I m \[Pi])+38988879445661728 E^(5 I m \[Pi])-6223192155784060448 E^(-7 I m \[Pi])-6223192155784060448 E^(7 I m \[Pi])+363774705994537152352 E^(-9 I m \[Pi])+363774705994537152352 E^(9 I m \[Pi])-10473415914950561061728 E^(-11 I m \[Pi])-10473415914950561061728 E^(11 I m \[Pi])+175131544378558179284272 E^(-13 I m \[Pi])+175131544378558179284272 E^(13 I m \[Pi])-1882569217830291061434992 E^(-15 I m \[Pi])-1882569217830291061434992 E^(15 I m \[Pi])+13906011478783501391765072 E^(-17 I m \[Pi])+13906011478783501391765072 E^(17 I m \[Pi])-73873697747138801927276944 E^(-19 I m \[Pi])-73873697747138801927276944 E^(19 I m \[Pi])+291307304221409145865626976 E^(-21 I m \[Pi])+291307304221409145865626976 E^(21 I m \[Pi])-871496595540463911780140384 E^(-23 I m \[Pi])-871496595540463911780140384 E^(23 I m \[Pi])+2006633898540980661316683424 E^(-25 I m \[Pi])+2006633898540980661316683424 E^(25 I m \[Pi])-3585139259070202259698561184 E^(-27 I m \[Pi])-3585139259070202259698561184 E^(27 I m \[Pi])+4982687320960827042344534776 E^(-29 I m \[Pi])+4982687320960827042344534776 E^(29 I m \[Pi])-5368963570749883546104720344 E^(-31 I m \[Pi])-5368963570749883546104720344 E^(31 I m \[Pi])+4439907637397011986035443432 E^(-33 I m \[Pi])+4439907637397011986035443432 E^(33 I m \[Pi])-2764095138898603535217873352 E^(-35 I m \[Pi])-2764095138898603535217873352 E^(35 I m \[Pi])+1253244012353802805104687632 E^(-37 I m \[Pi])+1253244012353802805104687632 E^(37 I m \[Pi])-390565358266341830356972880 E^(-39 I m \[Pi])-390565358266341830356972880 E^(39 I m \[Pi])+74802541903753894599211920 E^(-41 I m \[Pi])+74802541903753894599211920 E^(41 I m \[Pi])-6637553085023755473070800 E^(-43 I m \[Pi])-6637553085023755473070800 E^(43 I m \[Pi]),149184-794941610496 E^(-2 I m \[Pi])-794941610496 E^(2 I m \[Pi])+2923398142639788 E^(-4 I m \[Pi])+2923398142639788 E^(4 I m \[Pi])-1047154986392928512 E^(-6 I m \[Pi])-1047154986392928512 E^(6 I m \[Pi])+105358870386048566496 E^(-8 I m \[Pi])+105358870386048566496 E^(8 I m \[Pi])-4566718365731612623104 E^(-10 I m \[Pi])-4566718365731612623104 E^(10 I m \[Pi])+106366343388988256668996 E^(-12 I m \[Pi])+106366343388988256668996 E^(12 I m \[Pi])-1516575688568720003254272 E^(-14 I m \[Pi])-1516575688568720003254272 E^(14 I m \[Pi])+14387016584976539450588928 E^(-16 I m \[Pi])+14387016584976539450588928 E^(16 I m \[Pi])-96049117510864959812749312 E^(-18 I m \[Pi])-96049117510864959812749312 E^(18 I m \[Pi])+469189787498188789762199448 E^(-20 I m \[Pi])+469189787498188789762199448 E^(20 I m \[Pi])-1723454417701969310864666112 E^(-22 I m \[Pi])-1723454417701969310864666112 E^(22 I m \[Pi])+4851092536772706789412837440 E^(-24 I m \[Pi])+4851092536772706789412837440 E^(24 I m \[Pi])-10592284867389471146232245760 E^(-26 I m \[Pi])-10592284867389471146232245760 E^(26 I m \[Pi])+18060482275740041656327949448 E^(-28 I m \[Pi])+18060482275740041656327949448 E^(28 I m \[Pi])-24079209676212750765934728192 E^(-30 I m \[Pi])-24079209676212750765934728192 E^(30 I m \[Pi])+24997384601787805235397925248 E^(-32 I m \[Pi])+24997384601787805235397925248 E^(32 I m \[Pi])-19988264712890438902772093952 E^(-34 I m \[Pi])-19988264712890438902772093952 E^(34 I m \[Pi])+12069412427738578735188583448 E^(-36 I m \[Pi])+12069412427738578735188583448 E^(36 I m \[Pi])-5321642508401164402710274560 E^(-38 I m \[Pi])-5321642508401164402710274560 E^(38 I m \[Pi])+1616482885589372366120533440 E^(-40 I m \[Pi])+1616482885589372366120533440 E^(40 I m \[Pi])-302360065237788014020354560 E^(-42 I m \[Pi])-302360065237788014020354560 E^(42 I m \[Pi])+26248505381684851188961800 E^(-44 I m \[Pi])+26248505381684851188961800 E^(44 I m \[Pi]),-2535875100 E^(-I m \[Pi])-2535875100 E^(I m \[Pi])+111607835539500 E^(-3 I m \[Pi])+111607835539500 E^(3 I m \[Pi])-116456472540278244 E^(-5 I m \[Pi])-116456472540278244 E^(5 I m \[Pi])+22584744903958461000 E^(-7 I m \[Pi])+22584744903958461000 E^(7 I m \[Pi])-1565245909482083931000 E^(-9 I m \[Pi])-1565245909482083931000 E^(9 I m \[Pi])+52627226242039606854000 E^(-11 I m \[Pi])+52627226242039606854000 E^(11 I m \[Pi])-1017808348164735314005200 E^(-13 I m \[Pi])-1017808348164735314005200 E^(13 I m \[Pi])+12578426613223890956729496 E^(-15 I m \[Pi])+12578426613223890956729496 E^(15 I m \[Pi])-106466783303822381311732200 E^(-17 I m \[Pi])-106466783303822381311732200 E^(17 I m \[Pi])+647357573634129311659885800 E^(-19 I m \[Pi])+647357573634129311659885800 E^(19 I m \[Pi])-2924097470822213156581198200 E^(-21 I m \[Pi])-2924097470822213156581198200 E^(21 I m \[Pi])+10047239296177737567979544400 E^(-23 I m \[Pi])+10047239296177737567979544400 E^(23 I m \[Pi])-26693275309796689650143391024 E^(-25 I m \[Pi])-26693275309796689650143391024 E^(25 I m \[Pi])+55409732069140915765543819400 E^(-27 I m \[Pi])+55409732069140915765543819400 E^(27 I m \[Pi])-90342309853689783890861410200 E^(-29 I m \[Pi])-90342309853689783890861410200 E^(29 I m \[Pi])+115732632446036587283102803800 E^(-31 I m \[Pi])+115732632446036587283102803800 E^(31 I m \[Pi])-115905054929149202124450508200 E^(-33 I m \[Pi])-115905054929149202124450508200 E^(33 I m \[Pi])+89711747166481745864093387496 E^(-35 I m \[Pi])+89711747166481745864093387496 E^(35 I m \[Pi])-52587609401264687328601726200 E^(-37 I m \[Pi])-52587609401264687328601726200 E^(37 I m \[Pi])+22565684401899592933382634000 E^(-39 I m \[Pi])+22565684401899592933382634000 E^(39 I m \[Pi])-6685276668331739880777606000 E^(-41 I m \[Pi])-6685276668331739880777606000 E^(41 I m \[Pi])+1221913181561191348451670000 E^(-43 I m \[Pi])+1221913181561191348451670000 E^(43 I m \[Pi])-103827421287553411369671120 E^(-45 I m \[Pi])-103827421287553411369671120 E^(45 I m \[Pi]),-194688+1596229763472 E^(-2 I m \[Pi])+1596229763472 E^(2 I m \[Pi])-7719988118941440 E^(-4 I m \[Pi])-7719988118941440 E^(4 I m \[Pi])+3413709475153622480 E^(-6 I m \[Pi])+3413709475153622480 E^(6 I m \[Pi])-410716365348855782528 E^(-8 I m \[Pi])-410716365348855782528 E^(8 I m \[Pi])+20886242771769886061272 E^(-10 I m \[Pi])+20886242771769886061272 E^(10 I m \[Pi])-563869599338662051486528 E^(-12 I m \[Pi])-563869599338662051486528 E^(12 I m \[Pi])+9246415315171394319132840 E^(-14 I m \[Pi])+9246415315171394319132840 E^(14 I m \[Pi])-100402818025774991129230080 E^(-16 I m \[Pi])-100402818025774991129230080 E^(16 I m \[Pi])+765370614810549970800707376 E^(-18 I m \[Pi])+765370614810549970800707376 E^(18 I m \[Pi])-4266869693939026351327745664 E^(-20 I m \[Pi])-4266869693939026351327745664 E^(20 I m \[Pi])+17909799431298768641910623376 E^(-22 I m \[Pi])+17909799431298768641910623376 E^(22 I m \[Pi])-57779326476069101745979177472 E^(-24 I m \[Pi])-57779326476069101745979177472 E^(24 I m \[Pi])+145311670654828109316415671728 E^(-26 I m \[Pi])+145311670654828109316415671728 E^(26 I m \[Pi])-287416035936763840691766315904 E^(-28 I m \[Pi])-287416035936763840691766315904 E^(28 I m \[Pi])+448928044601082425162742186256 E^(-30 I m \[Pi])+448928044601082425162742186256 E^(30 I m \[Pi])-553399713459114812945744027904 E^(-32 I m \[Pi])-553399713459114812945744027904 E^(32 I m \[Pi])+535316074564475923586117107200 E^(-34 I m \[Pi])+535316074564475923586117107200 E^(34 I m \[Pi])-401480116849529085345229082880 E^(-36 I m \[Pi])-401480116849529085345229082880 E^(36 I m \[Pi])+228659913564174984909460175040 E^(-38 I m \[Pi])+228659913564174984909460175040 E^(38 I m \[Pi])-95559069311072475062068181760 E^(-40 I m \[Pi])-95559069311072475062068181760 E^(40 I m \[Pi])+27628214851988190170200013040 E^(-42 I m \[Pi])+27628214851988190170200013040 E^(42 I m \[Pi])-4937052212235124010297170560 E^(-44 I m \[Pi])-4937052212235124010297170560 E^(44 I m \[Pi])+410795449442059149332177040 E^(-46 I m \[Pi])+410795449442059149332177040 E^(46 I m \[Pi]),4172253168 E^(-I m \[Pi])+4172253168 E^(I m \[Pi])-257515126990640 E^(-3 I m \[Pi])-257515126990640 E^(3 I m \[Pi])+339306697647141360 E^(-5 I m \[Pi])+339306697647141360 E^(5 I m \[Pi])-79598414212709217840 E^(-7 I m \[Pi])-79598414212709217840 E^(7 I m \[Pi])+6513806449319989026784 E^(-9 I m \[Pi])+6513806449319989026784 E^(9 I m \[Pi])-254718630343701785743776 E^(-11 I m \[Pi])-254718630343701785743776 E^(11 I m \[Pi])+5673528602769055830068976 E^(-13 I m \[Pi])+5673528602769055830068976 E^(13 I m \[Pi])-80244261551334066560543216 E^(-15 I m \[Pi])-80244261551334066560543216 E^(15 I m \[Pi])+774417058475371106110766976 E^(-17 I m \[Pi])+774417058475371106110766976 E^(17 I m \[Pi])-5359581800333257274109082752 E^(-19 I m \[Pi])-5359581800333257274109082752 E^(19 I m \[Pi])+27556041704514696634336206208 E^(-21 I m \[Pi])+27556041704514696634336206208 E^(21 I m \[Pi])-107953120629069446020365508224 E^(-23 I m \[Pi])-107953120629069446020365508224 E^(23 I m \[Pi])+328098285951774465360733133424 E^(-25 I m \[Pi])+328098285951774465360733133424 E^(25 I m \[Pi])-783154776755329636086357015024 E^(-27 I m \[Pi])-783154776755329636086357015024 E^(27 I m \[Pi])+1479087594174844719244414935168 E^(-29 I m \[Pi])+1479087594174844719244414935168 E^(29 I m \[Pi])-2216943319678847934582545078592 E^(-31 I m \[Pi])-2216943319678847934582545078592 E^(31 I m \[Pi])+2633387145355152559481288637856 E^(-33 I m \[Pi])+2633387145355152559481288637856 E^(33 I m \[Pi])-2463241941879678582994227456096 E^(-35 I m \[Pi])-2463241941879678582994227456096 E^(35 I m \[Pi])+1791767523522854741081992946256 E^(-37 I m \[Pi])+1791767523522854741081992946256 E^(37 I m \[Pi])-992311187786117970626542018640 E^(-39 I m \[Pi])-992311187786117970626542018640 E^(39 I m \[Pi])+404149326457635790087741997760 E^(-41 I m \[Pi])+404149326457635790087741997760 E^(41 I m \[Pi])-114099778897125084718649709120 E^(-43 I m \[Pi])-114099778897125084718649709120 E^(43 I m \[Pi])+19943893358626563975269869920 E^(-45 I m \[Pi])+19943893358626563975269869920 E^(45 I m \[Pi])-1625701140345170250548615520 E^(-47 I m \[Pi])-1625701140345170250548615520 E^(47 I m \[Pi]),261184-3155088605184 E^(-2 I m \[Pi])-3155088605184 E^(2 I m \[Pi])+19950225954803712 E^(-4 I m \[Pi])+19950225954803712 E^(4 I m \[Pi])-10840620944748412928 E^(-6 I m \[Pi])-10840620944748412928 E^(6 I m \[Pi])+1553319313412409451008 E^(-8 I m \[Pi])+1553319313412409451008 E^(8 I m \[Pi])-92314457810713206816768 E^(-10 I m \[Pi])-92314457810713206816768 E^(10 I m \[Pi])+2877329951688862855794688 E^(-12 I m \[Pi])+2877329951688862855794688 E^(12 I m \[Pi])-54040016150041549673299968 E^(-14 I m \[Pi])-54040016150041549673299968 E^(14 I m \[Pi])+668664276722898241717977012 E^(-16 I m \[Pi])+668664276722898241717977012 E^(16 I m \[Pi])-5791549256927015663711535104 E^(-18 I m \[Pi])-5791549256927015663711535104 E^(18 I m \[Pi])+36644915408560878389204404224 E^(-20 I m \[Pi])+36644915408560878389204404224 E^(20 I m \[Pi])-174658637562626636276229586944 E^(-22 I m \[Pi])-174658637562626636276229586944 E^(22 I m \[Pi])+641131408318020398282944018688 E^(-24 I m \[Pi])+641131408318020398282944018688 E^(24 I m \[Pi])-1841291751907235776388578394112 E^(-26 I m \[Pi])-1841291751907235776388578394112 E^(26 I m \[Pi])+4181425894377756948724940998656 E^(-28 I m \[Pi])+4181425894377756948724940998656 E^(28 I m \[Pi])-7555147389832227075236528766976 E^(-30 I m \[Pi])-7555147389832227075236528766976 E^(30 I m \[Pi])+10883754540162875651774033899296 E^(-32 I m \[Pi])+10883754540162875651774033899296 E^(32 I m \[Pi])-12473729151888746322269135290368 E^(-34 I m \[Pi])-12473729151888746322269135290368 E^(34 I m \[Pi])+11294699147396676491737680381952 E^(-36 I m \[Pi])+11294699147396676491737680381952 E^(36 I m \[Pi])-7975538236203988314950943621120 E^(-38 I m \[Pi])-7975538236203988314950943621120 E^(38 I m \[Pi])+4298293030187155123066518616320 E^(-40 I m \[Pi])+4298293030187155123066518616320 E^(40 I m \[Pi])-1707192927165154065710534737920 E^(-42 I m \[Pi])-1707192927165154065710534737920 E^(42 I m \[Pi])+470899099538917767194363873280 E^(-44 I m \[Pi])+470899099538917767194363873280 E^(44 I m \[Pi])-80550869405489725962666885120 E^(-46 I m \[Pi])-80550869405489725962666885120 E^(46 I m \[Pi])+6435067013866298908421603100 E^(-48 I m \[Pi])+6435067013866298908421603100 E^(48 I m \[Pi]),-6787689426 E^(-I m \[Pi])-6787689426 E^(I m \[Pi])+583379799302944 E^(-3 I m \[Pi])+583379799302944 E^(3 I m \[Pi])-965884755078420576 E^(-5 I m \[Pi])-965884755078420576 E^(5 I m \[Pi])+272966975699026481544 E^(-7 I m \[Pi])+272966975699026481544 E^(7 I m \[Pi])-26275666132128160776504 E^(-9 I m \[Pi])-26275666132128160776504 E^(9 I m \[Pi])+1190564941206415328894256 E^(-11 I m \[Pi])+1190564941206415328894256 E^(11 I m \[Pi])-30423688420429273141505424 E^(-13 I m \[Pi])-30423688420429273141505424 E^(13 I m \[Pi])+490457364900987504584194456 E^(-15 I m \[Pi])+490457364900987504584194456 E^(15 I m \[Pi])-5372951921814803860967035624 E^(-17 I m \[Pi])-5372951921814803860967035624 E^(17 I m \[Pi])+42118227522084691335336284432 E^(-19 I m \[Pi])+42118227522084691335336284432 E^(19 I m \[Pi])-245136064442630730735163153328 E^(-21 I m \[Pi])-245136064442630730735163153328 E^(21 I m \[Pi])+1088090929695923165200978382544 E^(-23 I m \[Pi])+1088090929695923165200978382544 E^(23 I m \[Pi])-3755693139990129139447585373616 E^(-25 I m \[Pi])-3755693139990129139447585373616 E^(25 I m \[Pi])+10220463788671564285926386498544 E^(-27 I m \[Pi])+10220463788671564285926386498544 E^(27 I m \[Pi])-22130193520990608097714612075728 E^(-29 I m \[Pi])-22130193520990608097714612075728 E^(29 I m \[Pi])+38321957619907852501423642779432 E^(-31 I m \[Pi])+38321957619907852501423642779432 E^(31 I m \[Pi])-53136166156636819289067282345944 E^(-33 I m \[Pi])-53136166156636819289067282345944 E^(33 I m \[Pi])+58828514438723714297231529764336 E^(-35 I m \[Pi])+58828514438723714297231529764336 E^(35 I m \[Pi])-51616535777709000992894828977744 E^(-37 I m \[Pi])-51616535777709000992894828977744 E^(37 I m \[Pi])+35412189103206895159825935894520 E^(-39 I m \[Pi])+35412189103206895159825935894520 E^(39 I m \[Pi])-18585353196139832128378514884680 E^(-41 I m \[Pi])-18585353196139832128378514884680 E^(41 I m \[Pi])+7203056611277393239832938956960 E^(-43 I m \[Pi])+7203056611277393239832938956960 E^(43 I m \[Pi])-1942199200317431422907903593440 E^(-45 I m \[Pi])-1942199200317431422907903593440 E^(45 I m \[Pi])+325275703164062814297268821960 E^(-47 I m \[Pi])+325275703164062814297268821960 E^(47 I m \[Pi])-25477612258980856902730428600 E^(-49 I m \[Pi])-25477612258980856902730428600 E^(49 I m \[Pi]),-252016+6144977366250 E^(-2 I m \[Pi])+6144977366250 E^(2 I m \[Pi])-50521562688960000 E^(-4 I m \[Pi])-50521562688960000 E^(4 I m \[Pi])+33590506495715740000 E^(-6 I m \[Pi])+33590506495715740000 E^(6 I m \[Pi])-5710488673907343840000 E^(-8 I m \[Pi])-5710488673907343840000 E^(8 I m \[Pi])+395193596065092003986400 E^(-10 I m \[Pi])+395193596065092003986400 E^(10 I m \[Pi])-14169753916233335647040000 E^(-12 I m \[Pi])-14169753916233335647040000 E^(12 I m \[Pi])+303658855987104994887360000 E^(-14 I m \[Pi])+303658855987104994887360000 E^(14 I m \[Pi])-4264359173882435233746240000 E^(-16 I m \[Pi])-4264359173882435233746240000 E^(16 I m \[Pi])+41782919873011271023248640000 E^(-18 I m \[Pi])+41782919873011271023248640000 E^(18 I m \[Pi])-298599999506758935384174297600 E^(-20 I m \[Pi])-298599999506758935384174297600 E^(20 I m \[Pi])+1607260615983851577369359040000 E^(-22 I m \[Pi])+1607260615983851577369359040000 E^(22 I m \[Pi])-6671211090913547704836984640000 E^(-24 I m \[Pi])-6671211090913547704836984640000 E^(24 I m \[Pi])+21720940498886398101403179360000 E^(-26 I m \[Pi])+21720940498886398101403179360000 E^(26 I m \[Pi])-56150833056320255238973968960000 E^(-28 I m \[Pi])-56150833056320255238973968960000 E^(28 I m \[Pi])+116161119401202127707236362886900 E^(-30 I m \[Pi])+116161119401202127707236362886900 E^(30 I m \[Pi])-193099534328735782631154130860000 E^(-32 I m \[Pi])-193099534328735782631154130860000 E^(32 I m \[Pi])+258058714152116476796283736072500 E^(-34 I m \[Pi])+258058714152116476796283736072500 E^(34 I m \[Pi])-276302961080194546106639152640000 E^(-36 I m \[Pi])-276302961080194546106639152640000 E^(36 I m \[Pi])+235136398922450138917590600450000 E^(-38 I m \[Pi])+235136398922450138917590600450000 E^(38 I m \[Pi])-156859185207322308022677801072000 E^(-40 I m \[Pi])-156859185207322308022677801072000 E^(40 I m \[Pi])+80224360644650213198438289450000 E^(-42 I m \[Pi])+80224360644650213198438289450000 E^(42 I m \[Pi])-30357485102957641288867166400000 E^(-44 I m \[Pi])-30357485102957641288867166400000 E^(44 I m \[Pi])+8005594829219916001945551225000 E^(-46 I m \[Pi])+8005594829219916001945551225000 E^(46 I m \[Pi])-1313278982421693654779919000000 E^(-48 I m \[Pi])-1313278982421693654779919000000 E^(48 I m \[Pi])+100891344545564193334812497256 E^(-50 I m \[Pi])+100891344545564193334812497256 E^(50 I m \[Pi])};
FourierCp[m_]:=FourierCp[m]={0,(2 E^((I m \[Pi])/2))/(1+E^(I m \[Pi])),(8 E^((I m \[Pi])/2) (3-2 E^(I m \[Pi])+3 E^(2 I m \[Pi])))/(1+E^(I m \[Pi]))^3,(16 E^((I m \[Pi])/2) (9-44 E^(I m \[Pi])+86 E^(2 I m \[Pi])-44 E^(3 I m \[Pi])+9 E^(4 I m \[Pi])))/(1+E^(I m \[Pi]))^5,(32 E^((I m \[Pi])/2) (19-310 E^(I m \[Pi])+1277 E^(2 I m \[Pi])-1908 E^(3 I m \[Pi])+1277 E^(4 I m \[Pi])-310 E^(5 I m \[Pi])+19 E^(6 I m \[Pi])))/(1+E^(I m \[Pi]))^7,1/(1+E^(I m \[Pi]))^9 4 E^((I m \[Pi])/2) (525-21400 E^(I m \[Pi])+180588 E^(2 I m \[Pi])-551208 E^(3 I m \[Pi])+786318 E^(4 I m \[Pi])-551208 E^(5 I m \[Pi])+180588 E^(6 I m \[Pi])-21400 E^(7 I m \[Pi])+525 E^(8 I m \[Pi])),1/(1+E^(I m \[Pi]))^11 64 E^((I m \[Pi])/2) (99-8682 E^(I m \[Pi])+136743 E^(2 I m \[Pi])-751992 E^(3 I m \[Pi])+1927798 E^(4 I m \[Pi])-2606908 E^(5 I m \[Pi])+1927798 E^(6 I m \[Pi])-751992 E^(7 I m \[Pi])+136743 E^(8 I m \[Pi])-8682 E^(9 I m \[Pi])+99 E^(10 I m \[Pi])),1/(1+E^(I m \[Pi]))^13 32 E^((I m \[Pi])/2) (539-93100 E^(I m \[Pi])+2530262 E^(2 I m \[Pi])-23133436 E^(3 I m \[Pi])+97190069 E^(4 I m \[Pi])-218991448 E^(5 I m \[Pi])+285006516 E^(6 I m \[Pi])-218991448 E^(7 I m \[Pi])+97190069 E^(8 I m \[Pi])-23133436 E^(9 I m \[Pi])+2530262 E^(10 I m \[Pi])-93100 E^(11 I m \[Pi])+539 E^(12 I m \[Pi])),1/(1+E^(I m \[Pi]))^15 128 E^((I m \[Pi])/2) (339-108446 E^(I m \[Pi])+4793121 E^(2 I m \[Pi])-68637516 E^(3 I m \[Pi])+444153259 E^(4 I m \[Pi])-1541041218 E^(5 I m \[Pi])+3149770209 E^(6 I m \[Pi])-3977843112 E^(7 I m \[Pi])+3149770209 E^(8 I m \[Pi])-1541041218 E^(9 I m \[Pi])+444153259 E^(10 I m \[Pi])-68637516 E^(11 I m \[Pi])+4793121 E^(12 I m \[Pi])-108446 E^(13 I m \[Pi])+339 E^(14 I m \[Pi])),1/(1+E^(I m \[Pi]))^17 2 E^((I m \[Pi])/2) (51273-29002608 E^(I m \[Pi])+1992458808 E^(2 I m \[Pi])-42689400912 E^(3 I m \[Pi])+406188730108 E^(4 I m \[Pi])-2059800949360 E^(5 I m \[Pi])+6188060019592 E^(6 I m \[Pi])-11716536725968 E^(7 I m \[Pi])+14445634422262 E^(8 I m \[Pi])-11716536725968 E^(9 I m \[Pi])+6188060019592 E^(10 I m \[Pi])-2059800949360 E^(11 I m \[Pi])+406188730108 E^(12 I m \[Pi])-42689400912 E^(13 I m \[Pi])+1992458808 E^(14 I m \[Pi])-29002608 E^(15 I m \[Pi])+51273 E^(16 I m \[Pi])),1/(1+E^(I m \[Pi]))^19 16 E^((I m \[Pi])/2) (14375-13876650 E^(I m \[Pi])+1430341839 E^(2 I m \[Pi])-44232553584 E^(3 I m \[Pi])+596763971244 E^(4 I m \[Pi])-4255638735384 E^(5 I m \[Pi])+17970084887964 E^(6 I m \[Pi])-48219592213776 E^(7 I m \[Pi])+85842376602306 E^(8 I m \[Pi])-103782353468796 E^(9 I m \[Pi])+85842376602306 E^(10 I m \[Pi])-48219592213776 E^(11 I m \[Pi])+17970084887964 E^(12 I m \[Pi])-4255638735384 E^(13 I m \[Pi])+596763971244 E^(14 I m \[Pi])-44232553584 E^(15 I m \[Pi])+1430341839 E^(16 I m \[Pi])-13876650 E^(17 I m \[Pi])+14375 E^(18 I m \[Pi])),1/(1+E^(I m \[Pi]))^21 16 E^((I m \[Pi])/2) (30855-49383972 E^(I m \[Pi])+7424920338 E^(2 I m \[Pi])-321991528452 E^(3 I m \[Pi])+5982311835643 E^(4 I m \[Pi])-58204809771984 E^(5 I m \[Pi])+334212130353240 E^(6 I m \[Pi])-1222556154890000 E^(7 I m \[Pi])+2994438880764030 E^(8 I m \[Pi])-5067084812411640 E^(9 I m \[Pi])+6027054155892524 E^(10 I m \[Pi])-5067084812411640 E^(11 I m \[Pi])+2994438880764030 E^(12 I m \[Pi])-1222556154890000 E^(13 I m \[Pi])+334212130353240 E^(14 I m \[Pi])-58204809771984 E^(15 I m \[Pi])+5982311835643 E^(16 I m \[Pi])-321991528452 E^(17 I m \[Pi])+7424920338 E^(18 I m \[Pi])-49383972 E^(19 I m \[Pi])+30855 E^(20 I m \[Pi])),1/(1+E^(I m \[Pi]))^23 256 E^((I m \[Pi])/2) (3987-10325510 E^(I m \[Pi])+2212943621 E^(2 I m \[Pi])-131446129508 E^(3 I m \[Pi])+3283522685897 E^(4 I m \[Pi])-42530517095438 E^(5 I m \[Pi])+323558115682991 E^(6 I m \[Pi])-1567432055768752 E^(7 I m \[Pi])+5105286393063790 E^(8 I m \[Pi])-11596375952071660 E^(9 I m \[Pi])+18806328713083266 E^(10 I m \[Pi])-22063977947951064 E^(11 I m \[Pi])+18806328713083266 E^(12 I m \[Pi])-11596375952071660 E^(13 I m \[Pi])+5105286393063790 E^(14 I m \[Pi])-1567432055768752 E^(15 I m \[Pi])+323558115682991 E^(16 I m \[Pi])-42530517095438 E^(17 I m \[Pi])+3283522685897 E^(18 I m \[Pi])-131446129508 E^(19 I m \[Pi])+2212943621 E^(20 I m \[Pi])-10325510 E^(21 I m \[Pi])+3987 E^(22 I m \[Pi])),1/(1+E^(I m \[Pi]))^25 4 E^((I m \[Pi])/2) (510549-2095741960 E^(I m \[Pi])+628106103716 E^(2 I m \[Pi])-50112338131960 E^(3 I m \[Pi])+1649878310833194 E^(4 I m \[Pi])-27877034535496536 E^(5 I m \[Pi])+275100214374269460 E^(6 I m \[Pi])-1725063756991131304 E^(7 I m \[Pi])+7282303993155468635 E^(8 I m \[Pi])-21545575902663906256 E^(9 I m \[Pi])+45931536869613018696 E^(10 I m \[Pi])-71845464081266661936 E^(11 I m \[Pi])+83306878614087794764 E^(12 I m \[Pi])-71845464081266661936 E^(13 I m \[Pi])+45931536869613018696 E^(14 I m \[Pi])-21545575902663906256 E^(15 I m \[Pi])+7282303993155468635 E^(16 I m \[Pi])-1725063756991131304 E^(17 I m \[Pi])+275100214374269460 E^(18 I m \[Pi])-27877034535496536 E^(19 I m \[Pi])+1649878310833194 E^(20 I m \[Pi])-50112338131960 E^(21 I m \[Pi])+628106103716 E^(22 I m \[Pi])-2095741960 E^(23 I m \[Pi])+510549 E^(24 I m \[Pi])),1/(1+E^(I m \[Pi]))^27 128 E^((I m \[Pi])/2) (31017-198202158 E^(I m \[Pi])+81734601165 E^(2 I m \[Pi])-8615825894520 E^(3 I m \[Pi])+367637165367190 E^(4 I m \[Pi])-7965289581598620 E^(5 I m \[Pi])+100176118079770374 E^(6 I m \[Pi])-798102741874117272 E^(7 I m \[Pi])+4278830702958273243 E^(8 I m \[Pi])-16114872645064891122 E^(9 I m \[Pi])+43965388393223503375 E^(10 I m \[Pi])-88785476801469346032 E^(11 I m \[Pi])+134624527240441834884 E^(12 I m \[Pi])-154525728158977336456 E^(13 I m \[Pi])+134624527240441834884 E^(14 I m \[Pi])-88785476801469346032 E^(15 I m \[Pi])+43965388393223503375 E^(16 I m \[Pi])-16114872645064891122 E^(17 I m \[Pi])+4278830702958273243 E^(18 I m \[Pi])-798102741874117272 E^(19 I m \[Pi])+100176118079770374 E^(20 I m \[Pi])-7965289581598620 E^(21 I m \[Pi])+367637165367190 E^(22 I m \[Pi])-8615825894520 E^(23 I m \[Pi])+81734601165 E^(24 I m \[Pi])-198202158 E^(25 I m \[Pi])+31017 E^(26 I m \[Pi])),1/(1+E^(I m \[Pi]))^29 32 E^((I m \[Pi])/2) (235125-2301836100 E^(I m \[Pi])+1288062243042 E^(2 I m \[Pi])-176887988730132 E^(3 I m \[Pi])+9642798900969527 E^(4 I m \[Pi])-264007626918579032 E^(5 I m \[Pi])+4168506397510013492 E^(6 I m \[Pi])-41540427519017730072 E^(7 I m \[Pi])+278165415362180430237 E^(8 I m \[Pi])-1309394360545833112332 E^(9 I m \[Pi])+4478174507687018041662 E^(10 I m \[Pi])-11398459758880252780572 E^(11 I m \[Pi])+21965751331569070268343 E^(12 I m \[Pi])-32417617471379624774928 E^(13 I m \[Pi])+36882013595481880344408 E^(14 I m \[Pi])-32417617471379624774928 E^(15 I m \[Pi])+21965751331569070268343 E^(16 I m \[Pi])-11398459758880252780572 E^(17 I m \[Pi])+4478174507687018041662 E^(18 I m \[Pi])-1309394360545833112332 E^(19 I m \[Pi])+278165415362180430237 E^(20 I m \[Pi])-41540427519017730072 E^(21 I m \[Pi])+4168506397510013492 E^(22 I m \[Pi])-264007626918579032 E^(23 I m \[Pi])+9642798900969527 E^(24 I m \[Pi])-176887988730132 E^(25 I m \[Pi])+1288062243042 E^(26 I m \[Pi])-2301836100 E^(27 I m \[Pi])+235125 E^(28 I m \[Pi])),1/(1+E^(I m \[Pi]))^31 512 E^((I m \[Pi])/2) (27219-402476654 E^(I m \[Pi])+301938287401 E^(2 I m \[Pi])-53367092987324 E^(3 I m \[Pi])+3670972566717911 E^(4 I m \[Pi])-125412472873527322 E^(5 I m \[Pi])+2454135460748686781 E^(6 I m \[Pi])-30185316989010584344 E^(7 I m \[Pi])+248949787158125896391 E^(8 I m \[Pi])-1442799969305244620574 E^(9 I m \[Pi])+6083804570111142210965 E^(10 I m \[Pi])-19155536618519305731076 E^(11 I m \[Pi])+45911038619710267716139 E^(12 I m \[Pi])-84929746782853208511066 E^(13 I m \[Pi])+122400550533944344065129 E^(14 I m \[Pi])-138176814328302916597328 E^(15 I m \[Pi])+122400550533944344065129 E^(16 I m \[Pi])-84929746782853208511066 E^(17 I m \[Pi])+45911038619710267716139 E^(18 I m \[Pi])-19155536618519305731076 E^(19 I m \[Pi])+6083804570111142210965 E^(20 I m \[Pi])-1442799969305244620574 E^(21 I m \[Pi])+248949787158125896391 E^(22 I m \[Pi])-30185316989010584344 E^(23 I m \[Pi])+2454135460748686781 E^(24 I m \[Pi])-125412472873527322 E^(25 I m \[Pi])+3670972566717911 E^(26 I m \[Pi])-53367092987324 E^(27 I m \[Pi])+301938287401 E^(28 I m \[Pi])-402476654 E^(29 I m \[Pi])+27219 E^(30 I m \[Pi]))};
FourierCpp[m_]:=FourierCpp[m]={0, (-8*E^((I/2)*m*Pi))/(1 - E^(I*m*Pi)), 
 (32*E^((I/2)*m*Pi)*(3 + 2*E^(I*m*Pi) + 3*E^((2*I)*m*Pi)))/(1 - E^(I*m*Pi))^3, 
 (-64*E^((I/2)*m*Pi)*(9 + 44*E^(I*m*Pi) + 86*E^((2*I)*m*Pi) + 44*E^((3*I)*m*Pi) + 
    9*E^((4*I)*m*Pi)))/(1 - E^(I*m*Pi))^5, 
 (128*E^((I/2)*m*Pi)*(19 + 310*E^(I*m*Pi) + 1277*E^((2*I)*m*Pi) + 1908*E^((3*I)*m*Pi) + 
    1277*E^((4*I)*m*Pi) + 310*E^((5*I)*m*Pi) + 19*E^((6*I)*m*Pi)))/(1 - E^(I*m*Pi))^7, 
 (-16*E^((I/2)*m*Pi)*(525 + 21400*E^(I*m*Pi) + 180588*E^((2*I)*m*Pi) + 551208*E^((3*I)*m*Pi) + 
    786318*E^((4*I)*m*Pi) + 551208*E^((5*I)*m*Pi) + 180588*E^((6*I)*m*Pi) + 
    21400*E^((7*I)*m*Pi) + 525*E^((8*I)*m*Pi)))/(1 - E^(I*m*Pi))^9, 
 (256*E^((I/2)*m*Pi)*(99 + 8682*E^(I*m*Pi) + 136743*E^((2*I)*m*Pi) + 751992*E^((3*I)*m*Pi) + 
    1927798*E^((4*I)*m*Pi) + 2606908*E^((5*I)*m*Pi) + 1927798*E^((6*I)*m*Pi) + 
    751992*E^((7*I)*m*Pi) + 136743*E^((8*I)*m*Pi) + 8682*E^((9*I)*m*Pi) + 99*E^((10*I)*m*Pi)))/
  (1 - E^(I*m*Pi))^11, (-128*E^((I/2)*m*Pi)*(539 + 93100*E^(I*m*Pi) + 2530262*E^((2*I)*m*Pi) + 
    23133436*E^((3*I)*m*Pi) + 97190069*E^((4*I)*m*Pi) + 218991448*E^((5*I)*m*Pi) + 
    285006516*E^((6*I)*m*Pi) + 218991448*E^((7*I)*m*Pi) + 97190069*E^((8*I)*m*Pi) + 
    23133436*E^((9*I)*m*Pi) + 2530262*E^((10*I)*m*Pi) + 93100*E^((11*I)*m*Pi) + 
    539*E^((12*I)*m*Pi)))/(1 - E^(I*m*Pi))^13, 
 (512*E^((I/2)*m*Pi)*(339 + 108446*E^(I*m*Pi) + 4793121*E^((2*I)*m*Pi) + 
    68637516*E^((3*I)*m*Pi) + 444153259*E^((4*I)*m*Pi) + 1541041218*E^((5*I)*m*Pi) + 
    3149770209*E^((6*I)*m*Pi) + 3977843112*E^((7*I)*m*Pi) + 3149770209*E^((8*I)*m*Pi) + 
    1541041218*E^((9*I)*m*Pi) + 444153259*E^((10*I)*m*Pi) + 68637516*E^((11*I)*m*Pi) + 
    4793121*E^((12*I)*m*Pi) + 108446*E^((13*I)*m*Pi) + 339*E^((14*I)*m*Pi)))/
  (1 - E^(I*m*Pi))^15, (-8*E^((I/2)*m*Pi)*(51273 + 29002608*E^(I*m*Pi) + 
    1992458808*E^((2*I)*m*Pi) + 42689400912*E^((3*I)*m*Pi) + 406188730108*E^((4*I)*m*Pi) + 
    2059800949360*E^((5*I)*m*Pi) + 6188060019592*E^((6*I)*m*Pi) + 
    11716536725968*E^((7*I)*m*Pi) + 14445634422262*E^((8*I)*m*Pi) + 
    11716536725968*E^((9*I)*m*Pi) + 6188060019592*E^((10*I)*m*Pi) + 
    2059800949360*E^((11*I)*m*Pi) + 406188730108*E^((12*I)*m*Pi) + 
    42689400912*E^((13*I)*m*Pi) + 1992458808*E^((14*I)*m*Pi) + 29002608*E^((15*I)*m*Pi) + 
    51273*E^((16*I)*m*Pi)))/(1 - E^(I*m*Pi))^17, 
 (64*E^((I/2)*m*Pi)*(14375 + 13876650*E^(I*m*Pi) + 1430341839*E^((2*I)*m*Pi) + 
    44232553584*E^((3*I)*m*Pi) + 596763971244*E^((4*I)*m*Pi) + 4255638735384*E^((5*I)*m*Pi) + 
    17970084887964*E^((6*I)*m*Pi) + 48219592213776*E^((7*I)*m*Pi) + 
    85842376602306*E^((8*I)*m*Pi) + 103782353468796*E^((9*I)*m*Pi) + 
    85842376602306*E^((10*I)*m*Pi) + 48219592213776*E^((11*I)*m*Pi) + 
    17970084887964*E^((12*I)*m*Pi) + 4255638735384*E^((13*I)*m*Pi) + 
    596763971244*E^((14*I)*m*Pi) + 44232553584*E^((15*I)*m*Pi) + 1430341839*E^((16*I)*m*Pi) + 
    13876650*E^((17*I)*m*Pi) + 14375*E^((18*I)*m*Pi)))/(1 - E^(I*m*Pi))^19, 
 (-64*E^((I/2)*m*Pi)*(30855 + 49383972*E^(I*m*Pi) + 7424920338*E^((2*I)*m*Pi) + 
    321991528452*E^((3*I)*m*Pi) + 5982311835643*E^((4*I)*m*Pi) + 
    58204809771984*E^((5*I)*m*Pi) + 334212130353240*E^((6*I)*m*Pi) + 
    1222556154890000*E^((7*I)*m*Pi) + 2994438880764030*E^((8*I)*m*Pi) + 
    5067084812411640*E^((9*I)*m*Pi) + 6027054155892524*E^((10*I)*m*Pi) + 
    5067084812411640*E^((11*I)*m*Pi) + 2994438880764030*E^((12*I)*m*Pi) + 
    1222556154890000*E^((13*I)*m*Pi) + 334212130353240*E^((14*I)*m*Pi) + 
    58204809771984*E^((15*I)*m*Pi) + 5982311835643*E^((16*I)*m*Pi) + 
    321991528452*E^((17*I)*m*Pi) + 7424920338*E^((18*I)*m*Pi) + 49383972*E^((19*I)*m*Pi) + 
    30855*E^((20*I)*m*Pi)))/(1 - E^(I*m*Pi))^21, 
 (1024*E^((I/2)*m*Pi)*(3987 + 10325510*E^(I*m*Pi) + 2212943621*E^((2*I)*m*Pi) + 
    131446129508*E^((3*I)*m*Pi) + 3283522685897*E^((4*I)*m*Pi) + 
    42530517095438*E^((5*I)*m*Pi) + 323558115682991*E^((6*I)*m*Pi) + 
    1567432055768752*E^((7*I)*m*Pi) + 5105286393063790*E^((8*I)*m*Pi) + 
    11596375952071660*E^((9*I)*m*Pi) + 18806328713083266*E^((10*I)*m*Pi) + 
    22063977947951064*E^((11*I)*m*Pi) + 18806328713083266*E^((12*I)*m*Pi) + 
    11596375952071660*E^((13*I)*m*Pi) + 5105286393063790*E^((14*I)*m*Pi) + 
    1567432055768752*E^((15*I)*m*Pi) + 323558115682991*E^((16*I)*m*Pi) + 
    42530517095438*E^((17*I)*m*Pi) + 3283522685897*E^((18*I)*m*Pi) + 
    131446129508*E^((19*I)*m*Pi) + 2212943621*E^((20*I)*m*Pi) + 10325510*E^((21*I)*m*Pi) + 
    3987*E^((22*I)*m*Pi)))/(1 - E^(I*m*Pi))^23, 
 (-16*E^((I/2)*m*Pi)*(510549 + 2095741960*E^(I*m*Pi) + 628106103716*E^((2*I)*m*Pi) + 
    50112338131960*E^((3*I)*m*Pi) + 1649878310833194*E^((4*I)*m*Pi) + 
    27877034535496536*E^((5*I)*m*Pi) + 275100214374269460*E^((6*I)*m*Pi) + 
    1725063756991131304*E^((7*I)*m*Pi) + 7282303993155468635*E^((8*I)*m*Pi) + 
    21545575902663906256*E^((9*I)*m*Pi) + 45931536869613018696*E^((10*I)*m*Pi) + 
    71845464081266661936*E^((11*I)*m*Pi) + 83306878614087794764*E^((12*I)*m*Pi) + 
    71845464081266661936*E^((13*I)*m*Pi) + 45931536869613018696*E^((14*I)*m*Pi) + 
    21545575902663906256*E^((15*I)*m*Pi) + 7282303993155468635*E^((16*I)*m*Pi) + 
    1725063756991131304*E^((17*I)*m*Pi) + 275100214374269460*E^((18*I)*m*Pi) + 
    27877034535496536*E^((19*I)*m*Pi) + 1649878310833194*E^((20*I)*m*Pi) + 
    50112338131960*E^((21*I)*m*Pi) + 628106103716*E^((22*I)*m*Pi) + 
    2095741960*E^((23*I)*m*Pi) + 510549*E^((24*I)*m*Pi)))/(1 - E^(I*m*Pi))^25, 
 (512*E^((I/2)*m*Pi)*(31017 + 198202158*E^(I*m*Pi) + 81734601165*E^((2*I)*m*Pi) + 
    8615825894520*E^((3*I)*m*Pi) + 367637165367190*E^((4*I)*m*Pi) + 
    7965289581598620*E^((5*I)*m*Pi) + 100176118079770374*E^((6*I)*m*Pi) + 
    798102741874117272*E^((7*I)*m*Pi) + 4278830702958273243*E^((8*I)*m*Pi) + 
    16114872645064891122*E^((9*I)*m*Pi) + 43965388393223503375*E^((10*I)*m*Pi) + 
    88785476801469346032*E^((11*I)*m*Pi) + 134624527240441834884*E^((12*I)*m*Pi) + 
    154525728158977336456*E^((13*I)*m*Pi) + 134624527240441834884*E^((14*I)*m*Pi) + 
    88785476801469346032*E^((15*I)*m*Pi) + 43965388393223503375*E^((16*I)*m*Pi) + 
    16114872645064891122*E^((17*I)*m*Pi) + 4278830702958273243*E^((18*I)*m*Pi) + 
    798102741874117272*E^((19*I)*m*Pi) + 100176118079770374*E^((20*I)*m*Pi) + 
    7965289581598620*E^((21*I)*m*Pi) + 367637165367190*E^((22*I)*m*Pi) + 
    8615825894520*E^((23*I)*m*Pi) + 81734601165*E^((24*I)*m*Pi) + 198202158*E^((25*I)*m*Pi) + 
    31017*E^((26*I)*m*Pi)))/(1 - E^(I*m*Pi))^27, 
 (-128*E^((I/2)*m*Pi)*(235125 + 2301836100*E^(I*m*Pi) + 1288062243042*E^((2*I)*m*Pi) + 
    176887988730132*E^((3*I)*m*Pi) + 9642798900969527*E^((4*I)*m*Pi) + 
    264007626918579032*E^((5*I)*m*Pi) + 4168506397510013492*E^((6*I)*m*Pi) + 
    41540427519017730072*E^((7*I)*m*Pi) + 278165415362180430237*E^((8*I)*m*Pi) + 
    1309394360545833112332*E^((9*I)*m*Pi) + 4478174507687018041662*E^((10*I)*m*Pi) + 
    11398459758880252780572*E^((11*I)*m*Pi) + 21965751331569070268343*E^((12*I)*m*Pi) + 
    32417617471379624774928*E^((13*I)*m*Pi) + 36882013595481880344408*E^((14*I)*m*Pi) + 
    32417617471379624774928*E^((15*I)*m*Pi) + 21965751331569070268343*E^((16*I)*m*Pi) + 
    11398459758880252780572*E^((17*I)*m*Pi) + 4478174507687018041662*E^((18*I)*m*Pi) + 
    1309394360545833112332*E^((19*I)*m*Pi) + 278165415362180430237*E^((20*I)*m*Pi) + 
    41540427519017730072*E^((21*I)*m*Pi) + 4168506397510013492*E^((22*I)*m*Pi) + 
    264007626918579032*E^((23*I)*m*Pi) + 9642798900969527*E^((24*I)*m*Pi) + 
    176887988730132*E^((25*I)*m*Pi) + 1288062243042*E^((26*I)*m*Pi) + 
    2301836100*E^((27*I)*m*Pi) + 235125*E^((28*I)*m*Pi)))/(1 - E^(I*m*Pi))^29, 
 (2048*E^((I/2)*m*Pi)*(27219 + 402476654*E^(I*m*Pi) + 301938287401*E^((2*I)*m*Pi) + 
    53367092987324*E^((3*I)*m*Pi) + 3670972566717911*E^((4*I)*m*Pi) + 
    125412472873527322*E^((5*I)*m*Pi) + 2454135460748686781*E^((6*I)*m*Pi) + 
    30185316989010584344*E^((7*I)*m*Pi) + 248949787158125896391*E^((8*I)*m*Pi) + 
    1442799969305244620574*E^((9*I)*m*Pi) + 6083804570111142210965*E^((10*I)*m*Pi) + 
    19155536618519305731076*E^((11*I)*m*Pi) + 45911038619710267716139*E^((12*I)*m*Pi) + 
    84929746782853208511066*E^((13*I)*m*Pi) + 122400550533944344065129*E^((14*I)*m*Pi) + 
    138176814328302916597328*E^((15*I)*m*Pi) + 122400550533944344065129*E^((16*I)*m*Pi) + 
    84929746782853208511066*E^((17*I)*m*Pi) + 45911038619710267716139*E^((18*I)*m*Pi) + 
    19155536618519305731076*E^((19*I)*m*Pi) + 6083804570111142210965*E^((20*I)*m*Pi) + 
    1442799969305244620574*E^((21*I)*m*Pi) + 248949787158125896391*E^((22*I)*m*Pi) + 
    30185316989010584344*E^((23*I)*m*Pi) + 2454135460748686781*E^((24*I)*m*Pi) + 
    125412472873527322*E^((25*I)*m*Pi) + 3670972566717911*E^((26*I)*m*Pi) + 
    53367092987324*E^((27*I)*m*Pi) + 301938287401*E^((28*I)*m*Pi) + 402476654*E^((29*I)*m*Pi) + 
    27219*E^((30*I)*m*Pi)))/(1 - E^(I*m*Pi))^31};

QuantumVolume[m_]:=QuantumVolume[m]=2 I/Pi^2(PolyLog[2,-I Exp[I Pi m/2]]-PolyLog[2,I Exp[I Pi m/2]]);
QuantumVolumep[m_]:=QuantumVolumep[m]=-2 I/Pi^2(PolyLog[2,Exp[I Pi m/2]]-PolyLog[2,-Exp[I Pi m/2]]);


(** The definitions of EichlerTLV, EichlerTDLV, EichlerTp, EichlerTDp, EichlerTpp, EichlerTDpp are no longer used. **)
(** The definitions of Eichlerf1, Eichlerf2, Eichlerf1p, Eichlerf2p, Eichlerf1pp, Eichlerf2pp are no longer used. **)
(** The definitions of om2B, om3B, EichlerZ0, EichlerT0, EichlerTD0 are no longer used. **)
(** If you decide to delete them, also delete the ::usage **)

With[{nmax=50},
  Eichlerf1[tau_,m_]:=((Eichlerf1[ttau_,m]:=#[ttau])&[Compile[{{tttau,_Complex}},With[{qt=Exp[(2. I Pi) tttau+#1]}, #2]]&[2. Pi Im[tauB[m]], HornerForm[Sum[qt^n FourierCn[m,n]/n,{n,nmax,1,-1}]]]];  Eichlerf1[tau,m]);
  Eichlerf2[tau_,m_]:=((Eichlerf2[ttau_,m]:=#[ttau])&[Compile[{{tttau,_Complex}},With[{qt=Exp[(2. I Pi) tttau+#1]}, -#2]]&[2. Pi Im[tauB[m]], HornerForm[Sum[qt^n FourierCn[m,n]/n^2,{n,nmax,1,-1}]]]];  Eichlerf2[tau,m]);
  Eichlerf1p[taup_,m_]:=((Eichlerf1p[ttau_,m]:=#[ttau])&[Compile[{{tttau,_Complex}},With[{qt=Exp[(2. I Pi) tttau+#1]}, #2]]&[2. Pi Im[-.25/(tauB[m]-Boole[Re[tauB[m]]>=1/2])], HornerForm[Sum[qt^n FourierCpn[m,n]/n,{n,nmax,1,-1}]]]];  Eichlerf1p[taup,m]);
  Eichlerf2p[taup_,m_]:=((Eichlerf2p[ttau_,m]:=#[ttau])&[Compile[{{tttau,_Complex}},With[{qt=Exp[(2. I Pi) tttau+#1]}, (2. I Pi)tttau #2 - #3]]&[2. Pi Im[-.25/(tauB[m]-Boole[Re[tauB[m]]>=1/2])], HornerForm[Sum[qt^n FourierCpn[m,n]/n,{n,nmax,1,-1}]], HornerForm[Sum[qt^n FourierCpn[m,n]/n^2,{n,nmax,1,-1}]]]];  Eichlerf2p[taup,m]);
  Eichlerf1pp[taupp_,m_]:=((Eichlerf1pp[ttau_,m]:=#[ttau])&[Compile[{{tttau,_Complex}},With[{qt=Exp[(2. I Pi) tttau+#1]}, #2]]&[2. Pi Im[tauB[m]/(1-2tauB[m])], HornerForm[Sum[qt^n FourierCppn[m,n]/n,{n,nmax,1,-1}]]]];  Eichlerf1pp[taupp,m]);
  Eichlerf2pp[taupp_,m_]:=((Eichlerf2pp[ttau_,m]:=#[ttau])&[Compile[{{tttau,_Complex}},With[{qt=Exp[(2. I Pi) tttau+#1]}, 2. I Pi tttau #2 - #3]]&[2. Pi Im[tauB[m]/(1-2tauB[m])], HornerForm[Sum[qt^n FourierCppn[m,n]/n,{n,nmax,1,-1}]], HornerForm[Sum[qt^n FourierCppn[m,n]/n^2,{n,nmax,1,-1}]]]];  Eichlerf2pp[taupp,m]);
]

EichlerTLV[tau_,m_]:=tau-m/2+N[Eichlerf1[tau,m]/(2Pi I)];
EichlerTDLV[tau_,m_]:= 1/2tau^2 -1/8m^2+ 1/12+N[ tau Eichlerf1[tau,m]/(2Pi I)+Eichlerf2[tau,m]/(2Pi I)^2];
EichlerTp[tau_,m_]:=N[-8 I/(Pi^2 ) Eichlerf2p[-1/(4tau),m]+I*QuantumVolume[m]];
EichlerTDp[tau_,m_]:=N[-4/Pi Eichlerf1p[-1/(4tau),m]];
EichlerTpp[tau_,m_]:=N[-2/Pi^2 Eichlerf2pp[tau/(1-2tau),m]-2I/Pi Eichlerf1pp[tau/(1-2tau),m]+I*QuantumVolumep[m]];
EichlerTDpp[tau_,m_]:=N[ -1/Pi^2 Eichlerf2pp[tau/(1-2tau),m]+I/2 QuantumVolumep[m]];


om2B[tau_,m_,NMax_:20] :=Module[{an,anp1,anp2,om2,y,z, tb},
tb = tauB[m];
y=Exp[2Pi I m]; z=N[1/(Exp[I Pi m]ModularJ4[tau]+4(Exp[2Pi I m]+1))];
an=2EllipticK[1-y]/Pi;
anp1=(-2 EllipticE[1-y]+(1+y) EllipticK[1-y])/(12 \[Pi] (-1+y)^2);
om2=an+anp1/z;
Do[
anp2=(-(1+2 n)^3 an+32 (1+n)^2 (3+2 n) (1+y) anp1)/(64 (2+n) (5+7 n+2 n^2) (-1+y)^2);
om2+=anp2/z^(n+2);an=anp1; anp1=anp2;,{n,0,NMax}];om2 Sqrt[I(tau/(1-2 tau)-tb/(1-2 tb))]/Sqrt[I(tau/(1-2 tau)-tb/(1-2 tb))z]];

om3B[tau_,m_,NMax_:20] :=Module[{an,anp1,anp2,om3,y,z, tb},
tb = tauB[m];
y=Exp[2Pi I m]; z=N[1/(Exp[I Pi m]ModularJ4[tau]+4(Exp[2Pi I m]+1))];
an=-2 EllipticK[y]+(4 EllipticK[1-y] Log[4])/\[Pi];
anp1=(-2 \[Pi] EllipticE[y]+(\[Pi]-\[Pi] y) EllipticK[y]+2 (-2 EllipticE[1-y]+(1+y) EllipticK[1-y]) Log[4])/(12 \[Pi] (-1+y)^2);
om3=an+anp1/z;
Do[
anp2=(-(1+2 n)^3 an+32 (1+n)^2 (3+2 n) (1+y) anp1)/(64 (2+n) (5+7 n+2 n^2) (-1+y)^2);
om3+=anp2/z^(n+2);an=anp1; anp1=anp2;,{n,0,NMax}];om3 Sqrt[I(tau/(1-2 tau)-tb/(1-2 tb))]/Sqrt[I(tau/(1-2 tau)-tb/(1-2 tb))z]];




(* new definitions from Bruno, Aug 14, 2024 :*)

EichlerTTDLV[tau_,m_]:=(EichlerTTDLV[ttau_,m]:=#[ttau];#[tau])&[
  Compile[{{tttau,_Complex}}, With[{#1=Exp[(2I Pi)tttau+#2]}, {tttau-m/2+#3[[1]], tttau^2/2-m^2/8+1/12+#3[[1]]tttau-#3[[2]]}]]&[
    qt, 2. Pi Im[tauB[m]], HornerForm/@Sum[qt^n FourierCn[m,n] (2Pi I n)^{-1,-2},{n,50,1,-1}]]];
EichlerTTDp[tau_,m_]:=(EichlerTTDp[ttau_,m]:=#[ttau];#[tau])&[
  Compile[{{tttau,_Complex}}, With[{#1=Exp[(2I Pi)(-1/(4tttau))+#2]}, {#4+#3[[2]]+#3[[1]]/tttau,#3[[1]]}]]&[
    qt, 2. Pi Im[-.25/(tauB[m]-Boole[Re[tauB[m]]>=1/2])], HornerForm/@Sum[qt^n FourierCpn[m,n]{-4/(Pi n),8I/(Pi n)^2},{n,50,1,-1}], I QuantumVolume[m]]];
EichlerTTDpp[tau_,m_]:=(EichlerTTDpp[ttau_,m]:=#[ttau];#[tau])&[
  Compile[{{tttau,_Complex}}, With[{#1=tttau/(1-2tttau)}, With[{#2=Exp[(2I Pi)#1+#3]}, {1,0}#4[[1]]+{2,1}(#5+#1 #4[[1]]+#4[[2]])]]]&[
    taupp, qt, 2. Pi Im[tauB[m]/(1-2tauB[m])], HornerForm/@Sum[qt^n FourierCppn[m,n]{-2I/(Pi n),1/(Pi n)^2},{n,50,1,-1}], I QuantumVolumep[m]/2]];

EichlerTTDLVcut[tau_,m_]:=If[Im[#]>Abs[Re[#-1/2]]||Re[1/(1/2-tau)-1/(1/2-#)]Re[1/2-#]<0&[tauB[m]],
  EichlerTTDLV[tau,m], {1,(1+m)/2}-EichlerTTDLV[tau,m]];
EichlerTTDpcut[tau_,m_]:=If[Abs[#-(1+I)/4]<Sqrt[1/8]||Re[1/(1/2-tau)-1/(1/2-#)]>0&[tauB[m]],
  EichlerTTDp[tau,m], {1,(1+m)/2}-EichlerTTDp[tau,m]];
EichlerTTDppcut[tau_,m_]:=If[(2Sign[Im[m]]+Sign[Re[m]])Re[1/(1/2-tau)-1/(1/2-#)]>0&[tauB[m]],
  EichlerTTDpp[tau,m]-If[Im[m]<0&&Re[m]<=0,{2m,m},{0,0}], {1,(1+m)/2}-If[Im[m]<0&&Re[m]>0,{2m,m},{0,0}]-EichlerTTDpp[tau,m]];
EichlerTTDpppcut[tau_,m_]:=If[Abs[#-(3+I)/4]>Sqrt[1/8]&&Re[1/(1/2-tau)-1/(1/2-#)]>0&[tauB[m]],
  {-#1,-#1-#2}&@@EichlerTTDp[tau-1,m],{1+#1,(1+m)/2+#1+#2}&@@EichlerTTDp[tau-1,m]]

(* new definitions from Bruno, Sep 1st, 2024 :*)
(**New auxiliary functions that give the whole vector {1,m,T,TD} (and a variant for the expansion near tau=1/2)**)
Module[{nmax=50,prec=60,tau,taupp,j4tau,j4tauB,prefactor,scaledC,tauscaledC,tb,ismneg},
  Eichler1mTTDB[m_]:=Eichler1mTTDB[m]=(
    tb=SetPrecision[tauB[m],prec];
    tau=Simplify[taupp/(1+2taupp)/.taupp->tb/(1-2tb)-I EichlerTBtauR[tb]eps];
    prefactor=Series[Simplify[D[tau,eps]]DedekindEta[tau]^4 DedekindEta[2tau]^6/DedekindEta[4tau]^4,{eps,0,nmax}];
    j4tau=Series[ModularJ4[tau],{eps,0,nmax}];
    j4tauB=SeriesCoefficient[j4tau,0];
    scaledC=prefactor (eps(1+(8+j4tauB)/(j4tau-j4tauB)))^(1/2);
    tauscaledC=tau scaledC;
    ismneg=Im[m]<0||Im[m]==0&&Re[m]<0;
    Compile[{{ttau,_Complex}},With[{#eps=#scale(ttau/(1-2ttau)-#shift)},{1,m,#T,#TD}]]&[
      SetPrecision[#,$MachinePrecision]&[<|
          "eps"->eps,
          "scale"->I/EichlerTBtauR[tb],
          "shift"->tb/(1-2tb),
          "T"->HornerForm[Sum[Coefficient[scaledC,eps,i]eps^i/(i+1/2),{i,0,nmax-1}]]eps^(1/2)+.5-If[ismneg,m,0],
          "TD"->HornerForm[Sum[Coefficient[tauscaledC,eps,i]eps^i/(i+1/2),{i,0,nmax-1}]]eps^(1/2)+.25(1+If[ismneg,-m,m])|>]]);
];
Eichler1mTTDLV[m_]:=Eichler1mTTDLV[m]=
  Compile[{{tau,_Complex}}, With[{#qt=Exp[#scale tau+#shift]}, {1,m,tau-m/2+#F1, tau^2/2-m^2/8+1/12+#F1 tau-#F2}]]&[
  <|"qt"->qt, "scale"->N[2I Pi], "shift"->N[2 Pi Im[tauB[m]]],
  "F1"->HornerForm[Sum[qt^n FourierCn[m,n] (2Pi I n)^(-1),{n,50,1,-1}]],
  "F2"->HornerForm[Sum[qt^n FourierCn[m,n] (2Pi I n)^(-2),{n,50,1,-1}]]|>];
Eichler1mTTDp[m_]:=Eichler1mTTDp[m]=
  Compile[{{tau,_Complex}}, With[{#qt=Exp[#scale/tau+#shift]}, {1,m,#F2+#F1/tau,#F1}]]&[
  <|"qt"->qt, "scale"->-I Pi/2., "shift"->Im[-Pi/2./(tauB[m]-Boole[Re[tauB[m]]>=1/2])],
  "F1"->HornerForm[Sum[qt^n FourierCpn[m,n](-4/(Pi n)),{n,50,1,-1}]],
  "F2"->HornerForm[I QuantumVolume[m]+Sum[qt^n FourierCpn[m,n](8I/(Pi n)^2),{n,50,1,-1}]]|>];
Eichler1mF1TDpp[m_]:=Eichler1mF1TDpp[m]=
  Compile[{{tau,_Complex}}, With[{#taupp=tau/(1-2tau)}, With[{#qt=Exp[#scale #taupp+#shift]}, {1,m,#F1,#taupp #F1+#F2}]]]&[
  <|"taupp"->taupp, "qt"->qt, "scale"->2. I Pi, "shift"->2. Pi Im[tauB[m]/(1-2tauB[m])],
  "F1"->HornerForm[Sum[qt^n FourierCppn[m,n](-2I/(Pi n)),{n,50,1,-1}]],
  "F2"->HornerForm[I QuantumVolumep[m]/2+Sum[qt^n FourierCppn[m,n](1/(Pi n)^2),{n,50,1,-1}]]|>];
(* Eichler1mTTDpp[m_][tau_]:={{1,0,0,0},{0,1,0,0},{0,0,1,2},{0,0,0,1}} . Eichler1mF1TDpp[m][tau]; *)



EichlerTBtauR[tb_]:=EichlerTBtauR[tb]=
   With[{tbpp = (1-tb)/(1-2tb)},
     Min[{1, Abs[Which[Re[tbpp] < 1/4, 4(1-tb)^2/(6tb-5),
                      Re[tbpp] < 1/2, 1/(1+2tb),
                      Re[tbpp] < 3/4, 1/(3-2tb),
                      True, (4tb^2)/(6tb-1)]] / Abs[2tb-1]}]];
         
(* Probably obsolete: *)
Module[{nmax=50,prec=60,tau,taupp,j4tau,j4tauB,prefactor,scaledC,tb,sign,compiled},
EichlerTB[tau_,m_]:=Eichler1mTTDB[m][tau][[3]];
EichlerTDB[tau_,m_]:=Eichler1mTTDB[m][tau][[4]];
  EichlerTTDBcompile[m_]:=(tb=SetPrecision[tauB[m],prec]; ismneg=Im[m]<0||Im[m]==0&&Re[m]<0; tau=Simplify[taupp/(1+2taupp)/.taupp->tb/(1-2tb)-I EichlerTBtauR[tb]eps];
    prefactor=Series[Simplify[D[tau,eps]]DedekindEta[tau]^4 DedekindEta[2tau]^6/DedekindEta[4tau]^4,{eps,0,nmax}];
    j4tau=Series[ModularJ4[tau],{eps,0,nmax}]; j4tauB=SeriesCoefficient[j4tau,0]; scaledC=prefactor (eps(1+(8+j4tauB)/(j4tau-j4tauB)))^(1/2);
    compiled=(Compile[{{ttau,_Complex}},With[{eps=#1(ttau/(1-2ttau)-#2)},#3]]&@@SetPrecision[{I/EichlerTBtauR[tb],tb/(1-2tb),eps^(1/2) HornerForm[CoefficientList[#1 scaledC,eps] . (eps^(Range[nmax]-1)/(Range[nmax]-1/2))]+#2},$MachinePrecision])& @@@ {{1,.5-m If[ismneg,1,0]},{tau,.25(1+If[ismneg,-m,m])}};
    (EichlerTB[ttau_,m]:=#1[ttau];EichlerTDB[ttau_,m]:=#2[ttau];)&@@compiled; compiled);
];

SetAttributes[SwitchMin, HoldRest];
SwitchMin[l_List, code__] :=Hold[code][[Position[l,Min[l]][[1,1]]]];
    
(* EichlerTTDBcut[tau_,m_]:=If[False,
   {1-EichlerTB[tau,m],(1+m)/2-EichlerTDB[tau,m]},
   {EichlerTB[tau,m],EichlerTDB[tau,m]}]; *)
   
EichlerTTDBcut[tau_,m_]:={EichlerTB[tau,m],EichlerTDB[tau,m]};

ExpansionParametersLogA[m_]:=ExpansionParametersLogA[m]=Im[tauB[m]];
ExpansionParametersLogB[m_]:=ExpansionParametersLogB[m]=With[{tb=tauB[m]},Im[tb]/(4(Min[Re[tb]^2,(1-Re[tb])^2]+Im[tb]^2))];
ExpansionParametersLogC[m_]:=ExpansionParametersLogC[m]=1/(2-4tauB[m]);
ExpansionParametersLogD[m_]:=ExpansionParametersLogD[m]=-Log[EichlerTBtauR[tauB[m]]]/(2Pi);
ExpansionParametersLog=Compile[{{tau,_Complex},{m,_Complex}},
  With[{A=ExpansionParametersLogA[m], B=ExpansionParametersLogB[m], C=ExpansionParametersLogC[m], D=ExpansionParametersLogD[m]},
    With[{tauppshift=C+(4tau-2)^-1}, {A-Im[tau], B+.25Im[tau^-1], Im[tauppshift], B+.25Im[(tau-1)^-1], Log[Abs[tauppshift]]/(2.Pi)+D}]]];
      
EichlerZ[{r_,d1_,d2_,ch2_},tau_?NumericQ/;Im[tau]>0,m_?NumericQ/;-1<Re[m]<1]:=EichlerZaux[m][tau,r,d1,d2,ch2];
EichlerZaux[0]=EichlerZaux[0.]=EichlerZaux[0. I]=
  Compile[{{tau,_Complex},{r,_Real},{d1,_Real},{d2,_Real},{ch2,_Real}},
    Module[{tM=Table[0I,4,5], ttau,v,pLV,p0,p1},
      tM=#ToFundDomainAux[tau];
      ttau=tM[[1,1]];
      v={0,-ch2,d2,d1+d2,-2r} . tM;
      pLV=-Im[ttau];
      p0=+.25Im[ttau^-1];
      p1=.25Im[(ttau-1)^-1];
      Switch[Min[{pLV,p0,p1}],
        pLV, v . #ELV[ttau],
        p0, v . #Ep[ttau],
        p1, v . {{1,0,0,0},{0,1,0,0},{1,0,1,0},{1/2,1/2,1,1}} . #Ep[ttau-1],
             _, 0]]]&[
         <|
         "ToFundDomainAux"->ToFundDomainAux,
         "ELV"->Eichler1mTTDLV[0],
         "Ep"->Eichler1mTTDp[0]
         |>];
EichlerZaux[m_]:=EichlerZaux[m]=
  Compile[{{tau,_Complex},{r,_Real},{d1,_Real},{d2,_Real},{ch2,_Real}},
    Module[{tM=Table[0I,4,5], ttau,v,pLV,p0,p12,p1,pB,pB1,pB2,pB3,pB4,sign,tauppshift,M1shiftM,cutM1,cutM2,cutM},
      tM=#ToFundDomainAux[tau];
      ttau=tM[[1,1]];
      v={0,-ch2,d2,d1+d2,-2r} . tM;(*The entry 0 here kills the first row {ttau,0,0,0} of tM.*)
      pLV=#shiftLV-Im[ttau];
      p0=#shiftp+.25Im[ttau^-1];
      tauppshift=#shiftpp+(4ttau-2)^-1;
      p12=Im[tauppshift];
      p1=#shiftp+.25Im[(ttau-1)^-1];
      pB=Log[#tiny+Abs[tauppshift]]/(2.Pi)+#shiftB;
      pB1=Log[#tiny+Abs[#shiftpp+(4ttau+2)^-1]]/(2.Pi)+#shiftB;
      pB2=Log[#tiny+Abs[#shiftpp+(4ttau-6)^-1]]/(2.Pi)+#shiftB;
      pB3=Log[#tiny+Abs[tauppshift-1]]/(2.Pi)+#shiftB;
      pB4=Log[#tiny+Abs[tauppshift+1]]/(2.Pi)+#shiftB;
      sign=Sign[Re[1/(1/2-ttau)-1/(1/2-#tb)]];
      M1shiftM={{1,0,0,0},{0,1,0,0},{0,0,1,2},{0,0,0,1}};
      cutM1={{1,0,0,0},{0,1,0,0},{1,0,-1,0},{1/2,1/2,0,-1}};
      cutM2={{1,0,0,0},{0,1,0,0},{1,-2,-1,0},{1/2,-1/2,0,-1}};
      cutM=If[Im[m]>0||Im[m]==0&&Re[m]>0,cutM1,cutM2];
      Switch[Min[{pLV,p0,p12,p1,pB,pB1,pB2,pB3,pB4}],
        pLV, If[sign Re[#tb-1/2]>=0||Im[#tb]>=Abs[Re[#tb-1/2]], v, v . cutM] . #ELV[ttau],
        p0, If[sign>=0||Abs[#tb-(1+I)/4]<=Sqrt[1/8], v, v . cutM] . #Ep[ttau],
        p12, If[Im[m]<=0&&Re[m]<=0, If[sign<0,v . cutM2,v] . cutM1, If[Im[m]<0&&sign<=0||Im[m]>=0&&sign>=0, v, v . cutM]] . M1shiftM . #Epp[ttau],
        p1, If[sign<=0||Abs[#tb-(3+I)/4]<=Sqrt[1/8], v, v . cutM] . {{1,0,0,0},{0,1,0,0},{1,0,1,0},{1/2,1/2,1,1}} . #Ep[ttau-1],
        pB, v . #EB[ttau],
        pB1, If[sign>=0||4(1-Re[#tb])Im[#tb]^2-(1+Re[#tb])(1-2Re[#tb])^2>=0,v,v . cutM] . {{1,0,0,0},{0,1,0,0},{-1,0,1,0},{1/2,-1/2,-1,1}} . #EB[ttau+1],
        pB2, If[sign<=0||4 Re[#tb]Im[#tb]^2-(2-Re[#tb])(1-2Re[#tb])^2>=0,v,v . cutM] . {{1,0,0,0},{0,1,0,0},{1,0,1,0},{1/2,1/2,1,1}} . #EB[ttau-1],
        pB3, v . {{1,0,0,0},{0,1,0,0},{1,0,1,-4},{1/2,1/2,1,-3}} . #EB[(ttau-1)/(4ttau-3)],
        pB4, v . {{1,0,0,0},{0,1,0,0},{1,-2,-3,4},{1/2,-1/2,-1,1}} . #EB[(3ttau-1)/(4ttau-1)],
        _, 0]]]&[
         <|
         "tiny"->10^-300.,
         "tb"->tauB[m],
         "shiftLV"->Im[tauB[m]],
         "shiftp"->With[{tb=tauB[m]},Im[tb]/(4(Min[Re[tb]^2,(1-Re[tb])^2]+Im[tb]^2))],
         "shiftpp"->1/(2-4tauB[m]),
         "shiftB"->-Log[EichlerTBtauR[tauB[m]]]/(2Pi),
         "ToFundDomainAux"->ToFundDomainAux,
         "ELV"->Eichler1mTTDLV[m],
         "Ep"->Eichler1mTTDp[m],
         "Epp"->Eichler1mF1TDpp[m],
         "EB"->Eichler1mTTDB[m]
         |>];


Module[{ttau,logq,M,tb},
   EichlerZ0[{r_,d1_,d2_,ch2_},tau_?NumericQ/;Im[tau]>0,m_?NumericQ]:=(
     {ttau,M} = ToFundDomainApprox[tau];
     tb = tauB[m];
     logq = With[{sb = Im[tb]/(4(Min[Re[tb]^2,(1-Re[tb])^2]+Im[tb]^2))},
       {Im[-ttau+tb],
        Im[1/ttau]/4+sb,
        Im[1/(ttau-1/2)-1/(tb-1/2)]/4,
        Im[1/(ttau-1)]/4+sb}];
     {-ch2,d2,d1+d2,-2r} . M .
       Join[{1,m}, SwitchMin[logq,
         EichlerTTDLVcut[ttau,m],
         EichlerTTDpcut[ttau,m],
         EichlerTTDppcut[ttau,m],
         EichlerTTDpppcut[ttau,m]]]
   )];

EichlerZp[{r_,d1_,d2_,ch2_},taup_,m_]:=-2r EichlerTDp[-1/(4taup),m]+(d1+d2)EichlerTp[-1/4/taup,m]+d2 m-ch2;
EichlerZpp[{r_,d1_,d2_,ch2_},taupp_,m_]:=-2r EichlerTDpp[taupp/(1+2 taupp),m]+(d1+d2)EichlerTpp[taupp/(1+2 taupp),m]+d2 m-ch2;   
         
(* Boris: EichlerZ[{r_,d1_,d2_,ch2_},tau_/;Im[tau]>0,m_]:=Module[{ttau,M,tb},
{ttau,M}=ToFundDomainApprox[tau];tb=tauB[m];
(*Print[{ttau,tb}]; If[Im[ttau]>Im[tb],Print["LV"]]; *)
Which[Im[ttau]>Im[tb],
Tr[{{-ch2,d2,d1+d2,-2r}} . M . {{1},{m},{EichlerTLV[ttau,m]},{EichlerTDLV[ttau,m]}}],Re[ttau]<=1/2 && Im[-1/4/ttau]>Max[Im[-1/4/(tb-1)],Im[-1/4/(tb)]],
Tr[{{-ch2,d2,d1+d2,-2r}} . M . {{1},{m},{EichlerTp[ttau,m]},{EichlerTDp[ttau,m]}}],
Re[ttau]>1/2 && Im[-1/4/(ttau-1)]>Max[Im[-1/4/(tb-1)],Im[-1/4/(tb)]],
Tr[{{-ch2,d2,d1+d2,-2r}} . M . {{1},{m},{EichlerTp[ttau-1,m]+1},{EichlerTDp[ttau-1,m]+EichlerTp[ttau-1,m]+(m+1)/2}}],
True,
Tr[{{-ch2,d2,d1+d2,-2r}} . M . {{1},{m},{EichlerTpp[ttau,m]},{EichlerTDpp[ttau,m]}}]]];
*)

(* use classical version for testing: 
EichlerZ[{r_,d1_,d2_,ch2_},tau_,\[Lambda]_]:=-r tau(tau+m))+(d1+d2)tau+d2 m-ch2; *)

EichlerT[tau_,m_]:=EichlerZ[{0,1,0,0},tau,m];
EichlerTD[tau_,m_]:=EichlerZ[{-1/2,0,0,0},tau,m];
EichlerT0[tau_,m_]:=EichlerZ0[{0,1,0,0},tau,m];
EichlerTD0[tau_,m_]:=EichlerZ0[{-1/2,0,0,0},tau,m];

EichlerInt[tau0_,tau1_,m_]:=Module[{taupath,T,TD,s,so,eps,eps1},
taupath[s_]:=s tau1+(1-s)tau0;
so=NDSolve[{
eps[0]==1,eps1[0]==1,T[0]==0,TD[0]==0,
T'[s]==eps1[s] DtauT[taupath[s],m]taupath'[s],
TD'[s]==eps1[s]DtauT[taupath[s],m]taupath[s]taupath'[s],
WhenEvent[eps[s]Evaluate[Im[(ModularJ4[taupath[s]]+8)/(ModularJ4[taupath[s]]+4Exp[I Pi m]+4Exp[-I Pi m])]]>0 && eps[s]Re[(ModularJ4[taupath[s]]+8)/(ModularJ4[taupath[s]]+4Exp[I Pi m]+4Exp[-I Pi m])]<0,{
eps[s]->-eps[s],eps1[s]->eps1[s]Sign[Re[(ModularJ4[taupath[s]]+8)/(ModularJ4[taupath[s]]+4Exp[I Pi m]+4Exp[-I Pi m])]]}]
},
{T[s],TD[s],eps[s],eps1[s]},{s,0,1},DiscreteVariables->{Element[eps,{-1,1}],Element[eps1,{-1,1}]}];
{T[s],TD[s]}/.so[[1]]/.s->1
];

MeijerT0[tau_]:=(-1/(Pi)MeijerG[{{1/2,1/2},{1}},{{0,0},{0}},-16/(ModularJ4[tau]+8)]-I Pi)/(2Pi I);
MeijerTD0[tau_]:=(-1/(Pi)MeijerG[{{1/2,1/2,1},{}},{{0,0},{0}},16/(ModularJ4[tau]+8)]-2Pi^2/3)/(2Pi I)^2+1/12;
MeijerT[tau_]:=-1/(2Sqrt[2]Pi)MeijerG[{{1/4,3/4},{1}},{{0,0},{0}},-64/ModularJ4[tau]^2]/(2Pi I);
MeijerTD[tau_]:=(-1/(4 Sqrt[2] \[Pi])MeijerG[{{1/4,3/4,1},{}},{{0,0},{0}},64/ModularJ4[tau]^2]-Pi^2/6)/(2Pi I)^2+1/12;



(* Scattering diagram in Pi-stability slice *)




DtauT[tau_,m_]:=DedekindEta[tau]^4 DedekindEta[2tau]^6 /DedekindEta[4tau]^4 Sqrt[(ModularJ4[tau]+8)/(ModularJ4[tau]+4 Exp[I Pi m]+4 Exp[-I Pi m])];
DtauZ[{r_,d1_,d2_,ch2_},tau_,m_]:=(d1+d2-2r tau)DtauT[tau,m];
ArgDtauT[tau_,m_]:=Mod[4Arg[DedekindEta[tau]]+6Arg[DedekindEta[2tau]]-4Arg[DedekindEta[4tau]]+1/2 Arg[ModularJ4[tau]+8]-1/2Arg[ModularJ4[tau]+4 Exp[I Pi m]+4 Exp[-I Pi m]],2Pi,-Pi];
ArgDtauTD[tau_,m_]:=Mod[Arg[tau]+4Arg[DedekindEta[tau]]+6Arg[DedekindEta[2tau]]-4Arg[DedekindEta[4tau]]+1/2 Arg[ModularJ4[tau]+8]-1/2Arg[ModularJ4[tau]+4 Exp[I Pi m]+4 Exp[-I Pi m]],2Pi,-Pi];
UnitDtauT[tau_,m_]:=Normalize[DtauT[tau,m]];
UnitDtauTD[tau_,m_]:=Normalize[tau]Normalize[DtauT[tau,m]];
UnitDtauZ[{r_,d1_,d2_,ch2_},tau_,m_]:=Normalize[-2r tau+d1+d2]UnitDtauT[tau,m];

NormalizeFunctionDomain[fun_InterpolatingFunction]:=Function@@{fun[Rescale[#,{0,1},fun["Domain"][[1]]]]};
NormalizeApprox[z_,eps_:0.001]:=z/(eps+Abs[z]);

Module[{solution,tangentfunction,a,taunum,boundaries,eventstop},
IntegralCurve[tauinit_,tangent_,{ainit_,amin_,amax_},boundaries_List:{Im[tau]==0.01}]:=(eventstop=WhenEvent@@{boundaries/.tau->tau[a],"StopIntegration"};
tangentfunction[taunum_?NumericQ]=(If[Im[taunum]<=0,0,#]&[tangent/.tau->taunum]);
solution=tau/.Identity@@NDSolve[{tau[ainit]==tauinit,tau'[a]==tangentfunction[tau[a]],eventstop},tau,{a,amin,amax}];
NormalizeFunctionDomain[solution]);];

Module[{theta,tau1,tauinit,amin=-5,amax=5},
(* need to update argument tauexpr and start_List *)
RayCh[psi_,m_]:=RayCh[psi,m]=RayGeneral[{1,0,0,0},psi,m,1/((3(\[Pi]-theta))/(2\[Pi])-2.I),{theta,psi-Pi/2}];

RayGeneral[gamma:{_,_,_,_},psi_,m_,tauexpr_,start_List]:=(
tauinit=tauexpr/.FindRoot[Re[Exp[-I psi]EichlerZ[gamma,tauexpr,m]]==0,start];
IntegralCurve[tauinit,Normalize[EichlerZ[gamma,tau,m]]Conjugate[UnitDtauZ[gamma,tau,m]],{0,amin,amax},{Im[tau]==0.01(*,Norm[EichlerZch2[gamma,tau]]==10^-14*)}]);];

IntersectExactRaysLV[{r_,d1_,d2_,ch2_},{rr_,dd1_,dd2_,cch2_},psi_,m_]:=
(* returns (tau1,tau2) coordinate of intersection point of two rays, or {} if they are collinear *)
Module[{tau1,tau2,Exacttau1,Exacttau2},
If[(r (dd1+dd2)-rr (d1+d2)!=0) ,
If[TestBranch[{r,d1,d2,ch2},(cch2 r-dd2 m r-ch2 rr+d2 m rr)/(dd1 r+dd2 r-(d1+d2) rr),m]&&TestBranch[{rr,dd1,dd2,cch2},(cch2 r-dd2 m r-ch2 rr+d2 m rr)/(dd1 r+dd2 r-(d1+d2) rr),m],
            If[Length[Initialtau1]==0||Undefined[Initialtau1]||Length[Initialtau2]==0||Undefined[Initialtau2],{Initialtau1,Initialtau2}={1/(dd1 r+dd2 r-(d1+d2) rr) (cch2 r-dd2 m r-ch2 rr+d2 m rr-1/(2 Abs[dd1 r+dd2 r-(d1+d2) rr]) (dd1 r+dd2 r-(d1+d2) rr) \[Sqrt](4 (cch2-dd2 m)^2 r^2+2 (cch2-dd2 m) (-4 (ch2-d2 m) r rr+2 (d1+d2-m r)^2 rr-2 r (d1+d2-m r) (dd1+dd2-m rr))+(ch2-d2 m) (4 (ch2-d2 m) rr^2-4 (d1+d2-m r) rr (dd1+dd2-m rr)+4 r (dd1+dd2-m rr)^2)) Sin[psi]),1/(2 Abs[dd1 r+dd2 r-(d1+d2) rr]) \[Sqrt](4 (cch2-dd2 m)^2 r^2+2 (cch2-dd2 m) (-4 (ch2-d2 m) r rr+2 (d1+d2-m r)^2 rr-2 r (d1+d2-m r) (dd1+dd2-m rr))+(ch2-d2 m) (4 (ch2-d2 m) rr^2-4 (d1+d2-m r) rr (dd1+dd2-m rr)+4 r (dd1+dd2-m rr)^2)) Cos[psi]}];
{Exacttau1,Exacttau2}={tau1,tau2}/.FindRoot[{Re[Exp[-I psi](-r EichlerTDLV[tau1+I tau2,m]+(d1+d2) EichlerTLV[tau1+I tau2,m]+d2 m-ch2)],Re[Exp[-I psi](-rr EichlerTDLV[tau1+I tau2,m]+(dd1+dd2) EichlerTLV[tau1+I tau2,m]+dd2 m-cch2)]},{tau1,Initialtau1},{tau2,Initialtau2}];
{Initialtau1,Initialtau2}={Exacttau1,Exacttau2};If[Im[Exp[-I psi](-r EichlerTDLV[Exacttau1+I Exacttau2,m]+(d1+d2) EichlerTLV[Exacttau1+I Exacttau2,m]+d2 m-ch2)]>0&&Im[Exp[-I psi](-rr EichlerTDLV[Exacttau1+I Exacttau2,m]+(dd1+dd2) EichlerTLV[Exacttau1+I Exacttau2,m]+dd2 m-cch2)]>0,Exacttau1+I Exacttau2,0]
,0]]];

IntersectExactRaysC[{r_,d1_,d2_,ch2_},{rr_,dd1_,dd2_,cch2_},psi_,m_]:=
(* returns (tau1,tau2) coordinate of intersection point of two rays, or {} if they are collinear *)
Module[{tau1,tau2,Exacttau1,Exacttau2},
If[(r (dd1+dd2)-rr (d1+d2)!=0) ,
If[TestBranch[{r,d1,d2,ch2},(cch2 r-dd2 m r-ch2 rr+d2 m rr)/(dd1 r+dd2 r-(d1+d2) rr),m]&&TestBranch[{rr,dd1,dd2,cch2},(cch2 r-dd2 m r-ch2 rr+d2 m rr)/(dd1 r+dd2 r-(d1+d2) rr),m],
            If[Length[Initialtau1]==0||Undefined[Initialtau1]||Length[Initialtau2]==0||Undefined[Initialtau2],{Initialtau1,Initialtau2}={1/(dd1 r+dd2 r-(d1+d2) rr) (cch2 r-dd2 m r-ch2 rr+d2 m rr-1/(2 Abs[dd1 r+dd2 r-(d1+d2) rr]) (dd1 r+dd2 r-(d1+d2) rr) \[Sqrt](4 (cch2-dd2 m)^2 r^2+2 (cch2-dd2 m) (-4 (ch2-d2 m) r rr+2 (d1+d2-m r)^2 rr-2 r (d1+d2-m r) (dd1+dd2-m rr))+(ch2-d2 m) (4 (ch2-d2 m) rr^2-4 (d1+d2-m r) rr (dd1+dd2-m rr)+4 r (dd1+dd2-m rr)^2)) Sin[psi]),1/(2 Abs[dd1 r+dd2 r-(d1+d2) rr]) \[Sqrt](4 (cch2-dd2 m)^2 r^2+2 (cch2-dd2 m) (-4 (ch2-d2 m) r rr+2 (d1+d2-m r)^2 rr-2 r (d1+d2-m r) (dd1+dd2-m rr))+(ch2-d2 m) (4 (ch2-d2 m) rr^2-4 (d1+d2-m r) rr (dd1+dd2-m rr)+4 r (dd1+dd2-m rr)^2)) Cos[psi]}];
{Exacttau1,Exacttau2}={tau1,tau2}/.FindRoot[{Re[Exp[-I psi](-r EichlerTDp[tau1+I tau2,m]+(d1+d2) EichlerTp[tau1+I tau2,m]+d2 m-ch2)],Re[Exp[-I psi](-rr EichlerTDp[tau1+I tau2,m]+(dd1+dd2) EichlerTp[tau1+I tau2,m]+dd2 m-cch2)]},{tau1,Initialtau1},{tau2,Initialtau2}];
{Initialtau1,Initialtau2}={Exacttau1,Exacttau2};If[Im[Exp[-I psi](-r EichlerTDp[Exacttau1+I Exacttau2,m]+(d1+d2) EichlerTp[Exacttau1+I Exacttau2,m]+d2 m-ch2)]>0&&Im[Exp[-I psi](-rr EichlerTDp[Exacttau1+I Exacttau2,m]+(dd1+dd2) EichlerTp[Exacttau1+I Exacttau2,m]+dd2 m-cch2)]>0,Exacttau1+I Exacttau2,0]
,0]]];

CriticalPsi[m_]:=N[{{ArcTan[Im[QuantumVolume[m]-I m]/Re[QuantumVolume[m]-I m]],ArcTan[Im[QuantumVolume[m]]/Re[QuantumVolume[m]]]},{ArcTan[Im[QuantumVolumep[m]]/Re[QuantumVolumep[m]]],ArcTan[Im[QuantumVolumep[m]-I(m-1)]/Re[QuantumVolumep[m]-I(m-1)]]},If[Im[m]==0,Infinity,-ArcTan[Re[m]/Im[m]]]}];
XY[tau_,psi_,m_]:={Re[Exp[-I psi]EichlerT[tau,m]],-Re[Exp[-I psi]EichlerTD[tau,m]]}/Cos[psi];



EvaluateKronecker[f_]:=
f/.Subscript[Kr, kappa_][g1_,g2_]:>Simplify[F0HiggsBranchFormula[{{0,kappa},{-kappa,0}},{1/g1,-1/g2},{g1,g2}]];

(* taken in part from CoulombHiggs.m package *) 
F0HiggsBranchFormula[Mat_,Cvec_,Nvec_]:=Module[{Cvec0},
  If[Max[Nvec]<0,Print["F0HiggsBranchFormula: The dimension vector must be positive !"]];
  If[Plus@@Nvec==0,Return[0]];
  If[Plus@@Nvec==1,Return[1]];
  Cvec0=Cvec-(Plus@@(Nvec Cvec))/(Plus@@Nvec);
DivisorSum[GCD@@Nvec,(y-1/y)/(y^#-1/y^#)/# MoebiusMu[#] F0RationalInvariant[Mat,Cvec0,Nvec/#,y^#]&]];
F0RationalInvariant[Mat_,Cvec_,Nvec_,y_]:=Module[{Li,gcd},
	gcd=GCD@@Nvec;
	Li=Flatten[Map[Permutations,F0ListAllPartitions[{gcd}]],1];
	Sum[
	   Product[F0StackInvariant[Mat,Cvec,Nvec Li[[i,j,1]]/gcd,y],{j,Length[Li[[i]]]}]/Length[Li[[i]]]/(y-1/y)^(Length[Li[[i]]]-1),
	{i,Length[Li]}]];

F0QDeformedFactorial[n_,y_]:=If[n<0,Print["F0QDeformedFactorial[n,y] is defined only for n>=0"],
		If[n==0,1,(y^(2n)-1)/(y^2-1)F0QDeformedFactorial[n-1,y]]];
		
F0StackInvariant[Mat_,Cvec_,Nvec_,y_]:=Module[{m,JKListAllPermutations,pa,Cvec0},
  m=Length[Nvec];
  If[Max[Nvec]<0,Print["F0StackInvariant: The dimension vector must be positive !"]];
  If[Plus@@Nvec==0,Return[0]];
  If[Plus@@Nvec==1,Return[1]];
  Cvec0=Cvec-(Plus@@(Nvec Cvec))/(Plus@@Nvec);
  pa=Flatten[Map[Permutations,F0ListAllPartitions[Nvec]],1];
    (-y)^( Sum[-Max[Mat[[k,l]],0]Nvec[[k]]Nvec[[l]],{k,m},{l,m}]-1+Plus@@ Nvec)
	   (y^2-1)^(1-Plus@@Nvec)
	Sum[If[(Length[pa[[i]]]==1) ||And@@Table[Sum[Cvec0[[k]] pa[[i,a,k]],{a,b},{k,m}]>0,{b,Length[pa[[i]]]-1}],
      (-1)^(Length[pa[[i]]]-1)
       y^(2 Sum[Max[ Mat[[l,k]],0] pa[[i,a,k]]pa[[i,b,l]],
    {a,1,Length[pa[[i]]]},{b,a,Length[pa[[i]]]},{k,m},{l,m}])/
    Product[F0QDeformedFactorial[pa[[i,j,k]],y] ,{j,1,Length[pa[[i]]]},{k,m}],0],{i,Length[pa]}]
];
	
F0ListAllPartitions[gam_]:=Module[{m,unit,Li},
If[Plus@@gam==1, {{gam}},
		m=Max[Select[Range[Length[gam]],gam[[#]]>0&]];
        unit=Table[If[i==m,1,0],{i,Length[gam]}];        
	    Li=F0ListAllPartitions[gam-unit];
        Union[Map[Sort,
        Union[Flatten[
				Table[Union[Flatten[{{Flatten[{Li[[i]],{unit}},1]},
					    Table[
						  Table[If[j==k,Li[[i,j]]+unit,Li[[i,j]]]
						  ,{j,Length[Li[[i]]]}]
	                    ,{k,Length[Li[[i]]]}]}
                      ,1]]
				,{i,Length[Li]}],1]]
         ,1]]
	]];

F0BinarySplits[Nvec_]:=Module[{Li,Li1,Nl},
If[Plus@@Nvec==1,Li1=Nvec,
Li=Drop[Drop[Flatten[Table[Table[Nl[i],{i,Length[Nvec]}],Evaluate[Sequence@@Table[{Nl[i],0,Nvec[[i]]},{i,Length[Nvec]}]]],Length[Nvec]-1],1],-1];
Li1=Take[Li,Ceiling[Length[Li]/2]];
Li1]];



(* ::Section:: *)
(*Basic operations on trees and lists of trees*)


TreeCharge[arg : {r_, d1_, d2_, ch2_} /; FreeQ[arg, Ch]] := {r, d1, d2, ch2};
TreeCharge[trees_List] := Total[TreeCharge /@ trees];
TreeCharge[arg : Except[_List]] := arg /. repCh;

TreeConstituents[Tree_]:=(* gives the flattened list of constituents of a tree *)
If[!ListQ[Tree]||Length[Tree]>2,{Tree},
Join[TreeConstituents[Tree[[1]]],TreeConstituents[Tree[[2]]]]
];

FlipTree[Tree_]:=If[!ListQ[Tree],Tree/.Ch[m1_,m2_][1]:>Ch[-m1,-m2]/.Ch[m1_,m2_]:>Ch[-m1,-m2][1],If[Length[Tree]==4,{-Tree[[1]],Tree[[2]],Tree[[3]],-Tree[[4]]},{FlipTree[Tree[[2]]],FlipTree[Tree[[1]]]}]];

ShiftTree[{k1_,k2_},Tree_]:=Tree/.{Ch[m1_,m2_][1]:>Ch[m1+k1,m2+k2][1],Ch[m1_,m2_]:>Ch[m1+k1,m2+k2],GV[m1_,m2_,ch2_]:>GV[m1,m2,ch2+m2 k1+m1 k2],{r_Integer,d1_Integer,d2_Integer,ch2_Integer}:>SpecFlow[{r,d1,d2,ch2},{k1,k2}]};
TreeHue[n_,i_]:=Blend[{Brown,Green},i/n];



(* ::Section:: *)
(*Large volume scattering diagram*)


xytost[{x_,y_},m_:1/2]:={x,Sqrt[x^2+m x+2y]};
sttoxy[{s_,t_},m_:1/2]:={s,-1/2(s^2+m s-t^2)};

Rayxy[{r_,d1_,d2_,ch2_},xx_,L_,m_]:=Module[{x0,y0},
(* Ray starting from (xx,?) if r<>0, or (?,xx) if r=0, extending by distance L *)
If[r==0,{x0,y0}={(ch2 -m d2)/(d1+d2),xx},
{x0,y0}={xx,(ch2 -m d2-(d1+d2)xx)/(2 r)}];
{Point[{x0,y0}],Line[{{x0,y0},{x0-2r L,y0+L(d1+d2)}}]}
]

TreeConstituents[Tree_]:=(* gives the flattened list of constituents of a tree *)
If[!ListQ[Tree]||Length[Tree]>2,{Tree},
Join[TreeConstituents[Tree[[1]]],TreeConstituents[Tree[[2]]]]
];
TreeHue[n_,i_]:=Blend[{Brown,Green},i/n];

DiscF0[{r_,d1_,d2_,ch2_},m_:1/2]:=((d1-d2-m r)^2+4 d1 d2-4 r ch2)/(8 r^2);
TestBranch[{r_,d1_,d2_,ch2_},s_,m_:1/2]:=Module[{Di},
(* tests if (s,.) is on the branch with ImZ>0 *)
If[r==0,d1+d2-m r>0,
Di=Max[DiscF0[{r,d1,d2,ch2},m],0];
(s<=(d1+d2-m r)/(2r)-Sqrt[2Di] &&r>0)||(s>=(d1+d2-m r)/(2r)+Sqrt[2Di] &&r<0)
]];

InitialPosition[{r_,d1_,d2_,ch2_},m_:1/2]:=Module[{},
If[!IntegerQ[d1 d2-ch2],Print["Non integer second Chern class !"]];
If[r==0,(ch2-m d2)/(d1+d2-m r),If [r>0,(d1+d2-m r)/(2r)-Sqrt[2Max[DiscF0[{r,d1,d2,ch2},m],0]],(d1+d2-m r)/(2r)+Sqrt[2Max[DiscF0[{r,d1,d2,ch2},m],0]] ]]];

IntersectRays[{r_,d1_,d2_,ch2_},{rr_,dd1_,dd2_,cch2_},m_:1/2]:=
(* returns (x,y) coordinate of intersection point of two rays, or {} if they are collinear *)
If[(r (dd1+dd2)-rr (d1+d2)!=0) (*&&(r rr+d dd<=0)*),{-((-cch2 r+dd2 m r+ch2 rr-d2 m rr)/(dd1 r+dd2 r-d1 rr-d2 rr)),(-cch2 d1-cch2 d2+ch2 dd1+ch2 dd2-d2 dd1 m+d1 dd2 m)/(2 (dd1 r+dd2 r-d1 rr-d2 rr))},{}];

IntersectRays[{r_,d1_,d2_,ch2_},{rr_,dd1_,dd2_,cch2_},z_,zz_,m_:1/2]:=
(* returns (x,y) coordinate of intersection point of two rays, or {} if they don't intersect *)
Module[{zi},If[(r (dd1+dd2)-rr (d1+d2)!=0) ,zi={-((-cch2 r+dd2 m r+ch2 rr-d2 m rr)/(dd1 r+dd2 r-d1 rr-d2 rr)),(-cch2 d1-cch2 d2+ch2 dd1+ch2 dd2-d2 dd1 m+d1 dd2 m)/(2 (dd1 r+dd2 r-d1 rr-d2 rr))};
If[(zi-z) . {-2r,d1+d2-m r}>=0&&(zi-zz) . {-2rr,dd1+dd2-m rr}>=0,zi,{}],{}]];

Rays[{r_,d1_,d2_,ch2_},t_,psi_,m_:1/2]:=(* s(t) for the ray Re[E^(-I psi) Z_gam]\[Equal]0 *)
If[r!=0,(d1+d2-m r)/(2r)-t Tan[psi]-Sign[r]/Cos[psi]Sqrt[t^2+2DiscF0[{r,d1,d2,ch2},m]Cos[psi]^2],(ch2-m d2)/(d1+d2-m r)-t Tan[psi]];

IntersectRaysSt[{r_,d1_,d2_,ch2_},{rr_,dd1_,dd2_,cch2_},psi_,m_:1/2]:=
(* returns (s,t) coordinate of intersection point of two rays, or {} if they are collinear *)
If[(r (dd1+dd2)-rr( d1+d2)!=0) ,If[TestBranch[{r,d1,d2,ch2},(cch2 r-dd2 m r-ch2 rr+d2 m rr)/(dd1 r+dd2 r-(d1+d2) rr),m]&&TestBranch[{rr,dd1,dd2,cch2},(cch2 r-dd2 m r-ch2 rr+d2 m rr)/(dd1 r+dd2 r-(d1+d2) rr),m],{1/(dd1 r+dd2 r-(d1+d2) rr) (cch2 r-dd2 m r-ch2 rr+d2 m rr-1/(2 Abs[dd1 r+dd2 r-(d1+d2) rr]) (dd1 r+dd2 r-(d1+d2) rr) \[Sqrt](4 (cch2-dd2 m)^2 r^2+2 (cch2-dd2 m) (-4 (ch2-d2 m) r rr+2 (d1+d2-m r)^2 rr-2 r (d1+d2-m r) (dd1+dd2-m rr))+(ch2-d2 m) (4 (ch2-d2 m) rr^2-4 (d1+d2-m r) rr (dd1+dd2-m rr)+4 r (dd1+dd2-m rr)^2)) Sin[psi]),1/(2 Abs[dd1 r+dd2 r-(d1+d2) rr]) \[Sqrt](4 (cch2-dd2 m)^2 r^2+2 (cch2-dd2 m) (-4 (ch2-d2 m) r rr+2 (d1+d2-m r)^2 rr-2 r (d1+d2-m r) (dd1+dd2-m rr))+(ch2-d2 m) (4 (ch2-d2 m) rr^2-4 (d1+d2-m r) rr (dd1+dd2-m rr)+4 r (dd1+dd2-m rr)^2)) Cos[psi]},{}]];

WallCircle[{r_,d1_,d2_,ch2_},{rr_,dd1_,dd2_,cch2_},m_:1/2]:=Module[{R},
R=WallRadius[{r,d1,d2,ch2},{rr,dd1,dd2,cch2},m];
If[R>0,Circle[{(cch2 r-dd2 m r-ch2 rr+d2 m rr)/(dd1 r+dd2 r-(d1+d2) rr),0},Sqrt[R],{0,Pi}],{}]];

WallRadius[{r_,d1_,d2_,ch2_},{rr_,dd1_,dd2_,cch2_},m_:1/2]:=(cch2^2 r^2+cch2 r (-((d1+d2) (dd1+dd2))+(dd1-dd2) m r)+r ((dd1+dd2) (ch2 (dd1+dd2)-d2 dd1 m+d1 dd2 m)-dd1 dd2 m^2 r)-(d1+d2) (ch2 (dd1+dd2)-d2 dd1 m+d1 dd2 m) rr+m (-ch2 dd1+ch2 dd2+d2 dd1 m+d1 dd2 m) r rr+cch2 ((d1+d2)^2-2 ch2 r+(-d1+d2) m r) rr+(ch2+d1 m) (ch2-d2 m) rr^2)/((dd1+dd2) r-(d1+d2) rr)^2;

(* LV slice m=M1 + T M2 *)
ZLV[{r_,d1_,d2_,ch2_},{s_,t_},M1_,M2_]:=-r((1+M2)(s+I t)(s+I t)+M1(s+I t))+d1 (s+I t)+d2(s+I t+M1+M2(s+I t))-ch2;
Rays[{r_,d1_,d2_,ch2_},t_,psi_,M1_,M2_]:=(* s(t) for the ray Re[E^(-I psi) Z_gam]\[Equal]0 *)
If[r!=0,-t Tan[psi]+1/(2 (1+M2) r) (d1+d2+d2 M2-M1 r-Sec[psi] \[Sqrt](((d1+d2+d2 M2-M1 r) Cos[psi]-2 (1+M2) r t Sin[psi])^2+4 (1+M2) r Cos[psi] ((-ch2+d2 M1+(1+M2) r t^2) Cos[psi]+(d1+d2+d2 M2-M1 r) t Sin[psi]))),(ch2-d2 M1-(d1+d2+d2 M2) t Tan[psi])/(d1+d2+d2 M2)];
WallRadius[{r_,d1_,d2_,ch2_},{rr_,dd1_,dd2_,cch2_},M1_,M2_]:=(cch2^2 (1+M2) r^2+r ((dd1+dd2+dd2 M2) (-d2 dd1 M1+d1 dd2 M1+ch2 (dd1+dd2+dd2 M2))-dd1 dd2 M1^2 r)+cch2 r (-((d1+d2+d2 M2) (dd1+dd2+dd2 M2))+M1 (dd1-dd2 (1+M2)) r)-(d1+d2+d2 M2) (-d2 dd1 M1+d1 dd2 M1+ch2 (dd1+dd2+dd2 M2)) rr+M1 ((d2 dd1+d1 dd2) M1+ch2 (-dd1+dd2+dd2 M2)) r rr+cch2 ((d1+d2+d2 M2)^2-2 ch2 (1+M2) r+M1 (-d1+d2+d2 M2) r) rr+(ch2-d2 M1) (ch2+d1 M1+ch2 M2) rr^2)/((1+M2) ((dd1+dd2+dd2 M2) r-(d1+d2+d2 M2) rr)^2);
WallCircle[{r_,d1_,d2_,ch2_},{rr_,dd1_,dd2_,cch2_},M1_,M2_]:=Module[{R},
R=WallRadius[{r,d1,d2,ch2},{rr,dd1,dd2,cch2},M1,M2];
If[R>0,Circle[{(cch2 r-dd2 M1 r-ch2 rr+d2 M1 rr)/(dd1 r+dd2 (1+M2) r-(d1+d2+d2 M2) rr),0},Sqrt[R],{0,Pi}],{}]];
xytost[{x_,y_},M1_,M2_]:={x,Sqrt[M1 x+x^2+M2 x^2+2 y]/Sqrt[1+M2]};
sttoxy[{s_,t_},M1_,M2_]:={s,-1/2(1+M2)(s^2-t^2)-1/2 M1 s};
InitialPosition[{r_,d1_,d2_,ch2_},M1_,M2_]:=Module[{},
If[!IntegerQ[d1 d2-ch2],Print["Non integer second Chern class !"]];
If[r==0,(ch2-d2 M1)/(d1+d2+d2 M2),If [r>0,(d1+d2+d2 M2-M1 r-Sqrt[-4 (ch2-d2 M1) (1+M2) r+(d1+d2+d2 M2-M1 r)^2])/(2 (1+M2) r),(d1+d2+d2 M2-M1 r-Sqrt[-4 (ch2-d2 M1) (1+M2) r+(d1+d2+d2 M2-M1 r)^2])/(2 (1+M2) r) ]]];


ScattIndex[TreeList_]:=Table[
(* compute index for each tree in the list; do not trust the result if internal lines have non-primitive charges *)
Simplify[Times@@ScattIndexInternal[TreeList[[i]]][[2]]],{i,Length[TreeList]}];

ScattIndexInternal[Tree_]:=Module[{S1,S2,g1,g2,kappa,Li},
(* compute {total charge, list of Kronecker indices associated to each vertex *)
If[!ListQ[Tree]||Length[Tree]>2,{Tree,{1}},
S1=ScattIndexInternal[Tree[[1]]]/.repCh;
S2=ScattIndexInternal[Tree[[2]]]/.repCh;
g1=GCD@@S1[[1]];g2=GCD@@S2[[1]];
kappa=2Abs[(S1[[1,1]](S2[[1,2]]+S2[[1,3]])-(S1[[1,2]]+S1[[1,3]])S2[[1,1]])/g1/g2];
Li=Join[S1[[2]],S2[[2]]];
AppendTo[Li,Subscript[Kr, kappa][Min[g1,g2],Max[g1,g2]]];
If[GCD@@(S1[[1]]+S2[[1]])!=1,Print["Beware, non-primitive state"]];
{S1[[1]]+S2[[1]],Li}]];

(* more careful implementation taking care of internal non-primitive states *)
Options[ScattIndexImprovedInternal] = {"Debug"->False};

ScattIndexImproved[TreeList_, opt: OptionsPattern[]]:=Table[
	(* compute index for each tree in the list *)
	Simplify[FOmbToOm[Last@ScattIndexImprovedInternal[TreeList[[i]], opt][[2]]]],{i,Length[TreeList]}];

ScattIndexImprovedInternal[Tree_, opt: OptionsPattern[]]:=Module[{S1,S2,g1,g2,gFinal, kappa,Li, tem, repOmAttb, rrr},
(* compute {total charge, list of Kronecker indices associated to each vertex *)
	If[!ListQ[Tree]||Length[Tree]>2,{Tree,{Join[{1}, Table[(y-y^-1)/(j(y^j-y^-j)), {j, 2, GCD@@(Tree/.repCh)}]]}},
	If[OptionValue["Debug"], Print["Calling with args: ", Tree[[1]], "  |  ", Tree[[2]]]];
    S1=ScattIndexImprovedInternal[Tree[[1]], opt]/.repCh;
	S2=ScattIndexImprovedInternal[Tree[[2]], opt]/.repCh;
If[OptionValue["Debug"], Print["S1 is: ", S1, "   S2 is: ", S2]];
	g1=GCD@@S1[[1]];g2=GCD@@S2[[1]];
    gFinal = GCD@@(S1[[1]]+S2[[1]]);
	kappa=2Abs[(S1[[1,1]](S2[[1,2]]+S2[[1,3]])-(S1[[1,2]]+S1[[1,3]])S2[[1,1]])/g1/g2];
	Li=Join[S1[[2]],S2[[2]]];
If[OptionValue["Debug"], Print["Li is: ", Li, "  g1 is: ", g1, "  g2 is: ", g2, "  gFinal is: ", gFinal]];
AppendTo[Li,
repOmAttb = Join[
Table[CoulombHiggs`OmAttb[{P, 0}, y_]->Last[S1[[2]]][[P]], {P, 1, g1}],
Table[CoulombHiggs`OmAttb[{0, Q}, y_]->Last[S2[[2]]][[Q]], {Q, 1, g2}]
];
If[OptionValue["Debug"], Print["repOmAttb is: ", repOmAttb]];
tem = Table[
rrr = If[And@@(IntegerQ/@{P g1/gFinal, P g2/gFinal}),CoulombHiggs`FlowTreeFormulaRat[{{0, kappa}, {-kappa, 0}}, {g2, -g1}, {P g1/gFinal, P g2/gFinal}, y], 0];
Simplify[
rrr
/.repOmAttb
/.{CoulombHiggs`OmAttb[{p_, q_}, y___]:>0/;p>1||q>1||p q !=0}
],
{P, 1, gFinal}
];
If[OptionValue["Debug"], Print["tem is: ", tem]];tem
];
	(*If[GCD@@(S1[[1]]+S2[[1]])!=1,Print["Beware, non-primitive state"]];*)
	{S1[[1]]+S2[[1]],Li}]];

BinaryTreeIndex[TreeList_]:=Table[
(* compute index for each tree in the list; do not trust the result if internal lines have non-primitive charges *)
Simplify[Times@@BinaryTreeIndexInternal[TreeList[[i]]][[2]]],{i,Length[TreeList]}];

BinaryTreeIndexInternal[Tree_]:=Module[{S1,S2,kappa,Li},
If[!ListQ[Tree]||Length[Tree]>2,{Tree,{1}},
S1=BinaryTreeIndexInternal[Tree[[1]]]/.repCh;
S2=BinaryTreeIndexInternal[Tree[[2]]]/.repCh;
kappa=2Abs[(S1[[1,1]](S2[[1,2]]+S2[[1,3]])-(S1[[1,2]]+S1[[1,3]])S2[[1,1]])];
Li=Join[S1[[2]],S2[[2]]];
AppendTo[Li,(-1)^(kappa+1)(y^kappa-y^(-kappa))/(y-1/y)];
{S1[[1]]+S2[[1]],Li}]];

Module[{LiSlopes,LiHeights},
(* keep track of height of labels for each slope *)

InitialLabelPosition[m_]:=Module[{p},
p=Position[LiSlopes,m];
If[Length[p]==0,
AppendTo[LiSlopes,m];
AppendTo[LiHeights,0];{m,-.2},
LiHeights[[p[[1,1]]]]+=1;
{m,-.2-.28*LiHeights[[p[[1,1]]]]}]]; 

ScattDiagSt[TreeList_,m_:1/2]:=Module[{T,TNum,Diagram,xmin,xmax},
(* Draw overlayed diagrams in (s,t) plane for a LIST of trees *)
Diagram={};xmin={};xmax={};LiSlopes={};LiHeights={};
Do[T=ScattDiagInternalSt[TreeList[[i]],m];
TNum=T/.repCh;
AppendTo[Diagram,TreeHue[Length[TreeList],i]];
AppendTo[Diagram,Dashing[None]];
AppendTo[Diagram,TNum[[3]]];
AppendTo[Diagram,Dashed];
AppendTo[Diagram,Circle[{TNum[[2,1]],0},TNum[[2,2]],{0,Pi}]];
xmin=Min[xmin,TNum[[4,1]]];
xmax=Max[xmax,TNum[[4,2]]];
(*AppendTo[Diagram,Arrow[{T[[2]],T[[2]]+{0,1.5}}]]*)
,{i,Length[TreeList]}];
(*AppendTo[Diagram,Text[T[[1]],T[[2]]+{0,1}]]*);
AppendTo[Diagram,Dotted];
AppendTo[Diagram,Black];
(* AppendTo[Diagram,Table[{Line[{{k,0},{k+1/2,1/2},{k,1},{k-1/2,1/2},{k,0}}]},{k,-1+0Floor[xmin],-1+0Ceiling[xmax]}]]; *)
Graphics[Diagram]];

ScattDiagInternalSt[Tree_,m_:1/2]:=Module[{S1,S2,TreeNum,sInit,z,Li}(* construct total charge, coordinate of root and list of line segments in (s,t) coordinates, {min(x), max(x)} *), 
If[!ListQ[Tree]||Length[Tree]>2,
	TreeNum=Tree/.repCh;
sInit=InitialPosition[TreeNum,m];{Tree,{sInit,0},{Text[Tree/.repChO,InitialLabelPosition[sInit]]},{sInit,sInit}},
S1=ScattDiagInternalSt[Tree[[1]],m];
S2=ScattDiagInternalSt[Tree[[2]],m];
z=IntersectRays[S1[[1]]/.repCh,S2[[1]]/.repCh,m];
If[Length[z]==0,Print["Illegal tree"],
Li=S1[[3]];AppendTo[Li,S2[[3]]];AppendTo[Li,Arrow[{S1[[2]],xytost[z,m]}]];AppendTo[Li,Arrow[{S2[[2]],xytost[z,m]}]];
{S1[[1]]+S2[[1]],xytost[z,m],Li,{Min[S1[[4,1]],S2[[4,1]]],Max[S1[[4,2]],S2[[4,2]]]}}]]];

ScattDiagInternalLV[Tree_,psi_,Styl_,m_:1/2]:=Module[{S1,S2,TreeNum,sInit,z,Li},
(* construct total charge, coordinate of root and list of line segments in (s,t) coordinates, {min(x), max(x)} *) 
If[!ListQ[Tree]||Length[Tree]>2,
	TreeNum=Tree/.repCh;
sInit=InitialPosition[TreeNum,m];{Tree,{sInit,0},{Text[Tree/.repChO,InitialLabelPosition[sInit]]},{sInit,sInit}},
S1=ScattDiagInternalLV[Tree[[1]],psi,Styl,m];
S2=ScattDiagInternalLV[Tree[[2]],psi,Styl,m];
z=IntersectRaysSt[S1[[1]]/.repCh,S2[[1]]/.repCh,psi,m];
If[Length[z]==0,Print["Illegal tree"],
Li=S1[[3]];AppendTo[Li,S2[[3]]];
AppendTo[Li,If[(S1[[1]]/.repCh)[[1]]==0||S1[[2,2]]==z[[2]],Line[{S1[[2]],z}],ParametricPlot[{Rays[S1[[1]]/.repCh,t,psi,m],t},{t,S1[[2,2]],z[[2]]},PlotStyle->Styl][[1]]]];
AppendTo[Li,If[(S2[[1]]/.repCh)[[1]]==0||S2[[2,2]]==z[[2]],Line[{S2[[2]],z}],ParametricPlot[{Rays[S2[[1]]/.repCh,t,psi,m],t},{t,S2[[2,2]],z[[2]]},PlotStyle->Styl][[1]]]];
{S1[[1]]+S2[[1]],z,Li,{Min[S1[[4,1]],S2[[4,1]]],Max[S1[[4,2]],S2[[4,2]]]}}]]];

ScattDiagLV[TreeList_,psi_,m_:1/2]:=Module[{T,TNum,Diagram,xmin,xmax},
(* Draw overlayed diagrams in (s,t) plane for a LIST of trees *)
LiSlopes={};LiHeights={};
Diagram={AbsoluteThickness[1]};xmin={};xmax={};
Do[T=ScattDiagInternalLV[TreeList[[i]],psi,{TreeHue[Length[TreeList],i],AbsoluteThickness[1]},m];
TNum=T/.repCh;
AppendTo[Diagram,TreeHue[Length[TreeList],i]];
AppendTo[Diagram,Dashing[None]];
AppendTo[Diagram,TNum[[3]]];
 AppendTo[Diagram,Dashed];
AppendTo[Diagram,WallCircle[Plus@@TreeConstituents[TreeList[[i,1]]/.repCh],Plus@@TreeConstituents[TreeList[[i,2]]/.repCh],m]]; 
xmin=Min[xmin,TNum[[4,1]]];
xmax=Min[xmax,TNum[[4,2]]];AppendTo[Diagram,Dashing[None]];AppendTo[Diagram,Black];
AppendTo[Diagram,
If[(T[[1]]/.repCh)[[1]]==0, 
  Line[{T[[2]],T[[2]]+1.5{-Tan[psi],1}}],
  ParametricPlot[{Rays[T[[1]]/.repCh,t,psi,m],t},{t,T[[2,2]],T[[2,2]]+2},PlotStyle->{Black,AbsoluteThickness[1]}][[1]]]] 
,{i,Length[TreeList]}];
AppendTo[Diagram,Dotted];
AppendTo[Diagram,Black]; 
Graphics[Diagram]];

ScattDiag[TreeList_,psi_,m_:1/2]:=Module[{T,TNum,Diagram,xmin,xmax},
(* Draw overlayed diagrams in (s,t) plane for a LIST of trees *)
LiSlopes={};LiHeights={};
Diagram={AbsoluteThickness[1]};xmin={};xmax={};
Do[T=ScattDiagInternal[TreeList[[i]],psi,{TreeHue[Length[TreeList],i],AbsoluteThickness[1]},m];
TNum=T/.repCh;
AppendTo[Diagram,TreeHue[Length[TreeList],i]];
AppendTo[Diagram,Dashing[None]];
AppendTo[Diagram,TNum[[3]]];
 AppendTo[Diagram,Dashed];
AppendTo[Diagram,WallCircle[Plus@@TreeConstituents[TreeList[[i,1]]/.repCh],Plus@@TreeConstituents[TreeList[[i,2]]/.repCh],m]]; 
xmin=Min[xmin,TNum[[4,1]]];
xmax=Min[xmax,TNum[[4,2]]];AppendTo[Diagram,Dashing[None]];AppendTo[Diagram,Black];
AppendTo[Diagram,
If[(T[[1]]/.repCh)[[1]]==0, 
  Line[{T[[2]],T[[2]]+1.5{-Tan[psi],1}}],
  ParametricPlot[{Rays[T[[1]]/.repCh,t,psi,m],t},{t,T[[2,2]],T[[2,2]]+2},PlotStyle->{Black,AbsoluteThickness[1]}][[1]]]] 
,{i,Length[TreeList]}];
AppendTo[Diagram,Dotted];
AppendTo[Diagram,Black]; 
Graphics[Diagram]];


ScattDiagInternal[Tree_,m_:1/2]:=Module[{S1,S2,TreeNum,sInit,z,Li}(* construct total charge, coordinate of root and list of line segments in (s,t) coordinates, {min(x), max(x)} *), 
If[!ListQ[Tree]||Length[Tree]>2,
	TreeNum=Tree/.repCh;
sInit=InitialPosition[TreeNum,m];{Tree,{sInit,-1/2 sInit(sInit+m)},{Text[Tree/.repChO,InitialLabelPosition[sInit]]},{sInit,sInit}},
S1=ScattDiagInternal[Tree[[1]],m];
S2=ScattDiagInternal[Tree[[2]],m];
z=IntersectRays[S1[[1]]/.repCh,S2[[1]]/.repCh,m];
If[Length[z]==0,Print["Illegal tree"],
Li=S1[[3]];
AppendTo[Li,S2[[3]]];
AppendTo[Li,Arrow[{S1[[2]],xytost[z,m]}]];
AppendTo[Li,Arrow[{S2[[2]],xytost[z,m]}]];
{S1[[1]]+S2[[1]],z,Li,{Min[S1[[4,1]],S2[[4,1]]],Max[S1[[4,2]],S2[[4,2]]]}}]]];

ScattGraphInternal[Tree_,m_:1/2]:=Module[{S1,S2,TreeNum,sInit,z,Li}(* construct total charge, coordinate of root and list of line segments in (s,t) coordinates, {min(x), max(x)} *), 
If[!ListQ[Tree]||Length[Tree]>2,
	TreeNum=Tree/.repCh;
sInit=InitialPosition[TreeNum,m];{Tree,{sInit,-1/2 sInit(sInit+m)},{}},
S1=ScattGraphInternal[Tree[[1]],m];
S2=ScattGraphInternal[Tree[[2]],m];
z=IntersectRays[S1[[1]]/.repCh,S2[[1]]/.repCh,m];
If[Length[z]==0,Print["Illegal tree"],
Li={S1[[3]],S2[[3]],Arrow[{S1[[2]],z}],Arrow[{S2[[2]],z}]};
{S1[[1]]+S2[[1]],z,Li}]]];

ScattGraph[Tree_,m_:1/2]:=Module[{T,LiArrows,LiVertex},
(* extracts list of vertices in (x,y) plane and adjacency matrix *)
LiSlopes={};LiHeights={};
T=ScattGraphInternal[Tree,m];
LiArrows=Cases[Flatten[T[[3]]],x_Arrow]/.Arrow[x_]:>x;
LiVertex=Union[Flatten[LiArrows,1]];
{LiVertex,Table[If[i!=j,Sign[Count[LiArrows,{LiVertex[[i]],LiVertex[[j]]}]],0],{i,Length[LiVertex]},{j,Length[LiVertex]}]}];

ScattCheck[Tree_,m_:1/2]:=(* Check consistency of single tree, returns {charge,{xf,yf}} if tree is consistent, otherwise {charge,{}}; ignore whether leaves have Delta=0 or not *)
Module[{S1,S2,z,mf,mi,r,d1,d2,ch2},
mi=Floor[m];mf=m-mi;
If[!ListQ[Tree]||Length[Tree]>2,
(* tree consists of a single node *)
{r,d1,d2,ch2}=Tree/.repCh;
z=If[Disc[{r,d1,d2,ch2}]==0,
If[(d1==d2+mi r),{d2/r-mf/2,0},
If[(d1==d2+(1+mi)r),{d2/r+(1-mf)/2,0},{}],{}],{}];
If[Length[z]==0,Print["Illegal endpoint:",{r,d1,d2,ch2}],
z[[2]]=(ch2-d2 m-(d1+d2)z[[1]])/(2r);
(*Print["Initial pt:",z]; *)
{Tree/.repCh,z}],
(* otherwise, check each of the two branches *)
S1=ScattCheck[Tree[[1]],m];
S2=ScattCheck[Tree[[2]],m];
If[Length[S1[[2]]]>0&&Length[S2[[2]]]>0,
z=IntersectRays[S1[[1]]/.repCh,S2[[1]]/.repCh,S1[[2]],S2[[2]],m];(*Print[{S1[[1]],S2[[1]],S1[[2]],S2[[2]],z}];*){S1[[1]]+S2[[1]],z},
{S1[[1]]+S2[[1]],{}}]
]];

]; (* end of definition of LiSlopes, LiHeights *)

ScattSort[LiTree_,m_:1/2]:= (* sort trees by decreasing radius *)
Reverse[Map[#[[2]]&,SortBy[Table[{xytost[ScattCheck[LiTree[[i]],m][[2]],m][[2]],LiTree[[i]]},{i,Length[LiTree]}],N[First[#]]&]]];

Options[QuiverDomain]={"Style"->LightBlue};
QuiverDomain[Coll_,psi_,m_,OptionsPattern[]]:={RegionPlot[(And@@Table[Re[Exp[-I psi]ZLV[Coll[[i]],{s, t},m]]<0,{i,Length[Coll]}])&&t>0,{s,-1.5,1.5},{t,0,1},PlotPoints->100,AspectRatio->1,PlotStyle->Flatten[{OptionValue["Style"],Opacity[.5]}]],RegionPlot[SameHalfPlaneQ[Table[ZLV[Coll[[i]],{s, t},m],{i,Length[Coll]}]],{s,-1.5,1.5},{t,0,1},PlotPoints->100,AspectRatio->1,BoundaryStyle->Directive[Dashed],PlotStyle->Flatten[{OptionValue["Style"],Opacity[.3]}]]};

StabilityRegion[gam1_,gam2_,m_]:=Module[{},Print["Total charge:", gam1+gam2/.repCh,", DSZProduct:", DSZ[gam1,gam2]];RegionPlot[Im[ZLV[gam1/.repCh,{s,t},m]Conjugate[ZLV[gam2/.repCh,{s,t},m]]] DSZ[gam1,gam2]>0,{s,-2,2},{t,0,2}]];



(* ::Section:: *)
(*Constructing scattering trees in LV+epsilon*)


ListSubsetsAny=Module[{LiSub,LiSubHalf},
(* precompute all binary splits of Range[n] for n=2..10; used by ListStableTrees *)
Table[
LiSub=Drop[Subsets[Range[n],Floor[(n-1)/2]],1];
If[EvenQ[n],
LiSubHalf=Subsets[Range[n],{n/2}];
LiSub=Union[LiSub,Take[LiSubHalf,Length[LiSubHalf]/2]]];LiSub,{n,2,10}]];

ListStableTrees[LiCh_,{s0_,t0_},m_]:=
Module[{n,LiTree,LiCh1,LiCh2,gam1,gam2,s1,t1,z1,LiStable1,LiStable2},
(* for a given list of k_i Ch[m_i], construct consistent stable trees *)
n=Length[LiCh];
(*Print["**",LiCh,", (s,t)=",{s0,t0}]; *)
If[n==1, LiTree=LiCh,
(* construct possible subsets, avoiding double counting *) 
LiTree={};
Do[
LiCh1=ListSubsetsAny[[n-1,i]]/.n_Integer:>LiCh[[n]];
LiCh2=Complement[Range[n],ListSubsetsAny[[n-1,i]]]/.n_Integer:>LiCh[[n]];
(*Print[LiCh1,",",LiCh2];*)
gam1=Plus@@(LiCh1/.repCh);
gam2=Plus@@(LiCh2/.repCh);
(* reorder according to slope *)
If[GenSlope[gam1,m]>GenSlope[gam2,m],{gam1,gam2}={gam2,gam1};{LiCh1,LiCh2}={LiCh2,LiCh1}];
If[Wall[gam1,gam2,{s0,t0},m]DSZ[gam1,gam2]>0,
z1=IntersectRaysSt[gam1,gam2,0,m];
If [Length[z1]>0,{s1,t1}=z1;
(*Print[gam1,",",gam2,",",s1];*)
If[TestBranch[gam1,s1,m]&&TestBranch[gam2,s1,m],
(*Print[LiCh1,LiCh2,gam1,gam2,s1,t1]; *)
LiStable1=ListStableTrees[LiCh1,{s1,t1},m];
LiStable2=ListStableTrees[LiCh2,{s1,t1},m];
Do[
AppendTo[LiTree,{LiStable1[[j1]],LiStable2[[j2]]}],
{j1,Length[LiStable1]},{j2,Length[LiStable2]}];
]]];,
{i,Length[ListSubsetsAny[[n-1]]]}];
];LiTree
];

Wall[{r_,d1_,d2_,ch2_},{rr_,dd1_,dd2_,cch2_},{s_,t_},m_]:=d2 dd1 m-d1 dd2 m+dd2 m^2 r-d2 m^2 rr+2 dd2 m r s-2 d2 m rr s+dd1 r s^2+dd2 r s^2-d1 rr s^2-d2 rr s^2+cch2 (d1+d2-r (m+2 s))-ch2 (dd1+dd2-rr (m+2 s))+dd1 r t^2+dd2 r t^2-d1 rr t^2-d2 rr t^2;

ScanAllTrees[{r_,d1_,d2_,ch2_},{s0_,t0_},n_,m_]:=Module[
{Li,phimax},
phimax=d1+d2-r (2s0+m);
If[phimax<=0,Print["There are no trees at this point"];,
Li=ScanConstituents[{r,d1,d2,ch2},{s0-t0,s0+t0},{n,n},m,phimax];
Print[Li];
ScattSort[DeleteDuplicatesBy[Flatten[Select[Table[ListStableTrees[Li[[i]],{s0,t0},m],{i,Length[Li]}],Length[#]>0&],1],ScattGraph[#,m]&],m]]];


ScanConstituents[gam_,{smin_,smax_},{n_,np_},m_,phimax_]:=Module[
{eps=.001,Tabs,Tabsp,k,kp,Li,Lip,Lik,Likp,LiCompatible,LiConst,nc,i,ii,mi,mf},
mi=Floor[m];mf=m-mi;
(* allowed initial {position,charge,phi} *)
Tabs=Reverse[Sort[Flatten[{
Table[{k-mf+eps,Ch[k+mi,k],mf},{k,Ceiling[smin+mf],Floor[smax+mf]}],Table[{k+eps,Ch[k+1+mi,k],1-mf},{k,Ceiling[smin],Floor[smax]}]},1]]];
Tabsp=Sort[Flatten[{
Table[{k-eps,-Ch[k+mi,k],mf},{k,Ceiling[smin],Floor[smax]}],Table[{k+1-mf-eps,-Ch[k+1+mi,k],1-mf},{k,Ceiling[smin+m-1],Floor[smax+m-1]}]},1]];
Print["D-brane initial locations:",Tabs];
Print["Dbar-brane initial locations:",Tabsp];
(* choose up to n initial branes, and up to n initial antibranes *)
Li=Flatten[Table[Union[Map[Sort,Tuples[Range[Length[Tabs]],{k}]]],{k,1,n}],1];
Lip=Flatten[Table[Union[Map[Sort,Tuples[Range[Length[Tabsp]],{k}]]],{k,1,np}],1];
(*Print[Li];
Print[Lip]; *)
(* extract compatible choices of branes/antibranes *)
LiCompatible={};
Do[
If[
(* rightmost anti D-brane to the left of rightmost D-brane *)
Tabs[[Last[Li[[i]]],1]]>Tabsp[[First[Lip[[ii]]],1]]&&
(* leftmost anti D-brane to the left of lefttmost D-brane *)
Tabs[[First[Li[[i]]],1]]>Tabsp[[Last[Lip[[ii]]],1]]&&
(* condition on phi with unit multiplicity *)
Sum[Tabs[[Li[[i,j]],3]],{j,Length[Li[[i]]]}]+Sum[Tabsp[[Lip[[ii,jj]],3]],{jj,Length[Lip[[ii]]]}]<=phimax,AppendTo[LiCompatible,{i,ii}]],{i,Length[Li]},{ii,Length[Lip]}];
(*Print["Compatible pairs:", LiCompatible]; *)
nc=0;LiConst={};
Do[
(* Look for compatible choices leading to the correct total charge *)
i=LiCompatible[[l,1]];
ii=LiCompatible[[l,2]];
Lik=Table[{k[j],
If[j>1,If[Li[[i,j]]==Li[[i,j-1]],k[j-1],1],1],phimax-Sum[Tabs[[Li[[i,jj]],3]] k[jj],{jj,1,j-1}]},{j,Length[Li[[i]]]}];
Likp=Table[{kp[j],
If[j>1,If[Lip[[ii,j]]==Lip[[ii,j-1]],kp[j-1],1],1],phimax-Sum[Tabs[[Lip[[ii,jj]],3]]kp[jj],{jj,1,j-1}]},{j,Length[Lip[[ii]]]}];
(*Print[{i,ii,Lik,Likp}];*)
Do[Do[(* loop on kp *)
nc+=1;
If[(Sum[k[j]Tabs[[Li[[i,j]],2]],{j,Length[Li[[i]]]}]+Sum[kp[j]Tabsp[[Lip[[ii,j]],2]],{j,Length[Lip[[ii]]]}]/.repCh)==gam,AppendTo[LiConst,Join[
Table[k[j]Tabs[[Li[[i,j]],2]],{j,Length[Li[[i]]]}],Table[kp[j]Tabsp[[Lip[[ii,j]],2]],{j,Length[Lip[[ii]]]}]]
]],
##]&@@Likp,
##]&@@Lik;,{l,Length[LiCompatible]}];
Print[Length[LiConst]," possible sets of constituents out of ",nc," trials"];
LiConst
];

ContractInitialRays[Trees_,m_:1/2]:=Trees /.repCh/.{{r_Integer,d1_Integer,d2_Integer,ch2_Integer},{rr_Integer,dd1_Integer,dd2_Integer,cch2_Integer}}:>{r+rr,d1+dd1,d2+dd2,ch2+cch2}/;InitialPosition[{r,d1,d2,ch2},m]==InitialPosition[{rr,dd1,dd2,cch2},m]/.{{r_Integer,d1_Integer,d2_Integer,ch2_Integer}:>If[r==0,GV[d1,d2,ch2],If[Disc[{r,d1,d2,ch2}]==0,If[r>0,r Ch[d1/r,d2/r],-r Ch[d1/r,d2/r][1]],{r,d1,d2,ch2}]]};

DecontractInitialRays[Trees_,m_:1/2]:=Trees/.
Ch[d1_,d2_][1]:>Module[{s0,k},If[d1-d2==Floor[m]|| d1-d2==Floor[m]+1,Ch[d1,d2][1],
s0=InitialPosition[-{1,d1,d2,d1 d2},m]; k=Floor[s0+(m-Floor[m])/2];
If[s0+(m-Floor[m])/2<=k+1/2,
{(1-d1+k+Floor[m]) Ch[k+Floor[m],k][1],-(d1-k-Floor[m]) Ch[k+1+Floor[m],k]},
{-(d2+k+1) Ch[k+1+Floor[m],k][1],(-2-d2-k) Ch[k+1+Floor[m],k+1]}]]]/.
Ch[d1_,d2_]:>Module[{s0,k},If[d1-d2==Floor[m]|| d1-d2==Floor[m]+1,Ch[d1,d2],
s0=InitialPosition[{1,d1,d2,d1 d2},m]; k=Floor[s0+(m-Floor[m])/2];
If[s0+(m-Floor[m])/2<=k+1/2,
{(-1+d1-k-Floor[m]) Ch[k+Floor[m],k][1],(d1-k-Floor[m]) Ch[k+1+Floor[m],k]},
{(d2-k-1) Ch[k+1+Floor[m],k][1],(d2-k) Ch[k+1+Floor[m],k+1]}]]]/. GV[1,0,k_]->{Ch[k+1+Floor[m], k],Ch[k+Floor[m],k][1]} /.
GV[0,1,k_]->{Ch[k,k-Floor[m]],Ch[k, k-Floor[m]-1][1]};

ClearAll[KroneckerDims];
KroneckerDims[m_,Nn_]:=KroneckerDims[m,Nn]=Module[{Ta={}},
Do[If[m n1 n2-n1^2-n2^2+1>=0&&GCD[n1,n2]==1,AppendTo[Ta,{n1,n2}]],{n1,0,Nn},{n2,0,Nn-n1}];Drop[Ta,2]];

ConstructLVDiagram[smin_,smax_,phimax_,Nm_,m_,ListRays0_]:=Module[
{eps=.001,mi,mf,Inter,ListInter,ListRays,ListNewRays,kappa,KTab},
mi=Floor[m];mf=m-mi;
(* initial rays {charge, {x,y}, parent1, parent2,n1,n2 } *)
If[ListRays0=={},
         ListRays=Flatten[{
Table[{Ch[k+mi,k],{k-mf/2,-k^2/2-k/2mi+1/4 mf mi},0,0,0,0},{k,Ceiling[smin+mf/2],Floor[smax+mf/2]}],
Table[{Ch[k+mi,k][1],{k-mf/2,-k^2/2-k/2mi+1/4 mf mi},0,0,0,0},{k,Ceiling[smin+mf/2],Floor[smax+mf/2]}],Table[{Ch[k+1+mi,k],{k+(1-mf)/2,-k^2/2-k/2(1+mi)-1/4(1-mf)(1+mi)},0,0,0,0},{k,Ceiling[smin+(mf-1)/2],Floor[smax+(mf-1)/2]}],
Table[{Ch[k+1+mi,k][1],{k+(1-mf)/2,-k^2/2-k/2(1+mi)-1/4(1-mf)(1+mi)},0,0,0,0},{k,Ceiling[smin+(mf-1)/2],Floor[smax+(mf-1)/2]}]},1]/.repCh;
   ListInter={};,
(* If list of rays is already provided *)
ListRays=ListRays0;
ListInter=Select[Table[{ListRays[[i,3]],ListRays[[i,4]]},{i,Length[ListRays]}],First[#]>0&]];
While[True,
ListNewRays={};
       Monitor[ Do[
If[  !MemberQ[ListInter,{i,j}],
AppendTo[ListInter,{i,j}];
 kappa=DSZ[ListRays[[i,1]],ListRays[[j,1]]];
If[kappa!=0,Inter=IntersectRaysNoTest[ListRays[[i,1]],ListRays[[j,1]],ListRays[[i,2]],ListRays[[j,2]],m];
If[Inter!={},
KTab=KroneckerDims[Abs[kappa],Nm];
Do[If[CostPhi[KTab[[k,1]] ListRays[[i,1]]+KTab[[k,2]]ListRays[[j,1]],Inter[[1]],m]<=phimax,
AppendTo[ListNewRays,{KTab[[k,1]]ListRays[[i,1]]+KTab[[k,2]]ListRays[[j,1]],Inter,i,j,KTab[[k,1]],KTab[[k,2]]}]],{k,Length[KTab]}]
]]]
,{i,Length[ListRays]},{j,i+1,Length[ListRays]}],{i,j}];
If[ListNewRays=={},Break[],
Print["Adding ",Length[ListNewRays], " rays, "];
ListRays=Flatten[{ListRays,ListNewRays},1];
]];
Print[Length[ListRays], " in total."];
ListRays];

CostPhi[{r_,d1_,d2_,ch2_},s_,m_]:=d1+d2-r(2s+m);

IntersectRaysNoTest[{r_,d1_,d2_,ch2_},{rr_,dd1_,dd2_,cch2_},z_,zz_,m_:1/2]:=
(* returns (x,y) coordinate of intersection point of two rays, or {} if they don't intersect *)
(* here do not test if DSZ<>0, and require strictly in future of z and zz *)
Module[{zi},zi={-((-cch2 r+dd2 m r+ch2 rr-d2 m rr)/(dd1 r+dd2 r-d1 rr-d2 rr)),(-cch2 d1-cch2 d2+ch2 dd1+ch2 dd2-d2 dd1 m+d1 dd2 m)/(2 (dd1 r+dd2 r-d1 rr-d2 rr))};
If[(zi-z) . {-2r,d1+d2-m r}>0&&(zi-zz) . {-2rr,dd1+dd2-m rr}>0,zi,{}]];




(* Extract tree leading to k-th ray, internal *)
TreeFromListRays[ListRays_,k_]:=If[ListRays[[k,3]]==0,ListRays[[k,1]],{ListRays[[k,5]]TreeFromListRays[ListRays,ListRays[[k,3]]],ListRays[[k,6]]TreeFromListRays[ListRays,ListRays[[k,4]]]}];
(* Extract all trees leading to a ray with charge {r,d1,d2,ch2} *)

LVTreesFromListRays[ListRays_,{r_,d1_,d2_,ch2_},m_]:=Module[{Lipos,Div,LiTrees},
Div=Divisors[GCD@@{r,d1,d2,ch2}];
Lipos=Flatten[Join[Table[Position[ListRays,{r,d1,d2,ch2}/k],{k,Div}]],1];
If[Lipos=={},
Print["No such dimension vector in the list"],
LiTrees=(GCD[r,d1,d2,ch2]/GCD@@ListRays[[#,1]])TreeFromListRays[ListRays,#]&/@First[Transpose[Lipos]];
ScattSort[DeleteDuplicatesBy[SortBy[LiTrees,Length[TreeConstituents[#]]&],ScattGraph[#,m]&],m]
]];

IntersectRaysNoTest[{r_,d1_,d2_,ch2_},{rr_,dd1_,dd2_,cch2_},z_,zz_,m_:1/2]:=
(* returns (x,y) coordinate of intersection point of two rays, or {} if they don't intersect *)
(* here do not test if DSZ<>0, and require strictly in future of z and zz *)
Module[{zi},zi={-((-cch2 r+dd2 m r+ch2 rr-d2 m rr)/(dd1 r+dd2 r-d1 rr-d2 rr)),(-cch2 d1-cch2 d2+ch2 dd1+ch2 dd2-d2 dd1 m+d1 dd2 m)/(2 (dd1 r+dd2 r-d1 rr-d2 rr))};
If[(zi-z) . {-2r,d1+d2-m r}>0&&(zi-zz) . {-2rr,dd1+dd2-m rr}>0,zi,{}]];


LVTreesF0[{0,1,0,0}]={{Ch[1,0],Ch[0,0][1]}};
LVTreesF0[{0,0,1,0}]={{Ch[0,0],Ch[0,-1][1]}};
LVTreesF0[{0,1,1,0}]={{Ch[0,-1][1],Ch[1,0]}};
LVTreesF0[{0,1,2,0}]={{Ch[0,-2][1],Ch[1,0]}};
LVTreesF0[{0,2,2,0}]={{Ch[1,1],Ch[-1,-1][1]},{{2 Ch[1,0],Ch[0,0][1]},{Ch[0,0],2 Ch[0,-1][1]}},{2 Ch[1,0],2 Ch[0,-1][1]}};
LVTreesF0[{0,2,3,0}]={{{Ch[-1,-1][1],{Ch[0,0],Ch[0,-1][1]}},Ch[1,1]},{{3 Ch[0,-1][1],2 Ch[0,0]},{Ch[0,0][1],2 Ch[1,0]}},{{Ch[0,-1][1],2 Ch[1,0]},{2 Ch[0,-1][1],Ch[0,0]}}};
LVTreesF0[{1,0,0,-1}]={{Ch[0,-1],{Ch[-1,-1],Ch[-1,-2][1]}}};
LVTreesF0[{1,0,0,-2}]={{Ch[0,-1],{Ch[-2,-2],Ch[-2,-3][1]}},{{Ch[-1,-2][1],2 Ch[-1,-1]},{Ch[-1,-2],Ch[-2,-2][1]}},{Ch[-2,-2][1],2 Ch[-1,-1]}};
LVTreesF0[{1,0,0,-3}]={{Ch[0,-1],{Ch[-3,-3],Ch[-3,-4][1]}},{{Ch[-1,-2][1],2 Ch[-1,-1]},{Ch[-2,-3],Ch[-3,-3][1]}},{Ch[-1,-1],{Ch[-1,-2],Ch[-2,-3][1]}},{{Ch[-2,-2],Ch[-2,-3][1]},{{Ch[-1,-2],Ch[-2,-2][1]},Ch[-1,-1]}}};
LVTreesF0[{2,-1,0,-1}]={{Ch[0,-1],{2 Ch[-1,-2][1],3 Ch[-1,-1]}}};
LVTreesF0[{2,0,-1,-1}]={{2 Ch[0,-1],{Ch[-1,-1],Ch[-1,-2][1]}}};
LVTreesF0[{2,-1,-1,-1}]={{Ch[-2,-2][1],3 Ch[-1,-1]}};
LVTreesF0[{3,1,1,-1}]={{4 Ch[0,0],-Ch[-1,-1]}};
LVTreesF0[{3,2,2,0}]={{2 Ch[1,0],{2 Ch[0,-1][1],3 Ch[0,0]}}};
LVTreesF0[{5,3,4,0}]={{3 Ch[1,0],{4 Ch[0,-1][1],6 Ch[0,0]}}};
LVTreesF0[{5,1,2,-2}]={{6Ch[0,0],Ch[-1,-2][1]}};



(* ::Section:: *)
(*Common routines for quiver scattering diagram*)


repKr=Subscript[Kr, m_][p_,q_]:>1;

McKayrep={{n1_Integer,n2_Integer,n3_Integer,n4_Integer}->n1 Subscript[\[Gamma], 1]+n2 Subscript[\[Gamma], 2]+n3 Subscript[\[Gamma], 3]+n4 Subscript[\[Gamma], 4]};
McKayrepi={Subscript[\[Gamma], i_]:>IdentityMatrix[4][[i]]};
TreeToGamma[Tree_]:=Tree/.{Nvec_:>\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i\), \(Length[Nvec]\)]\(Nvec[\([i]\)]\ 
\*SubscriptBox[\(\[Gamma]\), \(i\)]\)\)/;Length[Nvec]>2};
McKayRay[Nvec_,{u_,v_},{k1_,k2_},tx_]:=(* produces from {u,v}+k1 vec to (u,v)+k2 vec, decorated with text at the target *)
{Arrow[{{u,v}+k1 McKayVec[Nvec],{u,v}+k2 McKayVec[Nvec]}],Text[tx,{u,v}+(k2+.1) McKayVec[Nvec]]};

McKayIntersectRays[Nvec_,NNvec_,z_,zz_,m_]:=Module[{zi},If[McKayDSZ[Nvec,NNvec]!=0,zi=McKayIntersectRays[Nvec,NNvec,m];If[(zi-z) . McKayVec[Nvec]>=0&&(zi-zz) . McKayVec[NNvec]>=0,zi,{}]]];

(*McKayInitialRays[L_]:=
Graphics[Table[McKayRay[IdentityMatrix[Length[Mat]][[i]],InitialRaysOrigin[[i]],L{-1,1},"\!\(\*SubscriptBox[\(\[Gamma]\), \("<>ToString[i]<>"\)]\)"],{i,Length[Mat]}]];
*)

McKayInitialRays[L_,m_]:=
Graphics[Table[{Thick,Red,McKayRay[IdentityMatrix[Length[Mat]][[i]],InitialRaysOrigin[m][[i]],L{-1,1},"\!\(\*SubscriptBox[\(\[Gamma]\), \("<>ToString[i]<>"\)]\)"]},{i,Length[Mat]}]];

McKayScattDiag[TreeList_,m_]:=Module[{T,TNum,Diagram},
(* Draws scattering diagram in (u,v) plane for each tree in Treelist *)
Diagram={};Do[
T=McKayScattDiagInternal[TreeList[[i]],m];
TNum=T/.repChn;
AppendTo[Diagram,TreeHue[Length[TreeList],i]];
AppendTo[Diagram,T[[3]]];
AppendTo[Diagram,Arrow[{TNum[[2]],TNum[[2]]+McKayVec[TNum[[1]]]/GCD@@(TNum[[1]])}]];
AppendTo[Diagram,Text[TNum[[1]]/.McKayrep,TNum[[2]]+1.2 McKayVec[TNum[[1]]]/GCD@@(TNum[[1]])]];
,{i,Length[TreeList]}];
Graphics[Diagram]];

McKayScattDiagInternal[Tree_,m_]:=Module[{S1,S2,TreeNum,z,Li},
(* construct total charge, coordinate of root and list of line segments in (x,y) coordinates *) 
If[!ListQ[Tree]||Length[Tree]>2,
(* tree consists of a single node with charge {0,0,..,p,0,0,..} *)
TreeNum=Tree/.repChn;z={};
If[Length[Position[TreeNum,n_Integer/;n>0]]==1,
	z=InitialRaysOrigin[m][[Tr[Position[TreeNum,n_Integer/;n>0]]]]-2 McKayVec[TreeNum/GCD@@TreeNum];,Print["Illegal endpoint"]];{Tree,z,{Text[Tree/.repChO/.McKayrep,z-.2 McKayVec[TreeNum]]}},
(* non-trivial tree *)
S1=McKayScattDiagInternal[Tree[[1]],m];
S2=McKayScattDiagInternal[Tree[[2]],m];
z=McKayIntersectRays[S1[[1]]/.repChn,S2[[1]]/.repChn,m];
If[Length[z]==0,Print["Illegal tree"],
Li=S1[[3]];AppendTo[Li,S2[[3]]];AppendTo[Li,Arrow[{S1[[2]],z}]];AppendTo[Li,Arrow[{S2[[2]],z}]];
{S1[[1]]+S2[[1]],z,Li}]]];

McKayScattCheck[Tree_,m_]:=(* Check consistency of single tree, returns {charge,{xf,yf}} if tree is consistent, otherwise {charge,{}} *)
Module[{S1,S2,TreeNum,z},
If[!ListQ[Tree]||Length[Tree]>2,
(* tree consists of a single node *)
TreeNum=Tree/.repChn;
If[Length[Position[TreeNum,n_Integer/;n>0]]==1,
	        z=InitialRaysOrigin[m][[Tr[Position[TreeNum,n_Integer/;n>0]]]]-20 McKayVec[TreeNum/GCD@@TreeNum],
z=-20 McKayVec[TreeNum]];
{Tree/.repChn,z},
(* otherwise, check each of the two branches *)
S1=McKayScattCheck[Tree[[1]],m];
S2=McKayScattCheck[Tree[[2]],m];
If[Length[S1[[2]]]>0&&Length[S2[[2]]]>0,{S1[[1]]+S2[[1]],McKayIntersectRays[S1[[1]]/.repChn,S2[[1]]/.repChn,S1[[2]],S2[[2]],m]},
{S1[[1]]+S2[[1]],{}}]
]];


McKayScattIndex[TreeList_]:=Table[
(* compute index for each tree in the list; do not trust the result if internal lines have non-primitive charges *)
Simplify[Times@@McKayScattIndexInternal[TreeList[[i]]][[2]]],{i,Length[TreeList]}];

McKayScattIndexInternal[Tree_]:=Module[{S1,S2,g1,g2,kappa,Li},
(* compute {total charge, list of Kronecker indices associated to each vertex *)
If[!ListQ[Tree]||Length[Tree]>2,{Tree,{1}},
S1=McKayScattIndexInternal[Tree[[1]]]/.repChn;
S2=McKayScattIndexInternal[Tree[[2]]]/.repChn;
g1=GCD@@S1[[1]];g2=GCD@@S2[[1]];
kappa=Abs[McKayDSZ[S1[[1]],S2[[1]]]]/g1/g2;
Li=Join[S1[[2]],S2[[2]]];
AppendTo[Li,Subscript[Kr, kappa][Min[g1,g2],Max[g1,g2]]];
(*If[GCD@@(S1[[1]]+S2[[1]])!=1,Print["Beware, non-primitive state"]];*)
{S1[[1]]+S2[[1]],Li}]];

(* more careful implementation taking care of internal non-primitive states *)
Options[McKayScattIndexImprovedInternal] = {"Debug"->False};
McKayScattIndexImproved[TreeList_, opt: OptionsPattern[]]:=Table[
	(* compute index for each tree in the list; do not trust the result if internal lines have non-primitive charges *)
	Simplify[FOmbToOm[Last@McKayScattIndexImprovedInternal[TreeList[[i]], opt][[2]]]],{i,Length[TreeList]}];

McKayScattIndexImprovedInternal[Tree_, opt: OptionsPattern[]]:=Module[{S1,S2,g1,g2,gFinal, kappa,Li, tem, repOmAttb, rrr},
(* compute {total charge, list of Kronecker indices associated to each vertex *)
	If[!ListQ[Tree]||Length[Tree]>2,{Tree,{Join[{1}, Table[(y-y^-1)/(j(y^j-y^-j)), {j, 2, GCD@@Tree}]]}},
	If[OptionValue["Debug"], Print["Calling with args: ", Tree[[1]], "  |  ", Tree[[2]]]];
S1=McKayScattIndexImprovedInternal[Tree[[1]], opt]/.repChn;
	S2=McKayScattIndexImprovedInternal[Tree[[2]], opt]/.repChn;
If[OptionValue["Debug"], Print["S1 is: ", S1, "   S2 is: ", S2]];
	g1=GCD@@S1[[1]];g2=GCD@@S2[[1]];
gFinal = GCD@@(S1[[1]]+S2[[1]]);
	kappa=Abs[McKayDSZ[S1[[1]],S2[[1]]]]/g1/g2;
	Li=Join[S1[[2]],S2[[2]]];
If[OptionValue["Debug"], Print["Li is: ", Li, "  g1 is: ", g1, "  g2 is: ", g2, "  gFinal is: ", gFinal]];
AppendTo[Li,
repOmAttb = Join[
Table[CoulombHiggs`OmAttb[{P, 0}, y_]->Last[S1[[2]]][[P]], {P, 1, g1}],
Table[CoulombHiggs`OmAttb[{0, Q}, y_]->Last[S2[[2]]][[Q]], {Q, 1, g2}]
];
If[OptionValue["Debug"], Print["repOmAttb is: ", repOmAttb]];
tem = Table[
rrr = If[And@@(IntegerQ/@{P g1/gFinal, P g2/gFinal}), CoulombHiggs`FlowTreeFormulaRat[{{0, kappa}, {-kappa, 0}}, {g2, -g1}, {P g1/gFinal, P g2/gFinal}, y], 0];
Simplify[
rrr
/.repOmAttb
/.{CoulombHiggs`OmAttb[{p_, q_}, y___]:>0/;p>1||q>1||p q !=0}
],
{P, 1, gFinal}
];
If[OptionValue["Debug"], Print["tem is: ", tem]];tem
];
	(*If[GCD@@(S1[[1]]+S2[[1]])!=1,Print["Beware, non-primitive state"]];*)
	{S1[[1]]+S2[[1]],Li}]];

FOmbToOm[OmbList_] := Module[{n},
If[Length[OmbList]<2, First@OmbList, 
n = Length[OmbList];
DivisorSum[n, (MoebiusMu[#] (y-y^-1)/(#(y^#-y^-#)) (OmbList[[n/#]]/.{y->y^#}))&]
(*Simplify[OmbList[[-1]] -Sum[ (y-y^-1)/(P(y^P-y^-P))(FOmbToOm[OmbList[[;;-P]]]/.{y->y^P}), {P, 2, n}]]*)
]];

McKayScattGraph[Tree_,m_]:=Module[{T,LiArrows,LiVertex},
(* extracts list of vertices and adjacency matrix *)
T=McKayScattDiagInternal[Tree,m];
LiArrows=Cases[Flatten[T[[3]]],x_Arrow]/.Arrow[x_]:>x;
LiVertex=Union[Flatten[LiArrows,1]];
{LiVertex,Table[If[i!=j,Sign[Count[LiArrows,{LiVertex[[i]],LiVertex[[j]]}]],0],{i,Length[LiVertex]},{j,Length[LiVertex]}]}];

McKayListAllTrees[Nvec_]:=Module[{LiTrees,LiTree1,LiTree2,Li},
(* generate all trees with leaves carrying charge {p,0,0,0}, {0,p,0,0}, etc and with non-zero DSZ pairing at each vertex *)
LiTrees={};
If[Count[Nvec,0]>=Length[Mat]-1,
LiTrees={Nvec};,
Li=Select[F0BinarySplits[Nvec],McKayDSZ[#,Nvec]!=0 &];
Do[
LiTree2=McKayListAllTrees[Li[[i]]];
LiTree1=McKayListAllTrees[Nvec-Li[[i]]];
Do[AppendTo[LiTrees,{LiTree1[[j]],LiTree2[[k]]}],{j,Length[LiTree1]},{k,Length[LiTree2]}];
,{i,Length[Li]}];
];LiTrees];

McKayScanAllTrees[Nvec_,m_]:=Module[{Li,Li2},
(* generate consistent scattering trees with leaves carrying charge {p,0,0,0}, {0,p,0,0}, etc with non-zero DSZ pairing at each vertex, with distinct support*)
Li=McKayListAllTrees[Nvec];
Li2=Select[Li,Length[McKayScattCheck[#,m][[2]]]>0&];
DeleteDuplicatesBy[ReverseSortBy[Li2,Length],McKayScattGraph[#,m]&]];

(* construct scattering diagram up to height phimax *)
ConstructMcKayDiagram[phimax_,m_,ListRays0_]:=Module[{ListRays,ListInter,kappa,Inter,KTab,ListNewRays,Nm},
If[ListRays0=={},
(* initial rays {charge, {x,y}, parent1, parent2,n1,n2,level } *)
ListRays=Table[{IdentityMatrix[4][[k]],InitialRaysOrigin[m][[k]]-5 McKayVec[IdentityMatrix[4][[k]]],0,0,0,0,1},{k,4}];
ListInter={};,
(* If list of rays is already provided *)
ListRays=ListRays0;
ListInter=Select[Table[{ListRays[[i,3]],ListRays[[i,4]]},{i,Length[ListRays]}],First[#]>0&]];
While[True,
ListNewRays={};
       Monitor[ Do[
If[  !MemberQ[ListInter,{i,j}],
AppendTo[ListInter,{i,j}];
 kappa=McKayDSZ[ListRays[[i,1]],ListRays[[j,1]]];
If[kappa!=0,Inter=McKayIntersectRaysNoTest[ListRays[[i,1]],ListRays[[j,1]],ListRays[[i,2]],ListRays[[j,2]],m];
If[Inter!={},
Nm=Floor[phimax/Min[Plus@@ListRays[[i,1]],Plus@@ListRays[[j,1]]]];
KTab=KroneckerDims[Abs[kappa],Nm];
Do[If[Plus@@(KTab[[k,1]] ListRays[[i,1]]+KTab[[k,2]]ListRays[[j,1]])<=phimax,
AppendTo[ListNewRays,{KTab[[k,1]]ListRays[[i,1]]+KTab[[k,2]]ListRays[[j,1]],Inter,i,j,KTab[[k,1]],KTab[[k,2]],ListRays[[i,7]]+ListRays[[j,7]]}]],{k,Length[KTab]}]
]]]
,{i,Length[ListRays]},{j,i+1,Length[ListRays]}],{i,j}];
If[ListNewRays=={},Break[],
Print["Adding ",Length[ListNewRays], " rays, "];
ListRays=Flatten[{ListRays,ListNewRays},1];
]];
Print[Length[ListRays], " in total."];
ListRays];

McKayIntersectRaysNoTest[Nvec_,NNvec_,z_,zz_,m_]:=
(* require strict inequality *)
Module[{zi},zi=McKayIntersectRays[Nvec,NNvec,m];If[(zi-z) . McKayVec[Nvec]>0&&(zi-zz) . McKayVec[NNvec]>0,zi,{}]];

(* Extract all trees leading up to a ray with dimension vector {n1,n2,n3,n4} *)
McKayTreesFromListRays[ListRays_,{n1_,n2_,n3_,n4_},m_]:=Module[{Lipos,Div,LiTrees},
Div=Divisors[GCD@@{n1,n2,n3,n4}];
Lipos=Flatten[Join[Table[Position[ListRays,{n1,n2,n3,n4}/k],{k,Div}]],1];
If[Lipos=={},
Print["No such dimension vector in the list"],
LiTrees=((n1+n2+n3+n4)/Plus@@ListRays[[#,1]])TreeFromListRays[ListRays,#]&/@First[Transpose[Lipos]];
DeleteDuplicatesBy[SortBy[LiTrees,Length[TreeConstituents[#]]&],McKayScattGraph[#,m]&]
]];



(* ::Section:: *)
(*Routines for quiver scattering diagram F0 phase I *)


gam1p=-{1,0,0,0};gam2p=-{-1,1,0,0};gam3p=-{-1,0,1,0};gam4p=-{1,-1,-1,1};
nIFromCh[{rk_,c11_,c12_,ch2_}]:=-{c11+c12+ch2+rk,c11+ch2,c12+ch2,ch2};
ChFromnI[{n1_,n2_,n3_,n4_}]:=-{n1-n2-n3+n4,n2-n4,n3-n4,n4};
(* 
{gamd[1],gamd[2],gamd[3],gamd[4]}={gam1p,gam2p,gam3p,gam3p};
repChO={Ch[m1_,m2_]:>"O("<>ToString[m1]<>","<>ToString[m2]<>")"};
nIFromCh[{rk_,c11_,c12_,ch2_}]:=-{c11+c12+ch2+rk,c11+ch2,c12+ch2,ch2};
ChFromnI[{n1_,n2_,n3_,n4_}]:=-{n1-n2-n3+n4,n2-n4,n3-n4,n4};
repChnI={Ch[m1_,m2_]:>nIFromCh[{1,m1,m2,m1 m2}]};
MatI={{0,2,2,-4},{-2,0,0,2},{-2,0,0,2},{4,-2,-2,0}};
{gamd[1],gamd[2],gamd[3],gamd[4]}={{1,0,0,0},{-1,1,0,0},{-1,0,1,0},{1,-1,-1,1}};

McKayDSZ[Nvec_,NNvec_]:=Tr[{Nvec}.MatI.Transpose[{NNvec}]];
McKayVec[Nvec_]:={-Nvec\[LeftDoubleBracket]1\[RightDoubleBracket]+Nvec\[LeftDoubleBracket]4\[RightDoubleBracket],Nvec\[LeftDoubleBracket]2\[RightDoubleBracket]+Nvec[[3]]-2Nvec\[LeftDoubleBracket]4\[RightDoubleBracket]};
McKayRayEq[Nvec_,{u_,v_}]:=(Nvec\[LeftDoubleBracket]2\[RightDoubleBracket]+Nvec\[LeftDoubleBracket]3\[RightDoubleBracket]-2Nvec\[LeftDoubleBracket]4\[RightDoubleBracket])u+(Nvec\[LeftDoubleBracket]1\[RightDoubleBracket]-Nvec\[LeftDoubleBracket]4\[RightDoubleBracket])v+(Nvec\[LeftDoubleBracket]4\[RightDoubleBracket]-Nvec\[LeftDoubleBracket]1\[RightDoubleBracket]-Nvec\[LeftDoubleBracket]2\[RightDoubleBracket]-3Nvec\[LeftDoubleBracket]3\[RightDoubleBracket]);
McKayIntersectRays[Nvec_,NNvec_]:={-(((NNvec\[LeftDoubleBracket]1\[RightDoubleBracket]-NNvec\[LeftDoubleBracket]4\[RightDoubleBracket]) (Nvec\[LeftDoubleBracket]2\[RightDoubleBracket]+3 Nvec\[LeftDoubleBracket]3\[RightDoubleBracket])-3 NNvec\[LeftDoubleBracket]3\[RightDoubleBracket] (Nvec\[LeftDoubleBracket]1\[RightDoubleBracket]-Nvec\[LeftDoubleBracket]4\[RightDoubleBracket])+NNvec\[LeftDoubleBracket]2\[RightDoubleBracket] (-Nvec\[LeftDoubleBracket]1\[RightDoubleBracket]+Nvec\[LeftDoubleBracket]4\[RightDoubleBracket]))/(-((NNvec\[LeftDoubleBracket]1\[RightDoubleBracket]-NNvec\[LeftDoubleBracket]4\[RightDoubleBracket]) (Nvec\[LeftDoubleBracket]2\[RightDoubleBracket]+Nvec\[LeftDoubleBracket]3\[RightDoubleBracket]-2 Nvec\[LeftDoubleBracket]4\[RightDoubleBracket]))+(NNvec\[LeftDoubleBracket]2\[RightDoubleBracket]+NNvec\[LeftDoubleBracket]3\[RightDoubleBracket]-2 NNvec\[LeftDoubleBracket]4\[RightDoubleBracket]) (Nvec\[LeftDoubleBracket]1\[RightDoubleBracket]-Nvec\[LeftDoubleBracket]4\[RightDoubleBracket]))),(2 NNvec\[LeftDoubleBracket]4\[RightDoubleBracket] Nvec\[LeftDoubleBracket]1\[RightDoubleBracket]+NNvec\[LeftDoubleBracket]1\[RightDoubleBracket] Nvec\[LeftDoubleBracket]2\[RightDoubleBracket]+NNvec\[LeftDoubleBracket]4\[RightDoubleBracket] Nvec\[LeftDoubleBracket]2\[RightDoubleBracket]+NNvec\[LeftDoubleBracket]1\[RightDoubleBracket] Nvec\[LeftDoubleBracket]3\[RightDoubleBracket]+5 NNvec\[LeftDoubleBracket]4\[RightDoubleBracket] Nvec\[LeftDoubleBracket]3\[RightDoubleBracket]-2 NNvec\[LeftDoubleBracket]1\[RightDoubleBracket] Nvec\[LeftDoubleBracket]4\[RightDoubleBracket]-NNvec\[LeftDoubleBracket]2\[RightDoubleBracket] (Nvec\[LeftDoubleBracket]1\[RightDoubleBracket]+2 Nvec\[LeftDoubleBracket]3\[RightDoubleBracket]+Nvec\[LeftDoubleBracket]4\[RightDoubleBracket])-NNvec\[LeftDoubleBracket]3\[RightDoubleBracket] (Nvec\[LeftDoubleBracket]1\[RightDoubleBracket]-2 Nvec\[LeftDoubleBracket]2\[RightDoubleBracket]+5 Nvec\[LeftDoubleBracket]4\[RightDoubleBracket]))/(2 NNvec\[LeftDoubleBracket]4\[RightDoubleBracket] Nvec\[LeftDoubleBracket]1\[RightDoubleBracket]+NNvec\[LeftDoubleBracket]1\[RightDoubleBracket] Nvec\[LeftDoubleBracket]2\[RightDoubleBracket]-NNvec\[LeftDoubleBracket]4\[RightDoubleBracket] Nvec\[LeftDoubleBracket]2\[RightDoubleBracket]+NNvec\[LeftDoubleBracket]1\[RightDoubleBracket] Nvec\[LeftDoubleBracket]3\[RightDoubleBracket]-NNvec\[LeftDoubleBracket]4\[RightDoubleBracket] Nvec\[LeftDoubleBracket]3\[RightDoubleBracket]-2 NNvec\[LeftDoubleBracket]1\[RightDoubleBracket] Nvec\[LeftDoubleBracket]4\[RightDoubleBracket]+NNvec\[LeftDoubleBracket]2\[RightDoubleBracket] (-Nvec\[LeftDoubleBracket]1\[RightDoubleBracket]+Nvec\[LeftDoubleBracket]4\[RightDoubleBracket])+NNvec\[LeftDoubleBracket]3\[RightDoubleBracket] (-Nvec\[LeftDoubleBracket]1\[RightDoubleBracket]+Nvec\[LeftDoubleBracket]4\[RightDoubleBracket]))};
InitialRaysOrigin={{3,1},{1,-1},{3,-5},{0,1}};

McKayIinitialRays[L_]:=
Graphics[{McKayIRay[{1,0,0,0},{1,1},L{-3,5},"\!\(\*SubscriptBox[\(\[Gamma]\), \(1\)]\)"],
McKayIRay[{0,1,0,0},{1,-1},L{-4,5.25},"\!\(\*SubscriptBox[\(\[Gamma]\), \(2\)]\)"],
McKayIRay[{0,0,1,0},{3,-5},L{-1,8},"\!\(\*SubscriptBox[\(\[Gamma]\), \(3\)]\)"],
McKayIRay[{0,0,0,1},{0,1},L{-1,2.5},"\!\(\*SubscriptBox[\(\[Gamma]\), \(4\)]\)"]}]; *)



(* ::Section:: *)
(*Routines for quiver scattering diagram F0 phase II*)


gam1=-{1,0,0,0};gam2=-{-1,1,0,0};gam3=-{-1,-1,1,1};gam4=-{1,0,-1,0};
{gamd[1],gamd[2],gamd[3],gamd[4]}={gam1,gam2,gam3,gam4};
repChO={Ch[m1_,m2_][1]:>"O("<>ToString[m1]<>","<>ToString[m2]<>")[1]",Ch[m1_,m2_]:>"O("<>ToString[m1]<>","<>ToString[m2]<>")",Ch1[m1_,m2_][1]:>"O("<>ToString[m1]<>","<>ToString[m2]<>")[1]",Ch1[m1_,m2_]:>"O("<>ToString[m1]<>","<>ToString[m2]<>")"};
nFromCh[{rk_,c11_,c12_,ch2_}]:=-{c11+c12+ch2+rk,c11+ch2,ch2,-c12+ch2};
ChFromn[{n1_,n2_,n3_,n4_}]:=-{n1-n2-n3+n4,n2-n3,n3-n4,n3};
repChn={Ch[m1_,m2_]:>nFromCh[{1,m1,m2,m1 m2}]};
Mat=\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{"0", "2", "0", 
RowBox[{"-", "2"}]},
{
RowBox[{"-", "2"}], "0", "2", "0"},
{"0", 
RowBox[{"-", "2"}], "0", "2"},
{"2", "0", 
RowBox[{"-", "2"}], "0"}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\);
(* {gamd[1],gamd[2],gamd[3],gamd[4]}={{1,0,0,0},{-1,1,0,0},{-1,-1,1,1},{1,0,-1,0}}; *)
McKayDSZ[Nvec_,NNvec_]:=Tr[{Nvec} . Mat . Transpose[{NNvec}]];
McKayVec[Nvec_]:={-Nvec[[1]]+Nvec[[3]],Nvec[[2]]-Nvec[[4]]};
McKayRayEq[Nvec_,{u_,v_}]:=(-1+v) Nvec[[1]]+(-1+u) Nvec[[2]]-(1+v) Nvec[[3]]-(1+u) Nvec[[4]];
McKayIntersectRays[Nvec_,NNvec_]:={-((-((NNvec[[1]]+NNvec[[2]]+NNvec[[3]]+NNvec[[4]]) (Nvec[[1]]-Nvec[[3]]))+(NNvec[[1]]-NNvec[[3]]) (Nvec[[1]]+Nvec[[2]]+Nvec[[3]]+Nvec[[4]]))/((NNvec[[2]]-NNvec[[4]]) (Nvec[[1]]-Nvec[[3]])-(NNvec[[1]]-NNvec[[3]]) (Nvec[[2]]-Nvec[[4]]))),(NNvec[[4]] (Nvec[[1]]+2 Nvec[[2]]+Nvec[[3]])+(NNvec[[1]]+NNvec[[3]]) (Nvec[[2]]-Nvec[[4]])-NNvec[[2]] (Nvec[[1]]+Nvec[[3]]+2 Nvec[[4]]))/(NNvec[[4]] (Nvec[[1]]-Nvec[[3]])+NNvec[[2]] (-Nvec[[1]]+Nvec[[3]])+(NNvec[[1]]-NNvec[[3]]) (Nvec[[2]]-Nvec[[4]]))};
InitialRaysOrigin={{1,1},{1,-1},{-1,-1},{-1,1}};

(* McKayIIinitialRays[L_]:=
Graphics[{McKayIIRay[{1,0,0,0},{1,1},L{-3,5},"\!\(\*SubscriptBox[\(\[Gamma]\), \(1\)]\)"],
McKayIIRay[{0,1,0,0},{1,-1},L{-3,5},"\!\(\*SubscriptBox[\(\[Gamma]\), \(2\)]\)"],
McKayIIRay[{0,0,1,0},{-1,-1},L{-3,5},"\!\(\*SubscriptBox[\(\[Gamma]\), \(3\)]\)"],
McKayIIRay[{0,0,0,1},{-1,1},L{-1,5},"\!\(\*SubscriptBox[\(\[Gamma]\), \(4\)]\)"]}]; *)


(* ::Section:: *)
(*Routines for quiver scattering diagram F0 phase II, m-dependent*)


ClearAll[InitialRaysOrigin];
InitialRaysOrigin[m_]:={{m/2,(1-m)/2},{m/2,1/2 (-1+m)},{-(m/2),1/2 (-1+m)},{-(m/2),(1-m)/2}};
McKayRayEq[{n1_,n2_,n3_,n4_},{u_,v_},m_]:=n1(v+1/2(m-1))+n2(u-m/2)+n3(-v+1/2(m-1))+n4(-u-m/2);
McKayIntersectRays[{n1_,n2_,n3_,n4_},{nn1_,nn2_,nn3_,nn4_},m_]:=If[(-2 n2+2 n4) nn1+(2 n1-2 n3) nn2+(2 n2-2 n4) nn3+(-2 n1+2 n3) nn4==0,{},{(-2 n3 nn1+2 n1 nn3+m (2 n3 nn1-n4 nn1+n1 nn2-n3 nn2-2 n1 nn3+n4 nn3+n2 (-nn1+nn3)+n1 nn4-n3 nn4))/(2 (n4 (nn1-nn3)+n2 (-nn1+nn3)+(n1-n3) (nn2-nn4))),(n4 (nn1-m nn1+2 m nn2+nn3-m nn3)-(-1+m) (n1+n3) (nn2-nn4)+n2 ((-1+m) nn1+(-1+m) nn3-2 m nn4))/(2 (n4 (nn1-nn3)+n2 (-nn1+nn3)+(n1-n3) (nn2-nn4)))}];


(* ::Section:: *)
(*Routines for quiver scattering diagram F0 phase I, m-dependent - comment out if not used !*)


(*ClearAll[InitialRaysOrigin];
McKayVec[Nvec_]:={-Nvec[[1]]+Nvec[[4]],Nvec[[2]]+Nvec[[3]]-2Nvec[[4]]};
McKayRayEq[{n1_,n2_,n3_,n4_},{u_,v_},m_]:=(-n3+2 n4) (-(m/2)-u)+n2 (-(m/2)+u)+n4 (1/2 (-1+m)-v)+n1 (1/2 (-1+m)+v);
McKayDSZ[{n1_,n2_,n3_,n4_},{nn1_,nn2_,nn3_,nn4_}]:=-2 n2 nn1-2 n3 nn1+4 n4 nn1+2 n1 nn2-2 n4 nn2+2 n1 nn3-2 n4 nn3-4 n1 nn4+2 n2 nn4+2 n3 nn4;
McKayIntersectRays[{n1_,n2_,n3_,n4_},{nn1_,nn2_,nn3_,nn4_},m_]:=If[-2 n2 nn1-2 n3 nn1+4 n4 nn1+2 n1 nn2-2 n4 nn2+2 n1 nn3-2 n4 nn3-4 n1 nn4+2 n2 nn4+2 n3 nn4==0,{},{(-2 n4 nn1+2 n1 nn4+m ((n1-n4) (nn2-nn3)+n3 (nn1-nn4)+n2 (-nn1+nn4)))/(2 (2 n4 nn1+n1 nn2-n4 nn2+n1 nn3-n4 nn3-2 n1 nn4+n2 (-nn1+nn4)+n3 (-nn1+nn4))),(-((-1+m) (n1+n4) (nn2+nn3-2 nn4))+n2 ((-1+m) nn1+2 m nn3-nn4-3 m nn4)+(-n3+2 n4) (nn1-m nn1+2 m nn2+nn4-m nn4))/(2 ((n1-n4) (nn2+nn3-2 nn4)-(n3-2 n4) (nn1-nn4)+n2 (-nn1+nn4)))}];
InitialRaysOrigin[m_]:={{m/2,(1-m)/2},{m/2,1/2 (-1-3 m)},{-(m/2),1/2 (-1+m)},{-(1/2),(1-m)/2}}; *)




(* ::Section:: *)
(*Mapping local F1 to local F0*)


F1ToF0[{r_,dH_,dC_,ch2_}]:={r,1/2(dH-dC),dH+dC,ch2};
F0ToF1[{r_,d1_,d2_,ch2_}]:={r,d1+1/2 d2,-d1+1/2 d2,ch2};
repCh1={Ch1[mH_,mC_][1]:>-{1,mH,mC,1/2 (mH^2-mC^2)},Ch1[mH_,mC_]:>{1,mH,mC,1/2 (mH^2-mC^2)},GV1[mH_,mC_,n_]:>{0,mH,mC,n}};
DSZ1[gam1_,gam2_]:=DSZ[F1ToF0[gam1/.repCh1],F1ToF0[gam2/.repCh1]];
ZLV1[gam_,{s_,t_},m_]:=ZLV[F1ToF0[gam/.repCh1],{s,t},m];
Euler1[gam1_,gam2_]:=Euler[F1ToF0[gam1/.repCh1],F1ToF0[gam2/.repCh1]];
Disc1[gam_]:=Disc[F1ToF0[gam/.repCh1]];
SpecFlow1[gam_, {mH_, mC_}] := Evaluate[F0ToF1@SpecFlow[F1ToF0[gam/.repCh1], (F1ToF0@{1, mH, mC, 1})[[2;;3]]]//Simplify];
DimGieseker1[{r_,dH_,dC_,ch2_}]:=DimGieseker[F1ToF0[{r,dH,dC,ch2}]];

InitialPosition1[gam_, m_] := InitialPosition[F1ToF0[gam/.repCh1], m];
IntersectRays1[gam1_,gam2_,m_]:=IntersectRays[F1ToF0[gam1/.repCh1],F1ToF0[gam2/.repCh1],m];
IntersectRays1[gam1_,gam2_,z_,zz_,m_]:=IntersectRays[F1ToF0[gam1/.repCh1],F1ToF0[gam2/.repCh1],z,zz,m];
IntersectRaysSt1[gam1_,gam2_,psi_,m_]:=IntersectRaysSt[F1ToF0[gam1/.repCh1],F1ToF0[gam2/.repCh1],psi,m];

Wall1[gam1_,gam2_,{s_,t_},m_]:=Wall[F1ToF0[gam1],F1ToF0[gam2],{s,t},m];
GenSlope1[gam_,m_]:=GenSlope[F1ToF0[gam],m];

MutateCollection1[Coll_,klist_]:=Module[{Coll0,k,eps},
  (* Coll is a list of Chern vectors, klist a list of {node,\pm 1} *)
  Coll0=If[Length[klist]>1, MutateCollection1[Coll,Drop[klist,-1]], Coll];
  k=Last[klist][[1]]; eps=Last[klist][[2]];
  Table[If[i==k,-Coll0[[k]],Coll0[[i]]+Max[0,eps DSZ1[Coll0[[i]],Coll0[[k]]]]Coll0[[k]]],{i,Length[Coll0]}]]

ExtFromStrong1[Coll_]:=Module[{S,Si},
S=Table[Euler1[Coll[[i]],Coll[[j]]],{i,Length[Coll]},{j,Length[Coll]}];
Si=Inverse[Transpose[S]];
Si . Coll
]
StrongFromExt1[Coll_]:=Module[{S,Si},
S=Table[Euler1[Coll[[j]],Coll[[i]]],{i,Length[Coll]},{j,Length[Coll]}];
Si=Inverse[Transpose[S]];
Si . Coll
]

Options[QuiverDomain1]={"Style"->LightBlue};
QuiverDomain1[Coll_,psi_,m_,OptionsPattern[]]:={RegionPlot[(And@@Table[Re[Exp[-I psi]ZLV1[Coll[[i]],{s, t},m]]<0,{i,Length[Coll]}])&&t>0,{s,-1.5,1.5},{t,0,1},PlotPoints->100,AspectRatio->1,PlotStyle->Flatten[{OptionValue["Style"],Opacity[.5]}]],RegionPlot[SameHalfPlaneQ[Table[ZLV1[Coll[[i]],{s, t},m],{i,Length[Coll]}]],{s,-1.5,1.5},{t,0,1},PlotPoints->100,AspectRatio->1,BoundaryStyle->Directive[Dashed],PlotStyle->Flatten[{OptionValue["Style"],Opacity[.3]}]]};

StabilityRegion1[gam1_,gam2_,m_]:=Module[{},Print["Total charge:", gam1+gam2/.repCh1,", DSZProduct:", DSZ1[gam1,gam2],", Discriminant",Disc1[gam1+gam2]];RegionPlot[Im[ZLV1[gam1/.repCh1,{s,t},m]Conjugate[ZLV1[gam2/.repCh1,{s,t},m]]] DSZ1[gam1,gam2]>0,{s,-2,2},{t,0,2}]];

ScanAllTrees1[{r_,dH_,dC_,ch2_},{s0_,t0_},n_,m_]:=Module[
{Li,phimax},
phimax=1/2(3dH+dC)-r (2s0+m);
If[phimax<=0,Print["There are no trees at this point"];,
Li=ScanConstituents1[{r,dH,dC,ch2},{s0-t0,s0+t0},{n,n},m,phimax];
Print[Li];
ScattSort1[DeleteDuplicatesBy[Flatten[Select[Table[ListStableTrees1[Li[[i]],{s0,t0},m],{i,Length[Li]}],Length[#]>0&],1],ScattGraph1[#,m]&],m]]];

ScanConstituents1[gam_,{smin_,smax_},{n_,np_},m_,phimax_]:=Module[
{eps=.001,Tabs,Tabsp,k,kp,Li,Lip,Lik,Likp,LiCompatible,LiConst,nc,i,ii,mi,mf},
mi=Floor[m];mf=m-mi;
If[m<0 || m>=1/4, Print["For now we assume 0<m<1/4 !"]];
(* allowed initial {position,charge,phi} *)
Tabs=Reverse[Sort[Flatten[{
	Table[{2k-1+eps,Ch1[3k-1,-k],1/2 m},{k,Ceiling[(smin+1)/2],Floor[(smax+1)/2]}],Table[{2k-m,Ch1[3k,-k],m},{k,Ceiling[(smin+m)/2],Floor[(smax+m)/2]}],
Table[{2k+eps,Ch1[3k+1,-k-1],1-m},{k,Ceiling[smin],Floor[smax]}],Table[{2k+1/2-m,Ch1[3k+1,-k],m+1/2},{k,Ceiling[(smin+m-1/2)/2],Floor[(smax+m-1/2)]}],
Table[{2k+1-2m,GV1[0,1,k+1/2],1/2},{k,Ceiling[(smin+2m-1)/2],Floor[(smax+2m-1)/2]}]
},1]]];
Tabsp=Sort[Flatten[{
                 Table[{2k-1/2-m,-Ch1[3k-1,-k],1/2-m},{k,Ceiling[(smin+m+1/2)/2],Floor[(smax+m+1/2)/2]}],Table[{2k-eps,-Ch1[3k,-k],m},{k,Ceiling[smin/2],Floor[smax/2]}],
Table[{2k+1-m,-Ch1[3k+1,-k-1],1-m},{k,Ceiling[(smin+m-1)/2],Floor[(smax+m-1)/2]}],
Table[{2k+1-eps,-Ch1[3k+1,-k],m+1/2},{k,Ceiling[(smin-1)/2],Floor[(smax-1)/2]}],Table[{2k+1-2m,GV1[0,1,k+1/2],1/2},{k,Ceiling[(smin+2m-1)/2],Floor[(smax+2m-1)/2]}]
},1]];
Print["D-brane initial locations:",Tabs];
Print["Dbar-brane initial locations:",Tabsp];
(* choose up to n initial branes, and up to n initial antibranes *)
Li=Flatten[Table[Union[Map[Sort,Tuples[Range[Length[Tabs]],{k}]]],{k,1,n}],1];
Lip=Flatten[Table[Union[Map[Sort,Tuples[Range[Length[Tabsp]],{k}]]],{k,1,np}],1];
(*Print[Li];
Print[Lip]; *)
(* extract compatible choices of branes/antibranes *)
LiCompatible={};
Do[
If[
(* rightmost anti D-brane to the left of rightmost D-brane *)
Tabs[[Last[Li[[i]]],1]]>Tabsp[[First[Lip[[ii]]],1]]&&
(* leftmost anti D-brane to the left of lefttmost D-brane *)
Tabs[[First[Li[[i]]],1]]>Tabsp[[Last[Lip[[ii]]],1]]&&
(* condition on phi with unit multiplicity *)
Sum[Tabs[[Li[[i,j]],3]],{j,Length[Li[[i]]]}]+Sum[Tabsp[[Lip[[ii,jj]],3]],{jj,Length[Lip[[ii]]]}]<=phimax,AppendTo[LiCompatible,{i,ii}]],{i,Length[Li]},{ii,Length[Lip]}];
(*Print["Compatible pairs:", LiCompatible]; *)
nc=0;LiConst={};
Do[
(* Look for compatible choices leading to the correct total charge *)
i=LiCompatible[[l,1]];
ii=LiCompatible[[l,2]];
Lik=Table[{k[j],
If[j>1,If[Li[[i,j]]==Li[[i,j-1]],k[j-1],1],1],phimax-Sum[Tabs[[Li[[i,jj]],3]] k[jj],{jj,1,j-1}]},{j,Length[Li[[i]]]}];
Likp=Table[{kp[j],
If[j>1,If[Lip[[ii,j]]==Lip[[ii,j-1]],kp[j-1],1],1],phimax-Sum[Tabs[[Lip[[ii,jj]],3]]kp[jj],{jj,1,j-1}]},{j,Length[Lip[[ii]]]}];
(*Print[{i,ii,Lik,Likp}];*)
Do[Do[(* loop on kp *)
nc+=1;
If[(Sum[k[j]Tabs[[Li[[i,j]],2]],{j,Length[Li[[i]]]}]+Sum[kp[j]Tabsp[[Lip[[ii,j]],2]],{j,Length[Lip[[ii]]]}]/.repCh1)==gam,AppendTo[LiConst,Join[
Table[k[j]Tabs[[Li[[i,j]],2]],{j,Length[Li[[i]]]}],Table[kp[j]Tabsp[[Lip[[ii,j]],2]],{j,Length[Lip[[ii]]]}]]
]],
##]&@@Likp,
##]&@@Lik;,{l,Length[LiCompatible]}];
Print[Length[LiConst]," possible sets of constituents out of ",nc," trials"];
LiConst
];

ListStableTrees1[LiCh_,{s0_,t0_},m_]:=
Module[{n,LiTree,LiCh1,LiCh2,gam1,gam2,s1,t1,z1,LiStable1,LiStable2},
(* for a given list of k_i Ch[m_i], construct consistent stable trees *)
n=Length[LiCh];
(*Print["**",LiCh,", (s,t)=",{s0,t0}]; *)
If[n==1, LiTree=LiCh,
(* construct possible subsets, avoiding double counting *) 
LiTree={};
Do[
LiCh1=ListSubsetsAny[[n-1,i]]/.n_Integer:>LiCh[[n]];
LiCh2=Complement[Range[n],ListSubsetsAny[[n-1,i]]]/.n_Integer:>LiCh[[n]];
(*Print[LiCh1,",",LiCh2];*)
gam1=Plus@@(LiCh1/.repCh1);
gam2=Plus@@(LiCh2/.repCh1);
(* reorder according to slope *)
If[GenSlope1[gam1,m]>GenSlope1[gam2,m],{gam1,gam2}={gam2,gam1};{LiCh1,LiCh2}={LiCh2,LiCh1}];
If[Wall1[gam1,gam2,{s0,t0},m]DSZ1[gam1,gam2]>0,
z1=IntersectRaysSt[F1ToF0[gam1],F1ToF0[gam2],0,m];
If [Length[z1]>0,{s1,t1}=z1;
(*Print[gam1,",",gam2,",",s1];*)
If[TestBranch[F1ToF0[gam1],s1,m]&&TestBranch[F1ToF0[gam2],s1,m],
(*Print[LiCh1,LiCh2,gam1,gam2,s1,t1]; *)
LiStable1=ListStableTrees1[LiCh1,{s1,t1},m];
LiStable2=ListStableTrees1[LiCh2,{s1,t1},m];
Do[
AppendTo[LiTree,{LiStable1[[j1]],LiStable2[[j2]]}],
{j1,Length[LiStable1]},{j2,Length[LiStable2]}];
]]];,
{i,Length[ListSubsetsAny[[n-1]]]}];
];LiTree
];

ScattCheck1[Tree_,m_]:=(* Check consistency of single tree, returns {charge,{xf,yf}} if tree is consistent, otherwise {charge,{}}; ignore whether leaves have Delta=0 or not *)
Module[{S1,S2,z, r,dH,dC,ch2},
If[!ListQ[Tree]||Length[Tree]>2,
(* tree consists of a single node *)
{r,dH,dC,ch2}=Tree/.repCh1;
z={(dC+3 dH-2 m r)/(4 r),-((5 dC^2+6 dC (dH+m r)+dH (5 dH+2 m r))/(16 r^2))};
(*Print["Initial pt:",z]; *)
{Tree/.repCh1,z},
(* otherwise, check each of the two branches *)
S1=ScattCheck1[Tree[[1]],m];
S2=ScattCheck1[Tree[[2]],m];
If[Length[S1[[2]]]>0&&Length[S2[[2]]]>0,
z=IntersectRays1[S1[[1]]/.repCh1,S2[[1]]/.repCh1,S1[[2]],S2[[2]],m];(*Print[{S1[[1]],S2[[1]],S1[[2]],S2[[2]],z}];*){S1[[1]]+S2[[1]],z},
{S1[[1]]+S2[[1]],{}}]
]];

ScattSort1[LiTree_,m_:1/2]:= (* sort trees by decreasing radius *)
Reverse[Map[#[[2]]&,SortBy[Table[{xytost[ScattCheck1[LiTree[[i]],m][[2]],m][[2]],LiTree[[i]]},{i,Length[LiTree]}],N[First[#]]&]]];


ScattGraphInternal1[Tree_,m_:1/2]:=Module[{S1,S2,TreeNum,sInit,z,Li}(* construct total charge, coordinate of root and list of line segments in (s,t) coordinates, {min(x), max(x)} *), 
If[!ListQ[Tree]||Length[Tree]>2,
	TreeNum=Tree/.repCh1;
sInit=InitialPosition1[TreeNum,m];{Tree,{sInit,-1/2 sInit(sInit+m)},{}},
S1=ScattGraphInternal1[Tree[[1]],m];
S2=ScattGraphInternal1[Tree[[2]],m];
z=IntersectRays1[S1[[1]]/.repCh1,S2[[1]]/.repCh1,m];
If[Length[z]==0,Print["Illegal tree"],
Li={S1[[3]],S2[[3]],Arrow[{S1[[2]],z}],Arrow[{S2[[2]],z}]};
{S1[[1]]+S2[[1]],z,Li}]]];

ScattGraph1[Tree_,m_:1/2]:=Module[{T,LiArrows,LiVertex},
(* extracts list of vertices in (x,y) plane and adjacency matrix *)
T=ScattGraphInternal1[Tree,m];
LiArrows=Cases[Flatten[T[[3]]],x_Arrow]/.Arrow[x_]:>x;
LiVertex=Union[Flatten[LiArrows,1]];
{LiVertex,Table[If[i!=j,Sign[Count[LiArrows,{LiVertex[[i]],LiVertex[[j]]}]],0],{i,Length[LiVertex]},{j,Length[LiVertex]}]}];

ContractInitialRays1[Trees_,m_:1/2]:=Trees /.repCh1/.{{r_Integer,dH_Integer,dC_Integer,ch2_},{rr_Integer,ddH_Integer,ddC_Integer,cch2_}}:>{r+rr,dH+ddH,dC+ddC,ch2+cch2}/;InitialPosition1[{r,dH,dC,ch2},m]==InitialPosition1[{rr,ddH,ddC,cch2},m]/.{{r_Integer,dH_Integer,dC_Integer,ch2_}:>If[r==0,GV[dH,dC,ch2],If[Disc1[{r,dH,dC,ch2}]==0,If[r>0,r Ch1[dH/r,dC/r],-r Ch1[dH/r,dC/r][1]],{r,dH,dC,ch2}]]};

ScattIndex1[TreeList_]:=Table[
(* compute index for each tree in the list; do not trust the result if internal lines have non-primitive charges *)
Simplify[Times@@ScattIndexInternal1[TreeList[[i]]][[2]]],{i,Length[TreeList]}];

GCD1[{r_,dH_,dC_,ch2_}]:=Module[{d},d=GCD[r,dH,dC];
If[EvenQ[(dH+dC-2ch2)/d],d,If[EvenQ[d],d/2,1]]]

ScattIndexInternal1[Tree_]:=Module[{S1,S2,g1,g2,kappa,Li},
(* compute {total charge, list of Kronecker indices associated to each vertex *)
If[!ListQ[Tree]||Length[Tree]>2,{Tree,{1}},
S1=ScattIndexInternal1[Tree[[1]]]/.repCh1;
S2=ScattIndexInternal1[Tree[[2]]]/.repCh1;
(*Print[S1[[1]],",",S2[[1]]];*)
g1=GCD1[S1[[1]]];g2=GCD1[S2[[1]]];
kappa=Abs[DSZ1[S1[[1]],S2[[1]]]]/g1/g2;
Li=Join[S1[[2]],S2[[2]]];
AppendTo[Li,Subscript[Kr, kappa][Min[g1,g2],Max[g1,g2]]];
If[GCD1[S1[[1]]+S2[[1]]]!=1,Print["Beware, non-primitive state"]];
{S1[[1]]+S2[[1]],Li}]];

Options[ScattIndexImprovedInternal1] = {"Debug"->False};

ScattIndexImproved1[TreeList_, opt: OptionsPattern[]]:=Table[
	(* compute index for each tree in the list *)
	Simplify[FOmbToOm[Last@ScattIndexImprovedInternal1[TreeList[[i]], opt][[2]]]],{i,Length[TreeList]}];

ScattIndexImprovedInternal1[Tree_, opt: OptionsPattern[]]:=Module[{S1,S2,g1,g2,gFinal, kappa,Li, tem, repOmAttb, rrr},
(* compute {total charge, list of Kronecker indices associated to each vertex *)
	If[!ListQ[Tree]||Length[Tree]>2,{Tree,{Join[{1}, Table[(y-y^-1)/(j(y^j-y^-j)), {j, 2, GCD1@(Tree/.repCh1)}]]}},
	If[OptionValue["Debug"], Print["Calling with args: ", Tree[[1]], "  |  ", Tree[[2]]]];
    S1=ScattIndexImprovedInternal1[Tree[[1]], opt]/.repCh1;
	S2=ScattIndexImprovedInternal1[Tree[[2]], opt]/.repCh1;
If[OptionValue["Debug"], Print["S1 is: ", S1, "   S2 is: ", S2]];
	g1=GCD1@S1[[1]];g2=GCD1@S2[[1]];
    gFinal = GCD1@(S1[[1]]+S2[[1]]);
	kappa=Abs[DSZ1[S1[[1]],S2[[1]]]]/g1/g2;
	Li=Join[S1[[2]],S2[[2]]];
If[OptionValue["Debug"], Print["Li is: ", Li, "  g1 is: ", g1, "  g2 is: ", g2, "  gFinal is: ", gFinal]];
AppendTo[Li,
repOmAttb = Join[
Table[CoulombHiggs`OmAttb[{P, 0}, y_]->Last[S1[[2]]][[P]], {P, 1, g1}],
Table[CoulombHiggs`OmAttb[{0, Q}, y_]->Last[S2[[2]]][[Q]], {Q, 1, g2}]
];
If[OptionValue["Debug"], Print["repOmAttb is: ", repOmAttb]];
tem = Table[
rrr = If[And@@(IntegerQ/@{P g1/gFinal, P g2/gFinal}),CoulombHiggs`FlowTreeFormulaRat[{{0, kappa}, {-kappa, 0}}, {g2, -g1}, {P g1/gFinal, P g2/gFinal}, y], 0];
Simplify[
rrr
/.repOmAttb
/.{CoulombHiggs`OmAttb[{p_, q_}, y___]:>0/;p>1||q>1||p q !=0}
],
{P, 1, gFinal}
];
If[OptionValue["Debug"], Print["tem is: ", tem]];tem
];
	(*If[GCD@@(S1[[1]]+S2[[1]])!=1,Print["Beware, non-primitive state"]];*)
	{S1[[1]]+S2[[1]],Li}]];
	
ConstructLVDiagram1[smin_,smax_,phimax_,Nm_,m_,ListRays0_]:=Module[
{eps=.001,mi,mf,Inter,ListInter,ListRays,ListNewRays,kappa,KTab},
mi=Floor[m];mf=m-mi;
(* initial rays {charge, {x,y}, parent1, parent2,n1,n2 } *)
If[ListRays0=={},
ListRays=Flatten[{
Table[{Ch1[3k,-k],{2k-m/2,-2k^2},0,0,0,0,0},{k,Ceiling[(smin+m/2)/2],Floor[(smax+m/2)/2]}],
Table[{Ch1[3k,-k][1],{2k-m/2,-2k^2},0,0,0,0,0},{k,Ceiling[(smin+m/2)/2],Floor[(smax+m/2)/2]}],
Table[{Ch1[3k+1,-k-1],{2k+(1-m)/2,(m-1)/4-2k^2-k m-2k (1-m)/2},0,0,0,0,0},{k,Ceiling[(smin+(m-1)/2)/2],Floor[(smax+(m-1)/2)/2]}],
Table[{Ch1[3k+1,-k-1][1],{2k+(1-m)/2,(m-1)/4-2 k^2-k m-2 k (1-m)/2},0,0,0,0,0},{k,Ceiling[(smin+(m-1)/2)/2],Floor[(smax+(m-1)/2)/2]}],
Table[{Ch1[3k+1,-k],{2k+3/4-m/2,-m/8-5/16-2k^2-k m-2 k(3/4-m/2)},0,0,0,0,0},{k,Ceiling[(smin-3/4+m/2)/2],Floor[(smax-3/4+m/2)/2]}],
Table[{Ch1[3k+1,-k][1],{2k+3/4-m/2,-m/8-5/16-2 k^2-k m-2 k(3/4-m/2)},0,0,0,0,0},{k,Ceiling[(smin-3/4+m/2)/2],Floor[(smax-3/4+m/2)/2]}],
Table[{Ch1[3k+2,-k-1],{2k+5/4-m/2,m/8-13/16-2 k^2-k m-2k(5/4-m/2)},0,0,0,0,0},{k,Ceiling[(smin-5/4+m/2)/2],Floor[(smax-5/4+m/2)/2]}],
Table[{Ch1[3k+2,-k-1][1],{2k+5/4-m/2,m/8-13/16-2k^2-k m-2 k(5/4-m/2)},0,0,0,0,0},{k,Ceiling[(smin-5/4+m/2)/2],Floor[(smax-5/4+m/2)/2]}]
},1]/.repCh1;
   ListInter={};,
(* If list of rays is already provided *)
	ListRays=ListRays0;
	ListInter=Select[Table[{ListRays[[i,3]],ListRays[[i,4]]},{i,Length[ListRays]}],First[#]>0&]];
While[True,
ListNewRays={};
       Monitor[ Do[
If[  !MemberQ[ListInter,{i,j}],
AppendTo[ListInter,{i,j}];
 kappa=DSZ1[ListRays[[i,1]],ListRays[[j,1]]];
If[kappa!=0,Inter=IntersectRaysNoTest1[ListRays[[i,1]],ListRays[[j,1]],ListRays[[i,2]],ListRays[[j,2]],m];
If[Inter!={},
KTab=KroneckerDims[Abs[kappa],Nm];
Do[If[CostPhi1[KTab[[k,1]] ListRays[[i,1]]+KTab[[k,2]]ListRays[[j,1]],Inter[[1]],m]<=phimax,
AppendTo[ListNewRays,{KTab[[k,1]]ListRays[[i,1]]+KTab[[k,2]]ListRays[[j,1]],Inter,i,j,KTab[[k,1]],KTab[[k,2]]}]],{k,Length[KTab]}]
]]]
,{i,Length[ListRays]},{j,i+1,Length[ListRays]}],{i,j}];
If[ListNewRays=={},Break[],
Print["Adding ",Length[ListNewRays], " rays, "];
ListRays=Flatten[{ListRays,ListNewRays},1];
]];
Print[Length[ListRays], " in total."];
ListRays];

CostPhi1[{r_,dH_,dC_,ch2_},s_,m_]:=(dC+3dH)/2-r(2s+m);

IntersectRaysNoTest1[{r_,dH_,dC_,ch2_},{rr_,ddH_,ddC_,cch2_},z_,zz_,m_]:=
(* returns (x,y) coordinate of intersection point of two rays, or {} if they don't intersect *)
(* here do not test if DSZ<>0, and require strictly in future of z and zz *)
Module[{zi},zi={(2 (cch2 r-ddC m r-ddH m r-ch2 rr+dC m rr+dH m rr))/(ddC r+3 ddH r-(dC+3 dH) rr),(ch2 (ddC+3 ddH)-cch2 (dC+3 dH)-2 dC ddH m+2 ddC dH m)/(2 ddC r+6 ddH r-2 (dC+3 dH) rr)};
If[(zi-z) . {-2r,(dC+3dH)/2-m r}>0&&(zi-zz) . {-2rr,(ddC+3ddH)/2-m rr}>0,zi,{}]];

(* Extract all trees leading to a ray with charge {r,d,chi} *)
LVTreesFromListRays1[ListRays_,{r_,dH_,dC_,ch2_},m_]:=Module[{Lipos,Div,LiTrees},
Div=Divisors[GCD1[{r,dH,dC,ch2}]];
Lipos=Flatten[Join[Table[Position[ListRays,{r,dH,dC,ch2}/k],{k,Div}]],1];
If[Lipos=={},
Print["No such dimension vector in the list"],
LiTrees=(GCD1[{r,dH,dC,ch2}]/GCD1[ListRays[[#,1]]])TreeFromListRays[ListRays,#]&/@First[Transpose[Lipos]];
ScattSort1[DeleteDuplicatesBy[SortBy[LiTrees,Length[TreeConstituents[#]]&],ScattGraph1[#,m]&],m]
]];

ChernToCh1[f_]:=f/.{{r_Integer,dH_Integer,dC_Integer,ch2_}:>If[r==0,GV1[dH,dC,ch2],If[Disc1[{r,dH,dC,ch2}]==0,If[r>0,r Ch1[dH/r,dC/r],-r Ch1[dH/r,dC/r][1]],{r,dH,dC,ch2}]]};	
	
	
(* precomputed sequences for m=1/4 *)
LVTreesF1[{1,2,1,3/2}]={{{Ch1[-2,0][1],Ch1[-2,1]},{Ch1[2,-1],{Ch1[1,-1][1],Ch1[1,0]}}}};
LVTreesF1[{1,2,0,2}]={{Ch1[2,-1],{Ch1[1,-1][1],Ch1[1,0]}}};
LVTreesF1[{1,0,1,-1/2}]={{Ch1[0,0],{Ch1[-2,0][1],Ch1[-2,1]}}};
LVTreesF1[{-1,0,1,1/2}]={{Ch1[0,0][1],{Ch1[1,-1][1],Ch1[1,0]}}};
LVTreesF1[{1,1,1,0}]={{Ch1[1,0],{Ch1[-2,0][1],Ch1[-2,1]}}};
LVTreesF1[{-1,2,1,-(3/2)}]={{Ch1[-2,0][1],{Ch1[1,-1][1],Ch1[1,0]}}};
LVTreesF1[{0,1,-1,0}]={{{-1,0,0,0},{1,1,-1,0}}};
LVTreesF1[{0,2,0,0}]={{Ch1[1,0],Ch1[-1,0][1]}};
LVTreesF1[{0,2,-1,1/2}]={{Ch1[1,0],Ch1[-1,1][1]}};
LVTreesF1[{0,1,0,1/2}]={{Ch1[0,0][1],Ch1[1,0]}};
LVTreesF1[{0,3,-1,0}]={{Ch1[-1,0][1],{Ch1[1,0],GV[1,-1,0]}}};
LVTreesF1[{0,4,-2,0}]={{Ch1[-2,1][1],Ch1[2,-1]},{Ch1[-1,0][1],{Ch1[1,0],{2 Ch1[0,0][1],2 Ch1[1,-1]}}}};
LVTreesF1[{1,0,0,-1}]={{{Ch1[-2,0][1],Ch1[-2,1]},Ch1[0,-1]},{Ch1[-2,0][1],2 Ch1[-1,0]}};
LVTreesF1[{1,0,0,-2}]={{{Ch1[-5,1][1],Ch1[-5,2]},{Ch1[-2,1][1],2 Ch1[-1,0]}},{{Ch1[-3,1][1],Ch1[-2,0]},{Ch1[-1,0],{Ch1[-2,0][1],Ch1[-2,1]}}},{Ch1[-1,0],{Ch1[-3,1][1],Ch1[-2,1]}}};
	


End[];
EndPackage[]
