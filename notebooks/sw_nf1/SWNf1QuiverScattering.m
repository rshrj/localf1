(* ::Package:: *)

(* ::Title:: *)
(*SWNf1QuiverScattering - Generic quiver scattering diagram machinery*)


Print["SWNf1QuiverScattering 1.0, June 2025 - Generic quiver scattering diagram machinery"];


BeginPackage["SWNf1QuiverScattering`"]
Needs["CoulombHiggs`"]


(* ===== Quiver-specific symbols: must be defined by the user ===== *)

Mat::usage = "DSZ matrix of the quiver (n\[Times]n antisymmetric integer matrix for an n-node quiver). \
Must be defined before using any McKay* or ConstructMcKayDiagram routines.";

McKayDSZ::usage = "McKayDSZ[Nvec_,NNvec_] computes the antisymmetric DSZ pairing between two \
dimension vectors, typically as Tr[{Nvec}.Mat.Transpose[{NNvec}]]. \
Must be defined for the specific quiver (2-arg form, no mass parameter).";

McKayVec::usage = "McKayVec[Nvec_] computes the positive direction vector along the scattering ray \
associated to dimension vector Nvec. Must be defined for the specific quiver.";

McKayRayEq::usage = "McKayRayEq[Nvec_,{u_,v_}] gives the linear form in (u,v) that vanishes on \
the scattering ray for dimension vector Nvec (no mass). \
McKayRayEq[Nvec_,{u_,v_},m_] is the mass-dependent version. \
At least the 3-arg version must be defined for the specific quiver.";

McKayIntersectRays::usage = "McKayIntersectRays[Nvec_,NNvec_,m_] (3-arg) returns the \
intersection point (u,v) of two rays for the given mass m, or {} if collinear. \
McKayIntersectRays[Nvec_,NNvec_,z_,zz_,m_] (5-arg) additionally tests that the intersection \
lies forward of z along Nvec and forward of zz along NNvec. \
The 3-arg form must be defined for the specific quiver; the 5-arg form is provided generically.";

InitialRaysOrigin::usage = "InitialRaysOrigin[m_] returns the list of initial (u,v) positions \
for the simple-root rays, one per quiver node, as a function of the mass parameter m. \
Must be defined for the specific quiver.";

repChn::usage = "Replacement rule converting symbolic charge labels (e.g. Ch[m1,m2]) to \
dimension vectors for the specific quiver. Must be defined for the specific quiver.";

repChO::usage = "Replacement rule converting symbolic charge labels (e.g. Ch[m1,m2]) to \
human-readable display strings. Must be defined for the specific quiver.";


(* ===== Symbolic helpers ===== *)

\[Gamma]::usage = "Subscript[\[Gamma],i] denotes the basis charge vector for the i-th node of the quiver.";

McKayrep::usage = "McKayrep is a replacement rule that converts a dimension vector {n1,n2,...} \
into the symbolic linear combination n1 \[Gamma]_1 + n2 \[Gamma]_2 + ... . \
Depends on Mat being defined (uses Length[Mat] for node count).";

McKayrepi::usage = "McKayrepi is a replacement rule that converts Subscript[\[Gamma],i] back to \
IdentityMatrix[Length[Mat]][[i]].";

TreeToGamma::usage = "TreeToGamma[Tree_] replaces dimension vectors of the correct length \
inside Tree by symbolic linear combinations of Subscript[\[Gamma],i].";


(* ===== Combinatorial utilities ===== *)

BinarySplits::usage = "BinarySplits[Nvec_] gives the list of all dimension vectors smaller \
than or equal to Nvec/2 (first half of all non-trivial binary splits of Nvec).";

ListAllPartitions::usage = "ListAllPartitions[gam_] returns the list of unordered partitions \
of the positive integer vector gam as a sum of positive integer vectors.";

QDeformedFactorial::usage = "QDeformedFactorial[n_,y_] gives the q-deformed factorial [n]_y!.";

StackInvariant::usage = "StackInvariant[Mat_,Cvec_,Nvec_,y_] gives the stack invariant of a \
quiver with DSZ matrix Mat, FI parameters Cvec, and dimension vector Nvec, \
computed using Reineke's formula. Assumes no oriented loops.";

RationalInvariant::usage = "RationalInvariant[Mat_,Cvec_,Nvec_,y_] computes the rational \
invariant of a quiver using Reineke's formula.";

HiggsBranchFormula::usage = "HiggsBranchFormula[Mat_,Cvec_,Nvec_] computes the refined Higgs \
branch index of a quiver with DSZ matrix Mat, FI parameters Cvec, and dimension vector Nvec, \
using Reineke's formula. Assumes no oriented loops.";

EvaluateKronecker::usage = "EvaluateKronecker[f_] replaces each Subscript[Kr,\[Kappa]][g1_,g2_] \
with the refined index of the Kronecker quiver with \[Kappa] arrows and dimension vector {g1,g2}, \
computed via HiggsBranchFormula.";

Kr::usage = "Subscript[Kr,m][p,q] denotes the index of the Kronecker quiver with m arrows and \
dimension vector (p,q). Use EvaluateKronecker to substitute numerical values.";

KroneckerDims::usage = "KroneckerDims[m_,Nn_] gives the list of coprime dimension vectors {n1,n2} \
that are populated (i.e. have non-negative discriminant m n1 n2 - n1^2 - n2^2 + 1 >= 0) \
for the Kronecker quiver with m arrows, with 0 <= n1, n2 <= Nn.";

TreeConstituents::usage = "TreeConstituents[Tree_] gives the flattened list of all leaf nodes \
(constituents) of Tree.";

TreeHue::usage = "TreeHue[n_,i_] returns the color (a Blend of Brown and Green) for the i-th \
tree among a list of n trees.";

TreeFromListRays::usage = "TreeFromListRays[ListRays_,k_] reconstructs the binary tree leading \
to the k-th ray in a ListRays output from ConstructMcKayDiagram.";

FOmbToOm::usage = "FOmbToOm[OmbList_] converts a list of rational (Ombar-type) indices \
{Ombar(1), Ombar(2), ...} to the integer DT index Om(n) via the plethystic inverse relation.";


(* ===== McKay generic machinery ===== *)

McKayRay::usage = "McKayRay[Nvec_,{u_,v_},{k1_,k2_},tx_] produces an Arrow graphics primitive \
from {u,v} + k1*McKayVec[Nvec] to {u,v} + k2*McKayVec[Nvec], with a Text label tx \
placed at the tip.";

McKayIntersectRaysNoTest::usage = "McKayIntersectRaysNoTest[Nvec_,NNvec_,z_,zz_,m_] returns \
the intersection point (u,v) of two rays if it lies strictly forward of z along Nvec and \
strictly forward of zz along NNvec, without testing non-vanishing of the DSZ product.";

McKayInitialRays::usage = "McKayInitialRays[L_,m_] draws the Length[Mat] initial rays in the \
(u,v) plane, each rescaled by factor L, starting from InitialRaysOrigin[m].";

McKayScattDiag::usage = "McKayScattDiag[TreeList_,m_] draws the McKay scattering diagram in \
the (u,v) plane for each tree in TreeList, using colors from TreeHue.";

McKayScattDiagInternal::usage = "McKayScattDiagInternal[Tree_,m_] computes \
{total charge, root coordinate {u,v}, list of Arrow segments}; used internally by McKayScattDiag.";

McKayScattCheck::usage = "McKayScattCheck[Tree_,m_] returns {total charge, {u,v}} of the root \
vertex if Tree is a consistent scattering tree, otherwise {total charge, {}}.";

McKayScattIndex::usage = "McKayScattIndex[TreeList_] computes the wall-crossing index (product of \
Kronecker quiver indices over all vertices) for each tree in TreeList. Do not trust the result \
if any internal line carries a non-primitive charge; use McKayScattIndexImproved instead.";

McKayScattIndexInternal::usage = "McKayScattIndexInternal[Tree_] returns \
{total charge, list of Subscript[Kr,\[Kappa]][g1,g2] factors at each vertex of Tree}.";

McKayScattIndexImproved::usage = "McKayScattIndexImproved[TreeList_,opt___] computes the \
wall-crossing index for each tree in TreeList, correctly handling non-primitive internal states \
using CoulombHiggs`FlowTreeFormulaRat.";

McKayScattIndexImprovedInternal::usage = "McKayScattIndexImprovedInternal[Tree_,opt___] \
computes the index for a single Tree, correctly handling non-primitive internal states. \
Supports option \"Debug\"->True for verbose output.";

McKayScattGraph::usage = "McKayScattGraph[Tree_,m_] extracts the geometric graph structure of \
Tree: returns {list of vertex coordinates, adjacency matrix}.";

McKayListAllTrees::usage = "McKayListAllTrees[Nvec_] generates all binary trees with leaves at \
simple-root charges and non-zero DSZ pairing at every internal vertex, with total dimension \
vector Nvec. Requires Mat and McKayDSZ to be defined.";

McKayScanAllTrees::usage = "McKayScanAllTrees[Nvec_,m_] generates all consistent scattering \
trees with total dimension vector Nvec and mass parameter m, removing duplicates by graph \
isomorphism, sorted by decreasing size.";

ConstructMcKayDiagram::usage = "ConstructMcKayDiagram[phimax_,Nm_,m_,ListRays_] constructs the \
quiver scattering diagram by iteratively intersecting rays up to total height phimax, with \
Kronecker multiplicity cutoff Nm. \
If ListRays is {}, starts from the Length[Mat] initial simple-root rays; otherwise resumes \
from the given list. Supports CheckAbort for graceful interruption. Returns a list of entries \
{charge, {u,v}, parent1Index, parent2Index, mult1, mult2, level}.";

McKayTreesFromListRays::usage = "McKayTreesFromListRays[ListRays_,Nvec_,m_] extracts the list \
of geometrically distinct binary trees with total dimension vector Nvec from a ray list \
produced by ConstructMcKayDiagram.";

McKayPlotDiagram::usage = "McKayPlotDiagram[ListRays_,m_] plots the quiver scattering diagram \
from a ray list produced by ConstructMcKayDiagram. \
Options: \"RayLength\"->5 (arrow scale factor), \"UseMaTeX\"->True (use MaTeX for labels), \
\"FontSize\"->11, \"LabelOffset\"->0.5, ImageSize->500.";


Begin["`Private`"]


(* ===== Symbolic helpers ===== *)

McKayrep := {Nvec_?(VectorQ[#,IntegerQ]&&Length[#]==Length[Mat]&) :>
  Sum[Nvec[[i]] Subscript[\[Gamma],i], {i,Length[Mat]}]};

McKayrepi := {Subscript[\[Gamma],i_] :> IdentityMatrix[Length[Mat]][[i]]};

TreeToGamma[Tree_] := Tree /. {Nvec_List :>
  Sum[Nvec[[i]] Subscript[\[Gamma],i], {i,Length[Nvec]}] /; Length[Nvec]==Length[Mat]};


(* ===== Combinatorial utilities ===== *)

BinarySplits[Nvec_] := Module[{Li,Li1,Nl},
  If[Plus@@Nvec==1, Li1=Nvec,
    Li=Drop[Drop[Flatten[
      Table[Table[Nl[i],{i,Length[Nvec]}],
        Evaluate[Sequence@@Table[{Nl[i],0,Nvec[[i]]},{i,Length[Nvec]}]]],
      Length[Nvec]-1],1],-1];
    Li1=Take[Li,Ceiling[Length[Li]/2]]];
  Li1];

ListAllPartitions[gam_] := Module[{m,unit,Li},
  If[Plus@@gam==1, {{gam}},
    m=Max[Select[Range[Length[gam]],gam[[#]]>0&]];
    unit=Table[If[i==m,1,0],{i,Length[gam]}];
    Li=ListAllPartitions[gam-unit];
    Union[Map[Sort,
      Union[Flatten[
        Table[Union[Flatten[{
          {Flatten[{Li[[i]],{unit}},1]},
          Table[Table[If[j==k,Li[[i,j]]+unit,Li[[i,j]]],
            {j,Length[Li[[i]]]}],
            {k,Length[Li[[i]]]}]},1]],
          {i,Length[Li]}],1]],
      1]]]];

QDeformedFactorial[n_,y_] :=
  If[n<0,
    Print["QDeformedFactorial[n,y] is defined only for n>=0"],
    If[n==0, 1, (y^(2n)-1)/(y^2-1) QDeformedFactorial[n-1,y]]];

StackInvariant[Mat_,Cvec_,Nvec_,y_] := Module[{m,pa,Cvec0},
  m=Length[Nvec];
  If[Max[Nvec]<0, Print["StackInvariant: The dimension vector must be positive!"]];
  If[Plus@@Nvec==0, Return[0]];
  If[Plus@@Nvec==1, Return[1]];
  Cvec0=Cvec-(Plus@@(Nvec Cvec))/(Plus@@Nvec);
  pa=Flatten[Map[Permutations,ListAllPartitions[Nvec]],1];
  (-y)^(Sum[-Max[Mat[[k,l]],0]Nvec[[k]]Nvec[[l]],{k,m},{l,m}]-1+Plus@@Nvec)
  (y^2-1)^(1-Plus@@Nvec)
  Sum[
    If[(Length[pa[[i]]]==1)||
        And@@Table[Sum[Cvec0[[k]]pa[[i,a,k]],{a,b},{k,m}]>0,
          {b,Length[pa[[i]]]-1}],
      (-1)^(Length[pa[[i]]]-1)
      y^(2 Sum[Max[Mat[[l,k]],0]pa[[i,a,k]]pa[[i,b,l]],
          {a,1,Length[pa[[i]]]},{b,a,Length[pa[[i]]]},{k,m},{l,m}])/
      Product[QDeformedFactorial[pa[[i,j,k]],y],
        {j,1,Length[pa[[i]]]},{k,m}],
      0],
    {i,Length[pa]}]];

RationalInvariant[Mat_,Cvec_,Nvec_,y_] := Module[{Li,gcd},
  gcd=GCD@@Nvec;
  Li=Flatten[Map[Permutations,ListAllPartitions[{gcd}]],1];
  Sum[
    Product[StackInvariant[Mat,Cvec,Nvec Li[[i,j,1]]/gcd,y],{j,Length[Li[[i]]]}]/
    Length[Li[[i]]]/(y-1/y)^(Length[Li[[i]]]-1),
    {i,Length[Li]}]];

HiggsBranchFormula[Mat_,Cvec_,Nvec_] := Module[{Cvec0},
  If[Max[Nvec]<0, Print["HiggsBranchFormula: The dimension vector must be positive!"]];
  If[Plus@@Nvec==0, Return[0]];
  If[Plus@@Nvec==1, Return[1]];
  Cvec0=Cvec-(Plus@@(Nvec Cvec))/(Plus@@Nvec);
  DivisorSum[GCD@@Nvec,
    (y-1/y)/(y^#-1/y^#)/# MoebiusMu[#]
    RationalInvariant[Mat,Cvec0,Nvec/#,y^#]&]];

EvaluateKronecker[f_] :=
  f /. Subscript[Kr,kappa_][g1_,g2_] :>
    Simplify[HiggsBranchFormula[{{0,kappa},{-kappa,0}},{1/g1,-1/g2},{g1,g2}]];

ClearAll[KroneckerDims];
KroneckerDims[m_,Nn_] := KroneckerDims[m,Nn] = Module[{Ta={}},
  Do[If[m n1 n2-n1^2-n2^2+1>=0&&GCD[n1,n2]==1,
    AppendTo[Ta,{n1,n2}]],
    {n1,0,Nn},{n2,0,Nn-n1}];
  Drop[Ta,2]];

TreeConstituents[Tree_] :=
  If[!ListQ[Tree]||Length[Tree]>2, {Tree},
    Join[TreeConstituents[Tree[[1]]],TreeConstituents[Tree[[2]]]]];

TreeHue[n_,i_] := Blend[{Brown,Green},i/n];

TreeFromListRays[ListRays_,k_] :=
  If[ListRays[[k,3]]==0,
    ListRays[[k,1]],
    {ListRays[[k,5]] TreeFromListRays[ListRays,ListRays[[k,3]]],
     ListRays[[k,6]] TreeFromListRays[ListRays,ListRays[[k,4]]]}];

FOmbToOm[OmbList_] := Module[{n},
  If[Length[OmbList]<2,
    First@OmbList,
    n=Length[OmbList];
    DivisorSum[n,
      (MoebiusMu[#](y-y^-1)/(#(y^#-y^-#))(OmbList[[n/#]]/.{y->y^#}))&]]];


(* ===== Generic McKay machinery ===== *)

McKayRay[Nvec_,{u_,v_},{k1_,k2_},tx_] :=
  {Arrow[{{u,v}+k1 McKayVec[Nvec],{u,v}+k2 McKayVec[Nvec]}],
   Text[tx,{u,v}+(k2+.1)McKayVec[Nvec]]};

(* 5-arg form: test forward from z and zz; requires McKayDSZ and 3-arg McKayIntersectRays *)
McKayIntersectRays[Nvec_,NNvec_,z_,zz_,m_] :=
  Module[{zi},
    If[McKayDSZ[Nvec,NNvec]!=0,
      zi=McKayIntersectRays[Nvec,NNvec,m];
      If[(zi-z).McKayVec[Nvec]>=0&&(zi-zz).McKayVec[NNvec]>=0,zi,{}]]];

McKayIntersectRaysNoTest[Nvec_,NNvec_,z_,zz_,m_] :=
  Module[{zi},
    zi=Quiet[McKayIntersectRays[Nvec,NNvec,m]];
    If[!FreeQ[zi,DirectedInfinity|Indeterminate],{},
      If[(zi-z).McKayVec[Nvec]>0&&(zi-zz).McKayVec[NNvec]>0,zi,{}]]];

McKayInitialRays[L_,m_] :=
  Graphics[Table[{Thick,Red,
    McKayRay[IdentityMatrix[Length[Mat]][[i]],
      InitialRaysOrigin[m][[i]],
      L{-1,1},
      "\!\(\*SubscriptBox[\(\[Gamma]\), \("<>ToString[i]<>"\)]\)"]},
    {i,Length[Mat]}]];

McKayScattDiag[TreeList_,m_] :=
  Module[{T,TNum,Diagram},
    Diagram={};
    Do[
      T=McKayScattDiagInternal[TreeList[[i]],m];
      TNum=T/.repChn;
      AppendTo[Diagram,TreeHue[Length[TreeList],i]];
      AppendTo[Diagram,T[[3]]];
      AppendTo[Diagram,
        Arrow[{TNum[[2]],TNum[[2]]+McKayVec[TNum[[1]]]/GCD@@(TNum[[1]])}]];
      AppendTo[Diagram,
        Text[TNum[[1]]/.McKayrep,TNum[[2]]+1.2 McKayVec[TNum[[1]]]/GCD@@(TNum[[1]])]];
      ,{i,Length[TreeList]}];
    Graphics[Diagram]];

McKayScattDiagInternal[Tree_,m_] :=
  Module[{S1,S2,TreeNum,z,Li},
    If[!ListQ[Tree]||Length[Tree]>2,
      (* leaf node: single simple-root charge *)
      TreeNum=Tree/.repChn; z={};
      If[Length[Position[TreeNum,n_Integer/;n>0]]==1,
        z=InitialRaysOrigin[m][[Tr[Position[TreeNum,n_Integer/;n>0]]]]-
          2 McKayVec[TreeNum/GCD@@TreeNum];,
        Print["Illegal endpoint"]];
      {Tree,z,{Text[Tree/.repChO/.McKayrep,z-.2 McKayVec[TreeNum]]}},
      (* internal vertex *)
      S1=McKayScattDiagInternal[Tree[[1]],m];
      S2=McKayScattDiagInternal[Tree[[2]],m];
      z=McKayIntersectRays[S1[[1]]/.repChn,S2[[1]]/.repChn,m];
      If[Length[z]==0,
        Print["Illegal tree"],
        Li=S1[[3]];
        AppendTo[Li,S2[[3]]];
        AppendTo[Li,Arrow[{S1[[2]],z}]];
        AppendTo[Li,Arrow[{S2[[2]],z}]];
        {S1[[1]]+S2[[1]],z,Li}]]];

McKayScattCheck[Tree_,m_] :=
  Module[{S1,S2,TreeNum,z},
    If[!ListQ[Tree]||Length[Tree]>2,
      TreeNum=Tree/.repChn;
      If[Length[Position[TreeNum,n_Integer/;n>0]]==1,
        z=InitialRaysOrigin[m][[Tr[Position[TreeNum,n_Integer/;n>0]]]]-
          20 McKayVec[TreeNum/GCD@@TreeNum],
        z=-20 McKayVec[TreeNum]];
      {Tree/.repChn,z},
      S1=McKayScattCheck[Tree[[1]],m];
      S2=McKayScattCheck[Tree[[2]],m];
      If[Length[S1[[2]]]>0&&Length[S2[[2]]]>0,
        {S1[[1]]+S2[[1]],
          McKayIntersectRays[S1[[1]]/.repChn,S2[[1]]/.repChn,S1[[2]],S2[[2]],m]},
        {S1[[1]]+S2[[1]],{}}]]];

McKayScattIndex[TreeList_] :=
  Table[Simplify[Times@@McKayScattIndexInternal[TreeList[[i]]][[2]]],
    {i,Length[TreeList]}];

McKayScattIndexInternal[Tree_] :=
  Module[{S1,S2,g1,g2,kappa,Li},
    If[!ListQ[Tree]||Length[Tree]>2,
      {Tree,{1}},
      S1=McKayScattIndexInternal[Tree[[1]]]/.repChn;
      S2=McKayScattIndexInternal[Tree[[2]]]/.repChn;
      g1=GCD@@S1[[1]]; g2=GCD@@S2[[1]];
      kappa=Abs[McKayDSZ[S1[[1]],S2[[1]]]]/g1/g2;
      Li=Join[S1[[2]],S2[[2]]];
      AppendTo[Li,Subscript[Kr,kappa][Min[g1,g2],Max[g1,g2]]];
      {S1[[1]]+S2[[1]],Li}]];

Options[McKayScattIndexImprovedInternal]={"Debug"->False};

McKayScattIndexImproved[TreeList_,opt:OptionsPattern[]] :=
  Table[
    Simplify[FOmbToOm[Last@McKayScattIndexImprovedInternal[TreeList[[i]],opt][[2]]]],
    {i,Length[TreeList]}];

McKayScattIndexImprovedInternal[Tree_,opt:OptionsPattern[]] :=
  Module[{S1,S2,g1,g2,gFinal,kappa,Li,tem,repOmAttb,rrr},
    If[!ListQ[Tree]||Length[Tree]>2,
      (* leaf: return list of Ombar values for multiplicities 1..GCD(charges) *)
      {Tree,{Join[{1},Table[(y-y^-1)/(j(y^j-y^-j)),{j,2,GCD@@Tree}]]}},
      If[OptionValue["Debug"],Print["Calling with args: ",Tree[[1]],"  |  ",Tree[[2]]]];
      S1=McKayScattIndexImprovedInternal[Tree[[1]],opt]/.repChn;
      S2=McKayScattIndexImprovedInternal[Tree[[2]],opt]/.repChn;
      If[OptionValue["Debug"],Print["S1 is: ",S1,"   S2 is: ",S2]];
      g1=GCD@@S1[[1]]; g2=GCD@@S2[[1]];
      gFinal=GCD@@(S1[[1]]+S2[[1]]);
      kappa=Abs[McKayDSZ[S1[[1]],S2[[1]]]]/g1/g2;
      Li=Join[S1[[2]],S2[[2]]];
      If[OptionValue["Debug"],
        Print["Li is: ",Li,"  g1 is: ",g1,"  g2 is: ",g2,"  gFinal is: ",gFinal]];
      AppendTo[Li,
        repOmAttb=Join[
          Table[CoulombHiggs`OmAttb[{P,0},y_]->Last[S1[[2]]][[P]],{P,1,g1}],
          Table[CoulombHiggs`OmAttb[{0,Q},y_]->Last[S2[[2]]][[Q]],{Q,1,g2}]];
        If[OptionValue["Debug"],Print["repOmAttb is: ",repOmAttb]];
        tem=Table[
          rrr=If[And@@(IntegerQ/@{P g1/gFinal,P g2/gFinal}),
            CoulombHiggs`FlowTreeFormulaRat[{{0,kappa},{-kappa,0}},{g2,-g1},
              {P g1/gFinal,P g2/gFinal},y],
            0];
          Simplify[rrr/.repOmAttb/.
            {CoulombHiggs`OmAttb[{p_,q_},y___]:>0/;p>1||q>1||p q!=0}],
          {P,1,gFinal}];
        If[OptionValue["Debug"],Print["tem is: ",tem]];
        tem];
      {S1[[1]]+S2[[1]],Li}]];

McKayScattGraph[Tree_,m_] :=
  Module[{T,LiArrows,LiVertex},
    T=McKayScattDiagInternal[Tree,m];
    LiArrows=Cases[Flatten[T[[3]]],x_Arrow]/.Arrow[x_]:>x;
    LiVertex=Union[Flatten[LiArrows,1]];
    {LiVertex,
      Table[If[i!=j,Sign[Count[LiArrows,{LiVertex[[i]],LiVertex[[j]]}]],0],
        {i,Length[LiVertex]},{j,Length[LiVertex]}]}];

McKayListAllTrees[Nvec_] :=
  Module[{LiTrees,LiTree1,LiTree2,Li},
    LiTrees={};
    If[Count[Nvec,0]>=Length[Mat]-1,
      (* Nvec is (a multiple of) a simple root: it is itself a leaf *)
      LiTrees={Nvec};,
      Li=Select[BinarySplits[Nvec],McKayDSZ[#,Nvec]!=0&];
      Do[
        LiTree2=McKayListAllTrees[Li[[i]]];
        LiTree1=McKayListAllTrees[Nvec-Li[[i]]];
        Do[AppendTo[LiTrees,{LiTree1[[j]],LiTree2[[k]]}],
          {j,Length[LiTree1]},{k,Length[LiTree2]}];
        ,{i,Length[Li]}];
      ];
    LiTrees];

McKayScanAllTrees[Nvec_,m_] :=
  Module[{Li,Li2},
    Li=McKayListAllTrees[Nvec];
    Li2=Select[Li,Length[McKayScattCheck[#,m][[2]]]>0&];
    DeleteDuplicatesBy[ReverseSortBy[Li2,Length],McKayScattGraph[#,m]&]];

ConstructMcKayDiagram[phimax_,Nm_,m_,ListRays0_] :=
  Module[{ListRays,seen,n,new,kappa,inter,ktab,nNodes},
    nNodes=Length[Mat];
    ListRays=If[ListRays0=={},
      Table[
        {IdentityMatrix[nNodes][[k]],
         InitialRaysOrigin[m][[k]]-5 McKayVec[IdentityMatrix[nNodes][[k]]],
         0,0,0,0,1},
        {k,nNodes}],
      ListRays0];
    seen=AssociationThread[
      Select[ListRays[[All,{3,4}]],First[#]>0&]->True];

    CheckAbort[
      While[True,
        n=Length[ListRays];

        new=Reap[
          Do[
            If[!KeyExistsQ[seen,{i,j}],
              seen[{i,j}]=True;
              kappa=McKayDSZ[ListRays[[i,1]],ListRays[[j,1]]];
              If[kappa=!=0,
                inter=McKayIntersectRaysNoTest[
                  ListRays[[i,1]],ListRays[[j,1]],
                  ListRays[[i,2]],ListRays[[j,2]],m];
                If[inter=!={},
                  ktab=KroneckerDims[Abs[kappa],Nm];
                  Do[
                    If[Plus@@(ktab[[k,1]]ListRays[[i,1]]+ktab[[k,2]]ListRays[[j,1]])<=phimax,
                      Sow[{
                        ktab[[k,1]]ListRays[[i,1]]+ktab[[k,2]]ListRays[[j,1]],
                        inter,i,j,ktab[[k,1]],ktab[[k,2]],
                        1+Max[ListRays[[i,7]],ListRays[[j,7]]]}]],
                    {k,Length[ktab]}]]]],
            {i,n},{j,i+1,n}]
          ][[2]];

        new=If[new==={},{},First[new]];
        If[new==={}||new===Null,Break[],
          Print["Adding ",Length[new]," rays, "];
          ListRays=Join[ListRays,new]]];
      Print[Length[ListRays]," in total."];
      ListRays,

      Print["Aborted -- returning partial ListRays (",Length[ListRays],")."];
      ListRays]];

McKayTreesFromListRays[ListRays_,Nvec_List,m_] :=
  Module[{Lipos,Div,LiTrees},
    Div=Divisors[GCD@@Nvec];
    Lipos=Flatten[Join[Table[Position[ListRays,Nvec/k],{k,Div}]],1];
    If[Lipos=={},
      Print["No such dimension vector in the list"],
      LiTrees=((Plus@@Nvec)/Plus@@ListRays[[#,1]])TreeFromListRays[ListRays,#]&/@
        First[Transpose[Lipos]];
      DeleteDuplicatesBy[
        SortBy[LiTrees,Length[TreeConstituents[#]]&],
        McKayScattGraph[#,m]&]]];



(* ===== Massive Nf=1 SU(2) SW quiver (3 nodes) ===== *)
(* Node basis: (electric, magnetic, flavor); simple roots: {0,1},{2,-1},{-1,0} *)

Mat = {{0,-2,1},{2,0,-1},{-1,1,0}};

McKayDSZ[Nvec_,NNvec_] := Tr[{Nvec}.Mat.Transpose[{NNvec}]];

McKayVec[{n1_,n2_,n3_}] := {Sqrt[3](-n1+n3), -n1+2n2-n3};

McKayRayEq[{n1_,n2_,n3_},{u_,v_}] :=
  n2(1/3+2u) + n3(1/3-u-Sqrt[3]v) + n1(1/3-u+Sqrt[3]v);

McKayRayEq[{n1_,n2_,n3_},{u_,v_},m_] := McKayRayEq[{n1,n2,n3},{u,v}];

McKayIntersectRays[{n1_,n2_,n3_},{nn1_,nn2_,nn3_},m_] :=
  {(n3(2nn1+nn2)+n2(nn1-nn3)-n1(nn2+2nn3)) /
     (6(n3(nn1-nn2)+n1(nn2-nn3)+n2(-nn1+nn3))),
   ((n1+n3)nn2-n2(nn1+nn3)) /
     (2Sqrt[3](n3(-nn1+nn2)+n2(nn1-nn3)+n1(-nn2+nn3)))};

InitialRaysOrigin[m_] := {{1/3,0},{-1/6,-1/(2Sqrt[3])},{-1/6,1/(2Sqrt[3])}};

(* Pass-through rules: dimension vectors are used directly as {n1,n2,n3} lists *)
repChn = {};
repChO = {};


(* ===== Plotting ===== *)

(* Private: convert integer charge vector to a LaTeX string *)
chargeToTeX[charge_] := Module[{terms, tex},
  terms = Select[Thread[{charge, Range@Length@charge}], #[[1]] != 0 &];
  If[terms === {}, Return["0"]];
  tex = Map[Function[t, Module[{c = t[[1]], i = t[[2]]},
    Which[
      c ==  1, "\\gamma_{" <> ToString[i] <> "}",
      c == -1, "-\\gamma_{" <> ToString[i] <> "}",
      c >   0, ToString[c] <> "\\,\\gamma_{" <> ToString[i] <> "}",
      True,    ToString[c] <> "\\,\\gamma_{" <> ToString[i] <> "}"]]],
    terms];
  StringReplace[StringRiffle[tex, " + "], " + -" -> " - "]];

(* Slab-method ray-box clip: returns {t_enter, t_exit} for the segment of
   orig + t*dir that lies inside box, with t >= 0.  Returns {} if ray misses. *)
rayBoxClip[orig_, dir_, {{xlo_,xhi_},{ylo_,yhi_}}] :=
  Module[{txlo,txhi,tylo,tyhi,tEnter,tExit},
    {txlo,txhi} = If[dir[[1]] != 0,
      Sort[{(xlo-orig[[1]])/dir[[1]], (xhi-orig[[1]])/dir[[1]]}],
      If[xlo <= orig[[1]] <= xhi, {-Infinity,Infinity}, {Infinity,-Infinity}]];
    {tylo,tyhi} = If[dir[[2]] != 0,
      Sort[{(ylo-orig[[2]])/dir[[2]], (yhi-orig[[2]])/dir[[2]]}],
      If[ylo <= orig[[2]] <= yhi, {-Infinity,Infinity}, {Infinity,-Infinity}]];
    tEnter = Max[txlo, tylo, 0];
    tExit  = Min[txhi, tyhi];
    If[tEnter < tExit, {tEnter, tExit}, {}]];

Options[McKayPlotDiagram] = {
  "RayLength"  -> 6,
  "PlotRange"  -> {{-1, 1}, {-1, 1}},
  "LabelFrac"  -> 0.72,
  "UseMaTeX"   -> True,
  "FontSize"   -> 11,
  ImageSize    -> 400};

McKayPlotDiagram[ListRays_, m_, opts : OptionsPattern[]] :=
  Module[{L, pRange, nLevels, levelColor, makeLabel, dir, tip, seg, labelPos},
    L      = OptionValue["RayLength"];
    pRange = OptionValue["PlotRange"];
    nLevels = Max[ListRays[[All, 7]]];

    levelColor[lev_] := Blend[
      {RGBColor[0.13, 0.47, 0.71], RGBColor[0.17, 0.63, 0.17],
       RGBColor[1.00, 0.50, 0.05], RGBColor[0.84, 0.15, 0.16]},
      (lev - 2)/Max[nLevels - 2, 1]];

    makeLabel[charge_, pos_] := If[TrueQ[OptionValue["UseMaTeX"]],
      Inset[MaTeX`MaTeX[chargeToTeX[charge], "FontSize" -> OptionValue["FontSize"]], pos],
      Text[Style[charge /. McKayrep, OptionValue["FontSize"]], pos]];

    Graphics[
      Table[
        dir = McKayVec[ListRays[[k, 1]]]/Max[1, GCD @@ ListRays[[k, 1]]];
        tip = ListRays[[k, 2]] + L dir;
        (* visible segment within PlotRange via slab clip *)
        seg = rayBoxClip[ListRays[[k, 2]], dir, pRange];
        labelPos = If[seg === {},
          ListRays[[k, 2]],  (* fallback: shouldn't occur *)
          ListRays[[k, 2]] + ((1 - OptionValue["LabelFrac"]) seg[[1]] +
                                  OptionValue["LabelFrac"]  seg[[2]]) dir];
        {If[ListRays[[k, 3]] == 0,
           {GrayLevel[0.30], AbsoluteThickness[1.2],
            AbsoluteDashing[{5, 3}], Arrowheads[{{0.020, 1}}],
            Arrow[{ListRays[[k, 2]], tip}]},
           {levelColor[ListRays[[k, 7]]], AbsoluteThickness[1.2],
            Arrowheads[{{0.020, 1}}],
            Arrow[{ListRays[[k, 2]], tip}]}],
         makeLabel[ListRays[[k, 1]], labelPos]},
        {k, Length[ListRays]}],
      PlotRange -> pRange,
      Frame -> True,
      FrameStyle -> Directive[GrayLevel[0.55], AbsoluteThickness[0.5]],
      Background -> White,
      ImageSize -> OptionValue[ImageSize]]];


End[]
EndPackage[]
