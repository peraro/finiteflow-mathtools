(* ::Package:: *)

<<LiteIBP`


(* Set space-time dimensions in LiteRed *)
SetDim[4-2 eps]


(* Kinematics *)
Internal = {k};
External = {p1,p2,p3,p4};
MomentumConservation = {mm[p4]->-mm[p1]-mm[p2]-mm[p3]};
Replacements = {
 mmp2[p1] -> 0,
 mmp2[p2] -> 0,
 mmp2[p3] -> 0,
 mmp[p1, p2] -> s/2, 
 mmp[p1, p3] -> (-s - t)/2, 
 mmp[p2, p3] -> t/2
};
Invariants = {
  s -> 2 mmp[p1,p2],
  t -> 2 mmp[p2,p3]
};


LoopMomenta = {
  mm[k],
  mm[k]-mm[p1],
  mm[k]-mm[p1]-mm[p2],
  mm[k]+mm[p4]
};


(* ::Text:: *)
(*We have three diagrams, which are the following permutations of the same box*)


permutations = {
  {1,2,3,4},
  {2,1,3,4},
  {4,1,3,2}
};


MomentaPermutations = Table[Inner[Rule,{p1,p2,p3,p4},{p1,p2,p3,p4}[[permutations[[ii]]]],List],{ii,Length[permutations]}]


InvariantsPermutations = Table[Inner[Rule,{s,t},{s,t}/.Invariants/.MomentaPermutations[[ii]]/.MomentumConservation/.Replacements//Expand,List],{ii,Length[permutations]}]


(* ::Text:: *)
(*For each permutation we define a new integral family in LiteIBP and LiteRed.  We also find mappings between different families.*)


(* Declare the internal fermion mass,
   and make LiteRed happy *)
Declare[mf2,Number];


(* permutation 1 *)


Propagators1 = (mp2[#]-mf2)&/@LoopMomenta/. MomentumConservation/.Replacements//Expand


LIBPFamilyLite[box1, 
  Propagators1,
  Internal,
  External,
  MomentumConservation,
  Replacements,
  Invariants]


NewBasis[box1,LIBPToLiteRed[Propagators1],LIBPLoopMomenta[box1],GenerateIBP->True]
AnalyzeSectors[box1,{_,_,_,_}];
LIBPFindSymmetries[box1,EMs->LIBPEMs[box1]];


(* permutation 2 *)


Propagators2 = (mp2[#]-mf2/.MomentaPermutations[[2]])&/@LoopMomenta/. MomentumConservation/.Replacements//Expand


LIBPFamilyLite[box2, 
  Propagators2,
  Internal,
  External,
  MomentumConservation,
  Replacements,
  Invariants]


NewBasis[box2,LIBPToLiteRed[Propagators2],LIBPLoopMomenta[box2],GenerateIBP->True];
AnalyzeSectors[box2,{_,_,_,_}];
FindExtSymmetries[box2,box1,EMs->LIBPEMs[box2]];
LIBPFindSymmetries[box2,EMs->LIBPEMs[box2]];


(* permutation 3 *)


Propagators3 = (mp2[#]-mf2/.MomentaPermutations[[3]])&/@LoopMomenta/. MomentumConservation/.Replacements//Expand


LIBPFamilyLite[box3, 
  Propagators3,
  Internal,
  External,
  MomentumConservation,
  Replacements,
  Invariants]


NewBasis[box3,LIBPToLiteRed[Propagators3],LIBPLoopMomenta[box3],GenerateIBP->True]
AnalyzeSectors[box3,{_,_,_,_}];
FindExtSymmetries[box3,box2,box1,EMs->LIBPEMs[box3]]
LIBPFindSymmetries[box3,EMs->LIBPEMs[box3]]



