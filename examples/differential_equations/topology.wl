(* ::Package:: *)

<<LiteIBP`


(* Set space-time dimensions in LiteRed *)
SetDim[4-2 eps]


(* Kinematics *)
Internal = {k1, k2};
External = {p1,p2,p3,p4};
MomentumConservation = {mm[p4]->-mm[p1]-mm[p2]-mm[p3]};
Replacements = {
 mmp2[p1] -> 0, 
 mmp2[p2] -> 0, 
 mmp2[p3] -> 0, 
 mmp[p1,p2] -> s12/2, 
 mmp[p1,p3] -> s13/2, 
 mmp[p2,p3] -> s23/2};
Invariants = {
  s12 -> 2 mmp[p1,p2],
  s13 -> 2 mmp[p1,p3],
  s23 -> 2 mmp[p2,p3]
};


(* Loop propagators *)
Propagators = {
-mp2[mm[k1]], (* = -mmp2[k1] *)
-mp2[mm[k1] + mm[p1]],
-mp2[mm[k1] - mm[p2] - mm[p3]],
-mp2[mm[k2]], (* = -mmp2[k2] *)
-mp2[mm[k2] - mm[p2] - mm[p3]],
-mp2[mm[k2] - mm[p3]],
-mp2[mm[k2] - mm[k1]],
-mp2[mm[k1] + mm[p2]],
-mp2[mm[k2] + mm[p1]]};


(* Define the family in LiteIBP *)
LIBPFamilyLite[family,
  Propagators /. MomentumConservation /. Replacements,
  Internal,
  External,
  MomentumConservation,
  Replacements,
  Invariants]


(* Define the family in LiteRed and find symmetries *)
NewBasis[family,LIBPToLiteRed[Propagators],LIBPLoopMomenta[family],GenerateIBP->True]
AnalyzeSectors[family,{_,_,_,_,_,_,_,0,0}];
LIBPFindSymmetries[family,EMs->LIBPEMs[family]];



