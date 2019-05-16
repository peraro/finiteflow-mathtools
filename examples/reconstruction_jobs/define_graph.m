(* ::Package:: *)

(* ::Text:: *)
(*This defines a graph.*)
(**)
(*We use, for simplicity, a list of analytic functions, but this is irrelevant for the purpose of the example in this directory.*)


<<FiniteFlow`


FFNewGraph[graph,in,{z1,z2,z3,z4}]


FFAlgRatFunEval[graph,out,{in},{z1,z2,z3,z4},
{(z1^3*z2^2)/(6123459*(1 + z2)^5), ((z1^3*z2^3)/6 - (z1^2*z2^2*z3)/3 - 
   (z1^2*z2^3*z3)/3 + (z1*z2^2*z3^2)/6 - (z1^2*z2^2*z3^2)/6 + 
   (z1*z2^3*z3^2)/6 - (z1^3*z2^3*z4)/3 + (z1^2*z2^2*z3*z4)/3 + 
   (z1^2*z2^3*z3*z4)/3 - (z1^2*z2^2*z4^2)/6 - (z1^2*z2^3*z4^2)/6)/
  (1346 + z2 + z1*z2)^5, (z1*z2*z3^7)/6, (z1^2*z2^2*z4^2)/6, 
 ((z1^4*z2^2)/6 - (z1^8*z2^2*z3)/3 + (z1^2*z2^2*z3^2)/6 + (z1^3*z2^2*z4)/3 - 
   (z1^2*z2^2*z3*z4)/3 + (z1^2*z2^2*z4^2)/6)/(116 + z1)^3}]


FFGraphOutput[graph,out]


(*FFReconstructFunction[graph,{z1,z2,z3,z4},"MaxPrimes"\[Rule]2,"PrintDebugInfo"\[Rule]1]*)
