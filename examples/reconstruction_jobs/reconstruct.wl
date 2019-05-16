(* ::Package:: *)

(* ::Text:: *)
(*In this file we attempt to perform the  functional reconstruction from the files with the evaluations.*)


<<FiniteFlow`


(* ::Text:: *)
(*Because in this file we don't perform any new evaluation of the graph, we don't need to load it, but we can use a dummy graph instead.  We only load the information about the degrees, similarly to what we did in "generate_points.wl".*)


{nparsin,nparsout} = {"NParsIn","NParsOut"}/.FFNParsFromDegreesFile["degrees.fflow"]


FFNewDummyGraph[graph,nparsin,nparsout]


FFLoadDegrees[graph,"degrees.fflow"]


(* ::Text:: *)
(*Now we load all the available evaluations*)


FFLoadEvaluations[graph,FileNames["evaluations_*.fflow"]]


(* ::Text:: *)
(*and finally attempt the reconstruction*)


Print[
  FFReconstructFromCurrentEvaluations[graph,{z1,z2,z3,z4},"MaxPrimes"->5]
]
