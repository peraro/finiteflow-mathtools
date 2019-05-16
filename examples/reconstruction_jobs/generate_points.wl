(* ::Package:: *)

(* ::Text:: *)
(*In this file we generate a list of sample points.  We need to decide how many primes are needed in total.*)


(* set the number of primes here *)
nprimes=1;

(* or from the command line *)
Module[{pos},
  pos = Position[$CommandLine,"nprimes"];
  If[Length[pos]>0, nprimes = ToExpression[$CommandLine[[pos[[1,1]]+1]]]; Return[]];
];


(* ::Text:: *)
(*We use the information about the degrees from "degrees.fflow".  We can also read the number of inputs and outputs from it.*)


<<FiniteFlow`


{nparsin,nparsout} = {"NParsIn","NParsOut"}/.FFNParsFromDegreesFile["degrees.fflow"]


(* ::Text:: *)
(*There's no need to load the actual graph, since we don't need to evaluate it in this file.  A dummy graph does the job.*)


FFNewDummyGraph[graph,nparsin,nparsout]


(* ::Text:: *)
(*We load the information about the total degrees.*)


FFLoadDegrees[graph,"degrees.fflow"]


(* ::Text:: *)
(*We also load the list of evaluations which have already completed (if any), so that we don't generate points again for them.*)


FFLoadEvaluations[graph,FileNames["evaluations_*.fflow"]]


(* ::Text:: *)
(*Hence, we generate and save a list of points.  These are only points for evaluations which have not been completed yet.*)


savefile = "points_"<>ToString[nprimes]<>".fflow"


FFDumpSamplePoints[graph,savefile,"MaxPrimes"->nprimes]


(* ::Text:: *)
(*We can read how many points we generated form the file.*)


FFSamplesFileSize[savefile]
