(* ::Package:: *)

(* ::Text:: *)
(*Here we choose which file to read the points from, and the subset of points we want to evaluate.*)


(* set the variables here *)
file="points_1.fflow";
start=0;
npoints=113;

(* or from the command line *)
Module[{pos},
  pos = Position[$CommandLine,"file"];
  If[Length[pos]>0, file = $CommandLine[[pos[[1,1]]+1]]; Return[]];
];
Module[{pos},
  pos = Position[$CommandLine,"start"];
  If[Length[pos]>0, start = ToExpression[$CommandLine[[pos[[1,1]]+1]]]; Return[]];
];
Module[{pos},
  pos = Position[$CommandLine,"npoints"];
  If[Length[pos]>0, npoints = ToExpression[$CommandLine[[pos[[1,1]]+1]]]; Return[]];
];


<<FiniteFlow`


(* ::Text:: *)
(*First we define the graph.*)


<<"define_graph.m"


(* ::Text:: *)
(*Now we evaluate the selection of points defined above*)


FFLoadDegrees[graph,"degrees.fflow"]


FFSampleFromPoints[graph,file,start,npoints]


(* ::Text:: *)
(*and save the results*)


outfile = "evaluations_"<>ToString[start]<>"_"<>ToString[npoints]<>"_"<>file


FFDumpEvaluations[graph,outfile]
