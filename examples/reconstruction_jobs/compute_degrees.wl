(* ::Package:: *)

<<FiniteFlow`


(* ::Text:: *)
(*First we define the graph.*)


<<"define_graph.m"


(* ::Text:: *)
(*Then we compute all the total and partial degrees. Note that only the total degrees are returned but the partial ones are also computed and stored internally.*)


FFAllDegrees[graph]


(* ::Text:: *)
(*Finally we save all the information about the computed degrees on file.*)


FFDumpDegrees[graph,"degrees.fflow"]
