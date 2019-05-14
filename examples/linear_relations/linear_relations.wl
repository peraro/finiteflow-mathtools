(* ::Package:: *)

(* ::Text:: *)
(*In this example we find linear relations between the elements of the output of a graph.   In particular, we try writing the most complicated elements in terms of simpler ones, before performing any analytic reconstruction.  We then proceed to reconstruct only an independent subset of the functions in the output.*)


<<FiniteFlow`


(* ::Text:: *)
(*First we must define the graph.*)
(**)
(*For simplicity, we define it from a list of rational functions, but it should be clear that the method presented here can be applied to any other arbitrarily complicated graph.  Indeed, it is significantly more useful where an analytic representation of the full output is not known.*)


FFNewGraph[graph,in,{z1,z2,z3,z4}]


FFAlgRatFunEval[graph,out,{in},{z1,z2,z3,z4},Import["full_output.m"]]


FFGraphOutput[graph,out]
FFGraphDraw[graph]


(* ::Text:: *)
(*The output is a list of 15 elements, as we can verify with*)


FFNParsOut[graph,out]


(* ::Text:: *)
(*Let's call these functions {f[1],f[2],...}.  We want to find linear relations between them.  This can be done by solving the linear fit problem*)
(**)
(*  c[1] f[1] + c[2] f[2] + ... == 0*)
(**)
(*for the unknown coefficients c[i].  Here we want the unknowns to be Q-numbers (it is also possible to find them such that they depend on a subset of the input variables).*)


(* ::Text:: *)
(*For later convenience, we define the lists of symbolic coefficients and functions*)


coefficients = c/@Range[FFNParsOut[graph,out]]

functions = f/@Range[FFNParsOut[graph,out]]


(* ::Text:: *)
(*Now we want to sort the output elements of the graph by their "complexity", so that the most complex elements are written in terms of the simplest ones. *)


(* ::Text:: *)
(*There are several ways of defining this complexity.  A simple and effective way is defining it based on the total degree of the final result.*)


degrees = FFTotalDegrees[graph];
complexity = Max@@#&/@degrees


(* ::Text:: *)
(*A more refined one is based of the number of evaluation points needed for the reconstruction of each function*)


complexity = FFNSamplePoints[graph][[2]]


(* ::Text:: *)
(*Hence we define our sorted functions and coefficients*)


sortedfunctions = SortBy[functions,complexity[[#[[1]]]]&]

sortedcoefficients = c@@#&/@sortedfunctions


(* ::Text:: *)
(*and we append to the graph a node which sorts the functions as above*)


FFAlgTake[graph,sorted,{out},{functions}->sortedfunctions]


(* ::Text:: *)
(*In order to do the fit, we create a node evaluating the r.h.s., i.e. zero, and we chain it with the current output.*)


FFAlgRatNumEval[graph,zero,{0}]

FFAlgChain[graph,fiteq,{sorted,zero}]

FFGraphOutput[graph,fiteq]


(* ::Text:: *)
(*Now we have modified the graph as follows*)


FFGraphDraw[graph]


(* ::Text:: *)
(*We now create a new graph with a SubgraphFit algorithm solving the linear fit equation above.*)


FFNewGraph[graphfit]

FFAlgSubgraphFit[graphfit,fit,{},graph,{z1,z2,z3,z4},sortedcoefficients]

FFGraphOutput[graphfit,fit]


(* ::Text:: *)
(*Like for all linear fits (and dense solvers in general) we need to run the learning phase*)


fitlearn=FFDenseSolverLearn[graphfit,sortedcoefficients]


(* ::Text:: *)
(*and then we can reconstruct the solution*)


fitrec = FFReconstructNumeric[graphfit];
fitsol = FFDenseSolverSol[fitrec,fitlearn];


(* ::Text:: *)
(*We can easily convert this into a set of linear relations between the functions f[i].  We have a dedicate utility function in the FFUtils package*)


<<FFUtils`


linrels=FFLinearRelationsFromFit[sortedfunctions,sortedcoefficients,fitsol]


(* ::Text:: *)
(*Now the graph "graphfit" is no longer needed.  Because it is using "graph" as a subgraph (which would prevent us from modifying "graph") we must delete it.*)


FFDeleteGraph[graphfit]


(* ::Text:: *)
(*Now we define our list of independent functions, which we want to reconstruct, based on the results of the fit.*)


independentfuncs = Complement[functions,First/@linrels]


(* ::Text:: *)
(*We can symbolically write the output of graph as*)
(**)
(*  f[1] e[1] + f[2]e[2] + ...*)
(**)
(*where e[i] is the unit vector in the i-th direction.  Hence, taking into account these relations, we write the output as:*)


symbolicout = Collect[Sum[e[i] f[i],{i,Length[functions]}]/.linrels,_f,Together]


(* ::Text:: *)
(*Now we finally reconstruct the independent functions appearing in the result.*)


(* ::Text:: *)
(*First we bring back "graph" as it was originally*)


FFGraphOutput[graph,out]

FFGraphPrune[graph]


FFGraphDraw[graph]


(* ::Text:: *)
(*Then we take only the independent elements of the output*)


FFAlgTake[graph,indep,{out},{functions}->independentfuncs]


(* ::Text:: *)
(*and reconstruct them*)


FFGraphOutput[graph,indep];
rec = FFReconstructFunction[graph,{z1,z2,z3,z4}];


(* ::Text:: *)
(*Here is their analytic expression, as a list of substitution rules*)


indepfrules=Inner[Rule,independentfuncs,rec,List]


(* ::Text:: *)
(*and the output of the graph can be represented as*)


symbolicout /. indepfrules
