(* ::Package:: *)

(* ::Text:: *)
(*In this file, we compute the 4-photon helicity amplitude with the -+-+ helicity configuration.*)
(**)
(*The input is a representation of the unreduced amplitude in terms of Feynman integrals in a standard representation (i.e. a representation typically used by IBP programs).  We directly reconstruct the coefficients of the polylogarithmic functions appearing in the final result.*)
(**)
(*Note that, by writing the kinematic dependence in terms of s and t only, we are neglecting an helicity-dependent overall prefactor, which is however easy to recover (see e.g. arXiv:1605.02172 or arXiv:1608.01902).*)
(**)
(*Although this is just a one-loop amplitude, the same setup can be used for multi-loop amplitudes as well.*)
(**)
(*A more complete (and more complicated) version of this calculation can be found in the file compute_amplitude.wl, where the input is an integrand written in terms of loop momenta and polarization vectors.  We recommend looking at the example in this file first, since it is significantly simpler.*)


<<FiniteFlow`


<<LiteMomentum`


<<"topology.wl"


(* ::Text:: *)
(*We load a representation of the unreduced amplitude (where we already set s=1).  In this example we assume such a representation to be known analytically.  One can also compute it numerically using a more complicated graph than the one presented here (an example is in compute_amplitude.wl).*)


invariants = {t,mf2}


unreducedamp=Import["unreduced_amplitude_-+-+.m"];


(* ::Text:: *)
(*We define our graph for the full amplitude, which depends on the invariants and epsilon.  Because eventually we want to expand in epsilon, the latter needs to be the first variable.*)


FFNewGraph[gredfull];
FFGraphInputVars[gredfull,in,Join[{eps},invariants]]


(* ::Text:: *)
(*The integrals we want to reduced are the one appearing in our expression for the unreduced amplitude*)


neededintegrals = LIBPIntegralsIn[unreducedamp];


(* ::Text:: *)
(*We load the IBP system of equations.  More comments of IBPs can be found in the example "examples/differential_equations/differential_equations.wl".*)


allintegrals = LIBPGetAllInts[FileNames["ibps/ints_*.mx"]];

(* create system.json, which collects info on the whole IBP system *)
LIBPWriteSystemJSON[
  FileNames["ibps/sids_*.json"],
  allintegrals,neededintegrals,
  {eps,t,mf2},
  "FileName"->"system.json"];
                         
FFAlgJSONSparseSolver[gredfull,ibps,{in},"system.json"];
FFSolverOnlyHomogeneous[gredfull,ibps];
FFGraphOutput[gredfull,ibps];


(* ::Text:: *)
(*Learning phase for IBP system.  This yields the list of master integrals:*)


ibplearn = FFSparseSolverLearn[gredfull,allintegrals];
nonmis = "DepVars"  /. ibplearn;
mis = "IndepVars" /. ibplearn;
Print["Master integrals = ",mis];


(* ::Text:: *)
(*We  now run the mark-and-sweep algorithm (highly recommended).*)


FFSparseSolverMarkAndSweepEqs[gredfull,ibps]
FFSparseSolverDeleteUnneededEqs[gredfull,ibps];


(* ::Text:: *)
(*Full list of relevant IBP integrals, sorted with the masters on the right.*)


ibpintegrals = Join[nonmis,mis];


(* ::Text:: *)
(*We want to implement the IBP reduction as a matrix multiplication.  Because in this example the master integrals also appear in the unreduced amplitude, we must complete the reduction returned by the IBP node with the trivial reduction of the masters to themselves, the latter defined by an Identity matrix.*)


FFAlgRatNumEval[gredfull,misred,Join@@IdentityMatrix[Length[mis]]];
FFAlgChain[gredfull,ibpsfull,{ibps,misred}]


(* ::Text:: *)
(*We then get the coefficients of the unreduced amplitude*)


unreducedcoeffs=FFLinearCoefficients[unreducedamp,ibpintegrals];


(* ::Text:: *)
(*and put them in an algorithm.*)


FFAlgRatFunEval[gredfull,unreduced,{in},
                Join[{eps},invariants],
                unreducedcoeffs]


(* ::Text:: *)
(*We finally get the reduced amplitude via a matrix multiplication between the coefficients of the unreduced amplitude and the ones of the IBP reduction.*)


FFAlgMatMul[gredfull,reduced,{unreduced,ibpsfull},1,Length[ibpintegrals],Length[mis]]


FFGraphOutput[gredfull,reduced]


(* ::Text:: *)
(*The graph so far*)


FFGraphDraw[gredfull]


(* ::Text:: *)
(*The following reconstructs the amplitude as a linear combination of master integrals.  This is however not our goal here, since we want to expand the amplitude into polylogarithms.*)


(*FFGraphOutput[gredfull,reduced];
reducedamp = FFReconstructFunction[gredfull,Join[{eps},invariants]].mis;




(* ::Text:: *)
(*Now we substitute analytic expressions for the master integrals, which are valid up to O (eps^0)*)
(**)
(*Analytic expressions for the master integrals, taken from https://arxiv.org/pdf/0804.0749.pdf*)
(*In particular:*)
(* - boxfun[X,Y,Z] is H[X,Y]+H[X,Z], with H defined in Eq. (C.28)*)
(* - trianglefun is Log[(1-rho)/(1+rho)], with rho = Sqrt[1-4 mf2/s]*)
(*  - bubblefun is f defined in Eq. (C.15-17)*)


(* the usual prefactor *)
rgamma[eps_] := Gamma[1+eps]Gamma[1-eps]^2/Gamma[1 - 2 eps];


boxint[s_,t_,mf2_]:=-1/(s t) ( boxfun[(s+t) mf2/(s t), mf2/s, mf2/t] );
trint[s_,mf2_]:=-1/(2 s) trianglefun[mf2/s];
bubble[s_,mf2_]:=Normal[Series[1/rgamma[eps](mf2^(- eps) Gamma[1+eps]/eps+2 + bubblefun[mf2/s]),{eps,0,0}]];
tadpole[mf2_]:=Normal[Series[-1/rgamma[eps] Gamma[1-(4-2 eps)/2]/Gamma[1] 1/(mf2)^(1-(4-2 eps)/2),{eps,0,0}]];


box1mis = {
 j[box1,1,1,1,1] -> boxint[s,t,mf2],
 j[box1,0,1,1,1] -> trint[t,mf2],
 j[box1,1,0,1,1] -> trint[s,mf2],
 j[box1,0,1,0,1] -> bubble[t,mf2],
 j[box1,1,0,1,0] -> bubble[s,mf2],
 j[box1,0,0,0,1] -> tadpole[mf2]
};


box1mis


(* All analytic master integrals, including permutations *)
analyticmis = Join[
     box1mis,
     box1mis /. box1->box2 /. InvariantsPermutations[[2]],
     box1mis /. box1->box3 /. InvariantsPermutations[[3]]];
analyticmis = mis /. analyticmis;


(* get list of special functions appearing in the master integrals *)
funcspattern = _Log|_boxfun|_trianglefun|_bubblefun;
allfuncs = Union[Cases[analyticmis,funcspattern,Infinity]]


(* ::Text:: *)
(**)
(*We create a list of rules for expressing master integrals in terms of *)
(*special functions*)


misrules  = Together[(Times@@(allfuncs^#[[1]]))->#[[2]]&/@CoefficientRules[#,allfuncs]&/@analyticmis];


(* ::Text:: *)
(*and obtain a function basis from the rules *)


functionbasis = Union@@(First/@#&/@misrules)


(* ::Text:: *)
(*We build a matrix converting the master integrals into the function basis *)


mis2basis = (functionbasis/.Join[#,(#->0)&/@Complement[functionbasis,First/@#]])&/@misrules


(* ::Text:: *)
(*put it in a graph*)


FFAlgRatFunEval[gredfull,mis2fs,{in},Join[{eps},invariants],Join@@mis2basis /. s -> 1]


(* ::Text:: *)
(*and multiply the reduced amplitude with it*)


FFAlgMatMul[gredfull,amp2fs,{reduced,mis2fs},1,Length[mis],Length[functionbasis]]


FFGraphOutput[gredfull,amp2fs]


(* ::Text:: *)
(*The final graph, which computes the coefficients of the polylog functions.  Now we want to expand these coefficients in epsilon*)


FFGraphDraw[gredfull]


(* ::Text:: *)
(*This graph defines the Laurent expansion*)


FFNewGraph[gexpansion,in,invariants];
FFAlgLaurent[gexpansion,laurent,{in},gredfull,0]


(* ::Text:: *)
(*The Laurent expansion has a learning step*)


FFGraphOutput[gexpansion,laurent];
explearn = FFLearn[gexpansion];


(* ::Text:: *)
(*We finally reconstruct analytically the coefficients of the expansion*)


rec = FFReconstructFunction[gexpansion,invariants];


(* ::Text:: *)
(*From this we build the final amplitude, up to the finite part.  Because the tree-level is vanishing, the result is finite.*)


amplitude=(Normal/@FFLaurentSol[Factor[rec],eps,explearn]).functionbasis


(* ::Text:: *)
(*In a more complex example, before reconstructing the final expression, we recommend finding linear relations between the output elements of the final graph.  This allows to reconstruct only a subset of linearly independent coefficients of the Laurent expansion of the output, which can simplify both the reconstruction and the result itself.  An explicit example of this can be found in the "examples/linear_relations" directory of this repository.*)
