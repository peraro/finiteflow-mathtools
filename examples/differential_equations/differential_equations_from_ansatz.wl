(* ::Package:: *)

(* ::Text:: *)
(*In this file we compute differential equations from an ansatz of the solution, which we infer from the alphabet of this topology, assuming the latter is known a priori.*)
(**)
(*We recommend looking at the example in differential_equations.wl before this one, since it is simpler and it shares part of the algorithm. The first steps are exactly the same as in the other example, so we refer to  the comments in differential_equations.wl for a description.*)


<<LiteIBP`


<<"topology.wl"


utints = Import["ibps/utints_def.m"];


integrals = Join[LIBPGetAllInts[FileNames["ibps/ints_family_*.mx"]], First/@utints];


FFNewGraph[graph,in,{eps,s12,s13}]


eqsfiles=FileNames["ibps/sids_*.json"];


neededintegrals = Select[integrals,LIBPIntS[#]<=3 && ((LIBPIntR[#]-LIBPIntT[#])<=1)&];


LIBPWriteSystemJSON[eqsfiles,integrals,neededintegrals,{eps,s12,s13},"FileName"->"system.json"]


FFAlgJSONSparseSolver[graph,ibps,{in},"system.json"]

FFSolverOnlyHomogeneous[graph,ibps]


FFGraphOutput[graph,ibps];
ibplearn = FFSparseSolverLearn[graph,integrals];
{nonmasters, masters} = {"DepVars", "IndepVars"}/.ibplearn;


mastersexpr = masters /. utints;


LIBPComputeDerivatives[family]


computederivative[x_]:=LIBPEliminateZeroSectors[Collect[LIBPDeriv[#,x]&/@mastersexpr,_j,Together]];
derivatives = computederivative[#]&/@{s12,s13,s23};


newneededints = LIBPIntegralsIn[derivatives];
SubsetQ[neededintegrals,newneededints]


FFSolverResetNeededVars[graph,ibps,integrals,newneededints]


newibplearn = FFSparseSolverLearn[graph,integrals];
masters == ("IndepVars"/.newibplearn)


{nonmasters, masters} = {"DepVars", "IndepVars"} /. newibplearn;


FFSparseSolverMarkAndSweepEqs[graph,ibps]

FFSparseSolverDeleteUnneededEqs[graph,ibps]


FFAlgRatFunEval[graph, unreduced, {in},
                {eps,s12,s13},
                Join@@(FFLinearCoefficients[#,nonmasters]&/@(Join@@derivatives))/.s23->1]


FFAlgMatMul[graph,de,{unreduced,ibps},
  Length[{s12,s13,s23}]*Length[masters], Length[nonmasters], Length[masters]
]


FFGraphOutput[graph,de]


FFAlgNonZeroes[graph,denz,{de}]


FFGraphOutput[graph,denz];
nzlearn = FFNonZeroesLearn[graph];


FFGraphDraw[graph]


(* ::Text:: *)
(*Until this point, everything was the same as in the file differential_equations.wl.  In order to make an appropriate ansatz for each element in the output, we must track the original matrix it is coming from.  This is easily done by assigning a symbolic label, say de[x,i,j] to the matrix element {i,j} of the differential equations w.r.t. x and see how it propagates in the node denz.*)


symbolicelements = Join@@Table[Join@@Table[de[x,ii,jj],{ii,Length[masters]},{jj,Length[masters]}],{x,{s12,s13,s23}}];
nonzeroentries = symbolicelements[["NonZero"/.nzlearn]];


nonzeroentries//Length


(* ::Text:: *)
(*We make the following ansatz for the solution:*)
(**)
(*  de[x,i,j]  = Sum[ eps * D[Log[w],x] coeff[w,x,i,j] , {w,alphabet}]*)
(**)
(*where the unknowns coeff[w,x,i,j] are Q-numbers.*)


(* ::Text:: *)
(*This is the alphabet we will use, and its derivatives (we only select the non-vanishing ones):*)


alphabet={w1->s12,w2->s23,w3->s13,w4->s23+s12,w5->s23+s13,w6->s12+s13,w7->s23+s12+s13};

Do[
  dlogwdx[x]=Select[Table[dlogwdx[w[[1]],x]->Together[eps D[Log[w[[2]]],x]],{w,alphabet}],!TrueQ[#[[2]]==0]&];
 ,{x,{s12,s13,s23}}];
 
alldlogdx = Join@@(dlogwdx/@{s12,s13,s23})


(* ::Text:: *)
(*We put the derivatives, which form our ansatz, inside a node*)


FFAlgRatFunEval[graph,ansatz,{in},{eps,s12,s13},(#[[2]]&/@alldlogdx)/.s23->1]


(* ::Text:: *)
(*and chain it to the previous output node*)


FFAlgChain[graph,final,{ansatz,denz}]

FFGraphOutput[graph,final]


finalelements = Join[First/@alldlogdx,nonzeroentries];


FFGraphDraw[graph]


(* ::Text:: *)
(*Now we create a second graph, with no input node, and use a subgraph multi-fit algorithm to perform the fit for the non-vanishing matrix elements, based on our ansatz.*)


FFNewGraph[graph2]


(* ::Text:: *)
(*Recall that the non-zero entries of the matrix elements have the symbolic form de[x,i,j].  We use that to define a linear fit for each of them:*)


fitentries = Table[
  Join[First/@dlogwdx[el[[1]]],{el}]
,{el,nonzeroentries}];


(* ::Text:: *)
(*and the  unknown coefficients coeff[w,x,i,j]*)


coefficients = Table[
  Table[coeff[w[[1]],el[[1]],el[[2]],el[[3]]],{w,First/@dlogwdx[el[[1]]]}]
,{el,nonzeroentries}];


FFAlgSubgraphMultiFit[graph2,multifit,{},graph,{eps,s12,s13},finalelements->fitentries]


(* ::Text:: *)
(*This will do several fits at once.  As an example, for the matrix element de[s12,1,1], it will perform the following fit*)


fitentries[[1]].Join[coefficients[[1]],{0}]==fitentries[[1,-1]]


(* ::Text:: *)
(*where we recall that*)


dlogwdx[s12]


(* ::Text:: *)
(*and similar for the other matrix elements.*)


(* ::Text:: *)
(*We now run the learning step*)


FFGraphOutput[graph2,multifit];
fitlearn = FFMultiFitLearn[graph2,coefficients];


(* verify the output of learn is a list of list, which implies no fit has failed *)
And@@(#[[0]]==List&/@fitlearn)


(* ::Text:: *)
(*We are now ready to reconstruct the numerical result*)


rec = FFReconstructNumeric[graph2];
sol = Join@@FFMultiFitSol[rec,fitlearn];


(* ::Text:: *)
(*insert it back in the ansatz*)


recnonzeroes =Table[fitentries[[i]].Join[coefficients[[i]],{0}],{i,Length[fitentries]}]/.Dispatch[sol]/.alldlogdx;


(* ::Text:: *)
(*and build the differential equation matrices from it*)


(* de[x] will contain the matrix of differential equations with respect to x *)
invariants = {s12,s13,s23};
deall = ArrayReshape[Normal[FFNonZeroesSol[recnonzeroes,nzlearn]],{Length[invariants],Length[masters],Length[masters]}];
Do[
  de[invariants[[i]]] = deall[[i]];
,{i,Length[{s12,s13,s23}]}];


(* Check some expected properties of the result *)
de[s12].de[s13]-de[s13].de[s12]//Together

D[de[s12],s13]-D[de[s13],s12]//Together



