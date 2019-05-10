(* ::Package:: *)

<<LiteIBP`


<<"topology.wl"


utints = Import["ibps/utints_def.m"];


integrals = Join[LIBPGetAllInts[FileNames["ibps/ints_family_*.mx"]], First/@utints];


(* ::Text:: *)
(*We create a dataflow graph for*)
(*1) finding a list of master integrals*)
(*2) writing differential equations for them*)


FFNewGraph[graph,in,{eps,s12,s13}]


eqsfiles=FileNames["ibps/sids_*.json"];


(* ::Text:: *)
(*In general, one cannot reduce all the integrals in a system to master.  Therefore, we need to specify the subset we are interested in.  Unfortunately, we still don't know what it the correct subset for computing differential equations, because we still have to determine the list of master integrals.  Therefore we specify one which we can reasonably guess to be a superset of them, i.e. integrals up to rank 3 and 1 dotted propagators.*)


neededintegrals = Select[integrals,LIBPIntS[#]<=3 && ((LIBPIntR[#]-LIBPIntT[#])<=1)&];


(* write system.json with all the info on the system *)
LIBPWriteSystemJSON[eqsfiles,integrals,neededintegrals,{eps,s12,s13},"FileName"->"system.json"]


(* create a linear solver form system.json *)
FFAlgJSONSparseSolver[graph,ibps,{in},"system.json"]

(* don't include constant term in the (numerical) solution *)
FFSolverOnlyHomogeneous[graph,ibps]


(* Learning phase: get the master integrals *)
FFGraphOutput[graph,ibps];
ibplearn = FFSparseSolverLearn[graph,integrals];
{nonmasters, masters} = {"DepVars", "IndepVars"}/.ibplearn;


(* ::Text:: *)
(*Now we have the list of master integrals*)


masters


(* ::Text:: *)
(*which is a subset of the UT integrals.  We get their explicit expression*)


mastersexpr = masters /. utints;


(* ::Text:: *)
(*We now want to compute their derivatives.  We use an utility for computing derivative operators.*)


LIBPComputeDerivatives[family]


(* ::Text:: *)
(*and now the derivative of a linear combination of integrals in LiteRed notation can be computed using the LIBPDeriv function.  Remember that vanishing sectors are not in the IBP system, and therefore we need to get rid of them manually, using LIBPEliminateZeroSectors.*)


(* Derivative of all the masters with respect to x *)
computederivative[x_]:=LIBPEliminateZeroSectors[Collect[LIBPDeriv[#,x]&/@mastersexpr,_j,Together]];
derivatives = computederivative[#]&/@{s12,s13,s23};


(* The new list of needed integrals must be a subset of the previous one *)
newneededints = LIBPIntegralsIn[derivatives];
SubsetQ[neededintegrals,newneededints]


(* ::Text:: *)
(*We reset the list of needed integrals in the solver*)


FFSolverResetNeededVars[graph,ibps,integrals,newneededints]


(* ::Text:: *)
(*and repeat the learning phase,  verifying the list of masters has not changed.*)


newibplearn = FFSparseSolverLearn[graph,integrals];
masters == ("IndepVars"/.newibplearn)


{nonmasters, masters} = {"DepVars", "IndepVars"} /. newibplearn;


(* ::Text:: *)
(*At this stage we strongly recommend  to also run the mark-and-sweep routine to eliminate unneeded equations from the system*)


FFSparseSolverMarkAndSweepEqs[graph,ibps]

FFSparseSolverDeleteUnneededEqs[graph,ibps] (* optional, but frees some memory *)


(* ::Text:: *)
(*We encode the coefficients of the unreduced derivatives of the masters w.r.t. the non-master integrals*)


FFAlgRatFunEval[graph, unreduced, {in},
                {eps,s12,s13},
                Join@@(FFLinearCoefficients[#,nonmasters]&/@(Join@@derivatives))/.s23->1]


(* ::Text:: *)
(*and apply the IBP reduction, via a matrix multiplication algorithm, to finally get the differential equations*)


FFAlgMatMul[graph,de,{unreduced,ibps},
  Length[{s12,s13,s23}]*Length[masters], Length[nonmasters], Length[masters]
]


FFGraphOutput[graph,de]


(* ::Text:: *)
(*Here is our graph*)


FFGraphDraw[graph]


(* ::Text:: *)
(*Since differential equations are sparse, we append a NonZero algorithm which only takes the non-zero elements.  This improves memory usage in the reconstruction (especially in more complicated applications).*)


FFAlgNonZeroes[graph,denz,{de}]


FFGraphOutput[graph,denz];
nzlearn = FFNonZeroesLearn[graph];


(* ::Text:: *)
(*Since we are using a UT basis, we expect the system to be in epsilon-form.  We can check this property before the analytic reconstruction, by multiplying the output by 1/eps and verifying the result does not depend on eps.*)


FFAlgRatFunEval[graph,1/eps,{in},{eps,s12,s13},{1/eps}]

FFAlgMatMul[graph,final,{1/eps,denz},1,1,FFNParsOut[graph,denz]]

FFGraphOutput[graph,final];

FFIndependentOf[graph,{eps,s12,s13},eps]


(* ::Text:: *)
(*This is the final graph which returns the (non-vanishing) matrix elements of the differential equations matrices, divided by eps.*)


FFGraphDraw[graph]


(* ::Text:: *)
(*Now we can reconstruct the matrix elements of the differential equations*)


rec0=FFReconstructFunction[graph,{eps,s12,s13}];


(* ::Text:: *)
(*we add back the removed zeroes*)


rec = Normal[FFNonZeroesSol[rec0,nzlearn]];


(* ::Text:: *)
(*and format this list of function as matrices.*)


(* de[x] will contain the matrix of differential equations with respect to x *)
invariants = {s12,s13,s23};
deall = ArrayReshape[rec,{Length[invariants],Length[masters],Length[masters]}];
Do[
  de[invariants[[i]]] = deall[[i]];
,{i,Length[{s12,s13,s23}]}];


(* Check some expected properties of the result *)
de[s12].de[s13]-de[s13].de[s12]//Together

D[de[s12],s13]-D[de[s13],s12]//Together



