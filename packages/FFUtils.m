(* ::Package:: *)

BeginPackage["FFUtils`",{"FiniteFlow`"}]


FFCorrelateFunctions::usage="FFCorrelateFunctions[{f1->expr1,f2->expr2,...},vars] finds linear relations between f1,f2,... assuming they have equal to expr1,expr2,... as rational functions of the variables vars with rational coefficients.  The functions must be listed from lower to higher weight."
FFCorrelateFunctionsWithRoots::usage="FFCorrelateFunctionsWithRoots[{f1->expr1,f2->expr2,...},vars] FFCorrelateFunctions[{f1->expr1,f2->expr2,...},vars] finds linear relations between f1,f2,... assuming they have equal to expr1,expr2,...  The expressions expr1,expr2,... must be rational functions in the variables vars, with rational coefficients, and in a number of square roots of rational functions of vars.  The functions must be listed from lower to higher weight."
FFLinearizeRoots::usage="FFLinearizeRoots[expr], where expr is an expression containing square roots, brings expr into a standard form which is multi-linear in the square roots.  Each square root in the output is written as FFSqrt[x]."
FFSqrt::usage="FFSqrt[x] represents the square root of x, in the output of FFLinearizeRoots."
FFLinearRelationsFromFit::usage = "FFFindLinearRelationsFromFit[{f1,f2,...},{c1,c2,...},sol], where sol is the solution of the linear fit problem c1*f1+c2*f2+...=0 with respect to the unknowns c1,c2,..., returns a list of linear relations satisfied by the functions {f1,f2,...}.  Note the the functions f1,f2,... passed as input can be symbolic entries."


Begin["`Private`"]


FFCorrelateFunctions[funcs_,vars_]:=Module[
  {nonzeroes,exprmap,fun,el,funmap,keys,cc,ccs,ls,indepccs,sol,cfun},
  nonzeroes = Select[funcs,!TrueQ[#[[2]]==0]&];
  funmap=Association[{}];
 
  (* Detect zeroes *)
  (funmap[#[[1]]]=0)&/@Select[funcs,TrueQ[#[[2]]==0]&];
  
  (* Detect identical functions up to a sign *)
  exprmap=Association[{}];
  Do[
    el = exprmap[fun[[2]]];
    If[!TrueQ[el[[0]]==Missing],
      funmap[fun[[1]]] = el;
      Continue[];
    ];
    el = exprmap[-fun[[2]]];
    If[!TrueQ[el[[0]]==Missing],
      funmap[fun[[1]]] = -el;
      Continue[];
    ];
    exprmap[fun[[2]]] = fun[[1]]; 
  ,{fun,nonzeroes}];
  
  (* find linear relations btw remaining functions *)
  keys = Keys[exprmap];
  ccs = (cc[exprmap[#]])&/@keys;
  cfun = First/@ccs;
  ls = FFLinearFit[{},Join[keys,{0}],{},vars,ccs];
  indepccs = Complement[ccs,First/@ls];
  ls = ccs /.Dispatch[ls];
  sol = Solve[((ls/.#->1/.cc[_]->0).cfun)==0,#[[1]]][[1,1]]&/@indepccs;
  
  (* update the other solutions *)
  SortBy[Join[sol, Normal@funmap /. Dispatch[sol]],Position[First/@funcs,#[[1]]]&]
];


FFLinearizeRoots[func_]:=Module[{roots,rootsyms,root,tmp,r},
  roots = Union[Cases[{func},_^(1/2),Infinity],1/Cases[{func},_^(-1/2),Infinity]];
  rootsyms = Unique[r]&/@Range[Length[roots]];
  tmp = func /. Join[Inner[Rule,roots,rootsyms,List],Inner[Rule,1/roots,1/rootsyms,List]];
  Do[
    tmp = PolynomialRemainder[tmp,roots[[ii]]^2-rootsyms[[ii]]^2,rootsyms[[ii]]];
  ,{ii,Length[roots]}];
  tmp /. Inner[Rule,rootsyms,FFSqrt/@(roots^2),List]
];


FFCorrelateFunctionsWithRoots[funcs_,vars_]:=Module[{roots,lrfuncs,rootrules},
  lrfuncs = FFLinearizeRoots/@(#[[2]]&/@funcs);
  roots = Union[Cases[{lrfuncs},_FFSqrt,Infinity]];
  FFCorrelateFunctions[Inner[Rule,First/@funcs,Together[lrfuncs],List],Join[vars,roots]]
];


FFLinearRelationsFromFit[funcs_,coeffsin_,sol_]:=Module[
  {indepv,dim,var,depcs,coeffs,deppos},
  coeffs = coeffsin /. Dispatch[sol];
  indepv = Complement[coeffsin,First/@sol];
  dim = Length[indepv];
  Table[
    deppos = Position[coeffsin,var][[1,1]];
    depcs = (coeffs /. var->1) /. Dispatch[#->0&/@indepv];
    depcs[[deppos]]=0;
    funcs[[deppos]] -> (-depcs).funcs
  ,{var,indepv}]
];


End[] (* "`Private`" *)


EndPackage[] (* "FFUtils`" *)
