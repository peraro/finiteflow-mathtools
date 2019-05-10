(* ::Package:: *)

BeginPackage["Symbols`",{"FiniteFlow`","FFUtils`"}]


Sym::usage="Sym[w1,w2,...] represents a symbol with letters w1,w2,..."
WLogDerivs::usage="WLogDerivs[alphabet,invariants,substitutions:{}] returns rules for the derivatives WDlogC of the letters listed in alphabet, in terms of the list of invariants.  The alphabet may be a list of symbolic letters, in which case the specified substitutions are used to cast it in terms of the invariants."
WDlogC::usage="WDlogC[W1,W2,x1,x2], where W1 and W2 are two letters and x1,x2 are Mandelstam invariants, represents the crossed-derivative D[Log[W1],x1]D[Log[W2],x2]-D[Log[W1],x2]D[Log[W2],x1]."
WDlogCLinearRelations::usage="WDlogCLinearRelations[alphabet,derivs,invariants], where derivs is the output of WLogDerivs, returns information about linear relations between the derivatives WLogDerivs built from the alphabet."
IntegrableSymbolsWeight::usage="IntegrableSymbolsWeight[n,alphabet,lower,invariants,linrels] returns a list of integrable symbols of weight n for the specified alphabet, where lower are a list of symbols of weight n-1, and linrels is the info returned by WDlogCLinearRelations about the linear relations between the derivatives WDlogC."


Begin["`Private`"]


WLogDerivs[w_,vars_,subst_:{}]:=Module[
  {dwdx,wi,xi,wpairs,xpairs,sol,zeroes},
  dwdx=Table[Together[(D[#,xi]/#)&[wi/.subst]],{wi,w},{xi,vars}];
  wpairs = Subsets[Range[Length[w]],{2}];
  xpairs = Subsets[Range[Length[vars]],{2}];
  sol = Join@@Table[WDlogC[w[[wi[[1]]]],w[[wi[[2]]]],vars[[xi[[1]]]],vars[[xi[[2]]]]]->Together[dwdx[[wi[[1]],xi[[1]]]]dwdx[[wi[[2]],xi[[2]]]]-dwdx[[wi[[2]],xi[[1]]]]dwdx[[wi[[1]],xi[[2]]]]],{wi,wpairs},{xi,xpairs}];
  sol
];


WDlogCLinearRelations[w_,dwlogc_,vars_]:=Module[
  {alldwlogs,correlations,indepdwlogc,wi,xi},
  alldwlogs=First/@dwlogc;
  correlations = FFCorrelateFunctionsWithRoots[dwlogc,vars];
  indepdwlogc=Complement[alldwlogs,First/@correlations];
  {indepdwlogc,
   Join[correlations,
        ((WDlogC[#[[1,2]],#[[1,1]],#[[1,3]],#[[1,4]]]->-#[[2]])&/@correlations),
        (WDlogC[#[[2]],#[[1]],#[[3]],#[[4]]]->-#)&/@indepdwlogc,
        Join@@Table[WDlogC[wi,wi,xi[[1]],xi[[2]]]->0,{wi,w},{xi,Subsets[vars,{2}]}]
   ]}
];


Options[IntegrableSymbolsWeight]={"Select"->Automatic,"Expand"->True};
IntegrableSymbolsWeight[w_,ls_,lowerw_,vars_,derivrel_,OptionsPattern[]]:=Module[
  {syms,ansatz,c,ccs,integrable,vpairs,eqs,sol,sols,indepcc,ccases,
   wpairs,dw,dwpairs,symv,eqv,ii,vv,ieq,ivar,taus,eqccs,
   dwindep,dwcorr,symsm2,mapsymv,cmap,odd,dwlist,makekprod},
  syms=Join@@Outer[#1/.SS_Sym:>Join[SS,Sym[#2]]&,lowerw,ls];
  If[!TrueQ[OptionValue["Select"]==Automatic], syms = Select[syms,OptionValue["Select"]];];
  wpairs = Union[#[[{w-1,w}]]&/@DeleteDuplicates[Cases[syms,_Sym,Infinity]]];
  ccs = (c/@Range[Length[syms]]);
  ansatz = ccs.syms;
  vpairs = vv@@#&/@Subsets[vars,{2}];
  {dwindep,dwcorr} = derivrel;
  mapsymv=Association[{}];
  symsm2=Union[#[[;;-2]]&/@DeleteDuplicates[Cases[lowerw,_Sym,Infinity]]];
  Table[mapsymv[symsm2[[symv]]]=symv,{symv,Length[symsm2]}];
  (*Table[mapsymv[symsm2[[symv]]]=SparseArray[{symv}\[Rule]1,{Length[symsm2]}],{symv,Length[symsm2]}];*)
  Table[mapsymv[dwindep[[ivar]]]=SparseArray[{ivar}->1,{Length[dwindep]}];,{ivar,Length[dwindep]}];
  dwlist[ll_]:=ll;
  dwlist[0]=SparseArray[{},{Length[dwindep]}];
  makekprod=Join[SparseArray[{},{(#1-1)*Length[dwindep]}],
                 #2,
                 SparseArray[{},{(Length[symsm2]-#1)*Length[dwindep]}]]&;
  (mapsymv[#[[1]]]=(#[[2]]/.WDlogC[a__]:>mapsymv[WDlogC[a]]);)&/@dwcorr;
  eqs = Normal@Select[Join@@(*Join@@*)Table[ansatz/.(SS_Sym:>makekprod(*Outer*)[mapsymv[SS[[;;-3]]],dwlist[mapsymv[WDlogC[SS[[w-1]],SS[[w]],ivar[[1]],ivar[[2]]]]]]),{ivar,vpairs}],!TrueQ[#==0]&];
  ClearAll[mapsymv,symsm2,dwindep,dwcorr,wpairs,syms];
  sol = FFSparseSolve[(#==0)&/@eqs,ccs,"VarsPattern"->(DeleteDuplicates[Cases[{#},_c,Infinity]]&),"Parameters"->{}];
  ClearAll[eqs];
  indepcc = Complement[ccs,#[[1]]&/@sol];
  (*ccases = Inner[Rule,indepcc,UnitVector[Length[indepcc],#],List]&/@Range[Length[indepcc]];*)
  cmap = ConstantArray[0,Length[ccs]];
  ccases = (cmap[[#]]&);
  (cmap[[indepcc[[#,1]]]]=SparseArray[{#}->1,Length[indepcc]];)&/@Range[Length[indepcc]];
  (cmap[[#[[1,1]]]]=(#[[2]]/.c->ccases))&/@sol;
  sol = Expand[ansatz /. c->ccases];
  ClearAll[syms,ansatz,c,ccs,integrable,vpairs,eqs,sols,indepcc,ccases,
           wpairs,dw,dwpairs,symv,eqv,ii,vv,ieq,ivar,taus,eqccs,
           dwindep,dwcorr,symsm2,mapsymv,cmap,odd,dwlist,makekprod];
  If[TrueQ[sol==0],Return[{}]];
  sol = Normal@Select[sol,!TrueQ[#==0]&];
  If[TrueQ[OptionValue["Expand"]], Scan[(sol[[#]] = Expand[sol[[#]]];)&,Range[Length[sol]]];];
  sol
  (*Select[(ansatz /. Dispatch[Join[#,sol/.Dispatch[#]]])&/@ccases,!TrueQ[#==0]&]*)  
];


End[] (* "`Private`" *)


EndPackage[] (* "Symbols`" *)
