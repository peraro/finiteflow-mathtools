(* ::Package:: *)

BeginPackage[ "LiteIBP`",{"LiteMomentum`","LiteRed`","Vectors`","Types`","FiniteFlow`","LiteIBPPerp`"}]


LIBPToLiteRed::usage = "LIBPToLiteRed[expr] converts expr to LiteRed format."
LIBPToLiteRedTopology::usage = "LIBPToLiteRedTopology[name] returns a topology in LiteRed format (note: it does not work if LIBPFamilyLite was used to define the family)."
LIBPFamily::usage = "LIBPFamily[family,propagators,loopmomenta,externalmomenta,momentumconservation,rules,invariants,invariantsdim] defines an integral family where:
* family: is a symbol for the family
* propagators: is a list of momenta flowing in the loop propagarors of the form {l,m2} indicating mp2[l]-m2, or {l,k,m2} indicating mp[l,k]-m2
* loopmomenta: is a list of symbols indicating the loop momenta
* externalmomenta: is a list of symbols indicating the external momenta
* momentumconservation: is a rule for momentum conservation
* rules: is a complete list of replacements rules for the scalar products among independent external momenta
* invariants: is a list of rules which defines the external invariants in terms of scalar products
* invariantsdim: is a list of rules replacing each invariant with its mass dimension"
LIBPEMs::usage = "LIBPEms[family] returns a list of rules which can be passed to LiteRed as EMs->LIBPEMs[family]"

LIBPPropagators::usage = "LIBPPropagators[family] returns the momenta defining the propagatos of family (note: it does not work if LIBPFamilyLite was used to define the family)."
LIBPDenoms::usage = "LIBPDenoms[family] returns rules defining the loop denominators of family."
LIBPLoopMomenta::usage = "LIBPLoopMomenta[family] returns a list of symbols indicating the loop momenta of family."
LIBPExtMomenta::usage = "LIBPExtMomenta[family] returns a list of symbols indicating the external momenta."
LIBPIndepExtMomenta::usage = "LIBPIndepExtMomenta[family] returns a list of symbols indicating the linearly independent external momenta with respect to momentum conservation."
LIBPMomCons::usage = "LIBPMomCons[family] returns the rule for momentum conservation of family."
LIBPIds::usage = "LIBPIds[family] returns the list of replacements rules for the scalar products among independent the external momenta of family."
LIBPInvariants::usage = "LIBPInvariants[family] returns the list of rules which defines the external invariants in terms of scalar products."
LIBPInvariantsDim::usage = "LIBPInvariantsDim[family] returns the mass dimension of the external invatiants of family (note: it does not work if LIBPFamilyLite was used to define the family)."
LIBPOrigFamily::usage = "For internal use.  LIBPOrigFamily[fam] returns the family which defined fam (it could be fam itself or its uncrossed version)."
LIBPIndepCrossings::usage = "LIBPIndepCrossings[family] returns a list of independent crossing between the external momenta of family."

LIBPSps::usage = "LIBPSps[family] returns a list of independent scalar products involving the loop momenta of family."
LIBPSpsToJ::usage = "LIBPSpsToJ[family] returns rules converting scalar products involvng the loop momenta of family into a linear combination of loop denominators."

LIBPSubsectorQ::usage = "LIBPSubsectorQ[js1,js2] returns True if js1 is a subsector of js2 with fewer loop propagators."
LIBPSectorOrSubsectorQ::usage = "LIBPSectorOrSubsectorQ[js1,js2] returns True if js1 == js2 or js1 is a subsector of js2."
LIBPMapSector::usage = "LIBPMapSector[sect], where sect is a mapped sector, returns the unique sector which sect is mapped to."
LIBPSectorDependencies::usage = "LIBPSectorDependencies[fam] lists all unique sectors and their dependencies."
LIBPFullSectorDependencies::usage = "\
LIBPFullSectorDependencies[family] recursively finds all sector dependencies for family.
LIBPFullSectorDependencies[{family1,family2,...}] recursively finds all sector dependencies for the families family, family2, etc..."

LIBPSectorOf::usage = "LIBPSectorOf[fam,a1,a2,...] or LIBPSectorOf[j[fam,a1,a2,...]] returns the sector js of the integral j[fam,a1,a2,...]."
LIBPGeneralizedSectorOf::usage = "LIBPGeneralizedSectorOf[fam,a1,a2,...] or LIBPGeneralizedSectorOf[j[fam,a1,a2,...]] returns the 'generalized' sector js of the specified integral with the same powers of propagators."

LIBPGenerateSeeds::usage = "LIBPGenerateSeeds[sect, deltar, smin, smax, constraints : {}] generates a list of seed integrals with deltar dots and rank from smin to smax, for sector sect.  In the last argument, one can define a list of additional constraints to be satisfied by the seed integrals."
LIBPUnsortedIntegralsIn::usage = "Returns a list of integrals in an expression, in no particular order."
LIBPIntegralsIn::usage = "Returns a list of integrals in an expression, sorted from higher to lower weight according to the default weight."

LIBPEliminateZeroSectors::usage = "Eliminate zero sectors in an expression."
LIBPZeroQ::usage = "LIBPZeroQ[sect] returns True if the argument is an integral of a vanishing sector."
LIBPMappedQ::usage = "LIBPMappedQ[sect] returns True if the argument is an integral of a mapped sector."

LIBPIntR::usage = "LIBPIntR[x], where x is an integral, returns the integer R for that integral."
LIBPIntT::usage = "LIBPIntT[x], where x is an integralor a sector, returns the integer T for that integral."
LIBPIntS::usage = "LIBPIntS[x], where x is an integral, returns the integer S for that integral."
LIBPIntDots::usage = "LIBPIntDots[x], where x is an integral, returns the number of dots of that integral, i.e. the integer R-T."

LIBPDefaultWeight::usage = "LIBPDefaultWeight[integral] returns the default weight for Feynman integrals.  SortBy[integrallist,LIBPDefaultWeight] sorts the integrals in integrallist from higher to lower weights."

LIBPProp::usage = "\
LIBPProp[k,m2] = mp[k]-m2
LIBPProp[k,q,m2] = mp[k,q]-m2"

LIBPPowers::usage = "\
LIBPPowers[tot,len] returns the lists of non-negative integers of length len whose total is equal to tot.
LIBPPowers[tot,len,vallist] is the same as LIBPPowers[tot,len], expect.
"
LIBPCombinePowers::usage = "For internal use."

LIBPRenormTotalRank::usage = "LIBPRenormTotalRank[sect] returns the maximum total rank for a numerator of the sector sect in a renormalizable theory."
LIBPSubLoopRanks::usage = "LIBPSubLoopRanks[sect] returns a list with the maximum ranks with respect to each each loop momentum, in a numerator of the sector sect in a renormalizable theory."
LIBPSelectRenormSubLoops::usage = "LIBPSelectRenormSubLoops[topo,{integrals...}] filters the list of integrals returning only the ones with renormalizable subloops."

LIBPSectorString::usage="LIBPSectorString[sect] returns a unique string identifying sector sect."

LIBPExtMappedQ::usage = "LIBPExtMappedQ[sect] returns True is the sector sect is mapped into another topology."

LIBPUnsortedSectorsIn::usage="LIBPUnsortedSectorsIn[expr] lists all the sectors for the integrals appearing in expr, in no particular order."
LIBPSectorsIn::usage="LIBPUnsortedSectorsIn[expr] lists all the sectors for the integrals appearing in expr, sorted by weight."
LIBPExtMapped::usage="LIBPExtMapped[sect], where sect is a sector mapped to an external topology, returns the unique sector which sect is mapped to."

LIBPComputeDerivatives::usage = "LIBPComputeDerivatives[family] computes the derivative operators for family.  This needs to be called once per family, before using LIBPDeriv."
LIBPDeriv::usage = "LIBPDeriv[expr,x], where expr is a linear combination of Feynman integrals and x is an invariant, computes the derivative of expr with respect to x.  Before using this function, you need to call LIBPComputeDerivatives[f] for each famly f appearing in expr."

LIBPFamilyLite::usage ="LIBPFamilyLite[family,loopdenominators,loopmomenta,externalmomenta,momentumconservation,rules,invariants] defines an integral family where:
* family: is a symbol for the family
* loopdenominators: is a list of loop denominators in terms of scalar products between loop momenta and external momenta
* loopmomenta: is a list of symbols indicating the loop momenta
* externalmomenta: is a list of symbols indicating the external momenta
* momentumconservation: is a rule for momentum conservation
* rules: is a complete list of replacements rules for the scalar products among independent external momenta
* invariants: is a list of rules which defines the external invariants in terms of scalar products"
LIBPEMs::usage = "LIBPEms[family] returns a list of rules which can be passed to LiteRed as EMs->LIBPEMs[family]"

LIBPFastGenIds::usage="LIBPFastGenIds[fam,GetSeeds] generates identities for the family fam, using the function GetSeeds in order to generate seed integrals.  Further options can be passed to control how the identities are generated.  The function GetSeeds will be called as GetSeeds[type][sector] for each identity type in {\"IBP\",\"LI\",\"SM\",\"Map\",\"ExtMap\"} and each sector for wich the identity type is relevant.  The identities are stored in files in the current directory."
LIBPSerializeFastIds::usage="LIBPSerializeFastIds[filelist,subst,vars] serializes the identities in the list filelist of mx files into JSON format.  The unknowns of the system are the Feynman integrals, and additional unknowns can be added via the \"ExtraInts\" option.  It applies the substitution subst to the identities (useful to set an invariant to 1).  The list vars are the free parameters in the identities (typically epsilon and the invarants which haven't been set to 1)."
LIBPGetAllInts::usage="LIBPGetAllInts[filelist] returns a list of integrals appearing in the expressions stored in the mathematica files listed in filelist.  The integrals are sorted by weight, from higher to lower."
LIBPWriteSystemJSON::usage="LIBPWriteSystemJSON[filelist,integrals,neededintegrals,params] writes a file \"system.json\" defining a system of equations stored in the files listed in filelist, with respected to the unknown integrals, where needeintegrals are the needed unknowns and params are the free parameters in the system.  This is best called after LIBPFastGenIds and LIBPSerializeFastIds."

LIBPFindSymmetries::usage = "LIBPFindSymmetries[args...] is equivalent to LiteRed's FindSymmetries[args...] except that it treats normal mappings and generalized mappings differently.  In particular, when the kinematics has non-trivial symmetries, this can yield a more complete set of identities for equivalent ranges of seed integrals."
LIBPUniqueSectorsNoEMs::usage = "LIBPUniqueSectorsNoEMs[family], after calling LIBPFindSymmetries[family,args...], yields a list of unique sectors without including symmetries of the kinematic invariants.  This is a superset of UniqueSectors[family]."

LIBPIdTypes::usage="A list of identity types."

LIBPUniqueUnmatchedSectors::usage="For internal use."
LIBPUniqueUnmatchedSectorsMaybeNoEms::usage="For internal use."

LIBPGenIds::usage="LIBPGenIds[type][sect,seeds] generates identities of the specified type, for the sector sect and using seed integrals seeds."
LIBPGenCutIds::usage="LIBPGenIds[type][sect,seeds] generates identities of the specified type, for the sector sect and using seed integrals seeds, on the cut of the sector sect."


Begin[ "`Private`"]


LIBPToLiteRed[expr_] := expr /.(mm[a_]:>a) /. (mp[a_,b_] :> sp[a,b]);


LIBPSectorOf[fam_,exponents___]:=js[fam,##]&@@(If[#>0,1,0]&/@{exponents});
LIBPSectorOf[j[fam_,exponents___]]:=LIBPSectorOf[fam,exponents];


LIBPGeneralizedSectorOf[fam_,exponents___]:=js[fam,##]&@@(If[#>0,#,0]&/@{exponents});
LIBPGeneralizedSectorOf[j[fam_,exponents___]]:=LIBPGeneralizedSectorOf[fam,exponents];


IdsToSPs[moms_,mmcons_,ids_]:=Module[{sps},
  Do[
    sp[moms[[i]],moms[[j]]]=mmp[moms[[i]],moms[[j]]]/.mmcons/.ids//Together;
  ,{i,Length[moms]},{j,Length[moms]}];
];


LIBPProp[p_,m_]:=mp2[p]-m;
LIBPProp[p_,q_,m_]:=mp[p,q]-m;


LIBPToLiteRedTopology[name_]:=(LIBPToLiteRed[LIBPProp@@#/.LIBPMomCons[name]/.LIBPIds[name]])&/@LIBPPropagators[name];
(*LIBPToLiteRedTopology[name_]:=(LIBPToLiteRed[#[[2]]&/@LIBPDenoms[name]]);*)


IsInvertibleInvMap[invmap_]:=Module[{inv2,sij,sol},
  inv2 = #[[2]]&/@invmap;
  sol = Solve[Table[sij[ii]==inv2[[ii]],{ii,Length[invmap]}],Variables[inv2]];
  TrueQ[Length[sol]==1]
]


LIBPFamilyLite[fam_Symbol, props_,
               loopmoms_, extmoms_, momconsin_,
               ids_, invsdefs_]:=Module[{momcons, allcrossings,invlist,extmasses},
  LIBPOrigFamily[fam] = fam;
  LIBPPropagators[fam] = $Failed;
  LIBPExtMomenta[fam] = extmoms;
  LIBPLoopMomenta[fam] = loopmoms;
  momcons = If[TrueQ[momconsin[[0]]==List],momconsin[[1]], momconsin];
  LIBPIndepExtMomenta[fam] = Select[extmoms,FreeQ[momcons[[1]],#]&];
  LIBPMomCons[fam] = momcons;
  LIBPIds[fam] = ids;
  LIBPInvariants[fam] = invsdefs;
  LIBPInvariantsDim[fam] = $Failed;
  LIBPDenoms[fam] = Table[j[fam,Sequence@@(-UnitVector[Length[props],i])]->Together[(props[[i]])/.ids],{i,1,Length[props]}];

  invlist = #[[1]]&/@invsdefs;
  extmasses = Together[(mmp2/@extmoms)/.momcons/.ids];
  allcrossings = MapThread[(#1->#2) &,{extmoms,#}]&/@Permutations[extmoms];
  allcrossings = Select[allcrossings,AllTrue[Together[((mmp2/@extmoms)/.#/.momcons/.ids)-extmasses],(#==0)&]&];
  LIBPIndepCrossings[fam] = Select[Normal@GroupBy[allcrossings, Inner[Rule,invlist,(Together[(invlist /. invsdefs /. # /. momcons /. ids)]),List]&],IsInvertibleInvMap[#[[1]]]&];
  LIBPEMs[fam] = Map[LIBPToLiteRed[mm[#[[1]]]->mm[#[[2]]]/.momcons]&/@#&,(Select[#,FreeQ[momcons[[1]],#[[1]]]&]&/@(((#->#)&/@invlist) /. LIBPIndepCrossings[fam]))];

  (* Some LiteRed stuff *)
  Declare[Evaluate[Join[extmoms,loopmoms]],Vector,Evaluate[invlist],Number];
  IdsToSPs[extmoms,momcons,ids];
];


LIBPFamily[fam_Symbol, props_,
           loopmoms_, extmoms_, momconsin_,
           ids_, invsdefs_, invsdims_]:=Module[{momcons, allcrossings,invlist,extmasses},
  LIBPOrigFamily[fam] = fam;
  LIBPPropagators[fam] = props;
  LIBPExtMomenta[fam] = extmoms;
  LIBPLoopMomenta[fam] = loopmoms;
  momcons = If[TrueQ[momconsin[[0]]==List],momconsin[[1]], momconsin];
  LIBPIndepExtMomenta[fam] = Select[extmoms,FreeQ[momcons[[1]],#]&];
  LIBPMomCons[fam] = momcons;
  LIBPIds[fam] = ids;
  LIBPInvariants[fam] = invsdefs;
  LIBPInvariantsDim[fam] = invsdims;
  LIBPDenoms[fam] = Table[j[fam,Sequence@@(-UnitVector[Length[props],i])]->Together[LIBPProp@@(props[[i]])/.ids],{i,1,Length[props]}];

  invlist = #[[1]]&/@invsdefs;
  extmasses = Together[(mmp2/@extmoms)/.momcons/.ids];
  allcrossings = MapThread[(#1->#2) &,{extmoms,#}]&/@Permutations[extmoms];
  allcrossings = Select[allcrossings,AllTrue[Together[((mmp2/@extmoms)/.#/.momcons/.ids)-extmasses],(#==0)&]&];
  LIBPIndepCrossings[fam] = Select[Normal@GroupBy[allcrossings, Inner[Rule,invlist,(Together[(invlist /. invsdefs /. # /. momcons /. ids)]),List]&],IsInvertibleInvMap[#[[1]]]&];
  LIBPEMs[fam] = Map[LIBPToLiteRed[mm[#[[1]]]->mm[#[[2]]]/.momcons]&/@#&,(Select[#,FreeQ[momcons[[1]],#[[1]]]&]&/@(((#->#)&/@invlist) /. LIBPIndepCrossings[fam]))];

  (* Some LiteRed stuff *)
  Declare[Evaluate[Join[extmoms,loopmoms]],Vector,Evaluate[invlist],Number];
  IdsToSPs[extmoms,momcons,ids];
];


LIBPSps[fam_]:= Join[DeleteDuplicates[Flatten[Outer[mmp[#1,#2]&,LIBPLoopMomenta[fam],LIBPLoopMomenta[fam]]]],
                     Flatten[Outer[mmp[#1,#2]&,LIBPLoopMomenta[fam],LIBPIndepExtMomenta[fam]]]];


LIBPSpsToJ[fam_]:=LIBPSpsToJ[fam,Sequence@@ConstantArray[1,Length[LIBPDenoms[fam]]]];
LIBPSpsToJ[fam_,dens___]:=Module[{densidx,thisdens},
  densidx = Flatten[Position[{dens},1]];
  thisdens = LIBPDenoms[fam][[densidx]];
  FFDenseSolve[(#[[2]]-#[[1]]==0)&/@thisdens,Join[LIBPSps[fam],#[[1]]&/@thisdens]]
]


LIBPSubsectorQ[js[fam1_,exp1___],js[fam1_,exp2___]]:=({exp1}!={exp2}) && AllTrue[{exp2}-{exp1},(#>=0)&];


LIBPSectorOrSubsectorQ[js[fam1_,exp1___],js[fam1_,exp2___]]:=AllTrue[{exp2}-{exp1},(#>=0)&];


LIBPMapSector[js[fam1_,exp1__]]:=If[MemberQ[UniqueSectors[fam1],js[fam1,exp1]],js[fam1,exp1], js@@Together[j[fam1,exp1]/.jRules[fam1,exp1]]];


LIBPSectorDependencies[fam_]:=Module[{unique,nonzero},
  unique = UniqueSectors[fam];
  nonzero = NonZeroSectors[fam];
  Table[unique[[i]]->DeleteDuplicates[LIBPMapSector/@Select[nonzero,LIBPSubsectorQ[#,unique[[i]]]&]],{i,Length[unique]}]
];


LIBPFullSectorDependencies[fams_List]:=Module[{nonzero,depof,sec,done},
  depof[_]:={};
  done=Association[{}];
  nonzero = SortBy[Join@@(NonZeroSectors/@fams),{LIBPIntT[#],LIBPExtMappedQ[#],LIBPMappedQ[#]}&];
  Do[
    depof[sec]=Union[#,Sequence@@(depof/@#)]&[Select[Keys[done],(LIBPSectorOrSubsectorQ[#,sec]||(LIBPExtMappedQ[sec]&&LIBPExtMapped[sec]==#)||LIBPMapSector[sec]==#)&]];
    done[sec] = True;
  ,{sec,nonzero}];
  #->depof[#]&/@nonzero
];
LIBPFullSectorDependencies[fam_]:=LIBPFullSectorDependencies[{fam}];


LIBPPowers[total_,n_,range_]:=Join@@(DeleteDuplicates[Permutations[#]]& /@ IntegerPartitions[total,{n},range]);
LIBPPowers[total_,n_]:=LIBPPowers[total,n,Range[0,total]];
LIBPCombinePowers[numidx_,numai_,denidx_,denai_]:=Module[
  {res},
  res = ConstantArray[0,Length[numai]+Length[denai]];
  res[[numidx]] = numai;
  res[[denidx]] = denai;
  res
];


LIBPGenerateSeeds[js[fam_Symbol,ai__],{drmin_,drmax_},smin_,smax_,filters_:{}]:=Module[{sectai,spsidx,denidx,spows,rpows,aipows,filter},
  sectai = {ai};
  spsidx = Flatten[Position[sectai,0]];
  denidx = Flatten[Position[sectai,1]];
  spows = -Join@@(LIBPPowers[#,Length[spsidx]]&/@Range[smin,smax]);
  rpows = Join@@(LIBPPowers[#,Length[denidx]]&/@Range[drmin,drmax]);
  aipows = (#+sectai)&/@(Join@@Table[LIBPCombinePowers[spsidx,sp,denidx,rp],{sp,spows},{rp,rpows}]);
  filter = filters;
  If[!TrueQ[filter[[0]]==List], filter = {filter}];
  aipows = Select[aipows,(And@@(Table[filter[[i]][IBPGI[fam,#]],{i,Length[filter]}]))&];
  j[fam,Sequence@@#]&/@aipows
];
LIBPGenerateSeeds[js[fam_Symbol,ai__],drmax_,smin_,smax_,filters_:{}]:=LIBPGenerateSeeds[js[fam,ai],{0,drmax},smin,smax,filters];


LIBPOrderedQ[j[fam1_,a1in__],j[fam2_,a2in__]]:=Module[{a1,a2,t1,t2,r1,r2,s1,s2},
  a1 = {a1in};
  a2 = {a2in};
  (* prefer subtopologies *)
  t1 = Length[a1];
  t2 = Length[a2];
  If[t1>t2,Return[True]];
  If[t2>t1,Return[False]];
  t1 = Length[Select[a1,(#>0)&]];
  t2 = Length[Select[a2,(#>0)&]];
  If[t1>t2,Return[True]];
  If[t2>t1,Return[False]];
  (* prefer no quadratic denoms *)
  r1 = Plus@@Select[a1,(#>0)&];
  r2 = Plus@@Select[a2,(#>0)&];
  If[r1>r2,Return[True]];
  If[r2>r1,Return[False]];
  (* prefer lower rank *)
  s1 = -Plus@@Select[a1,(#<0)&];
  s2 = -Plus@@Select[a2,(#<0)&];
  If[s1>s2,Return[True]];
  If[s2>s1,Return[False]];
  (* with equal rank, prefer distributing the rank across several scalar products *)
  s1 = Max[Select[-a1,(#>0)&]];
  s2 = Max[Select[-a2,(#>0)&]];
  If[s1>s2,Return[True]];
  If[s2>s1,Return[False]];
  (* otherwise let mathematica choose *)
  OrderedQ[j[fam1,a1in],j[fam2,a2in]]
];

LIBPDefaultWeight[j[fam1_,a1in__]]:=Module[{a1},
  a1 = {a1in};
  (* prefer subtopologies *)
  {(-Plus@@# + Length[#])&[Select[a1,(#>0)&]],
   -Length[a1],
   -Length[Select[a1,(#>0)&]],
   -Plus@@Select[a1,(#>0)&],
   Plus@@Select[a1,(#<0)&],
   !LIBPExtMappedQ[j[fam1,a1in]],
   !LIBPMappedQ[j[fam1,a1in]],
   -Max[Select[-a1,(#>0)&]],
   j[fam1,a1in]}
];

LIBPSortByWeight=LIBPDefaultWeight;


LIBPUniqueUnmatchedSectors[fam_]:=UniqueSectors[fam];


LIBPUniqueUnmatchedSectorsMaybeNoEms[fam_]:=LIBPUniqueSectorsNoEMs[fam];


LIBPUnsortedIntegralsIn[expr_]:=DeleteDuplicates[Cases[{expr},_j,Infinity]];
(*LIBPIntegralsIn[expr_]:=Sort[LIBPUnsortedIntegralsIn[expr],LIBPOrderedQ];*)
LIBPIntegralsIn[expr_]:=SortBy[LIBPUnsortedIntegralsIn[expr],LIBPSortByWeight];


LIBPZeroQ[j[fam_Symbol,a__]]:=MemberQ[ZeroSectors[fam],LIBPSectorOf[j[fam,a]]];
LIBPZeroQ[j[fam_Symbol[map_],a__]]:=MemberQ[ZeroSectors[fam],LIBPSectorOf[j[fam,a]]];
LIBPMappedQ[j[fam_Symbol,a__]]:=MemberQ[MappedSectors[fam],LIBPSectorOf[j[fam,a]]];
LIBPExtMappedQ[j[fam_Symbol,a__]]:=MemberQ[ExtMappedSectors[fam],LIBPSectorOf[j[fam,a]]];
LIBPZeroQ[js[fam_Symbol,a__]]:=MemberQ[ZeroSectors[fam],js[fam,a]];
LIBPMappedQ[js[fam_Symbol,a__]]:=MemberQ[MappedSectors[fam],js[fam,a]];
LIBPExtMappedQ[js[fam_Symbol,a__]]:=MemberQ[ExtMappedSectors[fam],js[fam,a]];
LIBPEliminateZeroSectors[expr_]:=expr /. Dispatch[(#->0)&/@Select[LIBPUnsortedIntegralsIn[expr],LIBPZeroQ]];


LIBPIntR[j[fam_,a__]]:=Plus@@(Select[{a},#>0&]);
LIBPIntT[j[fam_,a__]]:=Length[(Select[{a},#>0&])];
LIBPIntT[js[fam_,a__]]:=Length[(Select[{a},#>0&])];
LIBPIntS[j[fam_,a__]]:=-Plus@@(Select[{a},#<0&]);
LIBPIntDots[int_j]:=LIBPIntR[int]-LIBPIntT[int];


LIBPRenormTotalRank[js[fam_,dens__]]:=Plus@@(Select[{dens},#>0&])-Length[LIBPLoopMomenta[fam]]+1; (* <-- V = P-L+1 *)
LIBPSubLoopRanks[js[fam_,dens__]]:=Module[
    {densidx,numidx,numldep},
    numidx = Flatten[Position[{dens},0]];
    numldep = Table[Count[LIBPPropagators[fam][[ii]],kk,Infinity],{kk,LIBPLoopMomenta[fam]},{ii,numidx}];
    -numldep.List@@(#[[2;;]][[numidx]])&
];


LIBPSelectRenormSubLoops[topo_,seeds_List]:=Module[
  {llpos,subrenorm},
  llpos = Table[Flatten[Position[LIBPDenoms[topo],_?(!FreeQ[#,kkk]&),1]],{kkk,LIBPLoopMomenta[topo]}];
  subrenorm[j[topo,a__]]:=And@@((Plus@@({a}[[#]])>=0)&/@llpos);
  Select[seeds,subrenorm]
];


LIBPSectorString[sect_js]:=StringJoin@@(ToString/@sect);


(* Fast/efficient IBP generation:
   - sector mappings are added as identities
   - identities for sectors with many seeds are split into subsets
*)


LIBPMakeListPartitions[list_,length_]:=If[TrueQ[Mod[Length[list],length]==0],
                                          Partition[list,length],
                                          Join[Partition[list,length],{list[[-Mod[Length[list],length];;]]}]
                                       ];


LIBPIdTypes={"IBP","LI","SR","Map","ExtMap","OrtInt"} ;(* <-- types of identities *)


LIBPMaxIdsPerFile=100000;


LIBPGenIds["IBP"][sect_,seeds_]:=LIBPEliminateZeroSectors[Join@@((IBP[sect[[1]]]@@(#[[2;;]]))&/@seeds)];
LIBPGenIds["LI"][sect_,seeds_]:=LIBPEliminateZeroSectors[Join@@((LI[sect[[1]]]@@(#[[2;;]]))&/@seeds)];
LIBPGenIds["SR"][sect_,seeds_]:=LIBPEliminateZeroSectors[Join@@((SR[sect[[1]]]@@(#[[2;;]]))&/@seeds)];
LIBPGenIds["Map"][sect_,seeds_]:=LIBPEliminateZeroSectors[-seeds+(seeds/.(jRules@@sect))];
LIBPGenIds["ExtMap"][sect_,seeds_]:=LIBPEliminateZeroSectors[-seeds+(seeds/.(jExtRules@@sect))];
LIBPGenIds["OrtInt"][sect_,seeds_]:=LIBPEliminateZeroSectors[(#[[2]]-#[[1]])&/@LIBPIntegrateOut[seeds,sect]];


LIBPGenCutIds[type_][sector_,seeds_]:=(#/.Dispatch[(#->0)&/@Select[LIBPIntegralsIn[#],(LIBPIntT[#]<LIBPIntT[sector])&]])&@(LIBPGenIds[type][sector,seeds]);


LIBPLaunchKernels[n_]:=(LaunchKernels[2]; LIBPLaunchKernels[n-2];);
LIBPLaunchKernels[1]:=LaunchKernels[1];
LIBPLaunchKernels[0]:=Null;
LIBPLaunchKernels[Automatic]:=LaunchKernels[];


Options[LIBPFastGenIds]={"GetSectors"->Automatic,"MaxIdsPerFile"->LIBPMaxIdsPerFile,"GenIds"->LIBPGenIds,"Directory"->Automatic,"LaunchKernels"->0};
LIBPFastGenIds[fam_,GetSeeds_,OptionsPattern[]]:=Module[
  {GetSectors,NIds,IdSeeds,maxidsperfile,GenIds,idtype,dir},
  dir = If[TrueQ[#==Automatic],Directory[],Quiet[CreateDirectory[#],CreateDirectory::filex];#]&@(OptionValue["Directory"]);
  GetSectors = OptionValue["GetSectors"];
  GenIds=OptionValue["GenIds"];
  maxidsperfile=OptionValue["MaxIdsPerFile"];
  If[TrueQ[GetSectors==Automatic],
    Clear[GetSectors];
    GetSectors["IBP"]=LIBPUniqueUnmatchedSectorsMaybeNoEms[fam];
    GetSectors["LI"]=LIBPUniqueUnmatchedSectorsMaybeNoEms[fam];
    GetSectors["SR"]=LIBPUniqueUnmatchedSectorsMaybeNoEms[fam];
    GetSectors["Map"]=MappedSectors[fam];
    GetSectors["ExtMap"]=If[TrueQ[#[[0]]==List],#,{}]&[ExtMappedSectors[fam]];
    GetSectors["OrtInt"]=Complement[LIBPUniqueSectorsNoEMs[fam],LIBPUniqueUnmatchedSectorsMaybeNoEms[fam]];
  ];
  (NIds["IBP"][#]=Length[IBP[#[[1]]][[1]]];)&/@GetSectors["IBP"];
  (NIds["LI"][#]=Length[LI[#[[1]]][[1]]];)&/@GetSectors["LI"];
  (NIds["SR"][#]=Length[jSymmetries@@(#)];)&/@GetSectors["SR"];
  NIds["Map"][a_]:=1;
  NIds["ExtMap"][a_]:=1;
  NIds["OrtInt"][a_]:=1;
  LIBPLaunchKernels[OptionValue["LaunchKernels"]];
  Do[If[Length[GetSectors[idtype]]>0,
    (Print["Seeding id type ",idtype];);IdSeeds[idtype] = Join@@ParallelTable[If[NIds[idtype][sect]>0, Module[{seeds,filenames,ii},
      SetDirectory[dir];
      seeds={sect,#}&/@LIBPMakeListPartitions[GetSeeds[idtype][sect],Max[Quotient[maxidsperfile,NIds[idtype][sect]],1]];
      filenames = Table["seeds_"<>ToString[idtype]<>"_"<>LIBPSectorString[sect]<>"_"<>ToString[ii]<>".mx",{ii,Length[seeds]}];
      Do[Export[filenames[[ii]],seeds[[ii]],"MX"],{ii,Length[seeds]}];
      ResetDirectory[];
      filenames
     ],{}]
    ,{sect,GetSectors[idtype]},DistributedContexts->Automatic];
  ],{idtype,LIBPIdTypes}];
  Do[If[TrueQ[IdSeeds[idtype][[0]]==List] && Length[IdSeeds[idtype]]>0,
   Print["=========================="];
   Print["Id type: ",idtype, " # chunks = ",Length[IdSeeds[idtype]]];
   Print["=========================="];
   ParallelDo[
     SetDirectory[dir];
     Print[idtype," #",ii," of ",Length[IdSeeds[idtype]]];
     Print["-> ",idtype," #",ii," done in ", AbsoluteTiming[Module[{ids,sect,seeds,idseeds},
       idseeds = Import[IdSeeds[idtype][[ii]]];
       sect = idseeds[[1]];
       seeds = idseeds[[2]];
       ids = (GenIds[idtype][sect,seeds]);
       Export["ids_"<>ToString[sect[[1]]]<>"_"<>ToString[idtype]<>ToString[ii]<>".mx",ids,"MX"];
       Export["ints_"<>ToString[sect[[1]]]<>"_"<>ToString[idtype]<>ToString[ii]<>".mx",LIBPUnsortedIntegralsIn[ids],"MX"];
     ];]];
     ResetDirectory[];
   ,{ii,Length[IdSeeds[idtype]]},DistributedContexts->Automatic]];
  ,{idtype,LIBPIdTypes}];
];


Options[LIBPGetAllInts]={"IntegralWeight"->LIBPSortByWeight};
LIBPGetAllInts[intsfiles_,OptionsPattern[]]:=Module[{ints,file},
  (* This is okay, but uses too much memory
  SortBy[DeleteDuplicates[Join@@(Import/@intsfiles)],OptionValue["IntegralWeight"]];
  *)
  ints = Association[{}];
  Do[
    (ints[#]=1;)&/@Import[file];
  ,{file,intsfiles}];
  ints = Keys[ints];
  SortBy[ints,OptionValue["IntegralWeight"]]
];


Options[LIBPSerializeFastIds]:={"IntegralWeight"->LIBPSortByWeight,"ExtraInts"->{},"Overwrite"->True,"LaunchKernels"->0};
LIBPSerializeFastIds[todofilesin_,toonesubst_,pars_,OptionsPattern[]]:=Module[
  {todofiles,allints,IntegralWeight,position,ii,eqsfiles,extraints,nints},
  todofiles = todofilesin;
  If[!TrueQ[OptionValue["Overwrite"]], todofiles = Select[todofiles,!(FileExistsQ[StringReplace[#,{"ids_"->"posints_"}]]
                                                                   && FileExistsQ[StringReplace[#,{"ids_"->"sids_",".mx"->".json"}]])&]];
  IntegralWeight=OptionValue["IntegralWeight"];
  extraints = OptionValue["ExtraInts"];
  Print["Generating list of integrals..."];
  (*allints=Join@@(Import/@(StringReplace[#,{"ids_"->"ints_"}]&/@todofiles));
  allints=Join[SortBy[DeleteDuplicates[allints],IntegralWeight],extraints];*)
  allints = Join[LIBPGetAllInts[StringReplace[#,{"ids_"->"ints_"}]&/@todofilesin,
                                 "IntegralWeight"->IntegralWeight],extraints];
  Print["-> list of integrals generated"];
  position = Association[{}];
  Table[position[allints[[ii]]]=ii-1;,{ii,Length[allints]}];
  nints = Length[allints];
  ClearAll[allints];
  Print["Export position maps..."];
  Do[
  Export[StringReplace[file,{"ids_"->"posints_"}],
         ((#/.toonesubst)->position[#])&/@Join[Import[StringReplace[file,{"ids_"->"ints_"}]],extraints],
         "MX"];
  ,{file,todofiles}];
  Print["-> position maps exported"];
  ClearAll[position];
  LIBPLaunchKernels[OptionValue["LaunchKernels"]];
  Print["Total: ",AbsoluteTiming[ParallelDo[
  Print["File ",file];
  Print["File ",file," done in ",AbsoluteTiming[Module[{ids,integralsin,thisposition},
    thisposition = Association[{}];
    (thisposition[#[[1]]]=#[[2]])&/@Import[StringReplace[file,{"ids_"->"posints_"}]];
    integralsin = Join[LIBPUnsortedIntegralsIn[#],extraints]&;
    ids=((#/.toonesubst)==0)&/@Import[file];
    FFSparseEqsToJSON[StringReplace[file,{"ids_"->"sids_",".mx"->".json"}],
                      pars,ids,nints,integralsin,thisposition];
    Export[StringReplace[file,{"ids_"->"nids_",".mx"->".m"}],Length[ids]];
   ]]];,
   {file,todofiles},DistributedContexts->Automatic];]];
];


Options[LIBPWriteSystemJSON]={"FileName"->"system.json"}
LIBPWriteSystemJSON[eqsfiles_,allints_,needed_,pars_,OptionsPattern[]]:=Module[{neqs},
  neqs = Plus@@(Import/@(StringReplace[#,{"sids_"->"nids_",".json"->".m"}]&/@eqsfiles));
  Print["neqs = ",neqs];
  FFSparseSystemToJSON[OptionValue["FileName"],neqs,allints,pars,eqsfiles,"NeededVars"->needed]
];


LIBPUniqueSectorsNoEMs[fam_]:=UniqueSectors[fam];


Options[LIBPFindSymmetries]=Options[FindSymmetries];
LIBPFindSymmetries[fam_,opt:OptionsPattern[]]:=(
  If[TrueQ[Length[OptionValue[EMs]]>1],
    FindSymmetries[fam,Sequence@@FilterRules[{opt},Select[Options[FindSymmetries],FreeQ[#,EMs]&]]];
    LIBPUniqueSectorsNoEMs[fam]=UniqueSectors[fam];
  ];
  FindSymmetries[fam,opt];
);


LIBPUnsortedSectorsIn[expr_]:=First/@Normal[GroupBy[LIBPUnsortedIntegralsIn[expr],LIBPSectorOf]];
LIBPSectorsIn[expr_]:=SortBy[LIBPUnsortedSectorsIn[expr],LiteIBP`Private`LIBPSortByWeight[j@@#]&];


LIBPExtMapped[js[fam_,a__]]:=js@@(j[fam,a]/.jExtRules[fam,a]);
LIBPExtMapped[fam_,a__]:=LIBPExtMapped[js[fam,a]];


(* Derivatives of Feynman integrals *)


LIBPGetDerivatives[topo_]:=Module[
  {invs, invl, ids, pi,derivs,eqs,vars,pij,inveqs,mysys,learn},
  invs = LIBPInvariants[topo];
  invl = #[[1]]&/@invs;
  ids = LIBPIds[topo];
  pi = LIBPIndepExtMomenta[topo];
  derivs = Table[-mm[LIBPDerivV[pp]][mu]+ mm[pp][mu]LIBPDerivV[mmp2[pp]]+Sum[mm[qq][mu]LIBPDerivV[mmp[pp,qq]],{qq,pi}],{pp,pi}];
  eqs=Collect[Flatten[Table[LContract[mm[pp][mu] ddd],{ddd,derivs},{pp,pi}]],_LIBPDerivV,Expand[#/.ids]&];
  pij = Flatten[Table[mmp[pi[[i]],pi[[j]]],{i,Length[pi]},{j,i,Length[pi]}]];
  inveqs = Table[LIBPDerivV[xk]-Sum[D[ppij/.ids,xk] LIBPDerivV[ppij],{ppij,pij}],{xk,invl}];
  vars=Join[LIBPDerivV/@pij,Flatten[Table[mmp[LIBPDerivV[pp],qq],{pp,pi},{qq,pi}]]];
  FFDenseSolve[(#==0)&/@Join[inveqs,eqs],Join[LIBPDerivV/@invl,vars],"NeededVars"->(LIBPDerivV/@invl)](*;
  learn=FFDenseSolverLearn[mysys,Join[LIBPDerivV/@invl,vars]];
  FFDenseSolverSol[FFReconstructFunction[mysys,invl],learn]*)
  (*FFDenseSolve[(#\[Equal]0)&/@eqs,vars]*)
];


LIBPDerivivative[topo_,expr_,inv_,derivs_]:=Module[{mu},
  Together[((LIBPDerivV[inv]/.derivs)/.mp[mm[a_],mm[LIBPDerivV[b_]]]:>LContract[mm[a][mu]MomDerivative[expr/.LIBPInvariants[topo],mm[b][mu]]])/.LIBPIds[topo]]
];


LIBPComputeDerivatives[topo_]:=(LIBPDerivatives[topo]=LIBPGetDerivatives[topo];
                                Table[LIBPDenomsDeriv[topo,sss]=(Collect[LIBPDerivivative[topo,#,sss,LIBPDerivatives[topo]]/.LIBPSpsToJ[topo],_j,Together]&/@(#[[2]]&/@LIBPDenoms[topo]));,{sss,First/@LIBPInvariants[topo]}];
                                );


LIBPDeriv[j[t_Symbol,a__],s_]:=LIBPDenomsDeriv[t,s].((-{a}[[#]] j[t,Sequence@@(UnitVector[Length[{a}],#]+{a})])&/@Range[Length[{a}]]);


LIBPDeriv[j[t_[invs_List],a__],s_]:=(((LIBPDeriv[j[t,a],#]&/@(First/@LIBPInvariants[t]))/.Dispatch[Inner[Rule,First/@LIBPInvariants[t],invs,List]]).(D[#,s]&/@invs))/.t->t[invs];
LIBPDeriv[a_Plus,s_]:=LIBPDeriv[#,s]&/@a;
LIBPDeriv[a_Times,s_]:=Plus@@(MapAt[LIBPDeriv[#,s]&,a,#]&/@Range[Length[a]]);
LIBPDeriv[expr_,s_]:=D[expr,s];


End[] (* `Private` *)


EndPackage[]
